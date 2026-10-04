# lcm_ord2 — paramagnetic B-site order–disorder thermodynamics (triage item 1)

Job prefix `lcm/ord2/`; job scripts `infra/jobs_src/lcm_ord2/` (ord2.py, tasks_*.json) and `infra/jobs_src/lcm_ord2_vc/`;
analysis + results here (`tracks/lcm/lcm_ord2/`).

## Setup (2026-10-02)
- Parent YBaMnFeO5 cell: rock-salt P4/n, 110-Ry vc-relaxed (jobs/lcm/h110/YBMFO_P4n; a = 5.6649, c = 7.6888 A) -> `YBMFO_P4n_h110.json`.
- Protocol (both compounds): fixed parent cell; ions-only relax (FCONV 2e-3 Ry/bohr) in a random species-uncorrelated spin
  config r0; then SCF of r0..r3 (random, half of each species up, chosen to minimise species-resolved NN spin correlations),
  G-by-site (all NN AFM) and FM at the relaxed geometry. ECUT 90 Ry, k-spacing 0.28, gaussian 0.005, U Mn/Fe 4 (Cu 6),
  ortho-atomic, PseudoDojo SR NC, nosym (TR-only k set). Starting Hubbard occupations: Fe d5 (Hubbard_occ=5; QE default is
  d6 = Fe2+-like start, which biased earlier runs toward Mn3+/Fe2+), Cu d9.
- PM energy E_PM(sigma) = E0(sigma) of the Heisenberg mapping E = E0 - sum_k J_k C_k; estimated (i) as the plain mean of
  r0..r3, (ii) by a joint fit over all arrangements with shared species/direction-resolved J's (MnMn/MnFe/FeFe x ip/ap/Y).
- Arrangements (gen_tasks.py): 18 in 36-atom cells (A 1x1x2, B 2x1x1, C sqrt2-rotated) + 3 in 72-atom cells
  (SQS_D 2x2x1, SQS_F rot x 2c, AS1_D rock-salt + 1 NN antisite pair).
- YBaCuFeO5 benchmark parent: mp-505589 (P4mm, layered Cu/Fe = experimental bipyramid order) -> sqrt2 x sqrt2 x 2 cell,
  vc-relax at 110 Ry (lcm_hull110 settings, U Cu 6) in G-by-site: job lcm/ord2/vc_YBCFO_lay.

## 01:20 protocol fix: fixed total moment
- First submission (jobs lcm/ord2/YBMFO_*, cancelled after ~20 min): with free total moment, most non-rock-salt
  arrangements (LAY, STRd, MIX, RSapY, ...) sloshed for 20-36 SCF iterations with M drifting to 0.6-1.9 muB in
  balanced (M = 0 for d5/d5) configs = partially charge-transferred (Mn2+ -> Fe3+ spin-flip transfer), metallic solutions;
  rock-salt / columnar converged normally. Also confirmed QE's default starting Hubbard occupation for Fe is d6 (Fe2+-like).
- Resubmitted as lcm/ord2/v2/YBMFO_* with tot_magnetization = sum_i s_i m_i (m = 5 for Mn2+/Fe3+, 1 for Cu2+),
  Hubbard_occ Fe = 5 (d5 start), mixing_beta 0.2 (ndim 10). This is the "keep Mn2+/Fe3+" protocol required by the task;
  valence is checked per site from the ortho-atomic occupations; M-conserving charge transfer is still allowed and reported.
  A separate unconstrained probe (FIXM=0, d5 and QE-default starts) quantifies how much CT states could lower the
  non-rock-salt energies.

## 01:35 YBaCuFeO5 benchmark launched
- vc-relax (110 Ry, G-by-site, d5/d9 starts) of the sqrt2 x sqrt2 x 2 layered cell (from mp-505589): a_p = 3.893 A,
  c = 7.951 A (exp 3.87 / 7.66; MP PBE+U 3.90 / 7.95), P = 0.1 kbar, Cu m(ortho) 0.50 (n_d 9.39), Fe m 4.03 (n_d 5.90):
  Cu2+/Fe3+. -> `YBCFO_lay_h110.json`.
- 15 arrangements (same generator/protocol; BPR_F = 72-atom cell with every bipyramid Cu-O-Fe but random orientation,
  i.e. the experimentally proposed state) submitted as lcm/ord2/v2/YBCFO_*. Tasks: infra/jobs_src/lcm_ord2/tasks_ybcfo.json.

## 01:55 v3 relaunch (cost)
- The campaign has ~210 magdisc containers running; per-iteration cost in v2 was ~145 s for a 36-atom cell
  (2x the uncontended value), giving ~10 h per arrangement. v2 cancelled; relaunched as lcm/ord2/v3/* with
  the ions-only relax at 70 Ry and a k mesh halved in each direction (geometry only, energies are protocol SCFs at
  90 Ry / k-spacing 0.28), SCF conv_thr 1e-7 Ry/atom (<0.05 meV/cell). Fixed-moment / d5-start protocol unchanged.

## 02:15 v3 in-place restart (SCF robustness)
- Several non-rock-salt cells stagnated at 'estimated scf accuracy' ~1e-3 Ry for 60+ iterations (oscillation between
  nearly degenerate M-conserving charge-transfer patterns: e.g. COL_A r0 shows Mn m 4.04 / Fe m 3.79 pairs = partial
  Mn3+/Fe2+ even with the fixed total moment). Settings changed (jobs cancelled and resubmitted in place; the relax
  resumes from the last geometry): conv_thr 1e-6 Ry/atom (0.5 meV per 36-atom cell); relax electron_maxstep 80 with
  scf_must_converge = false, upscale 10; protocol SCFs electron_maxstep 150; an unconverged SCF keeps its last
  energy + accuracy (analysis accepts it only if accuracy < 2e-3 Ry/cell, flagged).

## 03:00 relax smearing / two-stage SCF (convergence of charge-transfer-prone cells)
- Non-rock-salt cells (LAY, LAYY, RSYap, SQS, ...) needed 80+ SCF iterations per relax step with 2-7x more Davidson work
  per iteration (sloshing between near-degenerate partial Mn3+/Fe2+ patterns, e.g. LAY_A r0: Fe m 3.2 / n_d 6.8).
- All jobs still relaxing were cancelled and resubmitted in place (relax resumes from the last geometry) with:
  relax smearing Marzari-Vanderbilt 0.02 Ry, mixing 0.1 (geometry only); protocol SCFs in two stages:
  (1) m-v 0.02 Ry, mixing 0.1, then (2) gaussian 0.005 Ry restarted from the stage-1 density (startingpot='file').
  The stage-2 (protocol) energy is used; if it fails, the converged stage-1 energy is used and flagged.
  YBMFO_RS_A, RS_C, LAYap_A had finished relaxing and continue with single-stage gaussian 0.005 SCFs (insulating).

## 03:20 orphan processes found -> clean v4 restart
- Diagnosis of the slow / erratic runs: `mrun.py cancel` (FunctionCall.cancel without terminate) interrupts the Modal
  function but leaves the `python ord2.py` + mpirun/pw.x process tree running; Modal then REUSES the container for new
  inputs. Scanning all magdisc containers (`modal container exec ... /proc`) found 57 ord2.py process trees, up to 5 x
  32-rank pw.x jobs per 32-core container (v1, v2, two v3 generations and the live one), all writing into the same
  job dirs (e.g. a crashed relax saved a premature 'relaxed' geometry into v3/YBMFO_STRp_C).
- All ord2 process trees were SIGKILLed via `modal container exec` (verified: 0 left). v1/v2/v3 data are discarded
  except as starting geometries: harvest_geom.py seeds each task with the best positions reached in v3 (clean
  converged relax or lowest-force last step) -> tasks_{ybmfo,ybcfo}_v4.json.
- ord2.py now writes owner.json (token) and polls it while pw.x runs (kills pw.x and exits if another submission
  claims the dir). Never use `mrun.py cancel` for these jobs; stop them by killing the process tree in the container.
- v4: 36 jobs lcm/ord2/v4/{YBMFO,YBCFO}_* (same physics settings: fixed moment, d5/d9 starts, m-v 0.02 relax at
  70 Ry/half k, two-stage SCF, energies at 90 Ry / k 0.28 / gaussian 0.005).

## 2026-10-02 10:45-12:00 analysis pass (session 66190d; PROVISIONAL, 13/21 YBMFO and 9/15 YBCFO arrangements complete)
- Status at 10:47: v4 jobs started 03:21; done = YBMFO RS_A, RS_C, COL, RSapY, RSYap, LAYap, CHKmix, MIX_A/B/C, STRd, STRdc, STRp,
  STRpc, HALF; YBCFO RS, LAY, RSapY, LAYap, MIXap, STRd, STRp, MIX_C. Running: YBMFO LAY (3/6 runs), LAYY (1), MIXap (2), SQS_D
  (relaxed, 0 runs), SQS_F (relaxing), AS1_D (1); YBCFO COL (2), LAYY (3), RSYap (0), MIX_A (1), BPR_F (1), SQS_D + CHKmix
  (relaxing); probes CTP2, KCHK2.
- Launched helper SCF jobs (helpers.py) for the relaxed 72-atom cells: lcm/ord2/v4/{YBMFO_SQS_D,YBMFO_AS1_D,YBCFO_BPR_F}_{h1,h2}
  (h1 = FM,G; h2 = r3,r2). fetch.py rewritten to read job dirs from infra/job_registry.jsonl and results/<ID>.json directly
  (no VolumeListFiles calls) and to merge helper runs.
- Data quality: all 96 YBMFO / 55 YBCFO accepted runs are converged protocol (stage-2 gaussian 0.005) SCFs; no unconverged or
  m-v-only energies were needed. 23/96 YBMFO runs are charge-transferred (any Mn3+/Fe2+ site): all runs of LAY_A, STRdc_B, STRpc_C,
  MIX_A; LAYap_A r2/r3/G; LAYY_A r0; AS1_D r0 (the NN antisite pair converts to Mn3+(Fe site)/Fe2+(Mn site)). Rock-salt, RSapY, STRp,
  CHKmix, HALF, STRd, RSYap, COL, MIX_B/C, MIXap are clean Mn2+/Fe3+.
  Probes: CTP1 (COL_A r0/G without fixed moment, d5 and QE starts) = constrained energy within 0.2 meV/cell; CTP2 MIX_A r0 (d5, free M)
  -33 meV/cell (-8 meV/f.u.) -> the fixed-moment constraint is not a significant bias. KCHK1: rock-salt in cells A/B/C differs by
  <= 0.6 meV/f.u. (protocol k meshes are consistent across cell shapes); RS_A vs RS_C (same order, different r0 relax geometry)
  differ by 8 meV/f.u. -> geometry noise sigma_geo = 6 meV/f.u.
- Why analyze.py failed: its joint Heisenberg fit mixed FM (half-metallic, +81 meV/f.u. above the PM line on average) and CT runs
  (rms 315 meV/cell). New analyze2.py: J fit on clean non-FM runs only (rms 26.6 meV/cell YBMFO, 11.2 YBCFO; species-resolved J;
  individual J's are poorly conditioned because random configs are near-uncorrelated by design, but E_PM is robust: fixing the
  Mn-Fe J's to the jfit110 values moves E_PM by <= 10 meV/f.u.; a species-independent J model (rms 75) moves near-RS E_PM by
  10-40 meV/f.u.). Per run E0 = E + x.J; per arrangement the low branch (runs within 60 meV/f.u. of the minimum) is used, which
  drops metastable CT solutions that do not match the r0-relaxed geometry (e.g. STRdc_B r0 63 vs r1-r3/G 440-560 meV/f.u.).
- Relax-level cross-check (relaxlevel.py): the r0 relax energy (70 Ry, half k, m-v 0.02), spin-corrected, reproduces the protocol
  r0 E0 to +-2-4 meV/f.u. for clean A/B/D cells (C cell offset -56 +- 3, k-inequivalent half mesh). It exposes metastable protocol
  states: LAY_A protocol 373 vs relax-level 258, LAYY_A 487 vs 95 meV/f.u. (YBCFO MIX_A 169 vs 90). 'best' = protocol low branch
  unless the relax-level state is > 30 meV/f.u. lower; SQS_D (no SCF yet) enters at relax level.
- PM ordering energies, YBaMnFeO5, meV/f.u. vs rock-salt RS_A (best; [] = source/flags):
  RS_C -8 | RSapY 18 | STRp_C 25 | AS1_D 27 (= 0.22 eV per NN antisite pair; CT; r0 only) | CHKmix 61 | STRdc_B 63 [CT, r0] |
  HALF 74 | MIX_A 85 [CT] | STRd 86 | LAYY 95 [relax-level, CT] | RSYap 99 | STRpc 107 [CT, r0] | SQS_D 116 +- 15 [relax-level] |
  COL 122 | MIX_C 166 | MIXap 180 [2 runs] | LAYap 200 | MIX_B 216 | LAY 258 [relax-level, CT].  (sigma 6-7 for 5-run clean entries.)
  Context: the earlier G-state antisite calc (lcm/s1, 75 Ry) gives 0.84 eV per NN pair in site-following G and 0.09 eV in FM;
  the PM value lies in between. The random-state energy (116) is half of the earlier unrelaxed pair-model estimate (230).
- YBaCuFeO5 (vs RS_A): LAY -107 (ground state = polar P4mm order) | LAYap -81 | MIXap -44 | STRd -42 | BPR_F -22 (r0 only) |
  RSapY 1 | STRp 21 | MIX_C 63 | MIX_A 90 [relax-level] | SQS_D 126 [relax running] | CHKmix 206 [relax running] | COL 419 |
  RSYap 421 [relax-level] | LAYY 501. Cu-Fe bipyramids (ap hetero) strongly favoured (ap-homo structures +420-500).
- Cluster expansion (analyze2.py): subset CEs are ill-conditioned (Y/Yd/c2 collinear in the A-cell data; CT states are not
  pair-additive): weighted LOO-CV 30-50 meV/f.u. Main model = Bayesian ridge CE over all 8 pair classes (shell priors 100/15/8 meV,
  scale 0.25 by CV), Boltzmann-weighted (kT_w = 100 meV) on the 'best' data: 
  YBMFO V (meV/bond, +- bootstrap): ip1 14.1+-6.6, ap 23.1+-17.4, Y 24.3+-16.3, ip2 4.3+-4.7, apd -7.2+-3.6, Yd 2.4+-3.3, ip3 -3.0+-2.5,
  c2 1.8+-1.3; wrms 22, CVw 39 meV/f.u.; predicts E(random) - E(RS) = 114 (SQS_D: 116).
  YBCFO V: ip1 -11.1+-6.4, ap 81.8+-15.9, Y 4.8+-11.9, ip2 -7.1+-3.2, apd 0.8, Yd 1.3, ip3 2.0, c2 -0.4; wrms 17, CVw 33.
- Canonical MC (torder.py, L = 8/12/16 + 30 bootstrap CEs; sro.py cooling runs, L = 16):
  YBMFO main: rock-salt ground state, C peak 915/965/915 K (L 8/12/16), eta_RS jumps between 1073 and 973 K; bootstrap 970 +- 250
  (16-84 %: 800-1140 K). CE variants with RS ground state: 836 (low-branch, frustrated ip1/ip2), 1057 (all runs incl. metastable
  CT), 1149 (subset CE), 1210 K (NN-only). The ionic-only (clean) CE predicts a NON-RS ground state q=(0,1,1/2) (in-plane stripes
  + Y-homo stacking) at 1500 K - an extrapolation outside the training set, flagged for a DFT check.
  -> T_order(YBaMnFeO5) = 950 (+250/-150) K, i.e. ~680 C (530-930 C).
  SRO at 900-1300 C (1173-1573 K): no LRO; hetero NN-pair fractions ip1/ap/Y = 0.56-0.59/0.60-0.66/0.58-0.61 (random 0.50).
  Equilibrium antisite fraction below T_order: 13 % (873 K), 6 % (773 K), 3 % (673 K).
  YBCFO main: LAY (polar) order, C peak 1118/1149/1118 K, bootstrap 1080 +- 280 (900-1200); alternative subset CE 1358 K.
  eta_LAY(eq.) ~0 at 1423 K, 0.06 at 1273 K, 0.19 at 1173 K, 0.71 at 1073 K; Cu-Fe bipyramid fraction 0.81 at 1423 K.
- Benchmark vs experiment (Morin et al., Nat. Commun. 7, 13758 (2016): sintered 1150 C/50 h, cooled 5-500 K/h or LN2-quenched;
  P4mm Cu/Fe occupation difference = eta_LAY = 0.10 quenched -> 0.16 at 5 K/h; LOG expectation: orientationally random Cu-Fe
  bipyramids = BPR_F). Protocol T_order(YBCFO) = 1120 K (900-1360) = 0.79 x the sintering T: consistent with the weak measured
  LRO only if B-site exchange freezes at >~1150 K on cooling. No gross overestimate of ordering tendency; a 10-25 % overestimate
  cannot be excluded. Practical lesson: ordering that requires equilibration below ~1150 K is not reached in these 112 ferrites.
- VERDICT (provisional): rock-salt Mn/Fe order in YBaMnFeO5 is NOT thermodynamically achievable by conventional annealing.
  T_order (950 K; 800-1210 K) lies below the 1173-1573 K synthesis/anneal window and below the benchmark's own T_order (1120 K),
  where 5 K/h cooling gives only eta ~ 0.16. >90 % order (x_AS < 5 %) needs equilibration at <~800 K; marginal at best at the top
  of the error bar. This explains the B-disordered GdBaMnFeO5/NdBaMnFeO5 samples thermodynamically (not only kinetically), and
  it supersedes the earlier lane estimates (Ising ~2000 K; PM-proxy ~1300 K), which ignored Mn2+/Fe3+ -> Mn3+/Fe2+ charge transfer
  in non-rock-salt environments and the near-degenerate in-plane stripe / Y-homo stacking variants.
- Still provisional / pending: SQS_D protocol SCFs (main + h1/h2), SQS_F (relaxing), AS1_D r1-r3/G (main + helpers), LAY/LAYY/MIXap
  remaining runs, YBCFO BPR_F r1-r3/G (helpers), YBCFO SQS_D + CHKmix (relaxing), RSYap/COL/LAYY/MIX_A SCFs, CTP2, KCHK2.
  Most consequential: the two YBMFO SQS energies and the AS1_D average (they anchor the disorder energy), and YBCFO BPR_F.
  Recommended extra DFT check (not launched): STRp with Y-homo stacking (q=(0,1,1/2), F-type cell), the ground state of the
  ionic-only CE; and a few more Y-homo/ap-hetero stackings to break the Y/Yd/c2 collinearity.
- Files: analyze2.py -> ce2_<TAG>.json; relaxlevel.py -> relaxlevel_<TAG>.json; torder.py -> torder_<TAG>_<KEY>[_CSEL].json
  (run_torder.sh launches the 8 variants; logs in mc_logs/); sro.py -> sro_<TAG>_<KEY>.json; make_summary.py -> ord2_summary.json
  (its docstring has the full rerun recipe). ce_YBMFO.json (old analyze.py) is superseded.
