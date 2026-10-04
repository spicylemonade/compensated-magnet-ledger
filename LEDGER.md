# Computational ledger

Every number in the blog post and the README is listed below with:

- the calculation it came from,
- the exact input and raw output files in this repository,
- a script that recomputes it.

Run one command to check all of them:

```bash
pip install numpy
python tools/verify.py
```

It re-reads the raw Quantum ESPRESSO outputs and recomputes each value: net magnetisation, band gap, spin windows and energy differences. It then compares each against the value the agents recorded. It needs only numpy and takes under a second. On 2026-10-04 it reported **55 pass, 0 fail**; three entries are experimental or literature values and cannot be computed.

## Tiers

| Tier | Meaning | How to check |
|---|---|---|
| **A** | Recomputed from raw `pw.x` output files in `materials/*/runs/` | `python tools/verify.py` |
| **L2** | Recomputed end to end by an included script from the raw energies (exchange fit and Monte Carlo; cluster-expansion Monte Carlo) | `materials/YBaMnFeO5/exchange_fit_and_Neel_T/rerun_fit_and_mc.py`, `materials/YBaMnFeO5/cation_order_cluster_expansion/rerun_torder.sh` |
| **B** | Read from a recorded analysis file in this repository; the script that produced it is included, but not every intermediate input is | open the cited JSON |
| **E** | Experimental or literature value; cannot be computed | see the citation |
| **R** | Re-run from scratch in a fresh cloud container from the inputs in this repository, with the pseudopotentials downloaded from the public source | `reproduce/`, then `python tools/verify.py --repro` |

## Definitions used throughout

- **Compensated (LCM) state:** the collinear magnetic order in which the two different magnetic sublattices point opposite ways. In this state the total spin moment of the cell is zero.
- **Spin window:** scanning the whole Brillouin zone, the energy range at a band edge that holds only one spin.
  - Hole window: `win_VB = VBM(edge spin) − VBM(other spin)`.
  - Electron window: `win_CB = CBM(other spin) − CBM(edge spin)`.
  - A state is occupied if its occupation is above 0.5.
  - These are the same definitions as the original analysis code; see `tools/qe_parse.py`.
- **Unipolar:** both band edges are in the same spin channel. **Bipolar:** they are in opposite channels.
- **Energies per magnetic ion:** `(E_state − E_ground) / (number of magnetic ions in the cell)`.
- **Moments:** QE's sphere-integrated "magn=" values. They are smaller than the nominal 3 μB (V²⁺, Cr³⁺) or 5 μB (Mn²⁺, Fe³⁺) because the spheres are small.

## Software and settings (common to all runs)

| Setting | Value |
|---|---|
| Code | Quantum ESPRESSO `pw.x` 7.5 (conda-forge build), MPI on cloud CPUs (Modal) |
| Pseudopotentials | PseudoDojo v0.4, PBE, scalar-relativistic, "standard", norm-conserving. md5 sums are in `reproduce/pseudo_md5.txt`; the public download matches them exactly (checked 2026-10-04) |
| DFT+U | Dudarev, ortho-atomic projectors. KV[Cr(CN)₆]: U(V 3d) = U(Cr 3d) = 3 eV, with a U grid {2, 3, 4}. YBaMnFeO₅: U(Mn 3d) = U(Fe 3d) = 4 eV, with scans over 0–6 eV |
| Hybrid functional | HSE06 without U, `ecutfock = 2 × ecutwfc`, Gygi–Baldereschi divergence treatment with `x_gamma_extrapolation` |
| Smearing | Gaussian, 0.005 Ry |
| Spin | Collinear (`nspin = 2`), no spin–orbit coupling unless stated |

Every input file states its own cutoff, k-mesh and q-mesh. `ledger/runs.csv` indexes all 128 input files (120 `pw.x` calculations and 8 `projwfc.x` post-processing inputs) with these settings and the final energy and magnetisation of each run.

## Run folders

The full index of inputs is in `ledger/runs.csv`. Folder descriptions are in `ledger/run_folder_descriptions.json`.

<details><summary>KV[Cr(CN)₆] run folders</summary>

| Folder | What it is |
|---|---|
| `A_pbeu_relax_scf_nscf_Ugrid` | vc-relax of the 15-atom primitive cell (F-43m). Then scf (5×5×5) and nscf (8×8×8) of the compensated state, plus the FM reference, at the base U and on the U grid |
| `B_pbeu_FM_fixed_moment` | FM with the total moment fixed at 6 μB (exchange-energy reference) |
| `C_hse06_ideal_LCM` | HSE06, ideal anhydrous crystal |
| `D_hse06_ideal_FM_fixed_moment` | HSE06 FM reference |
| `E_hse06_ideal_LCM_nq1` | HSE06, ideal crystal, q-mesh 1 (matched to the hydrate run) |
| `F_hse06_dihydrate_LCM` | HSE06, KV[Cr(CN)₆]·2H₂O. Stopped after the exact-exchange error plateaued; the last complete eigenvalue block is used |
| `G_pbeu_dihydrate` | PBE+U, same hydrate |
| `H_hse06_vacancy_with_water` | HSE06, 65-atom cell with a water-filled [Cr(CN)₆] vacancy |
| `I_pbeu_vacancy_with_water` | PBE+U, same cell |
| `J`–`N` | 60-atom defect cells: pristine, V/Cr antisite pair, flipped cyanide, K vacancy, dry [Cr(CN)₆] vacancy |
| `O_pbeu_reference_K_H2O` | Reference energies for bcc K and H₂O |
| `P_pbeu_bands_pdos` | Band structure and orbital-projected DOS |
| `Q_hse06_control_CrCr_same_element` | Control: Cr[Cr(CN)₆], with the same metal on both sites |
</details>

<details><summary>YBaMnFeO₅ run folders</summary>

| Folder | What it is |
|---|---|
| `A_pbeu_vcrelax_110Ry` | vc-relax of the 18-atom rock-salt-ordered P4/n cell (G-type) |
| `B_hse06_G_type_LCM`, `C_hse06_Yflip_competitor` | HSE06: the ground state and its closest competitor |
| `D_pbeu_72atom_G_FM_and_spin_flips` | 72-atom cell: G, FM, and spin-flip snapshots |
| `E`–`G` | Mn/Fe antisite study: reference, nearest-neighbour pair, far pair |
| `H_pbeu_stacking_vs_U` | G vs Y-layer flip vs Ba-layer flip at six U combinations |
| `I_pbeu_exchange_fit_configs` | 24 spin configurations for the exchange fit |
</details>

## Claims

The tables below are generated by `python tools/verify.py --write-ledger`. **recorded** is the value the agents wrote in their reports; **recomputed** is what the script gets from the files today.

<!-- CLAIMS-TABLE:BEGIN -->

#### KV[Cr(CN)6]

| id | tier | claim | recorded | recomputed | status | raw files |
|---|---|---|---|---|---|---|
| K01 | A | Relaxed cubic lattice constant of the ideal anhydrous crystal (exp. 10.55 A, hydrated powder) (A) | 10.680 | 10.680 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_vcrelax.out` |
| K02 | A | Net spin moment of the ideal crystal (zero = fully compensated) (Bohr mag/cell) | 0.000 | -0.000 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.out` |
| K03 | A | V moment (N-bound site), integrated in a small sphere: antiparallel to Cr (Bohr mag) | -2.020 | -2.016 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.out` |
| K04 | A | Cr moment (C-bound site), integrated in a small sphere (Bohr mag) | 2.320 | 2.318 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.out` |
| K05 | A | Band gap (eV) | 1.950 | 1.953 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.out` |
| K06 | A | Both band edges in the same spin channel (unipolar) (bool) | yes | yes | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.out` |
| K07 | A | Hole spin window (top of valence band that is 100 % one spin) (eV) | 2.020 | 2.023 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.out` |
| K08 | A | Electron spin window (bottom of conduction band that is 100 % one spin) (eV) | 1.150 | 1.146 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.out` |
| K09 | A | Robust to the Hubbard U choice: smallest hole window over the U grid (U_V, U_Cr in {2,3,4} eV) (eV) | 1.390 | 1.387 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_*_LCM_nscf.out` |
| K10 | A | Smallest electron window over the U grid (eV) | 0.870 | 0.873 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_*_LCM_nscf.out` |
| K11 | A | Unipolar at all five U points (bool) | yes | yes | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_*_LCM_nscf.out` |
| K12 | A | Compensated state is far below the ferromagnet: E(FM) - E(compensated) per magnetic ion (meV) | 152.2 | 152.2 | PASS | `KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.out`<br>`KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_FM.out` |
| K13 | A | HSE06: net spin moment (Bohr mag/cell) | 0.000 | -0.000 | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out` |
| K14 | A | HSE06: band gap (eV) | 2.090 | 2.091 | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out` |
| K15 | A | HSE06: unipolar band edges (bool) | yes | yes | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out` |
| K16 | A | HSE06: hole spin window (eV) | 2.640 | 2.636 | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out` |
| K17 | A | HSE06: electron spin window (eV) | 1.570 | 1.567 | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out` |
| K18 | A | HSE06: E(FM, moment fixed at 6) - E(compensated) per magnetic ion (meV) | 181.1 | 181.1 | PASS | `KV_Cr_CN6/runs/C_hse06_ideal_LCM/hse_lcm.out`<br>`KV_Cr_CN6/runs/D_hse06_ideal_FM_fixed_moment/KVCr_hse_FMfix.out` |
| K19 | A | HSE06, dihydrate KV[Cr(CN)6]*2H2O (one water arrangement; exact-exchange loop stopped before full convergence): hole window (eV) | 2.430 | 2.432 | PASS | `KV_Cr_CN6/runs/F_hse06_dihydrate_LCM/hse_lcm.out` |
| K20 | A | HSE06, dihydrate: electron window (still unipolar) (eV) | 1.420 | 1.415 | PASS | `KV_Cr_CN6/runs/F_hse06_dihydrate_LCM/hse_lcm.out` |
| K21 | A | PBE+U, same dihydrate: hole window collapses (methods disagree about water) (eV) | 0.930 | 0.929 | PASS | `KV_Cr_CN6/runs/G_pbeu_dihydrate/KVCr_2H2O_LCM_nscf.out` |
| K22 | A | HSE06, water-filled [Cr(CN)6] vacancy (65-atom cell): net moment = -3 (composition, not symmetry, sets the moment) (Bohr mag/cell) | -3.000 | -3.000 | PASS | `KV_Cr_CN6/runs/H_hse06_vacancy_with_water/hse_lcm.out` |
| K23 | A | HSE06, water-filled vacancy: electron window (in-cell) (eV) | 0.740 | 0.737 | PASS | `KV_Cr_CN6/runs/H_hse06_vacancy_with_water/hse_lcm.out` |
| K24 | A | HSE06, water-filled vacancy: hole window (in-cell) (eV) | 2.800 | 2.803 | PASS | `KV_Cr_CN6/runs/H_hse06_vacancy_with_water/hse_lcm.out` |
| K25 | A | PBE+U, V/Cr antisite pair (60-atom cell): moment stays 0 ... (Bohr mag/cell) | 0.000 | 0.000 | PASS | `KV_Cr_CN6/runs/K_pbeu_defect_antisite_pair/KVCr_AS_lcm.out` |
| K26 | A | ... but the hole window shrinks to (in-cell; 0.13 eV after alignment to the host) (eV) | 0.230 | 0.229 | PASS | `KV_Cr_CN6/runs/K_pbeu_defect_antisite_pair/KVCr_AS_lcm_nscf.out` |
| K27 | A | Control, same-element Cr[Cr(CN)6] (HSE06): edges in opposite spins (bipolar) ... (bool) | no | no | PASS | `KV_Cr_CN6/runs/Q_hse06_control_CrCr_same_element/hse_lcm.out` |
| K28 | A | ... with a tiny hole window (why two different metals matter) (eV) | 0.100 | 0.105 | PASS | `KV_Cr_CN6/runs/Q_hse06_control_CrCr_same_element/hse_lcm.out` |
| K29 | B | Spin space group of the ideal crystal (findspingroup): compensated ferrimagnet, Zeeman-type (s-wave) splitting (label) | 216.216.1.1.L | 216.216.1.1.L | PASS | `KV_Cr_CN6/analysis_data/pba_spingroup_check.txt` |
| K30 | E | Measured magnetic ordering temperature of KV[Cr(CN)6]*2H2O (Holmes & Girolami, JACS 121, 5593 (1999)) (K) | 376 | - | n/a | `literature` |
| K31 | E | Measured saturation moment at 5 K (expected 0 for perfect stoichiometry) (Bohr mag/f.u.) | 0.125 | - | n/a | `literature` |

#### YBaMnFeO5

| id | tier | claim | recorded | recomputed | status | raw files |
|---|---|---|---|---|---|---|
| Y01 | A | Relaxed cell, rock-salt-ordered P4/n: a (A) | 5.665 | 5.665 | PASS | `YBaMnFeO5/runs/A_pbeu_vcrelax_110Ry/vc110.out` |
| Y02 | A | Relaxed cell: c (A) | 7.689 | 7.689 | PASS | `YBaMnFeO5/runs/A_pbeu_vcrelax_110Ry/vc110.out` |
| Y03 | A | Net spin moment of the G-type compensated state (Bohr mag/cell) | 0.000 | 0.000 | PASS | `YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U44_G/single.out` |
| Y04 | A | Mn moment (sphere-integrated); Fe moment is about -3.7, antiparallel (Bohr mag) | 4.000 | 4.004 | PASS | `YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U44_G/single.out` |
| Y05 | A | Closest competing magnetic order (spin flip across the Y layer) above the ground state, U = 4 eV (meV/mag. ion) | 8.900 | 8.881 | PASS | `YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U44_G/single.out`<br>`YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U44_Yflip/single.out` |
| Y06 | A | Same competitor at U = 0 eV (meV/mag. ion) | 16.400 | 16.369 | PASS | `YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U00_*` |
| Y07 | A | Same competitor at U = 6 eV (smallest margin in the U scan) (meV/mag. ion) | 7.000 | 6.978 | PASS | `YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U66_*` |
| Y08 | A | Ferromagnet above the ground state, U = 4 eV (72-atom cell, 75 Ry) (meV/mag. ion) | 199.8 | 199.8 | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G.out`<br>`YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_FM.out` |
| Y09 | A | HSE06: competitor (Y-layer flip) above the ground state (meV/mag. ion) | 8.500 | 8.528 | PASS | `YBaMnFeO5/runs/B_hse06_G_type_LCM/hse_lcm.out`<br>`YBaMnFeO5/runs/C_hse06_Yflip_competitor/hse_lcm.out` |
| Y10 | A | HSE06: band gap (eV) | 2.350 | 2.351 | PASS | `YBaMnFeO5/runs/B_hse06_G_type_LCM/hse_lcm.out` |
| Y11 | A | HSE06: unipolar band edges (bool) | yes | yes | PASS | `YBaMnFeO5/runs/B_hse06_G_type_LCM/hse_lcm.out` |
| Y12 | A | HSE06: hole spin window (eV) | 1.000 | 0.995 | PASS | `YBaMnFeO5/runs/B_hse06_G_type_LCM/hse_lcm.out` |
| Y13 | A | HSE06: electron spin window (eV) | 1.400 | 1.398 | PASS | `YBaMnFeO5/runs/B_hse06_G_type_LCM/hse_lcm.out` |
| Y14 | A | PBE+U band gap (72-atom cell) (eV) | 1.400 | 1.401 | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out` |
| Y15 | A | PBE+U windows: hole (eV) | 0.450 | 0.445 | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out` |
| Y16 | A | PBE+U windows: electron (eV) | 1.430 | 1.432 | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out` |
| Y17 | A | Worst-case spin disorder (2 of 16 spins flipped): electron edge keeps its spin, window (eV) | 0.155 | 0.155 | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out` |
| Y18 | A | ... while the hole edge switches to the other spin (holes need good magnetic order) (bool) | no | no | PASS | `YBaMnFeO5/runs/D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out` |
| Y19 | A | One neighbouring Mn/Fe swap (12.5 % of sites, 72-atom cell): gap of the other spin channel collapses (dense grid) (eV) | 0.010 | 0.011 | PASS | `YBaMnFeO5/runs/F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite_nscf.out` |
| Y20 | A | ... same, on the coarser self-consistent grid (eV) | 0.300 | 0.300 | PASS | `YBaMnFeO5/runs/F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite.out` |
| Y21 | L2 | Exchange fit from 25 raw energies: in-plane Mn-O-Fe coupling J_ip (meV) | -32.600 | -32.609 | PASS | `YBaMnFeO5/exchange_fit_and_Neel_T/rerun_fit_and_mc.py`<br>`YBaMnFeO5/runs/I_pbeu_exchange_fit_configs/`<br>`YBaMnFeO5/runs/H_pbeu_stacking_vs_U/U44_*` |
| Y22 | L2 | Neel temperature, classical Monte Carlo on the fitted model (Binder crossings, L = 6/8/10) (K) | 417 | 414.3 | PASS | `YBaMnFeO5/exchange_fit_and_Neel_T/rerun_fit_and_mc.py`<br>`YBaMnFeO5/exchange_fit_and_Neel_T/jfit110/mc_clean_Mle10.5_10.json` |
| Y23 | B | Calibrated Neel temperature: Y22 x (measured/simulated ratio for YFeO3 with the same protocol, 1.17-1.19) (K) | 490 | 489.6 | PASS | `YBaMnFeO5/analysis_data/calibration_summary.json` |
| Y24 | B | Energy above the convex hull at 0 K (26 competing phases) (meV/atom) | 2.600 | 2.600 | PASS | `YBaMnFeO5/analysis_data/hull_QE110.json` |
| Y25 | L2 | DECISIVE: Mn/Fe order-disorder temperature from the paramagnetic cluster expansion (C peak, L = 16) (K) | 915 | 915.0 | PASS | `YBaMnFeO5/cation_order_cluster_expansion/ce2_YBMFO.json`<br>`YBaMnFeO5/cation_order_cluster_expansion/torder.py`<br>`YBaMnFeO5/cation_order_cluster_expansion/rerun_torder.sh` |
| Y26 | B | Headline order-disorder temperature with uncertainty (bootstrap + model variants) (K) | 950 (+250/-150) | 950 (+250/-150) | PASS | `YBaMnFeO5/cation_order_cluster_expansion/ord2_summary.json` |
| Y27 | E | Typical synthesis/anneal window for these layered perovskites (900-1300 C); B-site exchange inferred to freeze below ~1150 K (YBaCuFeO5 benchmark: Morin et al., Nat. Commun. 7, 13758 (2016)) (K) | 1173-1573 | - | n/a | `literature` |

<!-- CLAIMS-TABLE:END -->

## Independent re-runs (tier R)

On 2026-10-04 a fresh cloud container re-ran four sets of calculations from the input files in this repository, using:

- a clean conda-forge QE 7.5 image with no campaign code;
- pseudopotentials downloaded from pseudo-dojo.org, all md5-identical to the originals.

The procedure is `reproduce/modal_repro.py`, and the results are in `reproduce/results/` and `reproduce/RESULTS.md`. Check them with:

```bash
python tools/verify.py --repro
```

## Known caveats and method issues

Most of these come from the agents' own adversarial reviews and cross-checks. They are listed so that nobody has to rediscover them.

**Both materials**
1. **Ideal crystals.** Every DFT number is for a perfectly ordered crystal at 0 K. Finite-temperature spin disorder was tested only for YBaMnFeO₅ (spin-flip snapshots, Y17–Y18).
2. **DFT+U is a choice.** U values are literature-typical, not computed from first principles. Attempts to compute U by linear response (`hp.x`) never produced a result. That is why the U scans (K09–K11, Y05–Y07) and the independent HSE06 checks are part of the ledger.
3. **HSE06 windows come from the self-consistent k-grid, not a dense grid.** For KV[Cr(CN)₆] the 4×4×4 grid includes Γ, X and L. Band extrema off the grid could make windows slightly smaller.
4. **HSE06 energy differences between states of different magnetic symmetry can carry exact-exchange k-set offsets.** Gaps and windows are not affected. The KV[Cr(CN)₆] FM–LCM pair (K18) has the same crystal symmetry.
5. **Spin–orbit numbers were not used.** Earlier spin–orbit runs began from antiparallel non-collinear starting moments, which can trigger a known QE sign issue, so they are marked unverified.

**KV[Cr(CN)₆]**

6. **The real material is not the ideal crystal.** The only sample ever reported (1999) is a hydrated powder, never reproduced. Its saturation moment is 0.125 μB/f.u. instead of 0, and every missing [Cr(CN)₆] unit adds 3 μB.
7. **Water: the two methods disagree.**
   - PBE+U says the dihydrate's hole window drops from 2.02 to 0.93 eV (K21).
   - HSE06 says 2.64 → 2.43 eV (K19, matched comparison 2.68 → 2.43).
   - The agents judged HSE06 more reliable, because PBE+U misplaces the water levels. That is a judgement, not a measurement.
   - The HSE06 hydrate run was stopped before the exact-exchange loop fully converged. The valence-band maximum was still rising by about 16 meV per cycle, so the best estimate of the hole window is about 2.35 eV (range 2.25–2.45).
8. **Defect windows are "in-cell".** Aligning defect cells to the host bands (done by the original collectors, values in `runs/*/results/*.json`) gives smaller numbers in some cases. Examples: antisite hole window 0.23 eV in-cell vs 0.13 eV aligned; PBE+U water-filled-vacancy electron window 0.47 eV in-cell vs 0.07 eV aligned.
9. **Carriers will be heavy.** The relevant bands are 0.4–0.6 eV wide and a doped hole self-traps (polaron), so do not expect silicon-like transport. No conductivity, optical-gap or spin-resolved measurement exists for this compound.
10. **At room temperature the order is far from perfect.** The sublattice order is about 0.6 at 300 K (T_C = 376 K), so the ideal spin polarisation would be reduced.
11. **The physics is not new in general.** Compensated ferrimagnets are known to have spin-split bands. What is new is identifying this room-temperature compound as one and putting numbers on its band edges. The closest prior work is Schart et al., Inorg. Chem. 63, 22856 (2024), on Cr[Cr(CN)₆]; its windows are small and bipolar (K27–K28).

**YBaMnFeO₅**

12. **It has never been made.** The key result is negative: the rock-salt Mn/Fe order needed for the effect is predicted to disorder at about 950 K (Y25–Y26). That is below the temperatures at which cations move quickly during standard synthesis.
    - Every chemically similar compound whose B-site arrangement has been determined (GdBaMnFeO₅, NdBaMnFeO₅₊δ, YBaMnCoO₅) is B-site disordered. SmBaMnFeO₅₊δ has been made, but its B-site order is not reported.
    - A single nearest-neighbour Mn/Fe swap collapses the opposite-spin gap (Y19–Y20).
13. **The freeze-out temperature is an analogy.** The ~1150 K below which B-site exchange is taken to freeze comes from a different compound, YBaCuFeO₅ (Morin et al., Nat. Commun. 7, 13758 (2016)), computed with the same protocol.
    - That compound has Jahn–Teller physics that YBaMnFeO₅ lacks.
    - The highest cluster-expansion variant (1210 K) reaches into the synthesis window.
14. **The cluster expansion is provisional.** Weighted leave-one-out error is 39 meV/f.u.
    - A charge-transfer-free variant predicts a different, non-rock-salt ground state, which has not been checked by DFT. If it held, rock-salt order would not be the ground state either.
    - The raw outputs of its 96 DFT runs are not included because of their size. Their parsed energies, moments and charge-transfer flags are in `cation_order_cluster_expansion/results_YBMFO/` and `ce2_YBMFO.json`.
15. **The HSE06 YBaMnFeO₅ runs used a 75 Ry cutoff on the 110 Ry geometry.** An earlier HSE06 result (gap 2.32 eV, windows 1.10/1.25 eV) was on a superseded 75 Ry geometry; it is replaced by Y10–Y13.
16. **The Néel temperature is a model estimate.**
    - 417 K (raw classical Monte Carlo) and about 490 K (calibrated on YFeO₃) come from a fitted Heisenberg model.
    - At U = 6 eV the raw value drops to about 320–340 K.
    - PBE+U gets the ground state of the twin compound YBaMn₂O₅ wrong at U ≥ 2 eV. It keeps the right ground state for YBaMnFeO₅ at all U, and HSE06 agrees (Y09).
17. **Hull and phonons.**
    - The 0 K hull distance (+2.6 meV/atom) is against the competing phases that were computed; a few were still pending, and missing competitors can only raise it.
    - Single-cycle vc-relax energies should be re-relaxed once before precise reuse.
    - Phonons are real at every exactly computed q-point (Γ, Z, X, Y, M). Fourier-interpolated branches between them show imaginary values in the 1×1×2 supercell, which the agents attributed to interpolation; larger supercells were not run (`analysis_data/phonon110_summary.json`).

## What is deliberately not in this repository

- **The rest of the campaign.** This was one lane of a larger multi-agent search; the other lanes found nothing above their bars. Its other compounds, including several unpublished leads, are left out on purpose.
- **Very large raw files.** Wavefunctions, projections (`atomic_proj.xml`, up to 80 MB each) and the outputs of the 96 cluster-expansion and 26 hull runs. Their parsed results are included.
- **Literature PDFs and experimental crystal structures from databases with restrictive licences.**
