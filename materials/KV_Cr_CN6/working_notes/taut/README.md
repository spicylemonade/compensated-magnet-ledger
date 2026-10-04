# PBA taut: is the d3/d3 LCM the electronic ground state of KV[Cr(CN)6] and KV[Mo(CN)6]?

Track L (session 66190d), 2026-10-03. Folder `tracks/lcm/pba/taut/`. Gate 0 is in `../README.md`, the defect gate in `../defects/README.md`, the T_C model in `../tc/README.md`.

- Jobs: prefix `lcm/pba/taut/`. Call ids, commands, script md5s and fixed result paths are in `manifest.json`.
- Submit: `python tracks/lcm/pba/taut/submit_taut.py [--dry] [--only JOB]`. It never cancels anything.
- Collect: `python tracks/lcm/pba/taut/collect.py --fetch`. It reads fixed paths only (no listdir), caches them in `results/`, prints every start and writes `results/collect_snapshot.json`. `--peek JOB FILE` prints the tail of a running pw.x output.
- `summary.json` holds the numbers quoted below.

## 0. Status and bottom line (2026-10-03 10:50 PDT; snapshot in `summary.json`)

**Bottom line: no kill so far.** Every competitor that converged in PBE+U lies above the LCM.
- Neither compound has a metastable charge-transferred LCM-competitor at any U point tested so far. An orbitally polarised V(III)/M(II) start (V t2g↓², M t2g↓¹) relaxes back to the d3/d3 LCM within about 4 SCF iterations.
- The HSE06 kill test has **not run yet**. It is the last step of the four base-U chain jobs: vertical scan → relaxation of the lowest competitor → HSE06.

**PBE+U vertical energies** at the stage-1 LCM cell (nosym 4×4×4, 110 Ry, eV per f.u., relative to the LCM at identical settings). "→ LCM" means the converged state is the LCM (|ΔE| < 1 meV, same d moments and occupations).

| start (target) | KVMo 3/2 | KVMo 4/0 | KVMo 2/0 | KVCr 3/3 | KVCr 4/2 | KVCr 2/2 |
|---|---|---|---|---|---|---|
| FM6, d3/d3 FM with M = 6 | +0.552 (276 meV/ion = Gate 0) | +0.606 | +0.843 (gap −0.58) | +0.304 (152 meV/ion = Gate 0) | +0.291 | +0.393 |
| CT_AF, V(III)/M(II) LS AF, kicked | → LCM | → LCM | → LCM | pending | → LCM | → LCM |
| CT_AF_v, V side kicked only | → LCM | | | pending | | |
| CT_AF_ramp, converged at U = 6, then restarted at the target U | at U = 6 the CT state holds (V m_d −2.06, Mo +1.92, \|M\| 4.93); restart at 3/2 pending | pending | at U_V = 6 (no U on Mo) it already drifts back (\|M\| 6.6) | pending | pending | |
| CT_FM (M = 4), CT_HS_AF (M = 2), CT_HS_FM (M = 6) | pending | pending | pending | pending | pending | pending |
| MLS_AF: M(III) S = 1/2, M = −2 | **+0.685** (Mo m_d 0.72; gap 0.23) | pending | pending | +1.039 (Cr m_d 0.95; gap 0.94) | pending | |
| MLS_FM (M = 4) | +0.841 | pending | pending | +1.129 | pending | |
| VLS_AF: V(II) S = 1/2, M = +2 | +0.923 (V m_d −0.74) | | | **+0.889** (V m_d −0.85) | | |
| LSLS_AF (M = 0, fixed) | pending | | | → LCM | | |

**Kicked CT starts.**
- The kicked CT starts really were CT at first. For KVMo at 3/2 (`v_CT_AF.out`), the kick gave V Tr ns↓ 2.39 and Mo ↓ 1.58, with one Mo t2g↓ orbital at 0.99. By SCF iterations 4–5 the Mo t2g↓ occupation had fallen back to 0.10, and the run converged to E(LCM) + 5 μeV.
- At U_Mo = 0, where the Gate-0 free FM had charge-transferred (M_FM 4.91 / 5.91), the AF CT start also returns to the LCM.
- The CT drive seen in Gate 0 lives only in the FM arrangement: the fixed-M = 6 FM has gap −0.58 eV at 2/0, i.e. V t2g↑ above Mo t2g↓. It is not a competitor of the LCM ground state.

**Spin-state (non-CT) competitors.** The lowest is +0.685 eV/f.u. (Mo S = 1/2 ferrimagnet, KVMo) and +0.889 eV/f.u. (V S = 1/2, KVCr). Both are far above kT and above the FM. The d³ S = 3/2 configuration is robust on both sites.

**K variants**, LCM, base U, symmetry on, relaxed, ΔE per f.u.:
- **K is off-centre.** Started 0.3 Å along +[111], K relaxes *further*:
  - KVMo: to 0.67 Å toward the C-bound-metal corner of its cage, **−30.3 meV/f.u.**;
  - KVCr: to 0.57 Å, **−20.2 meV/f.u.**
- So the ideal F-43m Td K site is a saddle point. The LCM and its edges are untouched:

  | | M | gap (eV) | windows (eV) |
  |---|---|---|---|
  | KVMo, ideal | 0 | 1.39 | 0.97 / 1.87 |
  | KVMo, K off-centre | 0 | 1.34 | 0.98 / 1.89 |
  | KVCr, ideal | 0 | 2.02 | 2.03 / 1.11 |
  | KVCr, K off-centre | 0 | 1.98 | 2.01 / 1.19 |

  (SCF 5×5×5 grid; both edges V-majority.)
- K100, P2 and ALT2 (2-f.u. K-ordering) are still running.
- **Consequence for the phonon lane:** an unstable K Γ mode should appear at the F-43m geometry. Four equivalent ⟨111⟩ off-centre sites with a ~20–30 meV well suggest a K rattler / orientational glass at RT, not a structural problem for the LCM.

**Pending.**
- The remaining vertical starts.
- The relaxations: the lowest CT state per U point, and the lowest non-CT competitor at base U. At present these are KVMo MLS_AF and KVCr VLS_AF.
- **The HSE06 single points.** These carry the kill rule: a CT state below the LCM in HSE kills. If no CT start survives at a U point, the CT chain relaxes and HSE-tests the lowest fixed-M CT state (CT_FM / CT_HS_*), which cannot fall back to the M = 0 LCM.
- Collect with `python tracks/lcm/pba/taut/collect.py --fetch`. It rewrites `summary.json`; the HSE block and `kill_CT_below_LCM_in_HSE` fill in when the chains finish.
- **Expected finish:** vertical scans about 11:20–11:50; relax + HSE chains about 13:30–15:00 PDT. HSE with nosym fallback could take until about 17:00.

**Job history.**
- v1 (09:56) used a wrong `starting_ns_eigenvalue` mapping (§3.1); its 8 jobs were terminated at 10:12.
- v2 (10:12) crashed at the first kicked start on a parser bug.
- v3 (10:24) is the running set: `manifest.json` `jobs`, `terminated`, `superseded`.
- LCM, FM6 and CT_AF_free results carry over between versions. Kicked results do so only from v3 (`libver` 2).

## 1. Numbers in context

- **The LCM is an isolated minimum, robust to U.** The nearest competitors are the d3/d3 FM (exchange, +0.29–0.84 eV/f.u.) and the S = 1/2 spin states (≥ +0.68 eV/f.u.).
- **No V(III)/M(II) state was found as a local minimum** at U ≤ 4 on V with U_M ∈ {0, 2, 3}. The test was the spin-conserving V t2g↓ → M t2g↓ transfer that the Gate-0 band edges would favour.
- The Gate-0 statement "a charge-transferred ground state would be fatal" is therefore not triggered in PBE+U. The HSE check is in flight.


## 2. Question and kill rule

In the LCM, V(II) (N-bound, t2g³, S = 3/2) and M(III) (C-bound, t2g³, S = 3/2) are antiparallel. Both band edges sit in the V-majority (↓) channel: the VBM is V t2g↓ and the CBM is M t2g↓. So the Gate-0 gap (PBE+U 1.95 / 1.30 eV, HSE 2.09 / 1.52 eV for Cr / Mo) is itself the V→M charge-transfer excitation. A tautomer that moves this electron for good would compete with the LCM.

| competitor | V | M | M_cell (μB/f.u.) | why it matters |
|---|---|---|---|---|
| CT, low-spin M(II), AF | V(III) t2g↓², S = 1 | M(II) t2g↑³↓¹, S = 1 | 0 | spin-conserving V t2g↓ → M t2g↓ transfer; still compensated, but the windows change and moments become non-d3 |
| CT, low-spin M(II), FM | V(III) S = 1, ↑ | M(II) S = 1, ↑ | 4 | the free FM at U_Mo = 0 drifted toward it (M_FM 4.91) |
| CT, high-spin M(II), AF / FM | V(III) | M(II) t2g³eg¹, S = 2 | 2 / 6 | Jahn–Teller d4; implausible for C-bound M, tested anyway |
| M(III) low-spin, AF / FM | V(II) S = 3/2 | M(III) t2g↑²↓¹, S = 1/2 | −2 / 4 | "Mo low-spin S = 1/2" (task item ii) |
| V(II) low-spin, AF | V(II) t2g↓²↑¹, S = 1/2 | M(III) S = 3/2 | 2 | |
| both low-spin, AF | S = 1/2 | S = 1/2 | 0 | |
| FM d3/d3 (reference) | V(II) ↑ | M(III) ↑ | 6 | links this setup to the Gate-0 fixed-M FM numbers |

**Kill (task spec):** any charge-transfer state below the LCM in HSE06. In PBE+U, a CT state below the LCM at some U point is reported but is not by itself a kill (HSE decides, as in Gate 0).

## 3. Method

### 3.1 Why nosym, and how the starts are constrained

- Every competitor above is orbitally polarised: V(III) t2g², M(II) t2g⁴ or t2g³eg¹, and M(III)/V(II) S = 1/2 t2g³. QE symmetrises both ns and ρ with the crystal group. In the cubic cell the metal sites are Td, so any t2g polarisation would be averaged to 2/3 or 1/3 per orbital, giving a metallic, U-penalised (≈ U/3 per channel) and artificially high CT state.
- So **every PBE+U run uses `nosym = .true.`**, the LCM reference included, on the same 4×4×4 Γ-centred grid: 36 k with time reversal. It contains Γ, X, L and W.
- **Constrained starts:** site `starting_magnetization`, plus `starting_ns_eigenvalue` on all 5 d orbitals of both spins of every Hubbard site, plus `tot_magnetization` where M_cell ≠ 0.
  - **How the kick is applied (v2).** QE 7.5 applies `starting_ns_eigenvalue(m, s, type)` to the occupation matrix *after SCF iteration 1*, with m counted in ascending-eigenvalue order. So each constrained run is done in two passes:
    - pass 1: one SCF iteration with the same moments and `tot_magnetization`;
    - from pass 1's iteration-1 eigenvectors, each column is classified t2g or eg by its weight on z² + x²−y²;
    - pass 2: the real run, with the kick mapped onto those columns.
    - The spec is given per site and spin as (n_t2g, n_eg) occupied. All three t2g eigenvalues are set, the n_t2g most occupied to 1 and the rest to 0. Only n_eg eg eigenvalues are set (to 1); the others keep their covalent value.
  - **v1 was wrong.** v1 (`jobsrc/taut_v1.py`, md5 bbfc1d36, 09:56) assumed m = real harmonic m. Its CT_AF kick filled one eg and one t2g on V (`v_CT_AF.out`, iteration 2). Its 8 jobs were terminated at 10:12 (own jobs, `terminate_containers=True`; manifest `terminated`). The v2 jobs reuse the v1 LCM, CT_AF_free and FM6 runs, which carry no kick or are unaffected by it (FM6 reproduces Gate 0: 276.1 / 152.2 meV/ion vs 276.3 / 152.2).
  - **Variants.**
    - `CT_AF_v` kicks the V side only.
    - `CT_AF_free` uses moments only.
    - `CT_AF_ramp` converges the kicked CT state at U = 6 eV on every Hubbard element, then restarts at the target U from that density, its ns (occup.txt) and its wavefunctions, with no kick. This is the most forgiving test of whether a CT minimum exists at the target U.
  - The orbital *direction* inside the degenerate t2g set is whatever the iteration-1 eigenvectors are. With nosym the SCF is free to choose its own orbital order from there.
  - At U_Mo = 0 there is no Hubbard term on Mo, so only the V side is kicked.
- After SCF, a state is labelled by its ortho-atomic d moments m_d = Tr ns↑ − Tr ns↓ and its d counts N_d, both relative to the LCM at the same settings. It is CT if ΔN_d(V) < −0.15 e. It is "LCM-like" if both moments are within 0.25 μB of the LCM and |ΔN_d(V)| < 0.1. The labels come from the converged state, not from the start.

### 3.2 Settings

The Gate-0 stage-1 physics is kept: PseudoDojo SR NC, 110 Ry (ecutrho 440), ortho-atomic U, Gaussian 0.005 Ry, conv_thr 1e-7 × nat, mixing 0.3. A run that fails to converge is retried with mixing 0.1.

| job | U (V/M) | starts | then |
|---|---|---|---|
| `KVMo_U3_2_CT`, `KVCr_U3_3_CT` | base 3/2, 3/3 | LCM, CT_AF, CT_AF_v, CT_AF_ramp, CT_AF_free, CT_FM, CT_HS_AF, CT_HS_FM | relax the lowest CT state → HSE06 |
| `KVMo_U3_2_LS`, `KVCr_U3_3_LS` | base | LCM, FM6, MLS_AF, MLS_FM, VLS_AF, LSLS_AF | relax the lowest non-CT competitor → HSE06 |
| `KVMo_U4_0`, `KVMo_U2_0` | 4/0, 2/0 (no U on Mo) | LCM, FM6, CT_AF, CT_AF_ramp, CT_FM, CT_HS_AF, MLS_AF, MLS_FM | relax the lowest CT state |
| `KVCr_U4_2`, `KVCr_U2_2` | 4/2, 2/2 | same 8 | relax the lowest CT state |
| `KVMo_Kvar`, `KVCr_Kvar` | base, symmetry on, 5×5×5 | K off-centre and K ordering (§3.4) | |

All of these run on run_cpu64 with npool 16 (nosym) or 8 (HSE).

- **Vertical energies** are at the exact F-43m stage-1 LCM cell, i.e. the LCM's own relaxed geometry.
- **Relaxation:** an ionic relaxation at the fixed stage-1 cell, nosym, from the same constrained start. The JT and breathing distortions are internal coordinates, so the fixed cell still lets them develop.
  - Cell relaxation is not included. A ±1–2 % volume change with B ≈ 30–40 GPa is worth ≲ 15–30 meV/f.u.

### 3.3 HSE06 single points

**Geometry.** The relaxed competitor cell is symmetrised exactly with spglib/ase `refine_symmetry` (symprec 1e-3 Å), because QE HSE aborts with "lone vector" on ~1e-6 Å noise. Its space group, at several tolerances, is recorded.

**Steps.**
1. A PBE+U prep SCF runs at 90 Ry with the HSE k-grid, symmetry on and the same constrained start. In the distorted cell the orbital order is symmetry-allowed.
2. The prep must reproduce the relaxed d moments to within 0.3 μB. Otherwise it is retried with nosym, and the HSE then also runs with nosym.
3. HSE06 then runs with no U, nq 2³, ecutfock 180 Ry, Gygi–Baldereschi, x_gamma_extrapolation and conv_thr 1e-8 × nat (identical to `lcm_deep` / `pba_fmfix`). It starts from the prep wavefunctions and density (`startingwfc` / `startingpot = 'file'`), and keeps `tot_magnetization` where the state has M ≠ 0.

**Reference.** ΔE_HSE is taken against the existing HSE06 LCM at the cubic cell (`lcm/pba/hse/<m>_hse_LCM`, same 90 Ry / 4×4×4 / nq 2): −7840.41030 eV (KVCr) and −7335.94242 eV (KVMo).

**Caveat (smoke test).** Even with wavefunctions read from file, QE first runs a semilocal (plain PBE) SCF and only then the EXX loops. A CT state that is not a local minimum of plain PBE can relax back to the LCM in that first stage. The HSE sphere moments are therefore compared with the prep's (`m_sph` vs `m_sph_prep`). A collapse is reported as such, not as an HSE energy of the CT state.

### 3.4 K variants

All are LCM, base U, 110 Ry, symmetry on, ionic relaxation at the fixed cell:

- `K111` and `K100`: K displaced by 0.3 Å along [111] (R3m) or [100] (Imm2) in the 15-atom cell. Does K return to the Td site?
- `ALT2`: a 2-f.u. cell (a1 doubled) with the second K moved to the other tetrahedral-void set, giving R-3m K chains. It is compared with the ideal 2-f.u. cell `P2` on the same 3×5×5 grid.

## 4. Files

| path | content |
|---|---|
| `jobsrc/taut.py` | job script: vertical constrained starts, relaxation, symmetrisation, HSE chain |
| `jobsrc/kvar.py` | K variants |
| `jobsrc/smoke_hse_restart.py` | QE behaviour checks: hybrid restart from file; starting_ns + nosym + HUBBARD card |
| `submit_taut.py`, `collect.py` | submit and collect |
| `inputs/<job>/<job>.json` | job configs (exact stage-1 F-43m cells) |
| `results/` | cache of the fixed volume paths, plus `collect_snapshot.json` |
| `summary.json` | numbers quoted in this README |
