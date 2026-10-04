# PBA Gate 0: d3/d3 cyanide double perovskites (Prussian-blue analogues)

Track L (session 66190d), 2026-10-03. Spec: `tracks/lcm/trilemma/JUDGE.md` §2–3 (Gate 0), with input from `nonoxide_frameworks.md` §5 (D1–D4) and `kinetic.md` §5.

- Jobs: prefix `lcm/pba/`. Job list, call ids and fixed result paths are in `manifest.json`.
- Collect: `python tracks/lcm/pba/collect.py --fetch`. It reads fixed paths only (no listdir) and caches them in `results/`.
- Submit: `python tracks/lcm/pba/submit_pba.py s1`, then `... deep`. The deep phase reads the relaxed, symmetrized cells from `jobs/lcm/pba/s1/<id>/results/<id>_sym.json`.

## 1. What is being tested

In A[M_N][M_C(CN)6] the two magnetic sublattices are fixed by the cyanide linkage: M_C is C-bound (4b) and M_N is N-bound (4a). With d3/d3 (S = 3/2 on both ends) the cell is a Luttinger-compensated magnet (LCM) by construction. Gate 0 asks four questions:

- Is the computed ground state the d3/d3 LCM, with no charge transfer?
- Are the band edges unipolar, with windows of at least 0.3 eV, across U and in HSE?
- Is the V–Mo exchange at least 1.15× the V–Cr exchange, once calibrated on the measured V–Cr / Cr–Cr T_C ratio?
- Is the residual moment under SOC at most 0.1 μB/f.u.?

**Materials.** All cells are ideal: vacancy-free, anhydrous unless stated, primitive fcc.

| id | compound | role | builder bonds M_N–N / C–N / M_C–C (Å) |
|---|---|---|---|
| KVMo | KV(II)[Mo(III)(CN)6], F-43m, 15 atoms | **lead** (parent V1.37[Mo(CN)6]: 417 K, Bloch fit) | 2.08 / 1.16 / 2.17 |
| KVCr | KV(II)[Cr(III)(CN)6], F-43m, 15 atoms | experimental anchor (T_C 376 K, M_sat 0.125 μB/f.u.) | 2.08 / 1.16 / 2.04 |
| KCrV | KCr[V(CN)6], F-43m, 15 atoms | linkage isomer of KVCr (full-antisite limit) | 2.04 / 1.16 / 2.10 |
| CrCr | Cr(III)[Cr(III)(CN)6], Fm-3m, 14 atoms | T calibrant (T_C 240 K) | 2.04 / 1.16 / 2.01 |
| KVCr_2H2O | KV[Cr(CN)6]·2H2O, P1, 21 atoms | hydration robustness of the KVCr edges | KVCr builder plus 2 H2O |

**Calibrant composition.** Schart et al. (Inorg. Chem. 2024) report Cr[Cr(CN)6]0.996(7), made by a vacancy-suppressed NMF route: Fm-3m, a = 10.42 Å, T_C = 240 ± 10 K (neutron), non-compensation 0.038 μB/f.u. Charge balance at 0.996 requires N-bound Cr(III), so the cell is Cr(III)–NC–Cr(III), S = 3/2 on both ends, and it is the matched d3/d3 calibrant.

- Caveat: the classic Cr(II)3[Cr(III)(CN)6]2·xH2O (Gadet 1992) also orders at about 240 K.
- If the Schart assignment were wrong, the calibrant would not be d3/d3. The ratio test then loses meaning; it does not merely shift.

**Cells.** All cells come from `jobsrc/pba_geom.py`, which reproduces `trilemma/nonoxide_pba/*_prim.cif` to within 5e-4 Å; the differences come only from the CIF's rounded xN and xC.

- The cells are built at machine precision: a 60° rhombohedral cell from `ase.cellpar_to_cell`, with N and C coordinates exactly ±x.
- After the vc-relax, `pba_chain.py` projects the relaxed cell back onto exact F-43m / Fm-3m:
  - a is taken from the volume;
  - xN and xC are averaged over all 18 coordinates.
- That projection is the geometry for the U grid, HSE and SOC. QE HSE aborts with `sym_rho_init_shell: lone vector` on about 1e-6 cell noise (see TRACKL_LOG 2026-10-02 12:25).

**LCM convention.** The C-bound metal starts at +3 μB and the N-bound metal at −3 μB. FM is added automatically.

- In KVCr and KVMo the "V-majority channel" is spin index 1. `collect.py` labels the channels from the actual sign of the V moment.
- For the isomer, the N-bound metal is Cr.

## 2. Jobs

| job | fn | what it tests | results (volume, fixed) |
|---|---|---|---|
| `lcm/pba/s1/KVCr` | run_cpu32 | Steps 1–2. vc-relax in LCM at U_V = U_Cr = 3. LCM and FM SCF, LCM dense nscf. Then a 2-D U grid (U_V, U_Cr) ∈ {2,4}² at the symmetrized relaxed cell: LCM SCF + nscf, plus FM SCF | `results/KVCr.json`, `KVCr_sym.json`, `KVCr_U{2,4}_{2,4}.json` |
| `lcm/pba/s1/KVMo` | run_cpu32 | Same, at U_V 3 / U_Mo 2 (Mo 4d). U grid (U_V, U_Mo) ∈ {2,4} × {0,2}; U_Mo = 0 means no Hubbard term on Mo | `KVMo.json`, `KVMo_sym.json`, `KVMo_U{2,4}_{0,2}.json` |
| `lcm/pba/s1/KCrV` | run_cpu32 | Step 5 (isomer). vc-relax LCM at U 3/3, LCM + FM SCF, LCM nscf. Grid energies at (U_V, U_Cr) ∈ {2,4}², no nscf. Gives E(isomer) − E(KVCr) at five U points | `KCrV.json`, `KCrV_U*.json` |
| `lcm/pba/s1/CrCr` | run_cpu32 | Calibrant. vc-relax LCM at U_Cr = 3, LCM + FM + nscf. U_Cr ∈ {2,4}, so the calibration ratio also exists on the U diagonal | `CrCr.json`, `CrCr_sym.json`, `CrCr_U{2,4}.json` |
| `lcm/pba/s1/KVCr_2H2O` | run_cpu32 | JUDGE §3 hydrated cell. vc-relax LCM at U 3/3 (k 3×3×3; the P1 cell makes the 15-atom k-grid too expensive). LCM + FM SCF, nscf 6×6×6 | `results/KVCr_2H2O.json` |
| `lcm/pba/hse/{KVCr,KVMo,CrCr}_hse_LCM` and `CrCr_hse_FM` | run_cpu64 | Step 3. HSE06 with no U on the symmetrized relaxed cells. HSE is the arbiter for windows, gap, moments and the J ratios | `deep_hse.json` |
| `lcm/pba/hse/{KVCr,KVMo}_hse_FMfix` | run_cpu64 | HSE06 d3/d3 FM reference with M fixed at 6 μB (`tot_magnetization`, `jobsrc/pba_fmfix.py`) | `results/<id>.json` |
| `lcm/pba/s1a/FMfix` | run_cpu32 | PBE+U FM with M fixed at 6 at all five KVCr and five KVMo U points, at the stage-1 geometry, 110 Ry | `results/{KVCr,KVMo}_<U>_FMfix.json` |
| `lcm/pba/soc/KVMo_soc` | run_cpu32 | Step 4. `lcm_deep soc`: noncollinear SOC, dojo_fr, PBE+U 3/2, LCM with m ∥ z and m ∥ x. Gives the spin residue \|M_spin\| and the MAE. Then `pba_orbm.py`, a best-effort cell orbital magnetization (lorbm nscf) | `deep_soc.json`, `orbm.json` |

**Inside each s1 job** (`jobsrc/pba_chain.py`):
1. The unchanged `lcm_stage1.py` (pivot/garnet copy, md5 be16bdd5…) runs the base candidate.
2. The relaxed cell is symmetrized and written to `<id>_sym.json`.
3. `lcm_stage1.py` runs the U-grid candidates.

Every step resumes after preemption. In `results/<id>*.json`, `runs.{LCM,FM}` holds:
- `energy_eV`, `total_mag`, `abs_mag`, `site_moments`;
- `glob_nscf` / `glob_scf` with `vbm_spin`, `cbm_spin`, `win_VB`, `win_CB`, `gap_up`, `gap_dn`;
- `edges_*`.

## 3. Settings, and why

**Plane-wave cutoff: 110 Ry for all PBE+U runs.**
- JUDGE §2 specifies "vc-relax LCM at 110 Ry", and it is the lane's geometry rule: the 75/90-Ry P4 artefact of YBMFO disappeared at 110 Ry.
- The PseudoDojo "normal" hint for Cr is about 47 Ha (94 Ry), already above 90 Ry. The short C≡N bond and the vc-relax stress (Pulay) make 110 the safe choice.
- The cost is small for 15 atoms.
- The U-grid single points use the same cutoff, so every PBE+U ΔE is at one basis.

**HSE and SOC: 90 Ry (HSE ecutfock 180 Ry), at the 110-Ry geometry.** Energy differences between magnetic orders at a fixed geometry converge much faster than stress. 90 Ry keeps the EXX cost down.

**k-points**

| use | spacing | grid |
|---|---|---|
| SCF | KSCF 0.25 Å⁻¹ | 5×5×5 |
| dense nscf | **KNSCF 0.14** | **8×8×8** |
| HSE | 0.27 | 4×4×4 |
| SOC | 0.25 | 5×5×5 |
| orbital magnetization (lorbm) | 0.30 | 4×4×4, no symmetry |

- Dense nscf: the spec asks for ≤ 0.16. 0.16 gives a 7×7×7 grid, which contains neither X nor L. 8×8×8 contains Γ, X, L and W.
- HSE: the k-grid must be divisible by nq = 2. The grid is 4×4×4 for all three relaxed cells; `submit_pba.py` checks this and rescales if needed.

**U values.**
- Base: V 3, Cr 3, Mo 2 eV, ortho-atomic.
- Grid: (U_V, U_Cr) ∈ {2,4}² and (U_V, U_Mo) ∈ {2,4} × {0,2}. The grid is 2-D, not joint, because the windows follow the U_V − U_M offset.
- The isomer and the calibrant are added on their matching points.

**Scheduling.**
- All 15-atom PBE+U jobs run on run_cpu32 with NPOOL = 8 and NPOOL_NSCF = 16. The SCF IBZ has 10 k (20 with spin); the nscf IBZ has 29 (58).
- HSE runs on run_cpu64 with npool 8 (lcm_deep default NC // 8).

**Estimated cost:**

| block | 64-core node-h |
|---|---|
| s1 | ≈ 3 |
| hydrated | ≈ 1.5 |
| HSE | 6 × about 0.5 |
| SOC | ≈ 0.8 |
| **total** | **≈ 7–9** |

## 4. Pass / kill criteria (copied from JUDGE.md §2–3)

**Pass (all required):**
- Unipolar edges in the V-majority channel, with win_VB and win_CB ≥ 0.3 eV at every U-grid point and in HSE.
- HSE gap ≥ 1.0 eV.
- Integer moments.
- Calibration ΔE_FM−LCM(VCr) / ΔE(CrCr) = 1.57 ± 0.25 (measured 376/240 K). Otherwise the T scale is not trusted and no ≥ 400 K claim is made.
- ΔE(VMo) / ΔE(VCr) ≥ 1.15, in both PBE+U and HSE (⇒ T ≳ 430 K).
- SOC residue ≤ 0.1 μB/f.u.
- Isomer ≥ +0.8 eV/f.u.

**Kill:**
- Either HSE window < 0.2 eV, or the edges are bipolar with both windows < 0.3 eV.
- A charge-transferred or fractional-moment ground state in HSE, or an HSE gap < 0.5 eV.
- ΔE ratio < 1.0 (predicted T < 376 K).
- Any defect gate puts opposite-spin states inside a window. The defect gates are not part of Gate 0; they run only on a pass.

**KVCr·2H2O (JUDGE §3):**
- Pass: windows ≥ 0.3 eV on both edges, unipolar, in HSE and over the whole U grid; water states ≥ 0.5 eV away from both edges; predicted Δ_opt = win_VB + win_CB ≥ 0.6 eV (the MCD observable).
- Kill: as for item 1. If the hydrated cell shifts either window by > 50 %, kill the "existing sample" claim.

**Step-1 sanity checks (this README):**
- M_cell = 0.00 in the LCM.
- |m| 2.6–3.0 μB on both metals.
- The FM single point must give M = 6.00 μB, the d3 + d3 insulating value.
  - A V(III)/Cr(II) or V(III)/low-spin Mo(II) transfer would give 4.00 μB, with one moment near 2 μB.
  - So M_FM separates covalent moment reduction (which is fine) from charge transfer (which is fatal for the windows).
- `collect.py` flags any violation.
- **Which moment.** QE's "magnetic moment per site" is integrated on very small spheres (R ≈ 0.11 alat ≈ 0.86 Å), so it understates the d moments. At the first relax SCF of KVMo the sphere moments are V −1.90 / Mo +1.42, against ortho-atomic d-occupation moments (Tr ns up − down) of 2.58 / 2.44 μB with d counts 3.49 / 4.69, the covalency of σ donation and π back-bonding.
  - `collect.py` applies the 2.6–3.0 window to the Hubbard d moments. A value slightly below 2.6 alongside M_FM = 6.00 is covalency, not charge transfer.
  - HSE has no Hubbard projector, so it is judged on M_LCM = 0, M_FM = 6 and its sphere moments compared with PBE+U at the same geometry.

## 5. Adaptations and caveats

- **Fixed-moment FM reference (added 2026-10-03 ~07:10).** The unconstrained FM is not the d3/d3 Heisenberg reference here, so the J ratios use FM with M fixed at 6 (`tot_magnetization`, two Fermi levels).
  - In plain PBE (the first stage of an HSE run), the KVMo FM converged to M = 2.57 μB: V t2g↑ → Mo t2g↓ transfer, metallic. The KVCr FM start collapsed to M = 0.
  - In PBE+U at U_Mo = 0, M_FM = 4.91 / 5.91.
  - `KVCr_hse_FM` and `KVMo_hse_FM` were therefore terminated, with `terminate_containers=True` (see `manifest.terminated`), and replaced by the `*_hse_FMfix` jobs.
  - `CrCr_hse_FM` keeps the unconstrained FM, whose PBE start is at M = 6.00.
  - `collect.py` reports both the free-FM and the fixed-M ratios. The gate uses the fixed-M ratios.
- **Cell orientation (05:50).** The first four stage-1 jobs (`lcm/pba/s1/{KVCr,KVMo,KCrV,CrCr}`) were built in the `cellpar_to_cell` orientation. QE then found only 6/24 (12/48) symmetry operations, so they were terminated after about 6 minutes and rerun as `lcm/pba/s1a/*`, with cubic axes along the Cartesian axes. `lcm/pba/s1/KVCr_2H2O` is P1, so the orientation does not matter for it and it was kept.

- **Water states.** Gate 0 has no PDOS, so the hydrated cell's water states are judged indirectly: the gap and windows against anhydrous KVCr at the same U. If the change is borderline, a garnet_pdos-style projwfc run is the follow-up. The water model is ad hoc (P1, two H2O in the empty 4d void, O–H···N donors) and relaxed from the builder cell.
- **Isomer ΔE.** It compares two different cells (different a and FFT grids). Session 3 found cross-cell offsets of about 20 meV/cell at 75 Ry (inbox 2026-10-02), smaller at 110 Ry. That is negligible against the 0.8 eV bar.
- **SOC orbital moment.**
  - pw.x reports only the spin magnetization. `pba_orbm.py` tries the modern-theory `lorbm` nscf (NC pseudos, nosym), on a best-effort basis.
  - If it fails, the orbital part of the residue is bounded only by the ligand-field g-shift estimate (kinetic.md §4.2: ΔM ≈ 0.02–0.05 μB/f.u. for V–Mo) and by Mo XMCD M_L = −0.04 μB.
- **T predictions** use ratio scaling, T(VMo) = 376 K × ΔE(VMo)/ΔE(VCr). That is only valid if the VCr/CrCr calibration passes. Mean-field T ∝ J for equal S on the same lattice.
- **Not in Gate 0:** defect gates (K vacancy, capped [Mo(CN)6] vacancy, V↔Mo swap, CN flip; 8-f.u. cells) and Cr[Mo(CN)6]. These run only on a pass.

## 6. Startup verification, 2026-10-03 04:33–05:00 PDT (Track L, session 66190d)

Four checks between 04:33 and 05:00. For each job I read its status and stdout and the pw.x `.out` files at fixed paths (no listdir).

- No job crashed at startup.
- No QE input error and no NPOOL error. 8 pools on 32 or 64 ranks everywhere.
- Symmetry is right: 24/48 operations in the 15/14-atom cells, 8 under noncollinear SOC.
- Every LCM start has M_cell = 0.00, and every FM or FMfix start has M = 6.00.
- One sub-step failed: the best-effort lorbm nscf at the end of `KVMo_soc`. It is fixed and resubmitted as `KVMo_soc_r`.

### 6.1 Job status at 04:57

| job | state | check |
|---|---|---|
| s1a/KVCr | running, 1.1 h | Base and grid points (2,2), (2,4), (4,2) done. (4,4) running. |
| s1a/KVMo | **done**, rc 0, 0.73 h | Base and all 4 grid points done. |
| s1a/KCrV | **done**, rc 0, 0.52 h | Base and all 4 grid points done. |
| s1a/CrCr | **done**, rc 0, 0.58 h | Base, U 2 and U 4 done. |
| s1/KVCr_2H2O | running, 1.15 h | P1 cell (no symmetry), 14 k. vc-relax at BFGS step 21: M = 0.00, sphere moments V −2.04 / Cr +2.33, P ≈ −1.5 kbar. The energy is still falling by about 1e-3 Ry per step (RELAX_MAXSEC 4 h). |
| s1a/FMfix | **done**, rc 0, 0.44 h | All 10 points converged at M = 6.00. |
| hse/CrCr_hse_LCM | running, 0.70 h | 48 ops, 8 k, nq 2³. M = 0.00. EXX dexx 4e-7 Ry, so it is nearly finished. |
| hse/CrCr_hse_FM | running, 0.70 h | Unconstrained FM stays at M = 6.00. dexx 3e-7 Ry. |
| hse/KVMo_hse_LCM | running | 24 ops, 8 k. M = 0.00, \|M\| 6.77. Third EXX outer loop, dexx 1e-3. |
| hse/KVCr_hse_LCM | running | 24 ops. M = 0.00, \|M\| 6.99. dexx 8e-4. |
| hse/KVCr_hse_FMfix | running | M fixed at 6.00. PBE start: V +1.94 / Cr +2.21 (sphere), \|M\| 6.15. |
| hse/KVMo_hse_FMfix | running | M fixed at 6.00. PBE start: V +1.91 / Mo +1.41 (sphere), \|M\| 6.08. |
| soc/KVMo_soc | **done**, rc 0, 0.79 h | SOC z and x converged. The pba_orbm lorbm nscf failed after 1 s (see 6.3). |
| **soc/KVMo_soc_r** (new) | running | Rerun of the orbital step only (`fc-01M40SW1JE4EH8BHY2Z7FQ0VDS`, run_cpu32, submitted 04:54:55). The SCF is iterating: 8 ops, 24 k, M = 0. |

### 6.2 Numbers from finished parts (`collect.py --fetch`, 04:56)

PBE+U, 110 Ry. Windows come from the 8×8×8 nscf and are in the V-majority channel unless noted. ΔE = E(FM) − E(LCM) per magnetic ion.

**KVCr**, a = 10.680 Å. Unipolar V-majority, M_FM = 6.00 at every point.

| U point | gap (eV) | win_VB / win_CB (eV) | ΔE (meV) |
|---|---|---|---|
| 3/3 (base) | 1.95 | 2.02 / 1.15 | 152.2 |
| 2/2 | 1.42 | 1.95 / 1.41 | 196.6 |
| 2/4 | 1.84 | 2.63 / 0.91 | 159.8 |
| 4/2 | 2.04 | 1.39 / 1.40 | 145.7 |
| 4/4 | pending | | |

Hubbard d moments at base U: V 2.70 / Cr 2.93.

**KVMo**, a = 10.858 Å

| U point | gap (eV) | win_VB / win_CB (eV) | free FM: M_FM, ΔE (meV) | fixed FM (M = 6): ΔE (meV) |
|---|---|---|---|---|
| 3/2 (base) | 1.30 | 0.97 / 1.88 | 6.00, 276.3 | 276.3 |
| 2/0 | 0.66 | 0.41 / 2.13 | **4.91**, 356.4 | 423.2 |
| 2/2 | 1.04 | 1.19 / 1.70 | **5.97**, 321.3 | 321.8 |
| 4/0 | 1.11 | **0.05** / 2.54 | **5.91**, 301.8 | 303.5 |
| 4/2 | 1.60 | 0.71 / 2.01 | 6.00, 238.3 | 238.3 |

- Hubbard d moments at base U: V 2.61 / Mo 2.45.
- The (4,0) point fails the "≥ 0.3 eV at every U point" test in PBE+U, so HSE decides.
- In the fixed-M FM at (2,0), the minority edge lies 0.58 eV below the majority edge ("gap" −0.58 eV). So at U_Mo = 0 the d3/d3 FM is unstable to V→Mo transfer. Treat the 423 meV value as a constrained number.

**KCrV isomer**, E(KCr[V(CN)6]) − E(KVCr), each at its own relaxed cell:

| U point | ΔE_iso (eV/f.u.) |
|---|---|
| 3/3 (base) | **+0.626** |
| 2/2 | +0.540 |
| 2/4 | +0.442 |
| 4/2 | +0.810 |
| 4/4 | pending (needs KVCr 4/4) |

- The base-U value is below the 0.8 eV bar, so the isomer check **fails** in PBE+U.
- The isomer's own CB window is only 0.01–0.26 eV.

**CrCr**, a = 10.487 Å

| U_Cr | ΔE (meV) |
|---|---|
| 3 | 91.7 |
| 2 | 107.6 |
| 4 | 80.2 |

The edges are bipolar (win 0.10 / 0.26 eV), which is acceptable for a calibrant.

**Exchange ratios**

| ratio | value | note |
|---|---|---|
| VMo/VCr, base (free FM = fixed FM) | 1.815 | ratio scaling gives 682 K |
| VMo/VCr, all same-U_V pairs (fixed FM) | 1.64–2.65 | |
| VCr/CrCr, base | 1.660 | inside the 1.57 ± 0.25 target |
| VCr/CrCr, U = 2 diagonal | **1.828** | just above the 1.82 upper limit |
| VCr/CrCr, U = 4 diagonal | pending (KVCr 4/4) | |

**SOC (KVMo LCM, PBE+U 3/2, dojo_fr)**
- |M_spin| = 0.000 μB/f.u. with m ∥ z and with m ∥ x.
- Sphere moments: V −1.942 / Mo +1.409.
- Gap 1.37 eV.
- E(x) − E(z) = 0.0003 meV/cell. See 6.3 for why this is not a real MAE.

**HSE (preliminary, not converged)**
- CrCr: at dexx ~4e-7 Ry, E(FM) − E(LCM) ≈ 0.01441 Ry/cell, which is 98 meV/ion (PBE+U at U 3 gave 91.7).
- KVCr and KVMo are still in their EXX outer loops.

### 6.3 Problems found and actions

1. **lorbm step crashed at its startup (fixed and resubmitted).**
   - The `pba_orbm.py` nscf in `KVMo_soc` stopped after 1.2 s with `Error in routine iosys (1): Berry Phase/electric fields only for insulators!`. QE requires `occupations = 'fixed'` with lorbm, and the nscf had inherited Gaussian smearing.
   - Fix: the nscf `system` now sets `occupations 'fixed'`, `degauss 0`. The SCF is unchanged. pba_orbm.py md5 is now `41ac51c7…`; the v1 md5 `d556f22c…` is kept in the manifest.
   - I checked the generated namelist locally before submitting.
   - The orbital step alone was resubmitted as `lcm/pba/soc/KVMo_soc_r`, estimated at 0.7 h on cpu32. The SOC z/x runs are not repeated.
   - `collect.py` now prefers `soc/KVMo_soc_r/orbm.json`.
   - Risk: QE may also reject lorbm together with Hubbard U or noncollinear spin. That would show up as an error recorded in `orbm.json`, and the gate would fall back to the g-shift bound.
2. **`collect.py` bug (fixed).**
   - `lcm_deep.py` writes `deep_<mode>.json` flat while it runs. At exit it rewrites the file as `{"id", "<mode>": {...}}`.
   - After `KVMo_soc` finished, the SOC block read back as `None`. Every HSE result would have done the same once its job finished.
   - `load()` now unwraps the `hse` / `soc` key.
3. **The SOC MAE is zero by symmetry.**
   - The cell has its cubic axes along x, y, z, so m ∥ x and m ∥ z are equivalent ⟨100⟩ directions in F-43m. E(x) − E(z) = 0.0003 meV is the expected null, not a measured anisotropy.
   - The spin-residue gate is unaffected.
   - A real MAE (and a second, inequivalent residue direction) needs m ∥ [111] or [110]. I have not submitted this; it is not a Gate-0 criterion.
4. **Record inconsistencies, not corrected in the original entries.**
   - The `terminated` entries in the manifest (and §5 above) give termination times of "~05:50" and "~07:10". Those times had not yet happened.
   - The job registry shows the replacement submissions at 03:46 (s1a) and 04:28 (FMfix). So the terminations were at about 03:46 and about 04:28 PDT.
5. **Small corrections to the submitter's summary.**
   - Besides the U_Mo = 0 points, KVMo (2,2) also has a slightly charge-transferred free FM (M_FM 5.97, gap_FM 0.10). The fixed-M ΔE there is essentially the same (321.8 vs 321.3 meV).
   - The VCr/CrCr calibration passes at base U but sits just outside the target on the U = 2 diagonal (1.83).

**Not done.** No `mrun cancel` was issued, no other lane's files were edited, and TRACKL_LOG.md and LANES.md were not touched.
