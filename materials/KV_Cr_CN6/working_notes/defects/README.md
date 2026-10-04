# PBA defect Gate 1: 4-f.u. defect cells of KV[Cr(CN)6] and KV[Mo(CN)6]

Track L (session 66190d), 2026-10-03. Spec: `tracks/lcm/trilemma/JUDGE.md` §2 (defect gates), `kinetic.md` §5 (D1–D3), `nonoxide_frameworks.md` §5 (D4). Gate 0 is in `../README.md`.

- Jobs: prefix `lcm/pba/def/`. Call ids, commands and fixed result paths are in `manifest.json`.
- Build: `python tracks/lcm/pba/defects/build_defects.py` writes `inputs/<id>.json`, `structures/<id>.xyz` and `build_info.json`.
- Submit: `python tracks/lcm/pba/defects/submit_defects.py [--only ID] [--resubmit] [--dry]`. It never cancels anything.
- Collect: `python tracks/lcm/pba/defects/collect.py --fetch [--npz]`. It reads fixed paths only (no listdir) and caches them in `results/`.

## 1. Question

Gate 0 says the d3/d3 cyanide double perovskites are unipolar Luttinger-compensated magnets in the ideal cell. In PBE+U, both edges sit in the V-majority channel:

| | win_VB / win_CB (eV) |
|---|---|
| KVCr | 2.02 / 1.15 |
| KVMo | 0.97 / 1.88 |

Real PBAs are never ideal. Gate 1 asks whether the four common defects keep that edge structure:
- the compensation;
- the windows;
- the spin of the band-edge states.

**Kill criteria** (JUDGE §2 and the task spec). Any one defect kills the gate if it:
- puts an **opposite-spin** state (M-majority channel) within **0.3 eV** of either host edge, or inside the gap; or
- cuts either window by **more than 50 %**.

Same-spin (V-majority) gap states do not kill. They are reported, because they move the Fermi level and the carrier type.

## 2. Cells

All cells are 4-f.u. conventional cubic cells of the stage-1 relaxed geometry, with 59–62 atoms. The source is `jobs/lcm/pba/s1a/<m>/results/<m>.json`, field `relaxed`:

| | a (Å) | V–N (Å) | M–C (Å) |
|---|---|---|---|
| KVCr | 10.6799 | 2.112 | 2.060 |
| KVMo | 10.8579 | 2.097 | 2.161 |

**Geometry checks.**
- That geometry equals the exact F-43m projection (`_sym.json`) to 1e-15 Å.
- `build_defects.py` replicates the 15-atom primitive cell with ase `make_supercell` and matches all 60 atoms to the analytic conventional cell to < 1e-6 Å.
- Every defect cell is translated so that the defect centre sits at the origin. Its point operations then have no fractional translation, so QE keeps them on any FFT grid.
- Start spins are **site-based**: C-bound site +3 μB, N-bound site −3 μB. Water, K, C and N start at 0.

| id | cell | atoms | space group (spglib, magnetic types) | M_cell expected (μB) |
|---|---|---|---|---|
| `<m>_P` | pristine K4V4M4(CN)24 | 60 | F-43m (QE: 24 ops, supercell) | 0 |
| `<m>_VAC` | one [M(CN)6] vacancy at (½ 0 0), see note 1 | 62 | R3 (3 ops) | **−3** (one M³⁺ removed) |
| `<m>_AS` | V/M antisite pair, see note 2 | 60 | P-43m (24 ops) | 0 (isospin) |
| `<m>_FLIP` | one CN flipped in place (+x CN of V(0)): V–C and M–N | 60 | Cmm2 (4 ops) | 0 |
| `<m>_KVAC` | one K vacancy at (¼ ¼ ¼), one hole | 59 | P-43m (24 ops) | **+1** if the hole is V-majority (V(III)); −1 if it is M-majority |
| `KVMo_KVAC_polV` | K vacancy, see note 3 | 59 | R3m | +1 |
| `KVMo_KVAC_polM` | K vacancy, see note 4 | 59 | R3m | −1 (forced) |

**Notes on the defect cells**
1. **VAC.**
   - The six exposed V–N sites are capped by H2O: V–O 2.13 Å, O–H 0.97 Å, HOH 104.5°, H pointing into the cavity, in a T-compatible pinwheel.
   - Three of the four K are removed, giving KV4M3(CN)18(H2O)6. The formal charge is 1 + 8 + 9 − 18 = 0.
   - This is the ordered 25 % Prussian-blue vacancy lattice. In the 10.7 Å cell, the three capped V each carry two trans V–OH2 bonds; the fourth V is untouched.
   - Starting minimum distances: H···N 3.03 Å, K···H 3.28 Å, intermolecular H···H 3.29 Å.
2. **AS.** M sits on the N-bound site at the origin (spin −) and V on the C-bound site at (½ ½ ½) (spin +). The cyanide orientation is kept. The two antisites are not bonded (√3/2 a apart), so the pair keeps Td about the origin.
3. **polV.** The hole is pre-localised on V(0 0 0): start −2 μB, with its six V–N bonds shortened by 0.06 Å. M is free.
4. **polM.** The hole is forced into the Mo-majority channel: Mo(½ 0 0) starts at +2 μB (Mo(IV) d²) with Mo–C bonds 0.05 Å shorter, and `tot_magnetization` = −1. E(polM) − min(E(KVAC), E(polV)) is the cost of an opposite-spin hole.
   - This test is KVMo only. The Mo-majority VB top lies only 0.97 eV below the VBM, against 2.02 eV for Cr in KVCr, where a Cr(IV) hole is implausible.

The formal charge check (K⁺ V²⁺ M³⁺ CN⁻ H2O⁰) gives 0 for every cell except KVAC, which is −1, i.e. one hole.

## 3. Protocol (job `jobsrc/pba_defect.py`)

The physics is the Gate-0 stage-1 protocol:
- PseudoDojo SR NC, 110 Ry (ecutrho 440);
- ortho-atomic U: V 3, Cr 3, Mo 2 eV;
- Gaussian smearing 0.005 Ry, conv_thr 1e-7 × nat, mixing 0.3.

| step | what |
|---|---|
| relax | Ions only, at the fixed stage-1 cell. BFGS with forc 1e-3 Ry/bohr and etot 1e-5 × nat Ry, symmetry kept. k = 2×2×2 shifted (1 1 1), which has 1–2 irreducible k in every cell. Preemption restarts from the last ATOMIC_POSITIONS. RELAX_MAXSEC: P 1 h, KVAC 2 h, pol 3 h, AS/FLIP 4 h, VAC 5 h. |
| scf | Fresh SCF at the relaxed positions, same k. Gives E, M, \|M\|, sphere moments and ortho-atomic Tr[ns] (d occupations and moments). |
| nscf | 4×4×4 Γ-centred. It contains the folded Γ, X, L and W of the fcc BZ, where the Gate-0 edges are: KVCr VBM at X, CBM at L. nbnd = N_el/2 + 40. A failed run is retried with CG. |
| projwfc | Löwdin, lsym false. The raw files (`raw/<id>/`: schema xml, projwfc.out, atomic_proj.xml.gz) stay on the volume, so `ANALYZE_ONLY=1` can redo the analysis. |
| analysis | (1) Occupation-based per-spin edges and windows, using the lane definition. (2) A band table for every band within 1.5 eV of the gap: energy range, occupation, group weights (V(N), M(C), V(C), M(N), C, N, K, O_w, H_w), defect-region weight and top atoms. (3) Deep alignment references (below). (4) `<id>_bands.npz`: E, occ, w_k and per-atom s/p/d weights of every state. |

**Alignment.** Defect states are placed relative to the **pristine host edges**.
- Each cell is aligned to `<m>_P` by the k-weighted mean energy of the CN 3σ manifold (≈ VBM − 17.6 eV). That mean is projected on C/N atoms more than 4.5 Å from every defect centre (`E_ref_CN3s_far`).
- The K 3p semicore level on far K atoms is a cross-check.

**Pools.** NPOOL is the largest divisor of the rank count that is ≤ min(NPOOL_MAX, 2 × irreducible k), counted with spglib inside the job: 32 ranks give ≤ 8 pools, 64 ranks ≤ 16.
- With 1–2 irreducible k, the relax/scf runs on 2–4 pools.
- The nscf (10–21 irreducible k) runs on 8–16 pools.

**References** (`jobsrc/pba_refs.py`, job `lcm/pba/def/refs`, 110 Ry, no U):
- bcc K (vc-relax) gives μ_K, the K-rich limit. Result: E = −817.8394 eV/atom, a = 5.280 Å (PBE).
- Isolated H2O in a 14 Å box gives μ_H2O.

**Formation energies (where meaningful)**

| defect | quantity |
|---|---|
| AS, FLIP | E(d) − E(P), same composition. Compare the full linkage isomer, +0.44 to +0.81 eV/f.u. in Gate 0. |
| KVAC | E_f = E(KVAC) − E(P) + μ_K(bcc): neutral defect (vacancy plus hole), K-rich limit |
| VAC | Only the open quantity E(VAC) − E(P) + 3μ_K − 6μ_H2O. The missing μ[M(CN)6] is fixed by synthesis (the hexacyanometallate reservoir), so no absolute E_f is claimed. For PBAs the vacancy content is a synthesis variable, not an equilibrium one. |
| KVAC variants | Relative energies among themselves |

**Cost.**
- Estimated about 37 node-hours on 64-core nodes in total: P and KVAC on run_cpu32; VAC, AS, FLIP and pol on run_cpu64; refs on run_cpu16.
- JUDGE estimated about 10 node-hours for 120-atom cells; that underestimated the 110-Ry, 220-band size.

## 4. Reading the output

`collect.py` prints, per cell:
- status and stage. Running relaxes are summarised from the `.out` at its fixed path: BFGS steps, E, M, force, SCF iterations.
- M_cell against the counting value, |M|;
- metal sphere and d moments by role;
- the cell gap and edge spins (V-maj / M-maj), and the in-cell windows;
- after alignment: the opposite-spin VB top relative to the host VBM (must be ≤ −0.3 eV) and the opposite-spin CB bottom relative to the host CBM (must be ≥ +0.3 eV);
- aligned windows as a fraction of pristine (must be ≥ 50 %);
- every in-gap or near-edge opposite-spin state, with its character;
- formation energies;
- the K-vacancy variant energies and hole bookkeeping (ΔN_d of V and M against pristine).

## 5. Expectations (written before results)

- **VAC.** M = −3. The capped V(II) keeps t2g³. H2O σ-donation is weaker than N-bound cyanide, so the capped-V t2g levels should rise, which is same-spin and makes the VB edge V-majority again. The risk is the M-majority channel near the CBM: the V t2g↑ minority (empty) on capped V could drop.
- **AS.**
  - V on a C-site sees the strong field of C (low-lying occupied V t2g, a V(II)–C bond).
  - M on an N-site sees weak-field N.
  - For KVMo, an N-bound Mo(III) has its t2g raised, in the Mo-majority channel (up). That is exactly the opposite-spin VB-edge risk.
  - KCrV (the full isomer) had a CB window of only 0.01–0.26 eV, so the CB side is the likely casualty.
- **FLIP.** A milder version of AS: one V–C and one M–N bond.
- **KVAC.** The hole should sit in V t2g (the host VBM, V-majority), giving M = +1. The polM variant prices an opposite-spin Mo(IV) hole.

## 6. Status

The status log, with numbers as they arrive, is §7, appended by the session.

## 7. Log

- 2026-10-03 05:18 PDT. Submitted 13 jobs (`manifest.json`).
  - At the first check (05:21), all were running with the predicted symmetry: P/AS/KVAC 24 ops with 1 k; FLIP 4 ops with 2 k; VAC 3 ops with 2 k; pol 6 ops with 2 k.
  - First-iteration M matched the start: VAC −3, KVAC +1, polM −1 (fixed), polV +1.
  - refs: bcc K finished at a = 5.280 Å, E = −817.8394 eV/atom.
