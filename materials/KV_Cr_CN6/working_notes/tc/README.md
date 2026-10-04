# PBA T_C model: Heisenberg model of the d3/d3 cyanide double perovskites

Track L (session 66190d), 2026-10-03. Folder `tracks/lcm/pba/tc/`. Inputs are the Gate-0 DFT energies (`../collect.py`) and three new 2-f.u. supercell jobs (`lcm/pba/tc/J2_*`, §3). Every number below is in `summary.json`. Rerun with `python tracks/lcm/pba/tc/tc_model.py --fetch`.

## 0. Bottom line

| quantity | value |
|---|---|
| **A.** Stoichiometric KV[Mo(CN)6], DFT exchange ratio calibrated on KVCr and CrCr (PBE+U, NNN included) | **673 ± 77 K** (U pairings 609–750 K) |
| **B.** The same model calibrated directly on the V–Mo parent V1.37[Mo(CN)6] (417 K at p = 0.73) | **487 ± 49 K** (hard floor 377 K) |
| **C.** Combination of A and B (Birge-inflated; χ² = 4.1 for 1 dof) | **541 ± 85 K** |
| Calibration check T(KVCr)/T(CrCr), model vs measured 376/240 = 1.567 | −0.9 % (base U), +9 / −9 % (U = 2 / 4 diagonal), +10 % (HSE06) |
| f = T_exp / T_model, classical MC (S² scale), base U | 0.97 (CrCr), 0.98 (KVCr) |
| Model for the parent (p = 0.73) using A | 550 K, against 417 K measured (×0.76) |

- **Every route clears 400 K for the 1:1 phase.** That holds for A, for B and for C. B's lowest value (377 K) needs both the 340-K measured floor of the parent and the weakest dilution law.
- **The DFT ratio (A) is probably too high.** Run on the measured parent composition (27 % [Mo(CN)6] vacancies), A over-predicts T_C by ×1.32. For comparison, the crystalline V–Cr samples scatter around the same model by only ±13 % (§6). Read A as an upper scenario and B as the conservative one.
  - The experiment-implied V–Mo/V–Cr T ratio is 1.21–1.35.
  - The PBE+U ratio is 1.77 (1.60–1.97 across U).
  - Both clear the Gate-0 bar of ≥ 1.15.
- **[M(CN)6] vacancies lower T_C slowly** (§5):
  - The C-sublattice percolation threshold is p_c ≈ 0.135.
  - 5 % vacancies cost 3 %, and 10 % cost 6 %.
  - Each vacancy, however, costs the LCM its zero: M = 3v μB/f.u.
- **Beyond-NN exchange is small, antiferromagnetic and frustrating**, and it matters at the 3–10 % level:
  - J2/J1 is 0.01–0.05 across the cube faces;
  - J3 ≈ 0;
  - J4 (linear M–CN–M'–NC–M) ≈ 0.007 J1 (CrCr).
  - Including it brings the two calibrants into agreement: the VCr/CrCr ratio test moves from +6 % to −1 %, and A from KVCr and from CrCr agree to within 1 %.
- **The finite-T compensation is broken by J2C ≠ J2N, but only weakly.** The classical net moment stays below 0.02 μB/f.u. at every T < T_C (§7).

KVMo HSE06 (LCM) was still running at the time of writing; its job had restarted at 05:41. Once it lands, `tc_model.py` adds the HSE pairing automatically.

## 1. Model

**Lattice.** All metals of A[M_N][M_C(CN)6] sit on a simple-cubic net with spacing d = a/2 and rock-salt colouring. The N-bound metal is at even parity and the C-bound metal at odd parity.

| shell | vector (d) | n | pair | path |
|---|---|---|---|---|
| J1 | (1,0,0) | 6 | M_C–M_N | M–C≡N–M′, 5.3 Å |
| J2C, J2N | (1,1,0) | 12 | same sublattice | M–CN–M′–NC–M, 90° at M′, 7.5 Å (across the cube faces) |
| J3 | (1,1,1) | 8 | M_C–M_N | three bridges, 9.2 Å |
| J4C, J4N | (2,0,0) | 6 | same sublattice | linear M–CN–M′–NC–M, 10.6 Å |

**Hamiltonian.** H = Σ_<ij> J_ij S_i·S_j, with J > 0 antiferromagnetic and S = 3/2 on both sublattices.

**Mapping.** Collinear (broken-symmetry) DFT states are mapped with S² (Ising energies):

E(FM) − E(LCM) = (6 J1 + 8 J3) S² per magnetic ion.

- J2 and J4 pairs are parallel in both FM and LCM, so the 15-atom ΔE never contains them. The supercells in §3 resolve them.
- FM is the fixed-M = 6 d3/d3 reference (`dEfix`) wherever it exists.

**Temperature levels**

| level | definition | pure-J1 value |
|---|---|---|
| MC | classical spins of length S. Heat bath + over-relaxation + Wolff embedding clusters (valid for any coupling sign). Binder crossings at L = 12/16/24. | 1.4423 ± 0.006 J1S² (literature 1.4430) |
| Q | quantum-corrected classical, T_MC · S(S+1)/S² · θ_Q/θ_cl. θ_Q(3/2) = 1.395 ± 0.013, interpolated between the S = ½ sc-AF QMC value (T_N = 0.946 J, Sandvik 1998) and the classical limit | T_Q/T_MC = 1.61 |
| RPA | Tyablikov, quantum S = 3/2: kT = S(S+1)/3 · ⟨A/(A² − B²)⟩⁻¹, with A = J_AB(0) − J_AA(0) + J_AA(q) and B = J_AB(q). BZ sum normalised to the exact Watson integral. | 1.3189 J1 S(S+1) |
| MF | two-sublattice molecular field | 2 J1 S(S+1) |

All the compounds share S = 3/2 and the same lattice, so the S² vs S(S+1) and quantum factors cancel in the calibrated predictions. They change only the absolute f. After calibration, the levels differ only through how they weigh the frustrating J2, by ≤ 1 %.

## 2. DFT exchange (E(FM) − E(LCM), meV per magnetic ion → J1 with NNN ratios from §3)

| compound | point | ΔE | J1 (meV) | T_MC (K) | T_Q | T_RPA | T_MF |
|---|---|---|---|---|---|---|---|
| KVCr | 3/3 base | 152.2 | 11.28 | 391 | 629 | 593 | 927 |
| KVCr | 2/2, 2/4, 4/2, 4/4 | 196.6, 159.8, 145.7, 121.9 | 14.6, 11.8, 10.8, 9.0 | 505, 410, 374, 313 | | | |
| KVCr | HSE06 | 181.1 | 13.42 | 465 | 749 | 705 | 1103 |
| KVMo | 3/2 base | 276.3 | 20.41 | 693 | 1116 | 1045 | 1649 |
| KVMo | 2/2, 4/2 | 321.8, 238.3 | 23.8, 17.6 | 807, 597 | | | |
| KVMo | 2/0, 4/0 (U_Mo = 0) | 423.2, 303.5 | | flagged | | | |
| CrCr | U 3 base | 91.7 | 6.79 | 249 | 401 | 379 | 580 |
| CrCr | U 2, U 4 | 107.6, 80.2 | 8.0, 5.9 | 292, 218 | | | |
| CrCr | HSE06 | 98.0 | 7.26 | 266 | 429 | 405 | 620 |

**Why the U_Mo = 0 points are flagged.**
- At U_Mo = 0 the free FM is charge-transferred, with M = 4.91 and 5.91.
- The fixed-M FM has a gap of −0.58 / 0.01 eV.
- So the d3/d3 Heisenberg reference is unstable there. These points are reported but kept out of the error bars.

## 3. Beyond-nearest-neighbour exchange (2-f.u. supercells; jobs `lcm/pba/tc/J2_{KVCr,KVMo,CrCr}`)

**Method** (`jexch_submit.py`, runner `../jobsrc/pba_fmfix.py` unchanged):
- Two index-2 supercells of the relaxed, symmetrized stage-1 cell:
  - type I, tetragonal (1,−1,0)/(1,1,0)/(0,0,2) d;
  - type II, rhombohedral (2,1,1)/(1,2,1)/(1,1,2) d.
- Configurations: 6 in type I (FM, LCM, flipC, flipN, sc-A, sc-C) and 4 in type II (FM, LCM, alt, flipC).
- Every configuration has its total moment fixed at the d3/d3 count.
- PBE+U at each compound's base U, 90 Ry, k 0.30 Å⁻¹.
- `tc_model.py` recounts every pair coefficient from the actual geometry and solves least squares with one E0 per cell. It uses only within-cell differences, because the two k-meshes differ by about 1 meV/cell.

**Checks**
- FM − LCM in the supercell is 152.20 / 276.53 / 91.72 meV/ion, against 152.24 / 276.3 / 91.7 meV/ion in the 15-atom cell (110 Ry).
- The fit rms is 0.3–2 meV/cell, on a range of 600–1100 meV.

| | J1 (meV) | J2C/J1 (C-bound pair, via M_N) | J2N/J1 (N-bound pair, via M_C) | J3/J1 | J4/J1 (tied) | t = kT_MC/(J1S²) |
|---|---|---|---|---|---|---|
| KVCr | 11.28 | 0.019 (Cr–Cr) | 0.036 (V–V) | −0.000 | pending | 1.326 |
| KVMo | 20.43 | 0.046 (Mo–Mo) | 0.028 (V–V) | +0.002 | pending | 1.300 |
| CrCr | 6.79 | 0.010 | 0.008 | −0.000 | ≈ 0.007 | 1.403 |

- All J2 values are antiferromagnetic, which frustrates the LCM.
- A J2 of 0.05 J1 on both sublattices lowers the classical T_C by 14 %. MF predicts 10 % and RPA 15 %.
- J3 = 0.05 raises T_C by 10.5 %. J4 = 0.05 lowers it by 9 %.
- The NNN correction lowers T_C by 8.1 % (KVCr), 9.9 % (KVMo) and 2.7 % (CrCr). So it moves the VCr/CrCr and VMo/CrCr ratios down by 5.5 % and 7.4 %, respectively.

## 4. Calibrants: checked compositions

| calibrant | compound actually measured | T_exp | model composition |
|---|---|---|---|
| CrCr | **Cr(III)[Cr(III)(CN)6]0.996(7)**, NMF route, Fm-3m, a = 10.42 Å (Schart et al., Inorg. Chem. 2024; PMC11615938) | 240 ± 10 K. Neutron (111) magnetic intensity vanishes at 240 K; χT minimum at 265 K | p = 0.996 ± 0.007 → D = 0.998 |
| KVCr | **K V(II)[Cr(III)(CN)6]·2H2O·0.1 KO3SCF3**, sol–gel, crystalline, a = 10.55 Å (Holmes & Girolami, JACS 121, 5593 (1999); review Table 9.8) | 376 K (365 K after heating) | p = 0.979 ± 0.021 → D = 0.987 |

**CrCr: what is not the calibrant**
- The classic Cr(II)[Cr(III)(CN)6]2/3·10/3H2O (Mallah 1993) also orders at 240 K. It is high-spin Cr(II) (S = 2) with 1/3 vacancies, so it is not the d3/d3 calibrant.
- Schart's own numbers: HSE06 gives J/k_B = −100 K (convention H = −J Σ_<ij>), classical MC (S²) gives 323 K, and Θ = −836 K.

**KVCr: where the composition bracket comes from**
- z = 1 by analysis, and M_sat = 0.7 kG cm³ mol⁻¹ = 0.125 μB = |w − 3v|.
- In the all-vacancy limit this gives v ≤ 0.042, hence the p bracket. The alternative, about 12 % V(III), changes T_C by less.

**Results**
- Calibration factors at base U: f_MC = 0.967 (CrCr) and 0.975 (KVCr); f_Q ≈ 0.60; f_RPA ≈ 0.64; f_MF ≈ 0.41.
- So PBE+U (U = 3) energies with classical S² spins reproduce both T_C to 3 %. A quantum S = 3/2 treatment of the same J would overshoot by ×1.6. In other words, the Ising-mapped PBE+U J are about 1.65× too large for a quantum Heisenberg model.
- Note that f depends on U: f_MC(KVCr) runs from 0.75 to 1.22. Only ratios between compounds at a matched U are used.

## 5. Vacancy dilution ([M_C(CN)6] vacancies, C sublattice only; N-bound metal water-capped)

**Model T_C(p)** (classical MC, J1 only; `mc/dil_p*.json`):

| p | 1 | 0.95 | 0.90 | 0.85 | 0.80 | 0.73 | 2/3 | 0.60 | 0.50 |
|---|---|---|---|---|---|---|---|---|---|
| D(p) = T_C(p)/T_C(1) | 1 | 0.969 | 0.936 | 0.902 | 0.870 | 0.817 | 0.766 | 0.713 | 0.616 |

- The percolation threshold of the C-diluted rock-salt net is **p_c ≈ 0.135**. It equals fcc site percolation with 1st + 2nd neighbours, because N sites are linked only through C sites. It is not the sc site value 0.312.
- With the KVMo NNN set, D(0.73) = 0.811, against 0.817 for J1 only.

**KV[Mo(CN)6] versus vacancy fraction v = 1 − p.** In the table, M_uncomp = 3v μB/f.u. (Luttinger zero lost).
- A(p) = A·D(p).
- B(p) = 417 K rescaled from p = 0.73, using the mean of the MC and the empirical V–Cr dilution laws.
- C(p) = the Birge-weighted combination of A(p) and B(p).

| v | 0 | 0.03 | 0.05 | 0.10 | 0.15 | 0.20 | 0.27 (Magott) | 1/3 |
|---|---|---|---|---|---|---|---|---|
| A (K) | 673 ± 77 | 660 | 652 | 630 | 607 | 586 | 550 | 515 |
| B (K) | 487 ± 49 | 479 | 474 | 462 | 449 | 437 | 417 | 399 |
| C (K) | 535 ± 82 | 526 | 519 | 503 | 486 | 471 | 448 | 427 |
| M_uncomp (μB/f.u.) | 0 | 0.09 | 0.15 | 0.30 | 0.45 | 0.60 | 0.81 | 1.00 |

## 6. Consistency checks (vacancy-rich V–Cr; V–Mo parent)

**Reference point.** The stoichiometric V–Cr model, 376 K / D_sample = 381 K.

| sample | p, w | model | exp | exp/model |
|---|---|---|---|---|
| Cs0.82V[Cr(CN)6]0.92–0.94 (crystalline, Holmes–Girolami) | 0.93 | 364 K | 337 K | 0.93 |
| V(II)[Cr(CN)6]2/3 (crystalline, same paper) | 0.667 | 292 K | 330 K | 1.13 |
| **V(II)0.42V(III)0.58[Cr(CN)6]0.86·2.8H2O** (Ferlay 1995, amorphous) | 0.86, 58 % V(III) S = 1, cluster-DFT J(V(III)–Cr)/J(V(II)–Cr) = 0.41 | 207 K | 315 K | 1.52 |
| same, with V(III)–Cr coupling equal to V(II)–Cr | | 280 K | 315 K | 1.12 |
| **V(II)1.37[Mo(CN)6]** (Magott 2025, amorphous), with A | 0.73 | 550 K | 417 (340–454) K | **0.76** |

**V–Cr samples**
- The crystalline V–Cr samples scatter by ±13 % around the dilution model, with no systematic sign.
- **Ferlay 315 K is not reproduced with the assigned composition.** Even at a V(III) coupling equal to V(II), the model is 11 % low. With the cluster-DFT ratio it is 34 % low.
- Ferlay is therefore not used for calibration. It says the V(III)–Cr coupling cannot be much weaker than V(II)–Cr, or else the V(III) fraction is overstated.

**V–Mo parent: a 24 % over-prediction**
- The parent's T_C is an extrapolation: M(T) was measured only to 340 K, and the sample decomposes above that. The Bloch exponent 1.93 gives 417 K and the Bloch exponent 3/2 gives 454 K.
- Even against 454 K, the DFT-ratio model over-predicts by ×1.21.
- So either PBE+U overestimates J(V–Mo)/J(V–Cr), or the amorphous, MeCN-capped parent has a suppressed T_C.
  - The experiment-implied ratio is 1.21–1.35, against 1.77 in the model.
  - The amorphous Ferlay V–Cr sample sits 10–15 % below crystalline V–Cr, which can explain only part of the gap.
- The pending HSE06 KVMo ratio is the discriminating number.

## 7. Finite-temperature compensation (classical MC, L = 16, `mc/mt_*.json`)

With J2C ≠ J2N the two sublattices are no longer equivalent at T > 0. M_N(T) − M_C(T) then gives a net moment:

| | net moment (μB/f.u.) | grows from | peak | sign |
|---|---|---|---|---|
| KVCr | ≤ 0.018 | 0.003 at 0.1 T_C | 0.018 at 0.7–0.85 T_C | along Cr; the more frustrated V–V (J2N) sublattice disorders faster |
| KVMo | ≤ 0.020 | | | along V; Mo–Mo is more frustrated |

These are small next to the composition channel (3v) and the g-shift/SOC residue (≤ 0.05).

## 8. Caveats

- **Ratio transfer.** Prediction A assumes that f (or the DFT error in J) transfers from V–Cr/Cr–Cr to V–Mo. The parent check (§6) suggests it does not fully transfer.
- **NNN couplings.** They come from PBE+U at base U only, and the same ratios are applied at every U point and in HSE. J4 is resolved only via type-II alt; a tied J4C = J4N is used when flipC_II is missing.
- **S² mapping and classical MC.** The near-unity f_MC is an empirical fact, not a derivation. Only ratios are used.
- **Disorder model.** Vacancies are random and uncorrelated. Real PBAs have correlated vacancy networks (Simonov et al., Nature 2020), and the J of the water-capped metals is assumed unchanged.
- **Hydration.** The calibrant KVCr is hydrated (zeolitic 2 H2O). `tc_model.py` adds the hydrated-cell variant (KVCr·2H2O, U 3/3) once the parent job has finished.

## 9. Files

| file | content |
|---|---|
| `tc_model.py` | the model, calibration, predictions, error budget → `summary.json` |
| `sc_mc.py` | numba classical Heisenberg MC kernel |
| `tc_mc_runs.py` | Binder-crossing T_C driver. Modes: `dilution`, `nnn`, `model`, `fromsummary`, `mt`, `perc` |
| `jexch_submit.py` | supercell builder and Modal submit (`manifest_tc.json`) |
| `mc/` | MC tables |
| `results/` | cached supercell energies and preliminary HSE parse |
| `collect_snapshot.json` | collect.py output used |
| `inputs/` | supercell candidate JSONs |

**Compute used.** Three run_cpu64 jobs, `lcm/pba/tc/J2_{KVCr,KVMo,CrCr}`, at 10 SCFs each (about 1.5–2 h each). Local MC took about 20 core-minutes per T_C point. No `mrun cancel` was issued, and no other lane's files were touched.
