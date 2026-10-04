# Rock-salt-ordered YBaMnFeO5: an ideal d⁵/d⁵ Luttinger-compensated semiconductor whose cation order is predicted to be unreachable

**Track L package. Status 2026-10-02 11:30: DOWNGRADED. Not a discovery claim.** Session 66190d, campaign `magdisc`. All numbers come from Quantum ESPRESSO 7.5 runs on Modal unless a section says otherwise. The chronology is in `tracks/lcm/TRACKL_LOG.md`; the falsification evidence is in `tracks/lcm/lcm_ord2/`, `tracks/lcm/lcm_poly/` and `literature/lcm_scoop_check_2026-10-02.md`.

---

## V. Verdict after falsification (read this first)

Three independent adversarial referees (synthesizability, protocol and impact lenses; workflow wf_7a4c2237-88d) each concluded **"stopping condition not met; not big news"**. The decisive evidence comes from this lane's own calculations:

1. **The cation order cannot be made.**
   - The paramagnetic cluster expansion (21 YBaMnFeO5 arrangements, all relaxed; `lcm_ord2/ord2_summary.json`) gives **T_order ≈ 950 K (+250/−150 K)**. This is below the 1173–1573 K window in which Mn and Fe can actually exchange sites.
   - At 900–1300 °C there is only short-range order: unlike-neighbour fractions are 0.56–0.66 (0.5 is random).
   - The equilibrium antisite fraction is 13 % even at 873 K.
   - E(random SQS) − E(rock-salt) = 116 ± 15 meV/f.u. (paramagnetic). ΔE/(2k_B ln 2) ≈ 970 K, an independent cross-check.
   - The benchmark YBaCuFeO5, computed with the same protocol, has T_order = 1120 K. Experimentally it reaches only η = 0.10–0.16 after 50 h at 1150 °C and slow cooling (Morin et al., Nat. Commun. 7, 13758 (2016)). So the protocol does not overestimate disorder, and YBaMnFeO5 would freeze in near η ≈ 0.
   - Annealing below T_order is kinetically frozen (B-site hopping is 10³–10⁷ times slower than at the freeze-out temperature).
   - (111) layer-by-layer epitaxy cannot set Y/Ba [001] layering and Mn/Fe rock-salt order at the same time.
   - **Root cause (design principle):** exact d⁵/d⁵ compensation with 3d ions forces Mn²⁺/Fe³⁺. That charge difference, Δz = 1, is too small to drive rock-salt B-site order; Δz ≥ 2 is normally needed. Every mixed-B 112-O5 compound ever made (Gd-, Nd-, Sm-BaMnFeO5, YBaMnCoO5) is B-disordered.
2. **Disorder destroys the headline property.** Compensation survives disorder (all ions d⁵), but the s-wave splitting and unipolar edges need the order.
   - A nearest-neighbour Mn/Fe antisite pair charge-transfers to Mn³⁺/Fe²⁺ and closes the opposite-spin gap to 0.01 eV.
   - Rock-salt antiphase domains reverse the sign of the splitting.
3. **Corrections to earlier drafts (errors found by the referees):**
   - The HSE06 numbers (gap 2.32 eV, windows 1.10/1.25 eV, c_udud +51.7 meV/mag) were computed on the superseded **75-Ry** P4 cell (a = 5.618 Å), not the 110-Ry cell. The 110-Ry HSE runs (`lcm/deep/stack_hse_{G,Yflip}`) are pending.
   - "Every disordered arrangement is ≥ 0.10 eV/f.u. higher" is retracted. In the G state, RSapY (like pairs across Y) is only +9 meV/f.u. and in-plane stripes are +22. In the paramagnetic state they are +18 and +25–33. Both competitors are altermagnets, not LCMs.
   - "A far antisite pair costs about 1.9 eV" is retracted: that number was the Y/Ba antisite. A Mn/Fe pair costs 0.22 eV (paramagnetic) to 0.72–0.84 eV (G state).
   - The low-temperature-anneal and (111)-epitaxy synthesis routes are retracted.
   - P4/ncc rotations cost +6.7 meV/f.u., not +40. The +40 is the energy of the unrotated P4/nmm structure. A tilted P2₁/c minimum sits at +3.7 meV/f.u.
4. **What survives (for the ideal ordered cell only):**
   - G-type d⁵/d⁵ LCM ground state at every U from 0 to 6 eV.
   - Gap ≥ 1.4 eV and same-spin band edges.
   - T_N 417 K (raw MC) to about 490 K (calibrated).
   - Dynamically stable at Γ and Z.
   - No lower polytype: tilt variants are +3.7 meV/f.u. and up, brownmillerites +290 and up.
   - Hull: +2.6 meV/atom at 0 K, and the redox channel is closed (+36 meV/atom).
   - The twin YBaMn2O5 is described correctly by HSE06 at the experimental lattice: G-FiM is 181 meV/cell below A-type (at the PBE+U lattice the two are degenerate, −1 meV/cell). So HSE is reliable for this family, but the result is strongly geometry-sensitive.
   - Valid as a PRB-level design study of the ideal structure. **Not a realizable discovery.**
5. **Lesson carried into the next search:** the sublattice inequivalence that makes a magnet an LCM must be enforced by robust crystal chemistry, not by a weakly driven ordering. Examples: tetrahedral vs octahedral sites with strong site preferences, or Δz ≥ 2 ordering. Compensation should then be insensitive to swaps among the magnetic ions.

---

## 0. Original abstract (superseded; kept for the record — see §V for corrections)

Luttinger-compensated magnets (LCMs) are collinear ferrimagnets whose spin-up and spin-down sublattices are not related by any symmetry. In an insulator, Luttinger's theorem then pins the net moment to exactly zero, and the bands are still spin-split at Γ, as in a ferromagnet ("s-wave" splitting).

**No LCM semiconductor that orders above room temperature is known.** Experimental LCM insulators order at ≤ 48 K. The room-temperature compensated systems are metallic Heuslers (Mn2RuxGa), or 5d compounds whose spin–orbit coupling leaves a net moment (Sr2CrOsO6, about 0.7 μB).

We predict that **YBaMnFeO5 with rock-salt (3D checkerboard) order of Mn²⁺ and Fe³⁺ on the square-pyramidal B sites** of the oxygen-deficient "112" layered perovskite (P4/n) is such a material:

- **Magnetic order.** All six nearest-neighbour Mn–O–Fe superexchange paths are antiferromagnetic.
  - The ground state is G-type, Mn↑/Fe↓. This holds for PBE+U at U = 0–6 eV and for HSE06 (the competing c_udud order lies +52 meV/magnetic ion higher).
  - It was also reached by simulated annealing of a 10-parameter exchange model.
- **Compensation.** Both ions are d⁵ (S = 5/2, L = 0), so the moment cancels exactly.
  - The cancellation survives spin–orbit coupling (≤ 0.01 μB/cell) and Mn/Fe antisites (both d⁵).
  - At 300 K the thermal residue is about 0.05 μB/f.u.
- **Band gap.** 2.32 eV (HSE06; 1.3–1.4 eV in PBE+U).
- **Spin-polarised edges.** Both band edges lie in the same spin channel. Fully polarised windows are 1.10 eV (valence band) and 1.25 eV (conduction band) in HSE06.
- **Conduction electrons.** The conduction-band minimum is dispersive along c (m* ≈ 0.5 m₀). Its spin polarisation survives worst-case spin-flip disorder: a window ≥ 0.15 eV remains with 12.5 % of spins reversed.
- **Néel temperature.** Classical Monte Carlo on a 10-class DFT exchange model gives T_N = 417 ± 10 K. Calibration against YFeO3, computed with the same protocol, gives **≈ 490 K**. The range across calibrants is 490–640 K; the U-dependence is about −20 % at U = 6 eV.
- **Thermodynamics.** At 110 Ry against 25 competing phases, the ordered compound lies 2.6 meV/atom above the convex hull, which is on the hull within DFT accuracy.
  - Every B-site-disordered arrangement tested is higher in energy: columnar +0.10, layered +0.37 eV/f.u. (more in §5).

**The central risk is chemical order.** Every reported LnBaMnFeO5 (Ln = Gd, Nd) is B-site disordered. We give a neutron-diffraction test that separates chemical order from magnetic order (§9).

---

## 1. Material

| Item | Value |
|---|---|
| Formula | YBaMnFeO5 = Y³⁺Ba²⁺Mn²⁺Fe³⁺O²⁻₅ (ortho-atomic d occupations: Mn n_d 5.26, m 4.5 μB; Fe n_d 5.92, m 4.0 μB) |
| Structure type | Oxygen-deficient "112" double perovskite (YBaMn2O5 / YBaFe2O5 / YBaCuFeO5 family). Y/Ba A-site layering; O vacancies in the Y layer; square-pyramidal B sites |
| B-site order | Rock-salt: every Mn has only Fe as nearest magnetic neighbours (4 in-plane at 153°, 1 apical at 180° through the BaO layer, 1 across the O-free Y layer), and vice versa |
| Space group | **P4/n (No. 85)**, origin choice as in the CIF. At 110 Ry the P4 and P4/n relaxations converge to the same structure (ΔE = 0.05 meV/f.u.) |
| Cell (PBE+U, U = 4 eV, 110 Ry, G-type) | a = 5.665 Å, c = 7.689 Å (√2a_p × √2a_p × 2a_p), Z = 2, V = 246.74 Å³, ρ = 5.61 g cm⁻³ |
| Wyckoff (P4/n) | Ba 2a (0, 0, 0); Y 2b (0, 0, ½); Mn 2c (0, ½, 0.2609); Fe 2c (½, 0, 0.2583); O1 2c (½, 0, 0.0088); O2 8g (0.2144, 0.2287, 0.6858) |
| Bonds | Mn–O 2.074 (apical) / 2.113 Å (×4). Fe–O 1.918 / 2.005 Å (×4). Mn–O–Fe 153.1° in-plane, 180° apical |
| Files | `structures/YBaMnFeO5_rocksalt_P4n_PBEU110.cif`, `structures/POSCAR_YBaMnFeO5_rocksalt_P4n_PBEU110`, `structures/YBaMnFeO5_rocksalt_P4n_LCM_G_P1.mcif` (G-type; moments drawn along c, but the true easy axis is in-plane). Job-input JSON: `inputs/YBaMnFeO5_P4n_110relaxed.json`. The 75-Ry files are in `structures/superseded_75Ry/` |

## 2. Novelty (adversarial audit: `tracks/lcm/novelty_audit_YBaMnFeO5.json`; literature: `literature/lcm_prior_art.md`)

| Check | Result |
|---|---|
| Exact composition in databases | Alexandria agm011206190 and OQMD 1740299. Both use **layered** Mn/Fe order and were computed **FM only** (metallic in those calculations) |
| Rock-salt Y member | Absent from MP, Alexandria, GNoME, OQMD, AFLOW, JARVIS and COD |
| Rock-salt LnBaMnFeO5 / LnSrMnFeO5 | Enumerated in Alexandria for Ln = Pr, Nd, Ho (Ba) and Pr, Sm (Sr), FM only. No magnetic-ground-state, compensation or spin-splitting analysis |
| Synthesised relatives | GdBaMnFeO5 (Manabe 2016; García-Martín 2017), NdBaMnFeO5+δ (Hossain 2024), LnBaMnFeO5 cathode studies (Muñoz-Gil 2015; Gilev 2019). Mn²⁺/Fe³⁺ confirmed, but **B sites random**. YBaMn2O5 doped with ⁵⁷Fe (Rykov 2009) |
| LCM / compensated-FiM / spin-splitting discussion for any LnBaMnFeO5 | None found |

**Defensible novelty statement.** This is the first study of the magnetic ground state and spin-split band structure of rock-salt-ordered LnBaMnFeO5. It includes the prediction that the non-magnetic-A member, YBaMnFeO5, is a Luttinger-compensated semiconductor with fully spin-polarised band edges that orders above room temperature. The composition and the rock-salt prototype themselves are not new; the physics and the target state are.

## 3. Key predicted properties

| Property | Value | Method / file |
|---|---|---|
| Magnetic ground state | G-type, Mn↑ Fe↓: compensated FiM = LCM (findspingroup SSG 75.75.1.1.L, "Zeeman"/s-wave). Not altermagnetic, because no operation maps the Mn and Fe sublattices onto each other | 8 k = 0 orders × U = 0, 2, 4, 6 eV; k_z = ½ stackings; 25 supercell configurations; annealing of the fitted model |
| Margin to the next order (110-Ry geometry, meV per magnetic ion) | Y-layer flip (k_z = ½): +16.4 / +11.7 / +8.9 / +7.0 at U = 0 / 2 / 4 / 6 (+8.2 at U_Mn 4, U_Fe 6). Ba-layer flip: +35 to +76. FM: +160 to +330 | PBE+U, 1×1×2 cells, 90 Ry |
| HSE06 arbitration | E(c_udud) − E(G) = **+51.7 meV/mag** | HSE06, nq 2×2×1, k 0.25 Å⁻¹ |
| Net moment | 0 μB (collinear insulator); ≤ 0.01 μB/cell with SOC; ≈ 0.05 μB/f.u. at 300 K (sublattice MC) | |
| Band gap | **2.32 eV (HSE06, nq 2×2×1)**. PBE+U (U = 4): 1.32 (75-Ry cell) / 1.40 eV (110-Ry cell); 0.57–1.67 eV for U = 0–6 | |
| Edge spin polarisation | VBM and CBM both spin-up (unipolar). Fully polarised windows **1.10 eV (VB) / 1.25 eV (CB) in HSE06**; 0.45 / 1.43 eV in PBE+U (110 Ry) | global spin windows on a dense k grid plus the high-symmetry path |
| Band edges | VBM: flat Mn-d/O-p states on Z–R. CBM at Z: Fe-d (minority of Fe = global spin-up), m* ≈ 0.53 m₀ along c, ≈ 1.3 m₀ in-plane | PBE+U bands, `figures/bands_PBEU.png` |
| Exchange (meV; E = −Σ J s·s with \|S\|² included; 110-Ry geometry, U = 4) | J_ip −32.6, J_ap −53.7, J_Y −12.8. Next-nearest neighbours −0.2 to −3.2, all AFM, so not frustrating for the G state | 25 configs, 10 classes; rms 2.7, LOO 5.3 meV/cell (`data/jfit110_*`) |
| T_N | **417 ± 10 K** (classical Binder MC, L = 6/8/10). Calibrated **≈ 490 K** (YFeO3, exp/MC = 1.17–1.19); 490–640 K across d⁵ calibrants (MnFe2O4 1.46–1.53, weaker because experimental samples are about 20 % inverted). RPA cross-check ≈ 477 K | `tracks/lcm/lcm_cal/CAL_LOG.md` |
| U sensitivity of T_N | J scales roughly as 1/U; T_N(U = 6) ≈ 0.75–0.8 × T_N(U = 4) → ≈ 370–390 K calibrated. ⏳ linear-response U (hp.x) | |
| Magnetic anisotropy | Easy plane, 0.34 meV/magnetic ion | PBE+U+SOC (FR PseudoDojo) |
| Magnetic space group | Easy-plane moments (ground state): **P2'/c' (BNS 13.69, type III)**; moments along c: P4/n (BNS 85.59, type I). Both are ferromagnetic-compatible: a net in-plane (or c) magnetisation is symmetry-allowed, so the anomalous Hall effect, MOKE and XMCD are allowed at zero net spin moment | spglib magnetic symmetry on the 110-Ry cell |
| E_hull (0 K) | **+2.6 meV/atom**, decomposing to MnO + Y2O3 + MnFe2O4 + Ba2YFeO5. 26 competing phases at 110 Ry, each vc-relaxed with its own magnetic-order scan. Includes YBaFe2O5: the redox channel 2 YBaMnFeO5 → YBaMn2O5 + YBaFe2O5 is +36 meV/atom uphill. ⏳ Fe3O4, Mn2O3, Ba2Fe2O5, BaFe2O4, Ba6Y2Fe4O15 | `data/hull_QE110.json` |
| Phonons | 110-Ry P4/n, 1×1×2 supercell: **all real at the exact Γ and Z = (0, 0, ½)** (lowest optical mode at Z 0.82 THz), so the q_z ≠ 0 imaginary branches of the earlier 75-Ry 2×2×1 run were interpolation artefacts. 75-Ry 2×2×1: real at Γ, X, Y, M. 110-Ry 2×2×1: all real (exact X/Y 1.62, M 1.68 THz; in-plane mesh min +0.36 THz) | `data/phonon110_summary.json`, `figures/phonons_PBEU*.png` |

## 4. Why it matters

- **A missing material class.** No bulk compensated semiconductor (LCM or altermagnet) ordering above room temperature has experimentally verified spin-split band edges. An LCM has two properties that no other class combines:
  - zero net magnetisation, so no stray field and antiferromagnet-like robustness and dynamics;
  - **Γ-centred, isotropic (s-wave) spin splitting of the band edges**, as in a ferromagnetic semiconductor.
  - Altermagnets cannot do this. Their Γ point is spin-degenerate and their splitting is d/g/i-wave, so a polycrystal averages it out.
- **A "compensated half-semiconductor".** In YBaMnFeO5, electrons in the conduction band are 100 % spin-polarised over a 1.25 eV window (HSE). The polarisation survives strong spin disorder. This adds:
  - spin-polarised n-type transport and photocarriers at zero magnetisation;
  - spin filtering and TMR with compensated electrodes;
  - magneto-optics and the anomalous Hall/Nernst effect at zero net moment (the magnetic point group is ferromagnetic-like).
- **Chemically locked carrier spin.** Extra electrons reduce Fe³⁺ to Fe²⁺. The electron enters the Fe minority channel, and since every Fe is on the spin-down sublattice, that channel is global spin-up. Holes oxidise Mn²⁺ to Mn³⁺ and leave the Mn majority channel, which is also global spin-up.
  - Both carrier types therefore have the same spin even if they self-trap as small polarons, which is common in ferrites and is the usual failure of band-picture spin-polarisation claims in oxides.
  - Only mis-oriented spins or antisites (§5.5, §8.1) create opposite-spin carrier sites.
- **Protected compensation.** Previous near-room-temperature compensated insulators rely on 4d/5d ions. There, spin–orbit coupling leaves a sizeable net moment (Sr2CrOsO6, about 0.7 μB/f.u.). Here, L = 0 d⁵ ions keep the compensation exact. Mn/Fe antisites cannot unbalance it either, because both ions carry S = 5/2.
- **Chemistry.** Earth-abundant and non-toxic (Y, Ba, Mn, Fe, O). The parent family is already synthesised.

## 5. Computational evidence (all runs; job names in the log)

1. **Search.**
   - DB screen for exactly compensable collinear orders (MP 7,655 and Alexandria 13,616 hits).
   - A design panel independently proposed the rock-salt d⁵/d⁵ 112 design by two routes.
   - Stage-1 PBE+U multi-configuration screening of about 40 candidates.
   - This structure ranked first on (T_N proxy × gap × spin window × compensation robustness).
2. **Cation order (G-type, vc-relaxed).** Rock-salt 0; columnar +0.10; layered +0.37 eV/f.u.
   - ⏳ Full order–disorder thermodynamics (cluster expansion from 18 relaxed arrangements, plus a YBaCuFeO5 benchmark against its known partial disorder) in `tracks/lcm/lcm_ord2/`.
   - Non-rock-salt arrangements show partial Mn³⁺/Fe²⁺ charge transfer and smaller gaps (0.15–0.26 eV).
3. **Polytypes.** In `tracks/lcm/lcm_poly/`, all relative to the 110-Ry P4/n ground state:
   - inverted 112 (Ba/Y swapped): +2627 meV/f.u.;
   - octahedral-rotation variant P4/ncc: +40 meV/f.u.;
   - brownmillerite variants: +321 to +372 meV/f.u.;
   - antisite pairs: +479 meV/f.u. (far) and +1223 (stripe);
   - ⏳ tilt and stacking variants, MLIP screen of 20 Ba2Fe2O5-framework decorations.
4. **Magnetism.**
   - Ground-state search: U-scan, 110-Ry stackings, HSE arbitration.
   - Spin group (findspingroup and amcheck on a P1 mcif).
   - SOC, exchange fits (75-Ry 22-config and 110-Ry 25-config, direction-resolved pair classes).
   - Binder MC, RPA, annealing.
   - Sublattice-resolved finite-T magnetisation.
5. **Spin-flip robustness** (72-atom snapshots, 110-Ry geometry). Worst-case flips (flipFe, flipMn, 2 far flips = 12.5 %) keep the CBM spin-up with a window ≥ 0.15 eV and a gap > 0.85 eV. A flipped Mn puts occupied opposite-spin defect levels above the VBM. **The hole polarisation therefore needs good magnetic order; the electron polarisation does not.**
6. **Hull.** 110-Ry vc-relaxations for all phases. Magnetic order scanned for every competitor; the same PseudoDojo/U protocol throughout.
   - The 75-Ry estimates (−55 to −29 meV/atom) were Pulay artefacts and are retracted.
7. **Calibration of T.**
   - YFeO3 (same Fe³⁺ d⁵ pyramid/octahedron chemistry): exp/MC 1.17–1.19.
   - MnFe2O4: 1.46–1.53.
   - YBaMn2O5 (twin): invalid as a calibrant, because PBE+U with U ≥ 2 gives its wrong ground state (§8).

## 6. Methods (reproduction)
- **Code.** QE 7.5 (conda-forge). PBE with PseudoDojo v0.4 NC pseudopotentials, scalar-relativistic (FR for SOC).
- **Hubbard U.** Dudarev, on Mn/Fe 3d, ortho-atomic projectors. U = 4 eV production; scans over 0–6 eV, including U_Mn ≠ U_Fe.
- **Cutoffs.**
  - Geometries and hull: 110 Ry (Pulay stress converged: P4 and P4/n coincide; absolute energies changed by up to 66 meV/atom from 75 to 90 Ry for MnO, hence 110).
  - Magnetic energy differences: 90 Ry on 110-Ry geometries.
  - k-spacing 0.25–0.28 Å⁻¹; gaussian smearing 0.005 Ry.
- **HSE06.** No U, nq 2×2×1, ecutfock = 2 × ecutwfc, Gygi–Baldereschi with x_gamma_extrapolation; on the PBE+U geometry.
- **Exchange model.** Pairs classified by direction (in-plane / apical through Ba / across Y). Distance-only shells conflate the in-plane (3.967 Å) and apical (3.971 Å) pairs and are invalid. D-optimal configuration sets; LOO cross-validated model selection.
- **Monte Carlo.** Classical Heisenberg Metropolis plus over-relaxation (numba); L = 6–16; Binder crossings. RPA (Tyablikov) as cross-check.
- **Phonons.** phonopy finite displacements (0.01 Å), PBE+U forces at 110 Ry.
- **Scripts.**
  - `scripts/jobs/*`: Modal job scripts (lcm_stage1, lcm_deep, lcm_hull110, lcm_hull_scan110, lcm_phon, lcm_single).
  - `scripts/analysis/*`: screen, fits, MC, spin group, hull LP, phonon post-processing.
  - `scripts/lib_snapshot/*`: the qeutil/magtools/amscreen versions used.
  - `inputs/*`: exact JSON job inputs. `data/*`: results.

## 7. Comparison with known materials

| Material | Type | T_order (K) | Gap (eV) | Net moment | Band-edge spin polarisation | Status |
|---|---|---|---|---|---|---|
| **YBaMnFeO5 (rock-salt), this work** | LCM semiconductor | **417 raw / ≈ 490 calibrated** | **2.32 (HSE)** | **0 (Luttinger; SOC ≤ 0.01 μB)** | **both edges; 1.1 / 1.25 eV windows** | predicted |
| GdBaMnFeO5 / NdBaMnFeO5 (B-random) | AFM (+ 4f) | > 300 (exp.) | insulating | ≠ 0 (4f) | none (random sublattices) | exp. |
| Sr2CrOsO6 | near-compensated FiM | 725 (exp.) | ~0.3–0.8 | ~0.7 μB (Os SOC) | spin-split; strong 5d SOC | exp. |
| Mn2RuxGa | compensated FiM half-metal | ≤ 550 | metal | 0 at composition/T_comp | high (metal) | exp. |
| LaMn2SbO6 | LCM semiconductor | 48 | 1.57 | 0 | VB polarised | exp. |
| ε-Fe2O3 | compensated FiM (this work: LCM in DFT) | ~500 | 1.8 (PBE+U) | small | VB window 0.04 eV (this work) | metastable |
| GaFeO3 / AlFeO3 | LCM (ideal order) | 200–250 | ~2.3 | antisite-driven | 0.1–0.3 eV | exp. |
| MnTe | altermagnet | 307 | 1.3 | 0 | Γ spin-degenerate (g-wave) | exp. |
| CrVTiAl-type Heuslers | "FCF semiconductor" | > 400 | ≤ 0.2 with gapless channel | ~10⁻³ μB | — | exp. (disordered) |
| NiFe2O4 / YIG | ferrimagnetic insulators | 858 / 560 | 1.5–2.8 | 2 / 5 μB per f.u. | spin-split, large M | exp. |

## 8. Weaknesses and remaining uncertainty (scoped claims)

1. **B-site order is the main experimental risk.** Disorder leaves the compensation intact (d⁵/d⁵), but it destroys the inequivalence of the sublattices, and with it the s-wave splitting.
   - The ordered phase is favoured by ≥ 0.10 eV/f.u. at 0 K.
   - At synthesis temperature, the paramagnetic free energy of the ordered phase sits about +20 meV/atom above the hull. Configurational entropy of disorder (≤ kT ln 2 per B site ≈ 17 meV/atom at 1300 K) competes with this. The B-random analogues show that disorder wins under standard topotactic or high-T routes.
   - ⏳ The cluster-expansion T_order and the YBaCuFeO5 benchmark decide whether low-temperature annealing can order it.
   - **Antisites also damage the gap.** A nearest-neighbour Mn/Fe antisite pair (12.5 % of B sites, 72-atom cell) charge-transfers to Mn³⁺/Fe²⁺. M stays exactly 0, but the pair puts a donor–acceptor level in the opposite spin channel: that gap becomes 0.30 eV (SCF) / 0.01 eV (dense grid), while the host-channel gap stays 1.3–1.5 eV. Good chemical order is therefore needed for the semiconducting, unipolar edges, not only for the splitting. A single far antisite pair costs about 1.9 eV (110 Ry).
2. **DFT+U protocol caveat.** For the structural twin YBaMn2O5 (Mn²⁺/Mn³⁺, d⁵/d⁴), PBE+U with U ≥ 2 eV predicts the wrong ground state (A-type instead of the experimental G-FiM). We attribute this to Mn³⁺ e_g orbital physics, which is absent for d⁵–d⁵ bonds.
   - YBaMnFeO5 keeps G at every U from 0 to 6 eV, and HSE06 confirms it.
   - **HSE06 on the twin (PBE+U geometry): E(A) − E(G) = −1.0 meV/cell.** HSE removes most of the PBE+U error (−140 meV/cell at U = 4) but leaves the experimental order degenerate with A. ⏳ HSE at the experimental lattice.
   - ⏳ HSE of the closest non-LCM competitor of YBaMnFeO5 (Y-layer flip, +9 meV/mag in PBE+U).
3. **T_N carries DFT+U and calibration uncertainty.** The defensible range is 370–640 K; the central estimate is ≈ 490 K. The raw MC value without calibration is 417 K, above room temperature. At U = 6 the raw value drops to about 320–340 K.
4. **Hole polarisation is fragile.** It needs good magnetic order (spin-flip snapshots), so the claim is scoped to electrons, or to holes only well below T_N.
5. **Domains.** Time-reversed domains carry opposite s-wave splitting. Spin-ARPES and transport need a single domain: field-cooling through T_N with a small uncompensated moment, or exchange bias.
6. **Oxygen stoichiometry.** Each excess O (δ) oxidises Mn²⁺ and uncompensates about 2δ μB/f.u. Synthesis must hold δ ≲ 0.01; YBaMn2O5.00 shows this is achievable.
7. **Thermodynamics.** The compound is on the hull within DFT accuracy (+2.6 meV/atom), not comfortably below it. ⏳ Eight competitors are still running, including YBaFe2O5 for the redox channel 2 YBaMnFeO5 → YBaMn2O5 + YBaFe2O5.
8. ⏳ 110-Ry 2×2×1 phonons (in-plane zone boundary).
9. **Isolated vs paired antisites.** A far antisite pair keeps M = 0, a 1.0–1.2 eV gap and the CBM in the host spin channel (0.23 eV window); only the VBM turns into an opposite-spin defect level. NN pairs (charge transfer) are the damaging defect.

## 9. Proposed synthesis and decisive experiments

- **Synthesis.**
  - Ceramic route from Y2O3, BaCO3, MnO and Fe2O3 under controlled low pO2 (Ar/H2, or a sealed tube with a Ti or Zr getter) at 1100–1300 °C.
  - Then **long low-temperature annealing (≤ T_order)** to develop rock-salt order. Avoid topotactic reduction of B-disordered precursors.
  - Alternative: atomic-layer-controlled (111)-oriented epitaxy of the oxidised perovskite followed by topotactic reduction (precedent: rock-salt-ordered LaFeO3/LaCrO3 (111) superlattices, Ueda et al. *Science* 1998).
- **Detecting the order needs neutrons.**
  - Mn and Fe are indistinguishable by lab X-rays (Z = 25/26).
  - With neutrons (b_Mn = −3.73, b_Fe = +9.45 fm), the **(101) reflection of the √2×√2×2 cell** (2θ = 19.4° at λ = 1.54 Å) is the strongest line of the ordered pattern and vanishes for random B sites (`data/order_fingerprint_nuclear.json`).
  - G-type magnetic scattering falls on the same positions. Measure **above T_N** (or with polarisation analysis), because B-random GdBaMnFeO5 is also G-type.
  - Resonant X-ray diffraction at the Fe/Mn K edges is a lab-independent alternative.
- **Falsification tests.**
  1. Neutron diffraction: G-type with Mn and Fe moments of equal magnitude and opposite sign.
  2. M(T) ≈ 0 at low T, with no compensation point.
  3. **Spin-resolved photoemission or inverse photoemission on a single-domain sample: the CBM at Z must be single-spin over ≳ 1 eV.**
  - A spin-degenerate CBM in a chemically ordered, magnetically ordered sample would falsify the prediction.
