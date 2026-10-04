# PBA dynamic stability: phonons of KV[Cr(CN)6] and KV[Mo(CN)6] (PBE+U, LCM)

Track L (session 66190d), 2026-10-03. Gate 0 is in `../README.md`; defects are in `../defects/README.md`.

- Build: `python tracks/lcm/pba/phon/build_phon.py` writes `phonopy_disp_<m>_<tag>.yaml` and `inputs/<m>_<tag>.json`.
- Submit: `python tracks/lcm/pba/phon/submit_phon.py [--tags prim,conv,L30]`. It never cancels anything.
- Collect: `python tracks/lcm/pba/phon/collect.py --fetch`. It reads fixed paths from `manifest.json` (no listdir), caches them in `results/`, and rewrites `summary.json` and `results/phonons_full.json`. Run it again to pick up the jobs still running (§6).

## 0. Verdict

**The d3/d3 cyanide framework is dynamically stable. The ordered F-43m model with K at the centre of its cavity is not.**

- **Every framework mode is real at every exact q computed**, for both compounds. That covers Γ and X, and also L for KVCr; the KVMo L set is 20/23 done, §6.
  - The softest framework optic mode is the cooperative octahedral tilt (rigid rotation of the V(NC)6 and M(CN)6 units, with transverse cyanide motion):
    - Γ T1, the counter-rotating "a0a0c−" type: **+1.39 THz (46 cm⁻¹)** for KVCr, **+1.36 THz (45 cm⁻¹)** for KVMo;
    - X, the in-phase "a0a0c+" type: **+1.41 / +1.38 THz**.
  - Cyanide bending, libration, M–C / V–N stretching and C≡N stretching (2140–2150 cm⁻¹ for KVCr, 2093–2114 cm⁻¹ for KVMo) are all stiff.
- **The only imaginary branch is K⁺ rattling.** The modes are ≥ 87 % K by kinetic energy, and the branch is nearly flat:
  - Γ T2 (all K in phase, polar): **−1.29i THz** (KVCr), **−1.35i THz** (KVMo);
  - X (antiphase): **−1.04i / −1.14i THz**, doubly degenerate;
  - L (antiphase, KVCr): **−1.13i THz**, doubly degenerate.
  - The interpolated minimum from the 60-atom force constants is −1.27 / −1.34 THz, near Γ, on the same flat branch.
- **The K potential is a deep, multi-minimum well**, not a weak harmonic softening. This is the 15-atom cell with all K shifted in phase:
  - With the framework frozen, the minima lie about 0.9–1.1 Å from the 4c centre, along ⟨100⟩, towards the square windows into the empty 4d cavities. They sit **≈ 42 meV (KVCr) and ≈ 65 meV (KVMo) per K** below the centred site.
  - The −⟨111⟩ minimum (towards the N-bound V corner) is next, at 35 / 53 meV. The +⟨111⟩ direction (towards the C-bound M corner) is the weakest, at 17 / 28 meV.
  - Fixed-cell relaxations started at those minima are still running (§4). They can only deepen the wells. At 10:53 they stood at −46 (KVCr ⟨100⟩) and −64 meV (KVMo ⟨100⟩), with the K–N/C contacts shortened from 3.82–3.88 Å to about 3.2 Å.
- **For the LCM claim this is a spectator instability.**
  - The magnetism does not move with K: every displaced SCF and every K-scan point (K up to 1 Å off-centre) keeps M = 0.00 and |M| within 0.01 μB. The windows of the relaxed off-centred cells come with `koff2` (pending). The centred relax reproduces the Gate-0 SCF-grid windows: 2.02 / 1.14 eV for KVCr, 0.97 / 1.87 eV for KVMo.
  - Compensation is fixed by S = 3/2 on both ends of the cyanide, not by any symmetry that a K shift could break.
  - The real material is a disordered K / H2O cubic phase anyway: KV[Cr(CN)6]·2H2O, fcc, a = 10.55 Å. The fact sheet does not say how K is located inside the cage.
  - What the result does say: an ordered, centred-K F-43m crystal is not the DFT ground state. Expect static or dynamic K off-centring, polar if it orders (Γ T2 is the most unstable exact-q K mode). Hopping through the cube faces into the empty 4d cavities is possible but was not computed.
  - Any property that depends on the K arrangement (Raman/IR selection rules, pyroelectricity, the polar-phase band edges) must be computed in an off-centred cell. `results/koff*/` carries the edge windows of the relaxed off-centred cells once they finish.
- **Overall: pass for the framework; flag for the K sublattice.** No structural kill. One caveat must enter the lead's write-up: "cubic F-43m with K centred in its cavity" is a saddle point, about 40–65 meV per K, in PBE+U at the PBE+U lattice constant (+1.2 % on a, larger cage). The smaller experimental a should reduce it. Hydration and dispersion were not tested and could move it either way.

## 1. Cells and q-points

The unit cell is the stage-1 relaxed 15-atom F-43m primitive cell: PBE+U, 110 Ry, U_V 3 / U_Cr 3 / U_Mo 2, LCM vc-relax, taken from `jobs/lcm/pba/s1a/<m>/results/<m>.json` → `relaxed`. It equals the exact F-43m projection to 1e-15 Å, a = 10.680 Å (KVCr) / 10.858 Å (KVMo), with K on 4c.

| tag | supercell (phonopy S) | atoms | exact q | k-grid | displacements |
|---|---|---|---|---|---|
| `prim` | I | 15 | Γ | 5×5×5 Γ (the Gate-0 SCF grid) | 10 (5 ± pairs) + reference |
| `conv` | [[−1,1,1],[1,−1,1],[1,1,−1]] (conventional cube) | 60 | Γ, 3 X | 2×2×2 shifted (the defect-cell grid) | 10 + reference |
| `L30` | Sᵀ = [[1,−1,0],[0,1,−1],[1,1,0]] (index-2 sublattice in reduced basis, 7.7 / 7.7 / 13.3 Å) | 30 | Γ, L (½ ½ ½) | 3×3×2 Γ | 22 (11 ± pairs) + reference |
| `s222` | 2I | 120 | Γ, 3 X, 4 L | (2×2×2 shifted) | built only, not run |

- **conv + L30 = the 2×2×2-equivalent result.** Together they give exact frequencies at every commensurate q of the 120-atom 2×2×2 supercell (Γ, X, L). The cost is about 1/10 of running the 120-atom cell itself, which would be 22 SCFs of 120 atoms at about 1 h on 64 cores each.
- In a rock-salt double perovskite the parent-perovskite tilt instabilities fold to Γ (R-point, a⁻a⁻a⁻ / a0a0c−) and X (M-point, a0a0c+), so `conv` covers every simple tilt system.
- The L30 cell was first built in the rhombohedral basis a1+a2, a2+a3, a3+a1 (33.6°). That gave a 180³ FFT box and about 40 min per SCF, so those 16 jobs were terminated after about 1 min (`manifest.terminated`) and rebuilt in the reduced basis. It is the same lattice.

## 2. Protocol

- **Electronic settings.** PBE+U with the Gate-0 stage-1 protocol (`jobsrc/pba_phon.py`):
  - PseudoDojo SR NC, 110 Ry (ecutrho 440);
  - ortho-atomic U;
  - Gaussian smearing 0.005 Ry, mixing 0.3;
  - **conv_thr 2e-11 × nat Ry**;
  - LCM start: C-bound M +3 μB, N-bound V −3 μB.
- **Displacements.** phonopy 4.7 finite displacements of 0.02 Å, with `is_plusminus=True`, so every displacement has its mirror.
  - Each displaced cell is translated so the displaced atom sits at the origin. Its site operations then carry no fractional translation, and the atom order is unchanged.
  - The forces of the undisplaced reference are subtracted. They are at most 0.020 eV/Å on C/N for KVCr, a k-grid/convergence residue, and cancel in ± pairs anyway.
  - The force constants are then symmetrized (`symmetrize_force_constants`).
- **Jobs.** One job per displacement for conv, one or two per job for prim, two per job for L30, all on run_cpu32. NPOOL is the largest divisor of 32 that is ≤ 2 × the irreducible k count.
- **Mode character** (`collect.py`). Each eigenvector is taken as a real displacement pattern in the 120-atom 2×2×2 cell, which is commensurate with Γ, X and L. The pattern is decomposed into kinetic-energy fractions:
  - species weights;
  - K in-phase (polar) or antiphase;
  - rigid rotation of each V(N)6 / M(C)6 ligand shell (least-squares ω);
  - cyanide transverse centre-of-mass motion ("tilt-type bending");
  - libration, stretch, and longitudinal motion;
  - a bond-coherence index, which marks acoustic-like modes.
  - Γ irreps come from phonopy (F-43m).
- **Checks.**
  - **SCF and magnetism:** all SCFs converged, and every displaced cell kept M = 0.00. |M| is constant to ±0.01 μB (28.28–28.30 / 26.58 per 4 f.u.).
  - **± pairs:** asymmetry |F(+u) + F(−u)| / |F(+u) − F(−u)| is 1–5 %.
  - **Acoustic sum rule:** the raw rows sum to < 1e-4 eV/Å² (QE removes the net force).
  - **Γ from three cells:** prim (5³ k) and conv (2³ shifted k) agree to 0.015 THz (KVCr) / 0.027 THz (KVMo) over all 45 modes, so the 2×2×2-shifted grid is converged for the forces. L30 (3×3×2 k) repeats prim Γ to ≤ 0.005 THz for the lowest 9 modes.
  - **X-star:** all three X points of `conv` give identical spectra.

## 3. Frequencies at the exact q-points (THz; cm⁻¹ = THz × 33.36; i = imaginary)

### KV[Cr(CN)6]

| q | modes, lowest first (degeneracy) |
|---|---|
| Γ (prim = conv) | **−1.29i K T2 (polar K off-centring, K 87 %)** ×3; acoustic 0 ×3; **+1.39 T1 tilt** ×3; 5.59 T1 CN bend; 6.70 T2 CN bend; 7.11 T2, 7.89 T2 metal vs CN; 9.09 T1 CN libration; 11.68 T1 / 12.29 T2 / 12.86 T2 libration; 13.06 E, 13.27 A1, 15.08 T2 M–C / V–N stretch; 64.18 T2, 64.27 E, 64.43 A1 C≡N stretch |
| X (conv) | **−1.04i K antiphase** ×2; +0.83 K ×1; +1.29 TA-like ×2; **+1.41 in-phase tilt** ×1; 2.20; 3.58 tilt-type ×2; 4.60 tilt-type ×2; 4.72–6.74 CN bend; 8.65–15.18 libration / stretch; 64.1–64.5 C≡N |
| L (L30) | **−1.13i K antiphase** ×2; +0.59 K ×1; **+2.69** (TA-like / V(NC)6 rotation) ×2; 3.04 tilt-type; 3.29 CN bend ×2; 4.06–8.25 tilt / CN bend; 10.3–12.6 libration; 13.8–15.2 stretch; 64.15–64.39 C≡N. Γ from this cell repeats prim to ≤ 0.005 THz for the lowest 9 modes (K −1.291i, tilt +1.392; 0.07 THz at most over all 45) |

### KV[Mo(CN)6]

| q | modes, lowest first (degeneracy) |
|---|---|
| Γ (prim = conv) | **−1.35i K T2 (polar, K 88 %)** ×3; acoustic 0 ×3; **+1.36 T1 tilt** ×3; 5.43 T1 CN bend; 5.85 T2 metal vs CN; 6.52 T2 CN bend; 7.80 T2; 9.04 T1 libration; 11.86 libration ×6; 12.88 T2; 13.79 E, 13.99 T2, 14.16 A1 stretch; 62.76 T2, 62.96 E, 63.39 A1 C≡N stretch |
| X (conv) | **−1.14i K antiphase** ×2; +0.38 K ×1; +1.04 TA-like ×2; **+1.38 in-phase tilt** ×1; 1.79; 3.36 tilt-type ×2; 4.45; 4.57 tilt-type ×2; 4.73 CN bend ×2; 5.40; 6.17 … |
| L (L30) | 20/23 SCFs at 10:53; collect later (§6) |

**Mode character**
- The lowest Γ T1 tilt is 99 % transverse cyanide centre-of-mass motion. 100 % of it is rigid rotation of the ligand shells (V 61 %, M 38 % of the kinetic energy), with every V(NC)6 counter-rotating against every M(CN)6.
- The X tilt has the same local character, alternating along one cubic axis.
- The framework is soft, as PBAs are, but not unstable. The lowest framework optic modes are the tilts at 45–47 cm⁻¹. The X-point transverse-acoustic branch sits at 1.29 THz (KVCr) / 1.04 THz (KVMo), and the lowest non-K mode at L (KVCr) at 2.69 THz.

## 4. K off-centring: how deep is the well?

**Rigid scan.** `kscan`: 15-atom cell, all K moved in phase, framework frozen, same protocol. Values are E − E(centred) in meV per f.u. (= per K), with the K force along u in parentheses (eV/Å).

| | 0.40 Å | 0.55 Å | 0.70 Å | 0.85 Å | 1.00 Å | estimated minimum |
|---|---|---|---|---|---|---|
| KVCr ⟨100⟩ | −15.1 | −25.7 | −35.7 | −42.0 (+0.024) | −41.7 (−0.032) | **≈ −42.5 at ≈ 0.92 Å** |
| KVCr −⟨111⟩ (towards V) | −15.9 | −26.7 | −34.7 | −33.6 | −13.5 | ≈ −36 at ≈ 0.77 Å |
| KVCr +⟨111⟩ (towards Cr) | −12.5 | −17.5 | −15.1 | +3.1 | +49.2 | ≈ −18 at ≈ 0.58 Å |
| KVMo ⟨100⟩ | −17.6 | −30.9 | −44.8 | −57.0 | −64.5 (+0.028) | **≈ −66 at ≈ 1.1 Å** |
| KVMo −⟨111⟩ (towards V) | −18.5 | −32.2 | −45.4 | −52.8 | −46.3 | ≈ −53 at ≈ 0.87 Å |
| KVMo +⟨111⟩ (towards Mo) | −15.5 | −24.2 | −28.1 | −20.0 | +10.4 | ≈ −28 at ≈ 0.70 Å |

- **Harmonic check.** k_eff = m_K ω² from the Γ T2 mode gives −½ k u² ≈ −8 meV at 0.25 Å (k ≈ 0.27 eV/Å² for KVCr). The first SCFs of the relaxations started at 0.25 Å gave −5.6 to −7.5 meV there (`koff`, below), which is consistent once quartic terms are included.
- **Where K goes.** The ⟨100⟩ minimum moves K about 1 Å from the cage centre, 35–40 % of the way to the square V2M2(CN)4 window that opens into the empty 4d cavity. That shortens the four nearest K–N/C contacts from 3.82 to 3.24 Å (KVCr, 0.92 Å) and from 3.88 to 3.21 Å (KVMo, 1.1 Å).
- **Why the well is deeper in KVMo.** Plausibly because its cage is larger (a = 10.86 Å against 10.68 Å); not separated from the chemistry here.

**Relaxations** (`jobsrc/pba_koff.py`: fixed cell, ions, BFGS forc 2e-4 Ry/bohr, start symmetry kept).
- `koff/*_c` (centred, F-43m) is the energy reference. It sits 0.2 meV below the stage-1 positions (KVCr).
- `koff/*_{x,p111,m111}` started at 0.25 Å and crawled (K at 0.29–0.36 Å and −8 to −11 meV after about 10 BFGS steps). They were terminated and restarted from the scan minima as `koff2/*_{x (0.90 Å), m111 (0.80), p111 (0.55)}`.
- `koff/*_xy` (⟨110⟩, from 0.25 Å) is still running.
- `koff2` at 10:53, after 3–8 BFGS steps, relative to the centred relax, in meV per K:

  | | ⟨100⟩ (start 0.90 Å) | −⟨111⟩ (0.80 Å) | +⟨111⟩ (0.55 Å) |
  |---|---|---|---|
  | KVCr | −46.4 | −39.2 | −18.8 |
  | KVMo | −64.1 | −54.0 | −24.9 |

  The ordering ⟨100⟩ < −⟨111⟩ < +⟨111⟩ is unchanged from the rigid scan; framework relaxation adds a few meV so far.
- Final depths, K positions, the K–N/C contacts and the SCF-grid windows of the off-centred cells will appear in `summary.json` → `K_offcentring_relax_15atom` after `collect.py --fetch`.

## 5. Caveats

- **Harmonic, 0 K, fixed cell, PBE+U without dispersion, anhydrous, ordered K on 4c.**
  - All frequencies and K wells are computed at the PBE+U lattice constant, which is +1.2 % too large (KVCr: 10.68 Å against 10.55 Å measured for the dihydrate). A larger cage favours K off-centring.
  - Zeolitic water in the empty 4d cavities, as in the real KVCr·2H2O, would coordinate K and reshape its potential. Neither effect was computed.
- **No LO–TO term** (no Born charges / ε∞). The polar Γ T2 values are TO-like. Non-analytic corrections only raise LO branches and cannot remove an imaginary TO mode.
- **Collinear LCM only.** Spin–lattice coupling to other magnetic orders was not explored. The FM–LCM gap (≥ 152 meV/ion in PBE+U) makes it irrelevant for these modes.
- **Frequencies between the exact q, interpolated from the 60-atom conv force constants (8×8×8 mesh), are indicative only.** Their minimum (−1.27 / −1.34 THz near Γ) is the same flat K branch.

## 6. Jobs, status at 10:53, and how to collect the rest

All jobs are listed in `manifest.json` with call ids and result paths. Terminated jobs are under `manifest.terminated`; those were terminated with `terminate_containers=True` and are my own.

| block | jobs | status 10:53 |
|---|---|---|
| prim | `lcm/pba/phon/{KVCr,KVMo}_prim_g*` | done (11 + 11 SCFs) |
| conv | `lcm/pba/phon/{KVCr,KVMo}_conv_g*` | done (11 + 11) |
| L30 | `lcm/pba/phon/{KVCr,KVMo}_L30_g*` (12 jobs each) | KVCr done (23/23); KVMo 20/23, idx 18–20 running in `KVMo_L30_g7/g8/g9`, ETA about 11:10 |
| kscan | `lcm/pba/phon/kscan/{m}_{x,p111,m111}_g{0,1}` | done |
| koff | `lcm/pba/phon/koff/{m}_{c,xy}` | c done; xy running |
| koff2 | `lcm/pba/phon/koff2/{m}_{x,m111,p111}` | running (submitted 10:29; expect 0.5–1.5 h more) |

To finish: `source .venv/bin/activate && python tracks/lcm/pba/phon/collect.py --fetch`. This updates the L rows, the K relax table and the verdict string in `summary.json`.
- If the L30 sets are complete, `materials.<m>.exact["L[L30]"]` and `n_imag_exact` include L.
- `verdict` reads "UNSTABLE at an exact q" because it counts the K branch. `verdict_framework` and `n_imag_exact_framework` exclude modes labelled K rattling; that is the §0 statement.

## 7. Files

- `build_phon.py`, `submit_phon.py`, `collect.py`; `jobsrc/pba_phon.py` (forces), `jobsrc/pba_koff.py` (K relax).
- `phonopy_disp_<m>_<tag>.yaml`: phonopy displacement datasets.
- `inputs/<m>_kscan.json`: the rigid K-scan cells (idx = 10 × direction (1 ⟨100⟩, 2 +⟨111⟩, 3 −⟨111⟩) + amplitude index 0–4 for 0.40/0.55/0.70/0.85/1.00 Å). They were written at submission time and run with `pba_phon.py`; results are in `results/<m>_kscan/`.
- `results/phonopy_params_<m>_<tag>.yaml`: force constants.
- `results/phonons_full.json`: every mode at every exact q, with character.
- `summary.json`: compact.
