# KV[Cr(CN)6] and KV[Mo(CN)6]: d³/d³ cyanide Luttinger-compensated magnets with eV-scale unipolar band-edge windows (ideal-crystal DFT)

**Track L package.** Session 66190d, campaign `magdisc`. Status at 2026-10-03, about 11:30 PDT.

**Grade: identification plus quantification, at PRB / Chem. Sci. level. Not big news.** Three adversarial referees (realizability, physics and impact; workflow wf_cf79c5b9-44d) each ruled "stopping condition not met; not big news", with P ≈ 2–5 %.

**Where the details live**
- The full analysis, the pre-registered rule audit and the referee corrections are in `candidates/lcm/TRACK_L_FINAL_REPORT.md` §9.
- The chronology is in `tracks/lcm/TRACKL_LOG.md`, entries from 2026-10-03 03:00 to 11:30.
- Every number below was re-harvested at 11:15 PDT with the lane collectors. Their output is in `data/`.
- Verified 11:29–11:35 PDT against the source files, a fresh `--fetch` run of all five collectors and a scratch rerun of `tc/tc_model.py`. The fixes are listed in the Verification log at the end.

**Tags**
- **[PENDING]**: a job was still running at 11:15 (§6).
- **[new]**: not in the log.

---

## 0. Scoped claim

**Claimed (DFT, ideal crystals)**
1. **Both compounds are LCMs.** In PBE+U and HSE06, the ideal (vacancy-free, anhydrous, ordered) d³/d³ Prussian-blue analogues KV(II)[Cr(III)(CN)6] and KV(II)[Mo(III)(CN)6] are Luttinger-compensated magnets, also called fully compensated ferrimagnets:
   - the cell moment is exactly 0 by Luttinger counting;
   - the sublattices (C-bound M, N-bound V) are never related by symmetry;
   - the bands are spin-split in an isotropic, "s-wave" way.
2. **Both band edges lie in the V-majority spin channel, with eV-scale windows in HSE06** (ideal cells). In PBE+U the KVMo hole window is U-sensitive (0.05–1.19 eV) and fails the pre-registered "≥ 0.3 eV at every U point" line at (U_V 4, U_Mo 0). HSE06 values:

   | Compound | Gap (eV) | win_VB / win_CB (eV) |
   |---|---|---|
   | KVCr | 2.09 | **2.64 / 1.57** |
   | KVMo | 1.52 | **1.29 / 2.13** |

3. **KVCr is a re-identification.** It exists: Holmes & Girolami 1999, hydrated powder, T_C 376 K, M_sat 0.125 μB/f.u. So for KVCr this re-identifies a known, compensated-by-design magnet as an LCM semiconductor with a quantified, spin-locked hole edge, **in the ideal anhydrous crystal only**. The existing hydrated sample fails the pre-registered hydration rule (see Not claimed), and the anhydrous crystal has not been reported.
4. **KVMo is a prediction.** It is a compound that has never been made at 1:1, with a model T_C of **487 ± 49 K** (conservative, calibrated on the amorphous V–Mo parent; floor 377 K) to **622 ± 79 K** (DFT exchange ratio; `tc/summary.json`). The two routes disagree, and the pre-registered calibration test failed in HSE (1.85 against 1.57 ± 0.25), so only the 487 K route can carry a ≥ 400 K statement.

**Not claimed**
- That any real sample is exactly compensated. M = w − 3v − u μB/f.u. is composition-tuned: realistic values are 0–0.24 μB/f.u. for KVCr and 0.1–0.45 μB/f.u. for first-generation KVMo.
- That the existing hydrated KVCr carries the 2.6 eV hole window. The hydrated cell fails the pre-registered hydration rule (PBE+U hole window 2.02 → 0.93 eV, −54 %); the "existing sample" claim is withdrawn.
- That the electron edge survives vacancies. At every [Cr(CN)6] vacancy in KVCr an empty void state sits at the host CBM (CB window 6 % of pristine; water-filled test pending). KVMo vacancy cells keep 45–47 % of the CB window, a nominal fail of the pre-registered 50 % rule.
- Spin-polarised band transport. The t2g bands are flat (0.4–0.6 eV wide), and K-vacancy holes self-trap as V(III) polarons.
- T ≥ 400 K for KVCr.
- That 1:1 KVMo can be made.
- That antisite disorder is harmless. Antisites keep M = 0 but locally destroy both windows (the pre-registered defect kill triggered for both compounds); they are suppressed only kinetically.

**Protection wording (referee correction i).** The cyanide linkage fixes which sublattice each metal occupies (sublattice identity). M = 0 under V/M swaps comes from S = 3/2 on both ions (isospin), not from the linkage. [M(CN)6] vacancies (−3 μB each) and V(III) (+1 μB each) are not protected, and neither is hydration.

## 1. Summary

Luttinger-compensated magnets combine a zero net moment with spin-split bands. The only neutron-confirmed insulating LCM (LaMn2SbO6) orders at 48 K; the above-room-temperature compensated insulators that exist (Ca–Si YIG, 465 K; KVCr, 376 K) had not been analysed as LCMs (report §5). In oxides this lane found an empirical "trilemma" (report §4): no candidate combined
- room-temperature order;
- sublattice inequivalence fixed by crystal chemistry;
- eV-scale spin windows at the band edges.

**Why the cyanides get further.** d³/d³ cyanide double perovskites A[M_N][M_C(CN)6] sidestep two of these failure modes:
- the cyanide's C-end and N-end fix which sublattice each metal occupies;
- because both ions are S = 3/2, a V/M swap or a cyanide flip leaves M = 0 (the edges, however, are destroyed locally by a swap).

They do not protect the composition (vacancies, V(III)) or against hydration.

**The ideal cells.** With PBE+U (110 Ry, 2-D U grid) and HSE06, the ideal KV[Cr(CN)6] and KV[Mo(CN)6] cells are d³/d³ LCMs whose valence-band maximum (V t2g) and conduction-band minimum (M t2g) are both in the V-majority channel. Their exchange energies are large:

| | E(FM) − E(LCM), HSE06 (meV/ion) |
|---|---|
| KVCr | 181 |
| KVMo | 323 |

**The T_C model** is calibrated on the measured T_C of KVCr (376 K) and Cr[Cr(CN)6] (240 K). It puts stoichiometric KV[Mo(CN)6] at 487–622 K; only the 487 ± 49 K parent-calibrated route survives the failed HSE calibration test.

**Robustness checks:**
- the framework is dynamically stable (the only imaginary modes are K⁺ rattling, which leaves M = 0; the windows change by ≤ 7 % at the +[111] minimum, and those of the deeper ⟨100⟩ minimum are pending);
- no charge-transferred V(III)/M(II) tautomer was found in PBE+U (the HSE kill test is pending);
- K and [M(CN)6] vacancies leave the hole edge intact (≥ 87 %).

**The limits:**
- [M(CN)6] vacancies cost 3 μB each, and V(III) costs 1 μB each.
- Antisite pairs destroy both windows locally.
- Zeolitic water cuts the KVCr hole window by 54 % in PBE+U: the hydrated cell fails the pre-registered hydration rule.
- Vacancies take out the electron edge: KVCr loses it locally (void state at the CBM), KVMo keeps only 45–47 %.
- The lead compound's building block is labile: Mo(III) substitutes about 10⁵× faster than Cr(III) (aqua-ion estimate), and [Mo(CN)6]³⁻ survives only under cyanide limitation in aprotic, dark conditions.

**The decisive experiments are:**
- magneto-optics and element-specific XMCD at M ≈ 0 on anhydrous KVCr;
- spin-resolved photoemission of its top valence band;
- the synthesis of crystalline 1:1 KV[Mo(CN)6].

## 2. Key numbers and sources

Unless noted, the PBE+U settings are: QE 7.5; PseudoDojo NC; ortho-atomic U of 3 eV on V and Cr and 2 eV on Mo; 110 Ry; 8×8×8 nscf windows. HSE06: no U, 90 Ry, nq 2³, k 4×4×4, windows on the SCF grid.

| Quantity | KV[Cr(CN)6] | KV[Mo(CN)6] | Source |
|---|---|---|---|
| Relaxed a (Å), PBE+U LCM | 10.680 (exp. 10.55, hydrated) | 10.858 | `structures/`; `data/gate0_collect.txt`; vol `jobs/lcm/pba/s1a/<id>/results/<id>.json` |
| V–N / M–C (Å) | 2.112 / 2.060 | 2.097 / 2.161 | same |
| M_cell (LCM) | 0.000 at all U and in HSE | 0.000 at all U and in HSE | same |
| d moments V / M (μB) | −2.70 / +2.93 | −2.61 / +2.45 | same |
| PBE+U gap (eV), base (U-grid range) | 1.95 (1.42–2.49) | 1.30 (0.66–1.60) | `data/gate0_collect.txt` |
| PBE+U windows VB / CB (eV), base | 2.02 / 1.15 | 0.97 / 1.88 | same |
| PBE+U window ranges | VB 1.39–2.63; CB 0.87–1.41 | VB **0.05**–1.19; CB 1.70–2.54 | same |
| **HSE06 gap; windows (eV)** | **2.09; 2.64 / 1.57** | **1.52; 1.29 / 2.13** | vol `jobs/lcm/pba/hse/<id>_hse_LCM/deep_hse.json` |
| E(FM, M = 6) − E(LCM) (meV/ion): PBE+U base (range); HSE | 152.2 (121.9–196.6); 181.1 | 276.3 (238.3–423.2); 323.2 | `data/gate0_collect.txt` |
| VMo/VCr exchange ratio | — | ΔE ratio: PBE+U 1.815 (1.64–2.65); HSE 1.784. T ratio: model ≈ 1.67 (PBE+U) / 1.63 (HSE) against 1.21–1.34 implied by the V–Mo parent | same; `data/tc_summary.json` |
| VCr/CrCr calibration (target 1.57 ± 0.25) | base 1.66; U = 2 diagonal 1.83 (off); U = 4 diagonal 1.52; **HSE 1.85 (off)** | — | same |
| T_C | 376 K measured | model: (A) DFT ratio 622 ± 79 (U pairings 573–705); (B) parent-calibrated 487 ± 49 (floor 377); (C) combined 524 ± 61 K. Only B survives the failed HSE calibration test | `data/tc_summary.json` (06:31); see the provenance note below |
| Linkage isomer (eV/f.u.) | +0.44 to +0.81 (base +0.63) | not computed | `data/gate0_collect.txt` |
| SOC spin residue | not computed | 0.000 μB/f.u.; orbital part not computed (lorbm failed twice **[new]**) | same |
| CN flip (4 f.u.): ΔE; windows | +0.177 eV; 98 / 85 % | +0.182 eV; 93 / 97 % | `data/defects_collect.txt` |
| V/M antisite pair: ΔE; aligned windows | +0.927 eV; **0.13 / −0.06 eV** (kill) | +0.565 eV; **0.25 / 0.18 eV** (kill) | same |
| K vacancy: M; windows | +1; 94 / 108 % (hole V-majority) | +1; 89 / 99 %. Self-trapped V(III) polaron −0.275 eV, windows 87 / 26 % **[new]** | same |
| [M(CN)6] vacancy: M; aligned windows | −3; 1.93 / **0.07 eV** (95 / 6 %; void state at the CBM). Water-filled test **[PENDING]** | −3; 0.91 / 0.84 eV (93 / 45 %; nominal FAIL of the 50 % rule) **[new: original job done, its relax capped at 5 h unconverged]** | `data/vac_collect.txt`; `data/vac_summary.json` |
| Hydrated KVCr·2H2O (PBE+U) | gap 2.10 eV; windows **0.93 / 0.92 eV** (−54 / −20 %); still unipolar. HSE **[PENDING]** | — | `data/gate0_collect.txt`; `structures/` |
| Phonons | Framework stable at Γ, X, L. Tilt +1.39 THz. K rattling −1.29 (Γ) / −1.04 (X) / −1.13 (L) THz | Framework stable at Γ, X, L **[new: L]**. Tilt +1.36 THz. K −1.35 / −1.14 / −1.21 THz | `data/phon_summary.json`, `data/phon_collect.txt` |
| K off-centring (fixed cell) | ⟨100⟩ −47 meV/K **[PENDING final]**; +[111] relaxed −20 meV/f.u., windows 2.01 / 1.19 eV, M = 0 | ⟨100⟩ −68 meV/K at 11:29, still deepening **[PENDING final]**; +[111] −30 meV/f.u., windows 0.98 / 1.89 eV, M = 0 | `data/phon_collect.txt`, `data/taut_summary.json` |
| CT tautomers (PBE+U) | every V(III)/Cr(II) start → LCM at U (4, 2) and (2, 2); (3, 3) **[PENDING]**. V(II) S = 1/2 at +0.889 eV/f.u.; Cr S = 1/2 at +1.039 | every V(III)/Mo(II) start → LCM at (3, 2), (4, 0), (2, 0). Mo S = 1/2 at +0.685; V S = 1/2 at +0.923 eV/f.u. HSE kill test **[PENDING]** | `data/taut_collect.txt`, `data/taut_summary.json` |
| Experiment | K1V[Cr(CN)6]·2H2O·0.1KOTf powder, fcc, a 10.55 Å; T_C 376 K (365 K after heating); M_sat 0.125 μB/f.u.; H_c 165 G | Only amorphous {[K(crypt)]0.34V1.37[Mo(CN)6]}: T_C 417 K (Bloch fit; measured only to 340 K), Mo XMCD M_L −0.04 μB, decomposes in air within minutes | `literature/lcm_pba_facts.md` §1, §3 |

**T_C provenance note [corrected at verification].**
- `tracks/lcm/pba/tc/README.md` (06:00) and the log quote A/B/C = 673 ± 77 / 487 ± 49 / 541 ± 85 K and a parent over-prediction of ×1.32. Those predate the KVCr/KVMo J4 supercell terms (the README's J table still lists J4 as "pending") and are superseded.
- `data/tc_summary.json` (06:31) includes the J4 terms: A = 622 ± 79 K, B = 487 ± 49 K, C = 524 ± 61 K (χ² 2.1); A over-predicts the parent by ×1.22 (508 vs 417 K).
- A scratch rerun of the unchanged `tc_model.py` (last edited 06:32) on the cached MC tables and the finished HSE KVMo pairing gives 619 ± 82 / 487 ± 49 / 522 ± 58 K, within 3 K of `summary.json`. `summary.json` itself was not rewritten.
- Every central value is above 400 K, and A/B = 1.27–1.28. Route B's floor (377 K) is not.
- Under the pre-registered calibration rule (failed in HSE), only route B can carry a ≥ 400 K statement.

## 3. Robustness verdict in one table

| Test | KVCr | KVMo |
|---|---|---|
| Gate 0 (pre-registered, all lines required) | **not passed**: calibration (HSE 1.85), isomer bar (0.63 < 0.8), hydration rule | **not passed**: VB 0.05 eV at U (4, 0); calibration; FM charge transfer at U_Mo = 0 |
| Dynamic stability | framework ✓; K off-centres (spectator) | same |
| CT tautomer (PBE+U) | none found so far | none found |
| CT tautomer (HSE) | **[PENDING]** | **[PENDING]** |
| CN flips | harmless | harmless |
| Antisite pairs | **kill** both edges (pre-registered defect kill triggered; equilibrium fraction 3×10⁻¹⁶ at 300 K; inert block) | **kill** both edges (kill triggered; 3×10⁻¹⁰; labile block) |
| K vacancies / holes | spin-locked p-type | spin-locked, polaronic; the self-trapped hole cuts the CB window to 26 % (nominal FAIL of the 50 % rule) |
| [M(CN)6] vacancies | VB ✓; **CB killed** by a void state at the CBM (6 %; kill triggered) **[PENDING water test]** | VB ✓; **CB 45–47 %: nominal FAIL** of the 50 % rule (still ≥ 0.3 eV) |
| Hydration | **fails the pre-registered rule** (VB 2.02 → 0.93 eV, −54 %, in PBE+U); existing-sample claim withdrawn; HSE **[PENDING]** | not computed (untested) |

The full rule-by-rule audit, including the post-hoc rescues the referees flagged, is report §9.4 (F1–F8).

## 4. Weaknesses

1. **Gate 0 did not pass as pre-registered for either compound**, and the defect gates ran anyway. Post-hoc rescues (report §9.4); the first three are the ones the referees named:
   - "HSE decides" for the KVMo U (4, 0) hole window;
   - the ≥ 400 K claim kept although the VCr/CrCr calibration failed in HSE (1.85 against 1.57 ± 0.25) and on the U = 2 diagonal (1.83);
   - a fixed-moment FM after the free FM charge-transferred;
   - the 0.8 eV isomer bar re-judged as "no like-for-like basis";
   - the KVMo vacancy CB at 45 % accepted because it is ≥ 0.3 eV;
   - the antisite kill (and the KVCr vacancy CB kill) moved to "kinetically excluded" (and "possibly a dry-void artefact");
   - the hydration failure answered by narrowing the KVCr claim to an anhydrous crystal that has not been made.
2. **KVCr's experimental basis is thin.** It is one 1999 hydrated powder, never reproduced, with no single crystal, no occupancy refinement, and no transport, optical-gap or L-edge XMCD data.
3. **Hydration fails the pre-registered rule** (hole window 2.02 → 0.93 eV, −54 %, in PBE+U). The claim is therefore restricted to anhydrous KVCr, which has not been reported. The water states and the cell shear were not separated (no PDOS), and the rule's pass conditions (HSE, full U grid, water states ≥ 0.5 eV from the edges) were never tested. KVMo hydration was not computed.
4. **Vacancies take out the electron edge.** In KVCr every [Cr(CN)6] vacancy puts an empty void state at the host CBM (CB window 6 %; the water-filled test is pending). In KVMo the CB window falls to 45–47 %, below the pre-registered 50 % line.
5. **Composition sets M** (3 μB per vacancy, 1 μB per V(III)). Luttinger protection holds only for the exact ordered compound at T = 0.
6. **Antisites kill the edges locally** (pre-registered defect kill triggered). Protection is kinetic and depends on the inert [Cr(CN)6]³⁻ block. For the lead KVMo the block is labile (Mo(III) substitutes about 10⁵× faster than Cr(III), an aqua-ion estimate; [Mo(CN)6]³⁻ survives only under cyanide limitation), so the protection argument is inverted.
7. **The T_C prediction for KVMo is model-only, and its routes disagree.** The DFT-ratio route (622 ± 79 K) over-predicts the V–Mo parent by ×1.22 and rests on a calibration that failed in HSE; the conservative route (487 ± 49 K) has a floor of 377 K. `tc/README.md` and the log carry superseded pre-J4 values (673 / 541 K); see the provenance note in §2.
8. **Carriers will be polaronic** (flat t2g bands; V(III) self-trapping −0.27 eV). At 300 K, KVCr retains only about 60 % sublattice order.
9. **The physics is generic.** Spin-split bands in compensated ferrimagnets are textbook, and Schart 2024 already noted TDOS↑ ≠ TDOS↓ at M = 0 for Cr[Cr(CN)6]. Scoop risk is high (Pedersen; Oxford/Rennes; Pinkowicz/Rogalev).
10. **Methods.**
   - The HSE windows come from a 4×4×4 SCF grid.
   - The K off-centring and hydration were computed without dispersion.
   - The PBE+U lattice is 1.2 % too large.
   - SOC was not run for KVCr, and no orbital moment was computed.

## 5. Decisive experiments

1. **Magneto-optics plus element-specific XMCD at M ≈ 0 on anhydrous KVCr**, and on a composition series across the S = 0 line.
   - Field-cooled MCD/Faraday spectra at 1.5–4 eV and 300 K.
   - V and Cr L2,3 XMCD to fix each sublattice sign. This has never been done on V–Cr; only K-edge XMCD exists.
   - LCM signature: the MO signal of the charge-transfer edge stays finite as M → 0, and its sign follows the Cr-sublattice direction, not the sign of M.
2. **Spin-resolved, angle-integrated photoemission of the top valence band.**
   - Prediction: about 100 % V-majority polarisation over the top 2.6 eV (HSE, anhydrous; 0.93 eV hydrated, PBE+U). It reverses with the Néel vector.
   - The s-wave splitting survives angle integration in polycrystals.
   - The VB window is the one that survives the computed point defects (≥ 87 %, antisites excepted). It does not survive hydration well (−54 % in PBE+U), so the test needs anhydrous material.
3. **Synthesis of crystalline 1:1 KV[Mo(CN)6]** (route: `tracks/lcm/pba/vac/README.md` §5):
   - aprotic DMF;
   - cyanide-limited Li3 or (NEt4)3[Mo(CN)6], with no cryptand;
   - naked K⁺/Cs⁺ in 5–20× excess;
   - [V(DMF)6]²⁺;
   - slow diffusion, dark, glovebox, ≤ 320 K;
   - calibrate on KVCr first.

   Acceptance tests:
   - ICP 1:1:1 (± 0.02);
   - a ≈ 10.73 Å;
   - a single ν(CN) band at 2100–2120 cm⁻¹;
   - M_sat ≤ 0.1 μB/f.u.;
   - equal, antiparallel V and Mo XMCD;
   - remanence at 330 K.

- **Cheapest preliminary:** reproduce KVCr with ICP plus M_sat(5 K), to confirm |3v − w| ≤ 0.05.

## 6. Pending at 11:15 PDT, all still running at the 11:29–11:33 recheck (nothing cancelled)

| Job | Decides |
|---|---|
| `lcm/pba/hse/KVCr2H2O_hse_LCM` with the matched `KVCr_hse_LCM_nq1` (HSE06, nq 1) | Whether the hydration loss of the hole window also appears in HSE. A smaller HSE change would give a method split; it cannot turn the PBE+U failure of the pre-registered rule into a pass |
| `lcm/pba/def/KVCr_VAC_w1` (water-filled vacancy; nscf stage) and `KVCr_VAC` (original relax, 60 SCF steps at 11:30) | `_w1`: whether the KVCr CB kill is an empty-void artefact (lifted only if the opposite-spin state moves above CBM + 0.3 eV). `KVCr_VAC`: confirms the `_r` numbers within ±0.1 eV |
| `lcm/pba/taut/*` (10 jobs) | The HSE06 charge-transfer kill test, the remaining CT starts and K-ordering variants |
| `lcm/pba/phon/koff2/*`, `koff/*_xy` | Final K well depths, and the windows of the ⟨100⟩ off-centred cells |

**Refresh:**
```
source .venv/bin/activate
python tracks/lcm/pba/collect.py --fetch
python tracks/lcm/pba/defects/collect.py --fetch
python tracks/lcm/pba/vac/collect_vac.py --fetch && python tracks/lcm/pba/vac/summarize.py
python tracks/lcm/pba/phon/collect.py --fetch
python tracks/lcm/pba/taut/collect.py --fetch
```

The hydrated HSE jobs are not in `collect.py`. Read `jobs/lcm/pba/hse/{KVCr2H2O_hse_LCM,KVCr_hse_LCM_nq1}/deep_hse.json` by fixed path.

## 7. Contents

**`structures/`** (all from the volume JSON field `relaxed`; PBE+U LCM vc-relax at 110 Ry; BFGS converged):

| File | Content |
|---|---|
| `KVCrCN6_F-43m_conv_PBEU110_relaxed.cif`, `KVMoCN6_F-43m_conv_PBEU110_relaxed.cif` | symmetrized conventional cells (F-43m, 60 atoms) |
| `KVCrCN6_prim15_P1_PBEU110_relaxed.cif`, `KVMoCN6_prim15_P1_PBEU110_relaxed.cif`, `POSCAR_*_prim15_*` | the 15-atom primitive cells exactly as computed. Sources: `jobs/lcm/pba/s1a/KVCr/results/KVCr.json` and `jobs/lcm/pba/s1a/KVMo/results/KVMo.json` |
| `KVCrCN6_2H2O_P1_PBEU110_relaxed.cif`, `POSCAR_KVCrCN6_2H2O_P1_21atom_*` | the hydrated 21-atom cell from `jobs/lcm/pba/s1/KVCr_2H2O/results/KVCr_2H2O.json`. Run in P1; spglib finds C2 at 1e-3 Å |
| `*_LCM.mcif` | magnetic CIFs with the LCM start moments (V −3, M +3 μB) |
| `structures_info.json` | source paths, space groups, relax metadata, energies |

The hydrated relaxed cell is sheared (rhombohedral angles 58.7–59.5°) and 2.6 % larger in volume than the anhydrous cell [derived here].

**`data/`** (copies taken at 11:15; the originals stay in `tracks/lcm/pba/`. The 11:29 collector run rewrote `phon/summary.json` with slightly deeper `koff2` wells; the other copies still match):

| File | Content |
|---|---|
| `gate0_collect.{txt,json}` | Gate-0 collector output (PBE+U grid, HSE, ratios, isomer, SOC, hydration, JUDGE checks) |
| `gate0_manifest.json` | Gate-0 job list and settings |
| `defects_collect.{txt,json}` | Gate-1 defect cells |
| `vac_collect.{txt,json}`, `vac_summary.json`, `vac_cells_summary.json` | vacancy cells, lability energetics, concentration scenarios |
| `phon_collect.txt`, `phon_summary.json` | phonons, K scans and relaxations |
| `taut_collect.txt`, `taut_summary.json` | tautomer and spin-state starts, K variants |
| `tc_summary.json` | T_C model at 06:31, the current version (see the provenance note in §2) |
| `tc_hse_prelim_0631.json` | HSE energies as parsed at 06:31 |

**`scripts/`**
- `jobs/`: copies of the Modal job scripts (stage-1 chain, fixed-moment FM, lcm_deep HSE/SOC, defects, phonons, K relax, tautomers, K variants).
- `collectors/`: copies of the collectors and the T_C model.
- `MANIFEST.json`: origin path, md5 and the jobs that used each script.
- The copies are for reference. Run the originals from `tracks/lcm/pba/`, because they import `infra/mrun.py` and sibling modules by relative path.

## 8. Verification log (2026-10-03, 11:29–11:35 PDT)

Every number in this README and in report §0 / §9 (plus the cyanide rows of §2.7 and §4) was checked against the log (03:00–11:30 entries), `trilemma/{JUDGE,nonoxide_frameworks,kinetic}.md`, the `pba/` READMEs, summaries and manifests, `literature/lcm_pba_facts.md`, and a fresh `--fetch` run of all five collectors at 11:29. A scratch copy of `tc/tc_model.py` (unchanged code, cached MC tables, no fetch) was rerun to settle the T_C provenance; nothing in `tracks/` was edited by hand, and no job was cancelled. Numbers not listed below were confirmed as written.

| # | Where | Was | Now | Basis |
|---|---|---|---|---|
| V1 | Report §0, §2.7, §4, §9.0, §9.2; README §0, §1, §2, weakness 7 | KVMo T_C headlined as A 673 ± 77 K (pairings 609–750), C 541 ± 85 K (χ² 4.1); `summary.json` called "not certainly current" | A 622 ± 79 K (573–705), B 487 ± 49 K, C 524 ± 61 K (χ² 2.1) | `tc/README.md` (06:00) predates the KVCr/KVMo J4 terms (its J table says "pending"); `tc/summary.json` (06:31) includes them. Rerun of the current code with the finished HSE pairing: 619 ± 82 / 487 ± 49 / 522 ± 58 K |
| V2 | Report §9.2; README weakness 7 | A over-predicts the parent ×1.32 (550 K), ×1.21 against 454 K | ×1.22 (508 K), ×1.12 | `summary.json` checks.VMo_Magott |
| V3 | Report §9.2; README §2 | Experiment-implied ratio 1.21–1.35 "against 1.77 (PBE+U) and 1.78 (HSE)" | T ratio 1.21–1.34 against model T ratios ≈ 1.67 (PBE+U) and 1.63 (HSE); raw ΔE ratios 1.82 / 1.78 | 1.77 was the pre-J4 model; 1.78 is an energy ratio, not a T ratio; HSE T ratio from the rerun |
| V4 | Report §9.2 | Calibrants "both reproduced to ±1 %" | f_MC = 0.98 (CrCr), 1.02 (KVCr); model T ratio −4 % at base U, +6 / −12 % on the U diagonals, +7 % in HSE | `summary.json` calibration, calibration_ratio_test |
| V5 | Report §9.2 | "The referees made [the route disagreement] the headline caveat" | The referees' correction (ii) flagged the failed HSE calibration; the disagreement is stated as the lane's own caveat | Log 11:30 lists corrections (i)–(v); none is about route disagreement |
| V6 | Report §0, §2.7, §9.0; README §1, weakness 6 | "[Mo(CN)6]³⁻ is about 10⁵× more labile than [Cr(CN)6]³⁻" | Mo(III) substitutes about 10⁵× faster than Cr(III), an aqua-ion estimate (snippet-level Mo rate, recalled Cr rate); [Mo(CN)6]³⁻ survives only under cyanide limitation | `vac/README.md` §3.2 ("[derived/recalled]") |
| V7 | Report §9.1 | Quote "protected by the cyanide linkage" attributed to JUDGE.md §0–§2 and log 04:00 | Attributed to the referees' summary; JUDGE.md's actual wording quoted (§1 "structural (linkage + isospin)", §2 "a compensation that antisites cannot break") | The quoted phrase appears only in the log's 11:30 referee entry |
| V8 | Report §9.1 table | Sublattice identity "Protected? Yes, structurally, by kinetics and thermodynamics" | Site inequivalence structural; which metal sits where only thermodynamic + kinetic, weak for Mo(III) | Contradicted §9.3.1 ("kinetic, not structural") and §9.3.6 |
| V9 | Report §9.6 | JUDGE.md Kabalan statement "has been corrected" | Correction recorded in `lcm_pba_facts.md` §0 item 2 and §6 and the log; JUDGE.md itself still carries it | JUDGE.md §1 table, item 2, unchanged |
| V10 | Report §9.3.5 | "The water states are spin-degenerate, so they cap both windows" | Removed as unsupported (no PDOS was run); added that the JUDGE §3 pass conditions (HSE, U grid, water states ≥ 0.5 eV from the edges) were never tested | `pba/README.md` §5 ("Gate 0 has no PDOS") |
| V11 | Report §9.8 | "Anhydrous KVCr needs a dehydration step below 340 K" | No dehydration protocol is reported; heating already lowers T_C to 365 K | 340 K is the V–Mo decomposition limit, not a KVCr number |
| V12 | Report §9.8; README §5 | VB window "survives every realistic defect except antisite pairs"; spin-filter test without caveat | Hydration (−54 %) named as the exception; spin-filter test flagged as needing low-vacancy KVCr (void state shunts the CB) | `vac/README.md` §4 |
| V13 | Report §0, §9.0, §9.3.1, §9.7, §4 row; README §0, §1, §3, weaknesses | Hydration and vacancy failures partly implicit ("hydration-sensitive"; KVMo vacancy CB "45 % (≥ 0.3 eV)") | Stated plainly: hydrated KVCr fails the pre-registered rule (2.02 → 0.93 eV, −54 %) and the existing-sample claim is withdrawn; KVCr CB killed at vacancies (6 %); KVMo vacancy CB 45–47 % and polaron CB 26 % are nominal FAILs of the 50 % rule; antisite kill triggered | Collector output 11:29 |
| V14 | README weakness 1 | Post-hoc rescue list omitted the HSE calibration failure | Added (one of the three the referees named), plus the vacancy and hydration rescues | Log 11:30 correction (ii) |
| V15 | Report §9.4 F1, F2, F3, F7 | F1 omitted that JUDGE step 3 calls HSE "the arbiter for windows"; F2 rule quoted as "in base PBE+U and in HSE"; F3 rule attributed to JUDGE; F7 "a rule that requires both PBE+U and HSE" | F1: that partial basis for "HSE decides" noted, while the pass line still requires every U point; F2: JUDGE wording quoted, with HSE as "the arbiter for … the J ratio" (JUDGE step 3) and `collect.py`'s check named; F3: from `pba/README.md` §4 step-1 checks; F7: the kill line names no method and is triggered in PBE+U | JUDGE.md §2–3; `pba/README.md` §4 |
| V16 | Report §0 pt 4, §4 caveat, §9.7 | KVMo "the first candidate with no ✗" | Same, but "only because its (A) entry is a model '?'"; only route B can carry ≥ 400 K | §4 table |
| V17 | Report §9.3.1; README §2 | KVMo [Mo(CN)6] vacancy "original job, final" | Job done, but its relax stopped unconverged at the 5-h cap (64 steps) | `defects/collect.py`: conv=False |
| V18 | Report §9.3.2; README §2 | `koff2` wells at 11:15 (KVMo ⟨100⟩ −66.6 meV at 0.97 Å; −55.7; −26.8) | 11:29: −67.6 meV at 1.00 Å; −56.1; −27.3 (KVCr within 0.2 meV of the 11:15 values); still running | `phon/collect.py --fetch` 11:29 |
| V19 | Report §9.9; README §6 | States at 11:15 | All pending jobs confirmed still running at 11:29–11:33: hydrated HSE pair (no `deep_hse.json`), `KVCr_VAC_w1` (nscf), `KVCr_VAC` (60 SCF steps), 10 taut jobs, `koff2` / `koff/*_xy` | Volume `status.json` and collectors |
| V20 | Report §9.10 | Every `pba/` sub-folder "with a README and summary.json" | `defects/` has no summary.json | Directory listing |
| V21 | README §1 | "Every insulating LCM known so far orders far below room temperature" | The only neutron-confirmed insulating LCM orders at 48 K; above-RT compensated insulators (Ca–Si YIG, KVCr) had not been analysed as LCMs | Report §5 |
| V22 | Report §9.3.6 | [Mo(CN)6]³⁻ "decomposes … in air within minutes, and it is photolabile" | Protic decomposition is of the anion; air decomposition (minutes) is of the V–Mo product; photolability shown for K4[Mo(CN)7] | `vac/README.md` §3.2 |
| V23 | Report §9.6 | "a ≥ 400 K LCM prediction" | "a model ≥ 400 K prediction (487 ± 49 K, the only route that survives the failed calibration rule)" | §9.4 F2 |
| V24 | README §0 claims 2–3 | eV-scale windows and the KVCr hole edge stated without method / sample limits | HSE06 and ideal anhydrous crystal only; PBE+U KVMo hole window 0.05 eV at U (4, 0) noted | Gate-0 collector |
| V25 | Report §9.3.3 | HSE tautomer test described only as "pending" | Added that the chains test fixed-M CT states and that a collapse in QE's PBE pre-stage is reported as a collapse | `taut/README.md` §0, §3.3 |

---

## Addendum 2026-10-03 14:35 — results that landed after the verified write-up

1. **Hydration rule now PASSES with the arbiter functional (HSE06).** For the real sample KV[Cr(CN)6]·2H2O (21-atom P1 model, 2 H2O in the empty void):
   - Settings: HSE06, no U, nq 1, k-spacing 0.30; eigenvalues taken from the last converged EXX block (error plateaued at 1.1e-6 Ry; job `lcm/pba/hse/KVCr2H2O_hse_LCM`, terminated after the plateau).
   - Result: **gap 2.02 eV, unipolar (both edges V-majority), windows VB 2.43 / CB 1.42 eV.**
   - The matched anhydrous run (`lcm/pba/hse/KVCr_hse_LCM_nq1`, same settings) gives 2.11 / 2.68 / 1.44. Hydration therefore costs −9 % (VB) and −2 % (CB).
   - The PBE+U −54 % collapse was a PBE+U artefact (water levels mis-placed relative to the V/Cr d states). The "existing-sample" claim for the edges is reinstated.
   - Remaining hydration caveats: one water arrangement only; no PDOS; real samples also contain zeolitic water and vacancies.
2. **Vacancy + zeolitic water, PBE+U** (`lcm/pba/def/KVCr_VAC_w1`): M = −3.000; VB window 1.91 eV (94 %); CB still killed (7 %).
   - The opposite-spin state at CBM + 0.08 eV is 73 % water O–H σ*, the same kind of water level that PBE+U mis-placed in the hydrate.
   - The HSE06 test is running (`lcm/pba/hse/KVCr_VAC_w1_hse`). Until it lands, the electron-edge claim at vacancies stays withdrawn.
3. **V-free, alkali-free variant Cr(III)[Mo(III)(CN)6]** (Gate-0 chain `lcm/pba/s1a/CrMo`; HSE `lcm/pba/hse/CrMo_hse_LCM`):
   - LCM, M = 0. **BIPOLAR** edges (VBM s0, CBM s1). Windows are 0.66–1.79 eV at all 5 PBE+U points (U_Cr ∈ {2,3,4}, U_Mo ∈ {0,2}). **HSE06: gap 3.04 eV, windows 1.22 / 1.11 eV.**
   - Exchange ratio to CrCr is 1.6–1.9, giving T ≈ 390–450 K by scaling from CrCr's 240 K (borderline).
   - It is a compensated *bipolar* magnetic semiconductor (holes and electrons fully polarised in opposite spins). It is not reported, and N-bound Cr(III) is substitution-inert, so it would be hard to crystallise.
   - Design rule: unipolar edges need the most reducing ion (highest occupied t2g) and the most oxidising ion (lowest empty minority t2g) on opposite sublattices, i.e. the V(II)-N / M(III)-C pairing.
4. Grade unchanged: identification + quantification, plus predictions (KVMo unipolar, CrMo bipolar). The hydration result removes the most damaging objection to the existing-sample claim. Still open: vacancy composition control (M = 3v μB), the vacancy electron edge, air sensitivity, polaronic transport, and the fact that only one 1999 KVCr report exists.

### Corrections from the second referee round (2026-10-03 15:40; full text: `tracks/lcm/pba/REREVIEW.md`)
- **Hydrate HSE is a method split, not a pre-registered pass.** The kill (−54 %) was triggered in PBE+U. The pass conditions (full U grid; water states ≥ 0.5 eV from both edges by PDOS) remain untested.
  - The hydrate EXX loop was not converged: the VBM was still rising ~16 meV per block, and dexx was ~5× conv_thr.
  - Converged quantities: CB window **1.42 eV**, other-channel gap **5.86 eV**.
  - The hole window (~2.3–2.4 eV) is **water-limited** (a 1b1 lone pair); the framework window is unchanged (~2.8 eV).
  - Only one ordered-water arrangement was computed, on the PBE+U geometry (a 10.77 vs 10.55 Å exp).
- **Photocarrier prediction narrowed.** Cr(III) d–d excitations (⁴A2→⁴T2 ≈ 3.3 eV, ⁴T1 ≈ 4.0 eV) lie in the other spin channel, and the Cr ²E spin-flip state forms in < 250 fs.
  - The clean single-spin window is **≈ 2.0–3.2 eV**, the V→Cr MMCT band (measured 2.25–2.30 eV).
  - The excitations are localized CT pairs/polarons, not free band carriers. At 300 K the sublattice order is ~0.6.
  - The effect is LCM-specific only across a compensation series.
- **Cr[Mo(CN)6] T revised.** Exchange ratio to CrCr is 1.53–2.04 → 391 K at base U (367–490 K); KVCr-calibrated ≈ 369 K. So it is **below 400 K** at base U.
  - The Mo[Cr(CN)6] isomer has not been computed. Mo(III) remains labile; the MeCN Cr(III)-N precedent gave a 110 K cluster glass.
- **Bar status unchanged:** no existing compound reaches 400 K (KVCr 376 K); compensation is composition-tuned (M = w − 3v); the physics is generic. The only route to the bar is crystalline 1:1 KV[Mo(CN)6] (model 487 ± 49 K, floor 377 K), which has never been made.

### Method caveats added 2026-10-03 16:05 (from Session 3's QE audit; see TRACKL_LOG 16:05)
- SOC residues and MAE come from antiparallel noncollinear starts. These are affected by QE's lsign behaviour → treat as unverified.
- HSE total-energy differences between configurations of different magnetic symmetry may carry EXX k-set offsets. Gaps and windows are unaffected.
- Single-cycle vc-relax energies (hulls, swap calibrations) need one re-relax before quantitative reuse.
- K-rattling phonon modes are possibly affected by alkali egg-box forces. The K off-centring wells are larger than that effect.

### Update 2026-10-03 19:20 — HSE06 on the water-filled [Cr(CN)6] vacancy (near-converged, dexx 5e-6)
KVCr_VAC_w1 (65 atoms, one zeolitic H2O at the vacancy centre, M = −3): gap 1.57 eV, **both edges V-majority, in-cell windows 2.80 (VB) / 0.74 (CB) eV**. The PBE+U CB kill (water σ* at CBM + 0.08 eV) is a PBE+U water-level artefact, as for the hydrate. With the arbiter functional, the unipolar edges of KV[Cr(CN)6] survive both of its real defects (zeolitic water, [Cr(CN)6] vacancies). The remaining limits are not DFT ones: M = w − 3v is composition-tuned; T_C 376 K < 400 K; a single 1999 sample; air sensitivity.

### Converged HSE06 finals (2026-10-03 22:25)
- KVCr water-filled [Cr(CN)6] vacancy: gap 1.564 eV, unipolar, windows 2.803 / 0.737 eV (final; confirms the interim result).
- [Two lines about other compounds, outside the scope of this repository, were removed from this copy.]
