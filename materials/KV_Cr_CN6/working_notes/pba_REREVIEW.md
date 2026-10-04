# Track L re-review: hydrate HSE, water-filled vacancy, Cr[Mo(CN)6] and the photocarrier claim

Session 66190d (Track L), 2026-10-03, about 15:40 PDT.

This re-review synthesises the second three-referee round (realizability, physics, impact lenses). It covers the results that landed after the verified 11:30 write-up, as reported in:
- the README addendum and the report addendum, both 14:35;
- log entries 11:55–14:55.

Where the referees disagree, the disagreement is resolved below against the raw job outputs; the checks are listed in §7. The README, the report and the log were **not** edited. The corrections they need are in §5.

---

## 0. Verdict

| | |
|---|---|
| Stopping condition met | **No** (3 of 3 referees) |
| Grade | **Identification plus quantification.** Venue: PRB / PR Materials / Chem. Sci., or JACS only if paired with XMCD/MCD data. **Not big news.** |
| P(big news) | **≈ 4 %** (referees 4 / 4 / 5 %; range 3–6 %). At 11:30 it was 2–5 %, so it is effectively unchanged |
| What improved | Inside the same grade, the KVCr claim now plausibly covers the compound that actually exists: the reported dihydrate, the 1999 sample. It does so as a **method split** resolved in favour of the designated arbiter functional, **not** as a pre-registered pass |
| What did not move | The bar fails on independent grounds: <br>• no compound that exists reaches 400 K; <br>• no LCM-specific property has been measured; <br>• compensation is composition-tuned; <br>• the physics is generic, so the gain over prior art is quantitative, not conceptual |

**Why the hydrate result matters even though the grade stays.**
- The decisive experiment no longer needs an anhydrous crystal that nobody has made.
- A reproduction of the 1999 anaerobic aqueous synthesis is enough.
- That roughly doubles the feasibility of the KVCr demonstration path.

**Why it still does not reach the bar.**
- Even a perfect KVCr experiment would demonstrate a 376 K compensated magnet whose M is set by composition.
- That is below the lane's T ≥ 400 K bar.
- The only route to the bar as written is a crystalline 1:1 KV[Mo(CN)6] or another ≥ 400 K d³/d³ PBA. No new result touches that route.

---

## 1. The four new results and what each is worth

### 1.1 HSE06 on KV[Cr(CN)6]·2H2O: accepted in substance, but the numbers and the "PASS" wording must change

**Reparsed trajectories.** Both runs are reparsed below. Block 1 is QE's PBE pre-stage; blocks 2 onward are the EXX outer iterations. Each cell gives gap / win_VB / win_CB in eV.

| SCF block | Hydrate (`KVCr2H2O_hse_LCM`) | Matched anhydrous (`KVCr_hse_LCM_nq1`) |
|---|---|---|
| 1 (PBE) | 0.615 / 1.410 / 1.346 | 0.550 / 1.400 / 1.348 |
| 2 | 2.163 / 2.160 / 1.476 | 2.138 / 2.370 / 1.479 |
| 3 | 2.239 / 2.190 / 1.428 | 2.233 / 2.495 / 1.460 |
| 14 | **2.016 / 2.432 / 1.415** (last block; job aborted) | 1.930 / 2.875 / 1.447 |
| 15 | — | 1.926 / 2.877 / 1.446 (turning point) |
| 31 | — | **2.114 / 2.678 / 1.440** (converged, dexx 1.0e-7 Ry) |

**The hydrate EXX loop is not converged.**
- dexx sat at 1.06–1.15e-6 Ry, about 5× the conv_thr of 2.1e-7.
- The V-majority VBM was still rising by 16–17 meV per block when the job was aborted.
- The matched anhydrous run went through the same slow mode: it overshot by about 0.2 eV and then relaxed back.
- So the quoted 2.02 / 2.43 eV and −9 % compare a transient with a converged number. This resolves the disagreement between referee 1 ("plateau, does not matter") and referees 2 and 3 ("not converged") in favour of referees 2 and 3.
- What is converged: win_CB (1.415–1.417 eV), the other-channel gap (5.864 eV) and the water level relative to the CBM (−4.45 eV from block 8 onward).

**What sets the hydrate hole window.** Band counting at block 14, energies relative to the CBM:

| Level | Energy |
|---|---|
| Opposite-channel top: one water 1b1 lone pair, nearly spin-degenerate (≤ 0.05 eV splitting) | −4.45 eV |
| Cr t2g↑ | −4.84 eV |
| V t2g↓ VBM | −2.02 eV |

- In HSE the hole window is therefore **water-limited**. The framework window (Cr t2g↑ to VBM) is about 2.8 eV, unchanged from the anhydrous cell. Hydration costs the framework nothing; a water level enters the window.
- In plain PBE the hydrate loses nothing (1.41 against 1.40 eV), and Cr t2g↑ sets the window in both cells.
- In PBE+U the water level sits only 0.93 eV below the VBM, because +U pushes the d manifolds down and leaves the water level alone [inference from the band ordering; no PDOS yet].
- **PBE+U is the outlier of three functionals.** This supports the artefact explanation, which the 14:35 addendum asserted without evidence.

**Best converged estimate.**
- At block 3 the hydrate and anhydrous gaps agree to 6 meV. The anhydrous gap ended 0.12 eV below its block-3 value, which puts the hydrate gap at ≈ 2.1 eV (2.0–2.2).
- The VB window is then 4.45 − 2.1 ≈ **2.35 eV (2.25–2.45)**. The CB window is **1.42 eV**.
- Hydration therefore costs **≈ −12 % (−9 to −16 %) on the VB edge** and −2 % on the CB edge.
- The kill line (> 50 %, i.e. win_VB < 1.34 eV) is cleared by at least 0.8 eV in every EXX block, so the qualitative result is robust to the convergence problem.

**Status under the pre-registered rule: a method split, not a PASS.** The JUDGE §3 rule is quoted in full here:

> Pass: windows ≥ 0.3 eV on both edges, unipolar, in HSE and over the whole U grid; water states ≥ 0.5 eV away from both edges; predicted Δ_opt = win_VB + win_CB ≥ 0.6 eV (the MCD observable).
> Kill: as item 1. If the hydrated cell shifts either window by > 50 %, kill the "existing sample" claim.

- The kill line names no method. It **was triggered in PBE+U** (−54 %).
- The pass line is met in HSE and at base U. The **whole U grid was never run** for the hydrate, and the water separation comes from band counting, not from a PDOS.
- At 11:30, §9.4 F7 and §9.9 pre-committed that an HSE change ≤ 50 % "would give a method split, not a pass". The 14:35 "PASSES" reverses that after the fact, so it must be logged as post-hoc rescue **F9**.
- Correct wording: *"Method split. The hydration kill is triggered in PBE+U only (−54 %), not in PBE (≈ 0 %) or in HSE06 (−9 to −16 %), and HSE06 is the designated arbiter for windows. The existing-sample claim is reinstated with that caveat; the formal pass conditions (U grid, PDOS) are untested."*

**Limits of the hydrate model**, all noted by the referees:
- **One ordered-water arrangement.** The two waters' 1b1 levels already differ by 0.65 eV. The least-bound water in a disordered hydrate could lower the window by a few tenths of an eV.
- **Sheared cell.** The arrangement shears the cell (58.7–59.5°), while the powder is cubic by PXRD. It also halves the Cr t2g↓ CB width (0.39 against 0.78 eV).
- **Geometry and settings.** PBE+U geometry: a_eq 10.77 Å against 10.55 Å measured. nq 1, no dispersion.
- **Direction of the remaining functional error.** HSE probably still underbinds the water lone pair relative to the localized d levels: gas-phase or liquid water 1b1 lies roughly 4–5 eV below a V(II) t2g level [inference]. If so, the true window lies between the HSE value and the framework value of about 2.7–2.8 eV. Plausible overall range: **1.8–2.8 eV**.

**Calibration side effect** (referee 2):
- The hydrate's PBE+U exchange energy is 128 meV/ion, against 152 meV/ion anhydrous.
- So the KVCr T_C calibration point, which pairs the hydrated T_C with the anhydrous ΔE, is ambiguous at about the 15 % level for route A.
- Route B (calibrated on the V–Mo parent, 487 ± 49 K) is unaffected.

### 1.2 Water-filled [Cr(CN)6] vacancy: electron edge still unresolved, with a provisional HSE reading

**PBE+U (`KVCr_VAC_w1`), as reported and confirmed by referee 1:**
- M = −3.000; VB window 1.91 eV (94 %).
- The CB is killed: the lowest opposite-spin empty state (band 159) sits at CBM + 0.08 eV after alignment. Its weight is 53 % H_w, 20 % O_w and 12 % H_z, and it is 95 % localized at the vacancy.
- Water also mixes into the V-majority CBM in that cell (band 162: 21 % H_w).
- The relax stopped unconverged at its 45-min cap (last-step displacement 0.04 Å). The HSE single point therefore runs on an unconverged geometry.
- Referee 3's point: the state sits at CBM + 0.07 eV in the **empty** void cell and at + 0.08 eV with water. It is therefore a cap O–H / cavity state, not specifically a misplaced zeolitic-water level.
- The 14:35 claim that it is "the same kind of water level that PBE+U mis-placed in the hydrate" is not supported. The "PBE underbinds σ*" rationale in the log (12:05) cuts both ways.

**HSE06 (`KVCr_VAC_w1_hse`), provisional and not a harvest.**
- The job was still running at 15:35 PDT (1.27 h in). The output fetched read-only has the PBE stage plus the first EXX block.
- In-cell values relative to the host Cr t2g↓ CB manifold (bands 163–171, nine bands = 3 Cr × t2g, so the CBM is the host CBM):

| Stage | Gap | win_VB / win_CB | Opposite-spin water/void state (up, band 160) | Same-spin partner (dn, band 172) |
|---|---|---|---|---|
| PBE | 0.28 | 1.39 / 0.91 | CBM + 0.91 | + 0.85 |
| EXX block 1 | 1.33 | **2.73 / 0.76** | **CBM + 0.76** | + 0.68 |

- **How much further the state should move.** In the pristine and hydrate runs, the analogous empty diffuse band moved −0.52 eV relative to the CBM in the first EXX block and only −0.09 eV more by convergence (pristine +2.69 → +2.17 → +2.08). Extrapolating, the vacancy state should converge at about **CBM + 0.65–0.75 eV** in-cell.
- **Forecast against the pre-registered lines:**
  - **Kill lifted** (opposite-spin state ≥ CBM + 0.3 eV): P ≈ 0.85. The residual risk is the alignment: the in-cell value is not the aligned value, and the PBE+U alignments differed by up to 0.3 eV.
  - **50 % rule met:** P ≈ 0.3. There is no matched 4-f.u. pristine HSE run. The nearest references are the 15-atom nq 1 run (1.44 eV, giving a 0.72 eV threshold) and the nq 2 run (1.57 eV, giving 0.79 eV). The forecast value sits on or just under that band.
- **Likely outcome:** the electron edge survives water-filled vacancies at about half the pristine window (≈ 0.6–0.8 eV). That is the same status as KVMo's 45–47 %: a nominal fail of the 50 % rule, but ≥ 0.3 eV.
- **Open item:** the in-cell gap is 1.33 eV, against 2.14 eV for the pristine at the same EXX stage. Same-spin vacancy levels narrow it; they are not yet assigned to VB or CB. They do not violate the spin windows, but they are trap levels.

### 1.3 Cr(III)[Mo(III)(CN)6]: demote to a design-rule example

The numbers below were verified by referees 1 and 2 from the five PBE+U result JSONs and `deep_hse.json`.

| U point (Cr, Mo) | PBE+U windows (eV), bipolar | ΔE ratio to CrCr at matched U_Cr |
|---|---|---|
| (3, 2) base | 1.01 / 0.90 | 1.63 |
| (2, 0) | 1.23 / 1.26 | 1.87 |
| (2, 2) | 0.70 / 0.66 | 1.53 |
| (4, 0) | 1.79 / 1.75 | 2.04 |
| (4, 2) | 1.26 / 1.14 | 1.71 |

HSE06 gives M = 0, a gap of 3.04 eV and windows of 1.22 / 1.11 eV. That gap is spin-forbidden: the spin-allowed gaps are 4.15 and 4.26 eV.

**Temperature.**
- The ratio range is **1.53–2.04**, not the quoted "1.6–1.9".
- Scaling from CrCr (240 K) gives 367–490 K, or 391 K at base U. But that same scaling over-predicts KVCr itself (398 K at base, 444 K in HSE, against 376 K measured).
- Referenced instead to KVCr, which has almost the same ΔE (149.4 against 152.2 meV/ion), it gives **≈ 369 K**. The model's ratio compression (×0.92 for KVMo) gives ≈ 350–420 K.
- **Centre ≈ 370 K, below the 400 K bar.** There is no HSE FM run.

**Physics.**
- The edges are bipolar, so there is no spin-selective photoexcitation.
- It is conceptually close to the cyanide LCM of Guo et al., PRB 112, L180405 (2025).

**Makeability.** It is worse than the addendum says:
- Mo(III) is still protic-, air- and photo-labile.
- N-bound Cr(III) must come from a labile Cr(III) solvate. The only aprotic precedent (Nelson 2007, [Cr(NCMe)6]³⁺) gave a CrCr cluster glass at 110 K, against 240 K for the crystalline route.
- Schart's crystalline route uses NMF, which is protic and unlikely to tolerate [Mo(CN)6]³⁻.
- Both ions are C-philic, and the Mo[Cr(CN)6] linkage isomer was never computed.
- One plus: the 1:1 M(III)[M(III)(CN)6] framework is charge-balanced without vacancies.

**Treatment:** a footnote or design-rule example (unipolar edges need the most reducing ion and the most oxidising ion on opposite sublattices), not a candidate. No further compute is needed unless it is promoted (§4, K6).

### 1.4 The "100 % spin-polarized photocarriers at 2.0–5.8 eV" prediction: restate it narrowly

All three referees agree on three points:
- It is a corollary of unipolar edges, not an independent result.
- Its energy range is wrong.
- "Carriers" overstates what is excited.

**Energy range.** The 5.86 eV upper bound is a Kohn–Sham channel gap. Several real excitations sit inside it:

| Excitation | Energy | Spin channel and consequence |
|---|---|---|
| Cr(III) ⁴A2→⁴T2 | ≈ 3.3 eV [recalled, from [Cr(CN)6]³⁻] | Spin-allowed in the opposite (Cr-majority) channel |
| Cr(III) ⁴A2→⁴T1 | ≈ 4.0 eV [recalled] | Same |
| LMCT band of V–Cr films | ≈ 3.1 eV (Johansson 2016, `lcm_pba_facts.md` §1) | Channel unassigned |
| Formation of the Cr ²E spin-flip state | — | A 3.1 eV pump forms it in < 250 fs (Johansson 2016, per referee 1) |

- The HSE single-particle Cr t2g↑→eg↑ gap is about 6.7 eV, so these Frenkel excitations are invisible in the band picture.
- **Clean single-channel range: ≈ 2.0–3.1 eV.** This is the V→Cr MMCT band, measured at 2.25–2.30 eV in V(II)–Cr(III) PBAs. That agreement with the HSE V-majority gap is supporting evidence for the edge assignment.

**Nature of the excitation.**
- It is a metal-to-metal charge-transfer pair across 5.3 Å with ε∞ ≈ 2.3, giving a binding of about 1 eV (referee 2 estimate). It is an exciton, not a free pair.
- The bands are flat, and V(III) holes self-trap (−0.27 eV).
- Each pair has ΔS_z = 0, so a photo-spin-current will be small.

**Macroscopic polarization.**
- It scales with sublattice order, about 0.6 of saturation at 300 K.
- It averages to zero over Néel domains, and field-cooling selects domains less and less well as M → 0.
- It is LCM-specific only across a compensation series. At a single M = 0.125 μB/f.u. it is an ordinary ferrimagnet-semiconductor effect.

**Defensible version:** spin-selective (V-majority-only) interband absorption from ≈ 2.0–2.2 eV up to ≈ 3.1 eV. Its magneto-optical sign follows the Néel vector and its magnitude follows the sublattice order, so it stays finite through M = 0. Prefer L-edge XMCD, MCD and spin-resolved PES over photocurrent as probes.

---

## 2. Strongest claim Track L can defend now

> **KV[Cr(CN)6]·2H2O** is the compensated Prussian-blue analogue of Holmes & Girolami (1999): T_C 376 K, M_sat 0.125 μB/f.u., one unreproduced powder. In HSE06 it is a **Luttinger-compensated (M = 0 by counting), s-wave spin-split semiconductor whose valence- and conduction-band edges both lie in the V-majority spin channel**.
>
> | Cell | Gap (eV) | Windows VB / CB (eV) |
> |---|---|---|
> | Ideal anhydrous crystal | 2.09–2.11 | 2.64–2.68 / 1.44–1.57 |
> | One ordered-water model of the reported dihydrate | ≈ 2.1 (2.0–2.2; EXX loop not fully converged) | ≈ 2.35 (2.25–2.45) / 1.42 |
>
> **Hydrate.** The hole window is limited by a water lone-pair level, not by the framework. PBE, PBE+U and HSE06 all keep the hydrate unipolar with both windows ≥ 0.9 eV. Only PBE+U shows a large hydration loss (−54 %), so the pre-registered hydration test is a method split, resolved in favour of the designated arbiter functional.
>
> **Point defects.** The hole edge survives every computed point defect except antisite pairs: ≥ 87 % overall, and ≥ 94 % at K vacancies and water-filled [Cr(CN)6] vacancies. The electron edge at [Cr(CN)6] vacancies is unresolved: killed in PBE+U, with HSE pending and provisionally ≈ 0.7 eV in-cell.
>
> **Testable consequence.** Between the V→Cr charge-transfer onset (≈ 2.0–2.3 eV) and the LMCT and Cr(III) d–d bands (≈ 3.1–3.3 eV), interband absorption occurs only in the V-majority channel. Its magneto-optical sign follows the Néel vector and its magnitude follows the sublattice order. This can be tested on the as-made hydrate across a compensation series.
>
> **Temperature bar.** KVCr orders below the 400 K bar. KV[Mo(CN)6] (model T_C 487 ± 49 K, unipolar, never made at 1:1) and Cr[Mo(CN)6] (bipolar, ≈ 370 K, never made) remain ideal-crystal predictions.

**Do not claim:**
- a pre-registered hydration "PASS";
- −9 % or 2.43 eV as converged numbers;
- "100 % spin-polarized photocarriers at 2.0–5.8 eV";
- a ≥ 400 K or 390–450 K temperature for Cr[Mo(CN)6];
- that the vacancy CB kill is a demonstrated artefact.

---

## 3. Remaining fatal objections (to big news, not to the paper)

1. **Temperature bar.** No compound that exists reaches 400 K.
   - KVCr: 376 K.
   - 1:1 KVMo: never made; model only (487 ± 49 K, floor 377 K); labile, air-sensitive building block.
   - Cr[Mo(CN)6]: ≈ 370 K, never made, bipolar.
2. **No measurement of any LCM-specific property.** The experimental basis is one 1999 hydrated powder (elemental analysis plus PXRD):
   - never reproduced;
   - no single crystal or occupancy refinement;
   - no optical gap, transport or L-edge XMCD;
   - air-sensitive (oxidized V–Cr orders at 115 K).
3. **Compensation is composition-tuned, not protected.**
   - M = w − 3v − u μB/f.u. The measured 0.125 μB/f.u. needs v ≈ 0.04 or w ≈ 0.125.
   - The SOC orbital residue was never computed (estimated 0.02–0.1 μB/f.u.).
   - Antisites keep M = 0 but kill the edges locally, and are excluded only kinetically.
4. **Prior art caps novelty** (referee 3's list was not re-verified here):
   - spin-split compensated ferrimagnets in general (van Leuken & de Groot 1995; Pickett 1998);
   - CrVTiAl, an above-RT compensated ferrimagnetic semiconductor already proposed for spin filtering (PRB 97, 054407 (2018));
   - zero-moment MOKE and TMR in Mn2RuxGa;
   - Schart 2024 (Cr[Cr(CN)6], TDOS↑ ≠ TDOS↓ at M = 0);
   - Guo 2025 (cyanide LCM);
   - Cr(pyz)3 (Nat. Chem. 2026).

   What is new is quantitative: eV windows, a clean 2 eV gap and structurally fixed sublattices. Scoop risk is high.
5. **Excitations are polaronic or excitonic.** The bands are flat (0.4–0.6 eV; hydrate Cr t2g↓ 0.39 eV), V(III) self-traps (−0.27 eV), and the MMCT exciton is bound by about 1 eV. Spin-transport claims are not available.

**Not fatal, but must be fixed before any write-up:**
- the unconverged hydrate EXX loop;
- the "PASSES" wording;
- the pending vacancy electron edge;
- the overstated CrMo temperature and photocarrier range (§5).

---

## 4. Minimal decisive next steps

### 4.1 Computations (cheap, in this order)

**K1. Harvest `KVCr_VAC_w1_hse` against lines fixed now.** Do not cancel it.
- **Disclosure:** these lines were written after the PBE stage and the first EXX block had been seen. They were also set independently by referees 1 and 2 before that.
- **Lines:**
  - lowest opposite-spin empty state ≥ CBM + 0.3 eV after alignment → the F5 vacancy kill is lifted;
  - ≥ 0.72 eV → the 50 % rule is met.
- **Report** the 0.72–0.79 eV reference ambiguity, and the same-spin water-hybridized CBM.
- **When to harvest:** 3–5 EXX blocks suffice. In both reference runs, win_CB and the offset between the water state and the CBM converged within the first few blocks; only the V t2g VBM is slow.
- **If the result lands at 0.65–0.80 eV**, add a matched 4-f.u. pristine HSE (nq 1, k 2×2×2) before calling the 50 % rule.

**K2. Hydrate PBE+U U grid plus PDOS.**
- Four U points, plus projwfc on the PBE+U hydrate and on the PBE stage.
- This names the window-limiting state and formally tests the JUDGE §3 pass conditions.
- Expected: the water level stays put while the d manifolds move with U, so the PBE+U window tracks U_V and U_Cr. That would demonstrate the artefact explanation rather than infer it from band counting.

**K3. Converge the hydrate HSE.**
- Restart until dexx < conv_thr; the matched anhydrous run needed 31 blocks.
- Needed before publication, not for the verdict.
- Until then quote 2.0–2.2 / 2.25–2.45 / 1.42 eV and −9 to −16 %.

**K4. Second water arrangement.**
- Conventional 4-f.u. cell, cubic metric, antiparallel or disordered water dipoles.
- PBE+U relax, then HSE nq 1, compared at matched EXX block with the anhydrous run.
- This turns "one arrangement" into a bounded hole window (expected 1.8–2.8 eV).

**K5. Photocarrier claim.** Rewrite it as in §1.4; no compute is needed. Optional: spin-resolved independent-particle ε2↑(ω) and ε2↓(ω) in HSE.

**K6. Cr[Mo(CN)6].** Only if it is promoted beyond a footnote:
- the Mo[Cr(CN)6] linkage-isomer energy and its windows;
- HSE FM (M fixed at 6) for the exchange ratio;
- T referenced to KVCr.

Drop it if the isomer is downhill or if T stays below 400 K.

### 4.2 Experiments (the only steps that move P materially)

**X1 (cheapest, gating). Reproduce KV[Cr(CN)6]·2H2O by the Holmes–Girolami anaerobic aqueous route** (V(II) triflate + K3[Cr(CN)6]).
- Measure ICP, TGA (water count), PXRD, IR ν(CN) and M_sat(5 K), to confirm |3v − w| ≤ 0.05.
- If the HSE hydrate result holds, no dehydration step is needed.

**X2 (decisive for the LCM label). Sublattice-resolved magneto-optics across S = 0.**
- V and Cr L2,3 XMCD (sublattice signs), plus field-cooled MCD/Faraday at **2.0–3.1 eV** and 300 K.
- Run it across a K/Cs compensation series that spans S = 0.
- LCM signature: the MO signal of the V→Cr MMCT band stays finite as M → 0, and its sign follows the Cr-sublattice direction, not the sign of M.
- **Read first:** Ohkoshi et al., J. Phys. Chem. B 104, 9365 (2000), in full. Its reported "MO signal of the MMCT band ∝ M" in V–Cr films is either a temperature series (harmless) or a composition series. A composition series would be in tension with the prediction, or would reflect the loss of domain selection at small M.

**X3. Spin-resolved, angle-integrated photoemission of the top VB.**
- Field-cooled, air-free samples with glovebox-to-UHV transfer.
- Prediction: V-majority polarization over the top ≈ 2.3 eV (hydrate, HSE), reversing with the Néel vector.

**Route to the bar as written:** crystalline 1:1 KV[Mo(CN)6] with remanence ≥ 400 K (route in `vac/README.md` §5). Nothing in this round changes its odds.

### 4.3 What would move P(big news)

| Outcome | P after |
|---|---|
| K1 lands ≥ 0.72 eV and K3 converges inside 2.25–2.45 eV | ≈ 5 % |
| K1 lands < 0.3 eV: device and spin-filter claims restricted to the hole edge | ≈ 3 % |
| X1 + X2 show a finite, Néel-signed MMCT magneto-optical signal through M = 0 | ≈ 15–25 %. A strong JACS / Nat. Commun. result (an experimentally confirmed room-temperature LCM semiconductor), still below the 400 K bar |
| Crystalline 1:1 KVMo with remanence ≥ 400 K | The only event that meets the bar as written |
| Ohkoshi 2000 turns out to be a composition series with MO ∝ M and no domain explanation, or a scoop | ≤ 2 % |

---

## 5. Corrections required in the 14:35 addenda, the log and §9

| # | Where | Was | Should be |
|---|---|---|---|
| C1 | README addendum item 1; report addendum bullet 1; log 14:25 | "Hydration rule now PASSES with the arbiter functional" | "Method split: kill triggered in PBE+U only (−54 %), not in PBE (≈ 0 %) or HSE06 (−9 to −16 %); existing-sample claim reinstated with that caveat; pass conditions (U grid, PDOS) untested." Add **F9** to report §9.4 as a post-hoc reversal of the 11:30 pre-commitment |
| C2 | same | "gap 2.02, windows 2.43 / 1.42, −9 % / −2 %" | gap ≈ 2.1 (2.0–2.2), windows ≈ 2.35 (2.25–2.45) / 1.42 eV, −12 % (−9 to −16 %) / −2 %; matched-block comparison −15 % |
| C3 | same | "eigenvalues from the last converged EXX block (error plateaued at 1.1e-6 Ry)" | EXX loop not converged (dexx 1.06–1.15e-6 Ry against conv_thr 2.1e-7; VBM still moving 16–17 meV/block; job aborted). win_CB, the other-channel gap and the water level are converged |
| C4 | same; log 14:25 | "The PBE+U −54 % collapse was a PBE+U artefact (water levels mis-placed)" | "PBE+U is the outlier of three functionals; the HSE hole window is water-limited (water 1b1 at VBM − 2.4 eV; Cr t2g↑ at VBM − 2.8 eV); +U misplacing the water level is the likely cause (band counting; PDOS pending)" |
| C5 | README addendum item 2; log 12:05 | Vacancy state "73 % water O–H σ*, the same kind of water level that PBE+U mis-placed"; "PBE underbinds σ*" | The state is a cap O–H / cavity state present with and without the extra water (PBE+U +0.07 / +0.08 eV; PBE +0.91 eV; first EXX block +0.76 eV in-cell). The relax was unconverged at the 45-min cap. The artefact rationale cuts both ways |
| C6 | README addendum item 3; report addendum bullet 3; log 12:35, 14:05 | "Exchange ratio to CrCr 1.6–1.9 → T ≈ 390–450 K (borderline)" | Ratio 1.53–2.04 (base 1.63). T centre ≈ 370 K (KVCr-referenced 369 K; CrCr scaling 367–490 K, which over-predicts KVCr itself). Below the 400 K bar. No HSE FM. Also note the uncomputed Mo[Cr(CN)6] isomer, Mo(III) lability and the Nelson 2007 110 K cluster glass |
| C7 | log 14:55 | "≈ 2.0–5.8 eV … 100 % spin-polarized photocarriers from unpolarized light at M ≈ 0" | "V-majority-only interband absorption at ≈ 2.0–3.1 eV (below the ≈ 3.1 eV LMCT and the ≈ 3.3 eV Cr ⁴T2 bands); excitonic or polaronic; polarization ∝ sublattice order (≈ 0.6 at 300 K) and requires Néel-domain selection; a corollary of the unipolar edges" |
| C8 | README §0 claim 3, §5 item 1, "Not claimed"; report §0 pt 5, §9.0, §9.3.5, §9.7, §9.8 | Existing-sample claim withdrawn; experiments need anhydrous KVCr | Reinstated as a method split (C1). The experiments can use the as-made hydrate, with the caveats of §1.1 |
| C9 | README §6; report §9.9 | Hydrate HSE pair "[PENDING]" | Hydrate harvested (unconverged, aborted at block 14); anhydrous nq 1 converged (31 blocks); `KVCr_VAC_w1_hse` running (first EXX block at 15:35) |

---

## 6. Referee disagreements and how they were resolved

| Issue | Referee 1 (realizability) | Referee 2 (physics) | Referee 3 (impact) | Resolution |
|---|---|---|---|---|
| Hydrate EXX convergence | Adequate (plateau; 16 meV last step) | Not converged; −15 % at matched block; ≈ 2.26 eV | Not converged; −9 to −16 % | Reparse confirms referees 2 and 3. Best estimate ≈ 2.35 eV, −12 % (−9 to −16 %) |
| What caps the hydrate hole window | Unidentified; asks for PDOS | Water 1b1 (band counting) | Probably water 1b1 | Confirmed by band counting: water at −4.45 eV, Cr t2g↑ at −4.84 eV relative to the CBM. PDOS still needed (K2) |
| Status under the rule | Method split; F9 | Method split | Method split; HSE preferred on physical grounds | Unanimous: method split, not PASS |
| Vacancy HSE outcome | P(lift) ≈ 0.7 | Lift 0.6–0.7; ≥ 50 % a coin flip | P(lift) ≈ 0.4 (cap/cavity state) | First EXX block puts the state at CBM + 0.76 eV → P(lift) ≈ 0.85; P(≥ 50 %) ≈ 0.3 |
| CrMo T | 369 K (KVCr) / 391 K (CrCr) | ≈ 370 K | 350–420 K | Centre ≈ 370 K, below the bar |
| Clean photo-window | 2.0–3.2 eV | 2.2–3.3 eV | 2.0–3.2 eV | ≈ 2.0–3.1 eV (the LMCT at 3.1 eV, channel unassigned, also bounds it) |
| P(big news) | 4 % | 4 % | 5 % | ≈ 4 % |

---

## 7. Verification done for this synthesis

- **Fetched read-only from the Modal volume `magdisc-data`:**
  - `jobs/lcm/pba/hse/KVCr_hse_LCM_nq1/hse_lcm.out`;
  - `jobs/lcm/pba/hse/KVCr_VAC_w1_hse/{hse_lcm.out, status.json}` (state "running", 15:35 PDT).
- **Reused a referee's local copy** of `jobs/lcm/pba/hse/KVCr2H2O_hse_LCM/hse_lcm.out` (ends in "Abort"; input conv_thr 2.1d-07, nqx 1, ecutfock 180).
- **Parsers** (session scratchpad, temporary): `scratchpad/ref/parse.py` (block windows), plus `scratchpad/rr/edges.py` and `scratchpad/rr/track.py` (band-resolved levels relative to the CBM). Full path prefix: `/private/tmp/claude-501/-Users-geby-chemistry/0128b63c-3ab4-49af-9683-b3bd1a9fe21d/scratchpad/`.
- **Reproduced:**
  - every hydrate and anhydrous block in §1.1;
  - the dexx history (hydrate 1.06–1.15e-6 Ry at the end; anhydrous 1.0e-7 at convergence);
  - the water, Cr t2g↑ and VBM level positions;
  - the vacancy PBE stage and first EXX block in §1.2.
- **Taken from the referees without re-checking here:**
  - the CrMo PBE+U and HSE numbers and ratios (referees 1 and 2 independently agree);
  - the KVCr_VAC_w1 PBE+U band characters (referee 1);
  - the prior-art citations (referee 3);
  - the Johansson 2016 ²E timescale and the Ohkoshi 2000 statement (referee 1);
  - the recalled Cr(III) d–d energies.
- **Untouched:** no job was submitted or cancelled, and no lane file other than this one was written.
