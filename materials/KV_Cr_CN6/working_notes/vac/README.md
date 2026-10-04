# PBA vacancies and building-block lability: KV[Cr(CN)6] and 1:1 KV[Mo(CN)6]

Track L (session 66190d), 2026-10-03. Folder `tracks/lcm/pba/vac/`. This lane covers three things:
- the [M(CN)6]-vacancy cells of defect Gate 1 (`../defects/`);
- how labile the [Cr(CN)6]³⁻ and [Mo(CN)6]³⁻ building blocks are, and what defect concentrations that implies;
- a synthesis plan for crystalline 1:1 KV[Mo(CN)6].

The computed numbers below are in `summary.json`; literature values carry their sources in §3.

**Reproduce**
- `python tracks/lcm/pba/vac/collect_vac.py --fetch`: vacancy cells, original and `_r`. It reuses `../defects/collect.py analyse()` unchanged and reads fixed paths only.
- `python tracks/lcm/pba/defects/collect.py --fetch`: all Gate-1 cells, including the K-vacancy polaron variants.
- `python tracks/lcm/pba/vac/summarize.py`: writes `summary.json`.

**Spin convention.** As in the lane: C-bound M is +, N-bound V is −. "V-maj" is the spin channel of the V moment, which carries both host edges; "M-maj" is the opposite channel.

## 0. Bottom line

| question | answer |
|---|---|
| M_cell per [M(CN)6] vacancy | **−3.000 μB** in both KVCr and KVMo (counting value; V-sublattice sign). The capped V stays V(II) t2g³ (N_d 3.50–3.51, m_d −2.8). |
| gap (25 % ordered vacancies, cell) | 1.43 eV (KVCr, pristine 1.95); 1.04 eV (KVMo, pristine 1.30) |
| edge spins in the vacancy cell | both edges V-majority (unipolar); in-cell windows 2.07 / 0.45 eV (KVCr) and 0.98 / 1.08 eV (KVMo) |
| VB window vs host (aligned) | **survives**: 1.93 eV (95 %) and 0.91 eV (94 %) |
| CB window vs host (aligned) | KVCr: **lost locally** (an empty water-void state sits at the host CBM, +0.07 eV; K 3p alignment −0.25 eV). KVMo: **0.89 eV** (47 %; K 3p 0.45 eV), above the 0.3-eV bar but a nominal fail of the 50 % rule |
| defect states | Capped-V t2g lie inside the VB (same spin). The M t2g↓ CB is pulled 0.2–0.4 eV down at the vacancy (same spin). One empty, weakly spin-polarised void state (cap O–H, H_w 0.67), about 2.0–2.2 eV above the host VBM in both compounds, sets the opposite-spin CB bottom. |
| [Cr(CN)6]³⁻ lability | Inert in solution (k_aq(Cr³⁺) 2.4e-6 s⁻¹). Linkage flips in Fe–Cr PBAs run at 255–321 K and under 1 GPa, but only because LS Fe(II) is C-philic. In V–Cr the as-made linkage is thermodynamic (flip +0.18 eV, antisite pair +0.93 eV), so V–Cr aging is oxidation and dehydration, not isomerisation. |
| [Mo(CN)6]³⁻ lability | **Not inert.** It exists only under cyanide limitation (otherwise [Mo(CN)7]⁴⁻), decomposes in protic solvents, is photolabile, and Mo(III) substitutes about 10⁵× faster than Cr(III). The built framework resists flips as well as Cr (+0.18 eV), but antisite-like defects (antisite pair +0.57 eV) are a synthesis-stage risk. |
| realistic content, KVCr | v ≈ 0.04 (0–0.08), antisites < 10⁻⁴, flips ≤ 10⁻³ → **M(0) 0–0.24 μB/f.u.** (measured 0.125) |
| realistic content, 1:1 KVMo | v ≈ 0.06 (0.03–0.15), antisite-like 10⁻³–10⁻², flips ≤ 2 × 10⁻³ → **M(0) 0.1–0.45 μB/f.u.**, T_C −2 to −10 % |
| best synthesis | Aprotic (DMF), cyanide-limited [Mo(CN)6]³⁻ (Li or NEt4 salt, **no cryptand**), naked K⁺ (or Cs⁺) in 5–20× excess, [V(DMF)6]²⁺ "protected" V(II), slow diffusion or slow addition of V into Mo + K at 1–5 mM. Glovebox, dark, ≤ 320 K. Calibrate on KVCr first. |
| main failure modes | colloid or amorphous product; K template failure (V-rich, vacancies); [Mo(CN)6]³⁻ decomposition → [Mo(CN)7]⁴⁻ (the 110-K phase) or C-bound V; oxidation (minutes in air); > 340 K decomposition |

**Verdict for the lead**
- Vacancies do not touch the compensation mechanism beyond the counting rule. They keep the edges V-majority in the cell, and they leave the VB edge intact in both compounds.
- For KV[Mo(CN)6] the CB edge also survives (0.45–0.9 eV).
- For KV[Cr(CN)6], the experimental anchor, the electron edge is shunted at every vacancy by an empty void state. The VB edge, and with it the photoemission test, is unaffected; the MCD test is largely unaffected.
- That KVCr kill hangs on an empty-void model and on PBE+U, so `KVCr_VAC_w1` (zeolitic water at the vacancy centre) is running to test it.

## 1. The vacancy jobs: what happened and what was resubmitted

**Original jobs** (`../defects/manifest.json`). Both were submitted at 05:18 PDT on run_cpu64 with RELAX_MAXSEC = 5 h, as 62-atom R3 cells (3 operations, 2 k-points).

**Diagnosis at 09:40 PDT.** Neither job had died or got stuck in the SCF. Both were in slow BFGS descents (`results/` peeks; trajectories parsed from the relax `.out`).

| job | state at 09:40 | BFGS SCF steps | relaxation so far | last steps | max force component |
|---|---|---|---|---|---|
| KVCr_VAC | preempted once at about 07:10 (KeyboardInterrupt in stdout) and restarted at 07:11 from block 16 | 16 + 33 | −1.39 eV, then −0.89 eV | −5 to −20 meV/step | 0.004 Ry/bohr |
| KVMo_VAC | 4.3 h into the 5-h cap | 55 | −2.36 eV | −3 to −5 meV/step | 0.002–0.004 Ry/bohr |

**What is relaxing**
- The six V–OH2 caps rotate: H moves 0.5–0.8 Å, and V–O lengthens from 2.13 to 2.22–2.26 Å.
- The single remaining K slides 1.0–1.2 Å along the C3 axis onto a face of three N (K–N 2.93–2.97 Å).
- No water dissociates (O–H 0.97 Å), and M stays at −3.00 throughout.
- These are floppy modes with a flat landscape, which BFGS walks down slowly. At 0.004 Ry/bohr the electronic structure is already converged for our purpose.

**When the originals would finish**
- KVMo_VAC reaches its relax cap at about 10:19, then runs scf + nscf + projwfc (about 50 min).
- KVCr_VAC reaches its cap at about 12:10, so its full result lands near 13:00.
- Both land after this lane's ~60-min window. They were left running; nothing was cancelled.

**Resubmitted with faster settings (suffix `_r`, 09:42 PDT, `submit_vac_r.py`)**

| job | call id | positions | settings |
|---|---|---|---|
| `lcm/pba/def/KVCr_VAC_r` | fc-01M41AA0620YQXGFYMT0ZYDQ1B | last complete ATOMIC_POSITIONS block of the running relax (33 blocks; 49 steps in total) | SKIP_RELAX=1, otherwise identical: 110 Ry, U, scf 2×2×2s, nscf 4×4×4, projwfc, run_cpu64 |
| `lcm/pba/def/KVMo_VAC_r` | fc-01M41AA1T9Z2PQ21W2DPCMRSMY | same (55 blocks) | same |

- Both jobs use the same `pba_defect.py` (md5 in `manifest.json`) with the same 110-Ry / U / k protocol as every Gate-1 cell, so `../defects/collect.py analyse()` aligns them to the pristine 4 f.u. cell unchanged.
- The energies of the `_r` cells are not final, because the relax is still descending. Only the electronic quantities are used from them.

## 2. Vacancy cells: moments, gap, windows, defect states

Source: the `_r` cells. These are KV4M3(CN)18(H2O)6, i.e. one [M(CN)6] vacancy per 4 f.u. (v = 0.25, ordered), with six V–OH2 caps and three K removed. Alignment to the pristine 4 f.u. cell uses the far-field CN 3σ level (36 C/N atoms more than 4.5 Å from the vacancy); the K 3p cross-check is given in brackets.

| | KVCr_VAC_r | KVMo_VAC_r |
|---|---|---|
| M_cell (μB/cell), expected −3 | **−3.000** | **−3.000** |
| \|M\| (μB/cell); pristine 4 f.u. | 24.40; 28.29 | 23.16; 26.59 |
| capped V(II): m_d / N_d (pristine V) | −2.82 / 3.51 (−2.70 / 3.51) | −2.77 / 3.50 (−2.61 / 3.48) |
| uncapped V: m_d | −2.70 | −2.59 |
| M(C): m_d | +2.95 (pristine +2.93) | +2.47 (pristine +2.45) |
| cell gap (eV); pristine | **1.43**; 1.95 | **1.04**; 1.30 |
| edge spins in the cell | VBM V-maj, CBM V-maj | VBM V-maj, CBM V-maj |
| in-cell windows VB / CB (eV) | **2.07 / 0.45** | **0.98 / 1.08** |
| alignment shift: CN 3σ (K 3p) (eV) | −0.95 (−0.63) | −1.03 (−0.59) |
| opposite-spin VB top vs host VBM | −1.93 | −0.91 |
| opposite-spin CB bottom vs host CBM | **+0.07 (K 3p: −0.25)** | **+0.89 (K 3p: +0.45)** |
| aligned windows VB / CB (eV), % of pristine | 1.93 (95 %) / 0.07 (6 %) | 0.91 (94 %) / 0.89 (47 %) |
| Gate-1 flags (`analyse()`) | opposite-spin state near an edge: **yes** (CB); window cut: **yes** (CB) | opposite-spin state near an edge: no; window cut: **yes** (CB 47 % < 50 %) |

**Moments**
- Each vacancy removes exactly one M³⁺ moment, so **M = 3 μB per vacancy**, along the V sublattice.
- The capped V stays V(II) t2g³. N_d is unchanged, and the d moment is slightly larger (more ionic O donors). There is no charge transfer and no fractional moment.

**Where the defect states sit** (relative to the host VBM, CN 3σ alignment; band table in `results/jobs/.../<id>.json`)
- **Occupied V-maj states, harmless.**
  - The t2g of the three capped V lie inside the valence band: −0.91 to −0.44 eV (Cr) and −0.88 to −0.35 eV (Mo).
  - The one uncapped V sets the cell VBM at +0.04 to +0.14 eV (Cr) and −0.07 to +0.07 eV (Mo). Same spin as the host VBM.
  - The opposite-spin (M-maj) occupied t2g do not move: the VB windows keep 95 % and 94 %.
- **Empty V-maj states, harmless (same spin).** At the vacancy the M t2g↓ conduction band is pulled down by 0.38 eV (Cr), where it mixes 27 % with the void state below, and by 0.19 eV (Mo). This is what shrinks the cell gap.
- **The empty void ("cavity") state.** This one sets the opposite-spin CB bottom in both cells.
  - It is one band localised in the empty vacancy: region weight 0.95, with H_w 0.67, O_w 0.22 and V 0.06. It is built from the 12 inward-pointing cap O–H bonds and bound by the positive potential of the missing [M(CN)6]³⁻.
  - It is not a d state, and it sits at nearly the same energy above the V t2g VBM in both compounds:

| compound | M-maj copy (vs host VBM) | V-maj partner |
|---|---|---|
| KVCr | +2.02 eV | at +2.94 eV, partly mixed into the CB edge |
| KVMo | +2.19 eV | +2.11 eV, nearly spin-degenerate |

  - **KVCr:** the gap (1.95 eV) is larger than this energy, so the void state lands **at the host CBM**: +0.07 eV (K 3p alignment: −0.25 eV). That is inside the 0.3-eV kill zone for either alignment.
  - **KVMo:** the CB is the lower Mo t2g↓ band (VBM + 1.30 eV), so the void state stays **0.89 eV above the CBM** (K 3p: 0.45 eV). The CB window falls to 47 % of pristine, a nominal fail of the 50 % rule. It still clears the 0.3-eV absolute bar under either alignment.

**Caveats on the void state**
1. **Model.** The void is *empty*. Real Prussian-blue vacancies hold zeolitic water at the vacancy centre [recalled: Herren et al. 1980], and the aprotic route of §5 would cap V with MeCN or DMF instead of water. Either should push a void-bound state up. Test job `KVCr_VAC_w1` adds one central H2O (submitted 10:11, §6).
2. **Method.** PBE+U corrects only the d states. The position of an s/σ*-like void state relative to the U-shifted d edges is uncertain by several tenths of an eV. HSE moves the KVCr CB window from 1.15 to 1.57 eV in the pristine cell.
3. **Alignment.** The CN 3σ and K 3p references disagree by 0.32 eV (Cr) and 0.44 eV (Mo). Both probe atoms lie 4.5–9 Å from a +3 defect in a cell with 25 % vacancies. Take ±0.4 eV as the alignment error. This does not affect any VB conclusion, and it does not change the KVCr CB verdict.

**Relax status of these geometries.** At the positions used, the max force is 0.004 Ry/bohr and the energy is still descending at 3–20 meV per BFGS step; the floppy modes are cap-H rotation and K sliding. The originals (`KVCr_VAC`, `KVMo_VAC`) will give the final relaxed numbers (§6). Expect electronic changes of ≤ 0.1 eV, since the d levels are set by V–N/V–O and M–C bonds that no longer move.

**Harvest of the K-vacancy polaron variants (KVMo, Gate 1)**

| variant | E − E(KVAC delocalised hole) | state |
|---|---|---|
| polM: hole forced into the Mo-majority channel (Mo(IV)) | **+0.62 eV** (last relax step) / +0.73 eV (fresh SCF) | done 10:15. The relax stopped unconverged at 3 h (total force 0.05 Ry/bohr), so this is an upper bound. M = −1; the forced Mo-maj hole state lies 0.56 eV below the host CBM |
| polV: hole self-trapped on one V (V(III), M = +1) | **−0.27 eV** | relax still running (40 steps, F 0.009 Ry/bohr) |

- K-vacancy holes in KVMo self-trap as V(III) small polarons, which are V-majority, so the edge spin is kept.
- An opposite-spin Mo(IV) hole is not competitive.
- p-type KVMo would therefore be spin-locked but polaronic.

## 3. Building-block lability: literature and DFT reasoning

Read-status flags as in `literature/lcm_pba_facts.md`: [abstract], [snippet] (search-engine text only), [secondary], [recalled] (not re-read today), [derived] (our arithmetic), [inference].

### 3.1 [Cr(CN)6]³⁻: inert in solution, but the bridge can still flip in the solid if a flip runs downhill

**Inert as a complex**
- Cr(III) is the textbook inert ion: [Cr(H2O)6]³⁺ exchanges water at k = 2.4 × 10⁻⁶ s⁻¹ (Nelson et al., Inorg. Chem. 46, 10093 (2007) [abstract]).
- [Cr(CN)6]³⁻ aquation is slow and acid-assisted [snippet], and it photo-aquates [snippet].
- In flow-battery electrolytes, cyanide exchange of K3[Cr(CN)6] is the degradation channel. It is suppressed by added cyanide and base (Angew. Chem. 2025, doi 10.1002/anie.202507119 [abstract]).

**Linkage flips need a driving force**
- **Fe–Cr precedents.**
  - The as-made Cr–C≡N–Fe(II) linkage (red KFe[Cr(CN)6]) flips to Cr–N≡C–Fe(II) (green) on heating: Shriver, Shriver & Anderson, Inorg. Chem. 4, 725 (1965) [snippet].
  - In Fe(II)–Cr(III) nanocrystals the flip runs between 255 and 321 K, faster for smaller particles: Inorg. Chem. 2013, doi 10.1021/ic302764k [abstract].
  - In single crystals of K0.4Fe4[Cr(CN)6]2.8, 0–1.2 GPa turns the cyanide reversibly: JACS 127, 4580 (2005) [abstract]. After about 1.2 GPa a non-reverting phase forms: JACS 130, 15519 (2008) [snippet].
  - In CsFe[Cr(CN)6], the spin crossover at about 207 K occurs **without** isomerisation (J. Phys. Chem. B, doi 10.1021/jp068885+ [snippet]).
- **The common driver** is a C-philic partner. Low-spin Fe(II) gains a large ligand-field stabilisation on the C end, and Fe(II), Co(II) and BPh3 Lewis acids labilise [Cr(CN)6]³⁻ by linkage isomerism (Avendaño, Karadaş, Hilfiger, Shatruk & Dunbar, Inorg. Chem. 2010, doi 10.1021/ic901681e [abstract]).
  - So "[Cr(CN)6]³⁻ is inert" does **not** protect a linkage against flipping. In the solid, a downhill flip happens at room temperature.
- **V is N-philic.**
  - [V(III)(CN)6]³⁻ placed against Fe(II) or Ni(II) isomerises completely so that V ends up N-bound. With Co the isomerisation is partial; with Mn it does not occur (Nelson & Miller, Inorg. Chem. 47, 2526 (2008) [abstract]).
  - For V(II)/Cr(III), Gate 1 gives a CN flip at **+0.18 eV** and the full V/Cr antisite pair at **+0.93 eV** (§3.4). The as-assembled V–N≡C–Cr linkage is the thermodynamic one, so neither heat nor pressure can drive V–Cr towards antisites the way they drive Fe–Cr.

**What "aging" of V–Cr PBAs actually is**
- The reported aging is redox and loss of solvent, not isomerisation.
- **Oxidation.**
  - V–Cr PBAs are air-sensitive (V(II) → V(III) → vanadyl). The fully oxidised (VO)[Cr(CN)6]2/3 orders at 115 K (`literature/lcm_pba_facts.md` §1).
- **Dehydration.**
  - KVCr drops from 376 to 365 K after heat treatment, and V[Cr(CN)6]2/3 from 330 to 320 K after heating to 350 K (Holmes & Girolami 1999 via the review [secondary]).
- **Irreproducible T_C** in the early aqueous routes was cured by 1–4 % V(III) "catalyst". It gives "better organized" solids with reproducible T_C and an M_sat matching the V/Cr ratio (Garde, Villain & Verdaguer, JACS 124, 10531 (2002) [abstract]).
- No report of V–Cr linkage isomerism was found.

### 3.2 [Mo(CN)6]³⁻: not an inert building block

**Exists only under cyanide limitation**
- Mo(III) prefers seven-coordination: [Mo(CN)7]⁴⁻. Octahedral [Mo(CN)6]³⁻ was first made in 2002 by limiting the cyanide: Mo(OTf)3 + LiCN in DMF gives Li3[Mo(CN)6]·6DMF, alongside the cyano-bridged [Mo2(CN)11]⁵⁻ (Beauvais & Long, JACS 124, 2110 (2002) [abstract]).
- Any free CN⁻ during assembly therefore converts the building block to [Mo(CN)7]⁴⁻. The V(II)–[Mo(CN)7] network orders at only 110 K (Tomono et al., Inorg. Chem. 49, 1298 (2010) [title]).

**Decomposition and photolability**
- [Mo(CN)6]³⁻ decomposes in water, methanol and ethanol. The V–Mo product decomposes in air within minutes and irreversibly above 340 K (Magott et al., Adv. Sci. 13, e11285 (2026) [part, via PMC summary]).
- The Mo(III)–CN dative bond is photolabile: K4[Mo(CN)7]·2H2O switches reversibly between 6- and 7-coordination under visible light, inside the crystal (Nat. Commun. 2025, doi 10.1038/s41467-025-63523-x [abstract]).

**Substitution rate**
- Mo(III) substitutes associatively and fast. [Mo(H2O)6]³⁺ + NCS⁻ has k = 0.27 M⁻¹ s⁻¹ at 25 °C (Sasaki & Sykes [snippet]); Cr(III) anation is of order 10⁻⁶ M⁻¹ s⁻¹ [recalled]. Mo(III) is therefore about 10⁵× more labile than Cr(III) [derived/recalled].
- The reasons are the larger 4d radius and the low-lying seven-coordinate transition state.

**Evidence from the only V–Mo PBA**
- In V1.37[Mo(CN)6], ν(CN) goes from 2088 cm⁻¹ in the precursor to bands at 2114, 2098 and **2047** cm⁻¹. These are assigned to Mo–CN–V bridges plus non-bridging cyanide.
- EXAFS shows V N-coordinated. The bridge direction is not separately established [part].
- **Synthesis observations.**
  - Solution reactions in MeCN gave only colloids.
  - The amorphous product came from mechanochemistry with [K(crypt-222)]3[Mo(CN)6]·2MeCN.
  - K(crypt)⁺ is far too large for the PBA cavity, which explains the V-rich, 27 %-vacancy composition [inference].

**Consequence for KVMo [inference]**
- Inside a complete framework, all six CN of Mo are bridged to V and there is no room for a seventh ligand, so the built framework should be as robust as the Cr one against flips: E_flip is +0.18 eV in both.
- The danger is **during assembly**:
  - free cyanide → [Mo(CN)7]⁴⁻ inclusion, and C-first capture by V(II). C-bound V(II) is well known, e.g. [V(CN)6]⁴⁻ and Cs2Mn[V(CN)6] (Entley & Girolami 1995 [secondary]);
  - protic impurities → decomposition and vacancies;
  - light;
  - surfaces and vacancy sites, where Mo is under-coordinated.

### 3.3 V(II): the N-bound partner

- [V(H2O)6]²⁺ exchanges water at 87 s⁻¹ by associative interchange [snippet; Helm & Merbach, Chem. Rev. 105, 1923 (2005), recalled]. That is slow for an M(II) (10⁴–10⁹ s⁻¹ for most), but V–N≡C bond formation is still reversible on a seconds scale.
- V(III) exchanges faster (about 5 × 10² s⁻¹ [recalled]). V(II)/V(III) electron exchange therefore opens a redox-catalysed substitution path, which is a plausible reading of Garde's catalytic-V(III) effect [inference]. The price is a few % V(III) (w > 0).
- V(II) solvates usable in aprotic media include [V(NCMe)6]²⁺ (Magott) and V(OTf)2 (Holmes & Girolami, aqueous).

### 3.4 DFT energetics (Gate 1, PBE+U 110 Ry, 4 f.u. cells) and equilibrium populations

| defect | KVCr (eV) | KVMo (eV) | exp(−E/kT), 300 K | exp(−E/kT), 340 K | effect on M / windows (Gate 1) |
|---|---|---|---|---|---|
| CN flip (one V–C + one M–N bond) | +0.177 | +0.182 | 1.1e-3 / 8.8e-4 | 2.4e-3 / 2.0e-3 | M = 0; windows 98/85 % (Cr), 93/97 % (Mo): harmless |
| V/M antisite pair (non-bonded) | +0.927 | +0.565 | 2.7e-16 / 3.2e-10 | 1.8e-14 / 4.2e-9 | M = 0; VB window cut to 0.13 (Cr) / 0.25 eV (Mo); local opposite-spin levels (§4) |
| full linkage isomer KCr[V(CN)6] (per f.u.) | +0.44 to +0.81 (base +0.63) | not computed | 3e-11 (base) | | isomer CB window 0.01–0.26 eV |
| K-vacancy hole forced into the Mo-majority channel (polM) vs delocalised V-majority hole (KVAC) | – | +0.62 (last relax step) / +0.73 (fresh SCF); relax unconverged, so an upper bound | | | an opposite-spin hole is not competitive |
| K-vacancy hole self-trapped on one V (polV, V(III) small polaron) vs KVAC | – | −0.27 (relax still running) | | | the hole localises but stays V-majority (M = +1) |

**Reading.**
- The equilibrium antisite content is zero for practical purposes, and the flip content is at most about 10⁻³ per CN, which is harmless anyway.
- **The real antisite and flip content is therefore kinetic.** It is set by whether the building block survives assembly intact.
  - [Cr(CN)6]³⁻ does, in the dark, over minutes to hours.
  - [Mo(CN)6]³⁻ does only if the medium is aprotic, cyanide-limited and dark.
- The DFT sign matches experiment: V is N-philic (Nelson & Miller 2008). The KVMo antisite pair costs about 40 % less than KVCr's.

## 4. Realistic defect concentrations and their effect on M(0) and the windows

**Counting rule** (lane sign: C-bound M +, V −): **M = w − 3v − u** μB/f.u.
- v is the [M(CN)6] vacancy fraction, w the V(III) fraction, u the Mo(IV) fraction.
- CN flips and antisites keep M = 0 (both ions are S = 3/2).
- Charge balance: K_x = 1 − 3v − w + u. A K-rich synthesis therefore pushes v and w down together.
- Formula-unit density 3.3 × 10²¹ cm⁻³ (a ≈ 10.6–10.7 Å), so a 10⁻³ fraction per f.u. is 3 × 10¹⁸ cm⁻³.

| quantity | KVCr, best existing sample type (aqueous, K-rich; Holmes–Girolami) | KVCr, aprotic NMF-type route (not yet done) | 1:1 KVMo, first generation (route A, §5) | 1:1 KVMo, optimised |
|---|---|---|---|---|
| v ([M(CN)6] vacancies) | 0–0.08, best estimate 0.04 | ≤ 0.01 (Cr–Cr NMF gave 0.004) | 0.03–0.15, estimate 0.06 | 0.01–0.05 |
| w (V(III)), u (Mo(IV)) | w 0–0.12 (M_sat allows either v = 0.04 or w = 0.12) | ≤ 0.02 | ≤ 0.05 each (oxidation, catalyst) | ≤ 0.02 |
| full antisite pairs per f.u. | < 10⁻⁴ (inert block, uphill by 0.93 eV) | < 10⁻⁴ | 10⁻³–10⁻² antisite-like (C-bound V, N-bound / 7-coordinate Mo from decomposed blocks) | ≤ 10⁻³ |
| CN flips per CN | ≤ 10⁻³ (equilibrium bound at 300–340 K) | ≤ 10⁻³ | ≤ 2 × 10⁻³ | ≤ 10⁻³ |
| **M(0) = \|3v + u − w\|** (μB/f.u.) | 0–0.24; **measured 0.125** | ≤ 0.03 | **0.1–0.45** | 0.03–0.15 |
| T_C change from v (tc model, D(p)) | −0 to −5 % | ≈ 0 | −2 to −10 % (C: 535 → 480–525 K) | −1 to −3 % |

**Where these numbers come from**
- **KVCr.**
  - Holmes & Girolami's z = 1.00 is elemental-analysis precision. Its M_sat of 0.125 μB/f.u. needs |3v − w| = 0.125 (`literature/lcm_pba_facts.md` §3.1).
  - The best alkali V–Cr analogues have z = 0.92–0.95. Antisites are bounded by inertness plus thermodynamics (§3.4).
- **KVMo.**
  - The only data point is v = 0.27 with a template that cannot enter the framework.
  - K/Cs-templated V–Cr samples reach v = 0.05–0.08, and the NMF route reaches 0.004 for Cr–Cr.
  - The antisite-like range reflects the 10⁵× higher lability of Mo(III) and its photo- and protic decomposition. It is an [inference] bracket, not a measurement.

**Consequences for the windows at these concentrations**
- **VB (hole) edge: safe in both.**
  - Even at v = 0.25 the opposite-spin VB top stays 1.93 eV (Cr) and 0.91 eV (Mo) below the host VBM (95 % / 94 %).
  - Vacancy-related occupied states are V-majority and lie inside the VB.
  - V(III) holes (w, K deficiency) are V-majority, delocalised or self-trapped (KVAC: windows 89–108 %; polV −0.27 eV).
  - Spin-resolved photoemission of the top VB (`literature/lcm_pba_facts.md` §5B) is therefore robust to all realistic defects. Antisites at 10⁻³ add opposite-spin occupied levels at VBM −0.13 eV (Cr) / −0.25 eV (Mo), at a density of 10⁻³ per f.u. (about 3 × 10¹⁸ cm⁻³).
- **CB (electron) edge in KVMo: survives.** Vacancies put their void state 0.45–0.9 eV above the CBM. Antisites put Mo(N) empty levels at CBM +0.18 to +0.6 eV. Electrons injected into KVMo stay V-majority to within these margins.
- **CB edge in KVCr: compromised at vacancies.**
  - Each vacancy carries an empty, weakly spin-polarised void level at the host CBM (−0.25 to +0.07 eV). At v ≈ 0.04 that is about 10²⁰ cm⁻³ electron traps with mixed spin. Antisites add Cr(N) opposite-spin empty levels 0.02–0.06 eV below the CBM.
  - Unless the zeolitic-water test lifts the void state, the KVCr spin-filter/tunnelling use of the CB window (§5C of the fact sheet) needs low-vacancy material. The VB-based and magneto-optical tests are unaffected.
- **Compensation.** M(0) is set by v, w and u, not by antisites or flips. Realistic samples sit at 0.03–0.45 μB/f.u. That is small enough for the LCM claim (the bands are spin-split by the sublattice inequivalence, not by M). It is also what lets field-cooling select one Néel domain.
  - The composition zero (3v + u = w) can be crossed by oxidation in air. Each V(III) shifts M by +1 per f.u. fraction and each Mo(IV) by −1, so sealed handling is part of the measurement, not just of the synthesis.

## 5. The most realistic route to crystalline 1:1 KV[Mo(CN)6]

All of this section is [inference] built on the precedents cited in §3. No step has been tried for V–Mo.

**Design rules**
1. **Aprotic and cyanide-limited**, so that [Mo(CN)6]³⁻ survives.
2. **A small, "naked" alkali in large excess**, so that K⁺ templates the tetrahedral cavities and pushes v → 0. Each vacancy needs three K⁺ out (§2).
3. **Low supersaturation with reversible V–N bond formation**, which gives crystals rather than colloids.
4. **Dark, under 1 ppm O2/H2O, never above about 320 K.**

### 5.1 Route A: K-templated slow diffusion in DMF (recommended first attempt)

**Mo source**
- Li3[Mo(CN)6]·6DMF, made by the Beauvais–Long cyanide-limited route: Mo(OTf)3 + LiCN in DMF, strictly ≤ 6 CN per Mo.
- Option: metathesis to (NEt4)3 or (PPN)3[Mo(CN)6], so that Li⁺ cannot compete for the cavities and the salt dissolves in MeCN/DMF.
- **Do not** use [K(crypt-222)]⁺. It is the cation that gave Magott's V-rich, 27 %-vacancy solid.

**V source: "protected" V(II)**
- [V(NCMe)6](BF4)2 or anhydrous V(OTf)2, dissolved in DMF, gives [V(DMF)6]²⁺.
- The solvent itself is the protecting ligand. O-bound DMF is displaced by the N end of cyanide but competes with it, and it is too large to be built into an intact cavity.
- A milder moderator can be screened: a few equivalents of MeCN or a weakly chelating diamine such as tmeda. This is the aprotic analogue of the citrate-moderated growth used for low-vacancy battery PBAs (e.g. Materials 18, 3174 (2025) [abstract]; general practice [recalled]).
- Pyridine-type ligands strong enough to stay on V are to be avoided.

**Alkali source**
- KOTf, KPF6 or K[B(C6F5)4] at 5–20 equivalents per Mo (K-rich).
- CsOTf is the fallback template. Cs⁺ is the strongest vacancy suppressor in PBAs; the Cs V–Cr phase reached z = 0.92–0.94. CsV[Mo(CN)6] serves the LCM claim equally well.

**Assembly**
- H-tube or three-layer diffusion: V(II) in DMF / pure DMF buffer / [Mo(CN)6]³⁻ + K⁺ in DMF.
- 1–5 mM, room temperature or 4 °C, 1–3 weeks, in a glovebox, in the dark or under red light.
- Powder variant: syringe-pump addition of V(II) **into** the Mo + K solution over 24–48 h. The excess hexacyanometallate and alkali favour the A-rich, low-vacancy phase, as in aqueous K2Mn[Fe(CN)6]-type syntheses [recalled]. Then 1–3 days of ripening in the mother liquor.
- Optionally add 1–2 % V(III) (Garde's catalyst) if the product comes out poorly ordered. That buys order at the price of w ≈ 0.01–0.02 (M shift +0.01–0.02 μB/f.u.).

**Isolation**
- Wash with dry DMF, then MeCN; dry in vacuum at ≤ 300 K.
- Store and measure sealed (glass coverslip, glovebox-to-UHV transfer). The parent decomposes in air within minutes and above 340 K.

**Calibrate on Cr first.** Run the same protocol on KV[Cr(CN)6]. Cr is inert, so whatever vacancy content and crystallinity appear there are set by the template and the growth rate alone. Then switch to Mo.

### 5.2 Route B: electrodeposited films (for the MCD, photoemission and tunnelling tests)

- Reduce V(III) to V(II) at an FTO or metal electrode in an aprotic electrolyte (DMF or MeCN, KPF6) containing [Mo(CN)6]³⁻ and K⁺. This is the aprotic version of the V–Cr film route (Mizuno/Ohkoshi 2000; Johansson 2016).
- It gives thin, continuous films directly on the substrate that the experiments in `literature/lcm_pba_facts.md` §5 need.
- Risk: the reduction window. The [Mo(CN)6]³⁻/⁴⁻ potential has not been checked against V(III)/V(II) in the chosen solvent; measure it by cyclic voltammetry first.

### 5.3 Route C: fallback

Mechanochemistry as Magott did, but with a K⁺ or Cs⁺ salt of [Mo(CN)6]³⁻ that has no cryptand, followed by solvent-vapour annealing (DMF vapour, room temperature). This probably still gives poor crystallinity. Its advantage is that it is the only route already shown to keep Mo(III) intact.

### 5.4 Failure modes, ranked by likelihood

| # | failure | signature | consequence | mitigation |
|---|---|---|---|---|
| 1 | Colloidal or amorphous product (nucleation too fast; seen by Magott in MeCN) | broad PXRD, gel | no structure, so no window test on a real lattice | DMF, ≤ 5 mM, diffusion, moderators, ripening |
| 2 | Template failure: K⁺ not incorporated (crypt, Li⁺ competition, solvent in the cavities) | V-rich ICP (V:Mo > 1), K:Mo < 1 | v = 0.1–0.3, giving M(0) = 0.3–0.9 μB/f.u. and T_C −6 to −20 % | naked K⁺/Cs⁺ in excess; metathesise away Li⁺; reverse addition |
| 3 | [Mo(CN)6]³⁻ decomposition (protic impurity, light, free CN⁻) | IR bands of [Mo(CN)7]⁴⁻ or terminal CN; the 2047 cm⁻¹ band | vacancies; 7-coordinate Mo (the 110 K phase); C-bound V and antisite-like defects (opposite-spin levels near both edges, §4) | ≤ 1 ppm H2O, dark, ≤ 6 CN per Mo, no cyanide excess |
| 4 | Oxidation (V(II) → V(III)/VO²⁺; Mo(III) → Mo(IV)) | M drifts with exposure; T_C falls (vanadyl V–Cr: 115 K) | M = w − 3v − u moves through zero; V(III) holes stay V-majority (KVAC, §2); vanadyl destroys order | glovebox, sealed samples, no air step |
| 5 | Over-heating | irreversible M(T) above 340 K | decomposition | anneal ≤ 320 K; read T_C off XMCD or Arrott data below 330 K plus extrapolation |
| 6 | Wrong linkage (Mo N-bound, V C-bound) | IR, Mo L-edge / V K-edge XAS | windows lost | not expected: uphill by +0.57 eV per pair (§3.4) unless free CN⁻ is present |
| 7 | V→Mo charge transfer | M_FM ≠ 6, metallic | windows lost | not the PBE+U/HSE ground state (Gate 0); no action |

### 5.5 What a successful sample looks like (acceptance tests)

- **Composition and structure:** ICP K:V:Mo = 1.00 : 1.00 : 1.00 (± 0.02); single cubic phase.
- **Lattice constant:** a ≈ 10.73 Å. That is the PBE+U 10.86 Å scaled by the KVCr ratio exp/calc = 10.55/10.68 [derived].
- **IR:** one bridging ν(CN) band near 2100–2120 cm⁻¹, with no 2047 cm⁻¹ or [Mo(CN)7]⁴⁻ bands.
- **Magnetisation:** M_sat(5 K) ≤ 0.1 μB/f.u., i.e. |3v + u − w| ≤ 0.1.
- **Element-specific checks:** V and Mo XMCD give antiparallel sublattices of equal magnitude, within the covalency reduction.
- **Ordering:** remanence persisting at 330 K.

## 6. Files and how to collect the rest

**Job status (10:15 PDT)**

| job | call id | state | expected |
|---|---|---|---|
| `lcm/pba/def/KVCr_VAC_r` | fc-01M41AA0620YQXGFYMT0ZYDQ1B | **done** 10:07 (0.41 h) | – |
| `lcm/pba/def/KVMo_VAC_r` | fc-01M41AA1T9Z2PQ21W2DPCMRSMY | **done** 10:07 (0.42 h) | – |
| `lcm/pba/def/KVCr_VAC_w1` | fc-01M41BZ9NAZ1WFTJAQ36VQSD88 | running since 10:11 (P1, 4 k, npool 8) | relax cap 45 min, then scf + nscf (≈ 36 k) + projwfc: done ≈ 12:00–12:30 |
| `lcm/pba/def/KVMo_VAC` (original) | fc-01M40V7XG7GZ2WV91EBK4688DW | relax, 63 steps | relax cap at ≈ 10:19, then scf/nscf/projwfc: done ≈ 11:15 |
| `lcm/pba/def/KVCr_VAC` (original) | fc-01M40V7S5RP83C7YEKQVXWV9Q2 | relax, 16 + 41 steps (one preemption) | relax cap ≈ 12:11; done ≈ 13:00 |
| `lcm/pba/def/KVMo_KVAC_polM` | fc-01M40V7VNEZRRY5R5VAP070PH1 | **done** (3.95 h) | – |
| `lcm/pba/def/KVMo_KVAC_polV` | fc-01M40V7W9P19HH6NT71K523AGC | relax, 40 steps (restarted 07:21) | relax cap ≈ 10:21; done ≈ 11:00 |

No job was cancelled.

**Collect the rest**
```
source .venv/bin/activate
python tracks/lcm/pba/vac/collect_vac.py --fetch     # VAC, VAC_r, VAC_w1 (+ P, refs) -> vac/results/, vac_summary.json
python tracks/lcm/pba/defects/collect.py --fetch     # all Gate-1 cells incl. polM / polV and the original VAC jobs
python tracks/lcm/pba/vac/summarize.py               # rewrites summary.json (prefers the original VAC cell once done)
```

What to check when they land:
- **`KVCr_VAC_w1`:** the position of the lowest opposite-spin empty state, `opp_bot_rel_hostCBM`. If it moves above +0.3 eV, the KVCr CB kill is an empty-void artefact.
- **`KVCr_VAC` and `KVMo_VAC`:** confirm the `_r` numbers within ±0.1 eV.
- **polV:** confirm the −0.27 eV self-trapping and that its level stays V-majority.

**Files**

| file | purpose |
|---|---|
| `submit_vac_r.py` | builds `inputs/<id>_r.json` from the last relax positions (fixed volume path) and submits the SKIP_RELAX reruns |
| `submit_vac_w.py` | builds and submits `KVCr_VAC_w1` |
| `collect_vac.py` | collector; reuses `../defects/collect.py` `analyse()` with its own cache in `results/` |
| `summarize.py` | writes `summary.json`: cell numbers, void state, polaron harvest, lability energetics, concentration scenarios |
| `manifest.json` | jobs, call ids, md5 of `pba_defect.py` |
| `inputs/` | `KVCr_VAC_r.json`, `KVMo_VAC_r.json`, `KVCr_VAC_w1.json` |
| `results/` | cached volume files (fixed paths), `vac_summary.json` |

The job script is `../defects/jobsrc/pba_defect.py`, unchanged. No file outside `tracks/lcm/pba/vac/` was edited; `../defects/results/` is only refreshed by its own collector.

**Main sources** (details in §3)
- Magott et al., Adv. Sci. 13, e11285 (2026).
- Beauvais & Long, JACS 124, 2110 (2002).
- Nelson et al., Inorg. Chem. 46, 10093 (2007).
- Nelson & Miller, Inorg. Chem. 47, 2526 (2008).
- Avendaño et al., Inorg. Chem. 2010, doi 10.1021/ic901681e.
- Garde, Villain & Verdaguer, JACS 124, 10531 (2002).
- Fe–Cr linkage isomerism:
  - JACS 127, 4580 (2005);
  - JACS 130, 15519 (2008);
  - Inorg. Chem. 2013, doi 10.1021/ic302764k;
  - Shriver et al., Inorg. Chem. 4, 725 (1965).
- K4[Mo(CN)7] photoswitching, Nat. Commun. 2025, doi 10.1038/s41467-025-63523-x.
- Schart et al., Inorg. Chem. 63, 22856 (2024).

Abstracts were read via the Europe PMC API. No ACS/Wiley/PMC full text was readable (bot challenges, not bypassed).
