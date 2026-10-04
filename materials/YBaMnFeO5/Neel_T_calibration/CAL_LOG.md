# lcm_cal: matched T_N calibration (triage item 6) — log

Protocol: geometry 110-Ry vc-relaxed (lcm_hull110 settings); config energies ECUT 90 Ry, PseudoDojo SR NC, U(Mn)=U(Fe)=4 eV ortho-atomic,
k-spacing 0.28, gaussian 0.005 Ry, mixing 0.3. Mapping E = E0 - sum J s_i s_j (|S|^2 in J). Classical MC (Metropolis + over-relaxation,
unit spins), Binder crossings with jackknife. Multi-sublattice Tyablikov RPA with J_q = J/(S_a S_b).
Scripts: infra/jobs_src/lcm_cal/ (callib.py pair classes/fit/RPA; calmc.py MC kernel; cal_scf.py; cal_relax.py).

- 00:50 callib validated: jlib classes + fit reproduce J9 exactly (rms 2.82, LOO 6.35 meV/cell); RPA reproduces sc Watson-integral T_c (15.30 vs 15.305 K for J=1 meV classical).
- 00:50 YBaMnFeO5 J9 RPA (75-Ry geometry): quantum S=5/2 574 K; classical RPA 410 K (classical MC 467 K).
- 00:50 submitted lcm/cal/relax_YBaMn2O5 (mp-1182186 P4/nmm reordered to the YBaMnFeO5 template; Mn2+ +4.6 at idx 4,7, Mn3+ -3.8 at idx 5,6).
- 00:56 MC validation (lcm/cal/mc_YBMFO_J9b; new kernel Metropolis+2 OR, L=6..16, 5e4 sweeps): J9 Binder crossings 465-467 K (jk err 0.3-1.7 K) = previous 467 K. Pipeline validated.
- 00:56 YFeO3: 16 D-optimal configs (8 orbit classes to 6.9 A) submitted lcm/cal/YFeO3_g0..4. MnFe2O4 (110-Ry relaxed, a_conv 8.644 A, Neel M=-10): 18 configs (9 classes to 6.11 A) submitted lcm/cal/MnFe2O4_g0..3.
- 01:26 YBaMn2O5 110-Ry relax done: P4/nmm a=5.598, c=7.821 A (exp 5.536/7.615), M=2.0 muB/cell, Mn2+ 4.00 / Mn3+ -3.23 (CO kept). 22 configs (the YBMFO fit set mapped 1:1) submitted lcm/cal/YBaMn2O5_g0..6.
- 01:43 YFeO3 fit (16/16 configs clean, M = 5*sum s): best-LOO 6-class model J(meV) ap3.854 -32.22, ip3.877 -37.03, NNN 5.32 -1.21, 5.47 -0.79/-0.72, 5.64 -1.65; rms 5.0, LOO 11.1 meV/cell(4 Fe). NN-only: -32.22/-37.03 (rms 10.9, LOO 14.7). RPA S=5/2: 695 K (best), 758 K (NN); classical RPA 496/541 K. MC submitted lcm/cal/mc_YFeO3_{best,nn}. MnFe2O4 FM (c01) SCF stuck at 2e-6 Ry with M=28 (not 30: charge-transferred, non-d5) -> cancelled g0, resubmitted c02,c11,c14 (g4) and c01 with conv_thr 3e-7/atom (g5).
- 01:46 YFeO3 Binder MC (L=6..16, 5e4 sweeps): best model 545 K (crossings 541-546), NN-only 591 K. exp(640-648)/MC = 1.17-1.19 (best), 1.08-1.10 (NN). exp/RPA_q(695 K) = 0.92-0.93.
- 02:05 **YBaMn2O5 ground-state check (protocol failure):** at the G-relaxed 110-Ry geometry, E(A-type: in-plane Mn2+-Mn3+ FM, both interlayer AFM) - E(G FiM, exp. order) = +119 (U=0), -41 (U=2), -140 (U=4, protocol), -218 meV/cell(4 Mn) (U=6); FM - G = +580/+242/+41/-102. CO kept in all (Mn2+ 3.75-4.06 / Mn3+ 2.9-3.5, integer M). PBE+U >= 2 eV predicts the WRONG magnetic ground state (in-plane FM) for the structural twin; only U~0 gives the experimental ferrimagnet. Saved YBaMn2O5_Uscan.json.
- 02:35 MnFe2O4 fit (17/18 clean; FM c01 has M=28 not 30 -> charge-transferred, excluded; including it gives rms 121 meV/cell). Best-LOO 9-class: J_AB(Mn-Fe 3.58) -18.8, J_BB(3.06) -3.6, J_AA(3.74) -4.2, NNN -0.3..-0.6, Fe-Fe 6.11 (linear B chain) -1.5/-0.7 meV; rms 1.3, LOO 4.2 meV/cell. NN-only (3 classes): rms 15.4. Binder MC 385 K (best) / 401 K (NN); RPA_q 504/530 K.
- 02:50 YBaMn2O5 fit (22/22 clean, CO kept, integer M): J_ip(Mn2+-Mn3+) **+7.5 (FM)**, J_ap -24.7, J_Y -15.5 meV, NNN -1.3..+0.1; rms 1.3, LOO 2.3 meV/cell. Model GS (annealing) = A-type, as in DFT. Binder MC of the A-type order 191 K (best) / 180 K (NN); RPA_q 249/236 K. J_Y (-15.5) equals YBaMnFeO5's J_Y (-15.4).

## Calibration table (protocol: 110-Ry geometry, 90 Ry, U=4, k 0.28; Binder MC L up to 16 (12 for spinel), 5e4 sweeps; RPA S=5/2 (Mn3+ S=2))

| calibrant | bonds | model (classes, LOO meV/cell) | T_MC (K) | RPA_q (K) | T_exp (K) | exp/MC | exp/RPA_q |
|---|---|---|---|---|---|---|---|
| YFeO3 (Pnma, G) | Fe3+-O-Fe3+ d5-d5 | best (6, 11.1) | 545 +- 3 | 695 | 640-648 | 1.17-1.19 | 0.92-0.93 |
| | | NN (2, 14.7) | 592 +- 2 | 758 | | 1.08-1.10 | 0.84-0.86 |
| MnFe2O4 (normal spinel, Neel) | Mn2+-O-Fe3+ d5-d5 | best (9, 4.2) | 385 +- 1 | 504 | 560-590 (x~0.2 inverted) | 1.46-1.53 | 1.11-1.17 |
| | | NN (3, 20.4) | 401 +- 1 | 530 | | 1.40-1.47 | 1.06-1.11 |
| YBaMn2O5 (112 twin) | Mn2+-O-Mn3+ d5-d4 | best (9, 2.3) **wrong GS (A-type)** | 191 +- 1 | 249 | 160-170 (FiM) | 0.84-0.89 | 0.64-0.68 |
| | | NN (3, 9.3) | 180 | 236 | | 0.89-0.94 | 0.68-0.72 |
| YBaMnFeO5 J9 (75 Ry) | Mn2+-O-Fe3+ d5-d5 | 9 classes | 467 +- 2 | 574 (cl. RPA 410) | - | - | - |

Calibrated YBaMnFeO5 (J9 numbers): MC route 548-715 K (d5 calibrants, best models; YFeO3 548-555, MnFe2O4 679-715); RPA route 529-672 K.
Floors: 505 K (YFeO3 NN-only factor 1.08); 392 K (YBaMn2O5 MC factor 0.84, invalid calibrant); 369 K (YBaMn2O5 RPA factor).
Factors to apply to the 110-Ry T_MC: f_MC = 1.17-1.53 (central ~1.3); f_RPA = 0.92-1.17.
Kill criterion: d5 calibrants are UNDERestimated (exp/MC 1.17-1.53 > 1) -> no downward revision required by the matched d5 calibration; YBaMn2O5 is overestimated (0.84-0.89) but the U=4 protocol gets its ground state wrong, so it is not a valid T calibrant; it flags that the protocol is unreliable for e_g-orbital-active (Mn3+) bonds in this structure type.
