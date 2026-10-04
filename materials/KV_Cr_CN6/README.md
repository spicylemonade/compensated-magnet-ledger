# KV[Cr(CN)₆]: computed data

A Prussian-blue-type framework. Cr³⁺ (3 unpaired electrons) sits on the carbon end of each cyanide and V²⁺ (3 unpaired electrons) on the nitrogen end; K⁺ fills the cubic holes.

- **Made:** Holmes & Girolami, *J. Am. Chem. Soc.* 121, 5593 (1999), as a hydrated powder, KV[Cr(CN)₆]·2H₂O.
- **Measured:** orders magnetically at 376 K, with a saturation moment of 0.125 μB per formula unit (0 expected for perfect stoichiometry).

**What the calculations here say (ideal anhydrous crystal):**

- It is a Luttinger-compensated (fully compensated ferrimagnetic) magnet: V and Cr moments are antiparallel, and the net spin is exactly zero.
- Spin space group 216.216.1.1.L, with Zeeman-type (s-wave) splitting.
- Both band edges sit in the V-majority spin channel. HSE06 gives a 2.09 eV gap with spin windows of 2.64 eV (holes) and 1.57 eV (electrons); PBE+U gives smaller windows (1.4–2.6 / 0.9–1.4 eV over the U grid) and the same edge spins.

What is uncertain is spelled out in `../../LEDGER.md`: water, vacancies, composition, transport, and finite temperature.

| Folder | Contents |
|---|---|
| `runs/` | Raw `pw.x` inputs and outputs, one folder per study (see `../../ledger/runs.csv`). Large outputs are gzipped |
| `structures/` | Relaxed cells: CIF, POSCAR, and magnetic CIF with the compensated moments |
| `analysis_data/` | The agents' collector outputs (Gate-0 summary, defects, vacancies, phonons, tautomers, T_C model) and the spin-group check |
| `original_scripts/` | Job scripts and collectors exactly as run. They import the campaign's cloud tooling, so read them as provenance; use `tools/` to re-check numbers |
| `working_notes/` | The agents' working READMEs, including the second referee report (`pba_REREVIEW.md`) |
| `AGENT_DOSSIER_2026-10-03.md` | The agents' dossier with every correction they made. Its sections on KV[Mo(CN)₆] and Cr[Mo(CN)₆] are predictions not covered by the blog |

Settings and conventions are the same as in `../../LEDGER.md`.
