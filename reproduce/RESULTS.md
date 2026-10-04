# Independent re-runs: results

These were re-run on 2026-10-04 in a fresh cloud container (Modal) using `reproduce/modal_repro.py`:

- a clean conda-forge Quantum ESPRESSO 7.5 image, with no campaign code;
- pseudopotentials downloaded from their public sources;
- the unmodified input files from this repository (only `pseudo_dir` and `outdir` rewritten).

**Pseudopotential check.** Every PseudoDojo file used by these materials matched the md5 sum of the file used in the original runs:

- `kvcr_pbeu`: 12 files checked, all identical

## Reproduction (same inputs, same pseudopotentials)

| Calculation | Quantity | Original run | Independent re-run |
|---|---|---|---|
| KV[Cr(CN)6] PBE+U (U 3/3 eV) | net spin moment (Bohr mag/cell) | -0.00 | -0.00 |
|  | total energy (eV) | -7922.403400 | -7922.403406 |
|  | V / Cr sphere moments | -2.016 / 2.318 | -2.017 / 2.318 |
|  | band gap, nscf 8x8x8 (eV) | 1.9529 | 1.9509 |
|  | hole / electron window (eV) | 2.0230 / 1.1462 | 2.0251 / 1.1457 |
|  | both edges same spin | True | True |
|  | E(FM) - E(compensated) (meV/ion) | 152.24 | 152.24 |
| KV[Cr(CN)6] HSE06 | net spin moment | -0.00 | 0.00 |
|  | total energy (eV) | -7840.410303 | -7840.410303 |
|  | band gap (eV) | 2.0909 | 2.0909 |
|  | hole / electron window (eV) | 2.6365 / 1.5674 | 2.6365 / 1.5674 |
| YBaMnFeO5 PBE+U (U 4/4 eV) | G-type total energy (eV) | -42025.921471 | -42025.921471 |
|  | net spin moment (G) | 0.00 | 0.00 |
|  | E(Y-layer flip) - E(G) (meV/mag. ion) | 8.881 | 8.881 |

Differences in the last digits come from parallelisation and random starting wavefunctions. They are far below the precision any claim relies on.

Check it yourself: `python tools/verify.py --repro`.

## Robustness: a different pseudopotential family (not a reproduction)

This is the same KV[Cr(CN)6] PBE+U calculation with SSSP 1.3 efficiency pseudopotentials (GBRV ultrasoft for V and Cr, PSlibrary PAW for K and C, THEOS for N; 60/480 Ry) instead of PseudoDojo. Hubbard projectors depend on the pseudopotential, so the numbers are not expected to match exactly. The physical picture is unchanged.

| Quantity | PseudoDojo (original) | SSSP 1.3 (re-run) |
|---|---|---|
| net spin moment | 0.00 | -0.00 |
| V / Cr sphere moments | -2.016 / 2.318 | -2.001 / 2.325 |
| band gap (eV) | 1.953 | 1.795 |
| hole / electron window (eV) | 2.023 / 1.146 | 2.015 / 1.350 |
| both edges same spin | True | True |

