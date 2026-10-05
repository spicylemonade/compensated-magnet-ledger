# HSE06, KV[Cr(CN)₆]·2H₂O, carried to convergence

The agents' HSE06 run for the dihydrate (`../F_hse06_dihydrate_LCM/`) stopped after 12 exact-exchange (outer) cycles, before the loop converged. This folder is the same input run to completion on 2026-10-04:

- **Input:** identical to `../F_hse06_dihydrate_LCM/hse_lcm.in` apart from `pseudo_dir` and `outdir` (container paths).
- **Environment:** a fresh cloud container (Modal), Quantum ESPRESSO 7.5 from conda-forge, 64 MPI processes (`-nk 8`), PseudoDojo pseudopotentials downloaded from the public source and md5-checked (`run_summary.json`). Procedure: `reproduce/modal_repro.py --case kvcr_hydrate_hse`.
- **Reproduction:** the first 12 outer cycles reproduce the agents' run to every printed digit (all 14 total-energy lines and all 12 exchange-error values).
- **Convergence:** 41 outer cycles. The exchange error fell to 1.1e-6 Ry by cycle 10, rose again to 1.07e-5 Ry at cycle 28 (a slow oscillation; the agents' run stopped at cycle 12, inside it), then converged to 2.0e-7 Ry (threshold 2.1e-7). QE printed the converged `!!` total energy. Wall time 4.4 h.
- **Result** (claims K19b, K20b): gap 2.14 eV; both band edges in the same spin channel; hole window 2.31 eV, electron window 1.40 eV; net moment 0. The agents' stopped run gave 2.02 eV, 2.43 eV and 1.42 eV.
