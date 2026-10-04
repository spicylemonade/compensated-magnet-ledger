# YBaMnFeO₅ (rock-salt-ordered): computed data

A designed layered "112" perovskite, P4/n, in which Mn²⁺ and Fe³⁺ (both with 5 unpaired electrons) alternate in a three-dimensional checkerboard on the square-pyramidal sites.

**On paper:**

- G-type compensated ground state (Mn↑ Fe↓) at every U from 0 to 6 eV, confirmed by HSE06;
- net spin exactly zero;
- HSE06 gap 2.35 eV, with both band edges in one spin channel (windows 1.00 / 1.40 eV);
- Néel temperature 417 K (classical Monte Carlo), about 490 K after calibration;
- 2.6 meV/atom above the 0 K convex hull.

**The decisive negative result.** A paramagnetic cluster expansion puts the Mn/Fe order-disorder temperature at 950 (+250/−150) K. That is below the temperatures at which the cations can rearrange during synthesis. A neighbouring Mn/Fe swap collapses the spin-selective gap. The agents therefore downgraded the compound to a design study; see `AGENT_DOSSIER_2026-10-02.md`, section V.

| Folder | Contents |
|---|---|
| `runs/` | Raw `pw.x` inputs and outputs: relaxation, HSE06 ground state and competitor, U scans of the stackings, 72-atom spin-flip and antisite cells, the 24 exchange-fit configurations |
| `exchange_fit_and_Neel_T/` | `rerun_fit_and_mc.py` re-fits the exchange constants from the raw outputs and recomputes T_N (Level 2); the original fit and Monte Carlo files are alongside |
| `cation_order_cluster_expansion/` | The cluster-expansion data (`ce2_YBMFO.json`), per-arrangement DFT results, Monte Carlo outputs, the YBaCuFeO₅ benchmark, and `rerun_torder.sh` (Level 2) |
| `Neel_T_calibration/` | The same-protocol calibration on YFeO₃, MnFe₂O₄ and YBaMn₂O₅ |
| `structures/` | Relaxed cells, including superseded 75 Ry versions kept for the record |
| `analysis_data/`, `figures/` | Hull, phonons, magnons, calibration summaries; band and phonon plots |
| `original_job_inputs/` | The JSON job inputs as submitted |
| `AGENT_DOSSIER_2026-10-02.md` | The agents' dossier. Section V (verdict and corrections) supersedes the original abstract below it |
