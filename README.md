# Two magnets that add up to zero: a computational ledger

This repository holds the data behind a blog post about two candidate **room-temperature "Luttinger-compensated" magnetic semiconductors**. In these magnets two different kinds of atom carry exactly opposite spins, so the material has no net magnetisation, yet its electrons are still sorted by spin as they would be in a ferromagnet.

- **YBaMnFeO₅**: a new design. On paper it is nearly ideal. However, our own calculations predict that the atomic ordering it needs cannot be reached by normal synthesis.
- **KV[Cr(CN)₆]**: a Prussian-blue-type magnet first made in 1999 that orders at 376 K. Our calculations say it already has this electronic structure; as far as we could find, nobody had described it this way before.

Everything here was computed by a team of AI agents (Claude, by Anthropic), which ran density functional theory (DFT) calculations on rented cloud computers. The project was set up by someone who is **not** a chemist. That is why this repository exists: so that people who *are* chemists and physicists can check the work.

**Read the blog post:** [BLOG.md](BLOG.md)

## Status of the main claims

| | YBaMnFeO₅ (designed) | KV[Cr(CN)₆] (made in 1999) |
|---|---|---|
| Zero net magnetisation (ideal crystal, DFT) | yes | yes |
| Electrons spin-sorted at both band edges, same spin | yes: windows of 1.0 / 1.4 eV (HSE06) | yes: windows of 2.6 / 1.6 eV (HSE06) |
| Band gap (HSE06) | 2.35 eV | 2.09 eV |
| Orders above room temperature | predicted: about 420 K raw, about 490 K calibrated | **measured: 376 K** (Holmes & Girolami, 1999) |
| Can it be made? | **Probably not with known methods.** The required Mn/Fe checkerboard is predicted to scramble at about 950 K, below the temperatures where atoms can move during synthesis, and scrambling destroys the effect | **Already made**, but only once, as a hydrated powder with small imbalances (0.125 μB/f.u.) |
| Survives real-world imperfections? | No: one swapped Mn/Fe pair closes the spin-selective gap | Water and vacancies: HSE06 says yes, PBE+U says the hole window shrinks a lot. Unresolved |
| Measured spin polarisation, conductivity or band gap | none (not made) | none yet |
| Grade from the agents' own adversarial review | design study, not a realizable discovery | solid "identification plus numbers", not a breakthrough |

Each number above is traced in [LEDGER.md](LEDGER.md) to the exact input and raw output files and to a script that recomputes it.

## Check our work in three levels

**Level 1: check the arithmetic (under a second, laptop).** Recompute every number from the raw Quantum ESPRESSO outputs in this repository:

```bash
pip install numpy
python tools/verify.py            # 55 pass, 0 fail, 3 not computable (experiment/literature) as of 2026-10-04
```

**Level 2: re-run the models (minutes, laptop).** Re-derive the exchange constants and Néel temperature of YBaMnFeO₅ from the raw energies, and redo the cation-order simulation that rules it out:

```bash
pip install numpy pymatgen numba
python materials/YBaMnFeO5/exchange_fit_and_Neel_T/rerun_fit_and_mc.py              # recorded 417 K; re-run gives 414 ± 4 K
bash   materials/YBaMnFeO5/cation_order_cluster_expansion/rerun_torder.sh 3         # recorded and re-run: 915 / 965 / 915 K
```

**Level 3: re-run the quantum-mechanical calculations (CPU hours, cluster or cloud).** Download the same pseudopotentials (md5-checked) and run the original input files with Quantum ESPRESSO 7.5:

```bash
conda install -c conda-forge qe=7.5 openmpi
bash reproduce/get_pseudos.sh pseudo
NP=16 NK=8 PSEUDO=$PWD/pseudo bash reproduce/run_qe.sh \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.in \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.in
```

We did exactly this in a fresh cloud container for four sets of runs; see [reproduce/RESULTS.md](reproduce/RESULTS.md). There is also a one-command cloud version, `reproduce/modal_repro.py`.

## What is in here

```
BLOG.md                         the blog post (plain language)
LEDGER.md                       every claim -> files -> script -> status, plus all known caveats
ledger/claims.csv               the same claims, machine-readable
ledger/runs.csv                 index of all 120 included pw.x calculations + 8 post-processing inputs
ledger/harvest_manifest.json    provenance: original storage path and md5 of every raw file
tools/verify.py                 Level-1 checker (numpy only)
tools/qe_parse.py               small pw.x output parser (eigenvalues, moments, energies, cells)
tools/magtools.py               Heisenberg fit + classical Monte Carlo (from the campaign)
reproduce/                      Level-3: pseudopotential checksums, run scripts, independent re-run results
materials/KV_Cr_CN6/
    runs/                       raw inputs/outputs: PBE+U, HSE06, hydrate, vacancies, defects, bands, control
    structures/                 relaxed cells (CIF, POSCAR) and magnetic CIFs
    analysis_data/              collector outputs (defects, phonons, tautomers, T_C model), spin-group check
    original_scripts/           the agents' job scripts and collectors, as run
    working_notes/              the agents' working READMEs and second referee report
    AGENT_DOSSIER_2026-10-03.md the agents' own dossier, including corrections and retractions
materials/YBaMnFeO5/
    runs/                       raw inputs/outputs: relaxation, HSE06, U scans, spin flips, antisites, exchange fit
    exchange_fit_and_Neel_T/    fit + Monte Carlo (Level-2 script)
    cation_order_cluster_expansion/   the decisive negative result (Level-2 script)
    Neel_T_calibration/         YFeO3 / MnFe2O4 / YBaMn2O5 calibration runs and analysis
    structures/ analysis_data/ figures/ original_job_inputs/
    AGENT_DOSSIER_2026-10-02.md the agents' own dossier, including the downgrade verdict
```

About 34 MB in total; large text outputs are gzipped (the tools read them directly).

## How this was produced

- **When:** 1–4 October 2026.
- **Who:** several Claude agents working in parallel, each on its own "lane" of a broader search for unusual magnetic materials. This repository covers only the Luttinger-compensated-magnet lane. The other lanes did not find anything above their bars, and their results are not included.
- **How the agents worked:**
  - pre-registered pass/fail rules before running decisive calculations;
  - adversarial "referee" agents that tried to kill each claim;
  - literature searches for prior art.
  
  Several claims were retracted or corrected along the way. The dossiers in `materials/*/AGENT_DOSSIER_*.md` keep that record.
- **Compute:** Quantum ESPRESSO 7.5 on Modal cloud CPUs. The Track-L lane submitted about 750 jobs, most of them not part of this repository.
- **Human role:** setting goals, asking questions, and deciding what to publish. No calculation was checked by a human expert before publication, which is the point of publishing the ledger.

## Caveats in one paragraph

These are DFT predictions for ideal crystals at zero temperature.

- Density functional theory with a Hubbard U, and even HSE06, can misplace energy levels by tenths of an eV.
- The two methods used here disagree about how much water harms KV[Cr(CN)₆].
- Nobody has measured the spin polarisation, band gap or conductivity of either material.
- YBaMnFeO₅'s key ordering is predicted to be unreachable.
- KV[Cr(CN)₆]'s only sample is one hydrated powder from 1999.

The full list is in [LEDGER.md → Known caveats](LEDGER.md#known-caveats-and-method-issues).

## Experiments that would settle it

1. **Re-make KV[Cr(CN)₆]** and measure its composition and saturation magnetisation (the cheapest test).
2. **Element-specific X-ray magnetic circular dichroism** (V and Cr L-edges) and **magneto-optics** at near-zero magnetisation, across a series of compositions.
3. **Spin-resolved photoemission** of the top valence band of anhydrous KV[Cr(CN)₆]: the prediction is about 100 % one spin over about 2 eV, reversing with the magnetic order.

## License and citation

- Code (`tools/`, `reproduce/`, and scripts under `materials/`): MIT.
- Data and text: CC BY 4.0.

See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff).

The experimental compound KV[Cr(CN)₆]·2H₂O is from S. M. Holmes and G. S. Girolami, *J. Am. Chem. Soc.* **121**, 5593 (1999). Pseudopotentials are from PseudoDojo (van Setten et al., *Comput. Phys. Commun.* **226**, 39 (2018)), and calculations used Quantum ESPRESSO (Giannozzi et al., *J. Phys.: Condens. Matter* **29**, 465901 (2017)).
