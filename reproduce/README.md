# Re-running the calculations (Level 3)

Two ways to re-run the original Quantum ESPRESSO inputs from scratch.

## By hand, anywhere QE runs

```bash
conda install -c conda-forge qe=7.5 openmpi      # the version used for every run here
bash reproduce/get_pseudos.sh pseudo             # downloads PseudoDojo v0.4 PBE SR standard NC and checks md5 sums
NP=16 NK=8 PSEUDO=$PWD/pseudo bash reproduce/run_qe.sh \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.in \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.in \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_FM.in
```

`run_qe.sh` never edits the original input files. For each `X.in` it writes `X.repro.in`, which differs only in `pseudo_dir` and `outdir`, then runs it and writes `X.repro.out` next to the original.

Give an scf input before the nscf input that reads it. Each pair shares a prefix and a folder, so it shares a scratch directory.

Compare your results with the originals in two ways:

- run `python tools/verify.py` (originals);
- open your `.repro.out` files with `tools/qe_parse.py`.

```python
import sys; sys.path.insert(0, "tools")
from qe_parse import parse_pw, last_bands, spin_windows
print(spin_windows(last_bands("materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out")))
```

Approximate cost on 32 cores:

| Run | Time |
|---|---|
| KV[Cr(CN)₆] PBE+U scf | ~4 min |
| KV[Cr(CN)₆] PBE+U nscf | ~5 min |
| KV[Cr(CN)₆] HSE06 (64 cores) | ~1 h |
| YBaMnFeO₅ 36-atom PBE+U scf | 7–11 min |
| YBaMnFeO₅ HSE06 | ~5 h |

## In the cloud, in one command (Modal)

```bash
pip install modal && modal setup
modal run reproduce/modal_repro.py --case all        # or kvcr_pbeu, kvcr_hse, ybmfo_stack_U44, kvcr_pbeu_sssp
python tools/verify.py --repro
```

The script does four things:

- builds a fresh image with conda-forge QE 7.5, with no campaign code;
- downloads the pseudopotentials from their public sources;
- ships them to the container, where their md5 sums are checked;
- runs the inputs and writes the outputs to `reproduce/results/<case>/`.

**Pitfall we hit:** the container platform sets `OMP_NUM_THREADS` to the core count. With an OpenMP-enabled QE build, 32 MPI ranks then start 32 threads each, which makes runs about 50× slower. Both scripts force `OMP_NUM_THREADS=1`.

## What we got

See [RESULTS.md](RESULTS.md).
