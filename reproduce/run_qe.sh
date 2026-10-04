#!/usr/bin/env bash
# Re-run Quantum ESPRESSO inputs from this repository with your own pw.x.
#
#   NP=16 NK=8 PSEUDO=$PWD/pseudo SCRATCH=/tmp/qe bash reproduce/run_qe.sh path/to/a.in [path/to/b.in ...]
#
# The original inputs point at the paths of the cloud containers they ran in
# (pseudo_dir = '/data/pseudo/dojo_sr', outdir = '/tmp/scratch/<prefix>'). This script leaves the
# original files untouched: for each X.in it writes X.repro.in with only those two lines changed,
# runs it, and writes X.repro.out next to it. Inputs are run in the order given, so give an scf
# input before the nscf input that reads it (they share prefix and directory, hence outdir).
#
# Needs: pw.x from Quantum ESPRESSO 7.5 (e.g. `conda install -c conda-forge qe=7.5 openmpi`),
#        and the pseudopotentials from reproduce/get_pseudos.sh.
set -euo pipefail
NP="${NP:-4}"           # MPI ranks
NK="${NK:-1}"           # k-point pools (must divide NP)
PSEUDO="${PSEUDO:?set PSEUDO to the directory with the PseudoDojo .upf files}"
SCRATCH="${SCRATCH:-/tmp/qe_scratch}"
MPIRUN="${MPIRUN:-mpirun -np $NP}"
export OMP_NUM_THREADS=1   # MPI-only runs; with OpenMP-enabled builds, unset threads can oversubscribe the cores badly
for inp in "$@"; do
  dir="$(cd "$(dirname "$inp")" && pwd)"
  stem="$(basename "$inp" .in)"
  prefix="$(sed -n "s/^ *prefix *= *'\(.*\)'.*/\1/p" "$inp" | head -1)"
  tag="$(echo "$dir" | md5sum 2>/dev/null | cut -c1-8 || md5 -q -s "$dir" | cut -c1-8)"
  out_dir="$SCRATCH/${tag}/${prefix}"
  mkdir -p "$out_dir"
  rin="$dir/$stem.repro.in"
  sed -e "s|^ *pseudo_dir *=.*|  pseudo_dir = '$PSEUDO'|" \
      -e "s|^ *outdir *=.*|  outdir = '$out_dir'|" "$inp" > "$rin"
  echo "[$(date +%H:%M:%S)] pw.x  $inp  (np=$NP, nk=$NK)"
  t0=$(date +%s)
  $MPIRUN pw.x -nk "$NK" -in "$rin" > "$dir/$stem.repro.out" 2>&1 || echo "  pw.x returned non-zero, see $dir/$stem.repro.out"
  echo "  done in $(( $(date +%s) - t0 )) s -> $dir/$stem.repro.out"
done
