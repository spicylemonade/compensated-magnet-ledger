#!/usr/bin/env bash
# Download the PseudoDojo v0.4 PBE scalar-relativistic standard norm-conserving pseudopotentials
# and check them against the md5 sums of the files used in the original runs.
#
#   bash reproduce/get_pseudos.sh [target_dir]      (default: ./pseudo)
#
# Exit status is non-zero if any file needed by this repository is missing or differs.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST="${1:-pseudo}"
URL="https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz"
mkdir -p "$DEST"
cd "$DEST"
if [ ! -f .downloaded ]; then
  echo "downloading $URL"
  # Integrity is checked by md5 below, so a failed TLS chain check (seen with some wget builds) is retried without it.
  if command -v curl >/dev/null; then curl -sSL "$URL" -o dojo.tgz || curl -sSLk "$URL" -o dojo.tgz
  else wget -q "$URL" -O dojo.tgz || wget -q --no-check-certificate "$URL" -O dojo.tgz; fi
  tar xzf dojo.tgz && rm dojo.tgz
  # some tarball versions unpack into a sub-folder; flatten it
  for d in */; do [ -d "$d" ] && find "$d" -name '*.upf' -exec mv -n {} . \; ; done
  touch .downloaded
fi
md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
bad=0
while read -r el file sum zval; do
  case "$el" in ''|\#*) continue;; esac
  if [ ! -f "$file" ]; then echo "MISSING  $file"; bad=1; continue; fi
  got=$(md5 "$file")
  if [ "$got" = "$sum" ]; then echo "ok       $file  $got"; else echo "DIFFERS  $file  got $got expected $sum"; bad=1; fi
done < "$HERE/pseudo_md5.txt"
if [ $bad -ne 0 ]; then
  echo "WARNING: at least one pseudopotential differs from the files used in the original runs." >&2
  echo "         Results can still be close, but are no longer a like-for-like reproduction." >&2
  exit 1
fi
echo "all pseudopotentials match the original runs"
