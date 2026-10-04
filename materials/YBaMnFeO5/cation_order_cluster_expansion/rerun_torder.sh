#!/usr/bin/env bash
# Level-2 check of the decisive negative result for YBaMnFeO5: the Mn/Fe order-disorder temperature.
# Re-runs the canonical Monte Carlo of the fitted paramagnetic cluster expansion (ce2_YBMFO.json, headline
# variant E_PM_best_bw_bayes) in a scratch copy, so the recorded torder_*.json files are not overwritten.
#
#   bash rerun_torder.sh [n_bootstrap]     (default 3; the recorded run used 30)
#
# Recorded: T_order = 915 / 965 / 915 K (L = 8 / 12 / 16), headline 950 (+250/-150) K.
# Needs numpy and numba. Takes roughly 10-30 min on a laptop.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NBOOT="${1:-3}"
TMP="$(mktemp -d)"
cp "$HERE/torder.py" "$HERE/mc_ce.py" "$HERE/ce2_YBMFO.json" "$TMP/"
( cd "$TMP" && python torder.py ce2_YBMFO.json E_PM_best_bw_bayes best --boot "$NBOOT" )
cp "$TMP/torder_YBMFO_E_PM_best_bw_bayes.json" "$HERE/rerun_torder_YBMFO_E_PM_best_bw_bayes.json"
python - "$HERE/rerun_torder_YBMFO_E_PM_best_bw_bayes.json" "$HERE/torder_YBMFO_E_PM_best_bw_bayes.json" <<'EOF'
import json, sys
new, old = (json.load(open(p)) for p in sys.argv[1:3])
print("re-run  T_peak by L:", new["Tpeak"], " T_order(L=16) =", new["T_order_K"], "K")
print("recorded T_peak by L:", old["Tpeak"], " T_order(L=16) =", old["T_order_K"], "K")
EOF
rm -rf "$TMP"
