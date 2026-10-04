"""Equilibrium LRO/SRO of the B-site CE at synthesis / annealing temperatures (canonical MC, mc_ce.run, cooling sequence).
usage: python sro.py ce2_<TAG>.json KEY [CSEL=best] [L=16]   -> sro_<TAG>_<KEY>[_CSEL].json
Reports per T: eta(q) for the ground-state ordering vector (LRO; finite-size floor ~ 1/sqrt(N)), pair correlations <ss> per
class, the fraction of hetero (Mn-Fe / Cu-Fe) pairs per NN class, and for rock-salt (q = (1,1,1)) the antisite fraction
x_AS = (1 - eta_RS)/2 (meaningful only below T_order)."""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mc_ce

HERE = os.path.dirname(os.path.abspath(__file__))
cej = json.load(open(sys.argv[1])); key = sys.argv[2]
csel = sys.argv[3] if len(sys.argv) > 3 else "best"
L = int(sys.argv[4]) if len(sys.argv) > 4 else 16
ent = cej["ce"][key + "_best"] if csel == "best" else [e for e in cej["ce"][key] if "+".join(e["classes"]) == csel][0]
V = ent["V_meV_per_bond"]
TS = [2000, 1800, 1673, 1573, 1473, 1373, 1273, 1173, 1073, 973, 873, 773, 673, 573]
res, s = mc_ce.run(V, L, L, TS, nsweep=4000, nequil=3000, seed=11)
# ground-state ordering vector from the lowest-T point
q0 = res[-1]["q_max"]
out = []
print("CE", key, csel, {k: round(v, 1) for k, v in V.items()}, "ground q", q0)
print("%6s %8s %7s %7s %7s %7s  %s" % ("T(K)", "E", "eta_q0", "etaRS", "etaLAY", "x_AS", "hetero fraction ip1/ap/Y"))
for r in res:
    sr = r["sro"]
    het = {k: (1 - sr.get(k, 0.0)) / 2 for k in ("ip1", "ap", "Y")}
    eq0 = r["eta"].get(str(q0), 0.0)
    rec = {"T": r["T"], "E_meV_fu": r["E_meV_fu"], "eta_q0": eq0, "eta_RS": r.get("eta_RS"), "eta_LAY": r.get("eta_LAY"),
           "x_AS_RS": (1 - r["eta_RS"]) / 2 if r.get("eta_RS") is not None else None, "sro": sr, "hetero_frac": het, "q_max": r["q_max"]}
    out.append(rec)
    print("%6.0f %8.1f %7.3f %7.3f %7.3f %7.3f  %.2f/%.2f/%.2f" % (r["T"], r["E_meV_fu"], eq0, r.get("eta_RS", -1), r.get("eta_LAY", -1),
                                                               rec["x_AS_RS"] if rec["x_AS_RS"] is not None else -1, het["ip1"], het["ap"], het["Y"]))
json.dump({"key": key, "csel": csel, "V": V, "L": L, "q_ground": str(q0), "floor_eta": 1 / np.sqrt(L ** 3), "T": out},
          open(f"{HERE}/sro_{cej['tag']}_{key}{'_' + csel if csel != 'best' else ''}.json", "w"), indent=1, default=str)
