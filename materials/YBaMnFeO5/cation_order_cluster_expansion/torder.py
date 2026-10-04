"""Order-disorder temperature from canonical MC of the fitted PM cluster expansion.
usage: python torder.py ce_<TAG>.json [key=dE_PM_fit] [classes=ip1+ap+Y+...|best] [--boot N]
- MC (mc_ce.run) for L = 8, 12, 16 (Lz = L), fine T grid around the transition; T_order from (i) peak of C(T) for the
  largest L and (ii) crossing of the Binder cumulant U4 of the dominant order parameter (q of the T -> 0 ground state)
  for successive L.
- uncertainty: MC with L = 8 for N bootstrap CE parameter sets (CE refit on resampled training structures).
Writes torder_<TAG>_<key>.json"""
import json, sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mc_ce

HERE = os.path.dirname(os.path.abspath(__file__))
cej = json.load(open(sys.argv[1]))
key = sys.argv[2] if len(sys.argv) > 2 and not sys.argv[2].startswith("--") else "dE_PM_fit"
csel = sys.argv[3] if len(sys.argv) > 3 and not sys.argv[3].startswith("--") else "best"
NB = int(sys.argv[sys.argv.index("--boot") + 1]) if "--boot" in sys.argv else 30
if csel == "best":
    ent = cej["ce"][key + "_best"]
else:
    ent = [e for e in cej["ce"][key] if "+".join(e["classes"]) == csel][0]
V = ent["V_meV_per_bond"]
print("CE", key, "+".join(ent["classes"]), "CV %.1f meV/fu" % ent["cv"], {k: round(v, 2) for k, v in V.items()})


def ground_q(V):
    """dominant ordering vector at low T from a quick anneal (L=8)"""
    temps = np.geomspace(4000, 100, 26)
    res, s = mc_ce.run(V, 8, 8, temps, nsweep=400, nequil=300, seed=5)
    return res[-1]["q_max"], res[-1]


def scan(V, L, temps, q, nsw=3000, seed=1):
    res, s = mc_ce.run(V, L, L, temps, nsweep=nsw, nequil=nsw // 2, seed=seed)
    out = []
    for r in res:
        rr = {"T": r["T"], "E": r["E_meV_fu"], "C": r["C_kB_fu"], "sro": r["sro"]}
        rr["eta"] = r["eta"].get(str(q), 0.0)
        out.append(rr)
    # Binder cumulant for q: rerun quick sampling is costly; use the U4 stored when q is the max
    for rr, r in zip(out, res):
        rr["U4"] = r["U4_max"] if r["q_max"] == q else None
        rr["q_max"] = r["q_max"]; rr["eta_max"] = r["eta_max"]
    return out


q0, gs = ground_q(V)
print("ground-state ordering vector q =", q0, "eta =", round(gs["eta_max"], 3), "E =", round(gs["E_meV_fu"], 1), "meV/fu, SRO", {k: round(x, 2) for k, x in gs["sro"].items()})
# coarse scan to bracket
coarse = scan(V, 8, np.linspace(4000, 200, 39), q0, nsw=800)
Cc = np.array([r["C"] for r in coarse]); Tc0 = coarse[int(np.argmax(Cc))]["T"]
print("coarse C-peak (L=8): %.0f K" % Tc0)
temps = np.linspace(Tc0 * 1.35, max(100, Tc0 * 0.6), 28)
fine = {}
for L in (8, 12, 16):
    fine[L] = scan(V, L, temps, q0, nsw=4000 if L < 16 else 3000, seed=L)
    C = np.array([r["C"] for r in fine[L]]); i = int(np.argmax(C))
    print("L=%d C-peak %.0f K (C=%.2f)" % (L, fine[L][i]["T"], C[i]))
Tpk = {L: fine[L][int(np.argmax([r["C"] for r in fine[L]]))]["T"] for L in fine}
# eta(T) half-height for the largest L
L = 16
eta = np.array([r["eta"] for r in fine[L]]); Ts = np.array([r["T"] for r in fine[L]])
# bootstrap
boot = []
B = np.array(cej["ce"][key + "_best"]["boot"] if csel == "best" else ent["boot"])
cls = ent["classes"]
rng = np.random.default_rng(3)
for b in rng.choice(len(B), min(NB, len(B)), replace=False):
    Vb = dict(zip(cls, B[b]))
    qb, gsb = ground_q(Vb)
    cs = scan(Vb, 8, np.linspace(4000, 200, 39), qb, nsw=600, seed=int(b))
    boot.append({"V": Vb, "q": qb, "Tpeak": cs[int(np.argmax([r["C"] for r in cs]))]["T"]})
Tb = np.array([x["Tpeak"] for x in boot])
# finite-size correction factor from the central model: T(L=16)/T(L=8)
fs = Tpk[16] / Tpk[8]
res = {"key": key, "classes": cls, "V": V, "cv": ent["cv"], "q_ground": q0, "ground": gs, "Tpeak": Tpk, "fs_ratio_16_8": fs,
       "T_order_K": Tpk[16], "boot_Tpeak_L8": Tb.tolist(), "boot_T16_mean": float(np.mean(Tb) * fs), "boot_T16_std": float(np.std(Tb) * fs),
       "boot_T16_p16_p84": [float(np.percentile(Tb, 16) * fs), float(np.percentile(Tb, 84) * fs)], "boot_q": [str(x["q"]) for x in boot],
       "coarse": coarse, "fine": {str(k): v for k, v in fine.items()}}
json.dump(res, open(f"{HERE}/torder_{cej['tag']}_{key}{'_' + csel if csel != 'best' else ''}.json", "w"), indent=1, default=str)
print("T_order(C peak, L=16) = %.0f K; bootstrap (L=8 x fs %.3f): %.0f +- %.0f K (16-84%%: %.0f-%.0f)" % (
    Tpk[16], fs, res["boot_T16_mean"], res["boot_T16_std"], *res["boot_T16_p16_p84"]))
