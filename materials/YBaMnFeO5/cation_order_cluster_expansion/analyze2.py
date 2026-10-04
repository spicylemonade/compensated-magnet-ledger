"""PM cluster expansion for B-site order, v2 analysis (2026-10-02 pass; supersedes analyze.py's joint fit, which mixed in
FM and charge-transferred (CT) configurations and gave rms 315 meV/cell).
usage: python analyze2.py TAG parent.json tasks.json        -> ce2_<TAG>.json (+ printed tables)

1. Per run: protocol energy (stage-2 gaussian 0.005, converged only), species-resolved spin-correlation features, valence of
   every B site from the ortho-atomic Hubbard occupations; a run is 'CT' if any site is Mn3+/Fe2+/Cu+.
2. Exchange model: joint fit E = E0(arr) - sum_k J_k C_k on CLEAN (no CT), non-FM runs (12 J: ip1/ap/Y x AA/AB/BB + ip2/apd/Yd).
   FM is excluded everywhere (half-metallic; sits +35..50 meV/f.u. above the Heisenberg line in every clean arrangement).
3. Per-run PM estimate E0_i = E_i + x_i.J; per arrangement:
     low  = mean of E0 over non-FM runs within W = 60 meV/f.u. of the arrangement minimum (drops metastable CT solutions that do
            not match the relaxed geometry; for clean arrangements this is all runs)            [main]
     all  = mean over all non-FM runs (includes the metastable CT solutions)                     [upper]
     clean= mean over clean non-FM runs (arrangements without clean runs dropped)                [ionic d5/d5 landscape]
   sigma_arr^2 = 6^2 (geometry noise: RS_A vs RS_C) + s^2/n (s = spread, 25 meV/f.u. if n = 1).
4. Pair CE  E/f.u. = c0 + sum_k V_k (bonds_k per f.u.) Pi_k  (V per bond, sigma = +1 Mn/Cu, -1 Fe), weighted least squares,
   all supersets of {ip1, ap, Y} with up to 5 of {ip2, apd, Yd, ip3, c2}; weighted LOO-CV selection; 400 bootstrap refits.
"""
import json, sys, glob, itertools, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ord2lib import supercell, pair_table, features, spin_features, CLASSES

TAG, PAR, TASKS = sys.argv[1], sys.argv[2], sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
W = 60.0       # meV/f.u. window for the low branch
SGEO = 6.0     # meV/f.u.
S1 = 25.0      # meV/f.u. spread assumed for single-run arrangements
parent = json.load(open(PAR)); parent = parent.get("relaxed", parent)
tasks = {t["id"]: t for t in json.load(open(TASKS))["tasks"]}
SPC = ["ip1_AA", "ip1_AB", "ip1_BB", "ap_AA", "ap_AB", "ap_BB", "Y_AA", "Y_AB", "Y_BB"]
NNN = ["ip2", "apd", "Yd"]


def valence(el, nu, nd):
    m = abs(nu - nd)
    if el == "Mn":
        return "Mn2+" if m > 4.25 else ("Mn3+" if m > 3.3 else "Mn?")
    if el == "Fe":
        return "Fe3+" if m > 3.85 else ("Fe2+" if m > 3.3 else "Fe?")
    if el == "Cu":
        return "Cu2+" if 0.4 < m < 1.0 else ("Cu+" if m < 0.2 else "Cu?")
    return el


rows, info, bc = [], {}, {}
for f in sorted(glob.glob(f"{HERE}/results_{TAG}/*.json")):
    d = json.load(open(f)); t = tasks.get(d["id"])
    if t is None:
        continue
    key = json.dumps(t["cell"])
    if key not in bc:
        bc[key] = pair_table(supercell(parent, t["cell"]))
    bonds = bc[key]; sig = np.array(t["sigma"]); nfu = t["nB"] / 2
    S, N = features(sig, bonds)
    info[d["id"]] = {"nfu": nfu, "Pi": {k: S[k] / N[k] for k in CLASSES}, "bpf": {k: N[k] / nfu for k in CLASSES},
                     "relax": {k: (d.get("relax") or {}).get(k) for k in ("bfgs_converged", "total_force", "max_disp_A", "n_scf_cycles")},
                     "runs_done": sorted(k for k, r in d.get("runs", {}).items() if r.get("energy_eV") is not None)}
    for c in t["configs"]:
        r = d.get("runs", {}).get(c["name"])
        if not r or r.get("energy_eV") is None or not r.get("converged_scf"):
            continue
        bs = np.array(c["bspins"])
        sf = spin_features(bs, sig, bonds, classes=("ip1", "ap", "Y")); S2 = features(bs, bonds)[0]
        x = np.array([sf[k] for k in SPC] + [S2[k] for k in NNN])
        vl = [valence(t["species"][ai], *r["hub"][str(ai)]) if str(ai) in r.get("hub", {}) else "?" for ai in t["bidx"]]
        ct = sum(1 for v in vl if v in ("Mn3+", "Fe2+", "Cu+", "Mn?", "Fe?", "Cu?"))
        s1 = r.get("stage1_mv") or {}
        rows.append({"id": d["id"], "c": c["name"], "E": r["energy_eV"], "x": x, "nfu": nfu, "ct": ct,
                     "val": {v: vl.count(v) for v in sorted(set(vl))}, "gap": r.get("gap"),
                     "acc": r.get("last_acc_Ry"), "dmv_meV": 1000 * (r["energy_eV"] - s1["energy_eV"]) if s1.get("energy_eV") else None})

# ---------------- exchange fit on clean non-FM runs ----------------
fitrows = [r for r in rows if r["ct"] == 0 and r["c"] != "FM"]
ids_f = sorted(set(r["id"] for r in fitrows)); ia = {a: i for i, a in enumerate(ids_f)}
X = np.zeros((len(fitrows), len(ids_f) + 12)); y = np.zeros(len(fitrows))
for n, r in enumerate(fitrows):
    X[n, ia[r["id"]]] = 1.0; X[n, len(ids_f):] = -r["x"]; y[n] = r["E"]
ok = np.abs(X).sum(0) > 0
coef, *_ = np.linalg.lstsq(X[:, ok], y, rcond=None)
full = np.zeros(X.shape[1]); full[ok] = coef
res = y - X @ full
dof = max(1, len(y) - ok.sum())
cov = float(res @ res / dof) * np.linalg.pinv(X[:, ok].T @ X[:, ok])
err = np.zeros(X.shape[1]); err[ok] = np.sqrt(np.diag(cov))
J = full[len(ids_f):]
Jtab = {k: [round(1000 * J[i], 2), round(1000 * err[len(ids_f) + i], 2)] for i, k in enumerate(SPC + NNN)}
jrms = 1000 * float(np.sqrt(np.mean(res ** 2)))
print(f"{TAG}: {len(rows)} converged runs ({sum(1 for r in rows if r['ct'])} CT, {sum(1 for r in rows if r['c'] == 'FM')} FM)")
print("exchange fit (clean, non-FM): n=%d, %d arrangements, rms %.1f meV/cell" % (len(y), len(ids_f), jrms))
print("  J (meV, unit spins):", ", ".join(f"{k} {v[0]:.1f}+-{v[1]:.1f}" for k, v in Jtab.items()))

# ---------------- per-arrangement PM energies ----------------
ref = f"{TAG}_RS_A"
for r in rows:
    r["E0"] = (r["E"] + r["x"] @ J) / r["nfu"]
    r["fm_dev"] = None
e_ref = np.mean([r["E0"] for r in rows if r["id"] == ref and r["c"] != "FM"])
for r in rows:
    r["E0rel"] = 1000 * (r["E0"] - e_ref)
fmdev = []
arr = {}
for a in sorted(info):
    R = [r for r in rows if r["id"] == a and r["c"] != "FM"]
    if not R:
        continue
    e = np.array([r["E0rel"] for r in R])
    lowm = e <= e.min() + W
    cl = np.array([r["ct"] == 0 for r in R])
    ent = {"n": len(R), "configs": {r["c"]: {"E0": round(r["E0rel"], 1), "ct_sites": r["ct"], "val": r["val"],
                                              "gap": None if r["gap"] is None else round(r["gap"], 2)} for r in R}}
    for lab, m in (("low", lowm), ("all", np.ones(len(R), bool)), ("clean", cl)):
        if m.sum() == 0:
            ent[lab] = None; continue
        s = float(e[m].std(ddof=1)) if m.sum() > 1 else S1
        ent[lab] = {"E": float(e[m].mean()), "n": int(m.sum()), "spread": s,
                    "sigma": float(np.sqrt(SGEO ** 2 + s ** 2 / m.sum())), "ct": bool(any(R[i]["ct"] for i in np.where(m)[0]))}
    fm = [r for r in rows if r["id"] == a and r["c"] == "FM"]
    if fm and ent["clean"] and not fm[0]["ct"]:
        fmdev.append(fm[0]["E0rel"] - ent["clean"]["E"])
    ent.update(info[a]); arr[a] = ent
print("FM run above the clean PM line (Heisenberg violation), meV/f.u.: mean %.0f, range %.0f..%.0f over %d arrangements" % (
    np.mean(fmdev), np.min(fmdev), np.max(fmdev), len(fmdev)) if fmdev else "no FM check")
print("\n%-16s %4s %8s %6s %8s %8s  %-28s %s" % ("arrangement", "n", "low", "sig", "all", "clean", "Pi(ip1,ap,Y,ip2,apd,Yd,ip3,c2)", "CT runs"))
for a in sorted(arr, key=lambda a: arr[a]["low"]["E"]):
    e = arr[a]
    if not e.get("low"):
        continue
    f = lambda z: "%8.1f" % z["E"] if z else "%8s" % "-"
    pi = ",".join("%+.2g" % e["Pi"][k] for k in ("ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3", "c2"))
    print("%-16s %4d %s %6.1f %s %s  %-28s %s" % (a, e["n"], f(e["low"]), e["low"]["sigma"], f(e["all"]), f(e["clean"]), pi,
                                              ",".join(c for c, v in e["configs"].items() if v["ct_sites"]) or "-"))

# ---------------- best estimate: protocol low branch vs relax-level (variationally lower state wins) ----------------
RL = f"{HERE}/relaxlevel_{TAG}.json"
RLOFF = {"A": 0.0, "B": 0.0, "D": 0.0}          # relax-level offsets vs protocol (meV/f.u.); C/F only if calibrated below
if TAG == "YBMFO":
    RLOFF["C"] = -56.2                           # 4 clean C-cell arrangements, std 3.0
SRL = 15.0                                       # sigma of a relax-level single-config estimate (meV/f.u.)
if os.path.exists(RL):
    rl = json.load(open(RL))["arr"]
    for a, o in rl.items():
        if a not in arr:   # arrangement without any protocol SCF yet: Pi from the task
            t = tasks[a]
            key = json.dumps(t["cell"])
            if key not in bc:
                bc[key] = pair_table(supercell(parent, t["cell"]))
            S, N = features(np.array(t["sigma"]), bc[key]); nfu = t["nB"] / 2
            arr[a] = {"n": 0, "configs": {}, "low": None, "all": None, "clean": None, "nfu": nfu,
                      "Pi": {k: S[k] / N[k] for k in CLASSES}, "bpf": {k: N[k] / nfu for k in CLASSES}, "relax": {}, "runs_done": []}
    for a, e in arr.items():
        o = rl.get(a)
        rlv = None
        if o and o["cell"] in RLOFF:
            rlv = {"E": o["rel_relax"] - RLOFF[o["cell"]], "sigma": SRL if o.get("relaxed") else 25.0,
                   "relaxed": o.get("relaxed"), "cell": o["cell"]}
        e["relaxlevel"] = rlv
        cands = []
        if e["low"]:
            cands.append(("protocol", e["low"]["E"], e["low"]["sigma"]))
        if rlv:
            cands.append(("relax-level" + ("" if rlv["relaxed"] else " (relax running)"), rlv["E"], rlv["sigma"]))
        if not cands:
            e["best"] = None; continue
        src, E, s = min(cands, key=lambda c: c[1])
        if src == "relax-level" and e["low"] and e["low"]["E"] - E < 30:   # within noise: keep the protocol number
            src, E, s = cands[0]
        e["best"] = {"E": E, "sigma": s, "src": src, "ct": bool(e["low"] and e["low"]["ct"])}
    print("\nbest estimates (protocol low branch unless the relax-level state is >30 meV/f.u. lower or no protocol SCF yet):")
    for a in sorted(arr, key=lambda a: arr[a]["best"]["E"] if arr[a].get("best") else 1e9):
        b = arr[a].get("best")
        if b:
            print("  %-16s %7.1f +- %4.1f  %s" % (a, b["E"], b["sigma"], b["src"]))

# ---------------- cluster expansion ----------------
OPT = ["ip2", "apd", "Yd", "ip3", "c2"]
SETS = [["ip1", "ap", "Y"] + list(s) for n in range(0, 6) for s in itertools.combinations(OPT, n)]


KTW = 100.0   # meV: Boltzmann weight exp(-(E - Emin)/KTW) for the 'best_bw' dataset (low-energy structures dominate T_order)


def ce_fit(data, cls, nboot=400):
    yv = np.array([d["E"] for d in data]); wv = 1.0 / np.array([d["sigma"] for d in data]) ** 2
    wv = wv * np.array([d.get("bw", 1.0) for d in data])
    Xc = np.array([[1.0] + [d["Pi"][k] * d["bpf"][k] for k in cls] for d in data])
    sw = np.sqrt(wv)
    c, *_ = np.linalg.lstsq(Xc * sw[:, None], yv * sw, rcond=None)
    resid = yv - Xc @ c
    loo = []
    for i in range(len(data)):
        m = np.ones(len(data), bool); m[i] = False
        ci, *_ = np.linalg.lstsq(Xc[m] * sw[m, None], yv[m] * sw[m], rcond=None)
        loo.append(yv[i] - Xc[i] @ ci)
    loo = np.array(loo)
    cvw = float(np.sqrt(np.sum(wv * loo ** 2) / np.sum(wv)))
    cv = float(np.sqrt(np.mean(loo ** 2)))
    rng = np.random.default_rng(0); bs = []
    for _ in range(nboot):
        idx = rng.integers(0, len(data), len(data))
        if np.linalg.matrix_rank(Xc[idx]) < Xc.shape[1]:
            continue
        cb, *_ = np.linalg.lstsq(Xc[idx] * sw[idx, None], yv[idx] * sw[idx], rcond=None)
        bs.append(cb)
    bs = np.array(bs)
    return {"classes": cls, "V_meV_per_bond": dict(zip(cls, c[1:].tolist())), "c0": float(c[0]),
            "rms": float(np.sqrt(np.mean(resid ** 2))), "wrms": float(np.sqrt(np.sum(wv * resid ** 2) / np.sum(wv))), "cv": cvw, "cv_unw": cv,
            "V_boot_std": dict(zip(cls, bs[:, 1:].std(0).tolist())), "boot": bs[:, 1:].tolist(), "n": len(data),
            "resid": {d["id"]: round(float(r), 1) for d, r in zip(data, resid)}, "loo_resid": {d["id"]: round(float(l), 1) for d, l in zip(data, loo)}}


ALLC = ["ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3", "c2"]
PRIOR = {"ip1": 100.0, "ap": 100.0, "Y": 100.0, "ip2": 15.0, "apd": 15.0, "Yd": 15.0, "ip3": 8.0, "c2": 8.0}   # meV per bond


def ce_bayes(data, nboot=400):
    """Bayesian (ridge) pair CE with shell-dependent Gaussian priors V_k ~ N(0, (s*PRIOR_k)^2), all 8 classes, prior scale s
    chosen by weighted LOO-CV (Mueller & Ceder style); c0 unregularised."""
    yv = np.array([d["E"] for d in data]); wv = 1.0 / np.array([d["sigma"] for d in data]) ** 2
    wv = wv * np.array([d.get("bw", 1.0) for d in data])
    Xc = np.array([[1.0] + [d["Pi"][k] * d["bpf"][k] for k in ALLC] for d in data])

    def solve(Xs, ys, ws, sc):
        P = np.diag([0.0] + [1.0 / (sc * PRIOR[k]) ** 2 for k in ALLC])
        A = Xs.T @ (ws[:, None] * Xs) + P
        return np.linalg.solve(A, Xs.T @ (ws * ys))

    best = None
    for sc in (0.25, 0.5, 1.0, 2.0, 4.0):
        loo = np.array([yv[i] - Xc[i] @ solve(np.delete(Xc, i, 0), np.delete(yv, i), np.delete(wv, i), sc) for i in range(len(data))])
        cvw = float(np.sqrt(np.sum(wv * loo ** 2) / np.sum(wv)))
        if best is None or cvw < best[0]:
            best = (cvw, sc, loo)
    cvw, sc, loo = best
    c = solve(Xc, yv, wv, sc)
    resid = yv - Xc @ c
    rng = np.random.default_rng(0); bs = []
    for _ in range(nboot):
        idx = rng.integers(0, len(data), len(data))
        bs.append(solve(Xc[idx], yv[idx], wv[idx], sc))
    bs = np.array(bs)
    return {"classes": ALLC, "V_meV_per_bond": dict(zip(ALLC, c[1:].tolist())), "c0": float(c[0]), "prior_scale": sc,
            "rms": float(np.sqrt(np.mean(resid ** 2))), "wrms": float(np.sqrt(np.sum(wv * resid ** 2) / np.sum(wv))), "cv": cvw,
            "cv_unw": float(np.sqrt(np.mean(loo ** 2))), "V_boot_std": dict(zip(ALLC, bs[:, 1:].std(0).tolist())), "boot": bs[:, 1:].tolist(),
            "n": len(data), "resid": {d["id"]: round(float(r), 1) for d, r in zip(data, resid)},
            "loo_resid": {d["id"]: round(float(l), 1) for d, l in zip(data, loo)}}


ce = {}
for lab in ("best", "best_bw", "low", "all", "clean", "low_noCT"):
    if lab == "best_bw":
        emin = min(arr[a]["best"]["E"] for a in arr if arr[a].get("best"))
        data = [dict(id=a, E=arr[a]["best"]["E"], sigma=arr[a]["best"]["sigma"], Pi=arr[a]["Pi"], bpf=arr[a]["bpf"],
                     bw=float(np.exp(-(arr[a]["best"]["E"] - emin) / KTW))) for a in arr if arr[a].get("best")]
    elif lab == "low_noCT":
        data = [dict(id=a, E=arr[a]["low"]["E"], sigma=arr[a]["low"]["sigma"], Pi=arr[a]["Pi"], bpf=arr[a]["bpf"]) for a in arr
                if arr[a]["low"] and not arr[a]["low"]["ct"]]
    else:
        data = [dict(id=a, E=arr[a][lab]["E"], sigma=arr[a][lab]["sigma"], Pi=arr[a]["Pi"], bpf=arr[a]["bpf"]) for a in arr if arr[a].get(lab)]
    key = "E_PM_" + lab
    ents = []
    for cls in SETS:
        if len(cls) + 1 > len(data) - 2:
            continue
        ents.append(ce_fit(data, cls))
    if not ents:
        continue
    ents.sort(key=lambda e: e["cv"])
    ce[key] = ents; ce[key + "_best"] = ents[0]
    print(f"\nCE on {key} (meV/f.u., n={len(data)}): 5 best by weighted LOO-CV")
    for e in ents[:5]:
        print("  %-34s wrms %5.1f  CVw %5.1f (unw %5.1f)  V: %s" % ("+".join(e["classes"]), e["wrms"], e["cv"], e["cv_unw"],
                                                            ", ".join(f"{k} {v:.1f}+-{e['V_boot_std'][k]:.1f}" for k, v in e["V_meV_per_bond"].items())))
    b = ents[0]
    print("  best LOO residuals:", {k: v for k, v in b["loo_resid"].items() if abs(v) > 25})
    bb = ce_bayes(data)
    ce[key + "_bayes"] = [bb]; ce[key + "_bayes_best"] = bb
    print("  Bayesian CE (all 8 classes, prior scale %.2f): wrms %.1f CVw %.1f (unw %.1f) V: %s" % (bb["prior_scale"], bb["wrms"], bb["cv"], bb["cv_unw"],
          ", ".join(f"{k} {v:.1f}+-{bb['V_boot_std'][k]:.1f}" for k, v in bb["V_meV_per_bond"].items())))
    print("    LOO residuals > 25:", {k: v for k, v in bb["loo_resid"].items() if abs(v) > 25})

# keep bootstrap samples only for the entries used by torder.py (best per dataset, Bayesian, and variants named in run_torder.sh)
KEEP = {"E_PM_best": ["ip1+ap+Y"]}
for key, ents in list(ce.items()):
    if not isinstance(ents, list) or key + "_best" not in ce:
        continue
    for e in ents:
        if e is not ce.get(key + "_best") and "+".join(e["classes"]) not in KEEP.get(key, []):
            e.pop("boot", None)
out = {"tag": TAG, "window_meV_fu": W, "sigma_geo": SGEO, "J_clean_noFM_meV": Jtab, "J_fit_rms_meV_cell": jrms,
       "n_runs": len(rows), "n_ct_runs": sum(1 for r in rows if r["ct"]), "fm_heisenberg_dev_meV_fu": fmdev,
       "arrangements": arr, "ce": ce}
json.dump(out, open(f"{HERE}/ce2_{TAG}.json", "w"), indent=1, default=float)
