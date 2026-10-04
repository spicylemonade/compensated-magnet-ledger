"""PM cluster expansion for B-site order (lcm_ord2).
usage: python analyze.py TAG parent.json tasks.json [--noFM] [--mc]
 1. load results_<TAG>/*.json; valence check from ortho-atomic Hubbard occupations (moment m = n_up - n_dn, n_d)
 2. paramagnetic energy per arrangement:
      E_avg  = mean of the random configs r0..r3
      E0_fit = intercept of the joint Heisenberg fit  E(a,c) = E0(a) - sum_k J_k C_k(a,c)
               (C_k = sum over bonds of class k of s_i s_j, species-resolved AA/AB/BB x ip1/ap/Y + NNN ip2/apd/Yd)
 3. pair cluster expansion of E_PM per f.u. (sigma = +1 A, -1 Fe): E = c0 + sum_k V_k * (bonds_k per f.u.) * Pi_k,
    model selection by leave-one-out CV over nested class sets; bootstrap uncertainties.
Writes ce_<TAG>.json."""
import json, sys, glob, itertools, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ord2lib import supercell, pair_table, features, spin_features, CLASSES

TAG, PAR, TASKS = sys.argv[1], sys.argv[2], sys.argv[3]
NOFM = "--noFM" in sys.argv
ACCMAX = 2e-3   # max 'estimated scf accuracy' (Ry/cell) for accepting an unconverged SCF energy
NUNC = []
NMV, DMV = [], []
HERE = os.path.dirname(os.path.abspath(__file__))
parent = json.load(open(PAR)); parent = parent.get("relaxed", parent)
tasks = {t["id"]: t for t in json.load(open(TASKS))["tasks"]}
A_EL = "Cu" if "Cu" in json.dumps(parent["species"]) else "Mn"


def valence(el, nu, nd):
    m, n = abs(nu - nd), nu + nd
    if el == "Mn":
        return "Mn2+" if m > 4.25 else ("Mn3+" if m > 3.3 else "Mn?")
    if el == "Fe":
        return "Fe3+" if m > 3.85 else ("Fe2+" if m > 3.3 else "Fe?")
    if el == "Cu":
        return "Cu2+" if 0.4 < m < 1.0 else ("Cu+" if m < 0.2 else "Cu?")
    return el


SPC = ["ip1_AA", "ip1_AB", "ip1_BB", "ap_AA", "ap_AB", "ap_BB", "Y_AA", "Y_AB", "Y_BB"]
NNN = ["ip2", "apd", "Yd"]

rows, arr = [], {}
bondcache = {}
for f in sorted(glob.glob(f"{HERE}/results_{TAG}/*.json")):
    d = json.load(open(f))
    t = tasks.get(d["id"])
    if t is None:
        continue
    key = json.dumps(t["cell"])
    if key not in bondcache:
        bondcache[key] = pair_table(supercell(parent, t["cell"]))
    bonds = bondcache[key]
    sig = np.array(t["sigma"])
    nfu = t["nB"] / 2
    S, N = features(sig, bonds)
    a = {"id": d["id"], "nfu": nfu, "Pi": {k: S[k] / N[k] for k in CLASSES}, "bpf": {k: N[k] / nfu for k in CLASSES}, "runs": {}, "valence": {}}
    for c in t["configs"]:
        r = d.get("runs", {}).get(c["name"])
        if not r:
            continue
        E = r.get("energy_eV") if (r.get("energy_eV") is not None and r.get("converged_scf")) else None
        s1 = r.get("stage1_mv") or {}
        if E is None and s1.get("energy_eV") is not None and s1.get("converged_scf"):
            E = s1["energy_eV"]  # two-stage SCF: gaussian-0.005 restart did not converge -> m-v 0.02 Ry energy (flagged)
            NMV.append((d["id"], c["name"]))
        if E is not None and s1.get("energy_eV") is not None and s1.get("converged_scf") and r.get("converged_scf"):
            DMV.append(1000 * (r["energy_eV"] - s1["energy_eV"]))
        if E is None and r.get("energy_eV_unconv") is not None and (r.get("last_acc_Ry") or 1.0) < ACCMAX:
            E = r["energy_eV_unconv"]  # unconverged within 150 iterations but accuracy < ACCMAX Ry: flagged, kept
            NUNC.append((d["id"], c["name"], r.get("last_acc_Ry")))
        if E is None:
            continue
        bs = np.array(c["bspins"])
        sf = spin_features(bs, sig, bonds, classes=("ip1", "ap", "Y"))
        sf2 = {k: features(bs, bonds)[0][k] for k in NNN}
        x = [sf[k] for k in SPC] + [sf2[k] for k in NNN]
        # valence of each B site
        vl = []
        for b, ai in enumerate(t["bidx"]):
            h = r.get("hub", {}).get(str(ai))
            el = t["species"][ai]
            vl.append(valence(el, *h) if h else "?")
        cnt = {v: vl.count(v) for v in set(vl)}
        a["valence"][c["name"]] = cnt
        a["runs"][c["name"]] = {"E": E, "x": x, "gap": r.get("gap"), "M": r.get("total_mag"), "flip": r.get("n_flipped")}
        rows.append((d["id"], c["name"], E, x))
    if a["runs"]:
        arr[d["id"]] = a
    rr = d.get("relax", {})
    a["relax"] = {k: rr.get(k) for k in ("energy_eV", "bfgs_converged", "pressure_kbar", "max_disp_A", "n_scf_cycles", "total_force")}

MINC = 4  # arrangements need >= MINC configurations (incl. >= 2 random) for the joint fit
ids = sorted(a for a in arr if len(arr[a]["runs"]) >= MINC and sum(1 for c in arr[a]["runs"] if c.startswith("r")) >= 2)
print("incomplete (excluded):", {a: sorted(arr[a]["runs"]) for a in arr if a not in ids})
arr = {a: arr[a] for a in ids}
rows = [r for r in rows if r[0] in arr]
print(f"{TAG}: {len(ids)} arrangements, {len(rows)} energies; unconverged-but-accepted: {NUNC}; m-v-only: {NMV}")
if DMV:
    print("  E(gauss 0.005) - E(m-v 0.02) over %d two-stage runs: mean %.1f, max |.| %.1f meV/cell" % (len(DMV), np.mean(DMV), np.max(np.abs(DMV))))
# ---------------- joint Heisenberg fit ----------------
use = [r for r in rows if not (NOFM and r[1] == "FM")]
ia = {a: i for i, a in enumerate(ids)}
nJ = len(SPC) + len(NNN)
X = np.zeros((len(use), len(ids) + nJ)); y = np.zeros(len(use))
for n, (aid, cn, E, x) in enumerate(use):
    X[n, ia[aid]] = 1.0
    X[n, len(ids):] = -np.array(x)
    y[n] = E
# drop J columns with no data
colok = np.abs(X).sum(0) > 0
coef, *_ = np.linalg.lstsq(X[:, colok], y, rcond=None)
full = np.zeros(X.shape[1]); full[colok] = coef
res = y - X @ full
dof = max(1, len(y) - colok.sum())
s2 = float(res @ res / dof)
cov = s2 * np.linalg.pinv(X[:, colok].T @ X[:, colok])
err = np.zeros(X.shape[1]); err[colok] = np.sqrt(np.diag(cov))
J = {k: (1000 * full[len(ids) + i], 1000 * err[len(ids) + i]) for i, k in enumerate(SPC + NNN)}
print("joint fit: rms %.1f meV/cell over %d energies; J (meV, E = E0 - J s.s, unit spins):" % (1000 * np.sqrt(np.mean(res ** 2)), len(y)))
print("  " + ", ".join(f"{k} {v[0]:.1f}+-{v[1]:.1f}" for k, v in J.items()))
for aid in ids:
    a = arr[aid]
    a["E0_fit"] = full[ia[aid]]; a["E0_err"] = err[ia[aid]]
    rnd = [a["runs"][c]["E"] for c in a["runs"] if c.startswith("r")]
    a["E_avg"] = float(np.mean(rnd)) if rnd else None
    a["n_rand"] = len(rnd)
    a["resid_meV"] = {cn: round(1000 * float(res[n]), 1) for n, (aid2, cn, E, x) in enumerate(use) if aid2 == aid}
# reference: rock-salt (YBMFO) or the experimental layered order (YBCFO) -> use RS_A as common reference if present
ref = f"{TAG}_RS_A"
e_ref_fit = arr[ref]["E0_fit"] / arr[ref]["nfu"]
e_ref_avg = arr[ref]["E_avg"] / arr[ref]["nfu"]
tab = []
for aid in ids:
    a = arr[aid]
    eG = a["runs"].get("G", {}).get("E"); eF = a["runs"].get("FM", {}).get("E")
    row = {"id": aid, "nfu": a["nfu"], "Pi": a["Pi"], "n_rand": a["n_rand"],
           "dE_PM_fit": 1000 * (a["E0_fit"] / a["nfu"] - e_ref_fit), "dE_PM_err": 1000 * a["E0_err"] / a["nfu"],
           "dE_PM_avg": 1000 * (a["E_avg"] / a["nfu"] - e_ref_avg) if a["E_avg"] is not None else None,
           "dE_G": 1000 * (eG / a["nfu"] - arr[ref]["runs"]["G"]["E"] / arr[ref]["nfu"]) if eG and "G" in arr[ref]["runs"] else None,
           "dE_FM": 1000 * (eF / a["nfu"] - arr[ref]["runs"]["FM"]["E"] / arr[ref]["nfu"]) if eF and "FM" in arr[ref]["runs"] else None,
           "valence": a["valence"], "gaps": {c: (round(v["gap"], 2) if v.get("gap") is not None else None) for c, v in a["runs"].items()},
           "resid_meV": a["resid_meV"], "relax": a["relax"]}
    tab.append(row)
print("\n%-18s %6s %8s %8s %8s %8s  valence(r-configs)" % ("arrangement", "nrand", "PM_fit", "PM_avg", "G", "FM"))
for r in sorted(tab, key=lambda r: r["dE_PM_fit"]):
    vv = {}
    for c, cnt in r["valence"].items():
        if c.startswith("r"):
            for k, n in cnt.items():
                vv[k] = vv.get(k, 0) + n
    f = lambda x: "%8.1f" % x if x is not None else "%8s" % "-"
    print("%-18s %6d %s %s %s %s  %s" % (r["id"], r["n_rand"], f(r["dE_PM_fit"]), f(r["dE_PM_avg"]), f(r["dE_G"]), f(r["dE_FM"]), vv))

# ---------------- cluster expansion ----------------
def design(rows_, cls):
    Xc = np.array([[1.0] + [r["Pi"][k] * arrbpf[r["id"]][k] for k in cls] for r in rows_])
    return Xc


arrbpf = {aid: arr[aid]["bpf"] for aid in ids}
SETS = [["ip1", "ap", "Y"], ["ip1", "ap", "Y", "ip2"], ["ip1", "ap", "Y", "ip2", "apd", "Yd"], ["ip1", "ap", "Y", "ip2", "ip3"],
        ["ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3"], ["ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3", "c2"], ["ip1", "ap", "Y", "apd", "Yd"],
        ["ip1", "ap", "Y", "c2"], ["ip1", "ap", "Y", "ip2", "c2"]]
ce = {}
for key in ("dE_PM_fit", "dE_PM_avg", "dE_G", "dE_FM"):
    data = [r for r in tab if r[key] is not None]
    if len(data) < 6:
        continue
    yv = np.array([r[key] for r in data])
    best = None
    for cls in SETS:
        if len(cls) + 1 >= len(data) - 1:
            continue
        Xc = design(data, cls)
        c, *_ = np.linalg.lstsq(Xc, yv, rcond=None)
        rms = float(np.sqrt(np.mean((yv - Xc @ c) ** 2)))
        loo = []
        for i in range(len(data)):
            m = np.ones(len(data), bool); m[i] = False
            ci, *_ = np.linalg.lstsq(Xc[m], yv[m], rcond=None)
            loo.append(yv[i] - Xc[i] @ ci)
        cv = float(np.sqrt(np.mean(np.square(loo))))
        # bootstrap
        rng = np.random.default_rng(0); bs = []
        for _ in range(400):
            idx = rng.integers(0, len(data), len(data))
            cb, *_ = np.linalg.lstsq(Xc[idx], yv[idx], rcond=None)
            bs.append(cb)
        bs = np.array(bs)
        entry = {"classes": cls, "V_meV_per_bond": dict(zip(cls, c[1:].tolist())), "c0": float(c[0]), "rms": rms, "cv": cv,
                 "V_boot_std": dict(zip(cls, bs[:, 1:].std(0).tolist())), "boot": bs[:, 1:].tolist(), "n": len(data),
                 "loo_resid": {r["id"]: float(l) for r, l in zip(data, loo)}}
        ce.setdefault(key, []).append(entry)
        if best is None or cv < best["cv"]:
            best = entry
    print(f"\nCE on {key} (meV/f.u., n={len(data)}):")
    for e in ce[key]:
        print("  %-40s rms %6.1f  CV %6.1f  V: %s" % ("+".join(e["classes"]), e["rms"], e["cv"],
                                                     ", ".join(f"{k} {v:.1f}+-{e['V_boot_std'][k]:.1f}" for k, v in e["V_meV_per_bond"].items())))
    ce[key + "_best"] = best
out = {"tag": TAG, "noFM": NOFM, "J_joint_meV": J, "joint_rms_meV": float(1000 * np.sqrt(np.mean(res ** 2))), "table": tab, "ce": ce}
json.dump(out, open(f"{HERE}/ce_{TAG}{'_noFM' if NOFM else ''}.json", "w"), indent=1, default=float)
