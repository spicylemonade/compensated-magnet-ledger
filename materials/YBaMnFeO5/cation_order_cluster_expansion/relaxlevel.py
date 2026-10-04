"""Relax-level PM energy estimates (cross-check + provisional values for arrangements whose protocol SCFs are not done).
The ions-only relax runs in spin config r0 at 70 Ry / half k mesh / m-v 0.02 Ry. Its final (or, if still running, latest) '!'
energy, spin-corrected with the clean exchange fit of ce2_<TAG>.json (E0 = E + x_r0.J), is compared with the protocol r0 E0.
Cell-dependent offsets (half k meshes are not k-equivalent between cell shapes) are calibrated on clean arrangements that have
both numbers (relative to RS_A); offsets are applied per cell type when available.
usage: python relaxlevel.py TAG parent.json tasks.json [RUNNING_ID ...]   -> relaxlevel_<TAG>.json"""
import json, sys, os, re, glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ord2lib import supercell, pair_table, features, spin_features

HERE = os.path.dirname(os.path.abspath(__file__))
TAG, PAR, TASKS = sys.argv[1], sys.argv[2], sys.argv[3]
running = sys.argv[4:]
parent = json.load(open(PAR)); parent = parent.get("relaxed", parent)
tasks = {t["id"]: t for t in json.load(open(TASKS))["tasks"]}
ce = json.load(open(f"{HERE}/ce2_{TAG}.json"))
J = np.array([v[0] for v in ce["J_clean_noFM_meV"].values()]) / 1000
SPC = ["ip1_AA", "ip1_AB", "ip1_BB", "ap_AA", "ap_AB", "ap_BB", "Y_AA", "Y_AB", "Y_BB"]; NNN = ["ip2", "apd", "Yd"]
CELLNAME = {json.dumps(v): k for k, v in {"A": [[1, 0, 0], [0, 1, 0], [0, 0, 2]], "B": [[2, 0, 0], [0, 1, 0], [0, 0, 1]],
                                          "C": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]], "D": [[2, 0, 0], [0, 2, 0], [0, 0, 1]],
                                          "F": [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]}.items()}


def xr0(t):
    bonds = pair_table(supercell(parent, t["cell"])); sig = np.array(t["sigma"])
    c = [c for c in t["configs"] if c["name"] == t.get("relax_config", "r0")][0]; bs = np.array(c["bspins"])
    sf = spin_features(bs, sig, bonds, classes=("ip1", "ap", "Y")); S2 = features(bs, bonds)[0]
    return np.array([sf[k] for k in SPC] + [S2[k] for k in NNN])


out = {}
for f in sorted(glob.glob(f"{HERE}/results_{TAG}/*.json")):
    d = json.load(open(f)); t = tasks[d["id"]]; rel = d.get("relax") or {}
    if rel.get("energy_eV") is None:
        continue
    x = xr0(t); nfu = t["nB"] / 2
    r0 = d.get("runs", {}).get("r0") or {}
    out[d["id"]] = {"cell": CELLNAME[json.dumps(t["cell"])], "E_relax": (rel["energy_eV"] + x @ J) / nfu,
                    "E_prot": (r0["energy_eV"] + x @ J) / nfu if r0.get("energy_eV") is not None and r0.get("converged_scf") else None,
                    "relaxed": True, "bfgs": rel.get("bfgs_converged"), "F": rel.get("total_force"),
                    "ct_r0": None}
# still-relaxing jobs: parse the live relax output (in memory only)
if running:
    import modal
    v = modal.Volume.from_name("magdisc-data")
    for tid in running:
        t = tasks[tid]; nfu = t["nB"] / 2
        try:
            txt = b"".join(v.read_file(f"jobs/lcm/ord2/v4/{tid}/{tid}_relax.out")).decode(errors="replace")
        except Exception as ex:
            print(tid, "no relax out", ex); continue
        en = re.findall(r"^!\s+total energy\s+=\s+([-0-9.]+)\s+Ry", txt, re.M)
        fo = re.findall(r"Total force =\s+([0-9.]+)", txt)
        if not en:
            print(tid, "no ! energy yet"); continue
        x = xr0(t)
        out[tid] = {"cell": CELLNAME[json.dumps(t["cell"])], "E_relax": (float(en[-1]) * 13.605693122994 + x @ J) / nfu, "E_prot": None,
                    "relaxed": False, "bfgs": False, "F": float(fo[-1]) if fo else None, "nsteps_current_run": len(en),
                    "E_first_last_meV_fu": 1000 * (float(en[0]) - float(en[-1])) * 13.605693122994 / nfu}
ref = out[f"{TAG}_RS_A"]
for k, o in out.items():
    o["rel_relax"] = 1000 * (o["E_relax"] - ref["E_relax"])
    o["rel_prot"] = 1000 * (o["E_prot"] - ref["E_prot"]) if o["E_prot"] is not None else None
# offsets per cell from clean arrangements (protocol low-branch clean, no CT)
arrs = ce["arrangements"]
off = {}
for cell in sorted(set(o["cell"] for o in out.values())):
    dd = [o["rel_relax"] - o["rel_prot"] for k, o in out.items() if o["cell"] == cell and o["rel_prot"] is not None
          and arrs.get(k, {}).get("clean") and arrs[k]["clean"]["n"] == arrs[k]["n"]]
    if dd:
        off[cell] = (float(np.mean(dd)), float(np.std(dd)), len(dd))
print("cell offsets relax-level minus protocol (meV/f.u., clean arrangements):", {c: tuple(round(z, 1) for z in v) for c, v in off.items()})
for k in sorted(out, key=lambda k: out[k]["rel_relax"]):
    o = out[k]
    oc = off.get(o["cell"], (None,))[0]
    o["rel_relax_corr"] = o["rel_relax"] - oc if oc is not None else None
    print("%-16s cell %s relax %7.1f corr %s protocol(r0) %s %s F=%s" % (k, o["cell"], o["rel_relax"],
          "%7.1f" % o["rel_relax_corr"] if o["rel_relax_corr"] is not None else "    n/a",
          "%7.1f" % o["rel_prot"] if o["rel_prot"] is not None else "    n/a", "" if o["relaxed"] else "(RELAX RUNNING)", o["F"]))
json.dump({"offsets": off, "arr": out}, open(f"{HERE}/relaxlevel_{TAG}.json", "w"), indent=1, default=float)
