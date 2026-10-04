"""Direction-resolved exchange fit on the 110-Ry P4/n geometry (90 Ry energies): 24 configs (lcm/jfit110) + 1x1x2 G/Yflip/Baflip (lcm/jy U44).
Fits with all configs and with 'clean' configs only (|M| <= 10.5 muB per cell-normalized, i.e. excluding FM-like states)."""
import json, sys, numpy as np
sys.path.insert(0, "infra/lib"); sys.path.insert(0, "tracks/lcm")
import magtools as mt, jlib
from pymatgen.core import Structure, Lattice
idx = json.load(open("tracks/lcm/jfit110/index.json"))
P = Structure(Lattice(idx["prim"]["lattice"]), idx["prim"]["species"], idx["prim"]["frac"])
res = json.load(open("tracks/lcm/jfit110/results.json"))
st = json.load(open("tracks/lcm/jfit110/stack_U44_U64.json"))["U44"]
st_in = json.load(open("tracks/lcm/conv/ybmfo110_stack.json"))
C = []
for o in idx["configs"]:
    r = res.get(o["id"])
    if r and r["energy_eV"] and r["converged"] and r["n_flipped"] == 0:
        n = round(o["natoms"] / len(P)); C.append({"matrix": o["matrix"], "spins": o["spins"], "E": r["energy_eV"], "n": n, "M": abs(r.get("total_mag") or 0) / n, "tag": o["id"]})
for cf in st_in["configs"]:
    r = st.get(cf["name"])
    if r and r["energy_eV"] and r["converged"] and r["n_flipped"] == 0:
        C.append({"matrix": [[1,0,0],[0,1,0],[0,0,2]], "spins": cf["spins"], "E": r["energy_eV"], "n": 2, "M": abs(r.get("total_mag") or 0) / 2, "tag": "st_" + cf["name"]})
seen, U = set(), []
for c in C:
    k = (json.dumps(c["matrix"]), tuple(c["spins"]))
    if k not in seen: seen.add(k); U.append(c)
print(len(U), "unique configs; per-cell |M| values:", sorted(set(round(c["M"], 1) for c in U)))
# discover pair classes on THIS geometry (rounded-distance names depend on the cell)
found = {}
for i, s in enumerate(P):
    if s.specie.symbol not in ("Mn", "Fe"): continue
    for nb in P.get_neighbors(s, 7.0):
        if P[nb.index].specie.symbol not in ("Mn", "Fe"): continue
        k = jlib.cls(s.specie.symbol, P[nb.index].specie.symbol, nb.nn_distance, abs(nb.coords[2] - s.coords[2]))
        if k: found[k] = found.get(k, 0) + 1
CL = ["ip", "ap", "Y"] + sorted([k for k in found if k not in ("ip", "ap", "Y")], key=lambda k: (float(k.split("_")[-1]), k))
print("classes:", CL)
F = []
for c in U:
    sc = P.copy(); sc.make_supercell(c["matrix"]); F.append(jlib.features(sc, c["spins"], CL))
F = np.array(F); E = np.array([c["E"] for c in U]); N = np.array([c["n"] for c in U]); M = np.array([c["M"] for c in U])
out = {}
for label, mask in [("all", np.ones(len(U), bool)), ("clean_M<=10.5", M <= 10.5)]:
    print(f"--- {label}: {mask.sum()} configs")
    for n in range(3, len(CL) + 1):
        fit = mt.fit_J(F[mask][:, :n], E[mask], N[mask])
        print(n, dict(zip(CL[:n], [round(x, 2) for x in fit["J_meV"]])), "rms %.2f loo %.2f" % (fit["rms_meV_cell"], fit["loo_rms_meV_cell"] or -1))
        out[f"{label}_{n}"] = {"classes": CL[:n], "J_meV": fit["J_meV"], "rms": fit["rms_meV_cell"], "loo": fit["loo_rms_meV_cell"], "nconf": int(mask.sum())}
json.dump(out, open("tracks/lcm/jfit110/fits.json", "w"), indent=1)
