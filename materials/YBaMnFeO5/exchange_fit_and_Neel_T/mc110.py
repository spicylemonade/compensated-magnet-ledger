"""Binder-crossing classical MC (and sublattice M(T)) for the 110-Ry-geometry J model of YBaMnFeO5.
usage: python tracks/lcm/mc110.py <fit_key> e.g. clean_M<=10.5_10"""
import json, sys, math, numpy as np
sys.path.insert(0, "infra/lib"); sys.path.insert(0, "tracks/lcm")
import magtools as mt, jlib
from pymatgen.core import Structure, Lattice
key = sys.argv[1]
fit = json.load(open("tracks/lcm/jfit110/fits.json"))[key]
J = dict(zip(fit["classes"], fit["J_meV"]))
idx = json.load(open("tracks/lcm/jfit110/index.json"))
P = Structure(Lattice(idx["prim"]["lattice"]), idx["prim"]["species"], idx["prim"]["frac"])
def lattice(L):
    sc = P.copy(); sc.make_supercell([L, L, max(2, L // 2)])
    ii = [i for i, s in enumerate(sc) if s.specie.symbol in ("Mn", "Fe")]; pos = {j: k for k, j in enumerate(ii)}
    nb = [[] for _ in ii]
    for i in ii:
        for n in sc.get_neighbors(sc[i], 7.0):
            if n.index not in pos: continue
            k = jlib.cls(sc[i].specie.symbol, sc[n.index].specie.symbol, n.nn_distance, abs(n.coords[2] - sc[i].coords[2]))
            if k in J and abs(J[k]) > 1e-9: nb[pos[i]].append((pos[n.index], J[k]))
    m = max(len(x) for x in nb); NB = -np.ones((len(ii), m), dtype=np.int64); JJ = np.zeros((len(ii), m))
    for a, l in enumerate(nb):
        for b, (j, v) in enumerate(l): NB[a, b] = j; JJ[a, b] = v
    return NB, JJ, np.array([1.0 if sc[i].specie.symbol == "Mn" else -1.0 for i in ii])
NB, JJ, pat = lattice(6)
Eg = -0.5 * sum(JJ[a, b] * pat[a] * pat[NB[a, b]] for a in range(len(pat)) for b in range(NB.shape[1]) if NB[a, b] >= 0) / len(pat)
Tmf = -Eg * 2 / 3 / 0.08617333
print("J:", {k: round(v, 2) for k, v in J.items()}, "T_MF %.0f K" % Tmf, flush=True)
temps = np.arange(round(0.5 * Tmf / 10) * 10, round(0.72 * Tmf / 10) * 10 + 1, 10.0)
out = {}
for L in (6, 8, 10):
    NB, JJ, pat = lattice(L)
    sc = mt.mc_scan(NB, JJ, pat, temps, nsweep=12000, nequil=5000, seed=300 + L)
    out[L] = [(s["T"], s["M"], s["U4"], s["chi"], s["C"]) for s in sc]
    print("L", L, [(int(t), round(u, 3)) for t, m, u, c, cc in out[L]], flush=True)
cr = []
for a, b in ((6, 8), (8, 10), (6, 10)):
    ta = np.array([x[0] for x in out[a]]); d = np.array([x[2] for x in out[a]]) - np.array([x[2] for x in out[b]])
    ii = np.where(np.diff(np.sign(d)) != 0)[0]
    c = [ta[i] + (ta[i+1] - ta[i]) * d[i] / (d[i] - d[i+1]) for i in ii]
    print("Binder crossing", a, b, [round(x) for x in c]); cr += c[:1]
json.dump({"key": key, "J": J, "T_MF": Tmf, "scan": out, "binder_first_crossings": cr}, open(f"tracks/lcm/jfit110/mc_{key.replace('<=','le')}.json", "w"))
