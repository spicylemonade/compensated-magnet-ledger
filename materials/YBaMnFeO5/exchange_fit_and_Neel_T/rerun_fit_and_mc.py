"""Level-2 check for YBaMnFeO5: raw DFT energies -> Heisenberg exchange fit -> classical Monte Carlo Neel temperature.

Re-derives the exchange constants (J_ip, J_ap, J_Y, ...) and the Binder-crossing T_N (recorded: 417 +- 10 K)
directly from the raw pw.x outputs stored in this repository:
    ../runs/I_pbeu_exchange_fit_configs/j110_*/single.out    24 spin configurations (PBE+U, U = 4 eV, 90 Ry)
    ../runs/H_pbeu_stacking_vs_U/U44_*/single.out            G / Yflip / Baflip stackings (same settings)
This is the same procedure as the campaign scripts fit110.py and mc110.py (kept next to this file), with the
energies read from the raw outputs instead of an intermediate results file.

usage:  python rerun_fit_and_mc.py            (full: L = 6, 8, 10; 12000 sweeps; a few minutes with numba)
        python rerun_fit_and_mc.py --quick    (fewer sweeps; T_N to within ~10-20 K)
needs:  numpy, pymatgen, numba
"""
import json
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(HERE))
from qe_parse import parse_pw  # noqa: E402
import magtools as mt  # noqa: E402
import jlib  # noqa: E402
from pymatgen.core import Lattice, Structure  # noqa: E402

QUICK = "--quick" in sys.argv
RUNS = HERE.parent / "runs"
idx = json.load(open(HERE / "jfit110" / "index.json"))
P = Structure(Lattice(idx["prim"]["lattice"]), idx["prim"]["species"], idx["prim"]["frac"])


def n_flipped(spins, moments):
    """Magnetic sites whose converged moment sign differs from the intended spin."""
    s = np.array(spins, float); m = np.array(moments[: len(s)], float)
    mag = (s != 0) & (np.abs(m) > 1.0)
    return int(np.sum(np.sign(m[mag]) != np.sign(s[mag])))


C = []
for o in idx["configs"]:
    r = parse_pw(RUNS / "I_pbeu_exchange_fit_configs" / o["id"] / "single.out")
    if r["energy_eV"] is None or not r["converged"] or n_flipped(o["spins"], r["site_moments"]):
        print("skip", o["id"]); continue
    n = round(o["natoms"] / len(P))
    C.append({"matrix": o["matrix"], "spins": o["spins"], "E": r["energy_eV"], "n": n, "M": abs(r["total_mag"] or 0) / n, "tag": o["id"]})
st_in = json.load(open(RUNS / "H_pbeu_stacking_vs_U" / "U44_G" / "ybmfo110_stack.json"))
for cf in st_in["configs"]:
    r = parse_pw(RUNS / "H_pbeu_stacking_vs_U" / f"U44_{cf['name']}" / "single.out")
    if r["energy_eV"] is None or not r["converged"] or n_flipped(cf["spins"], r["site_moments"]):
        print("skip stack", cf["name"]); continue
    C.append({"matrix": [[1, 0, 0], [0, 1, 0], [0, 0, 2]], "spins": cf["spins"], "E": r["energy_eV"], "n": 2,
              "M": abs(r["total_mag"] or 0) / 2, "tag": "st_" + cf["name"]})
seen, U = set(), []
for c in C:
    k = (json.dumps(c["matrix"]), tuple(c["spins"]))
    if k not in seen:
        seen.add(k); U.append(c)
print(len(U), "unique configurations from raw outputs")

found = {}
for i, s in enumerate(P):
    if s.specie.symbol not in ("Mn", "Fe"):
        continue
    for nb in P.get_neighbors(s, 7.0):
        if P[nb.index].specie.symbol not in ("Mn", "Fe"):
            continue
        k = jlib.cls(s.specie.symbol, P[nb.index].specie.symbol, nb.nn_distance, abs(nb.coords[2] - s.coords[2]))
        if k:
            found[k] = found.get(k, 0) + 1
CL = ["ip", "ap", "Y"] + sorted([k for k in found if k not in ("ip", "ap", "Y")], key=lambda k: (float(k.split("_")[-1]), k))
F = []
for c in U:
    sc = P.copy(); sc.make_supercell(c["matrix"]); F.append(jlib.features(sc, c["spins"], CL))
F = np.array(F); E = np.array([c["E"] for c in U]); N = np.array([c["n"] for c in U]); M = np.array([c["M"] for c in U])
mask = M <= 10.5  # 'clean' set: exclude ferromagnet-like, partly metallic states (as in the original fit)
nfit = 10
fit = mt.fit_J(F[mask][:, :nfit], E[mask], N[mask])
J = dict(zip(CL[:nfit], fit["J_meV"]))
rec = json.load(open(HERE / "jfit110" / "fits.json"))["clean_M<=10.5_10"]
print("re-fitted J (meV, unit spins, E = E0 - sum J s_i.s_j):")
for k in CL[:nfit]:
    print(f"   {k:14s} {J[k]:8.2f}    recorded {dict(zip(rec['classes'], rec['J_meV']))[k]:8.2f}")
print(f"   rms {fit['rms_meV_cell']:.2f} meV/cell (recorded {rec['rms']:.2f}); LOO {fit['loo_rms_meV_cell']:.2f} (recorded {rec['loo']:.2f}); {int(mask.sum())} configs")


def lattice(L):
    sc = P.copy(); sc.make_supercell([L, L, max(2, L // 2)])
    ii = [i for i, s in enumerate(sc) if s.specie.symbol in ("Mn", "Fe")]
    pos = {j: k for k, j in enumerate(ii)}
    nb = [[] for _ in ii]
    for i in ii:
        for n in sc.get_neighbors(sc[i], 7.0):
            if n.index not in pos:
                continue
            k = jlib.cls(sc[i].specie.symbol, sc[n.index].specie.symbol, n.nn_distance, abs(n.coords[2] - sc[i].coords[2]))
            if k in J and abs(J[k]) > 1e-9:
                nb[pos[i]].append((pos[n.index], J[k]))
    m = max(len(x) for x in nb)
    NB = -np.ones((len(ii), m), dtype=np.int64); JJ = np.zeros((len(ii), m))
    for a, l in enumerate(nb):
        for b, (j, v) in enumerate(l):
            NB[a, b] = j; JJ[a, b] = v
    return NB, JJ, np.array([1.0 if sc[i].specie.symbol == "Mn" else -1.0 for i in ii])


NB, JJ, pat = lattice(6)
Eg = -0.5 * sum(JJ[a, b] * pat[a] * pat[NB[a, b]] for a in range(len(pat)) for b in range(NB.shape[1]) if NB[a, b] >= 0) / len(pat)
Tmf = -Eg * 2 / 3 / 0.08617333
temps = np.arange(round(0.5 * Tmf / 10) * 10, round(0.72 * Tmf / 10) * 10 + 1, 10.0)
nsw, neq = (3000, 1500) if QUICK else (12000, 5000)
print(f"mean-field T = {Tmf:.0f} K; Monte Carlo {temps[0]:.0f}-{temps[-1]:.0f} K, {nsw} sweeps")
out = {}
for L in (6, 8, 10):
    NB, JJ, pat = lattice(L)
    sc = mt.mc_scan(NB, JJ, pat, temps, nsweep=nsw, nequil=neq, seed=300 + L)
    out[L] = [(s["T"], s["U4"]) for s in sc]
crossings = []
for a, b in ((6, 8), (8, 10), (6, 10)):
    ta = np.array([x[0] for x in out[a]]); d = np.array([x[1] for x in out[a]]) - np.array([x[1] for x in out[b]])
    ii = np.where(np.diff(np.sign(d)) != 0)[0]
    c = [ta[i] + (ta[i + 1] - ta[i]) * d[i] / (d[i] - d[i + 1]) for i in ii]
    if c:
        crossings.append(c[0])
    print(f"Binder crossing L={a}/{b}: {[round(x) for x in c]}")
if crossings:
    print(f"T_N (classical MC, raw) = {np.mean(crossings):.0f} +- {np.std(crossings):.0f} K   (recorded 417 +- 10 K)")
json.dump({"J_meV": J, "rms": fit["rms_meV_cell"], "loo": fit["loo_rms_meV_cell"], "T_MF": Tmf,
           "binder_crossings_K": crossings, "quick": QUICK},
          open(HERE / ("rerun_result_quick.json" if QUICK else "rerun_result.json"), "w"), indent=1)
