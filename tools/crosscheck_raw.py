"""Cross-check recorded energies against the raw pw.x outputs of the cation-ordering and hull runs.

The cluster expansion (claims Y25, Y26) and the hull distance (Y24) were computed by the agents from parsed
energies stored in JSON files. This script re-reads the raw outputs included in
materials/YBaMnFeO5/runs/J_cation_order_cluster_expansion_raw/ and K_hull_competing_phases_raw/, checks that
every recorded energy matches the output it came from, and recomputes the hull distance from the raw energies
(lowest-energy mix of competing phases at fixed Y:Ba:Mn:Fe:O, as in the agents' hull_post.py). Only numpy is
needed.

    python tools/crosscheck_raw.py
"""
import glob
import json
import pathlib
import sys
from collections import Counter
from functools import lru_cache
from itertools import combinations

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from qe_parse import parse_pw  # noqa: E402

Y = ROOT / "materials" / "YBaMnFeO5"
J = Y / "runs" / "J_cation_order_cluster_expansion_raw"
K = Y / "runs" / "K_hull_competing_phases_raw"
ELEMENTS = ["Y", "Ba", "Mn", "Fe", "O"]
TARGET = np.array([1, 1, 1, 1, 5]) / 9  # YBaMnFeO5, atomic fractions in ELEMENTS order
TOL_EV = 1e-4  # recorded values were parsed from the same files; anything above this is a mismatch
# The agents wrote analysis_data/hull_QE110.json at 2026-10-02 10:43:49 PDT (timestamp of the file in the campaign
# archive). Competing phases whose job finished after that were not in their hull.
HULL_FILE_WRITTEN = 1790963029  # 2026-10-02 10:43:49 PDT as a Unix time


def _find_out(jobs, fname):
    for job in jobs:
        d = J / job.split("/ord2/", 1)[-1]
        for cand in (d / (fname + ".gz"), d / fname):
            if cand.exists():
                return cand
    return None


def check_ordering():
    """Every energy in cation_order_cluster_expansion/results_*/ vs the raw output it was parsed from."""
    out = {"n": 0, "ok": 0, "worst_eV": 0.0, "truncated": [], "missing": [], "bad": []}
    for res in sorted(glob.glob(str(Y / "cation_order_cluster_expansion" / "results_*" / "*.json"))):
        d = json.load(open(res))
        recs = [(name, r, f"{d['id']}_{name}.out") for name, r in d.get("runs", {}).items()]
        if d.get("relax") and d["relax"].get("energy_eV") is not None:
            recs.append(("relax", d["relax"], f"{d['id']}_relax.out"))
        for name, r, fname in recs:
            if r.get("energy_eV") is None:
                continue
            out["n"] += 1
            f = _find_out(d.get("jobs", []), fname)
            if f is None:
                out["missing"].append(f"{d['id']} {name}")
                continue
            e = parse_pw(f)["energy_eV"]
            if e is None:  # the saved output stops before the final energy line
                out["truncated"].append(f"{d['id']} {name} ({f.relative_to(ROOT)})")
                continue
            dE = abs(e - r["energy_eV"])
            if dE <= TOL_EV:
                out["ok"] += 1
                out["worst_eV"] = max(out["worst_eV"], dE)
            else:
                out["bad"].append(f"{d['id']} {name}: recorded {r['energy_eV']:.6f} eV, raw {e:.6f} eV")
    return out


@lru_cache(maxsize=None)
def hull_phases():
    """Final vc-relax energy per atom of every completed hull run, re-read from its raw output."""
    rows, bad = {}, []
    for jd in sorted(p for p in K.iterdir() if p.is_dir()):
        rec = jd / f"h110_{jd.name}.json"
        if not rec.exists():
            continue  # job failed before a final energy (see its status.json)
        d = json.load(open(rec))
        fin = d.get("final") or {}
        if not fin.get("E_eV"):
            continue
        f = jd / "vc110.out.gz" if (jd / "vc110.out.gz").exists() else jd / "vc110.out"
        e = parse_pw(f)["energy_eV"]
        if e is None or abs(e - fin["E_eV"]) > TOL_EV:
            bad.append(f"{jd.name}: recorded {fin['E_eV']:.6f} eV, raw {e}")
            continue
        cnt = Counter(d["relaxed"]["species"])
        n = sum(cnt.values())
        status = json.load(open(jd / "status.json")) if (jd / "status.json").exists() else {}
        rows[jd.name] = {"formula": d.get("formula", jd.name), "frac": np.array([cnt.get(el, 0) / n for el in ELEMENTS]),
                         "E_atom": e / n, "finished": status.get("end")}
    return rows, tuple(bad)


def lowest_mix(phases, target=TARGET):
    """Minimise sum(x_i E_i) subject to sum(x_i frac_i) = target, x >= 0 (a small linear program).

    The optimum sits at a vertex that uses at most len(ELEMENTS) phases, so every such subset is tried
    (exact; a few hundred thousand 5x5 solves, about a second with numpy)."""
    names = list(phases)
    A = np.array([phases[k]["frac"] for k in names]).T  # elements x phases
    c = np.array([phases[k]["E_atom"] for k in names])
    m = len(ELEMENTS)
    idx = np.array(list(combinations(range(len(names)), m)))
    M = A[:, idx].transpose(1, 0, 2)  # subsets x elements x m
    keep = np.abs(np.linalg.det(M)) > 1e-12
    idx, M = idx[keep], M[keep]
    X = np.linalg.solve(M, np.broadcast_to(target, (len(M), m))[..., None])[..., 0]
    feasible = (X > -1e-10).all(axis=1)
    cost = (c[idx] * X).sum(axis=1)
    cost[~feasible] = np.inf
    i = int(np.argmin(cost))
    return float(cost[i]), {phases[names[j]]["formula"]: round(float(x), 4) for j, x in zip(idx[i], X[i]) if x > 1e-8}


def hull_distance(target_id="YBMFO_P4n", before=None):
    """Energy above the hull (meV/atom) of a YBaMnFeO5 run, from raw energies; `before` limits the competing
    phases to jobs that finished before that Unix time."""
    rows, _ = hull_phases()
    comps = {k: v for k, v in rows.items() if not np.allclose(v["frac"], TARGET)
             and (before is None or (v["finished"] or 0) < before)}
    e_mix, dec = lowest_mix(comps)
    return 1000 * (rows[target_id]["E_atom"] - e_mix), dec, len(comps)


if __name__ == "__main__":
    o = check_ordering()
    print(f"cation-ordering runs: {o['ok']}/{o['n']} recorded energies match their raw outputs "
          f"(largest difference {o['worst_eV']:.1e} eV); {len(o['truncated'])} saved outputs end before the final energy, "
          f"{len(o['missing'])} not found, {len(o['bad'])} mismatches")
    for t in o["truncated"]:
        print("  truncated:", t)
    for t in o["missing"] + o["bad"]:
        print("  PROBLEM:", t)
    rows, hbad = hull_phases()
    ncomp = sum(not np.allclose(v["frac"], TARGET) for v in rows.values())
    print(f"hull runs: {len(rows)} recorded final energies match their raw vc-relax outputs "
          f"({ncomp} competing phases, {len(rows) - ncomp} YBaMnFeO5 cells); {len(hbad)} mismatches")
    for t in hbad:
        print("  PROBLEM:", t)
    rec = json.load(open(Y / "analysis_data" / "hull_QE110.json"))["YBMFO_P4n"]["E_hull_meV_atom"]
    eh0, dec0, n0 = hull_distance(before=HULL_FILE_WRITTEN)
    eh, dec, n = hull_distance()
    print(f"YBaMnFeO5 above the hull, competing phases finished before the agents' hull file ({n0}): {eh0:+.1f} meV/atom "
          f"(agents recorded {rec}); decomposes to {dec0}")
    print(f"YBaMnFeO5 above the hull, all {n} completed competing phases: {eh:+.1f} meV/atom; decomposes to {dec}")
    sys.exit(1 if (o["bad"] or o["missing"] or hbad) else 0)
