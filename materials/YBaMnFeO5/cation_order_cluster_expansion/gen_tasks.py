"""Generate B-site arrangements + species-uncorrelated random collinear spin configurations for the PM cluster expansion.
usage: python gen_tasks.py parent.json(with key 'relaxed' or plain) TAG A_ELEMENT out_tasks.json
Arrangement set (identical for both compounds; sigma=+1 on A = Mn or Cu, -1 on Fe):
  36-atom cells  A = 1x1x2, B = 2x1x1, C = [[1,1,0],[-1,1,0],[0,0,1]] (8 B sites): all symmetry/feature-distinct
                 50/50 arrangements are enumerated; a fixed, feature-diverse subset is used (see SEL below)
  72-atom cells  D = 2x2x1: SQS-like (NN correlations ~0) and rock-salt with one in-plane NN antisite pair;
                 F = [[1,1,0],[-1,1,0],[0,0,2]]: SQS-like.
Spin configs per arrangement: r0..r3 (half of each species up; chosen among all such configs to minimise the
species-resolved NN spin correlations, individually and on average), G (G-by-site: all NN AFM), FM.  r0 = relax config."""
import itertools, json, sys
import numpy as np
sys.path.insert(0, ".")
from ord2lib import supercell, pair_table, features, spin_features, CLASSES, ideal_site

CELLS = {"A": [[1, 0, 0], [0, 1, 0], [0, 0, 2]], "B": [[2, 0, 0], [0, 1, 0], [0, 0, 1]], "C": [[1, 1, 0], [-1, 1, 0], [0, 0, 1]],
         "D": [[2, 0, 0], [0, 2, 0], [0, 0, 1]], "F": [[1, 1, 0], [-1, 1, 0], [0, 0, 2]]}
SCLS = ("ip1", "ap", "Y")


def parity(sc):
    out = []
    for _, (X, Y, K) in sc["bsites"]:
        out.append((int(round(2 * Y)) + int(K)) % 2)
    return np.array(out)


def fkey(sig, bonds):
    S, N = features(sig, bonds)
    return tuple(round(S[k] / N[k], 3) if N[k] else 0.0 for k in CLASSES)


def spin_sets(sig, bonds, nset=4, seed=0):
    sig = np.array(sig)
    A = np.where(sig == 1)[0]; B = np.where(sig == -1)[0]
    nb = len(sig)
    cands = []
    combsA = list(itertools.combinations(A, len(A) // 2)); combsB = list(itertools.combinations(B, len(B) // 2))
    rng = np.random.default_rng(seed)
    allc = [(a, b) for a in combsA for b in combsB]
    if len(allc) > 6000:
        allc = [allc[i] for i in rng.choice(len(allc), 6000, replace=False)]
    seen = set()
    for a, b in allc:
        s = -np.ones(nb, int); s[list(a)] = 1; s[list(b)] = 1
        key = tuple(s) if s[0] == 1 else tuple(-s)
        if key in seen:
            continue
        seen.add(key)
        f = spin_features(s, sig, bonds, classes=SCLS)
        # normalise by the number of bonds of each species-resolved class
        n = spin_features(np.ones(nb), sig, bonds, classes=SCLS)
        v = np.array([f[k] / n[k] if n[k] else 0.0 for k in sorted(f)])
        cands.append((np.array(key), v))
    # individual score, keep the 60 best, then exhaustive/greedy set selection minimising |mean| and mean|.|
    cands.sort(key=lambda c: float(np.sum(c[1] ** 2)))
    pool = cands[:60]
    best, bsc = None, 1e9
    if len(pool) <= 22:
        it = itertools.combinations(range(len(pool)), nset)
    else:
        it = (tuple(rng.choice(len(pool), nset, replace=False)) for _ in range(20000))
    for comb in it:
        V = np.array([pool[i][1] for i in comb])
        sc = float(np.sum(V.mean(0) ** 2) + 0.25 * np.mean(np.sum(V ** 2, 1)))
        if sc < bsc:
            best, bsc = comb, sc
    return [pool[i][0].tolist() for i in best], bsc


def build(parent, tag, Ael, M, sig, name):
    sc = supercell(parent, M)
    bonds = pair_table(sc)
    species = list(sc["species"])
    for b, (ai, _) in enumerate(sc["bsites"]):
        species[ai] = Ael if sig[b] == 1 else "Fe"
    # reorder atoms: Ba, Y, B..., O (keep B-site order); remember mapping
    order = sorted(range(len(species)), key=lambda i: (["Ba", "Y", Ael, "Fe", "O"].index(species[i]) if species[i] != "Fe" else 2, i))
    # keep A and Fe interleaved in B-site order: sort B atoms by their bsite index
    bpos = {ai: b for b, (ai, _) in enumerate(sc["bsites"])}
    order = sorted(range(len(species)), key=lambda i: ({"Ba": 0, "Y": 1, "O": 3}.get(species[i], 2), bpos.get(i, 0), i))
    sp = [species[i] for i in order]; fr = [sc["frac"][i] for i in order]
    nb = len(sc["bsites"])
    bidx = [order.index(ai) for ai, _ in sc["bsites"]]  # position of each bsite in the reordered list
    sets, score = spin_sets(sig, bonds)
    par = parity(sc)
    G = [1 if p == 0 else -1 for p in par]
    cfgs = [(f"r{i}", s) for i, s in enumerate(sets)] + [("G", G), ("FM", [1] * nb)]
    configs = []
    for n, s in cfgs:
        full = [0] * len(sp)
        for b in range(nb):
            full[bidx[b]] = int(s[b])
        configs.append({"name": n, "spins": full, "bspins": [int(x) for x in s]})
    S, N = features(sig, bonds)
    return {"id": f"{tag}_{name}", "compound": tag, "cell": M, "natoms": len(sp), "nB": nb, "lattice": sc["lattice"], "species": sp, "frac": fr,
            "sigma": [int(x) for x in sig], "bidx": bidx, "parity": [int(p) for p in par],
            "S": S, "N": N, "Pi": {k: (S[k] / N[k] if N[k] else 0.0) for k in CLASSES}, "spin_score": score,
            "configs": configs, "relax_config": "r0"}


def sqs_search(parent, M, nsamp=200000, seed=1, weights=None, target=None):
    sc = supercell(parent, M); bonds = pair_table(sc); nb = len(sc["bsites"])
    rng = np.random.default_rng(seed)
    w = weights or {"ip1": 4, "ap": 2, "Y": 2, "ip2": 1, "apd": 1, "Yd": 1, "ip3": 0.5, "c2": 0.5}
    best, bs = None, 1e9
    for _ in range(nsamp):
        s = -np.ones(nb, int); s[rng.choice(nb, nb // 2, replace=False)] = 1
        S, N = features(s, bonds)
        sc_ = sum(w[k] * (S[k] / N[k] - (target or {}).get(k, 0.0)) ** 2 for k in w if N[k])
        if sc_ < bs:
            best, bs = s.copy(), sc_
    return best, bs


if __name__ == "__main__":
    pj = json.load(open(sys.argv[1]))
    parent = pj.get("relaxed", pj)
    tag, Ael, outp = sys.argv[2], sys.argv[3], sys.argv[4]
    tasks = []
    # ---- 36-atom: enumerate per cell
    cand = {}
    for cn in ("A", "B", "C"):
        sc = supercell(parent, CELLS[cn]); bonds = pair_table(sc); nb = len(sc["bsites"])
        par = parity(sc)
        for comb in itertools.combinations(range(nb), nb // 2):
            sig = -np.ones(nb, int); sig[list(comb)] = 1
            k = fkey(sig, bonds)
            cand.setdefault((cn, k), sig)
    # selection: (cell, Pi vector [ip1, ap, Y, ip2, apd, Yd, ip3, c2, ip4]) -> name
    SEL = {
        ("A", (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0)): "RS_A",
        ("C", (-1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0)): "RS_C",
        ("A", (-1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0)): "COL_A",
        ("A", (1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, 1.0)): "LAY_A",
        ("A", (-1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0)): "RSapY_A",     # ip checker, ap hetero, Y homo
        ("A", (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0)): "RSYap_A",     # ip checker, ap homo, Y hetero
        ("A", (1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0)): "LAYap_A",       # ip homo, ap hetero, Y homo
        ("A", (1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0)): "LAYY_A",        # ip homo, ap homo, Y hetero
        ("A", (-1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, -1.0)): "CHKmix_A",      # ip checker, stacking half RS/half COL
        ("A", (0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0)): "MIX_A",           # NN-uncorrelated (mixed layers)
        ("A", (0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0)): "MIXap_A",
        ("B", (0.0, -1.0, -1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0)): "STRd_B",       # in-plane stripes along the sqrt2 axis, c hetero
        ("B", (0.0, 1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 1.0, 0.0)): "STRdc_B",        # same stripes, c homo
        ("B", (0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 1.0, 0.0)): "MIX_B",
        ("C", (0.0, -1.0, -1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0)): "STRp_C",       # stripes along pseudo-cubic axis, c hetero
        ("C", (0.0, 1.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0)): "STRpc_C",        # same stripes, c homo
        ("C", (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0)): "MIX_C",
        ("C", (-0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, -0.5)): "HALF_C",
    }
    for (cn, k), name in SEL.items():
        assert (cn, k) in cand, (cn, k)
        tasks.append(build(parent, tag, Ael, CELLS[cn], cand[(cn, k)], name))
    # ---- 72-atom
    s, sc_ = sqs_search(parent, CELLS["D"], nsamp=60000, seed=11)
    tasks.append(build(parent, tag, Ael, CELLS["D"], s, "SQS_D"))
    s, sc_ = sqs_search(parent, CELLS["F"], nsamp=60000, seed=12)
    tasks.append(build(parent, tag, Ael, CELLS["F"], s, "SQS_F"))
    if Ael == "Mn":
        # rock-salt with one in-plane NN antisite pair (2x2x1)
        scD = supercell(parent, CELLS["D"]); bD = pair_table(scD); parD = parity(scD)
        sig = np.where(parD == 0, 1, -1)
        i0 = 0; j0 = [j for i, j, k in bD if i == i0 and k == "ip1"][0]
        sig[i0], sig[j0] = sig[j0], sig[i0]
        tasks.append(build(parent, tag, Ael, CELLS["D"], sig, "AS1_D"))
    else:
        # YBaCuFeO5 benchmark: every bipyramid Cu-O-Fe (ap hetero) with random bipyramid orientation (ip1, Y ~ 0)
        w = {"ip1": 4, "ap": 20, "Y": 2, "ip2": 1, "apd": 1, "Yd": 1, "ip3": 0.5, "c2": 0.5}
        s, sc_ = sqs_search(parent, CELLS["F"], nsamp=200000, seed=13, weights=w, target={"ap": -1.0, "apd": 0.0})
        tasks.append(build(parent, tag, Ael, CELLS["F"], s, "BPR_F"))
    json.dump({"tasks": tasks}, open(outp, "w"))
    for t in tasks:
        print(t["id"], t["natoms"], {k: round(v, 3) for k, v in t["Pi"].items()}, "spin_score=%.3f" % t["spin_score"])
