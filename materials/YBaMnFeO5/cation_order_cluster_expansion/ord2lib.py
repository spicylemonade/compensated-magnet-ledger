"""B-site (Mn/Fe or Cu/Fe) order-disorder library for the 112 layered perovskite YBaBB'O5 (Track L, lcm_ord2).

Parent cell: sqrt2 x sqrt2 x 2 (18 atoms, 4 B sites). B sites in fractional coords of the parent: in-plane (0,1/2) or (1/2,0),
z ~ 0.26 (lower pyramid layer, apical O in the Ba plane at z=0) or z ~ 0.74. Pseudo-cubic in-plane coordinates:
u = x + y, v = y - x (units of a_p).
Pair classes on the ideal B lattice (direction-resolved, as in tracks/lcm/jlib.py):
  ip1  in-plane NN (a_p ~ 4.0 A)          ap   along c through the Ba layer (apical O, ~4.0 A)
  Y    along c across the O-free Y layer (~3.7 A)
  ip2  in-plane (1,1) a_p (5.66 A)        apd  Ba-side diagonal (in-plane 1 a_p + apical)     Yd   Y-side diagonal
  ip3  in-plane (2,0) a_p (8.0 A)         c2   same column, next bilayer (dz = c)            ip4 in-plane (2,1)
Ising variable sigma = +1 (Mn or Cu), -1 (Fe).  Pair energy convention: E = E0 + sum_k V_k sum_{bonds in k} sigma_i sigma_j.
"""
import itertools
import numpy as np

CLASSES = ["ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3", "c2", "ip4"]


def parent_B(parent):
    """indices and frac coords of the B sites (Mn/Fe/Cu) in the parent dict {lattice, species, frac}."""
    idx = [i for i, s in enumerate(parent["species"]) if s in ("Mn", "Fe", "Cu")]
    return idx


def ideal_site(f):
    """ideal (x, y, layer) for a parent-cell fractional coordinate of a B site; layer 0: z~0.26, 1: z~0.74"""
    x = round(2 * (f[0] % 1.0)) / 2 % 1.0
    y = round(2 * (f[1] % 1.0)) / 2 % 1.0
    z = f[2] % 1.0
    lay = 0 if z < 0.5 else 1
    return x, y, lay


def supercell(parent, M):
    """Supercell of the parent with integer matrix M (rows = new vectors in parent basis). Returns dict with lattice,
    species, frac (positions), and B-site info: list of (atom index, ideal parent-frame coords (X, Y, Zlayer-int))."""
    M = np.array(M, int)
    L = np.array(parent["lattice"], float)
    nL = M @ L
    det = int(round(abs(np.linalg.det(M))))
    # lattice points of parent inside the supercell
    rng = range(-4, 5)
    Minv = np.linalg.inv(M)
    pts = []
    for t in itertools.product(rng, rng, rng):
        f = np.array(t, float) @ Minv
        if np.all(f > -1e-8) and np.all(f < 1 - 1e-8):
            pts.append(np.array(t))
    assert len(pts) == det, (len(pts), det)
    species, frac, bsites = [], [], []
    pfr = np.array(parent["frac"], float)
    for i, s in enumerate(parent["species"]):
        for t in pts:
            pf = pfr[i] + t
            sf = pf @ Minv
            sf = sf - np.floor(sf + 1e-9)
            species.append(s)
            frac.append(sf.tolist())
            if s in ("Mn", "Fe", "Cu"):
                x, y, lay = ideal_site(pfr[i])
                # absolute ideal parent-frame coordinate: in-plane (x+tx, y+ty), vertical index 2*tz + lay
                bsites.append((len(species) - 1, (x + t[0], y + t[1], 2 * t[2] + lay)))
    # sort atoms by species order of parent (keep as is: grouped by parent atom); fine
    return {"lattice": nL.tolist(), "species": species, "frac": frac, "M": M.tolist(), "bsites": bsites}


def pair_table(sc, rmax_class=CLASSES):
    """Bond list on the B sublattice of a supercell: list of (i, j, class) with i, j indices into sc['bsites'] (0..NB-1),
    each bond counted once (including bonds to periodic images; i may equal j)."""
    M = np.array(sc["M"], int)
    Minv = np.linalg.inv(M)
    B = sc["bsites"]
    nb = len(B)
    # ideal coordinates as array (X, Y, K) with K vertical layer index (pyramid layer), parent periodicity in K is 2
    P = np.array([b[1] for b in B], float)
    bonds = []
    # search translations of the supercell
    T = [np.array(t) @ M for t in itertools.product(range(-2, 3), repeat=3)]
    for i in range(nb):
        for j in range(nb):
            for t in T:
                d = P[j] + np.array([t[0], t[1], 2 * t[2]]) - P[i]
                dx, dy, dk = d
                du, dv = dx + dy, dy - dx  # pseudo-cubic in-plane displacement (integers)
                du, dv = int(round(du)), int(round(dv))
                dk = int(round(dk))
                k = classify(du, dv, dk, int(P[i][2]) % 2)
                if k is None:
                    continue
                if (j, i, k) == (i, j, k) and np.allclose(d, 0):
                    continue
                bonds.append((i, j, k))
    # each bond counted twice (i->j and j->i); keep as directed list and weight 1/2
    return bonds


def classify(du, dv, dk, lay_i):
    """class for displacement (du, dv) in a_p, dk in pyramid layers, from a site in layer parity lay_i (0: z~0.26)."""
    a, b = sorted((abs(du), abs(dv)))
    if dk == 0:
        return {(0, 1): "ip1", (1, 1): "ip2", (0, 2): "ip3", (1, 2): "ip4"}.get((a, b))
    if abs(dk) == 1:
        # from layer 0 (z~0.26): dk=-1 goes through the Ba layer (apical), dk=+1 across Y; from layer 1 the reverse
        thruBa = (lay_i == 0 and dk == -1) or (lay_i == 1 and dk == 1)
        if (a, b) == (0, 0):
            return "ap" if thruBa else "Y"
        if (a, b) == (0, 1):
            return "apd" if thruBa else "Yd"
        return None
    if abs(dk) == 2 and (a, b) == (0, 0):
        return "c2"
    return None


def features(sigma, bonds, classes=CLASSES):
    """correlation sums S_k = sum_{bonds in k} s_i s_j (each bond once) and bond counts N_k."""
    S = {k: 0.0 for k in classes}
    N = {k: 0.0 for k in classes}
    for i, j, k in bonds:
        if k in S:
            S[k] += 0.5 * sigma[i] * sigma[j]
            N[k] += 0.5
    return S, N


def spin_features(spins, sigma, bonds, classes=("ip1", "ap", "Y", "ip2", "apd", "Yd")):
    """species-resolved spin correlation sums: key f'{k}_{pair}' with pair in AA (sigma=+1 both), AB, BB."""
    out = {}
    for k in classes:
        for p in ("AA", "AB", "BB"):
            out[f"{k}_{p}"] = 0.0
    for i, j, k in bonds:
        if k not in classes:
            continue
        p = "AA" if sigma[i] + sigma[j] == 2 else ("BB" if sigma[i] + sigma[j] == -2 else "AB")
        out[f"{k}_{p}"] += 0.5 * spins[i] * spins[j]
    return out
