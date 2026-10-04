"""Exact (machine-precision) primitive cells for d3/d3 cyanide double perovskites A[M_N][M_C(CN)6] (Track L, PBA Gate 0).

Primitive fcc cell a/2 (0,1,1), (1,0,1), (1,1,0) (cubic axes = Cartesian axes), atom order
  [A (if any), M_N, M_C, C x6, N x6]       (same order as tracks/lcm/trilemma/nonoxide_pba/*_prim.cif)
Conventional frame: M_N at 0, M_C at (a/2,0,0), N at +-xN*a e_i, C at +-xC*a e_i (around M_N), A at (a/4,a/4,a/4) (4c).
Cartesian +e_x maps to primitive fractions (-1, 1, 1); so N(+e_x) = (-xN, xN, xN) etc. (mod 1).
symmetrize(): projects a (vc-)relaxed cell back onto the exact F-43m / Fm-3m form (a from the volume, xN / xC averaged
over all 18 coordinates each). QE HSE aborts with 'sym_rho_init_shell: lone vector' on ~1e-6 cell noise (TRACKL_LOG).
pure numpy + ase (no pymatgen / spglib needed in the Modal image)."""
import numpy as np

DIRS = [(-1, 1, 1), (1, -1, 1), (1, 1, -1)]  # primitive fractions of conventional +e_x, +e_y, +e_z (per unit a)
IMG = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]


def pba_cell(a, xN, xC, MN, MC, A=None):
    # cubic axes along Cartesian x, y, z (ase/QE fcc primitive): QE's ibrav=0 symmetry finder only tests rotations about
    # Cartesian axes, so the ase cellpar_to_cell orientation found 6 of 24 (F-43m) / 12 of 48 (Fm-3m) operations
    lat = 0.5 * a * np.array([[0.0, 1.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 0.0]])
    species, frac = [], []
    if A:
        species.append(A); frac.append([0.25, 0.25, 0.25])
    species += [MN, MC]; frac += [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]
    for el, x in (("C", xC), ("N", xN)):
        for d in DIRS:
            for sg in (1, -1):
                species.append(el)
                frac.append([(sg * x * c) % 1.0 for c in d])
    return {"lattice": lat.tolist(), "species": species, "frac": frac, "a": float(a), "xN": float(xN), "xC": float(xC)}


def spins_for(species, MN_index, MC_index, order="LCM"):
    """LCM: C-bound metal +1, N-bound metal -1 (Gate-0 spec); FM: both +1."""
    s = [0] * len(species)
    s[MC_index] = 1
    s[MN_index] = -1 if order == "LCM" else 1
    return s


def bonds(a, xN, xC):
    return {"M_N-N": xN * a, "C-N": (xC - xN) * a, "M_C-C": (0.5 - xC) * a}


def symmetrize(lattice, frac, species, MN, MC, A=None):
    """Exact projection of a relaxed PBA primitive cell (same atom order as pba_cell). Returns (cell_dict, info)."""
    lat = np.array(lattice, float); fr = np.array(frac, float)
    a = (4.0 * abs(np.linalg.det(lat))) ** (1.0 / 3.0)
    sp = list(species)
    iN = [i for i, s in enumerate(sp) if s == "N"]; iC = [i for i, s in enumerate(sp) if s == "C"]
    fN = fr[iN] % 1.0; fC = fr[iC] % 1.0
    xN = float(np.mean(np.minimum(fN, 1 - fN))); xC = float(np.mean(np.minimum(fC, 1 - fC)))
    cell = pba_cell(a, xN, xC, MN, MC, A)
    assert cell["species"] == sp, (cell["species"], sp)
    # deviation of the relaxed structure from the projection (Cartesian, minimum image in the relaxed cell)
    d = fr - np.array(cell["frac"]); d -= np.round(d)
    dev = float(np.max(np.linalg.norm(d @ lat, axis=1)))
    L = np.linalg.norm(lat, axis=1)
    ang = [float(np.degrees(np.arccos(np.dot(lat[i], lat[j]) / (L[i] * L[j])))) for i, j in ((1, 2), (0, 2), (0, 1))]
    info = {"a": a, "xN": xN, "xC": xC, "bonds": bonds(a, xN, xC), "max_atom_dev_A": dev,
            "spread_xN": float(np.ptp(np.minimum(fN, 1 - fN))), "spread_xC": float(np.ptp(np.minimum(fC, 1 - fC))),
            "relaxed_lengths": L.tolist(), "relaxed_angles": ang}
    return cell, info


def add_two_waters(cell, d_OO=2.80, d_OH=0.97):
    """KV[Cr(CN)6].2H2O model: two H2O in the EMPTY tetrahedral void (4d, 3/4 3/4 3/4; K sits in 4c), O-O along
    conventional e_x through the void centre (O 1.40 A from the centre, ~2.9-3.0 A from the nearest N/C of the cube
    faces). Each water donates its two H toward its two nearest framework N (O-H...N); P1 starting guess for a relax."""
    lat = np.array(cell["lattice"]); fr = np.array(cell["frac"]); sp = list(cell["species"])
    a = cell["a"]
    ex = np.array(DIRS[0], float) @ lat / a  # Cartesian unit vector of conventional +e_x
    c = np.array([0.75, 0.75, 0.75]) @ lat
    cart = fr @ lat
    O = [c + 0.5 * d_OO * ex, c - 0.5 * d_OO * ex]
    inv = np.linalg.inv(lat)
    newsp, newcart = [], []
    for o in O:
        # minimum-image vectors to framework N atoms
        vec = []
        for i, s in enumerate(sp):
            if s != "N":
                continue
            d = (cart[i] - o) @ inv; d -= np.round(d)
            vec += [(d + np.array(t)) @ lat for t in IMG]
        vec = sorted(vec, key=lambda v: np.linalg.norm(v))
        u1 = vec[0] / np.linalg.norm(vec[0])
        u2 = None
        for v in vec[1:]:
            u = v / np.linalg.norm(v)
            ang = np.degrees(np.arccos(np.clip(np.dot(u1, u), -1, 1)))
            if 75 <= ang <= 135:
                u2 = u; break
        if u2 is None:  # fall back to an ideal 104.5 deg pair in the plane of u1 and e_x
            w = ex - np.dot(ex, u1) * u1; w /= np.linalg.norm(w)
            th = np.radians(104.5); u2 = np.cos(th) * u1 + np.sin(th) * w
        newsp += ["O", "H", "H"]; newcart += [o, o + d_OH * u1, o + d_OH * u2]
    allc = np.vstack([cart, np.array(newcart)])
    out = dict(cell)
    out["species"] = sp + newsp
    out["frac"] = ((allc @ inv) % 1.0).tolist()
    return out


def min_dists(cell):
    """Shortest interatomic distance per element pair (27 images)."""
    lat = np.array(cell["lattice"]); fr = np.array(cell["frac"]); sp = cell["species"]
    out = {}
    for i in range(len(sp)):
        for j in range(i + 1, len(sp)):
            d = fr[j] - fr[i]; d -= np.round(d)
            r = float(min(np.linalg.norm((d + np.array(t)) @ lat) for t in IMG))
            k = "-".join(sorted((sp[i], sp[j])))
            out[k] = min(out.get(k, 99.0), r)
    return out
