# Copied unchanged from the campaign library (infra/lib/magtools.py) so that the exchange fit and Monte Carlo
# used for YBaMnFeO5 can be re-run from this repository. Only fit_J, mc_scan and their helpers are used here.
"""
magtools: collinear magnetic-configuration enumeration, Heisenberg fitting, and classical Monte Carlo.

Hamiltonian convention (used everywhere in this project):
    E = E0 - sum_{<ij>} J_ij  e_i . e_j          (each pair counted once, e_i unit vectors)
so J > 0 is ferromagnetic, and J already contains |S_i||S_j| (energy-mapping convention).
Classical MC with unit spins and these J gives T_c (classical). Quantum/classical rescaling
(e.g. (S+1)/S) is NOT applied; report both if needed.
"""
from __future__ import annotations

import itertools
import json
import math

import numpy as np

KB_MEV = 0.08617333262  # meV/K


# ----------------------------------------------------------------------------------------------
# enumeration
# ----------------------------------------------------------------------------------------------
def supercell_matrices(max_det=2):
    """A small set of HNF-like supercells up to det=max_det (diagonal + simple shears)."""
    mats = [np.eye(3, dtype=int)]
    if max_det >= 2:
        for i in range(3):
            m = np.eye(3, dtype=int); m[i, i] = 2; mats.append(m)
        for i, j in itertools.permutations(range(3), 2):
            m = np.eye(3, dtype=int); m[i, j] = 1; m[i, i] = 1; m[j, j] = 2
            mats.append(m)
            m2 = np.eye(3, dtype=int); m2[i, i] = 1; m2[j, i] = 1; m2[j, j] = 2
            mats.append(m2)
    if max_det >= 4:
        for i, j in itertools.combinations(range(3), 2):
            m = np.eye(3, dtype=int); m[i, i] = 2; m[j, j] = 2; mats.append(m)
            m = np.eye(3, dtype=int); m[i, i] = 1; m[i, j] = 1; m[j, i] = -1; m[j, j] = 1; mats.append(m)  # sqrt2 x sqrt2
    return mats


def enumerate_configs(structure, mag_sites_fn, max_det=2, max_per_cell=40, max_sites=16, include_fm=True, seed=0):
    """Enumerate symmetry-distinct collinear configurations over a few supercells.
    structure: pymatgen Structure (primitive). mag_sites_fn(site)->bool selects magnetic sites.
    Returns list of dicts {"matrix", "structure", "spins"} with spins in {+1,-1,0}.
    Dedup by StructureMatcher on species-decorated supercells (spin up/down -> distinct dummy labels)."""
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.core import Structure
    rng = np.random.default_rng(seed)
    out = []
    sm = StructureMatcher(primitive_cell=True, scale=False, attempt_supercell=True, ltol=0.1, stol=0.2, angle_tol=3)
    seen = []
    def deco(s, spins):
        sp = []
        for site, sg in zip(s, spins):
            el = site.specie.symbol
            sp.append({1: "Po", -1: "At", 0: el}[sg] if sg != 0 else el)
        return Structure(s.lattice, sp, s.frac_coords)
    for M in supercell_matrices(max_det):
        sc = structure.copy(); sc.make_supercell(M)
        idx = [i for i, site in enumerate(sc) if mag_sites_fn(site)]
        n = len(idx)
        if n == 0 or n > max_sites:
            continue
        combos = list(itertools.product([1, -1], repeat=n - 1))
        if len(combos) > 4096:
            combos = [tuple(rng.choice([1, -1], size=n - 1)) for _ in range(4096)]
        cnt = 0
        for c in combos:
            sp = np.zeros(len(sc), dtype=int)
            sp[idx] = (1,) + tuple(c)
            if not include_fm and np.all(sp[idx] == 1):
                continue
            d = deco(sc, sp)
            # also compare with the global spin flip
            flip = deco(sc, -sp)
            dup = False
            for q in seen:
                if len(q) == len(d) or True:
                    if sm.fit(q, d) or sm.fit(q, flip):
                        dup = True
                        break
            if dup:
                continue
            seen.append(d)
            out.append({"matrix": M.tolist(), "structure": sc, "spins": sp.tolist()})
            cnt += 1
            if cnt >= max_per_cell:
                break
    return out


# ----------------------------------------------------------------------------------------------
# Heisenberg fit
# ----------------------------------------------------------------------------------------------
def shells(structure, mag_sites_fn, rmax=6.0, tol=0.03):
    """Distinct (species pair, distance) shells among magnetic sites in a primitive structure."""
    idx = [i for i, s in enumerate(structure) if mag_sites_fn(s)]
    ds = []
    for i in idx:
        for nb in structure.get_neighbors(structure[i], rmax):
            if nb.index in idx:
                pair = tuple(sorted((structure[i].specie.symbol, structure[nb.index].specie.symbol)))
                ds.append((pair, nb.nn_distance))
    ds.sort(key=lambda x: x[1])
    sh = []
    for pair, d in ds:
        for s in sh:
            if s["pair"] == pair and abs(s["d"] - d) < tol:
                break
        else:
            sh.append({"pair": pair, "d": d})
    sh.sort(key=lambda s: s["d"])
    return sh


def pair_counts(sc_structure, spins, mag_sites_fn, shell_list, tol=0.03):
    """For a configuration, return sum over pairs (each once) of s_i s_j for each shell."""
    idx = [i for i, s in enumerate(sc_structure) if mag_sites_fn(s)]
    rmax = max(s["d"] for s in shell_list) + 0.1
    feats = np.zeros(len(shell_list))
    for i in idx:
        for nb in sc_structure.get_neighbors(sc_structure[i], rmax):
            j = nb.index
            if j not in idx:
                continue
            pair = tuple(sorted((sc_structure[i].specie.symbol, sc_structure[j].specie.symbol)))
            for k, s in enumerate(shell_list):
                if s["pair"] == pair and abs(s["d"] - nb.nn_distance) < tol:
                    feats[k] += 0.5 * spins[i] * spins[j]  # each pair seen twice
                    break
    return feats


def fit_J(features, energies, n_prim_cells):
    """Least squares E/ncell = E0 - sum_k J_k f_k/ncell. Returns J (meV), E0, rms (meV/cell), LOO rms."""
    F = np.asarray(features) / np.asarray(n_prim_cells)[:, None]
    E = np.asarray(energies) / np.asarray(n_prim_cells) * 1000.0  # meV per primitive cell
    X = np.hstack([np.ones((len(E), 1)), -F])
    coef, *_ = np.linalg.lstsq(X, E, rcond=None)
    pred = X @ coef
    rms = float(np.sqrt(np.mean((pred - E) ** 2)))
    loo = []
    for i in range(len(E)):
        m = np.ones(len(E), bool); m[i] = False
        if m.sum() <= X.shape[1]:
            continue
        c, *_ = np.linalg.lstsq(X[m], E[m], rcond=None)
        loo.append((X[i] @ c - E[i]) ** 2)
    return {"E0_meV": float(coef[0]), "J_meV": coef[1:].tolist(), "rms_meV_cell": rms,
            "loo_rms_meV_cell": float(np.sqrt(np.mean(loo))) if loo else None}


# ----------------------------------------------------------------------------------------------
# Monte Carlo (classical Heisenberg, unit spins), numba-accelerated
# ----------------------------------------------------------------------------------------------
def build_mc_lattice(structure, mag_sites_fn, shell_list, J_meV, L, tol=0.03, Jcut=1e-6):
    """Neighbor table for an L x L x L supercell of the primitive structure."""
    sc = structure.copy(); sc.make_supercell([L, L, L])
    idx = [i for i, s in enumerate(sc) if mag_sites_fn(s)]
    pos = {j: k for k, j in enumerate(idx)}
    rmax = max(s["d"] for s, J in zip(shell_list, J_meV) if abs(J) > Jcut) + 0.1
    nbrs = [[] for _ in idx]
    for i in idx:
        for nb in sc.get_neighbors(sc[i], rmax):
            j = nb.index
            if j not in pos:
                continue
            pair = tuple(sorted((sc[i].specie.symbol, sc[j].specie.symbol)))
            for k, s in enumerate(shell_list):
                if s["pair"] == pair and abs(s["d"] - nb.nn_distance) < tol and abs(J_meV[k]) > Jcut:
                    nbrs[pos[i]].append((pos[j], J_meV[k]))
                    break
    maxn = max(len(x) for x in nbrs)
    NB = -np.ones((len(idx), maxn), dtype=np.int64)
    JJ = np.zeros((len(idx), maxn))
    for a, lst in enumerate(nbrs):
        for b, (j, J) in enumerate(lst):
            NB[a, b] = j; JJ[a, b] = J
    return NB, JJ, sc, idx


try:
    import numba

    @numba.njit(cache=True)
    def _rand_unit():
        z = 2.0 * np.random.random() - 1.0
        phi = 2.0 * np.pi * np.random.random()
        r = math.sqrt(1.0 - z * z)
        return r * math.cos(phi), r * math.sin(phi), z

    @numba.njit(cache=True)
    def _mc_run(S, NB, JJ, Dz, T, nsweep, nequil, pattern, seed):
        np.random.seed(seed)
        N = S.shape[0]
        beta = 1.0 / (0.08617333262 * T) if T > 0 else 1e12
        Es = 0.0; E2s = 0.0; Ms = 0.0; M2s = 0.0; M4s = 0.0; cnt = 0
        for sweep in range(nsweep + nequil):
            for _ in range(N):
                i = np.random.randint(N)
                hx = 0.0; hy = 0.0; hz = 0.0
                for b in range(NB.shape[1]):
                    j = NB[i, b]
                    if j < 0:
                        break
                    hx += JJ[i, b] * S[j, 0]; hy += JJ[i, b] * S[j, 1]; hz += JJ[i, b] * S[j, 2]
                nx, ny, nz = _rand_unit()
                # E_i = -h.S_i - Dz*(S_iz)^2   (Dz>0 easy axis z)
                dE = -(hx * (nx - S[i, 0]) + hy * (ny - S[i, 1]) + hz * (nz - S[i, 2])) - Dz * (nz * nz - S[i, 2] * S[i, 2])
                if dE <= 0.0 or np.random.random() < math.exp(-beta * dE):
                    S[i, 0] = nx; S[i, 1] = ny; S[i, 2] = nz
            if sweep >= nequil:
                E = 0.0
                mx = 0.0; my = 0.0; mz = 0.0
                for i in range(N):
                    for b in range(NB.shape[1]):
                        j = NB[i, b]
                        if j < 0:
                            break
                        E -= 0.5 * JJ[i, b] * (S[i, 0] * S[j, 0] + S[i, 1] * S[j, 1] + S[i, 2] * S[j, 2])
                    E -= Dz * S[i, 2] * S[i, 2]
                    mx += pattern[i] * S[i, 0]; my += pattern[i] * S[i, 1]; mz += pattern[i] * S[i, 2]
                m = math.sqrt(mx * mx + my * my + mz * mz) / N
                Es += E; E2s += E * E; Ms += m; M2s += m * m; M4s += m ** 4; cnt += 1
        return Es / cnt, E2s / cnt, Ms / cnt, M2s / cnt, M4s / cnt

    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False


def mc_scan(NB, JJ, pattern, temps, nsweep=4000, nequil=2000, Dz=0.0, seed=1):
    """Temperature scan, returns list of dicts (T, E, C, M, chi, U4). pattern: +-1 ground-state signs."""
    N = NB.shape[0]
    S = np.zeros((N, 3)); S[:, 2] = np.asarray(pattern, float)
    out = []
    for k, T in enumerate(sorted(temps, reverse=True)[::-1]):
        pass
    S = np.zeros((N, 3)); S[:, 2] = np.asarray(pattern, float)
    for k, T in enumerate(sorted(temps)):  # heat from ordered state
        E, E2, M, M2, M4 = _mc_run(S, NB, JJ, Dz, T, nsweep, nequil, np.asarray(pattern, float), seed + k)
        C = (E2 - E * E) / (KB_MEV * T * T) / N
        chi = N * (M2 - M * M) / (KB_MEV * T)
        U4 = 1 - M4 / (3 * M2 * M2) if M2 > 0 else 0
        out.append({"T": float(T), "E": E / N, "C": C, "M": M, "chi": chi, "U4": U4})
    return out


def tc_from_scan(scan):
    T = np.array([s["T"] for s in scan]); chi = np.array([s["chi"] for s in scan]); C = np.array([s["C"] for s in scan])
    return {"Tc_chi": float(T[np.argmax(chi)]), "Tc_C": float(T[np.argmax(C)])}
