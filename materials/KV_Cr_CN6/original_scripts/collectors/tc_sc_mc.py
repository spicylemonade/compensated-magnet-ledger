"""Classical Heisenberg Monte Carlo for the PBA magnetic lattice (Track L, PBA T_C model).

Lattice: all metals of A[M_N][M_C(CN)6] form a simple-cubic net with spacing d = a/2 (rock-salt colouring):
  N-bound metal at even parity (x+y+z even), C-bound metal at odd parity.
Hamiltonian (J > 0 antiferromagnetic, each pair once):  H = sum_<ij> J_ij  S_i . S_j , classical spins |S_i| = S_i.
Shells (vector in units of d): J1 (1,0,0) inter, 6 nb | J2C/J2N (1,1,0) intra, 12 nb | J3 (1,1,1) inter, 8 nb |
  J4C/J4N (2,0,0) intra, 6 nb (linear M-CN-M'-NC-M).
Disorder: [M_C(CN)6] vacancies = C-sublattice sites removed with probability 1 - p (N-bound metals stay, water-capped);
optional fraction w of N sites carrying a different spin (V(III), S = 1) with its own J scale.
Units: the kernel works with unit vectors and couplings K_ij = J_ij * (S_i/S0) * (S_j/S0); temperatures are returned in
units of J1 * S0^2 (S0 = 3/2), i.e. the classical (S^2) energy scale. For the pure J1 model T_C = 1.4430 J1 S0^2.
Updates per MC step: 1 heat-bath sweep + 1 over-relaxation sweep + Wolff embedding clusters (valid for any coupling
signs; reflections about a random plane) until ~N spins have been flipped.
Observables: staggered magnetisation m_s = (sum_N s - sum_C s)/N_occ (unit-vector weights), Binder ratio <m^4>/<m^2>^2,
and the two sublattice projections on the instantaneous staggered axis (finite-T compensation, M_N - M_C).
usage (library): from sc_mc import tc_scan; or  python sc_mc.py selftest"""
import math
import os
import sys

import numba
import numpy as np

SHELLS = [("J1", (1, 0, 0), "x"), ("J2", (1, 1, 0), "s"), ("J3", (1, 1, 1), "x"), ("J4", (2, 0, 0), "s")]


def shell_vectors(v):
    import itertools
    out = set()
    for p in itertools.permutations(v):
        for sg in itertools.product((1, -1), repeat=3):
            out.add(tuple(a * b for a, b in zip(p, sg)))
    return sorted(out)


def make_lattice(L, J, p=1.0, w=0.0, S_w=1.0, rJ_w=1.0, seed=0, S0=1.5):
    """J: dict with J1, J2C, J2N, J3, J4C, J4N (absolute or relative; only ratios to J1 matter for T/J1).
    Returns (nbr[N, nb], K[N, nb], sub[N] (+1 N-site, -1 C-site), occ[N] bool, slen[N])."""
    assert L % 2 == 0
    rng = np.random.default_rng(seed)
    N = L ** 3
    xs = np.arange(N) % L; ys = (np.arange(N) // L) % L; zs = np.arange(N) // (L * L)
    sub = np.where((xs + ys + zs) % 2 == 0, 1, -1).astype(np.int64)
    occ = np.ones(N, np.bool_)
    cidx = np.where(sub == -1)[0]
    nvac = int(round((1.0 - p) * len(cidx)))
    if nvac > 0:
        occ[rng.choice(cidx, nvac, replace=False)] = False
    slen = np.ones(N)
    jsc = np.ones(N)  # per-site exchange scale (V(III) on N sites)
    nidx = np.where(sub == 1)[0]
    nw = int(round(w * len(nidx)))
    if nw > 0:
        pick = rng.choice(nidx, nw, replace=False)
        slen[pick] = S_w / S0
        jsc[pick] = rJ_w
    slen[~occ] = 0.0
    vecs, keys = [], []
    for name, v, kind in SHELLS:
        for t in shell_vectors(v):
            vecs.append(t); keys.append((name, kind))
    nb = len(vecs)
    nbr = np.zeros((N, nb), np.int64)
    K = np.zeros((N, nb))
    for k, (t, (name, kind)) in enumerate(zip(vecs, keys)):
        j = ((xs + t[0]) % L) + L * (((ys + t[1]) % L) + L * ((zs + t[2]) % L))
        nbr[:, k] = j
        if kind == "x":
            Jv = np.full(N, J.get(name, 0.0))
        else:
            Jv = np.where(sub == 1, J.get(name + "N", 0.0), J.get(name + "C", 0.0))
        # bond scale: inter-sublattice bonds touching a V(III) get rJ_w (geometric mean for intra pairs)
        sc = np.sqrt(jsc * jsc[j]) if kind == "s" else jsc * jsc[j]
        K[:, k] = Jv * sc * slen * slen[j]
    # drop zero-coupling columns (unused shells) to speed up the kernel
    keep = np.where(np.any(K != 0.0, axis=0))[0]
    return nbr[:, keep].copy(), K[:, keep].copy(), sub, occ, slen


@numba.njit(cache=True)
def _rand_unit():
    z = 2.0 * np.random.random() - 1.0
    ph = 2.0 * math.pi * np.random.random()
    r = math.sqrt(max(0.0, 1.0 - z * z))
    return r * math.cos(ph), r * math.sin(ph), z


@numba.njit(cache=True)
def _field(i, s, nbr, K):
    hx = 0.0; hy = 0.0; hz = 0.0
    for k in range(nbr.shape[1]):
        kk = K[i, k]
        if kk != 0.0:
            j = nbr[i, k]
            hx -= kk * s[j, 0]; hy -= kk * s[j, 1]; hz -= kk * s[j, 2]
    return hx, hy, hz


@numba.njit(cache=True)
def _heatbath(s, nbr, K, act, beta):
    for ii in range(act.shape[0]):
        i = act[ii]
        hx, hy, hz = _field(i, s, nbr, K)
        H = math.sqrt(hx * hx + hy * hy + hz * hz)
        if H * beta < 1e-10:
            a, b, c = _rand_unit()
            s[i, 0] = a; s[i, 1] = b; s[i, 2] = c
            continue
        bh = beta * H
        u = np.random.random()
        ct = 1.0 + math.log(1.0 - u * (1.0 - math.exp(-2.0 * bh))) / bh
        if ct > 1.0:
            ct = 1.0
        if ct < -1.0:
            ct = -1.0
        st = math.sqrt(max(0.0, 1.0 - ct * ct))
        ph = 2.0 * math.pi * np.random.random()
        ex = hx / H; ey = hy / H; ez = hz / H
        # orthonormal frame around e
        if abs(ez) < 0.9:
            ux = -ey; uy = ex; uz = 0.0
        else:
            ux = 0.0; uy = -ez; uz = ey
        nu = math.sqrt(ux * ux + uy * uy + uz * uz)
        ux /= nu; uy /= nu; uz /= nu
        vx = ey * uz - ez * uy; vy = ez * ux - ex * uz; vz = ex * uy - ey * ux
        cp = math.cos(ph) * st; sp = math.sin(ph) * st
        s[i, 0] = ct * ex + cp * ux + sp * vx
        s[i, 1] = ct * ey + cp * uy + sp * vy
        s[i, 2] = ct * ez + cp * uz + sp * vz


@numba.njit(cache=True)
def _overrelax(s, nbr, K, act):
    for ii in range(act.shape[0]):
        i = act[ii]
        hx, hy, hz = _field(i, s, nbr, K)
        H2 = hx * hx + hy * hy + hz * hz
        if H2 < 1e-24:
            continue
        d = (s[i, 0] * hx + s[i, 1] * hy + s[i, 2] * hz) / H2
        s[i, 0] = 2.0 * d * hx - s[i, 0]
        s[i, 1] = 2.0 * d * hy - s[i, 1]
        s[i, 2] = 2.0 * d * hz - s[i, 2]


@numba.njit(cache=True)
def _wolff(s, nbr, K, act, beta, inclu, stack):
    """one Wolff embedding cluster; returns its size"""
    rx, ry, rz = _rand_unit()
    i0 = act[np.random.randint(act.shape[0])]
    top = 0
    stack[top] = i0; top += 1
    inclu[i0] = True
    size = 0
    while top > 0:
        top -= 1
        i = stack[top]
        pi = s[i, 0] * rx + s[i, 1] * ry + s[i, 2] * rz
        # flip i
        s[i, 0] -= 2.0 * pi * rx; s[i, 1] -= 2.0 * pi * ry; s[i, 2] -= 2.0 * pi * rz
        size += 1
        for k in range(nbr.shape[1]):
            kk = K[i, k]
            if kk == 0.0:
                continue
            j = nbr[i, k]
            if inclu[j]:
                continue
            pj = s[j, 0] * rx + s[j, 1] * ry + s[j, 2] * rz
            # embedded Ising bond energy before the flip of i: kk * pi * pj ; satisfied if < 0
            e = kk * pi * pj
            if e < 0.0:
                if np.random.random() < 1.0 - math.exp(2.0 * beta * e):
                    inclu[j] = True
                    stack[top] = j; top += 1
    # reset flags
    return size


@numba.njit(cache=True)
def _clear(inclu, s, act):
    for ii in range(act.shape[0]):
        inclu[act[ii]] = False


@numba.njit(cache=True)
def _measure(s, sub, act):
    mx = 0.0; my = 0.0; mz = 0.0
    ax = 0.0; ay = 0.0; az = 0.0; bx = 0.0; by = 0.0; bz = 0.0
    nA = 0; nB = 0
    for ii in range(act.shape[0]):
        i = act[ii]
        if sub[i] == 1:
            ax += s[i, 0]; ay += s[i, 1]; az += s[i, 2]; nA += 1
        else:
            bx += s[i, 0]; by += s[i, 1]; bz += s[i, 2]; nB += 1
    n = nA + nB
    mx = (ax - bx) / n; my = (ay - by) / n; mz = (az - bz) / n
    m = math.sqrt(mx * mx + my * my + mz * mz)
    if m > 0:
        ex = mx / m; ey = my / m; ez = mz / m
    else:
        ex = 0.0; ey = 0.0; ez = 1.0
    MA = (ax * ex + ay * ey + az * ez) / max(nA, 1)
    MB = -(bx * ex + by * ey + bz * ez) / max(nB, 1)
    return m, MA, MB


@numba.njit(cache=True)
def _energy(s, nbr, K, act):
    E = 0.0
    for ii in range(act.shape[0]):
        i = act[ii]
        for k in range(nbr.shape[1]):
            kk = K[i, k]
            if kk != 0.0:
                j = nbr[i, k]
                E += 0.5 * kk * (s[i, 0] * s[j, 0] + s[i, 1] * s[j, 1] + s[i, 2] * s[j, 2])
    return E


@numba.njit(cache=True)
def _run(s, nbr, K, sub, act, beta, nequil, nmeas, seed):
    np.random.seed(seed)
    N = s.shape[0]
    inclu = np.zeros(N, np.bool_)
    stack = np.zeros(N, np.int64)
    nact = act.shape[0]
    out = np.zeros((nmeas, 4))
    for it in range(nequil + nmeas):
        _heatbath(s, nbr, K, act, beta)
        _overrelax(s, nbr, K, act)
        flipped = 0
        ncl = 0
        while flipped < nact and ncl < 4 * nact:
            c = _wolff(s, nbr, K, act, beta, inclu, stack)
            _clear(inclu, s, act)
            flipped += c
            ncl += 1
        if it >= nequil:
            m, MA, MB = _measure(s, sub, act)
            out[it - nequil, 0] = m
            out[it - nequil, 1] = MA
            out[it - nequil, 2] = MB
            out[it - nequil, 3] = _energy(s, nbr, K, act) / nact
    return out


def simulate(L, T, J, p=1.0, w=0.0, S_w=1.0, rJ_w=1.0, dis_seed=0, mc_seed=1, nequil=400, nmeas=2000):
    """T in units of J1*S0^2 (J['J1'] should be 1 for that convention). Returns dict of moments."""
    nbr, K, sub, occ, slen = make_lattice(L, J, p, w, S_w, rJ_w, dis_seed)
    # unit spins on occupied sites; vacancies carry no coupling (K = 0 rows/cols) and are excluded from act
    rng = np.random.default_rng(mc_seed)
    s = rng.normal(size=(L ** 3, 3)); s /= np.linalg.norm(s, axis=1)[:, None]
    act = np.where(occ)[0].astype(np.int64)
    # order the active sites so N and C alternate (sequential sweep)
    o = _run(s, nbr, K, sub, act, 1.0 / T, nequil, nmeas, mc_seed)
    m = o[:, 0]
    return {"L": L, "T": T, "m1": float(m.mean()), "m2": float((m ** 2).mean()), "m4": float((m ** 4).mean()),
            "MA": float(o[:, 1].mean()), "MB": float(o[:, 2].mean()), "dM": float((o[:, 1] - o[:, 2]).mean()),
            "dM_err": float((o[:, 1] - o[:, 2]).std() / math.sqrt(len(m) / 20.0)),
            "E": float(o[:, 3].mean()), "nact": int(occ.sum())}


def selftest():
    import time
    J = {"J1": 1.0}
    for L in (8, 12, 16):
        for T in (1.40, 1.443, 1.48):
            t0 = time.time()
            r = simulate(L, T, J, nequil=200, nmeas=1000, mc_seed=L * 100 + int(T * 1000))
            print(L, T, round(r["m2"] ** 0.5, 4), round(r["m4"] / r["m2"] ** 2, 4), round(r["E"], 4), f"{time.time() - t0:.1f}s",
                  flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
