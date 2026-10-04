"""Canonical (fixed 50/50 composition) Monte Carlo of the B-site pair cluster expansion on the 112 B sublattice.
B sublattice = simple tetragonal grid (u, v, K): u, v pseudo-cubic in-plane indices, K pyramid-layer index; vertical bond
(K, K+1) crosses the Y layer for even K and the Ba layer (apical O) for odd K (see ord2lib.classify).
E = sum_k V_k sum_{bonds in k} sigma_i sigma_j  (V in meV per bond, sigma = +1 A (Mn/Cu), -1 Fe).
Moves: Kawasaki exchanges of an unlike pair (50 % nearest-neighbour, 50 % arbitrary).
Order parameters: Fourier amplitudes eta(q) = |sum_i sigma_i exp(i q.r_i)| / N for q in {0, pi/2, pi}^3 (u, v, K);
reported: eta_RS = eta(pi,pi,pi), eta_COL = eta(pi,pi,0), eta_LAY = eta(0,0,pi), eta_max = max over q != 0, plus pair SRO
<sigma sigma> per class.
usage: python mc_ce.py '{"ip1":40,"ap":20,"Y":10}' [L] [Lz] [Tmin] [Tmax] [nT] [nsweep] [seed]"""
import json, sys, math
import numpy as np
from numba import njit

KB = 0.08617333262  # meV/K
CLS = ["ip1", "ap", "Y", "ip2", "apd", "Yd", "ip3", "c2", "ip4"]


def offsets():
    """list of (du, dv, dK, class_even, class_odd): class for a site in an even-K / odd-K layer"""
    out = []
    for du, dv in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        out.append((du, dv, 0, "ip1", "ip1"))
    for du, dv in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        out.append((du, dv, 0, "ip2", "ip2"))
    for du, dv in [(2, 0), (-2, 0), (0, 2), (0, -2)]:
        out.append((du, dv, 0, "ip3", "ip3"))
    for du, dv in [(2, 1), (2, -1), (-2, 1), (-2, -1), (1, 2), (1, -2), (-1, 2), (-1, -2)]:
        out.append((du, dv, 0, "ip4", "ip4"))
    # vertical: from even K, dK=+1 across Y, dK=-1 through Ba; from odd K the reverse
    out.append((0, 0, 1, "Y", "ap")); out.append((0, 0, -1, "ap", "Y"))
    for du, dv in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        out.append((du, dv, 1, "Yd", "apd")); out.append((du, dv, -1, "apd", "Yd"))
    out.append((0, 0, 2, "c2", "c2")); out.append((0, 0, -2, "c2", "c2"))
    return out


def build(L, Lz, V):
    N = L * L * Lz
    offs = offsets()
    nb, vv, cl = [], [], []
    for K in range(Lz):
        for v in range(L):
            for u in range(L):
                lst, lv, lc = [], [], []
                for du, dv, dK, ce, co in offs:
                    c = ce if K % 2 == 0 else co
                    j = (((K + dK) % Lz) * L + (v + dv) % L) * L + (u + du) % L
                    lst.append(j); lv.append(V.get(c, 0.0)); lc.append(CLS.index(c))
                nb.append(lst); vv.append(lv); cl.append(lc)
    return np.array(nb, np.int64), np.array(vv, float), np.array(cl, np.int64)


@njit(cache=True)
def local_field(s, NB, VV, i):
    h = 0.0
    for b in range(NB.shape[1]):
        h += VV[i, b] * s[NB[i, b]]
    return h


@njit(cache=True)
def total_energy(s, NB, VV):
    e = 0.0
    for i in range(s.shape[0]):
        e += s[i] * local_field(s, NB, VV, i)
    return 0.5 * e


@njit(cache=True)
def sweeps(s, NB, VV, beta, nsw, seed):
    np.random.seed(seed)
    N = s.shape[0]
    acc = 0
    for _ in range(nsw * N):
        i = np.random.randint(N)
        if np.random.random() < 0.5:
            j = NB[i, np.random.randint(4)]
        else:
            j = np.random.randint(N)
        if s[i] == s[j]:
            continue
        hi = local_field(s, NB, VV, i); hj = local_field(s, NB, VV, j)
        vij = 0.0
        for b in range(NB.shape[1]):
            if NB[i, b] == j:
                vij += VV[i, b]
        dE = -2.0 * s[i] * hi - 2.0 * s[j] * hj - 4.0 * vij
        if dE <= 0.0 or np.random.random() < math.exp(-beta * dE):
            s[i] = -s[i]; s[j] = -s[j]; acc += 1
    return acc


def order_params(s, L, Lz, phases):
    return {k: float(abs(np.sum(s * ph)) / len(s)) for k, ph in phases.items()}


def make_phases(L, Lz):
    K, v, u = np.meshgrid(np.arange(Lz), np.arange(L), np.arange(L), indexing="ij")
    u = u.ravel(); v = v.ravel(); K = K.ravel()
    ph = {}
    qs = [0, 0.5, 1.0]
    for qu in qs:
        for qv in qs:
            for qk in qs:
                if qu == qv == qk == 0:
                    continue
                if (qu == 0.5 and L % 4) or (qv == 0.5 and L % 4) or (qk == 0.5 and Lz % 4):
                    continue
                ph[(qu, qv, qk)] = np.exp(1j * np.pi * (qu * u + qv * v + qk * K))
    return ph


def sro(s, NB, CL):
    out = np.zeros(len(CLS)); cnt = np.zeros(len(CLS))
    pr = s[:, None] * s[NB]
    for c in range(len(CLS)):
        m = CL == c
        if m.any():
            out[c] = pr[m].mean(); cnt[c] = m.sum()
    return {CLS[c]: float(out[c]) for c in range(len(CLS)) if cnt[c]}


def run(V, L=12, Lz=12, temps=None, nsweep=2000, nequil=1000, seed=0, start=None):
    NB, VV, CL = build(L, Lz, V)
    N = NB.shape[0]
    rng = np.random.default_rng(seed)
    s = np.ones(N, np.int64); s[rng.choice(N, N // 2, replace=False)] = -1
    if start is not None:
        s = start.copy()
    phases = make_phases(L, Lz)
    keys = list(phases)
    P = np.array([phases[k] for k in keys])
    out = []
    for it, T in enumerate(temps):
        beta = 1.0 / (KB * T)
        sweeps(s, NB, VV, beta, nequil, seed + 7 * it)
        Es, eta, eta2, eta4 = [], [], [], []
        sr = {k: 0.0 for k in CLS}; nsr = 0
        nblk = 50
        per = max(1, nsweep // nblk)
        etaq = np.zeros(len(keys))
        for b in range(nblk):
            sweeps(s, NB, VV, beta, per, seed + 1000 * it + b + 1)
            Es.append(total_energy(s, NB, VV))
            a = np.abs(P @ s) / N
            etaq += a
            eta.append(a)
            if b % 5 == 0:
                for k, x in sro(s, NB, CL).items():
                    sr[k] += x
                nsr += 1
        Es = np.array(Es); eta = np.array(eta)
        e_fu = 2 * Es / N
        C = np.var(Es) / (KB * T) ** 2 / N * 2  # per f.u., in k_B units
        etaq /= nblk
        imax = int(np.argmax(etaq))
        rec = {"T": float(T), "E_meV_fu": float(e_fu.mean()), "C_kB_fu": float(C), "q_max": keys[imax], "eta_max": float(etaq[imax]),
               "U4_max": float(1 - np.mean(eta[:, imax] ** 4) / (3 * np.mean(eta[:, imax] ** 2) ** 2)),
               "eta": {str(k): float(etaq[i]) for i, k in enumerate(keys) if etaq[i] > 0.05}, "sro": {k: v / nsr for k, v in sr.items() if v != 0}}
        for name, q in (("RS", (1.0, 1.0, 1.0)), ("COL", (1.0, 1.0, 0)), ("LAY", (0, 0, 1.0))):
            if q in keys:
                i = keys.index(q)
                rec[f"eta_{name}"] = float(etaq[i])
                rec[f"U4_{name}"] = float(1 - np.mean(eta[:, i] ** 4) / (3 * np.mean(eta[:, i] ** 2) ** 2))
        out.append(rec)
    return out, s


if __name__ == "__main__":
    V = json.loads(sys.argv[1])
    L = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    Lz = int(sys.argv[3]) if len(sys.argv) > 3 else L
    Tmin = float(sys.argv[4]) if len(sys.argv) > 4 else 300
    Tmax = float(sys.argv[5]) if len(sys.argv) > 5 else 4000
    nT = int(sys.argv[6]) if len(sys.argv) > 6 else 30
    nsw = int(sys.argv[7]) if len(sys.argv) > 7 else 2000
    seed = int(sys.argv[8]) if len(sys.argv) > 8 else 0
    temps = np.linspace(Tmax, Tmin, nT)
    res, s = run(V, L, Lz, temps, nsw, nsw // 2, seed)
    for r in res:
        print("%6.0f E %8.2f C %6.3f etaRS %.3f etaLAY %.3f etaCOL %.3f max %s %.3f U4 %.3f" % (
            r["T"], r["E_meV_fu"], r["C_kB_fu"], r.get("eta_RS", -1), r.get("eta_LAY", -1), r.get("eta_COL", -1), r["q_max"], r["eta_max"], r["U4_max"]))
