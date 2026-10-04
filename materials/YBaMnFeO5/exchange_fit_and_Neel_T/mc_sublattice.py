"""Finite-T sublattice magnetizations for the LCM order (classical Heisenberg, Metropolis, numba).
Tracks m_Mn(T), m_Fe(T) (projected on the instantaneous staggered axis) and the net moment
M_net = 5 muB * (N_Mn m_Mn - N_Fe m_Fe)/N_pairs (spin-only, S=5/2, g=2 classical-spin mapping).
usage: python tracks/lcm/mc_sublattice.py '<J json>' L"""
import json, sys, math
import numpy as np, numba
sys.path.insert(0, "tracks/lcm")
from pymatgen.core import Structure, Lattice
exec(open("tracks/lcm/tn_mc_dir.py").read().split("rl = json.load")[0].split("from pymatgen.core import Structure, Lattice")[1]) if False else None

def cls(si, sj, d, dz):
    pair = "".join(sorted((si, sj)))
    if pair == "FeMn":
        if dz < 0.6 and d < 4.3: return "ip"
        if dz > 3.8 and d < 4.3: return "ap"
        if 3.3 < dz < 3.8 and d < 3.9: return "Y"
        if d < 7.0: return f"FeMn_{round(d, 1)}"
        return None
    if d < 6.0:
        return f"{pair}_{'ip' if dz < 0.6 else 'oop'}_{round(d, 1)}"
    return None

@numba.njit(cache=True)
def run(S, NB, JJ, sub, T, nsweep, nequil, seed):
    np.random.seed(seed)
    N = S.shape[0]; beta = 1.0 / (0.08617333262 * T)
    mA = 0.0; mB = 0.0; mst = 0.0; cnt = 0
    for sw in range(nsweep + nequil):
        for _ in range(N):
            i = np.random.randint(N)
            hx = 0.0; hy = 0.0; hz = 0.0
            for b in range(NB.shape[1]):
                j = NB[i, b]
                if j < 0: break
                hx += JJ[i, b] * S[j, 0]; hy += JJ[i, b] * S[j, 1]; hz += JJ[i, b] * S[j, 2]
            z = 2.0 * np.random.random() - 1.0; ph = 2.0 * np.pi * np.random.random(); r = math.sqrt(1 - z * z)
            nx = r * math.cos(ph); ny = r * math.sin(ph); nz = z
            dE = -(hx * (nx - S[i, 0]) + hy * (ny - S[i, 1]) + hz * (nz - S[i, 2]))
            if dE <= 0.0 or np.random.random() < math.exp(-beta * dE):
                S[i, 0] = nx; S[i, 1] = ny; S[i, 2] = nz
        if sw >= nequil:
            ax = 0.0; ay = 0.0; az = 0.0; bx = 0.0; by = 0.0; bz = 0.0; na = 0; nb = 0
            for i in range(N):
                if sub[i] > 0:
                    ax += S[i, 0]; ay += S[i, 1]; az += S[i, 2]; na += 1
                else:
                    bx += S[i, 0]; by += S[i, 1]; bz += S[i, 2]; nb += 1
            # staggered axis
            sx = ax - bx; sy = ay - by; sz = az - bz
            sn = math.sqrt(sx * sx + sy * sy + sz * sz) + 1e-12
            ux = sx / sn; uy = sy / sn; uz = sz / sn
            mA += (ax * ux + ay * uy + az * uz) / na
            mB += -(bx * ux + by * uy + bz * uz) / nb
            mst += sn / N
            cnt += 1
    return mA / cnt, mB / cnt, mst / cnt

if __name__ == "__main__":
    J = json.loads(sys.argv[1]); L = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    rl = json.load(open("tracks/lcm/batch4/YBaMnFeO5_rs_relaxed.json"))
    prim = Structure(Lattice(rl["lattice"]), rl["species"], rl["frac"])
    sc = prim.copy(); sc.make_supercell([L, L, max(2, L // 2)])
    idx = [i for i, s in enumerate(sc) if s.specie.symbol in ("Mn", "Fe")]
    pos = {j: k for k, j in enumerate(idx)}
    nbrs = [[] for _ in idx]
    for i in idx:
        for nb in sc.get_neighbors(sc[i], 7.0):
            j = nb.index
            if j not in pos: continue
            k = cls(sc[i].specie.symbol, sc[j].specie.symbol, nb.nn_distance, abs(nb.coords[2] - sc[i].coords[2]))
            if k in J and abs(J[k]) > 1e-9: nbrs[pos[i]].append((pos[j], J[k]))
    maxn = max(len(x) for x in nbrs)
    NB = -np.ones((len(idx), maxn), dtype=np.int64); JJ = np.zeros((len(idx), maxn))
    for a, lst in enumerate(nbrs):
        for b, (j, v) in enumerate(lst): NB[a, b] = j; JJ[a, b] = v
    sub = np.array([1 if sc[i].specie.symbol == "Mn" else -1 for i in idx], dtype=np.int64)
    S = np.zeros((len(idx), 3)); S[:, 2] = sub.astype(float)
    temps = [float(t) for t in (sys.argv[3].split(",") if len(sys.argv) > 3 else "50,150,250,300,350,400,450,500")]
    out = []
    for k, T in enumerate(temps):
        mA, mB, ms = run(S, NB, JJ, sub, T, 4000, 2000, 11 + k)
        out.append({"T": T, "m_Mn": round(mA, 4), "m_Fe": round(mB, 4), "m_stag": round(ms, 4), "M_net_muB_per_fu": round(5.0 * (mA - mB), 4)})
        print(out[-1], flush=True)
    json.dump({"J": J, "L": L, "rows": out}, open(f"tracks/lcm/tn/mc_sublattice_{abs(hash(sys.argv[1])) % 10000}.json", "w"))
