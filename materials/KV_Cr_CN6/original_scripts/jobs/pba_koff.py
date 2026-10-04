"""Track L PBA: K off-centring test in the 15-atom F-43m primitive cell (all K shifted in phase = Gamma T2 polar mode).

usage: python pba_koff.py <cfg.json> <tag>          cfg = inputs/<m>_prim.json of the phonon build (uses supercell idx -1)
  tag: x | xy | p111 | m111 | c (centred reference relax)  -> results/koff_<m>_<tag>.json
K starts displaced by DISP0 (0.25 A) along [100], [110], +[111] (toward the C-bound M corner) or -[111] (toward V);
'c' relaxes the centred cell (symmetry F-43m, so only C/N can move: a consistency check on the reference energy).
Fixed-cell ionic relax (BFGS, symmetry of the start kept: C2v / Cs / C3v), same protocol as the phonon forces
(PBE+U ortho-atomic, 110 Ry, gaussian 0.005 Ry, k 5x5x5, conv_thr 2e-11*nat), forc_conv_thr FORC (2e-4 Ry/bohr).
Reports E_final, the K displacement from the 4c centre, K-N/C shortest distances, M, |M|, and the per-spin band edges
on the SCF grid (windows, lane definition) from the final data-file-schema.xml.
env: DISP0 (0.25) FORC (2e-4) NPOOL (16) MAXSEC (5 h)
"""
import json
import os
import shutil
import sys
import time

import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, SCRATCH

NC = int(os.environ.get("NCORES", 32))
DISP0 = float(os.environ.get("DISP0", 0.25))
FORC = float(os.environ.get("FORC", 2e-4))
NPOOL = int(os.environ.get("NPOOL", 16))
MAXSEC = float(os.environ.get("MAXSEC", 5 * 3600))
DIRS = {"x": [1, 0, 0], "xy": [1, 1, 0], "p111": [1, 1, 1], "m111": [-1, -1, -1], "c": [0, 0, 0]}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def last_positions(txt, nat):
    lines = txt.splitlines()
    best = None
    for i, ln in enumerate(lines):
        if ln.startswith("ATOMIC_POSITIONS") and "crystal" in ln:
            try:
                fr = np.array([[float(x) for x in lines[i + 1 + k].split()[1:4]] for k in range(nat)])
                best = fr
            except (IndexError, ValueError):
                pass
    return best


def windows(xmlf):
    b = parse_xml_bands(xmlf)
    E, O = b["eig"], b["occ"]  # [nk, 2, nb]
    out = {}
    for s in (0, 1):
        occ = E[:, s, :][O[:, s, :] > 0.5]; emp = E[:, s, :][O[:, s, :] <= 0.5]
        out[s] = (float(occ.max()), float(emp.min()))
    vbm_s = 0 if out[0][0] >= out[1][0] else 1
    cbm_s = 0 if out[0][1] <= out[1][1] else 1
    vbm, cbm = out[vbm_s][0], out[cbm_s][1]
    return {"vbm_spin": vbm_s, "cbm_spin": cbm_s, "gap": cbm - vbm, "win_VB": vbm - out[1 - vbm_s][0],
            "win_CB": out[1 - cbm_s][1] - cbm, "edges": {str(s): out[s] for s in out}}


def main():
    cfg = json.load(open(sys.argv[1])); tag = sys.argv[2]
    m = cfg["m"]
    sc = [c for c in cfg["supercells"] if c["idx"] == -1][0]
    lat = np.array(sc["lattice"]); fr = np.array(sc["frac"]); sp = sc["species"]
    iK = sp.index("K")
    d = np.array(DIRS[tag], float)
    if np.linalg.norm(d) > 0:
        d = d / np.linalg.norm(d) * DISP0
    cart = fr @ lat
    centre = cart[iK].copy()
    cart[iK] += d
    at = Atoms(sp, cell=lat, positions=cart, pbc=True)
    mom = [s * cfg["mom"].get(e, 0.0) for e, s in zip(sp, cfg["spins"])]
    U = {k: tuple(v) for k, v in cfg["hubbard"].items()}
    kp = cfg["kpts"]
    pref = f"koff_{m}_{tag}"
    os.makedirs("results", exist_ok=True)
    if os.path.exists(f"{pref}.out"):  # preemption restart from the last positions
        g = last_positions(open(f"{pref}.out", errors="replace").read(), len(at))
        n = len([f for f in os.listdir(".") if f.startswith(f"{pref}_prev")])
        shutil.move(f"{pref}.out", f"{pref}_prev{n}.out")
        if g is not None:
            at = Atoms(sp, cell=lat, scaled_positions=g, pbc=True)
            log("restart from previous positions")
    write_pw(f"{pref}.in", at, calculation="relax", pset="dojo_sr", ecutwfc=cfg.get("ecut", 110.0), kpts=kp[:3],
             kshift=tuple(kp[3:]), moments=mom, smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref,
             conv_thr=2e-11 * len(at), tstress=False, tprnfor=True,
             control={"forc_conv_thr": FORC, "etot_conv_thr": 1e-6 * len(at), "nstep": 200},
             electrons={"mixing_beta": 0.3, "electron_maxstep": 400})
    log("relax", m, tag, "K start disp", d.round(3).tolist())
    r = run_pw(f"{pref}.in", ncores=NC, npool=NPOOL, timeout=MAXSEC)
    txt = open(r["out"], errors="replace").read()
    p = parse_pw(r["out"])
    g = last_positions(txt, len(at))
    res = {"m": m, "tag": tag, "disp0_A": d.tolist(), "rc": r["rc"], "sec": r["seconds"], "energy_eV": p.get("energy_eV"),
           "bfgs_converged": "bfgs converged" in txt, "n_scf": p.get("n_scf_cycles"), "M": p.get("total_mag"),
           "absM": p.get("abs_mag"), "site_moments": p.get("site_moments")}
    if g is not None:
        c2 = g @ lat
        dK = c2[iK] - centre
        res["K_disp_A"] = dK.round(4).tolist(); res["K_disp_norm_A"] = float(np.linalg.norm(dK))
        dist = []
        for j, s in enumerate(sp):
            if s in ("C", "N"):
                v = (g[j] - g[iK]); v -= np.round(v); dist.append((float(np.linalg.norm(v @ lat)), s))
        dist.sort()
        res["K_nearest"] = [(round(x, 3), s) for x, s in dist[:6]]
        res["final_frac"] = g.tolist()
    xmlf = SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"
    try:
        res["edges_scf"] = windows(xmlf)
    except Exception as e:
        res["edges_err"] = repr(e)[:200]
    json.dump(res, open(f"results/{pref}.json", "w"), indent=1)
    shutil.rmtree(SCRATCH / pref, ignore_errors=True)
    log("done", json.dumps({k: res.get(k) for k in ("energy_eV", "bfgs_converged", "K_disp_norm_A", "M", "absM")}))


if __name__ == "__main__":
    main()
