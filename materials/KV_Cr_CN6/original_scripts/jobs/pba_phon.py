"""Track L PBA phonons: PBE+U (Gate-0 stage-1 protocol) forces on phonopy displaced supercells, LCM spins.

usage: python pba_phon.py <cfg.json> <idx> [<idx> ...]     (idx -1 = undisplaced reference)
  -> forces_<idx>.json in the job dir (skipped if it already exists: resumes after preemption)
Protocol: PseudoDojo SR NC, ECUT 110 Ry (ecutrho 440), ortho-atomic U from cfg, gaussian 0.005 Ry, mixing 0.3 (local-TF),
conv_thr CONV_PER_ATOM (2e-11) x nat Ry, k-grid cfg['kpts'] (mesh + shift), LCM start moments spin x |m0|.
NPOOL: largest divisor of NCORES <= min(NPOOL_MAX, 2 x irreducible k) (spglib, species = element + start moment).
env: NCORES (set by the backend) NPOOL_MAX (16) CONV_PER_ATOM (2e-11) MAXSEC (5 h per SCF)
"""
import json
import os
import re
import shutil
import sys
import time

import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, SCRATCH

NC = int(os.environ.get("NCORES", 32))
NPOOL_MAX = int(os.environ.get("NPOOL_MAX", 16))
CONV = float(os.environ.get("CONV_PER_ATOM", 2e-11))
MAXSEC = float(os.environ.get("MAXSEC", 5 * 3600))
RY_BOHR_TO_EV_A = 13.605693122994 / 0.529177210903


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def nk_irr(at, mom, mesh, shift):
    try:
        import spglib
        keys = {}
        types = [keys.setdefault((s, round(float(m), 2)), len(keys) + 1) for s, m in zip(at.get_chemical_symbols(), mom)]
        mp, _ = spglib.get_ir_reciprocal_mesh(mesh, (at.cell.array, at.get_scaled_positions(), types), is_shift=shift,
                                              is_time_reversal=True, symprec=1e-5)
        return int(len(np.unique(mp)))
    except Exception as e:
        log("spglib k count failed:", repr(e))
        return None


def pick_npool(nk):
    if nk is None:
        return 2
    best = 1
    for p in range(1, NC + 1):
        if NC % p == 0 and p <= 2 * nk and p <= NPOOL_MAX:
            best = p
    return best


def parse_forces(txt, nat):
    """last 'Forces acting on atoms' block (total forces, Ry/bohr) -> eV/A"""
    blocks = txt.split("Forces acting on atoms")
    if len(blocks) < 2:
        return None
    F = re.findall(r"atom\s+(\d+)\s+type\s+\d+\s+force\s+=\s+([-0-9.Ee+]+)\s+([-0-9.Ee+]+)\s+([-0-9.Ee+]+)", blocks[-1])
    if len(F) < nat:
        return None
    F = F[:nat]  # the first nat lines are the total force; later blocks are contributions
    arr = np.zeros((nat, 3))
    for i, fx, fy, fz in F:
        arr[int(i) - 1] = [float(fx), float(fy), float(fz)]
    return (arr * RY_BOHR_TO_EV_A).tolist()


def parse_hub(txt):
    out = {}
    for i, u, d, t in re.findall(r"Tr\[ns\(\s*(\d+)\)\]\s*\(up, down, total\)\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", txt):
        out[int(i)] = [float(u), float(d), float(t)]
    return out


def main():
    cfg = json.load(open(sys.argv[1]))
    idxs = [int(x) for x in sys.argv[2:]]
    U = {k: tuple(v) for k, v in cfg["hubbard"].items()}
    kp = cfg["kpts"]
    by_idx = {sc["idx"]: sc for sc in cfg["supercells"]}
    for idx in idxs:
        fout = f"forces_{idx}.json"
        if os.path.exists(fout):
            log("skip", idx, "(done)")
            continue
        sc = by_idx[idx]
        at = Atoms(sc["species"], cell=sc["lattice"], scaled_positions=sc["frac"], pbc=True)
        mom = [s * cfg["mom"].get(e, 0.0) for e, s in zip(sc["species"], cfg["spins"])]
        nk = nk_irr(at, mom, kp[:3], kp[3:])
        npool = pick_npool(nk)
        pref = f"{cfg['id']}_d{idx}".replace("-", "m")
        write_pw(f"{pref}.in", at, calculation="scf", pset="dojo_sr", ecutwfc=cfg.get("ecut", 110.0), kpts=kp[:3],
                 kshift=tuple(kp[3:]), moments=mom, smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref,
                 conv_thr=CONV * len(at), tstress=False, tprnfor=True,
                 electrons={"mixing_beta": 0.3, "electron_maxstep": 400})
        log("start", idx, "nat", len(at), "nk_irr(spglib)", nk, "npool", npool)
        r = run_pw(f"{pref}.in", ncores=NC, npool=npool, timeout=MAXSEC)
        txt = open(r["out"], errors="replace").read()
        p = parse_pw(r["out"])
        F = parse_forces(txt, len(at))
        hub = parse_hub(txt)
        nks = re.findall(r"number of k points=\s+(\d+)", txt)
        nsym = re.findall(r"(\d+) Sym\. Ops\.", txt)
        res = {"idx": idx, "forces_eV_A": F, "energy_eV": p.get("energy_eV"), "conv": p.get("converged_scf"),
               "M": p.get("total_mag"), "absM": p.get("abs_mag"), "site_moments": p.get("site_moments"),
               "hub_tr_ns": hub, "n_scf_iter": len(re.findall(r"iteration #", txt)), "rc": r["rc"],
               "sec": r["seconds"], "nk_qe": int(nks[-1]) if nks else None, "nsym_qe": int(nsym[-1]) if nsym else None,
               "npool": npool}
        if F is None or not p.get("converged_scf"):
            log("FAILED", idx, "rc", r["rc"], "conv", p.get("converged_scf"))
            json.dump(res, open(f"failed_{idx}.json", "w"))
        else:
            json.dump(res, open(fout, "w"))
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
        log("done", idx, "E", p.get("energy_eV"), "M", p.get("total_mag"), "|M|", p.get("abs_mag"), "sec", round(r["seconds"]),
            "nk", res["nk_qe"], "nsym", res["nsym_qe"], "maxF", None if F is None else float(np.abs(F).max()))


if __name__ == "__main__":
    main()
