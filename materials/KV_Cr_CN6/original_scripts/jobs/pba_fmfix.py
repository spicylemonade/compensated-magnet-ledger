"""Track L PBA Gate 0: fixed-moment FM reference SCFs (tot_magnetization = 6, i.e. the d3/d3 Heisenberg FM).
Why: the unconstrained FM is not the Heisenberg reference in these cyanides.
  * plain PBE (initial stage of the HSE runs): KVMo FM converges to M = 2.57 (V t2g up -> Mo t2g down charge transfer,
    metallic) and the KVCr FM start collapses to M = 0
  * PBE+U at U_Mo = 0: M_FM = 4.91 / 5.91
With M fixed at 6 muB/cell (two Fermi levels) the transfer is forbidden; E(FMfix) - E(LCM) is the spin-flip energy
used for the J ratios.
usage: python pba_fmfix.py cands.json   (list of {id, lattice, species, frac, spins, tot_mag, ecut, kspacing,
       hubbard (PBE+U; {} = none) | hse: true + hse_nq, moment_mag})
output: results/<id>.json  {energy_eV, converged_scf, total_mag, abs_mag, site_moments, glob, homo_lumo, sec, settings}
Settings mirror lcm_stage1 (PBE+U: conv_thr 1e-7*nat, mixing 0.3) and lcm_deep hse (HSE06: ecutfock 2*ecut, gygi-
baldereschi, x_gamma_extrapolation, npool NC//8)."""
import json, os, sys, shutil, traceback
import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, save_json, SCRATCH

NC = int(os.environ.get("NCORES", 32))
os.makedirs("results", exist_ok=True)


def glob_windows(b):
    E = b["eig"]; occ = b["occ"]
    occd = occ > 0.5
    vb = np.where(occd, E, -np.inf); cb = np.where(~occd, E, np.inf)
    if E.shape[1] != 2:
        return {"gap": float(cb.min() - vb.max())}
    tv = vb.max(axis=(0, 2)); bc = cb.min(axis=(0, 2))
    sv = int(np.argmax(tv)); sc = int(np.argmin(bc))
    return {"gap": float(bc.min() - tv.max()), "vbm_spin": sv, "cbm_spin": sc, "win_VB": float(tv[sv] - tv[1 - sv]),
            "win_CB": float(bc[1 - sc] - bc[sc]), "gap_up": float(bc[0] - tv[0]), "gap_dn": float(bc[1] - tv[1])}


for c in json.load(open(sys.argv[1])):
    path = f"results/{c['id']}.json"
    if os.path.exists(path) and json.load(open(path)).get("converged_scf"):
        continue
    try:
        at = Atoms(c["species"], cell=c["lattice"], scaled_positions=c["frac"], pbc=True)
        mom = [s * c.get("moment_mag", 3.0) for s in c["spins"]]
        sysx = {"tot_magnetization": float(c["tot_mag"])}
        if c.get("hse"):
            nq = c["hse_nq"]
            sysx.update({"input_dft": "hse", "nqx1": nq[0], "nqx2": nq[1], "nqx3": nq[2], "ecutfock": 2 * c["ecut"],
                         "exxdiv_treatment": "gygi-baldereschi", "x_gamma_extrapolation": True})
            U, conv, npool, elec, tmo = None, 1e-8 * len(at), max(1, NC // 8), {"mixing_beta": 0.3, "electron_maxstep": 200}, 22 * 3600
        else:
            U = {k: tuple(v) for k, v in (c.get("hubbard") or {}).items()} or None
            conv, npool, elec, tmo = 1e-7 * len(at), int(os.environ.get("NPOOL", max(1, NC // 4))), {"mixing_beta": 0.3}, 5 * 3600
        pref = c["id"]
        write_pw(f"{pref}.in", at, pset="dojo_sr", ecutwfc=c["ecut"], kspacing=c["kspacing"], moments=mom, smearing="gaussian",
                 degauss=0.005, hubbard=U, prefix=pref, conv_thr=conv, tstress=False, tprnfor=False, system=sysx, electrons=elec)
        r = run_pw(f"{pref}.in", ncores=NC, npool=npool, timeout=tmo)
        p = parse_pw(r["out"]); p["sec"] = r["seconds"]
        try:
            b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
            p["glob"] = glob_windows(b)
        except Exception as e:
            p["glob_err"] = repr(e)
        out = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag", "abs_mag", "site_moments", "glob", "homo_lumo", "sec", "errors")}
        out["settings"] = {k: c.get(k) for k in ("tot_mag", "ecut", "kspacing", "hubbard", "hse", "hse_nq")}
        save_json(out, path)
        print(c["id"], out.get("energy_eV"), out.get("total_mag"), out.get("converged_scf"), (out.get("glob") or {}).get("gap"), flush=True)
    except Exception:
        traceback.print_exc()
    shutil.rmtree(SCRATCH / c["id"], ignore_errors=True)
