"""Track L PBA taut (iii): K off-centre / K-ordering variants of the LCM, PBE+U base U, 110 Ry, symmetry ON.
usage: python kvar.py cfg.json   -> results/<id>.json
cfg: {id, MN, MC, cell (exact F-43m stage-1 cell), U {el: u}, ecut, kspacing, disp_A}
runs (ionic relax at the fixed stage-1 cell, LCM start V -3 / M +3; QE keeps the symmetry of the start):
  P1      15-atom ideal cell (reference; forces ~0)
  K111    K displaced by disp_A along [111] (toward the C-bound metal side of its void); C3v site, R3m
  K100    K displaced by disp_A along [100]; C2v site
  P2      2-f.u. supercell (a1 doubled), both K in the same tetrahedral void set (= ideal, reference)
  ALT2    2-f.u. supercell, second K moved to the other void set (T- = T+ + (1/2,1/2,1/2) primitive), i.e. K chains along x
E per f.u. relative to P1 / P2, M, gap and spin windows (SCF grid)."""
import json, os, sys, time, shutil, traceback, re
import numpy as np
from ase import Atoms
from ase.build import make_supercell
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, save_json, SCRATCH

NC = int(os.environ.get("NCORES", 64))
NPOOL = int(os.environ.get("NPOOL", 8))
cfg = json.load(open(sys.argv[1]))
os.makedirs("results", exist_ok=True)
RES = f"results/{cfg['id']}.json"
ORB = {"V": "3d", "Cr": "3d", "Mo": "4d"}
U = {el: (ORB[el], float(u)) for el, u in cfg["U"].items() if float(u) > 0}
out = json.load(open(RES)) if os.path.exists(RES) else {"id": cfg["id"], "runs": {}}
c = cfg["cell"]
P1 = Atoms(c["species"], cell=c["lattice"], scaled_positions=c["frac"], pbc=True)
iK = c["species"].index("K")


def glob_windows(b):
    E = b["eig"]; occ = b["occ"]
    occd = occ > 0.5
    vb = np.where(occd, E, -np.inf); cb = np.where(~occd, E, np.inf)
    tv = vb.max(axis=(0, 2)); bc = cb.min(axis=(0, 2))
    sv = int(np.argmax(tv)); sc = int(np.argmin(bc))
    return {"gap": float(bc.min() - tv.max()), "vbm_spin": sv + 1, "cbm_spin": sc + 1, "win_VB": float(tv[sv] - tv[1 - sv]),
            "win_CB": float(bc[1 - sc] - bc[sc])}


def variants():
    d = cfg.get("disp_A", 0.3)
    v = {"P1": P1.copy()}
    a = P1.copy(); p = a.get_positions(); p[iK] += d * np.array([1, 1, 1]) / np.sqrt(3); a.set_positions(p); v["K111"] = a
    a = P1.copy(); p = a.get_positions(); p[iK] += d * np.array([1, 0, 0]); a.set_positions(p); v["K100"] = a
    sc = make_supercell(P1, np.diag([2, 1, 1]))
    v["P2"] = sc.copy()
    alt = sc.copy()
    ks = [i for i, s in enumerate(alt.get_chemical_symbols()) if s == "K"]
    f = alt.get_scaled_positions()
    # K #2 (second image): shift by (1/2,1/2,1/2) primitive = (1/4, 1/2, 1/2) in supercell fractional coordinates
    f[ks[1]] = (f[ks[1]] + np.array([0.25, 0.5, 0.5])) % 1.0
    alt.set_scaled_positions(f)
    v["ALT2"] = alt
    return v


for name, at in variants().items():
    if name in out["runs"] and out["runs"][name].get("energy_eV") is not None:
        continue
    t0 = time.time()
    sym = at.get_chemical_symbols()
    mom = [(-3.0 if s == cfg["MN"] else (3.0 if s == cfg["MC"] else 0.0)) for s in sym]
    pref = f"k_{name}"
    write_pw(f"{pref}.in", at, calculation="relax", pset="dojo_sr", ecutwfc=cfg["ecut"], kspacing=cfg["kspacing"], moments=mom,
             smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref, conv_thr=1e-7 * len(at), tstress=False, tprnfor=True,
             electrons={"mixing_beta": 0.3}, control={"nstep": 60, "max_seconds": 3.0 * 3600})
    r = run_pw(f"{pref}.in", ncores=NC, npool=NPOOL, timeout=4 * 3600)
    p = parse_pw(r["out"])
    txt = open(r["out"], errors="replace").read()
    q = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag", "abs_mag", "bfgs_converged", "total_force", "errors", "n_scf_cycles")}
    q["nsym"] = (re.findall(r"(\d+) Sym\. Ops\.", txt) or [None])[0]
    q["nks"] = (re.findall(r"number of k points=\s*(\d+)", txt) or [None])[0]
    q["nfu"] = sym.count("K")
    try:
        from ase.io import read
        af = read(r["out"], format="espresso-out", index=-1)
        a0 = at.get_positions(); a1 = af.get_positions()
        if name in ("P1", "K111", "K100"):  # K offset from the ideal Td site after relaxation
            q["K_offset_final_A"] = float(np.linalg.norm(a1[iK] - P1.get_positions()[iK]))
            q["K_offset_vec_A"] = (a1[iK] - P1.get_positions()[iK]).tolist()
        q["K_move_A"] = [float(np.linalg.norm(a1[i] - a0[i])) for i, s in enumerate(sym) if s == "K"]
        q["max_move_A"] = float(np.abs(a1 - a0).max())
    except Exception as e:
        q["geo_err"] = repr(e)
    try:
        b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
        q["glob"] = glob_windows(b)
    except Exception as e:
        q["glob_err"] = repr(e)
    q["sec"] = time.time() - t0
    out["runs"][name] = q
    for ref, group in (("P1", ("P1", "K111", "K100")), ("P2", ("P2", "ALT2"))):
        if out["runs"].get(ref, {}).get("energy_eV") is not None:
            for g in group:
                rr = out["runs"].get(g)
                if rr and rr.get("energy_eV") is not None:
                    rr["dE_eV_per_fu"] = (rr["energy_eV"] / rr["nfu"]) - out["runs"][ref]["energy_eV"] / out["runs"][ref]["nfu"]
    save_json(out, RES)
    print(cfg["id"], name, q, flush=True)
    shutil.rmtree(SCRATCH / pref, ignore_errors=True)
out["t_end"] = time.time()
save_json(out, RES)
