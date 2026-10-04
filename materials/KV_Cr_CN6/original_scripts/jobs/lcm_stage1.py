"""Track L stage-1: QE PBE+U (PseudoDojo SR NC) on LCM candidates.
For each candidate: SCF for every config in c["configs"] ({name, spins}); FM always included.
Records energies, site moments, total/abs magnetization, gap and band-edge spin analysis on the SCF grid and a
denser nscf grid (for gapped compensated configs), incl. GLOBAL spin-polarization windows:
  win_VB = VBM - max_k(top valence band of the opposite spin)   (>0: the top win_VB of the VB is 100% one spin)
  win_CB = min_k(bottom conduction band of opposite spin) - CBM
and Gamma-point splitting of the edge bands.
usage: python lcm_stage1.py cands.json      (results/<id>.json)

pivot/ferrite copy (identical physics/defaults to infra/jobs_src/lcm_stage1/lcm_stage1.py md5 61529d21...), adds env knobs:
  NPOOL (scf/vc-relax pools, default NC//2), NPOOL_NSCF (default NC)  -- QE only warns when npool > nks and then idles
      cores; large cells with few irreducible k need fewer pools
  RELAX_MAXSEC (QE max_seconds for the vc-relax: pw.x stops cleanly and the last geometry is used; subprocess timeout is
      raised to RELAX_MAXSEC + 1 h), SCF_MAXSEC (same for each scf; default none)

pivot/garnet copy (= pivot/ferrite copy md5 b73a0c7a + two workflow-only changes; physics, QE inputs and defaults identical):
  * candidate field "nscf_configs" (list of config names): run the dense nscf only for those configs (default: every gapped
      non-FM config, as before). Used to skip nscf on the extra exchange-fit configurations.
  * vc-relax resume after container preemption: if <id>_vcrelax.out from an interrupted attempt exists in the job dir and
      no relaxed geometry was saved, it is renamed to <id>_vcrelax_prev<n>.out and the new vc-relax starts from its last
      complete CELL_PARAMETERS + ATOMIC_POSITIONS block; recorded as out["relax_restarts"]. Otherwise every preemption restarted the
      relax from the input geometry (seen on 4 of 5 pivot/ferrite run_cpu64/32 jobs, 2026-10-02)."""
import glob
import json, os, sys, time, traceback, shutil
import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, band_edge_analysis, save_json, SCRATCH

NC = int(os.environ.get("NCORES", 16))
ECUT = float(os.environ.get("ECUT", 75))
KSCF = float(os.environ.get("KSCF", 0.28))
KNSCF = float(os.environ.get("KNSCF", 0.16))
NPOOL = int(os.environ.get("NPOOL", max(1, NC // 2)))
NPOOL_NSCF = int(os.environ.get("NPOOL_NSCF", NC))
RELAX_MAXSEC = os.environ.get("RELAX_MAXSEC")
SCF_MAXSEC = os.environ.get("SCF_MAXSEC")
os.makedirs("results", exist_ok=True)
cands = json.load(open(sys.argv[1]))
UDEF = {"Fe": 4.0, "Mn": 4.0, "Co": 4.0, "Ni": 5.0, "Cr": 3.0, "V": 3.0, "Cu": 6.0, "Ti": 3.0,
        "Mo": 2.0, "Ru": 2.0, "Os": 2.0, "Re": 2.0, "Ir": 2.0, "W": 2.0, "Nb": 2.0, "Rh": 2.0}
ORB = {el: ("3d" if el in "Ti V Cr Mn Fe Co Ni Cu".split() else ("4d" if el in "Nb Mo Ru Rh".split() else "5d")) for el in UDEF}


def global_windows(b):
    E = b["eig"]; occ = b["occ"]
    if E.shape[1] != 2:
        return {}
    occd = occ > 0.5
    vb = np.where(occd, E, -np.inf); cb = np.where(~occd, E, np.inf)
    tv = vb.max(axis=(0, 2)); bc = cb.min(axis=(0, 2))  # per spin
    sv = int(np.argmax(tv)); sc = int(np.argmin(bc))
    out = {"vbm_spin": sv, "cbm_spin": sc, "win_VB": float(tv[sv] - tv[1 - sv]), "win_CB": float(bc[1 - sc] - bc[sc]),
           "gap_up": float(bc[0] - tv[0]), "gap_dn": float(bc[1] - tv[1])}
    # Gamma point (k = 0) splitting of top VB / bottom CB
    k = b["k"]
    ig = np.where(np.linalg.norm(k, axis=1) < 1e-6)[0]
    if len(ig):
        g = ig[0]
        out["gamma_top_vb"] = [float(vb[g, 0].max()), float(vb[g, 1].max())]
        out["gamma_bot_cb"] = [float(cb[g, 0].min()), float(cb[g, 1].min())]
    return out


def last_geometry(txt, nat):
    """Last COMPLETE (CELL_PARAMETERS, ATOMIC_POSITIONS crystal) block of a (possibly truncated) vc-relax output.
    Returns (cell[3x3] Angstrom, frac[nat x 3], n_blocks) or None if no complete block exists."""
    BOHR = 0.529177210903
    lines = txt.splitlines()
    best, nblk = None, 0
    for i, ln in enumerate(lines):
        if not ln.startswith("CELL_PARAMETERS"):
            continue
        try:
            cell = np.array([[float(x) for x in lines[i + k].split()[:3]] for k in (1, 2, 3)])
            if "angstrom" in ln:
                fac = 1.0
            elif "alat" in ln:
                fac = float(ln.split("=")[1].strip(" )")) * BOHR
            else:  # bohr
                fac = BOHR
            j = i + 4
            while j < len(lines) and not lines[j].startswith("ATOMIC_POSITIONS"):
                j += 1
            if j >= len(lines) or "crystal" not in lines[j]:
                continue
            fr = np.array([[float(x) for x in lines[j + 1 + k].split()[1:4]] for k in range(nat)])
            if fr.shape != (nat, 3):
                continue
            best, nblk = (cell * fac, fr), nblk + 1
        except (IndexError, ValueError):
            continue
    return None if best is None else (best[0], best[1], nblk)


def run_one(c):
    out = {"id": c["id"], "formula": c["formula"], "t0": time.time(), "runs": {}}
    if os.path.exists(f"results/{c['id']}.json"):
        try:
            prev = json.load(open(f"results/{c['id']}.json"))
            out["runs"] = prev.get("runs", {})
            for k in ("relaxed", "relax_parse"):
                if k in prev:
                    out[k] = prev[k]
        except Exception:
            pass
    at = Atoms(c["species"], cell=c["lattice"], scaled_positions=c["frac"], pbc=True)
    mags = np.abs(np.array(c.get("magmoms") or [0] * len(at), float))
    U = c.get("hubbard")
    if U is None:
        U = {el: (ORB[el], UDEF[el]) for el in set(c["species"]) if el in UDEF}
    out["hubbard"] = U
    configs = c["configs"]
    if out.get("relaxed"):
        at = Atoms(out["relaxed"]["species"], cell=out["relaxed"]["lattice"], scaled_positions=out["relaxed"]["frac"], pbc=True)
    elif c.get("relax_config"):
        # optional vc-relax (PBE+U, same settings) in the given magnetic configuration before the order scan
        rc = [x for x in configs if x["name"] == c["relax_config"]][0]
        mom = [sg * max(m, 3.0) if sg != 0 else 0.0 for sg, m in zip(rc["spins"], mags)]
        pref = f"{c['id']}_vcrelax"
        if os.path.exists(f"{pref}.out"):
            # interrupted attempt (preemption / timeout kill): continue from its last geometry
            try:
                g = last_geometry(open(f"{pref}.out", errors="replace").read(), len(at))
                if g is not None:
                    nprev = len(glob.glob(f"{pref}_prev*.out"))
                    shutil.move(f"{pref}.out", f"{pref}_prev{nprev}.out")
                    at = Atoms(c["species"], cell=g[0], scaled_positions=g[1], pbc=True)
                    out["relax_restarts"] = nprev + 1
                    out["relax_restart_nsteps"] = g[2]
            except Exception as e:
                out["relax_restart_err"] = repr(e)
        write_pw(f"{pref}.in", at, calculation="vc-relax", pset="dojo_sr", ecutwfc=ECUT, kspacing=KSCF, moments=mom,
                 smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref, conv_thr=1e-7 * len(at),
                 electrons={"mixing_beta": 0.3},
                 control=dict({"nstep": 120}, **({"max_seconds": float(RELAX_MAXSEC)} if RELAX_MAXSEC else {})))
        r = run_pw(f"{pref}.in", ncores=NC, npool=NPOOL,
                   timeout=(float(RELAX_MAXSEC) + 3600) if RELAX_MAXSEC else 8 * 3600)
        from qeutil import final_atoms
        try:
            at = final_atoms(r["out"])
            out["relaxed"] = {"lattice": at.cell.array.tolist(), "species": at.get_chemical_symbols(), "frac": at.get_scaled_positions().tolist()}
            out["relax_parse"] = {k: v for k, v in parse_pw(r["out"]).items() if k in ("energy_eV", "pressure_kbar", "bfgs_converged", "total_mag")}
            out["relax_parse"]["sec"] = r["seconds"]
        except Exception as e:
            out["relax_err"] = repr(e)
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
        save_json(out, f"results/{c['id']}.json")
    if not any(x["name"] == "FM" for x in configs):
        fm = [1 if s != 0 else 0 for s in configs[0]["spins"]]
        configs = [{"name": "FM", "spins": fm}] + configs
    for cf in configs:
        name, spins = cf["name"], cf["spins"]
        if name in out["runs"] and out["runs"][name].get("energy_eV") is not None:
            continue
        mom = [sg * max(m, 3.0) if sg != 0 else 0.0 for sg, m in zip(spins, mags)]
        pref = f"{c['id']}_{name}"
        inp = f"{pref}.in"
        write_pw(inp, at, calculation="scf", pset="dojo_sr", ecutwfc=ECUT, kspacing=KSCF, moments=mom,
                 smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref, conv_thr=1e-7 * len(at),
                 tstress=False, tprnfor=False, electrons={"mixing_beta": 0.3},
                 control=({"max_seconds": float(SCF_MAXSEC)} if SCF_MAXSEC else None))
        r = run_pw(inp, ncores=NC, npool=NPOOL, timeout=(float(SCF_MAXSEC) + 3600) if SCF_MAXSEC else 5 * 3600)
        p = parse_pw(r["out"]); p["sec"] = r["seconds"]
        try:
            b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
            p["edges_scf"] = band_edge_analysis(b)
            p["glob_scf"] = global_windows(b)
        except Exception as e:
            p["edges_scf_err"] = repr(e)
        sm = p.get("site_moments") or []
        p["n_flipped"] = int(sum(1 for i, sg in enumerate(spins) if sg != 0 and i < len(sm) and abs(sm[i]) > 0.5 and np.sign(sm[i]) != sg))
        p["n_collapsed"] = int(sum(1 for i, sg in enumerate(spins) if sg != 0 and i < len(sm) and abs(sm[i]) < 0.5))
        want_nscf = (name in c["nscf_configs"]) if c.get("nscf_configs") is not None else (name != "FM")
        if want_nscf and p.get("edges_scf", {}).get("gap", 0) > 0.05 and p.get("converged_scf"):
            ninp = f"{pref}_nscf.in"
            write_pw(ninp, at, calculation="nscf", pset="dojo_sr", ecutwfc=ECUT, kspacing=KNSCF, moments=mom,
                     smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref, tstress=False, tprnfor=False)
            r2 = run_pw(ninp, ncores=NC, npool=NPOOL_NSCF, timeout=4 * 3600)
            p["nscf_sec"] = r2["seconds"]
            try:
                b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
                p["edges_nscf"] = band_edge_analysis(b)
                p["glob_nscf"] = global_windows(b)
            except Exception as e:
                p["edges_nscf_err"] = repr(e)
        p["spins"] = spins
        out["runs"][name] = p
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
        save_json(out, f"results/{c['id']}.json")
    nmag = sum(1 for s in configs[0]["spins"] if s)
    eFM = out["runs"]["FM"].get("energy_eV")
    for k, v in out["runs"].items():
        if v.get("energy_eV") is not None and eFM is not None:
            v["dE_vs_FM_meV_per_mag"] = (v["energy_eV"] - eFM) * 1000 / max(nmag, 1)
    out["nmag"] = nmag
    out["t1"] = time.time()
    save_json(out, f"results/{c['id']}.json")
    return out


for c in cands:
    if os.path.exists(f"results/{c['id']}.json") and json.load(open(f"results/{c['id']}.json")).get("t1"):
        continue
    try:
        o = run_one(c)
        for k, v in o["runs"].items():
            g = v.get("glob_nscf", v.get("glob_scf", {}))
            print(c["id"], c["formula"], k, "dE=%.1f" % v.get("dE_vs_FM_meV_per_mag", 0), "M=%s" % v.get("total_mag"),
                  "gap=%.2f" % v.get("edges_scf", {}).get("gap", -1), "winVB=%.2f winCB=%.2f" % (g.get("win_VB", 0), g.get("win_CB", 0)),
                  "flip=%d col=%d" % (v["n_flipped"], v["n_collapsed"]), flush=True)
    except Exception:
        traceback.print_exc()
