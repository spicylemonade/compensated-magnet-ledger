"""Track L deep validation for a finalist (QE 7.5, PseudoDojo NC).
input JSON: {id, formula, lattice, species, frac, spins (LCM order, +-1/0), competitors: [{name, spins}], hubbard: {el: [orb, U]}}
modes (argv[2], comma-separated):
  uscan  : PBE+U SCF for LCM + competitors at U scale factors {0.5, 1.0, 1.5} (and U=0) -> dE, gap, spin windows, moments
  soc    : noncollinear + SOC SCF (FR pseudos) for LCM with moments along z and along x -> total spin+orbital moment, MAE
  hse    : HSE06 SCF (nq=1-2, ecutfock = 2*ecutwfc) for LCM -> gap, global spin windows (from SCF k-grid)
  bands  : PBE+U nscf on a spin-resolved high-symmetry path (seekpath) + dense uniform grid -> edges, effective masses
outputs deep_<mode>.json"""
import json, os, sys, time, shutil
import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, band_edge_analysis, save_json, SCRATCH, kmesh

NC = int(os.environ.get("NCORES", 32))
cfg = json.load(open(sys.argv[1]))
modes = sys.argv[2].split(",")
ECUT = cfg.get("ecut", 75.0); KSP = cfg.get("kspacing", 0.25)
at = Atoms(cfg["species"], cell=cfg["lattice"], scaled_positions=cfg["frac"], pbc=True)
U0 = {k: tuple(v) for k, v in cfg["hubbard"].items()}
MOM = cfg.get("moment_mag", 4.0)


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


def scf(pref, spins, U, extra_sys=None, pset="dojo_sr", noncolin=False, lspinorb=False, angles=(0.0, 0.0), kspacing=None, nbnd=None, elec=None, timeout=12 * 3600, npool=None):
    if noncolin:
        mom = [(s * MOM, angles[0], angles[1]) if s != 0 else (0.0, 0.0, 0.0) for s in spins]
    else:
        mom = [s * MOM for s in spins]
    write_pw(f"{pref}.in", at, pset=pset, ecutwfc=ECUT, kspacing=kspacing or KSP, moments=mom, smearing="gaussian", degauss=0.005,
             hubbard=U, prefix=pref, conv_thr=1e-8 * len(at), tstress=False, tprnfor=False, noncolin=noncolin, lspinorb=lspinorb,
             system=extra_sys, nbnd=nbnd, electrons=elec or {"mixing_beta": 0.3})
    r = run_pw(f"{pref}.in", ncores=NC, npool=npool or max(1, NC // 4), timeout=timeout)
    p = parse_pw(r["out"]); p["sec"] = r["seconds"]
    try:
        b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
        p["glob"] = glob_windows(b)
    except Exception as e:
        p["glob_err"] = repr(e)
    return p


res = {"id": cfg["id"]}
if "uscan" in modes:
    out = json.load(open("deep_uscan.json")) if os.path.exists("deep_uscan.json") else {}
    for fac in cfg.get("u_factors", [0.0, 0.5, 1.0, 1.5]):
        if str(fac) in out and len(out[str(fac)]) == 1 + len(cfg.get("competitors", [])):
            continue
        U = {k: (v[0], round(v[1] * fac, 3)) for k, v in U0.items()} if fac > 0 else None
        row = {}
        for nm, sp in [("LCM", cfg["spins"])] + [(c["name"], c["spins"]) for c in cfg.get("competitors", [])]:
            pref = f"u{fac}_{nm}"
            p = scf(pref, sp, U)
            row[nm] = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag", "abs_mag", "site_moments", "glob", "sec")}
            shutil.rmtree(SCRATCH / pref, ignore_errors=True)
            out[str(fac)] = row
            save_json(out, "deep_uscan.json")
        print("U factor", fac, {k: (round(v["energy_eV"] - row["LCM"]["energy_eV"], 4) if v.get("energy_eV") and row["LCM"].get("energy_eV") else None, (v.get("glob") or {}).get("gap")) for k, v in row.items()}, flush=True)
    res["uscan"] = out
if "soc" in modes:
    out = {}
    for tag, ang in [("z", (0.0, 0.0)), ("x", (90.0, 0.0))]:
        pref = f"soc_{tag}"
        p = scf(pref, cfg["spins"], U0, pset="dojo_fr", noncolin=True, lspinorb=True, angles=ang, elec={"mixing_beta": 0.2, "electron_maxstep": 400})
        out[tag] = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag_vec", "abs_mag", "site_moments", "glob", "sec")}
        # orbital moments are printed only with verbosity high in some versions; grep them
        try:
            txt = open(f"{pref}.out").read()
            import re
            orb = re.findall(r"orbital magnetization.*?=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", txt)
            out[tag]["orbital_mag_raw"] = orb[-3:]
        except Exception:
            pass
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
        save_json(out, "deep_soc.json")
        print("SOC", tag, out[tag].get("energy_eV"), out[tag].get("total_mag_vec"), flush=True)
    if out["z"].get("energy_eV") and out["x"].get("energy_eV"):
        out["MAE_meV_cell_x_minus_z"] = (out["x"]["energy_eV"] - out["z"]["energy_eV"]) * 1000
    res["soc"] = out
if "hse" in modes:
    nq = cfg.get("hse_nq", [1, 1, 1])
    sysx = {"input_dft": "hse", "nqx1": nq[0], "nqx2": nq[1], "nqx3": nq[2], "ecutfock": cfg.get("ecutfock", 2 * ECUT), "exxdiv_treatment": "gygi-baldereschi", "x_gamma_extrapolation": True}
    pref = "hse_lcm"
    p = scf(pref, cfg["spins"], None if cfg.get("hse_noU", True) else U0, extra_sys=sysx, kspacing=cfg.get("hse_kspacing", 0.30), elec={"mixing_beta": 0.3, "electron_maxstep": 200}, timeout=22 * 3600, npool=max(1, NC // 8))
    res["hse"] = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag", "abs_mag", "site_moments", "glob", "sec", "homo_lumo")}
    save_json(res["hse"], "deep_hse.json")
    print("HSE", res["hse"], flush=True)
if "bands" in modes:
    import seekpath
    struct = (at.cell.array, at.get_scaled_positions(), at.get_atomic_numbers())
    sp = seekpath.get_explicit_k_path(struct, reference_distance=0.025)
    # keep the ORIGINAL cell (no standardization) by transforming k-points back; simpler: use seekpath on original cell via
    # its 'primitive_transformation_matrix'; to stay safe we compute bands on seekpath's primitive cell with spins mapped.
    prim_lat = np.array(sp["primitive_lattice"]); prim_pos = np.array(sp["primitive_positions"]); prim_types = sp["primitive_types"]
    # map spins by matching positions
    cart_o = at.get_positions(); cart_p = prim_pos @ prim_lat
    spins_p = []
    from ase.data import chemical_symbols
    for cp in cart_p:
        fr = np.linalg.solve(at.cell.array.T, cp) % 1.0
        d = at.get_scaled_positions() - fr; d -= np.round(d)
        j = int(np.argmin(np.linalg.norm(d @ at.cell.array, axis=1)))
        spins_p.append(cfg["spins"][j])
    atp = Atoms([chemical_symbols[t] for t in prim_types], cell=prim_lat, scaled_positions=prim_pos, pbc=True)
    global_at = at
    at = atp
    p = scf("bnd", spins_p, U0)
    kpts = sp["explicit_kpoints_rel"]; labels = sp["explicit_kpoints_labels"]
    klines = "K_POINTS crystal\n%d\n" % len(kpts) + "\n".join("%.8f %.8f %.8f 1.0" % tuple(k) for k in kpts)
    mom = [s * MOM for s in spins_p]
    write_pw("bnd_bands.in", at, calculation="bands", pset="dojo_sr", ecutwfc=ECUT, kpts=klines, moments=mom, smearing="gaussian", degauss=0.005,
             hubbard=U0, prefix="bnd", tstress=False, tprnfor=False, nbnd=cfg.get("nbnd_bands"))
    r = run_pw("bnd_bands.in", ncores=NC, npool=max(1, NC // 4), timeout=8 * 3600)
    b = parse_xml_bands(str(SCRATCH / "bnd" / "bnd.save" / "data-file-schema.xml"))
    res["bands"] = {"scf": {k: p.get(k) for k in ("energy_eV", "total_mag", "glob")}, "labels": labels, "kpts": np.array(kpts).tolist(),
                    "eig": b["eig"].tolist(), "occ_scf_ef": p.get("fermi_eV"), "homo_scf": (p.get("homo_lumo") or [None])[0],
                    "prim_lattice": prim_lat.tolist(), "prim_species": [chemical_symbols[t] for t in prim_types], "prim_frac": prim_pos.tolist(), "prim_spins": spins_p}
    save_json(res["bands"], "deep_bands.json")
    at = global_at
save_json(res, f"deep_{'_'.join(modes)}.json")
print("done", modes)
