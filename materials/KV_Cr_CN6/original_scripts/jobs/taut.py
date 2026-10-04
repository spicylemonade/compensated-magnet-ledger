"""Track L PBA "taut": electronic competitors of the d3/d3 LCM in KV[M(CN)6] (M = Cr, Mo), 15-atom primitive cell.
usage: python taut.py cfg.json        (results/<id>.json, rewritten after every step; every step resumes)

cfg: {id, compound, MN ("V"), MC ("Cr"/"Mo"), cell {lattice, species, frac} (exact F-43m stage-1 projection),
      U {el: value} (value 0 -> no Hubbard term on that element), ecut (110), kpts [4,4,4], npool,
      starts [names from LIB], relax: {"classes": [...], "n": 1, "maxsec": s} | null,
      hse: {"ecut": 90, "nq": [2,2,2], "kpts": [4,4,4], "E_LCM": <eV, existing HSE LCM at the cubic cell>} | null}

Physics / protocol
  * PseudoDojo SR NC, PBE+U ortho-atomic, Gaussian 0.005 Ry, conv_thr 1e-7*nat, mixing 0.3 (= Gate-0 stage 1).
  * nosym = .true. for every PBE+U run (the LCM reference included, same k set): the competitors are orbitally
    polarised (V(III) t2g^2, M(II) t2g^4 / t2g^3 eg^1, M(III) S=1/2 t2g^3), and QE symmetrises ns and rho with the
    crystal group, so in the cubic cell (Td site symmetry) any t2g orbital polarisation would be averaged away.
  * Constrained starts: site starting_magnetization + starting_ns_eigenvalue + (optional) tot_magnetization.
    QE 7.5 applies starting_ns_eigenvalue to the STARTING ns, which is diagonal and degenerate per spin, so index m is
    the real harmonic m in QE order (1 z2, 2 xz, 3 yz, 4 x2-y2, 5 xy) in the cubic axes of the cell (smoke test
    lcm/pba/taut/smoke1). Each start fixes all 5 orbitals per spin on each Hubbard site (1 = occupied, 0 = empty).
    Default pattern: the D2d "complementary" one (V hole in xy, M(II) extra electron in xy; S=1/2 d3 = up xz,yz + down xy).
  * relax: for the lowest converged start whose final state is not the LCM, in each class of cfg.relax.classes:
    ionic relaxation at the fixed stage-1 cell (nosym, same start). Then the relaxed cell is symmetrised exactly
    (spglib, symprec 1e-3 A) and, if cfg.hse: a PBE+U prep SCF at the HSE cutoff with symmetry ON (the distorted cell
    keeps the orbital order allowed) and the same start, then HSE06 (no U, nq 2, ecutfock 2*ecut, gygi-baldereschi,
    x_gamma_extrapolation, conv_thr 1e-8*nat = lcm_deep/pba_fmfix settings) started from the prep wavefunctions and
    density (startingwfc/startingpot = 'file'; QE still runs a semilocal SCF first, then the EXX loops). If the prep
    does not reproduce the relaxed state, it is retried with nosym and the HSE then runs with nosym too.
"""
import json, os, re, sys, time, shutil, traceback
import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, save_json, SCRATCH

NC = int(os.environ.get("NCORES", 64))
cfg = json.load(open(sys.argv[1]))
CID = cfg["id"]
NPOOL = int(cfg.get("npool", 16))
os.makedirs("results", exist_ok=True)
RES = f"results/{CID}.json"
MN, MC = cfg["MN"], cfg["MC"]
ORB = {"V": "3d", "Cr": "3d", "Mo": "4d"}
Z2, XZ, YZ, X2, XY = 1, 2, 3, 4, 5
T2G = (XZ, YZ, XY)


LIBVER = 2  # v1 (bbfc1d36) mapped starting_ns_eigenvalue indices to orbitals; wrong in QE 7.5 (see docstring)
RAMP_U = 6.0  # eV, pre-stage U for the *_ramp starts


def site(sign, maj=(0, 0), mino=(0, 0)):
    """sign of the moment; (n_t2g, n_eg) occupied in the majority / minority channel -> {ispin: (n_t2g, n_eg)}"""
    smaj, smin = (1, 2) if sign > 0 else (2, 1)
    return {smaj: tuple(maj), smin: tuple(mino)}


# name: (class, (mom_N, mom_C), tot_mag, ns spec N-bound site, ns spec C-bound site[, "ramp"])
LIB = {
    "LCM":        ("LCM", (-3, 3), None, None, None),
    "FM6":        ("FM", (3, 3), 6.0, site(+1, (3, 0)), site(+1, (3, 0))),
    # charge transfer V(II)M(III) -> V(III) t2g^2 (S=1) / M(II) d4
    "CT_AF":      ("CT", (-2, 2), None, site(-1, (2, 0)), site(+1, (3, 0), (1, 0))),
    "CT_AF_v":    ("CT", (-2, 2), None, site(-1, (2, 0)), None),
    "CT_AF_ramp": ("CT", (-2, 2), None, site(-1, (2, 0)), site(+1, (3, 0), (1, 0)), "ramp"),
    "CT_AF_free": ("CT", (-2, 2), None, None, None),
    "CT_FM":      ("CT", (2, 2), 4.0, site(+1, (2, 0)), site(+1, (3, 0), (1, 0))),
    "CT_HS_AF":   ("CT", (-2, 4), 2.0, site(-1, (2, 0)), site(+1, (3, 1))),
    "CT_HS_FM":   ("CT", (2, 4), 6.0, site(+1, (2, 0)), site(+1, (3, 1))),
    # no charge transfer, other spin states: M(III) / V(II) low spin S=1/2 (t2g up^2 down^1)
    "MLS_AF":     ("LS", (-3, 1), -2.0, site(-1, (3, 0)), site(+1, (2, 0), (1, 0))),
    "MLS_FM":     ("LS", (3, 1), 4.0, site(+1, (3, 0)), site(+1, (2, 0), (1, 0))),
    "VLS_AF":     ("LS", (-1, 3), 2.0, site(-1, (2, 0), (1, 0)), site(+1, (3, 0))),
    "LSLS_AF":    ("LS", (-1, 1), 0.0, site(-1, (2, 0), (1, 0)), site(+1, (2, 0), (1, 0))),
}
NS_FREE = {k for k, v in LIB.items() if v[3] is None and v[4] is None}


cell = cfg["cell"]
AT0 = Atoms(cell["species"], cell=cell["lattice"], scaled_positions=cell["frac"], pbc=True)
IN = cell["species"].index(MN)
IC = cell["species"].index(MC)
U = {el: (ORB[el], float(u)) for el, u in cfg["U"].items() if float(u) > 0}

out = json.load(open(RES)) if os.path.exists(RES) else {}
out.update(id=CID, compound=cfg["compound"], U=cfg["U"], ecut=cfg["ecut"], kpts=cfg["kpts"], idx={"N": IN + 1, "C": IC + 1})
out.setdefault("runs", {})
out.setdefault("relax", {})
out.setdefault("hse", {})
# relax / HSE records written by taut.py v1 (killed 10:12) are stale: keep only current-version records
out["relax"] = {k: v for k, v in out["relax"].items() if v.get("libver") == LIBVER and v.get("run")}
out["hse"] = {k: v for k, v in out["hse"].items() if v.get("libver") == LIBVER}


def save():
    save_json(out, RES)


# ---------------------------------------------------------------------------------------------------- parsing
def parse_hub(txt):
    """last HUBBARD OCCUPATIONS block -> {atom(1-based): {"tr": [up, dn, tot], "eig": {s: [5]}, "diag": {s: [5]}}}"""
    i = txt.rfind("=================== HUBBARD OCCUPATIONS")
    if i < 0:
        return {}
    blk = txt[i:]
    j = blk.find("Number of occupied Hubbard levels")
    blk = blk[:j] if j > 0 else blk[:20000]
    res = {}
    for part in re.split(r"-{10,} ATOM\s+", blk)[1:]:
        try:
            na = int(part.split()[0])
            tr = re.search(r"Tr\[ns\(\s*\d+\)\] \(up, down, total\) =\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", part)
            d = {"tr": [float(x) for x in tr.groups()] if tr else None, "eig": {}, "diag": {}}
            for sp in re.split(r"SPIN\s+", part)[1:]:
                s = int(sp.split()[0])
                lines = sp.splitlines()
                for k, ln in enumerate(lines):
                    if "eigenvalues:" in ln:
                        d["eig"][s] = [float(x) for x in lines[k + 1].split()]
                    if "occupation matrix ns" in ln:
                        mat = [[float(x) for x in lines[k + 1 + r].split()] for r in range(5)]
                        d["diag"][s] = [round(mat[r][r], 3) for r in range(5)]
            res[na] = d
        except Exception:
            continue
    return res


def glob_windows(b):
    E = b["eig"]; occ = b["occ"]
    occd = occ > 0.5
    vb = np.where(occd, E, -np.inf); cb = np.where(~occd, E, np.inf)
    if E.shape[1] != 2:
        return {"gap": float(cb.min() - vb.max())}
    tv = vb.max(axis=(0, 2)); bc = cb.min(axis=(0, 2))
    sv = int(np.argmax(tv)); sc = int(np.argmin(bc))
    # fractional occupations (metal / two-Fermi-level smearing): count states with 0.05 < occ < 0.95
    nfrac = int(((occ > 0.05) & (occ < 0.95)).sum())
    return {"gap": float(bc.min() - tv.max()), "vbm_spin": sv + 1, "cbm_spin": sc + 1, "win_VB": float(tv[sv] - tv[1 - sv]),
            "win_CB": float(bc[1 - sc] - bc[sc]), "gap_up": float(bc[0] - tv[0]), "gap_dn": float(bc[1] - tv[1]), "n_frac_occ": nfrac}


def bonds(at):
    """M-ligand bond lengths along +-x,y,z for both metals (min image)"""
    pos = at.get_positions(); cellv = at.cell.array; inv = np.linalg.inv(cellv)
    sym = at.get_chemical_symbols()
    res = {}
    for i, lig in ((IN, "N"), (IC, "C")):
        bl = []
        for j, s in enumerate(sym):
            if s != lig:
                continue
            d = pos[j] - pos[i]; f = d @ inv; f -= np.round(f); d = f @ cellv
            r = np.linalg.norm(d)
            if r < 2.6:
                ax = "xyz"[int(np.argmax(np.abs(d)))]
                bl.append((ax + ("+" if d[np.argmax(np.abs(d))] > 0 else "-"), round(float(r), 4)))
        res[sym[i]] = dict(sorted(bl))
    return res


def summarize(p, sname):
    hub = p.get("hub") or {}
    d = {}
    for role, ia in (("N", IN + 1), ("C", IC + 1)):
        h = hub.get(ia) or hub.get(str(ia))
        if h and h.get("tr"):
            d[role] = {"m_d": round(h["tr"][0] - h["tr"][1], 3), "N_d": round(h["tr"][2], 3), "diag": h.get("diag")}
        sm = p.get("site_moments") or []
        if len(sm) > ia - 1:
            d.setdefault(role, {})["m_sph"] = sm[ia - 1]
    return d


def classify(s, ref):
    """label of the converged state from d occupations relative to the LCM at the same settings"""
    try:
        sN, sC = s["sites"]["N"], s["sites"]["C"]
        rN, rC = ref["sites"]["N"], ref["sites"]["C"]
        dNN = sN["N_d"] - rN["N_d"] if "N_d" in sN else None
        dNC = sC["N_d"] - rC["N_d"] if ("N_d" in sC and "N_d" in rC) else None
        mN = sN.get("m_d", sN.get("m_sph")); mC = sC.get("m_d", sC.get("m_sph"))
        rmN = rN.get("m_d", rN.get("m_sph")); rmC = rC.get("m_d", rC.get("m_sph"))
        lab = []
        ct = dNN is not None and dNN < -0.15
        lab.append("CT" if ct else "noCT")
        lab.append("N:%+.2f" % mN); lab.append("C:%+.2f" % mC)
        same = (abs(mN - rmN) < 0.25 and abs(mC - rmC) < 0.25 and (dNN is None or abs(dNN) < 0.1))
        return {"ct": bool(ct), "dN_d_N": dNN, "dN_d_C": dNC, "lcm_like": bool(same), "label": " ".join(lab),
                "aligned": bool(np.sign(mN) == np.sign(mC))}
    except Exception as e:
        return {"err": repr(e)}


# ---------------------------------------------------------------------------------------------------- runners
def lib(sname):
    e = LIB[sname]
    return e[0], e[1], e[2], e[3], e[4], (len(e) > 5 and e[5] == "ramp")


def species_index(at, mom, hubbard):
    """QE species index of the two metals (write_pw groups atoms by (element, moment) in order of first appearance)"""
    info = write_pw("probe.in", at, moments=mom, hubbard=hubbard)
    labs = list(info["species"].keys())
    return {MN: labs.index(info["labels"][IN]) + 1, MC: labs.index(info["labels"][IC]) + 1}


def iter1_ns(txt):
    """ns after SCF iteration 1 (2nd HUBBARD OCCUPATIONS block): {atom: {spin: (eigenvalues[5], eigvec columns [5x5])}}"""
    blocks = txt.split("=================== HUBBARD OCCUPATIONS")
    if len(blocks) < 3:
        return {}
    blk = blocks[2]
    j = blk.find("Number of occupied Hubbard levels")
    blk = blk[:j] if j > 0 else blk[:20000]
    res = {}
    for part in re.split(r"-{10,} ATOM\s+", blk)[1:]:
        na = int(part.split()[0])
        res[na] = {}
        for sp in re.split(r"SPIN\s+", part)[1:]:
            if not sp.split() or not sp.split()[0].isdigit():
                continue
            s = int(sp.split()[0])
            L = sp.splitlines()
            ev = vec = None
            for k, ln in enumerate(L):
                if "eigenvalues:" in ln:
                    ev = [float(x) for x in L[k + 1].split()]
                if "eigenvectors" in ln:
                    vec = np.array([[float(x) for x in L[k + 1 + r].split()] for r in range(5)])
            if ev is not None and vec is not None:
                res[na][s] = (ev, vec)
    return res


def ns_map(spec, ev, vec):
    """(n_t2g, n_eg) occupied -> {m: value} for starting_ns_eigenvalue. Columns are classified t2g / eg by their weight on
    z2 + x2-y2 (QE real-harmonic rows 1 and 4); within each group the most occupied eigenvectors get 1.
    t2g: all three set (n -> 1, rest -> 0); eg: only the n_eg most occupied set to 1 (the rest keep their covalent value)."""
    nt, ne = spec
    egw = vec[0] ** 2 + vec[3] ** 2
    cols = list(range(5))
    t2g = sorted([c for c in cols if egw[c] < 0.5], key=lambda c: -ev[c])
    eg = sorted([c for c in cols if egw[c] >= 0.5], key=lambda c: -ev[c])
    out_ = {}
    if len(t2g) == 3:
        for i, c in enumerate(t2g):
            out_[c + 1] = 1.0 if i < nt else 0.0
    for i, c in enumerate(eg[:ne]):
        out_[c + 1] = 1.0
    return out_, {"t2g_cols": [c + 1 for c in t2g], "eg_cols": [c + 1 for c in eg], "ev": ev}


def scf_run(key, sname, at, nosym=True, ecut=None, kpts=None, hubbard=U, extra_sys=None, extra_el=None, calc="scf",
            npool=None, conv=None, timeout=4 * 3600, keep=False, control=None, ramp=None):
    """One constrained run. Starts with an ns spec use two passes: pass 1 = 1 SCF iteration with the same moments /
    tot_magnetization to read the iteration-1 ns eigenvectors (QE 7.5 applies starting_ns_eigenvalue to those, in
    ascending-eigenvalue order); pass 2 = the real run with starting_ns_eigenvalue mapped onto t2g / eg eigenvectors.
    ramp (start flag or arg): pass 2 is an SCF at U = RAMP_U on every Hubbard element, then the real run (scf or relax)
    restarts from its density, ns (occup.txt) and wavefunctions at the target U, with no ns kick."""
    ecut = ecut or cfg["ecut"]; kpts = tuple(kpts or cfg["kpts"])
    cls, (mN, mC), tot, nsN, nsC, rflag = lib(sname)
    ramp = rflag if ramp is None else ramp
    mom = [0.0] * len(at); mom[IN] = float(mN); mom[IC] = float(mC)
    base_sys = {"nosym": True} if nosym else {}
    if tot is not None:
        base_sys["tot_magnetization"] = float(tot)
    base_sys.update(extra_sys or {})
    el0 = {"mixing_beta": 0.3, "electron_maxstep": 250}
    el0.update(extra_el or {})
    pref = key
    info = {}
    hub_kick = ({el: (o, RAMP_U) for el, (o, u) in hubbard.items()} if (ramp and hubbard) else hubbard)

    def write_run(tag, calc_, hub, sysx, el, ctrl=None, conv_=None):
        write_pw(f"{pref}{tag}.in", at, calculation=calc_, pset="dojo_sr", ecutwfc=ecut, kpts=kpts, moments=mom, smearing="gaussian",
                 degauss=0.005, hubbard=hub, prefix=pref, conv_thr=conv_ or conv or 1e-7 * len(at), tstress=False,
                 tprnfor=(calc_ == "scf"), system=sysx, electrons=el, control=ctrl)
        return run_pw(f"{pref}{tag}.in", ncores=NC, npool=npool or NPOOL, timeout=timeout)

    kick = {}
    if hub_kick and (nsN or nsC):
        r1 = write_run("_p1", "scf", hub_kick, dict(base_sys), dict(el0, electron_maxstep=1, scf_must_converge=False))
        t1 = open(r1["out"], errors="replace").read()
        it1 = iter1_ns(t1)
        ity = species_index(at, mom, hub_kick)
        for el_, ns, ia in ((MN, nsN, IN + 1), (MC, nsC, IC + 1)):
            if ns is None or el_ not in hub_kick or ia not in it1:
                continue
            for s, spec in ns.items():
                if s not in it1[ia]:
                    continue
                mp, inf = ns_map(spec, *it1[ia][s])
                info[f"{el_}_s{s}"] = dict(inf, set=mp, spec=spec)
                for m, v in mp.items():
                    kick[f"starting_ns_eigenvalue({m},{s},{ity[el_]})"] = v
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
    pre = None
    if ramp:
        rr = write_run("_ramp", "scf", hub_kick, dict(base_sys, **kick), dict(el0))
        pre = parse_pw(rr["out"])
        pre_hub = parse_hub(open(rr["out"], errors="replace").read())
        pre = {"energy_eV": pre.get("energy_eV"), "converged_scf": pre.get("converged_scf"), "total_mag": pre.get("total_mag"),
               "sites": summarize({"hub": pre_hub, "site_moments": pre.get("site_moments")}, sname), "sec": rr["seconds"]}
        r = write_run("", calc, hubbard, dict(base_sys), dict(el0, startingwfc="file", startingpot="file"), control)
    else:
        r = write_run("", calc, hubbard, dict(base_sys, **kick), dict(el0), control)
    p = parse_pw(r["out"]); p["sec"] = r["seconds"]; p["rc"] = r["rc"]
    txt = open(r["out"], errors="replace").read()
    p["hub"] = parse_hub(txt)
    p["nks"] = (re.findall(r"number of k points=\s*(\d+)", txt) or [None])[0]
    p["nsym"] = (re.findall(r"(\d+) Sym\. Ops\.", txt) or (["1 (no symmetry)"] if "No symmetry found" in txt else [None]))[0]
    p["ns_modified"] = "Modifying starting ns" in txt
    try:
        b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
        p["glob"] = glob_windows(b)
    except Exception as e:
        p["glob_err"] = repr(e)
    keepk = ("energy_eV", "converged_scf", "total_mag", "abs_mag", "site_moments", "glob", "sec", "rc", "errors", "nks", "nsym",
             "ns_modified", "total_force", "n_scf_cycles", "bfgs_converged", "hub")
    q = {k: p.get(k) for k in keepk}
    q["start"] = sname
    q["libver"] = LIBVER
    q["kick"] = info
    if pre:
        q["ramp_stage"] = pre
    q["sites"] = summarize(p, sname)
    if not keep:
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
    return q, r["out"]


def run_vertical():
    for sname in cfg["starts"]:
        prev = out["runs"].get(sname)
        if prev and prev.get("energy_eV") is not None and (prev.get("converged_scf") or prev.get("retried")) and \
                (prev.get("libver") == LIBVER or sname in NS_FREE or sname == "FM6"):
            continue
        print(f"[{CID}] vertical {sname}", flush=True)
        t0 = time.time()
        q, _ = scf_run(f"v_{sname}", sname, AT0, timeout=3 * 3600)
        if not q.get("converged_scf"):
            q2, _ = scf_run(f"v_{sname}", sname, AT0, timeout=3 * 3600, extra_el={"mixing_beta": 0.1, "electron_maxstep": 400})
            q2["retried"] = True; q2["first_attempt"] = {k: q.get(k) for k in ("energy_eV", "total_mag", "n_scf_cycles")}
            q = q2
        q["cls_start"] = lib(sname)[0]
        out["runs"][sname] = q
        save()
        print(f"[{CID}] {sname}: E={q.get('energy_eV')} M={q.get('total_mag')} conv={q.get('converged_scf')} "
              f"gap={(q.get('glob') or {}).get('gap')} sites={q['sites']} {time.time() - t0:.0f}s", flush=True)
    ref = out["runs"].get("LCM")
    if ref and ref.get("energy_eV") is not None:
        for k, q in out["runs"].items():
            if q.get("energy_eV") is not None:
                q["dE_eV"] = q["energy_eV"] - ref["energy_eV"]
                q["state"] = classify(q, ref)
        save()


def relaxed_atoms(outfile):
    from ase.io import read
    return read(outfile, format="espresso-out", index=-1)


def symmetrize(at, prec=1e-3):
    import spglib
    from ase.spacegroup.symmetrize import refine_symmetry
    cellt = (at.cell.array, at.get_scaled_positions(), at.numbers)
    sg = {str(p): spglib.get_spacegroup(cellt, symprec=p) for p in (1e-4, 1e-3, 1e-2, 3e-2)}
    a2 = at.copy()
    refine_symmetry(a2, symprec=prec)
    dev = float(np.abs(a2.get_positions() - at.get_positions()).max())
    return a2, sg, dev


def run_relax_chain():
    rc = cfg.get("relax")
    if not rc:
        return
    missing = [k for k in cfg["starts"] if (out["runs"].get(k) or {}).get("energy_eV") is None
               or not ((out["runs"][k].get("libver") == LIBVER) or k in NS_FREE or k == "FM6")]
    if missing:
        print(f"[{CID}] vertical scan incomplete ({missing}); no relax", flush=True)
        return
    ref = out["runs"].get("LCM")
    for cls in rc["classes"]:
        if out["relax"].get(cls, {}).get("done"):
            continue
        cands = [(q["energy_eV"], k) for k, q in out["runs"].items()
                 if q.get("energy_eV") is not None and q.get("cls_start") == cls and q.get("converged_scf")
                 and not (q.get("state") or {}).get("lcm_like")]
        if cls == "CT":  # the CT class must actually be charge transferred
            cands = [(e, k) for e, k in cands if (out["runs"][k].get("state") or {}).get("ct")]
        cands.sort()
        rec = {"candidates": [[k, e - ref["energy_eV"]] for e, k in cands], "libver": LIBVER}
        if not cands:
            rec.update(done=True, note="no converged non-LCM state in this class")
            out["relax"][cls] = rec; save(); continue
        sname = cands[0][1]
        rec["target"] = sname
        print(f"[{CID}] relax {cls}: {sname}", flush=True)
        t0 = time.time()
        q, of = scf_run(f"r_{sname}", sname, AT0, calc="relax", keep=True, timeout=int(rc.get("maxsec", 9000)) + 3600,
                        control={"max_seconds": float(rc.get("maxsec", 9000)), "nstep": 80})
        try:
            at = relaxed_atoms(of)
            q["relaxed"] = {"lattice": at.cell.array.tolist(), "species": at.get_chemical_symbols(), "frac": at.get_scaled_positions().tolist()}
            q["max_disp_A"] = float(np.abs(at.get_positions() - AT0.get_positions()).max())
            q["bonds"] = bonds(at)
            q["bonds0"] = bonds(AT0)
        except Exception as e:
            q["relax_parse_err"] = repr(e)
        q["dE_eV"] = q["energy_eV"] - ref["energy_eV"] if q.get("energy_eV") is not None else None
        q["E_relax_gain_eV"] = (q["energy_eV"] - out["runs"][sname]["energy_eV"]) if q.get("energy_eV") is not None else None
        q["state"] = classify(q, ref)
        q["sec"] = time.time() - t0
        rec["run"] = q
        out["relax"][cls] = rec
        save()
        print(f"[{CID}] relaxed {sname}: dE={q.get('dE_eV')} gain={q.get('E_relax_gain_eV')} state={q.get('state')} bonds={q.get('bonds')}", flush=True)
        if cfg.get("hse") and q.get("relaxed"):
            try:
                run_hse(cls, sname, q)
            except Exception:
                traceback.print_exc()
        rec["done"] = True
        save()
        shutil.rmtree(SCRATCH / f"r_{sname}", ignore_errors=True)


def run_hse(cls, sname, rq):
    h = cfg["hse"]
    if out["hse"].get(cls, {}).get("done"):
        return
    rl = rq["relaxed"]
    at = Atoms(rl["species"], cell=rl["lattice"], scaled_positions=rl["frac"], pbc=True)
    at_s, sg, dev = symmetrize(at)
    rec = {"target": sname, "spacegroups": sg, "sym_dev_A": dev, "libver": LIBVER,
           "cell": {"lattice": at_s.cell.array.tolist(), "species": at_s.get_chemical_symbols(), "frac": at_s.get_scaled_positions().tolist()}}
    out["hse"][cls] = rec; save()
    want = rq["sites"]
    nosym = False
    prep = None
    for attempt, ns in enumerate((False, True)):
        pref = f"h_{sname}"
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
        q, _ = scf_run(pref, sname, at_s, nosym=ns, ecut=h["ecut"], kpts=h["kpts"], conv=1e-8 * len(at_s), keep=True,
                       npool=max(1, NC // 8), timeout=6 * 3600)
        ok = q.get("converged_scf") and all(abs(q["sites"][r].get("m_d", 0) - want[r].get("m_d", 0)) < 0.3 for r in ("N", "C")
                                            if "m_d" in want[r])
        rec[f"prep{attempt}"] = q
        save()
        print(f"[{CID}] HSE prep ({'nosym' if ns else 'sym'}) {sname}: ok={ok} nsym={q.get('nsym')} nks={q.get('nks')} sites={q['sites']}", flush=True)
        if ok:
            nosym, prep = ns, q
            break
    if prep is None:
        rec.update(done=True, note="PBE+U prep at the HSE cutoff did not reproduce the relaxed state; no HSE")
        save(); return
    cls_, (mN, mC), tot, _, _, _ = lib(sname)
    mom = [0.0] * len(at_s); mom[IN] = float(mN); mom[IC] = float(mC)
    nq = h["nq"]
    sysx = {"input_dft": "hse", "nqx1": nq[0], "nqx2": nq[1], "nqx3": nq[2], "ecutfock": 2 * h["ecut"],
            "exxdiv_treatment": "gygi-baldereschi", "x_gamma_extrapolation": True}
    if nosym:
        sysx["nosym"] = True
    if tot is not None:
        sysx["tot_magnetization"] = float(tot)
    pref = f"h_{sname}"
    write_pw(f"{pref}_hse.in", at_s, pset="dojo_sr", ecutwfc=h["ecut"], kpts=tuple(h["kpts"]), moments=mom, smearing="gaussian",
             degauss=0.005, hubbard=None, prefix=pref, conv_thr=1e-8 * len(at_s), tstress=False, tprnfor=False, system=sysx,
             electrons={"mixing_beta": 0.3, "electron_maxstep": 200, "startingwfc": "file", "startingpot": "file"})
    t0 = time.time()
    r = run_pw(f"{pref}_hse.in", ncores=NC, npool=max(1, NC // 8), timeout=22 * 3600)
    p = parse_pw(r["out"])
    txt = open(r["out"], errors="replace").read()
    hq = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag", "abs_mag", "site_moments", "errors")}
    hq["sec"] = time.time() - t0
    hq["nsym"] = (re.findall(r"(\d+) Sym\. Ops\.", txt) or [None])[0]
    hq["nks"] = (re.findall(r"number of k points=\s*(\d+)", txt) or [None])[0]
    hq["nosym"] = nosym
    try:
        b = parse_xml_bands(str(SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml"))
        hq["glob"] = glob_windows(b)
    except Exception as e:
        hq["glob_err"] = repr(e)
    sm = hq.get("site_moments") or []
    hq["m_sph"] = {"N": sm[IN] if len(sm) > IN else None, "C": sm[IC] if len(sm) > IC else None}
    hq["m_sph_prep"] = {r_: prep["sites"][r_].get("m_sph") for r_ in ("N", "C")}
    if h.get("E_LCM") is not None and hq.get("energy_eV") is not None:
        hq["dE_eV"] = hq["energy_eV"] - h["E_LCM"]
    rec["run"] = hq
    rec["done"] = True
    save()
    print(f"[{CID}] HSE {sname}: {hq}", flush=True)
    shutil.rmtree(SCRATCH / pref, ignore_errors=True)


try:
    run_vertical()
except Exception:
    traceback.print_exc()
try:
    run_relax_chain()
except Exception:
    traceback.print_exc()
out["t_end"] = time.time()
save()
print(f"[{CID}] done", flush=True)
