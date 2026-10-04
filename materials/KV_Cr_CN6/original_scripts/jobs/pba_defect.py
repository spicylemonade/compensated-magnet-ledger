"""Track L PBA defect Gate 1 job: one 59-62-atom defect (or pristine 4 f.u.) cell of KV[Cr(CN)6] / KV[Mo(CN)6].

usage: python pba_defect.py <cfg.json>       (cfg written by tracks/lcm/pba/defects/build_defects.py)

Steps (every step resumes after container preemption; state in results/<id>.json):
  1. relax   PBE+U ions-only BFGS at the fixed stage-1 cell, LCM start (site spins x |m0|), k 2x2x2 shifted (1 1 1).
             An interrupted relax restarts from the last complete ATOMIC_POSITIONS block of its .out (renamed _prev<n>).
  2. scf     fresh SCF at the relaxed positions (same k) -> energy, M, |M|, sphere moments, ortho-atomic d occupations.
  3. nscf    4x4x4 Gamma-centred grid (contains the folded Gamma, X, L, W of the fcc primitive BZ), nbnd = Nel/2 + NB_EXTRA.
             A failed Davidson nscf is retried once with CG.
  4. projwfc Lowdin projections (lsym false) on the nscf -> per-atom s/p/d weights of every state.
  5. analysis (in the container; raw files kept in raw/<id>/ so ANALYZE_ONLY=1 can redo it):
       occupation-based edges per spin, windows (lane definition), Fermi level;
       band table near the gap (energy range, mean occupation, group weights, defect-region weight, top atoms);
       deep potential references for cell-to-cell alignment: k-weighted mean energy of the CN 3sigma manifold
       (C+N weight > 0.6, E < VBM - 14 eV) projected on C/N atoms farther than R_FAR from every defect centre
       ("E_ref_CN3s_far"), the same on all C/N, and the K 3p semicore level on far K.
     -> results/<id>.json (summary), results/<id>_bands.npz (E, occ, wk, k, per-atom s/p/d weights).
Protocol = PBA Gate 0 stage 1: PseudoDojo SR NC, 110 Ry (ecutrho 440), ortho-atomic U from cfg, gaussian 0.005 Ry,
conv_thr 1e-7*nat, mixing_beta 0.3; BFGS forc_conv_thr 1e-3 Ry/bohr, etot_conv_thr 1e-5*nat Ry.
NPOOL: largest divisor of NCORES <= min(NPOOL_MAX, 2 * irreducible k) (spglib count with the QE species = element +
starting moment, time reversal on; QE only warns if it finds fewer k).
env: ECUT(110) KRELAX("2 2 2 1 1 1") KNSCF("4 4 4 0 0 0") NB_EXTRA(40) RELAX_MAXSEC(14400) SCF_MAXSEC(7200)
     NPOOL_MAX(NCORES//4) R_FAR(4.5 A) ANALYZE_ONLY(0) SKIP_RELAX(0)"""
import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
import traceback
import xml.etree.ElementTree as ET

import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, parse_xml_bands, save_json, zvalence, SCRATCH

NC = int(os.environ.get("NCORES", 32))
ECUT = float(os.environ.get("ECUT", 110))
KRELAX = [int(x) for x in os.environ.get("KRELAX", "2 2 2 1 1 1").split()]
KNSCF = [int(x) for x in os.environ.get("KNSCF", "4 4 4 0 0 0").split()]
NB_EXTRA = int(os.environ.get("NB_EXTRA", 40))
RELAX_MAXSEC = float(os.environ.get("RELAX_MAXSEC", 14400))
SCF_MAXSEC = float(os.environ.get("SCF_MAXSEC", 7200))
NPOOL_MAX = int(os.environ.get("NPOOL_MAX", max(1, NC // 4)))
R_FAR = float(os.environ.get("R_FAR", 4.5))
ANALYZE_ONLY = os.environ.get("ANALYZE_ONLY", "0") == "1"
SKIP_RELAX = os.environ.get("SKIP_RELAX", "0") == "1"
RY_EV = 13.605693122994
os.makedirs("results", exist_ok=True)
os.makedirs("raw", exist_ok=True)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------------------------------------ helpers
def nk_irr(at, mom, mesh, shift):
    try:
        import spglib
        keys = {}
        types = [keys.setdefault((s, round(float(m), 2)), len(keys) + 1) for s, m in zip(at.get_chemical_symbols(), mom)]
        mp, _ = spglib.get_ir_reciprocal_mesh(mesh, (at.cell.array, at.get_scaled_positions(), types), is_shift=shift,
                                              is_time_reversal=True, symprec=1e-3)
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


def last_positions(txt, nat):
    """last complete 'ATOMIC_POSITIONS (crystal)' block -> (frac[nat,3], n_blocks) or None"""
    lines = txt.splitlines()
    best, n = None, 0
    for i, ln in enumerate(lines):
        if not ln.startswith("ATOMIC_POSITIONS") or "crystal" not in ln:
            continue
        try:
            fr = np.array([[float(x) for x in lines[i + 1 + k].split()[1:4]] for k in range(nat)])
            if fr.shape == (nat, 3):
                best, n = fr, n + 1
        except (IndexError, ValueError):
            continue
    return None if best is None else (best, n)


def parse_hub(txt):
    out = {}
    for i, u, d, t in re.findall(r"Tr\[ns\(\s*(\d+)\)\]\s*\(up, down, total\)\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", txt):
        out[int(i)] = [float(u), float(d), float(t)]
    return out


def write(path, at, calc, mom, U, kp, pref, extra_sys, nbnd=None, electrons=None, control=None):
    sysx = dict(extra_sys or {})
    if nbnd:
        sysx["nbnd"] = nbnd
    el = {"mixing_beta": 0.3}
    el.update(electrons or {})
    return write_pw(path, at, calculation=calc, pset="dojo_sr", ecutwfc=ECUT, kpts=kp[:3], kshift=tuple(kp[3:]),
                    moments=mom, smearing="gaussian", degauss=0.005, hubbard=U, prefix=pref, conv_thr=1e-7 * len(at),
                    tstress=False, tprnfor=(calc == "relax"), electrons=el, control=control, system=sysx or None)


def save(out):
    save_json(out, f"results/{out['id']}.json")


# ------------------------------------------------------------------------------------------------ QE steps
def do_relax(c, out, at, mom, U):
    st = out.setdefault("relax", {})
    if st.get("done"):
        fr = np.array(st["frac"])
        return Atoms(at.get_chemical_symbols(), cell=at.cell, scaled_positions=fr, pbc=True)
    pref = f"{c['id']}_relax"
    if os.path.exists(f"{pref}.out"):
        g = last_positions(open(f"{pref}.out", errors="replace").read(), len(at))
        nprev = len(glob.glob(f"{pref}_prev*.out"))
        shutil.move(f"{pref}.out", f"{pref}_prev{nprev}.out")
        if g is not None:
            at = Atoms(at.get_chemical_symbols(), cell=at.cell, scaled_positions=g[0], pbc=True)
            st.setdefault("restarts", []).append({"prev": f"{pref}_prev{nprev}.out", "n_blocks": g[1]})
            log("relax restart from", f"{pref}_prev{nprev}.out", "blocks", g[1])
    nk = nk_irr(at, mom, KRELAX[:3], KRELAX[3:])
    npool = pick_npool(nk)
    write(f"{pref}.in", at, "relax", mom, U, KRELAX, pref, c.get("system_extra"),
          control={"nstep": 120, "max_seconds": RELAX_MAXSEC})
    log("relax start nk_irr", nk, "npool", npool)
    r = run_pw(f"{pref}.in", ncores=NC, npool=npool, timeout=RELAX_MAXSEC + 3600)
    txt = open(f"{pref}.out", errors="replace").read()
    p = parse_pw(f"{pref}.out")
    g = last_positions(txt, len(at))
    nsteps = len(re.findall(r"^!\s+total energy", txt, re.M))
    if "JOB DONE" not in txt:
        # crash: keep the .out so a resubmission restarts from its last positions; never mark the relax done
        st.update(rc=r["rc"], errors=p.get("errors"), n_scf=nsteps, sec=st.get("sec", 0) + r["seconds"])
        save(out)
        raise RuntimeError(f"relax crashed rc={r['rc']} errors={p.get('errors')}")
    fr = g[0] if g is not None else at.get_scaled_positions()
    st.update(done=True, frac=fr.tolist(), sec=st.get("sec", 0) + r["seconds"], rc=r["rc"], nk_irr=nk, npool=npool,
              bfgs_converged=p.get("bfgs_converged"), n_scf=nsteps, energy_eV=p.get("energy_eV"),
              total_force=p.get("total_force"), total_mag=p.get("total_mag"), abs_mag=p.get("abs_mag"),
              errors=p.get("errors"), max_seconds_hit=("Maximum CPU time exceeded" in txt),
              stopped_unconverged=not p.get("bfgs_converged"))
    a0 = c["_at0"]
    d = fr - a0.get_scaled_positions(); d -= np.round(d)
    disp = np.linalg.norm(d @ at.cell.array, axis=1)
    st["max_disp_A"] = float(disp.max())
    st["top_disp"] = [[int(i), a0.get_chemical_symbols()[i], round(float(disp[i]), 4)] for i in np.argsort(-disp)[:8]]
    shutil.rmtree(SCRATCH / pref, ignore_errors=True)
    save(out)
    log("relax done", st["bfgs_converged"], st["energy_eV"], st["total_mag"], "nscf", nsteps, "maxdisp", st["max_disp_A"])
    return Atoms(at.get_chemical_symbols(), cell=at.cell, scaled_positions=fr, pbc=True)


def do_scf(c, out, at, mom, U, pref):
    nk = nk_irr(at, mom, KRELAX[:3], KRELAX[3:])
    npool = pick_npool(nk)
    write(f"{pref}.in", at, "scf", mom, U, KRELAX, pref, c.get("system_extra"), control={"max_seconds": SCF_MAXSEC})
    r = run_pw(f"{pref}.in", ncores=NC, npool=npool, timeout=SCF_MAXSEC + 3600)
    txt = open(f"{pref}.out", errors="replace").read()
    p = parse_pw(f"{pref}.out")
    hub = parse_hub(txt)
    s = {k: p.get(k) for k in ("energy_eV", "total_mag", "abs_mag", "converged_scf", "n_scf_cycles", "site_moments",
                                "fermi_eV", "homo_lumo", "errors", "wall")}
    s["n_iter"] = int(re.findall(r"convergence has been achieved in\s+(\d+)", txt)[-1]) if "convergence has been achieved" in txt else None
    s.update(sec=r["seconds"], rc=r["rc"], nk_irr=nk, npool=npool,
             hub={int(k): v for k, v in hub.items()},
             d_moment={int(k): round(v[0] - v[1], 4) for k, v in hub.items()})
    ef2 = re.findall(r"the spin up/dw Fermi energies are\s+([-0-9.]+)\s+([-0-9.]+)", txt)
    if ef2:
        s["fermi_updw_eV"] = [float(x) for x in ef2[-1]]
    out["scf"] = s
    out["scf_done_t"] = time.time()
    save(out)
    log("scf", s["energy_eV"], "M", s["total_mag"], "|M|", s["abs_mag"], "conv", s["converged_scf"], "iter", s["n_iter"], round(r["seconds"]))
    if not p.get("converged_scf"):
        raise RuntimeError("SCF not converged")


def do_nscf(c, out, at, mom, U, pref):
    nel = sum(zvalence(s, "dojo_sr") for s in at.get_chemical_symbols()) - float((c.get("system_extra") or {}).get("tot_charge", 0.0))
    nbnd = int(np.ceil(nel / 2)) + NB_EXTRA
    nk = nk_irr(at, mom, KNSCF[:3], KNSCF[3:])
    npool = pick_npool(nk)
    write(f"{pref}_nscf.in", at, "nscf", mom, U, KNSCF, pref, c.get("system_extra"), nbnd=nbnd)
    r = run_pw(f"{pref}_nscf.in", ncores=NC, npool=npool, timeout=6 * 3600)
    rec = {"sec": r["seconds"], "rc": r["rc"], "nbnd": nbnd, "nelec": nel, "nk_irr": nk, "npool": npool, "grid": KNSCF}
    if r["rc"] != 0 or "JOB DONE" not in open(f"{pref}_nscf.out", errors="replace").read():
        shutil.move(f"{pref}_nscf.out", f"{pref}_nscf_try1.out")
        write(f"{pref}_nscf.in", at, "nscf", mom, U, KNSCF, pref, c.get("system_extra"), nbnd=nbnd, electrons={"diagonalization": "cg"})
        r = run_pw(f"{pref}_nscf.in", ncores=NC, npool=npool, timeout=10 * 3600)
        rec["retry_cg"] = {"sec": r["seconds"], "rc": r["rc"]}
    out["nscf"] = rec
    save(out)
    log("nscf rc", r["rc"], round(r["seconds"]), "nk", nk, "npool", npool)
    if r["rc"] != 0:
        raise RuntimeError("nscf failed")
    return npool


def do_projwfc(c, out, pref, npool, raw):
    sdir = SCRATCH / pref
    open(f"{pref}_proj.in", "w").write(
        "&PROJWFC\n"
        f"  prefix = '{pref}'\n  outdir = '{sdir}'\n  filpdos = '{sdir}/pdos'\n  filproj = '{sdir}/proj'\n"
        "  lsym = .false.\n  DeltaE = 0.02\n  degauss = 0.0037\n/\n")
    cmd = (f"OMP_NUM_THREADS=1 PMIX_MCA_gds=hash mpirun --allow-run-as-root -np {NC} -x OMP_NUM_THREADS --bind-to none "
           f"--mca pml ob1 --mca btl self,sm --mca osc ^ucx projwfc.x -nk {npool} -in {pref}_proj.in > {pref}_proj.out 2>&1")
    t0 = time.time()
    rc = subprocess.run(cmd, shell=True, timeout=4 * 3600).returncode
    out["projwfc"] = {"rc": rc, "sec": time.time() - t0, "npool": npool}
    save(out)
    log("projwfc rc", rc, round(time.time() - t0))
    shutil.copy(sdir / f"{pref}.save" / "data-file-schema.xml", f"{raw}/data-file-schema.xml")
    shutil.copy(f"{pref}_proj.out", f"{raw}/projwfc.out")
    xmls = glob.glob(str(sdir / "**" / "atomic_proj.xml"), recursive=True)
    if not xmls:
        raise RuntimeError("atomic_proj.xml not found")
    with open(xmls[0], "rb") as fi, gzip.open(f"{raw}/atomic_proj.xml.gz", "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    open(f"{raw}/DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S"))


# ------------------------------------------------------------------------------------------------ projections
def parse_states(projout):
    txt = open(projout, errors="replace").read()
    st = re.findall(r"state #\s*(\d+):\s+atom\s+(\d+)\s+\((\S+?)\s*\)\s*,\s+wfc\s+(\d+)\s+\(l=\s*(\d+)\s+m=\s*(\d+)\)", txt)
    return [(int(a) - 1, int(b) - 1, lab, int(w), int(l), int(m)) for a, b, lab, w, l, m in st]


def parse_atomic_proj(path):
    opener = gzip.open if path.endswith(".gz") else open
    recs, cur, header = [], None, {}
    with opener(path, "rb") as fh:
        for ev, el in ET.iterparse(fh, events=("start", "end")):
            tag = el.tag
            if ev == "start":
                if tag == "K-POINT":
                    cur = {"proj": {}}
                    recs.append(cur)
                continue
            if tag == "HEADER":
                header = dict(el.attrib)
            elif tag == "K-POINT":
                cur["k"] = [float(x) for x in el.text.split()]
            elif tag == "E":
                cur["E"] = np.array(el.text.split(), float) * RY_EV
            elif tag == "ATOMIC_WFC":
                v = np.array(el.text.split(), float).reshape(-1, 2)
                cur["proj"][(int(el.attrib.get("spin", 1)), int(el.attrib["index"]) - 1)] = (v[:, 0] ** 2 + v[:, 1] ** 2).astype(np.float32)
                el.clear()
            elif tag == "PROJS":
                el.clear()
    return header, recs


def analyse(c, out, at, raw):
    nat = len(at)
    sym = at.get_chemical_symbols()
    role = c["role"]
    b = parse_xml_bands(f"{raw}/data-file-schema.xml")
    E, occ = b["eig"], b["occ"]
    nk, ns, nb = E.shape
    wk = np.asarray(b["w"], float); wk = wk / wk.sum()
    states = parse_states(f"{raw}/projwfc.out")
    nw = len(states)
    hdr, recs = parse_atomic_proj(f"{raw}/atomic_proj.xml.gz")
    if int(hdr.get("NUMBER_OF_ATOMIC_WFC", nw)) != nw:
        raise RuntimeError(f"state list {nw} != header {hdr.get('NUMBER_OF_ATOMIC_WFC')}")
    P = np.zeros((nk, ns, nb, nw), np.float32)
    filled = np.zeros((nk, ns), bool)
    for r in recs:
        spins_here = sorted({s for s, _ in r["proj"]})
        Er = r["E"]
        for si, s in enumerate(spins_here):
            e = Er[si * nb:(si + 1) * nb] if len(Er) == len(spins_here) * nb else Er[:nb]
            best = None
            for kk in range(nk):
                for ss in range(ns):
                    if filled[kk, ss]:
                        continue
                    dev = np.max(np.abs(E[kk, ss, :len(e)] - e))
                    if best is None or dev < best[0]:
                        best = (dev, kk, ss)
            dev, kk, ss = best
            if dev > 2e-3:
                raise RuntimeError(f"atomic_proj record does not match schema eigenvalues (dev {dev:.4f} eV)")
            for iw in range(nw):
                P[kk, ss, :, iw] = r["proj"][(s, iw)][:nb]
            filled[kk, ss] = True
    if not filled.all():
        raise RuntimeError(f"projections missing for {np.argwhere(~filled).tolist()}")
    # per-atom s/p/d weights  A[nk, ns, nb, nat, 3]
    A = np.zeros((nk, ns, nb, nat, 3), np.float32)
    for iw, (_, a, _, w, l, m) in enumerate(states):
        if l <= 2:
            A[..., a, l] += P[..., iw]
    Aat = A.sum(-1)
    tot = np.maximum(Aat.sum(-1), 1e-6)
    # groups by role (+ defect region)
    groups = {}
    for i, r in enumerate(role):
        groups.setdefault(r, []).append(i)
    region = c.get("region") or []
    gnames = sorted(groups)
    Gw = np.stack([Aat[..., groups[g]].sum(-1) for g in gnames], -1) / tot[..., None]      # fractions
    Rw = (Aat[..., region].sum(-1) / tot) if region else np.zeros(tot.shape)
    # edges (occupation based)
    occd = occ > 0.5
    vb = np.where(occd, E, -np.inf); cb = np.where(~occd, E, np.inf)
    tv = vb.max(axis=(0, 2)); bc = cb.min(axis=(0, 2))
    VBM, CBM = float(tv.max()), float(bc.min())
    sv, sc = int(np.argmax(tv)), int(np.argmin(bc))
    edges = {"VBM_eV": VBM, "CBM_eV": CBM, "gap": CBM - VBM, "vbm_spin": sv, "cbm_spin": sc,
             "win_VB": float(tv[sv] - tv[1 - sv]), "win_CB": float(bc[1 - sc] - bc[sc]),
             "top_occ_per_spin": tv.tolist(), "bot_emp_per_spin": bc.tolist(), "gap_per_spin": (bc - tv).tolist(),
             "n_occ_per_spin": [float((wk[:, None] * occ[:, s]).sum()) for s in range(ns)],
             "ef_xml_eV": b.get("ef_eV")}
    out["edges"] = edges
    # deep alignment references
    pos = at.get_positions(); cell = at.cell.array; inv = np.linalg.inv(cell)
    cen = [np.asarray(f) @ cell for f in (c.get("centers_frac") or [])]

    def dist_to_centres(i):
        if not cen:
            return 99.0
        best = 99.0
        for p0 in cen:
            d = (pos[i] - p0) @ inv; d -= np.round(d)
            best = min(best, float(np.linalg.norm(d @ cell)))
        return best
    dcen = np.array([dist_to_centres(i) for i in range(nat)])
    cnN = [i for i in range(nat) if sym[i] in ("C", "N")]
    far = [i for i in cnN if dcen[i] > R_FAR]
    kfar = [i for i in range(nat) if sym[i] == "K" and dcen[i] > R_FAR]
    refs = {"R_FAR": R_FAR, "n_far_CN_atoms": len(far), "n_far_K": len(kfar)}
    cn_frac = Aat[..., cnN].sum(-1) / tot
    for name, atoms_, sel_fn in (("CN3s_far", far, lambda s: (cn_frac[:, s] > 0.6) & (E[:, s] < VBM - 14.0)),
                                 ("CN3s_all", cnN, lambda s: (cn_frac[:, s] > 0.6) & (E[:, s] < VBM - 14.0)),
                                 ("K3p_far", kfar, lambda s: ((Aat[:, s][..., [i for i in range(nat) if sym[i] == 'K']].sum(-1) / tot[:, s]) > 0.6)
                                  & (E[:, s] < VBM - 9.0) & (E[:, s] > VBM - 20.0))):
        if not atoms_:
            continue
        vals = []
        for s in range(ns):
            m = sel_fn(s)
            W = wk[:, None] * m * A[:, s][..., atoms_, :].sum((-1, -2))
            vals.append(float((W * E[:, s]).sum() / W.sum()) if W.sum() > 1e-8 else None)
        refs[f"E_ref_{name}"] = vals
        refs[f"E_ref_{name}_mean"] = float(np.mean([v for v in vals if v is not None])) if any(v is not None for v in vals) else None
    out["refs"] = refs
    # band table near the gap: bands whose k-range intersects [VBM - 1.5, CBM + 1.5]
    lo, hi = VBM - 1.5, CBM + 1.5
    table = {}
    for s in range(ns):
        rows = []
        for ib in range(nb):
            e = E[:, s, ib]
            if e.max() < lo or e.min() > hi:
                continue
            gw = (wk[:, None] * Gw[:, s, ib]).sum(0)
            aw = (wk[:, None] * (Aat[:, s, ib] / tot[:, s, ib][:, None])).sum(0)
            ta = np.argsort(-aw)[:4]
            rows.append({"band": ib, "Emin": round(float(e.min()), 4), "Emax": round(float(e.max()), 4),
                         "occ": round(float((wk * occ[:, s, ib]).sum()), 3),
                         "groups": {g: round(float(x), 3) for g, x in zip(gnames, gw) if x >= 0.02},
                         "region": round(float((wk * Rw[:, s, ib]).sum()), 3),
                         "d_frac_metals": round(float((wk * (A[:, s, ib][:, [i for i in range(nat) if sym[i] in ('V', 'Cr', 'Mo')], 2].sum(-1) / tot[:, s, ib])).sum()), 3),
                         "top_atoms": [[int(a), sym[a], role[a], round(float(aw[a]), 3)] for a in ta]})
        table[f"s{s}"] = rows
    out["band_table"] = table
    out["groups"] = {g: len(v) for g, v in groups.items()}
    out["dist_to_centre"] = [round(float(x), 3) for x in dcen]
    np.savez_compressed(f"results/{out['id']}_bands.npz", E=E.astype(np.float32), occ=occ.astype(np.float16), wk=wk,
                        k=b["k"], A=A.astype(np.float16), species=np.array(sym), role=np.array(role))
    return out


# ------------------------------------------------------------------------------------------------ main
def main(path):
    c = json.load(open(path))
    cid = c["id"]
    at0 = Atoms(c["species"], cell=c["lattice"], scaled_positions=c["frac"], pbc=True)
    c["_at0"] = at0
    mom = [s * m for s, m in zip(c["spins"], c["mom"])]
    U = {k: tuple(v) for k, v in c["hubbard"].items()}
    rp = f"results/{cid}.json"
    out = json.load(open(rp)) if os.path.exists(rp) else {}
    if out.get("done") and not ANALYZE_ONLY:
        log(cid, "already done"); return
    out.update(id=cid, material=c["material"], defect=c["defect"], formula=c["formula"], hubbard=c["hubbard"],
               M_expected=c.get("M_expected"), system_extra=c.get("system_extra"), region=c.get("region"),
               settings={"ECUT": ECUT, "KRELAX": KRELAX, "KNSCF": KNSCF, "NB_EXTRA": NB_EXTRA, "RELAX_MAXSEC": RELAX_MAXSEC,
                         "NC": NC, "NPOOL_MAX": NPOOL_MAX, "R_FAR": R_FAR}, t_last=time.time())
    out.setdefault("t0", time.time())
    save(out)
    raw = f"raw/{cid}"
    os.makedirs(raw, exist_ok=True)
    if not os.path.exists(f"{raw}/DONE"):
        if ANALYZE_ONLY:
            raise RuntimeError("ANALYZE_ONLY but raw files missing")
        if SKIP_RELAX and not out.get("relax", {}).get("done"):
            out["relax"] = {"done": True, "skipped": True, "frac": at0.get_scaled_positions().tolist()}
        at = do_relax(c, out, at0, mom, U)
        out["relaxed"] = {"lattice": at.cell.array.tolist(), "species": at.get_chemical_symbols(), "frac": at.get_scaled_positions().tolist()}
        pref = f"{cid}_lcm"
        # scf (rerun if the save dir was lost to a preemption before the nscf finished)
        if not (out.get("scf", {}).get("converged_scf") and (SCRATCH / pref / f"{pref}.save" / "data-file-schema.xml").exists()):
            do_scf(c, out, at, mom, U, pref)
        npool = do_nscf(c, out, at, mom, U, pref)
        do_projwfc(c, out, pref, npool, raw)
        shutil.rmtree(SCRATCH / pref, ignore_errors=True)
    else:
        rl = out["relaxed"]
        at = Atoms(rl["species"], cell=rl["lattice"], scaled_positions=rl["frac"], pbc=True)
    analyse(c, out, at, raw)
    out["done"] = True
    out["t1"] = time.time()
    save(out)
    e = out["edges"]; s = out["scf"]
    log(cid, "E %.4f M %.3f |M| %.3f gap %.3f vbm_s %d win_VB %.3f cbm_s %d win_CB %.3f" %
        (s["energy_eV"], s["total_mag"] if s.get("total_mag") is not None else float("nan"), s["abs_mag"], e["gap"],
         e["vbm_spin"], e["win_VB"], e["cbm_spin"], e["win_CB"]))


if __name__ == "__main__":
    try:
        main(sys.argv[1])
    except Exception:
        traceback.print_exc()
        sys.exit(1)
