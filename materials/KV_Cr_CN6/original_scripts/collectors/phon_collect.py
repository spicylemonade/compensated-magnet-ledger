"""Collect + analyse the PBA phonons (Track L): KV[Cr(CN)6] / KV[Mo(CN)6], PBE+U 110 Ry, LCM, F-43m 15-atom cell.

usage: python tracks/lcm/pba/phon/collect.py [--fetch] [--fake]
  --fetch  read jobs/lcm/pba/phon/<m>_<tag>_g<k>/{forces,failed}_<idx>.json (fixed paths from manifest.json, no listdir)
           into results/<m>_<tag>/
  --fake   replace DFT forces by a nearest-neighbour spring model (pipeline test only; writes nothing to summary.json)
Writes summary.json and results/phonons_full.json (+ results/phonopy_params_<m>_<tag>.yaml force constants).

Force constants: phonopy, forces minus the undisplaced reference (index -1), symmetrized (ASR + index symmetry).
Exact frequencies only at the supercell's commensurate q: prim -> Gamma; conv (60 atoms) -> Gamma + 3 X; s222 -> + 4 L.
Mode character (real displacement pattern in the 2x2x2 (120-atom) supercell; kinetic-energy fractions, sum m|u|^2 = 1):
  K           K weight; 'polar' if the 4 K move in phase (Gamma T2 off-centring), else 'antiphase'
  rot_V/rot_M rigid-rotation part of the six ligands about each V(N6) / M(C6) centre (least-squares omega,
              delta_i - delta_centre ~ omega x d_i), KE of omega x d_i on the ligands  -> octahedral rotation / tilt
  CN_bend_tr  transverse centre-of-mass motion of the C#N units (the 'tilt-type' bending of M-C#N-V)
  CN_lib      relative transverse C vs N motion (C#N libration)
  CN_str      C#N stretch;  CN_long  longitudinal C#N translation
"""
import argparse
import json
import pathlib

import numpy as np
import phonopy

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE / "results"
THZ_CM = 33.35641
MATS = ("KVCr", "KVMo")
TAGS = ("prim", "conv", "L30", "s222")
S_MOD = 2 * np.eye(3, dtype=int)  # 120-atom fcc supercell: commensurate with Gamma, X and L (used for mode characters)
QNAMES = {(0.0, 0.0, 0.0): "G", (0.0, 0.5, 0.5): "X", (0.5, 0.0, 0.5): "X", (0.5, 0.5, 0.0): "X",
          (0.5, 0.0, 0.0): "L", (0.0, 0.5, 0.0): "L", (0.0, 0.0, 0.5): "L", (0.5, 0.5, 0.5): "L"}


def split_acoustic(modes, name):
    """Gamma: the 3 modes with the largest bond coherence (uniform translation, coherence -> 1) are acoustic, whatever
    their sign; sorting by frequency would mislabel imaginary optic modes as acoustic. Returns (acoustic, optic) sorted."""
    so = sorted(modes, key=lambda x: x["THz"])
    if name != "G":
        return [], so
    ac = sorted(so, key=lambda x: -(x.get("com_frac") or 0))[:3]
    ids = {x["band"] for x in ac}
    return sorted(ac, key=lambda x: x["THz"]), [x for x in so if x["band"] not in ids]


def qname(q):
    k = tuple(float(abs(round(x * 2)) / 2 % 1) for x in q)
    return QNAMES.get(k, "?")


# ------------------------------------------------------------------------------------------------ fetch / load
def fetch():
    import modal
    vol = modal.Volume.from_name("magdisc-data")
    man = json.loads((HERE / "manifest.json").read_text())
    got = 0
    for job, rec in man["jobs"].items():
        if rec["tag"].startswith("koff"):
            continue
        d = RES / f"{rec['m']}_{rec['tag']}"
        d.mkdir(parents=True, exist_ok=True)
        for i in rec["idx"]:
            if (d / f"forces_{i}.json").exists():
                continue
            for kind in ("forces", "failed"):
                try:
                    b = b"".join(vol.read_file(f"jobs/{job}/{kind}_{i}.json"))
                    (d / f"{kind}_{i}.json").write_bytes(b)
                    got += kind == "forces"
                    break
                except Exception:
                    pass
    print("fetched", got, "new force files")
    for job, rec in man["jobs"].items():
        for path in rec.get("results", []) if rec["tag"].startswith("koff_") else []:
            dst = RES / rec["tag"].split("_")[0] / pathlib.Path(path).name  # results/koff or results/koff2
            if dst.exists():
                continue
            try:
                b = b"".join(vol.read_file(path))
                dst.parent.mkdir(parents=True, exist_ok=True); dst.write_bytes(b)
                print("fetched", dst.name)
            except Exception:
                pass


def kscan_table():
    """rigid K-only scans in the 15-atom cell (all K in phase, framework frozen): E(u) - E(centred), K force along u."""
    out = {}
    for m in MATS:
        ref = RES / f"{m}_prim" / "forces_-1.json"
        cfgp = HERE / "inputs" / f"{m}_kscan.json"
        if not ref.exists() or not cfgp.exists():
            continue
        E0 = json.loads(ref.read_text())["energy_eV"]
        cfg = json.loads(cfgp.read_text())
        meta = {c["idx"]: c for c in cfg["supercells"]}
        iK = cfg["supercells"][0]["species"].index("K")
        rows = {}
        for idx, c in sorted(meta.items()):
            f = RES / f"{m}_kscan" / f"forces_{idx}.json"
            if not f.exists():
                continue
            r = json.loads(f.read_text())
            u = {"x": [1, 0, 0], "p111": [1, 1, 1], "m111": [-1, -1, -1]}[c["dir"]]
            u = np.array(u, float) / np.linalg.norm(u)
            FK = float(np.array(r["forces_eV_A"])[iK] @ u)
            rows.setdefault(c["dir"], []).append({"amp_A": c["amp_A"], "dE_meV": round(1000 * (r["energy_eV"] - E0), 2),
                                                 "F_K_along_u_eV_A": round(FK, 4), "M": r.get("M"), "absM": r.get("absM")})
        if rows:
            out[m] = rows
            print(f"rigid K scan {m} (E - E_centred, meV/f.u.; framework frozen):")
            for dname, rr in rows.items():
                print("   ", dname, [(x["amp_A"], x["dE_meV"], x["F_K_along_u_eV_A"]) for x in rr])
    return out


def koff_peek():
    """unfinished K off-centring relaxes (koff: start 0.25 A; koff2: start at the rigid-scan minimum): last energy vs the
    centred relax, K position, read from the .out on the volume. koff x/p111/m111 were terminated (state 'terminated')."""
    import re
    import modal
    vol = modal.Volume.from_name("magdisc-data")
    man = json.loads((HERE / "manifest.json").read_text())
    out = {}
    for m in MATS:
        fc = RES / "koff" / f"koff_{m}_c.json"
        if not fc.exists():
            continue
        Ec = json.loads(fc.read_text())["energy_eV"]
        cfg = json.loads((HERE / "inputs" / f"{m}_prim.json").read_text())
        sc0 = [c for c in cfg["supercells"] if c["idx"] == -1][0]
        lat = np.array(sc0["lattice"]); iK = sc0["species"].index("K"); k0 = np.array(sc0["frac"][iK])
        for kd, tag in [("koff", t) for t in ("x", "xy", "p111", "m111")] + [("koff2", t) for t in ("x", "p111", "m111")]:
            if (RES / kd / f"koff_{m}_{tag}.json").exists():
                continue
            job = f"lcm/pba/phon/{kd}/{m}_{tag}"
            try:
                txt = b"".join(vol.read_file(f"jobs/{job}/koff_{m}_{tag}.out")).decode(errors="replace")
            except Exception:
                continue
            en = [float(e) * 13.605693122994 for e in re.findall(r"^!\s+total energy\s+=\s+([-0-9.]+)", txt, re.M)]
            L = txt.splitlines()
            bl = [i for i, ln in enumerate(L) if ln.startswith("ATOMIC_POSITIONS")]
            dK = None
            if bl:
                fk = np.array([float(x) for x in L[bl[-1] + 1 + iK].split()[1:4]])
                d = fk - k0; d -= np.round(d); dK = d @ lat
            tf = re.findall(r"Total force =\s+([0-9.]+)", txt)
            state = "running" if job in man["jobs"] else "terminated"
            out.setdefault(m, {})[f"{kd}_{tag}"] = {"n_scf": len(en), "dE_last_meV": round(1000 * (en[-1] - Ec), 2) if en else None,
                                                    "dE_min_meV": round(1000 * (min(en) - Ec), 2) if en else None,
                                                    "K_disp_A": None if dK is None else np.round(dK, 3).tolist(),
                                                    "K_disp_norm_A": None if dK is None else round(float(np.linalg.norm(dK)), 3),
                                                    "total_force_Ry_bohr": float(tf[-1]) if tf else None, "state": state}
            r_ = out[m][f"{kd}_{tag}"]
            print(f"   {kd} {state} {m} {tag}: n_scf {len(en)} dE {r_['dE_last_meV']} meV |dK| {r_['K_disp_norm_A']} A F {r_['total_force_Ry_bohr']}")
    return out


def koff_table():
    """K off-centring relaxes (15-atom cell, all K in phase): E(tag) - E(centred relax), K shift, edges."""
    out = {}
    for m in MATS:
        rows = {}
        for kd, tag in [("koff", t) for t in ("c", "x", "xy", "p111", "m111")] + [("koff2", t) for t in ("x", "p111", "m111")]:
            f = RES / kd / f"koff_{m}_{tag}.json"
            if f.exists():
                rows[tag if kd == "koff" else f"{kd}_{tag}"] = json.loads(f.read_text())
        if not rows:
            continue
        ref = rows.get("c", {}).get("energy_eV")
        tab = {}
        for tag, r in rows.items():
            e = r.get("edges_scf", {})
            tab[tag] = {"dE_meV_per_fu": None if (ref is None or r.get("energy_eV") is None) else round(1000 * (r["energy_eV"] - ref), 2),
                        "K_disp_A": r.get("K_disp_A"), "K_disp_norm_A": r.get("K_disp_norm_A"), "K_nearest": r.get("K_nearest"),
                        "bfgs_converged": r.get("bfgs_converged"), "M": r.get("M"), "absM": r.get("absM"),
                        "gap": e.get("gap"), "win_VB": e.get("win_VB"), "win_CB": e.get("win_CB"),
                        "vbm_spin": e.get("vbm_spin"), "cbm_spin": e.get("cbm_spin")}
        out[m] = tab
        print(f"K off-centring {m}:")
        for tag, t in tab.items():
            print(f"   {tag:5s} dE {t['dE_meV_per_fu']} meV/f.u.  |dK| {t['K_disp_norm_A']}  conv {t['bfgs_converged']}  M {t['M']} |M| {t['absM']}  "
                  f"gap {t['gap']} win {t['win_VB']}/{t['win_CB']} spins {t['vbm_spin']}/{t['cbm_spin']}  nearest {t['K_nearest'][:2] if t['K_nearest'] else None}")
    return out


def load_forces(m, tag):
    d = RES / f"{m}_{tag}"
    out = {}
    if d.exists():
        for p in d.glob("forces_*.json"):
            r = json.loads(p.read_text())
            if r.get("forces_eV_A"):
                out[int(r["idx"])] = r
    return out


def fake_forces(ph):
    """spring model on the perfect supercell (pipeline test)."""
    sc = ph.supercell
    lat = np.array(sc.cell); fr0 = np.array(sc.scaled_positions); sym = list(sc.symbols)
    n = len(sym)
    pairs = []
    for i in range(n):
        d = fr0 - fr0[i]; d -= np.round(d); c = d @ lat; r = np.linalg.norm(c, axis=1)
        for j in np.where((r > 0.1) & (r < 3.9))[0]:
            if j > i:
                k = 30.0 if r[j] < 1.3 else (8.0 if r[j] < 2.3 else 0.4)
                pairs.append((i, j, k, c[j] / r[j]))
    out = {}
    for idx, dsc in enumerate(ph.supercells_with_displacements):
        u = np.array(dsc.scaled_positions) - fr0; u -= np.round(u); u = u @ lat
        F = np.zeros((n, 3))
        for i, j, k, e in pairs:
            du = u[j] - u[i]
            f = k * (du @ e) * e + 0.05 * k * (du - (du @ e) * e)
            F[i] += f; F[j] -= f
        out[idx] = {"forces_eV_A": F.tolist(), "M": 0.0, "absM": 0.0, "conv": True}
    out[-1] = {"forces_eV_A": np.zeros((n, 3)).tolist(), "M": 0.0, "absM": 0.0, "conv": True}
    return out


# ------------------------------------------------------------------------------------------------ geometry / character
def supercell_topology(sc):
    lat = np.array(sc.cell); fr = np.array(sc.scaled_positions); sym = list(sc.symbols)
    n = len(sym)

    def vec(i, j):
        d = fr[j] - fr[i]; d -= np.round(d)
        return d @ lat
    metals = [i for i, s in enumerate(sym) if s in ("V", "Cr", "Mo")]
    octs = []
    for c in metals:
        lig_el = "N" if sym[c] == "V" else "C"
        ligs = [(j, vec(c, j)) for j in range(n) if sym[j] == lig_el and np.linalg.norm(vec(c, j)) < 2.4]
        assert len(ligs) == 6, (sym[c], len(ligs))
        octs.append({"c": c, "el": sym[c], "ligs": ligs})
    cn = []
    for i in range(n):
        if sym[i] != "C":
            continue
        js = [j for j in range(n) if sym[j] == "N" and np.linalg.norm(vec(i, j)) < 1.4]
        assert len(js) == 1
        e = vec(i, js[0]); cn.append((i, js[0], e / np.linalg.norm(e)))
    K = [i for i, s in enumerate(sym) if s == "K"]
    bonds = [(o["c"], j) for o in octs for j, _ in o["ligs"]] + [(i, j) for i, j, _ in cn]
    return {"octs": octs, "cn": cn, "K": K, "bonds": bonds, "sym": sym, "masses": np.array(sc.masses), "frac": fr, "lat": lat}


def real_pattern(u):
    """complex supercell modulation -> real pattern with the phase that maximises its norm"""
    s = np.sum(u * u)
    ph = np.exp(-0.5j * np.angle(s)) if abs(s) > 1e-14 else 1.0
    r = np.real(u * ph)
    if np.linalg.norm(r) < 1e-8:
        r = np.imag(u * ph)
    return r


def character(r, top):
    m = top["masses"]
    ke = float(np.sum(m[:, None] * r ** 2))
    r = r / np.sqrt(ke)
    w = {}
    for el in sorted(set(top["sym"])):
        idx = [i for i, s in enumerate(top["sym"]) if s == el]
        w[el] = float(np.sum(m[idx, None] * r[idx] ** 2))
    # K
    rK = r[top["K"]]
    nK = np.linalg.norm(rK, axis=1)
    if nK.max() > 1e-6:
        mean = rK.mean(axis=0)
        polar = float(np.linalg.norm(mean) ** 2 / np.mean(nK ** 2))  # 1 = all K parallel (in phase), 0 = net zero
        kdir = (mean / np.linalg.norm(mean)).round(3).tolist() if np.linalg.norm(mean) > 1e-6 else None
    else:
        polar, kdir = 0.0, None
    # octahedral rotations
    rot = {"V": 0.0, "M": 0.0}
    omegas = []
    for o in top["octs"]:
        A, b = [], []
        for j, d in o["ligs"]:
            dl = r[j] - r[o["c"]]
            A.append(np.array([[0, d[2], -d[1]], [-d[2], 0, d[0]], [d[1], -d[0], 0]]))  # omega x d = A @ omega
            b.append(dl)
        A = np.vstack(A); b = np.concatenate(b)
        om = np.linalg.lstsq(A, b, rcond=None)[0]
        part = (A @ om).reshape(6, 3)
        ke_rot = float(sum(m[j] * np.sum(p ** 2) for (j, _), p in zip(o["ligs"], part)))
        rot["V" if o["el"] == "V" else "M"] += ke_rot
        omegas.append({"el": o["el"], "pos": (top["frac"][o["c"]] % 1).round(3).tolist(), "omega": om.round(4).tolist()})
    # cyanide
    cn = {"CN_bend_tr": 0.0, "CN_lib": 0.0, "CN_str": 0.0, "CN_long": 0.0}
    for i, j, e in top["cn"]:
        mC, mN = m[i], m[j]; M = mC + mN; mu = mC * mN / M
        com = (mC * r[i] + mN * r[j]) / M; rel = r[i] - r[j]
        cl, rl = com @ e, rel @ e
        cn["CN_long"] += M * cl ** 2; cn["CN_bend_tr"] += M * float(np.sum((com - cl * e) ** 2))
        cn["CN_str"] += mu * rl ** 2; cn["CN_lib"] += mu * float(np.sum((rel - rl * e) ** 2))
    cn = {k: float(v) for k, v in cn.items()}
    out = {"w": {k: round(v, 3) for k, v in w.items()}, "K_polar": round(polar, 3), "K_dir": kdir,
           "rot_V": round(rot["V"], 3), "rot_M": round(rot["M"], 3), **{k: round(v, 3) for k, v in cn.items()}}
    # dominant rotation axis pattern
    oms = np.array([o["omega"] for o in omegas])
    if np.abs(oms).max() > 1e-6 and rot["V"] + rot["M"] > 0.05:
        ax = int(np.argmax(np.sum(oms ** 2, axis=0)))
        out["rot_axis"] = "xyz"[ax]
        out["rot_signs"] = {f"{o['el']}@{o['pos']}": float(np.sign(o["omega"][ax]) * (abs(o["omega"][ax]) > 0.1 * np.abs(oms).max()))
                            for o in omegas}
    # bond coherence: 1 = bonded framework atoms move together (acoustic-like), ~0 random, <0 antiphase
    num = sum(float(np.sum((r[i] - r[j]) ** 2)) for i, j in top["bonds"])
    den = sum(float(np.sum(r[i] ** 2) + np.sum(r[j] ** 2)) for i, j in top["bonds"])
    out["coherence"] = round(1.0 - num / den, 3) if den > 0 else None
    # centre-of-mass fraction: 1 for a uniform translation (acoustic at Gamma), 0 for any optic mode at Gamma
    P = np.sum(m[:, None] * r, axis=0)
    out["com_frac"] = round(float(P @ P / (np.sum(m) * np.sum(m[:, None] * r ** 2))), 4)
    out["label"] = label(out)
    return out


def label(c):
    w = c["w"]
    rot = c["rot_V"] + c["rot_M"]
    if w.get("K", 0) > 0.5:
        return "K rattling (" + ("in-phase / polar off-centring" if c["K_polar"] > 0.8 else "antiphase") + ")"
    if (c.get("coherence") or 0) > 0.85:
        return "acoustic-like (rigid framework translation)"
    if rot > 0.35 and c["CN_bend_tr"] >= c["CN_lib"]:
        return "octahedral tilt (cooperative rotation, CN transverse)"
    if rot > 0.35:
        return "CN libration (ligand shells rotate, C and N opposite)"
    if c["CN_bend_tr"] + c["CN_lib"] > 0.5:
        return "CN bending (" + ("transverse translation" if c["CN_bend_tr"] >= c["CN_lib"] else "libration") + ")"
    if c["CN_str"] > 0.5:
        return "CN stretch"
    if c["CN_long"] + w.get("V", 0) + w.get("Cr", 0) + w.get("Mo", 0) > 0.5:
        return "M-C / V-N stretch or metal vs cyanide translation"
    return "mixed framework"


# ------------------------------------------------------------------------------------------------ analysis
def analyse(m, tag, fake=False):
    ph = phonopy.load(HERE / f"phonopy_disp_{m}_{tag}.yaml", produce_fc=False, log_level=0)
    nd = len(ph.dataset["first_atoms"])
    F = fake_forces(ph) if fake else load_forces(m, tag)
    missing = [i for i in range(-1, nd) if i not in F]
    if missing:
        return {"missing": missing, "n_have": len(F), "n_need": nd + 1}
    ref = np.array(F[-1]["forces_eV_A"])
    Fd = np.array([np.array(F[i]["forces_eV_A"]) - ref for i in range(nd)])
    diag = {"ref_residual_max_eV_A": float(np.abs(ref).max())}
    # magnetic sanity of every displaced SCF
    diag["M_range"] = [float(min(F[i].get("M", 0) for i in F)), float(max(F[i].get("M", 0) for i in F))]
    diag["absM_range"] = [float(min(F[i].get("absM", 0) for i in F)), float(max(F[i].get("absM", 0) for i in F))]
    diag["all_conv"] = bool(all(F[i].get("conv", True) for i in F))
    # +/- pair asymmetry (second-order anharmonicity + noise): |F(+u)+F(-u)| / |F(+u)-F(-u)|
    fa = ph.dataset["first_atoms"]
    asym = []
    for i in range(nd):
        for j in range(i + 1, nd):
            if fa[i]["number"] == fa[j]["number"] and np.allclose(fa[i]["displacement"], -np.array(fa[j]["displacement"])):
                asym.append(float(np.linalg.norm(Fd[i] + Fd[j]) / np.linalg.norm(Fd[i] - Fd[j])))
    diag["pm_asymmetry"] = [round(a, 4) for a in asym]
    ph.forces = Fd
    ph.produce_force_constants()
    # before symmetrisation: lowest Gamma modes and the acoustic-sum-rule residual of the raw force constants
    ph.run_qpoints([[0, 0, 0]]); f_raw = np.sort(ph.get_qpoints_dict()["frequencies"][0])
    diag["gamma_lowest3_rawFC_THz"] = [round(float(x), 3) for x in f_raw[:3]]
    fc = ph.force_constants
    diag["rawFC_ASR_max_eV_A2"] = {"row": float(np.abs(fc.sum(axis=1)).max()), "col": float(np.abs(fc.sum(axis=0)).max())}
    ph.symmetrize_force_constants()
    ph.save(HERE / "results" / f"phonopy_params_{m}_{tag}{'_fake' if fake else ''}.yaml", settings={"force_constants": True})
    try:
        from phonopy.harmonic.dynmat_to_fc import get_commensurate_points
        comm = get_commensurate_points(ph.supercell_matrix)
    except Exception:
        comm = np.array([[0, 0, 0]])
    nb = 3 * len(ph.primitive)
    out = {"supercell_matrix": np.array(ph.supercell_matrix).tolist(), "nat_super": len(ph.supercell), "diag": diag, "q": {}}
    for q in comm:
        q = np.array(q, float); q = ((q + 0.5) % 1.0) - 0.5; q[np.isclose(q, -0.5)] = 0.5
        name = qname(q)
        ph.run_qpoints([q], with_eigenvectors=True)
        f = ph.get_qpoints_dict()["frequencies"][0]
        rec = {"q_prim_frac": q.round(3).tolist(), "name": name, "freq_THz": [round(float(x), 4) for x in f]}
        # irreps
        labels = [None] * nb
        try:
            irr = ph.set_irreps(q, degeneracy_tolerance=0.02)
            if getattr(irr, "_ir_labels", None):
                for bi, lab in zip(irr.band_indices, irr._ir_labels):
                    for b in bi:
                        labels[b] = lab
        except Exception as e:
            rec["irreps_error"] = str(e)[:120]
        # characters via the 2x2x2 modulation (commensurate with Gamma, X, L)
        ph.run_modulations(S_MOD, [[q, b, 1.0, 0.0] for b in range(nb)])
        U, sc = ph.get_modulations_and_supercell()
        top = supercell_topology(sc)
        modes = []
        for b in range(nb):
            c = character(real_pattern(U[b]), top)
            modes.append({"band": b, "THz": round(float(f[b]), 3), "cm1": round(float(f[b]) * THZ_CM, 1), "irrep": labels[b], **c})
        rec["modes"] = modes
        out["q"].setdefault(name, []).append(rec)
    # one representative per star; check the X's agree
    summ = {}
    for name, recs in out["q"].items():
        fr = np.array([r["freq_THz"] for r in recs])
        summ[name] = {"n_points": len(recs), "star_spread_THz": float(np.ptp(fr, axis=0).max()) if len(recs) > 1 else 0.0,
                      "min_THz": float(fr.min()), "lowest_optic": None}
        ac, op = split_acoustic(recs[0]["modes"], name)
        summ[name]["acoustic_THz"] = [mm["THz"] for mm in ac] if name == "G" else None
        summ[name]["lowest_optic"] = [{k: mm[k] for k in ("THz", "cm1", "irrep", "label", "w", "K_polar", "rot_V", "rot_M",
                                                           "CN_bend_tr", "CN_lib", "CN_str", "coherence")} for mm in op[:9]]
        summ[name]["n_imag_below_-0.1THz"] = int(sum(1 for mm in op if mm["THz"] < -0.1))
        summ[name]["min_optic_THz"] = float(op[0]["THz"])
    out["summary"] = summ
    # interpolated (indicative only) mesh / path minimum from this supercell's force constants
    if tag == "conv":
        ph.run_mesh([8, 8, 8], with_eigenvectors=False)
        mfreq, mq = ph.mesh.frequencies, ph.mesh.qpoints
        out["interp_mesh8_min_THz"] = float(mfreq.min())
        fmin_q = mq[int(np.argmin(mfreq.min(axis=1)))]
        out["interp_mesh8_argmin_q"] = np.round(fmin_q, 3).tolist()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--fake", action="store_true")
    ap.add_argument("--tags", default="prim,conv,L30,s222")
    a = ap.parse_args()
    if a.fetch:
        fetch()
    res = {}
    for m in MATS:
        for tag in a.tags.split(","):
            if not (HERE / f"phonopy_disp_{m}_{tag}.yaml").exists():
                continue
            if tag == "s222" and not a.fake and not (RES / f"{m}_{tag}").exists():
                continue
            r = analyse(m, tag, fake=a.fake)
            res[f"{m}_{tag}"] = r
            if "missing" in r:
                print(f"{m}_{tag}: forces {r['n_have']}/{r['n_need']}, missing {r['missing']}")
                continue
            print(f"== {m}_{tag} ({r['nat_super']} atoms) diag {json.dumps(r['diag'])}")
            for name, s in r["summary"].items():
                print(f"  {name}: min optic {s['min_optic_THz']:+.3f} THz  n_imag {s['n_imag_below_-0.1THz']}  star spread {s['star_spread_THz']:.4f}"
                      + (f"  acoustic {s['acoustic_THz']}" if s["acoustic_THz"] else ""))
                for mm in s["lowest_optic"][:6]:
                    print(f"     {mm['THz']:+7.3f} THz {mm['cm1']:7.1f} cm-1 {str(mm['irrep']):>4s}  {mm['label']:<42s} w {mm['w']} "
                          f"rotV {mm['rot_V']} rotM {mm['rot_M']} CNtr {mm['CN_bend_tr']} CNlib {mm['CN_lib']} Kpol {mm['K_polar']} coh {mm['coherence']}")
            if "interp_mesh8_min_THz" in r:
                print(f"  interpolated 8x8x8 mesh min {r['interp_mesh8_min_THz']:+.3f} THz at {r['interp_mesh8_argmin_q']} (indicative)")
    if a.fake:
        return
    # Gamma consistency prim vs conv
    for m in MATS:
        p = res.get(f"{m}_prim", {})
        for t in ("conv", "L30"):
            c = res.get(f"{m}_{t}", {})
            if "q" in p and "q" in c:
                fp = np.sort(p["q"]["G"][0]["freq_THz"]); fc = np.sort(c["q"]["G"][0]["freq_THz"])
                res[f"{m}_gamma_prim_vs_{t}_maxdiff_THz"] = float(np.abs(fp - fc).max())
                res[f"{m}_gamma_prim_vs_{t}_low9_diff_THz"] = [round(float(x), 3) for x in (fc - fp)[3:12]]
                print(f"{m}: Gamma prim (k 5^3) vs {t} max |df| = {np.abs(fp - fc).max():.3f} THz; lowest optic diffs "
                      f"{res[f'{m}_gamma_prim_vs_{t}_low9_diff_THz']}")
    res["koff"] = koff_table()
    res["kscan"] = kscan_table()
    if a.fetch:
        try:
            res["koff_running"] = koff_peek()
        except Exception as e:
            print("koff peek failed", e)
    (HERE / "results" / "phonons_full.json").write_text(json.dumps(res, indent=1))
    make_summary(res)


CLASSES = {"K": lambda mm: mm["label"].startswith("K rattling"),
           "tilt": lambda mm: mm["label"].startswith("octahedral tilt"),
           "CN_bend": lambda mm: mm["label"].startswith("CN bending")}


def make_summary(res, tol=-0.1):
    out = {"protocol": "PBE+U (ortho-atomic U V 3 / Cr 3 / Mo 2), PseudoDojo SR NC, 110 Ry, LCM (V down / M up), phonopy finite "
                       "displacements 0.02 A, plus-minus pairs, forces minus undisplaced reference, FC symmetrized",
           "cells": {"prim": "15-atom F-43m primitive, k 5x5x5 Gamma -> Gamma",
                     "conv": "60-atom conventional cube, k 2x2x2 shifted -> Gamma + X",
                     "L30": "30-atom index-2 cell (reduced basis 7.6/7.6/13.1 A), k 3x3x2 Gamma -> Gamma + L"},
           "imag_threshold_THz": tol, "materials": {}}
    for m in MATS:
        d = {"exact": {}, "lowest_by_class": {}, "diag": {}, "missing": {}}
        nimag, nimag_fw, fmin = 0, 0, None
        for tag in ("prim", "conv", "L30"):
            r = res.get(f"{m}_{tag}")
            if not r:
                continue
            if "missing" in r:
                d["missing"][tag] = r["missing"]
                continue
            d["diag"][tag] = r["diag"]
            for name, recs in r["q"].items():
                rec = recs[0]
                ac, op = split_acoustic(rec["modes"], name)
                modes = sorted(rec["modes"], key=lambda x: x["THz"])
                key = f"{name}[{tag}]"
                d["exact"][key] = {"q_prim_frac": rec["q_prim_frac"], "n_star_points": len(recs),
                                   "star_spread_THz": r["summary"][name]["star_spread_THz"],
                                   "acoustic_THz": [x["THz"] for x in ac] if name == "G" else None,
                                   "freq_THz": [x["THz"] for x in modes],
                                   "lowest_optic": [{k: x[k] for k in ("THz", "cm1", "irrep", "label", "w", "rot_V", "rot_M",
                                                                        "CN_bend_tr", "CN_lib", "K_polar", "coherence")}
                                                    for x in op[:8]]}
                n = sum(1 for x in op if x["THz"] < tol)
                nf = sum(1 for x in op if x["THz"] < tol and not x["label"].startswith("K rattling"))
                nimag += n
                nimag_fw += nf
                d["exact"][key]["n_imag"] = n
                d["exact"][key]["n_imag_framework"] = nf
                fw = [x for x in op if not x["label"].startswith("K rattling")]
                d["exact"][key]["lowest_framework_THz"] = fw[0]["THz"]
                d["exact"][key]["lowest_framework_label"] = fw[0]["label"]
                lo = op[0]["THz"]
                fmin = lo if fmin is None else min(fmin, lo)
                for cl, fn in CLASSES.items():
                    sel = [x for x in op if fn(x)]
                    if sel:
                        d["lowest_by_class"].setdefault(cl, {})[key] = {"THz": sel[0]["THz"], "cm1": sel[0]["cm1"],
                                                                         "irrep": sel[0]["irrep"], "label": sel[0]["label"]}
            if "interp_mesh8_min_THz" in r:
                d["interp_conv_mesh8_min_THz"] = r["interp_mesh8_min_THz"]
                d["interp_conv_mesh8_argmin_q"] = r["interp_mesh8_argmin_q"]
        for t in ("conv", "L30"):
            if f"{m}_gamma_prim_vs_{t}_maxdiff_THz" in res:
                d[f"gamma_prim_vs_{t}_maxdiff_THz"] = res[f"{m}_gamma_prim_vs_{t}_maxdiff_THz"]
        d["n_imag_exact"] = nimag
        d["n_imag_exact_framework"] = nimag_fw
        d["min_optic_or_zone_boundary_THz"] = fmin
        done = [t for t in ("prim", "conv", "L30") if t in d["diag"]]
        d["verdict"] = ("incomplete" if not done else
                        ("UNSTABLE at an exact q" if nimag else "dynamically stable at all exact q computed") +
                        f" ({', '.join(done)})")
        if done:
            d["verdict_framework"] = (("framework UNSTABLE" if nimag_fw else "framework (all non-K modes) stable") +
                                      f" at {', '.join(k for k in d['exact'])}; " +
                                      (f"K rattling imaginary ({nimag - nimag_fw} modes): ordered centred-K F-43m is a saddle"
                                       if nimag - nimag_fw else "K rattling real"))
        if m in res.get("koff", {}):
            d["K_offcentring_relax_15atom"] = res["koff"][m]
        if m in res.get("koff_running", {}):
            d["K_offcentring_relax_15atom_running"] = res["koff_running"][m]
        if m in res.get("kscan", {}):
            d["K_rigid_scan_15atom"] = res["kscan"][m]
        out["materials"][m] = d
    (HERE / "summary.json").write_text(json.dumps(out, indent=1))
    print("wrote summary.json:", {m: v.get("verdict_framework", v["verdict"]) for m, v in out["materials"].items()})


if __name__ == "__main__":
    main()
