"""Collect Track L PBA Gate 0 (tracks/lcm/pba/README.md; spec tracks/lcm/trilemma/JUDGE.md section 2-3).
usage: python tracks/lcm/pba/collect.py [--fetch] [--json out.json]
  --fetch : download the FIXED result paths below from the Modal volume (no listdir) into tracks/lcm/pba/results/
            (mirrors the volume path); otherwise only the local cache is read. Missing files are reported, never fatal.
Prints: relaxed geometry; per (material, U point / HSE): M_cell, |M|, metal moments (+ flags), gap, global spin
windows (dense nscf when present) with edge channel labels, E(FM) - E(LCM) per magnetic ion; J ratios
dE(VMo)/dE(VCr) and calibration dE(VCr)/dE(CrCr) per method with ratio-scaled T; isomer energy; SOC residue;
hydrated-vs-anhydrous windows; and the JUDGE pass / kill checks."""
import argparse
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
CACHE = HERE / "results"
P = "jobs/lcm/pba"

# --- fixed result paths (same as manifest.json) ---------------------------------------------------------------------
GRID = {"KVCr": ["U2_2", "U2_4", "U4_2", "U4_4"], "KVMo": ["U2_0", "U2_2", "U4_0", "U4_2"],
        "KCrV": ["U2_2", "U2_4", "U4_2", "U4_4"], "CrCr": ["U2", "U4"]}
BASE_TAG = {"KVCr": "U3_3", "KVMo": "U3_2", "KCrV": "U3_3", "CrCr": "U3"}
# (M_N index, M_C index, M_N element, M_C element)
SITES = {"KVCr": (1, 2, "V", "Cr"), "KVMo": (1, 2, "V", "Mo"), "KCrV": (1, 2, "Cr", "V"), "CrCr": (0, 1, "Cr(N)", "Cr(C)"),
         "KVCr_2H2O": (1, 2, "V", "Cr")}
PATHS = {}
for m, tags in GRID.items():
    PATHS[(m, "s1", BASE_TAG[m])] = f"{P}/s1a/{m}/results/{m}.json"
    PATHS[(m, "sym", None)] = f"{P}/s1a/{m}/results/{m}_sym.json"
    for t in tags:
        PATHS[(m, "s1", t)] = f"{P}/s1a/{m}/results/{m}_{t}.json"
PATHS[("KVCr_2H2O", "s1", "U3_3")] = f"{P}/s1/KVCr_2H2O/results/KVCr_2H2O.json"
for m in ("KVCr", "KVMo", "CrCr"):
    for o in ("LCM", "FM"):
        PATHS[(m, "hse", o)] = f"{P}/hse/{m}_hse_{o}/deep_hse.json"
# fixed-moment (tot_magnetization = 6) d3/d3 FM references (pba_fmfix.py): the unconstrained FM is charge-transferred /
# collapsed in plain PBE (HSE start) and in PBE+U at U_Mo = 0; KVCr_hse_FM / KVMo_hse_FM were terminated (manifest)
for m in ("KVCr", "KVMo"):
    PATHS[(m, "hse", "FMfix")] = f"{P}/hse/{m}_hse_FMfix/results/{m}_hse_FMfix.json"
    for t in [BASE_TAG[m]] + GRID[m]:
        PATHS[(m, "fmfix", t)] = f"{P}/s1a/FMfix/results/{m}_{t}_FMfix.json"
# SCF outputs of lcm_stage1 (fixed names <id>_<config>.out) -> ortho-atomic Hubbard d occupations Tr[ns] (up, down);
# parsed on --fetch and cached as small *.hub.json (the .out files themselves are not kept)
OUTS = {}
for m, tags in GRID.items():
    for t in [BASE_TAG[m]] + tags:
        pre = m if t == BASE_TAG[m] else f"{m}_{t}"
        for o in ("LCM", "FM"):
            OUTS[(m, t, o)] = f"{P}/s1a/{m}/{pre}_{o}.out"
for o in ("LCM", "FM"):
    OUTS[("KVCr_2H2O", "U3_3", o)] = f"{P}/s1/KVCr_2H2O/KVCr_2H2O_{o}.out"
SOCOUT = {t: f"{P}/soc/KVMo_soc/soc_{t}.out" for t in ("z", "x")}  # noncollinear per-site moments (parse_pw misses them)
PATHS[("KVMo", "soc", "deep")] = f"{P}/soc/KVMo_soc/deep_soc.json"
PATHS[("KVMo", "soc", "orbm")] = f"{P}/soc/KVMo_soc/orbm.json"
# 2026-10-03 04:55: the lorbm nscf inside KVMo_soc failed (iosys: "Berry Phase/electric fields only for insulators!",
# smeared occupations); the orbital step was rerun with occupations fixed as KVMo_soc_r, which is preferred when present
PATHS[("KVMo", "soc", "orbm_r")] = f"{P}/soc/KVMo_soc_r/orbm.json"
STATUS = [f"{P}/s1a/{m}/status.json" for m in GRID] + [f"{P}/s1/KVCr_2H2O/status.json"] + \
         [f"{P}/hse/{m}_hse_LCM/status.json" for m in ("KVCr", "KVMo", "CrCr")] + [f"{P}/hse/CrCr_hse_FM/status.json"] + \
         [f"{P}/hse/{m}_hse_FMfix/status.json" for m in ("KVCr", "KVMo")] + [f"{P}/s1a/FMfix/status.json"] + \
         [f"{P}/soc/KVMo_soc/status.json", f"{P}/soc/KVMo_soc_r/status.json"]

T_VCR, T_CRCR = 376.0, 240.0  # K, measured (Holmes & Girolami 1999; Schart 2024)
MLO, MHI = 2.6, 3.0           # Gate-0 window for |m| of both metals (muB)


def parse_hub(txt):
    """last 'Tr[ns(i)] (up, down, total)' per Hubbard atom (QE 1-based index) -> {i: [up, down, total]}"""
    import re
    out = {}
    for i, u, d, t in re.findall(r"Tr\[ns\(\s*(\d+)\)\]\s*\(up, down, total\)\s*=\s*([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", txt):
        out[int(i)] = [float(u), float(d), float(t)]
    return out


def parse_nc_sites(txt):
    """last noncollinear per-atom block: 'atom number N ...' + 'magnetization : mx my mz' -> {N: [mx, my, mz]}"""
    import re
    out = {}
    for n, x, y, z in re.findall(r"atom number\s+(\d+)\s+relative position.*?\n.*?\n\s+magnetization :\s+([-0-9.]+)\s+([-0-9.]+)\s+([-0-9.]+)", txt):
        out[int(n)] = [float(x), float(y), float(z)]
    return out


def fetch():
    import importlib.util
    spec = importlib.util.spec_from_file_location("mrun", HERE.parents[2] / "infra" / "mrun.py")
    mrun = importlib.util.module_from_spec(spec); spec.loader.exec_module(mrun)
    got, miss = 0, []
    for p in list(PATHS.values()) + STATUS:
        try:
            b = mrun._read(p)
        except Exception:
            miss.append(p); continue
        dst = CACHE / p; dst.parent.mkdir(parents=True, exist_ok=True); dst.write_bytes(b); got += 1
    for p in OUTS.values():
        dst = CACHE / (p + ".hub.json")
        if dst.exists():
            continue
        try:
            txt = mrun._read(p).decode(errors="replace")
        except Exception:
            miss.append(p); continue
        if "JOB DONE" not in txt:
            continue  # still running: parse again next time
        dst.parent.mkdir(parents=True, exist_ok=True); dst.write_text(json.dumps(parse_hub(txt))); got += 1
    for p in SOCOUT.values():
        dst = CACHE / (p + ".sites.json")
        if dst.exists():
            continue
        try:
            txt = mrun._read(p).decode(errors="replace")
        except Exception:
            miss.append(p); continue
        if "JOB DONE" in txt:
            dst.parent.mkdir(parents=True, exist_ok=True); dst.write_text(json.dumps(parse_nc_sites(txt))); got += 1
    print(f"fetched {got} files; missing {len(miss)}")
    return miss


def hub(m, t, o):
    try:
        return {int(k): v for k, v in json.loads((CACHE / (OUTS[(m, t, o)] + ".hub.json")).read_text()).items()}
    except Exception:
        return None


def load(key):
    f = CACHE / PATHS[key]
    try:
        d = json.loads(f.read_text())
    except Exception:
        return None
    # lcm_deep.py writes deep_<mode>.json flat while running, then rewrites it at exit as {"id": ..., "<mode>": {...}}
    # (seen 2026-10-03 on KVMo_soc): unwrap so finished and running files read the same
    if isinstance(d, dict) and "id" in d:
        for mode in ("hse", "soc"):
            if isinstance(d.get(mode), dict):
                return d[mode]
    return d


def f2(x, n=2):
    return "  -  " if x is None else (f"{x:.{n}f}" if isinstance(x, (int, float)) else str(x))


def chan(m, run):
    """spin-index -> label (from the actual sign of the M_N moment)."""
    iN, iC, eN, eC = SITES[m]
    sm = run.get("site_moments") or []
    if len(sm) > max(iN, iC) and isinstance(sm[iN], (int, float)):
        nmaj = 0 if sm[iN] > 0 else 1
    else:
        nmaj = 1
    return {nmaj: f"{eN}-maj", 1 - nmaj: f"{eC}-maj"}, nmaj


def summarize(m, run, label="", hubocc=None):
    """one row dict for a collinear run (stage-1 run dict or HSE dict). Moment check |m| 2.6-3.0 uses the ortho-atomic
    Hubbard d-occupation moments (Tr ns up - down) when available; QE sphere moments (R ~ 0.11 alat ~ 0.86 A) understate
    them (KVMo first SCF: sphere V -1.90 / Mo 1.42 vs d-occupation 2.58 / 2.44), so HSE (no U) is judged on M_FM = 6,
    M_LCM = 0 and sphere moments relative to PBE+U."""
    if not run:
        return None
    iN, iC, eN, eC = SITES[m]
    sm = run.get("site_moments") or []
    mN = sm[iN] if len(sm) > iN else None
    mC = sm[iC] if len(sm) > iC else None
    g = run.get("glob_nscf") or run.get("glob_scf") or run.get("glob") or {}
    src = "nscf" if run.get("glob_nscf") else ("scf" if (run.get("glob_scf") or run.get("glob")) else "-")
    gap = (run.get("edges_nscf") or run.get("edges_scf") or {}).get("gap", g.get("gap"))
    if gap is None and g.get("gap_up") is not None:
        gap = min(g["gap_up"], g["gap_dn"])
    lab, nmaj = chan(m, run)
    flags = []
    dN = dC = None
    if hubocc and (iN + 1) in hubocc and (iC + 1) in hubocc:
        dN = hubocc[iN + 1][0] - hubocc[iN + 1][1]; dC = hubocc[iC + 1][0] - hubocc[iC + 1][1]
        for el, x in ((eN, dN), (eC, dC)):
            if not (MLO <= abs(x) <= MHI):
                flags.append(f"|m_d({el})|={abs(x):.2f}")
    if run.get("converged_scf") is False:
        flags.append("SCF-NOT-CONVERGED")
    if run.get("n_flipped"):
        flags.append(f"flipped={run['n_flipped']}")
    if run.get("n_collapsed"):
        flags.append(f"collapsed={run['n_collapsed']}")
    return {"label": label, "E": run.get("energy_eV"), "M": run.get("total_mag"), "absM": run.get("abs_mag"), "mN": mN, "mC": mC,
            "mdN": dN, "mdC": dC, "nd": ([hubocc[iN + 1][2], hubocc[iC + 1][2]] if dN is not None else None),
            "gap": gap, "vb": lab.get(g.get("vbm_spin")), "cb": lab.get(g.get("cbm_spin")), "vbm_spin": g.get("vbm_spin"),
            "cbm_spin": g.get("cbm_spin"), "nmaj": nmaj, "winVB": g.get("win_VB"), "winCB": g.get("win_CB"),
            "gap_up": g.get("gap_up"), "gap_dn": g.get("gap_dn"), "src": src, "flags": flags}


def pair(m, lcm, fm, label, hl=None, hf=None):
    a = summarize(m, lcm, label, hl); b = summarize(m, fm, label, hf)
    dE = None
    if a and b and a["E"] is not None and b["E"] is not None:
        dE = (b["E"] - a["E"]) * 1000 / 2.0  # meV per magnetic ion (2 per primitive cell)
    if a is not None:
        if a["M"] is not None and abs(a["M"]) > 0.02:
            a["flags"].append(f"M_LCM={a['M']:.3f}")
        if b and b["M"] is not None and abs(abs(b["M"]) - 6.0) > 0.05:
            a["flags"].append(f"M_FM={b['M']:.2f} (6 = d3/d3; 4 = charge transfer)")
    return {"lcm": a, "fm": b, "dE": dE, "label": label}


def print_row(m, r):
    a, b = r["lcm"], r["fm"]
    if a is None:
        print(f"  {m:5s} {r.get('label', ''):8s} (no LCM result)"); return
    print(f"  {m:5s} {a['label']:8s} M={f2(a['M'], 3):>6s} |M|={f2(a['absM']):>5s} m_N={f2(a['mN']):>6s} m_C={f2(a['mC']):>6s} "
          f"md_N={f2(a['mdN']):>6s} md_C={f2(a['mdC']):>6s} "
          f"gap={f2(a['gap']):>5s} VBM:{a['vb'] or '-':8s} CBM:{a['cb'] or '-':8s} winVB={f2(a['winVB']):>5s} "
          f"winCB={f2(a['winCB']):>5s} ({a['src']}) | M_FM={f2(b['M'] if b else None):>5s} gap_FM={f2(b['gap'] if b else None):>5s} "
          f"| dE(FM-LCM)={f2(r['dE'], 1):>6s} meV/ion  {' '.join(a['flags'])}")


def unipolar_ok(a, thr=0.3):
    """both edges in the N-bound-metal (V) majority channel, both windows >= thr."""
    if not a or a["vbm_spin"] is None:
        return None
    return a["vbm_spin"] == a["nmaj"] and a["cbm_spin"] == a["nmaj"] and a["winVB"] >= thr and a["winCB"] >= thr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()
    if a.fetch:
        fetch()
    out = {"pbeU": {}, "hse": {}, "relax": {}}

    print("== status ==")
    for s in STATUS:
        try:
            st = json.loads((CACHE / s).read_text())
            print(f"  {s.split('lcm/pba/')[1].rsplit('/', 1)[0]:22s} {st.get('state')} rc={st.get('rc')} h={f2(st.get('hours'))}")
        except Exception:
            print(f"  {s.split('lcm/pba/')[1].rsplit('/', 1)[0]:22s} (no status cached)")

    print("\n== relaxed geometry (PBE+U base U, 110 Ry, LCM; exact F-43m / Fm-3m projection) ==")
    for m in GRID:
        sym = load((m, "sym", None)); base = load((m, "s1", BASE_TAG[m])) or {}
        if not sym:
            print(f"  {m}: not yet"); continue
        i = sym["info"]; rp = sym.get("relax_parse") or {}
        print(f"  {m:5s} src={sym['source']:8s} a={i['a']:.4f} A  M_N-N={i['bonds']['M_N-N']:.3f} C-N={i['bonds']['C-N']:.3f} "
              f"M_C-C={i['bonds']['M_C-C']:.3f}  P={f2(rp.get('pressure_kbar'), 1)} kbar bfgs={rp.get('bfgs_converged')} "
              f"dev={i['max_atom_dev_A']:.1e} A restarts={base.get('relax_restarts', 0)}")
        out["relax"][m] = {"a": i["a"], **i["bonds"], "source": sym["source"]}

    print("\n== PBE+U: LCM (C-end +, N-end -) and FM at fixed geometry; windows from the dense nscf (8x8x8) ==")
    print("   m_N/m_C = QE sphere moments (R ~ 0.86 A); md_N/md_C = ortho-atomic Hubbard d moments Tr ns (up - down)")
    for m in GRID:
        out["pbeU"][m] = {}
        for t in [BASE_TAG[m]] + GRID[m]:
            d = load((m, "s1", t))
            if not d:
                print(f"  {m:5s} {t:8s} (missing)"); continue
            r = pair(m, d["runs"].get("LCM"), d["runs"].get("FM"), t, hub(m, t, "LCM"), hub(m, t, "FM"))
            print_row(m, r)
            fx = load((m, "fmfix", t)) if (m, "fmfix", t) in PATHS else None
            if fx and fx.get("energy_eV") is not None and r["lcm"] and r["lcm"]["E"] is not None:
                r["dEfix"] = (fx["energy_eV"] - r["lcm"]["E"]) * 1000 / 2.0
                r["fmfix"] = {"M": fx.get("total_mag"), "conv": fx.get("converged_scf"), "gap": (fx.get("glob") or {}).get("gap")}
                print(f"  {'':5s} {'':8s} FM fixed M=6: dE(FMfix-LCM)={r['dEfix']:.1f} meV/ion (M={fx.get('total_mag')}, "
                      f"conv={fx.get('converged_scf')}, gap={f2((fx.get('glob') or {}).get('gap'))})")
            out["pbeU"][m][t] = r
    h = load(("KVCr_2H2O", "s1", "U3_3"))
    if h:
        r = pair("KVCr_2H2O", h["runs"].get("LCM"), h["runs"].get("FM"), "U3_3", hub("KVCr_2H2O", "U3_3", "LCM"),
                 hub("KVCr_2H2O", "U3_3", "FM"))
        print("  hydrated KVCr.2H2O (P1, SCF 3x3x3, nscf 6x6x6):"); print_row("KVCr_2H2O", r)
        out["pbeU"]["KVCr_2H2O"] = {"U3_3": r}

    print("\n== HSE06 (no U, nq 2x2x2, k 4x4x4, 90 Ry) on the exactly symmetrized relaxed cells; windows on the SCF grid ==")
    print("   FM reference: KVCr / KVMo = fixed M = 6 (tot_magnetization); CrCr = unconstrained (M = 6.00)")
    for m in ("KVCr", "KVMo", "CrCr"):
        L = load((m, "hse", "LCM"))
        F = load((m, "hse", "FMfix")) if m in ("KVCr", "KVMo") else load((m, "hse", "FM"))
        if not L and not F:
            print(f"  {m}: not yet"); continue
        r = pair(m, L, F, "HSE")
        print_row(m, r)
        out["hse"][m] = r

    # ---------------------------------------------------------------- ratios
    def dE(meth, m, t=None):
        r = (out["hse"].get(m) if meth == "hse" else out["pbeU"].get(m, {}).get(t)) or {}
        if meth == "fix":  # fixed-moment FM reference (KVCr / KVMo); CrCr free FM (always M = 6.00)
            return r.get("dEfix") if m in ("KVCr", "KVMo") else r.get("dE")
        return r.get("dE")

    def ratio(x, y):
        return None if (x is None or y is None or y == 0) else x / y

    print("\n== exchange-energy ratios (E(FM)-E(LCM) per magnetic ion) and ratio-scaled T ==")
    rows = []
    # VMo / VCr
    pairs = [("PBE+U base", ("KVMo", "U3_2"), ("KVCr", "U3_3"))]
    pairs += [(f"PBE+U U_V={uv} M-low", ("KVMo", f"U{uv}_0"), ("KVCr", f"U{uv}_2")) for uv in (2, 4)]
    pairs += [(f"PBE+U U_V={uv} M-high", ("KVMo", f"U{uv}_2"), ("KVCr", f"U{uv}_4")) for uv in (2, 4)]
    pairs += [(f"PBE+U U_V={uv} same U_M=2", ("KVMo", f"U{uv}_2"), ("KVCr", f"U{uv}_2")) for uv in (2, 4)]
    for meth, tagm in (("pbe", "free FM"), ("fix", "FM M=6")):
        for lab, (m1, t1), (m2, t2) in pairs:
            q = ratio(dE(meth, m1, t1), dE(meth, m2, t2))
            rows.append(("VMo/VCr", f"{lab} [{tagm}]", q))
    allq = [ratio(dE("fix", "KVMo", t1), dE("fix", "KVCr", t2)) for t1 in GRID["KVMo"] + ["U3_2"] for t2 in GRID["KVCr"] + ["U3_3"]
            if t1[1] == t2[1]]  # same U_V, fixed-moment FM
    allq = [q for q in allq if q is not None]
    rows.append(("VMo/VCr", "HSE06", ratio(dE("hse", "KVMo"), dE("hse", "KVCr"))))
    # calibration VCr / CrCr
    for meth, tagm in (("pbe", "free FM"), ("fix", "FM M=6")):
        for lab, t1, t2 in (("PBE+U base (3/3 vs 3)", "U3_3", "U3"), ("PBE+U U=2 diagonal", "U2_2", "U2"), ("PBE+U U=4 diagonal", "U4_4", "U4")):
            rows.append(("VCr/CrCr", f"{lab} [{tagm}]", ratio(dE(meth, "KVCr", t1), dE(meth, "CrCr", t2))))
    rows.append(("VCr/CrCr", "HSE06", ratio(dE("hse", "KVCr"), dE("hse", "CrCr"))))
    for kind, lab, q in rows:
        extra = ""
        if q is not None and kind == "VMo/VCr":
            extra = f"T_pred(VMo) = {T_VCR * q:.0f} K (x 376 K)"
        if q is not None and kind == "VCr/CrCr":
            extra = f"target 1.57 +- 0.25 -> {'OK' if abs(q - 1.57) <= 0.25 else 'OFF'}; T_pred(VCr) from CrCr = {T_CRCR * q:.0f} K"
        print(f"  {kind:8s} {lab:40s} {f2(q, 3):>6s}  {extra}")
    if allq:
        print(f"  VMo/VCr  all same-U_V pairs [FM M=6]: min {min(allq):.3f} max {max(allq):.3f} (n={len(allq)})")
    out["ratios"] = [{"kind": k, "method": l, "ratio": q} for k, l, q in rows]

    # ---------------------------------------------------------------- isomer
    print("\n== linkage isomer: E(KCr[V(CN)6]) - E(KV[Cr(CN)6] LCM) per f.u. (each at its own relaxed cell) ==")
    iso = {}
    for t in ["U3_3"] + GRID["KCrV"]:
        a_, b_ = load(("KCrV", "s1", t)), load(("KVCr", "s1", t))
        if not a_ or not b_:
            continue
        ei = [v.get("energy_eV") for v in a_["runs"].values() if v.get("energy_eV") is not None]
        e0 = (b_["runs"].get("LCM") or {}).get("energy_eV")
        if ei and e0 is not None:
            best = min(a_["runs"], key=lambda k: a_["runs"][k].get("energy_eV") or 1e9)
            iso[t] = min(ei) - e0
            print(f"  {t:6s} dE_iso = {iso[t]:+.3f} eV/f.u. (isomer ground config {best})  {'OK' if iso[t] >= 0.8 else 'BELOW 0.8'}")
    out["isomer_eV"] = iso

    # ---------------------------------------------------------------- SOC
    print("\n== SOC (KVMo LCM, dojo_fr, PBE+U 3/2, noncollinear) ==")
    soc, orbm = load(("KVMo", "soc", "deep")), (load(("KVMo", "soc", "orbm_r")) or load(("KVMo", "soc", "orbm")))
    soc_res = None
    if soc:
        for tag in ("z", "x"):
            s = soc.get(tag) or {}
            mv = s.get("total_mag_vec")
            try:
                sm = {int(k): v for k, v in json.loads((CACHE / (SOCOUT[tag] + ".sites.json")).read_text()).items()}
            except Exception:
                sm = {}
            mV = sm.get(2); mMo = sm.get(3)  # QE 1-based: 2 = V (N-end), 3 = Mo (C-end); sphere moments
            nrm = (sum(x * x for x in mv) ** 0.5) if mv else None
            soc_res = max(soc_res or 0, nrm or 0) if nrm is not None else soc_res
            print(f"  m||{tag}: M_spin vec={mv} |M_spin|={f2(nrm, 3)} muB/f.u.  m_V={mV} m_Mo={mMo} conv={s.get('converged_scf')} "
                  f"gap={f2((s.get('glob') or {}).get('gap'))}")
        print(f"  MAE E(x)-E(z) = {f2(soc.get('MAE_meV_cell_x_minus_z'), 3)} meV/f.u.")
    else:
        print("  not yet")
    if orbm:
        n = orbm.get("nscf") or {}
        print(f"  lorbm nscf: rc={n.get('rc')} done={n.get('job_done')} errors={n.get('errors')} exception={orbm.get('exception')}")
        for ln in (n.get("orbm_lines") or [])[:20]:
            print("    " + ln)
    out["soc_spin_residue"] = soc_res

    # ---------------------------------------------------------------- hydrated
    if h and out["pbeU"].get("KVCr", {}).get("U3_3"):
        a0 = out["pbeU"]["KVCr"]["U3_3"]["lcm"]; a1 = out["pbeU"]["KVCr_2H2O"]["U3_3"]["lcm"]
        if a0 and a1 and a0["winVB"] is not None and a1["winVB"] is not None:
            ch = [abs(a1[k] - a0[k]) / max(abs(a0[k]), 1e-3) for k in ("winVB", "winCB")]
            print(f"\n== hydration: KVCr.2H2O vs KVCr (U 3/3): winVB {a0['winVB']:.2f}->{a1['winVB']:.2f}, winCB {a0['winCB']:.2f}->"
                  f"{a1['winCB']:.2f} (rel. change {ch[0]:.0%}/{ch[1]:.0%}; kill 'existing sample' claim if > 50 %), gap "
                  f"{f2(a0['gap'])}->{f2(a1['gap'])}, edges {a1['vb']}/{a1['cb']}")
            out["hydration_rel_change"] = ch

    # ---------------------------------------------------------------- gate checks
    print("\n== JUDGE Gate-0 checks (None = not computable yet) ==")
    chk = {}
    for m in ("KVCr", "KVMo"):
        pts = [r["lcm"] for r in out["pbeU"].get(m, {}).values() if r.get("lcm")]
        oks = [unipolar_ok(p) for p in pts]
        chk[f"{m}: unipolar V-maj, both windows >= 0.3 eV at every PBE+U point ({len(pts)} pts)"] = \
            (all(oks) if oks and None not in oks else None)
        hs = (out["hse"].get(m) or {}).get("lcm")
        chk[f"{m}: HSE unipolar V-maj, windows >= 0.3 eV"] = unipolar_ok(hs)
        chk[f"{m}: HSE gap >= 1.0 eV"] = (hs["gap"] >= 1.0) if hs and hs["gap"] is not None else None
        lcm_ok = [p["M"] is not None and abs(p["M"]) <= 0.02 and (p["gap"] or 0) > 0.1 and not any("NOT-CONV" in f for f in p["flags"])
                  for p in pts + ([hs] if hs else [])]
        chk[f"{m}: LCM integer: M_cell = 0 and gapped at every PBE+U point and in HSE"] = all(lcm_ok) if lcm_ok else None
        fms = [r["fm"] for r in list(out["pbeU"].get(m, {}).values()) + ([out["hse"][m]] if out["hse"].get(m) else []) if r.get("fm")]
        chk[f"{m}: FM reference d3/d3 (M_FM = 6.00) at every point (else dE / J ratio contaminated)"] = \
            all(f["M"] is not None and abs(abs(f["M"]) - 6) <= 0.05 for f in fms) if fms else None
        dfl = [f for p in pts for f in p["flags"] if f.startswith("|m_d")]
        chk[f"{m}: Hubbard d moments 2.6-3.0 muB (covalency flag; not a kill if M_FM = 6)"] = (not dfl) if pts else None
        # kills
        chk[f"{m}: KILL HSE window < 0.2 eV"] = (min(hs["winVB"], hs["winCB"]) < 0.2) if hs and hs["winVB"] is not None else None
        chk[f"{m}: KILL HSE gap < 0.5 eV"] = (hs["gap"] < 0.5) if hs and hs["gap"] is not None else None
    cal = [q for k, l, q in rows if k == "VCr/CrCr" and l in ("PBE+U base (3/3 vs 3) [FM M=6]", "HSE06")]
    chk["calibration dE(VCr)/dE(CrCr) = 1.57 +- 0.25 (base PBE+U and HSE)"] = \
        (all(abs(q - 1.57) <= 0.25 for q in cal) if cal and None not in cal else None)
    vm = {l: q for k, l, q in rows if k == "VMo/VCr"}
    vb = vm.get("PBE+U base [FM M=6]")
    chk["dE(VMo)/dE(VCr) >= 1.15 in PBE+U (base, FM M=6) and HSE"] = \
        ((vb or 0) >= 1.15 and (vm.get("HSE06") or 0) >= 1.15) if vb and vm.get("HSE06") else None
    chk["KILL dE(VMo)/dE(VCr) < 1.0 (base PBE+U FM M=6, or HSE)"] = \
        any(q < 1.0 for q in (vb, vm.get("HSE06")) if q is not None) if (vb or vm.get("HSE06")) else None
    chk["SOC spin residue <= 0.1 muB/f.u. (orbital part: see lorbm lines / g-shift bound)"] = \
        (soc_res <= 0.1) if soc_res is not None else None
    chk["isomer >= +0.8 eV/f.u. (base U)"] = (iso["U3_3"] >= 0.8) if "U3_3" in iso else None
    if out.get("hydration_rel_change"):
        chk["hydrated KVCr windows change <= 50 %"] = max(out["hydration_rel_change"]) <= 0.5
    for k, v in chk.items():
        print(f"  [{'PASS' if v is True and 'KILL' not in k else ('KILL' if v is True else ('ok' if v is False and 'KILL' in k else ('FAIL' if v is False else '....')))}] {k}")
    out["checks"] = chk
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1, default=str))


if __name__ == "__main__":
    main()
