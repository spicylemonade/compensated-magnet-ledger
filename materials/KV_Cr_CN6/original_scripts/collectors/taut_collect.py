"""Collect Track L PBA "taut" (electronic competitors of the d3/d3 LCM). See README.md.
usage: python tracks/lcm/pba/taut/collect.py [--fetch] [--peek JOB FILE]
  --fetch : download the FIXED result paths (jobs/lcm/pba/taut/<job>/results/<job>.json, status.json) into results/
            (mirrors the volume path; no listdir). Missing files are reported, never fatal.
  --peek JOB FILE : print the tail of a file in a job dir (e.g. a running pw.x .out)
Prints, per job: every start (vertical, at the stage-1 LCM cell, nosym 4x4x4, 110 Ry): dE vs LCM (eV/f.u.), M, |M|,
d moments / d counts of V (N-bound) and M (C-bound), state label, gap and spin windows; the relaxed lowest competitor;
HSE06 single points; K variants. Writes summary.json."""
import argparse, json, pathlib, importlib.util

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[3]
CACHE = HERE / "results"
P = "jobs/lcm/pba/taut"
TAUT = ["KVMo_U3_2_CT", "KVMo_U3_2_LS", "KVCr_U3_3_CT", "KVCr_U3_3_LS", "KVMo_U4_0", "KVMo_U2_0", "KVCr_U4_2", "KVCr_U2_2"]
KVAR = ["KVMo_Kvar", "KVCr_Kvar"]


def mrun():
    spec = importlib.util.spec_from_file_location("mrun", ROOT / "infra" / "mrun.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def fetch():
    m = mrun(); miss = []
    for j in TAUT + KVAR:
        for p in (f"{P}/{j}/results/{j}.json", f"{P}/{j}/status.json"):
            try:
                b = m._read(p)
            except Exception:
                miss.append(p); continue
            dst = CACHE / p; dst.parent.mkdir(parents=True, exist_ok=True); dst.write_bytes(b)
    print(f"missing {len(miss)}: {miss[:6]}")


def load(j, f="results"):
    p = CACHE / (f"{P}/{j}/results/{j}.json" if f == "results" else f"{P}/{j}/status.json")
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


NS_FREE = {"LCM", "CT_AF_free", "FM6"}  # runs whose v1 result is valid (no kick, or kick irrelevant: FM6 = Gate-0 FM)


def classify(q, ref):
    """same rule as jobsrc/taut.py: CT if dN_d(V) < -0.15; LCM-like if both moments within 0.25 and |dN_d(V)| < 0.1"""
    try:
        sN, sC, rN, rC = q["sites"]["N"], q["sites"]["C"], ref["sites"]["N"], ref["sites"]["C"]
        dNN = sN["N_d"] - rN["N_d"] if "N_d" in sN else None
        mN = sN.get("m_d", sN.get("m_sph")); mC = sC.get("m_d", sC.get("m_sph"))
        rmN = rN.get("m_d", rN.get("m_sph")); rmC = rC.get("m_d", rC.get("m_sph"))
        ct = dNN is not None and dNN < -0.15
        same = abs(mN - rmN) < 0.25 and abs(mC - rmC) < 0.25 and (dNN is None or abs(dNN) < 0.1)
        return {"ct": bool(ct), "dN_d_N": dNN, "lcm_like": bool(same), "label": f"{'CT' if ct else 'noCT'} V:{mN:+.2f} M:{mC:+.2f}"}
    except Exception as e:
        return {"err": repr(e)}


def fmt_sites(s):
    o = []
    for r, el in (("N", "V"), ("C", "M")):
        d = (s or {}).get(r, {})
        md = d.get("m_d"); nd = d.get("N_d"); ms = d.get("m_sph")
        o.append(f"{el} m_d={md:+.2f} N_d={nd:.2f}" if md is not None else f"{el} m_sph={ms:+.2f}" if ms is not None else f"{el} -")
    return " | ".join(o)


def orb(s, r, spin):
    d = (s or {}).get(r, {}).get("diag") or {}
    v = d.get(str(spin)) or d.get(spin)
    return "[" + " ".join(f"{x:.2f}" for x in v) + "]" if v else "-"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--peek", nargs=2)
    a = ap.parse_args()
    if a.peek:
        t = mrun()._read(f"{P}/{a.peek[0]}/{a.peek[1]}").decode(errors="replace")
        print(t[-4000:]); return
    if a.fetch:
        fetch()
    summ = {"vertical": {}, "relaxed": {}, "hse": {}, "kvar": {}, "status": {}}
    for j in TAUT + KVAR:
        st = load(j, "status") or {}
        summ["status"][j] = {k: st.get(k) for k in ("state", "rc", "hours")}
    print("== status ==")
    for j, s in summ["status"].items():
        print(f"  {j:16s} {s.get('state')} rc={s.get('rc')} h={s.get('hours') and round(s['hours'], 2)}")
    print("\n== vertical constrained starts at the stage-1 LCM cell (PBE+U, nosym 4x4x4, 110 Ry); dE = E - E(LCM) per f.u. ==")
    print("   orbital diag order (QE real harmonics): z2 xz yz x2-y2 xy; spin 1 = up (M-majority in the AF starts), 2 = down")
    for j in TAUT:
        d = load(j)
        if not d:
            print(f"  {j}: no results yet"); continue
        print(f"  -- {j}  U={d.get('U')}")
        ref = d.get("runs", {}).get("LCM")
        for k, q in d.get("runs", {}).items():
            if q.get("libver") != 2 and k not in NS_FREE:
                print(f"    {k:11s} (v1 run with the wrong ns kick; superseded when the v3 job reaches it)"); continue
            g = q.get("glob") or {}
            st = classify(q, ref) if ref else {}
            e = (q["energy_eV"] - ref["energy_eV"]) if (ref and q.get("energy_eV") is not None) else None
            print(f"    {k:11s} dE={'%+.3f' % e if e is not None else '   -  '} M={q.get('total_mag')} |M|={q.get('abs_mag')} conv={q.get('converged_scf')} "
                  f"{fmt_sites(q.get('sites'))} [{st.get('label', '')}{' LCM-like' if st.get('lcm_like') else ''}] "
                  f"gap={g.get('gap', float('nan')):.2f} vbm_s={g.get('vbm_spin')} cbm_s={g.get('cbm_spin')} "
                  f"win={g.get('win_VB', float('nan')):.2f}/{g.get('win_CB', float('nan')):.2f} nfrac={g.get('n_frac_occ')} "
                  f"nks={q.get('nks')} {q.get('sec') and round(q['sec'])}s")
            if k != "LCM" and q.get("sites"):
                print(f"                V(dn,up) {orb(q['sites'], 'N', 2)} {orb(q['sites'], 'N', 1)}  M(up,dn) {orb(q['sites'], 'C', 1)} {orb(q['sites'], 'C', 2)}")
            summ["vertical"].setdefault(j, {})[k] = {"dE_eV": e, "M": q.get("total_mag"), "absM": q.get("abs_mag"), "conv": q.get("converged_scf"),
                                                    "sites": q.get("sites"), "state": st, "gap": g.get("gap"), "win_VB": g.get("win_VB"),
                                                    "win_CB": g.get("win_CB"), "vbm_spin": g.get("vbm_spin"), "cbm_spin": g.get("cbm_spin")}
        for cls, rec in (d.get("relax") or {}).items():
            q = rec.get("run") or {}
            g = q.get("glob") or {}
            print(f"    RELAX[{cls}] target={rec.get('target')} cands={rec.get('candidates')} {rec.get('note', '')}")
            if q:
                print(f"       dE_relaxed={q.get('dE_eV')} gain={q.get('E_relax_gain_eV')} M={q.get('total_mag')} {fmt_sites(q.get('sites'))} "
                      f"[{(q.get('state') or {}).get('label')}] bfgs={q.get('bfgs_converged')} max_disp={q.get('max_disp_A')} gap={g.get('gap')} "
                      f"win={g.get('win_VB')}/{g.get('win_CB')}")
                print(f"       bonds={q.get('bonds')}")
            summ["relaxed"].setdefault(j, {})[cls] = {"target": rec.get("target"), "candidates": rec.get("candidates"),
                                                     "dE_eV": q.get("dE_eV"), "gain_eV": q.get("E_relax_gain_eV"), "M": q.get("total_mag"),
                                                     "sites": q.get("sites"), "state": q.get("state"), "gap": g.get("gap"),
                                                     "bonds": q.get("bonds"), "bfgs": q.get("bfgs_converged"), "note": rec.get("note")}
        for cls, rec in (d.get("hse") or {}).items():
            q = rec.get("run") or {}
            g = q.get("glob") or {}
            preps = {k: (v.get("nsym"), v.get("nks"), fmt_sites(v.get("sites"))) for k, v in rec.items() if k.startswith("prep")}
            print(f"    HSE[{cls}] target={rec.get('target')} spg={rec.get('spacegroups')} preps={preps} {rec.get('note', '')}")
            if q:
                print(f"       dE_HSE={q.get('dE_eV')} M={q.get('total_mag')} |M|={q.get('abs_mag')} m_sph={q.get('m_sph')} (prep {q.get('m_sph_prep')}) "
                      f"conv={q.get('converged_scf')} gap={g.get('gap')} win={g.get('win_VB')}/{g.get('win_CB')} vbm_s={g.get('vbm_spin')} cbm_s={g.get('cbm_spin')} "
                      f"nsym={q.get('nsym')} nks={q.get('nks')}")
            summ["hse"].setdefault(j, {})[cls] = {"target": rec.get("target"), "spacegroups": rec.get("spacegroups"), "dE_eV": q.get("dE_eV"),
                                                 "M": q.get("total_mag"), "absM": q.get("abs_mag"), "m_sph": q.get("m_sph"),
                                                 "m_sph_prep": q.get("m_sph_prep"), "conv": q.get("converged_scf"), "glob": g,
                                                 "note": rec.get("note")}
    print("\n== K variants (LCM, PBE+U base U, 110 Ry, symmetry on, ionic relax at fixed cell) ==")
    for j in KVAR:
        d = load(j)
        if not d:
            print(f"  {j}: no results yet"); continue
        for k, q in d.get("runs", {}).items():
            g = q.get("glob") or {}
            print(f"  {j} {k:5s} dE/f.u.={q.get('dE_eV_per_fu')} M={q.get('total_mag')} bfgs={q.get('bfgs_converged')} nsym={q.get('nsym')} "
                  f"K_offset={q.get('K_offset_final_A')} gap={g.get('gap')} win={g.get('win_VB')}/{g.get('win_CB')}")
            summ["kvar"].setdefault(j, {})[k] = {"dE_eV_per_fu": q.get("dE_eV_per_fu"), "M": q.get("total_mag"), "bfgs": q.get("bfgs_converged"),
                                                "K_offset_A": q.get("K_offset_final_A"), "glob": g}
    (HERE / "results" / "collect_snapshot.json").write_text(json.dumps(summ, indent=1, default=str))
    verdict(summ)


# HSE06 LCM sphere moments at the cubic cell (deep_hse.json): V(N) / M(C)
HSE_LCM_SPH = {"KVCr": (-2.01, 2.33), "KVMo": (-1.93, 1.43)}


def verdict(summ):
    """per compound: lowest competitor and lowest CT state per U point (vertical, relaxed), HSE status, kill flag"""
    V = {}
    for j, runs in summ["vertical"].items():
        m = j.split("_")[0]; u = "_".join(j.split("_")[1:3])
        row = V.setdefault(m, {}).setdefault(u, {"lowest": None, "lowest_CT": None, "n_starts": 0, "relaxed_CT": None, "relaxed_LS": None})
        for k, q in runs.items():
            if k == "LCM" or q.get("dE_eV") is None or not q.get("conv"):
                continue
            st = q.get("state") or {}
            row["n_starts"] += 1
            if st.get("lcm_like"):
                continue
            ent = {"start": k, "dE_eV": round(q["dE_eV"], 4), "M": q.get("M"), "label": st.get("label"), "gap": q.get("gap"),
                   "win_VB": q.get("win_VB"), "win_CB": q.get("win_CB")}
            if row["lowest"] is None or q["dE_eV"] < row["lowest"]["dE_eV"]:
                row["lowest"] = ent
            if st.get("ct") and (row["lowest_CT"] is None or q["dE_eV"] < row["lowest_CT"]["dE_eV"]):
                row["lowest_CT"] = ent
    for j, rel in summ["relaxed"].items():
        m = j.split("_")[0]; u = "_".join(j.split("_")[1:3])
        for cls, r in rel.items():
            if r.get("dE_eV") is not None:
                V[m][u][f"relaxed_{cls}"] = {"target": r.get("target"), "dE_eV": round(r["dE_eV"], 4), "gain_eV": r.get("gain_eV"),
                                            "label": (r.get("state") or {}).get("label"), "lcm_like": (r.get("state") or {}).get("lcm_like"),
                                            "gap": r.get("gap")}
    H = {}
    for j, hh in summ["hse"].items():
        m = j.split("_")[0]
        for cls, h in hh.items():
            ent = {"target": h.get("target"), "dE_eV": h.get("dE_eV"), "M": h.get("M"), "m_sph": h.get("m_sph"), "m_sph_prep": h.get("m_sph_prep"),
                   "gap": (h.get("glob") or {}).get("gap"), "conv": h.get("conv"), "note": h.get("note"), "spacegroups": h.get("spacegroups")}
            ms = h.get("m_sph") or {}
            if ms.get("N") is not None:
                ref = HSE_LCM_SPH[m]
                ent["collapsed_to_LCM"] = bool(abs(ms["N"] - ref[0]) < 0.12 and abs(ms["C"] - ref[1]) < 0.12 and abs(h.get("M") or 0) < 0.05)
            H.setdefault(m, {})[cls] = ent
    kill = {}
    for m, hh in H.items():
        c = hh.get("CT") or {}
        kill[m] = (None if c.get("dE_eV") is None else bool(c["dE_eV"] < 0 and not c.get("collapsed_to_LCM")))
    out = {"what": "PBA taut: electronic competitors of the d3/d3 LCM (tracks/lcm/pba/taut/README.md)",
           "units": "dE per f.u. (15-atom cell) in eV, relative to the LCM at identical settings (PBE+U: nosym 4x4x4 110 Ry; "
                    "HSE: existing HSE06 LCM at the cubic cell, 90 Ry)",
           "pbeu": V, "hse": H, "kill_CT_below_LCM_in_HSE": kill, "kvar": summ["kvar"], "status": summ["status"],
           "vertical_all": summ["vertical"], "relaxed_all": summ["relaxed"],
           "notes": ["pbeu 'lowest' includes FM6 (the d3/d3 FM exchange reference, M fixed 6)",
                     "kicked CT_AF / CT_AF_v starts that converge LCM-like are excluded from 'lowest' (lcm_like = True in vertical_all)",
                     "kill_CT_below_LCM_in_HSE is null until the HSE06 chains finish (see README section 0 'Pending')",
                     "v1 runs with the wrong ns kick (libver != 2, kicked starts) are excluded"],
           "jobs": json.loads((HERE / "manifest.json").read_text()).get("jobs") if (HERE / "manifest.json").exists() else None}
    (HERE / "summary.json").write_text(json.dumps(out, indent=1, default=str))
    print("\n== verdict ==")
    for m, us in V.items():
        for u, r in us.items():
            lo, ct = r.get("lowest") or {}, r.get("lowest_CT") or {}
            print(f"  {m} {u}: lowest competitor {lo.get('start')} {lo.get('dE_eV')} eV [{lo.get('label')}]; lowest CT {ct.get('start')} "
                  f"{ct.get('dE_eV')} eV; relaxed CT {r.get('relaxed_CT')}; relaxed LS {r.get('relaxed_LS')}")
    for m, hh in H.items():
        print(f"  {m} HSE: {hh}")
    print(f"  kill (CT below LCM in HSE): {kill}")


if __name__ == "__main__":
    main()
