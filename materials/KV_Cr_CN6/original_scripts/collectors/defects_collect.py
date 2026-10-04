"""Collect Track L PBA defect Gate 1 (tracks/lcm/pba/defects/README.md; jobs lcm/pba/def/*).
usage: python tracks/lcm/pba/defects/collect.py [--fetch] [--npz] [--json out.json]
  --fetch : download the FIXED paths in manifest.json (results/<id>.json, status.json; refs.json) from the Modal volume
            into tracks/lcm/pba/defects/results/ (mirrors the volume path; never lists the volume). For jobs still in
            their relax, the relax .out is read at its fixed path and summarised (BFGS steps, E, M, force) - not kept.
  --npz   : also fetch results/<id>_bands.npz of finished jobs (1-5 MB each)
Per cell: M_cell (vs counting value), |M|, metal d moments, gap, edge spins (V-maj / M-maj), windows, and - after
aligning each cell to the pristine 4 f.u. cell with the far-field CN 3sigma level - the positions of the opposite-spin
(M-majority) states and of any in-gap states relative to the PRISTINE host edges, with their character.
Kill (JUDGE.md section 2, task spec): any opposite-spin state within 0.3 eV of a host edge (or inside the gap), or a
window cut by > 50 % (aligned window < 0.5 x pristine). Formation energies: E(AS) - E(P), E(FLIP) - E(P);
K vacancy E(KVAC) - E(P) + mu_K(bcc K); vacancy: E(VAC) - E(P) + 3 mu_K - 6 mu_H2O (+ mu[M(CN)6], not defined here)."""
import argparse
import json
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
CACHE = HERE / "results"
MAN = json.loads((HERE / "manifest.json").read_text())
JOBS = {j["id"]: j for j in MAN["jobs"]}
EDGE_TOL, WIN_FRAC = 0.3, 0.5
RY = 13.605693122994


def mrun():
    import importlib.util
    spec = importlib.util.spec_from_file_location("mrun", HERE.parents[3] / "infra" / "mrun.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def relax_peek(txt):
    en = [float(x) for x in re.findall(r"^!\s+total energy\s+=\s+([-0-9.]+)\s+Ry", txt, re.M)]
    tm = re.findall(r"total magnetization\s+=\s+([-0-9.]+)\s+Bohr", txt)
    am = re.findall(r"absolute magnetization\s+=\s+([-0-9.]+)\s+Bohr", txt)
    fz = re.findall(r"Total force =\s+([0-9.]+)", txt)
    nk = re.findall(r"number of k points=\s+(\d+)", txt)
    sym = re.findall(r"(\d+) Sym. Ops.|(No symmetry found)", txt)
    it = re.findall(r"convergence has been achieved in\s+(\d+)", txt)
    tt = re.findall(r"total cpu time spent up to now is\s+([0-9.]+)", txt)
    return {"n_scf": len(en), "E_first_eV": en[0] * RY if en else None, "E_last_eV": en[-1] * RY if en else None,
            "dE_relax_meV": (en[-1] - en[0]) * RY * 1000 if len(en) > 1 else None,
            "M": float(tm[-1]) if tm else None, "absM": float(am[-1]) if am else None, "force": float(fz[-1]) if fz else None,
            "nk": int(nk[0]) if nk else None, "symops": (sym[0][0] or sym[0][1]) if sym else None,
            "scf_iters": [int(x) for x in it][-5:], "cpu_s": float(tt[-1]) if tt else None,
            "bfgs_converged": "bfgs converged" in txt, "errors": re.findall(r"Error in routine.*\n.*", txt)[:2],
            "done": "JOB DONE" in txt}


def fetch(npz=False):
    m = mrun()
    got, miss, peeks = 0, [], {}
    for cid, j in JOBS.items():
        base = j["job"]
        paths = [f"jobs/{base}/status.json"] + [p for p in j["results"] if p.endswith(".json")]
        for p in paths:
            try:
                b = m._read(p)
            except Exception:
                miss.append(p); continue
            dst = CACHE / p; dst.parent.mkdir(parents=True, exist_ok=True); dst.write_bytes(b); got += 1
        res = CACHE / f"jobs/{base}/results/{cid}.json"
        done = False
        if res.exists():
            try:
                done = json.loads(res.read_text()).get("done", False)
            except Exception:
                pass
        if npz and done:
            for p in [p for p in j["results"] if p.endswith(".npz")]:
                dst = CACHE / p
                if dst.exists():
                    continue
                try:
                    dst.write_bytes(m._read(p)); got += 1
                except Exception:
                    miss.append(p)
        if cid != "refs" and not done:
            try:
                txt = m._read(f"jobs/{base}/{cid}_relax.out").decode(errors="replace")
                peeks[cid] = relax_peek(txt)
            except Exception:
                pass
    (CACHE / "relax_peek.json").write_text(json.dumps(peeks, indent=1))
    print(f"fetched {got} files; missing {len(miss)}")
    return peeks


def load(cid):
    j = JOBS[cid]
    p = CACHE / f"jobs/{j['job']}/results/{cid}.json"
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def status(cid):
    try:
        return json.loads((CACHE / f"jobs/{JOBS[cid]['job']}/status.json").read_text())
    except Exception:
        return {}


def vmaj_spin(r):
    """spin index of the V-majority channel = channel of the occupied V t2g (V moments negative -> 1)."""
    sm = (r.get("scf") or {}).get("site_moments") or []
    sp = (r.get("relaxed") or {}).get("species") or []
    v = [sm[i] for i, s in enumerate(sp) if s == "V" and i < len(sm)]
    if not v:
        return 1
    return 1 if sum(v) < 0 else 0


def lab(s, sv, M):
    return "V-maj" if s == sv else f"{M}-maj"


def metal_moments(r):
    sp = (r.get("relaxed") or {}).get("species") or []
    sm = (r.get("scf") or {}).get("site_moments") or []
    dm = (r.get("scf") or {}).get("d_moment") or {}
    role = None
    try:
        role = json.loads((HERE / "inputs" / f"{r['id']}.json").read_text())["role"]
    except Exception:
        pass
    out = []
    for i, s in enumerate(sp):
        if s in ("V", "Cr", "Mo"):
            out.append((i, s, role[i] if role else s, sm[i] if i < len(sm) else None, dm.get(str(i + 1), dm.get(i + 1))))
    return out


def ref(r):
    rf = r.get("refs") or {}
    return rf.get("E_ref_CN3s_far_mean") if rf.get("E_ref_CN3s_far_mean") is not None else rf.get("E_ref_CN3s_all_mean")


def analyse(cid, P, refs):
    r = load(cid)
    if not r or not r.get("done"):
        return None
    M = "Cr" if r["material"] == "KVCr" else "Mo"
    sv = vmaj_spin(r); so = 1 - sv
    e = r["edges"]; s = r["scf"]
    out = {"id": cid, "defect": r["defect"], "E_eV": s.get("energy_eV"), "M": s.get("total_mag"), "absM": s.get("abs_mag"),
           "M_expected": r.get("M_expected"), "gap": e["gap"], "vbm": lab(e["vbm_spin"], sv, M), "cbm": lab(e["cbm_spin"], sv, M),
           "win_VB_cell": e["win_VB"], "win_CB_cell": e["win_CB"], "sv": sv,
           "relax": {k: (r.get("relax") or {}).get(k) for k in ("bfgs_converged", "n_scf", "max_disp_A", "stopped_unconverged", "total_force")},
           "E_relax_vs_scf_meV": (((r.get("relax") or {}).get("energy_eV") or 0) - (s.get("energy_eV") or 0)) * 1000 if (r.get("relax") or {}).get("energy_eV") else None}
    out["metal_moments"] = [{"i": i, "el": el, "role": rl, "m_sphere": m_, "m_d": dm} for i, el, rl, m_, dm in metal_moments(r)]
    # partially occupied bands (holes / electrons in a metallic defect cell): spin, sum of (1 - occ), mean character
    hb = []
    for skey, rows in (r.get("band_table") or {}).items():
        for row in rows:
            if 0.02 < row["occ"] < 0.98:
                hb.append({"spin": lab(int(skey[1:]), sv, M), "band": row["band"], "occ": row["occ"],
                           "E_range": [row["Emin"], row["Emax"]], "groups": row["groups"], "top_atoms": row["top_atoms"][:3]})
    if hb:
        out["partial_bands"] = hb
        out["holes_per_spin"] = {sl: round(sum(1 - x["occ"] for x in hb if x["spin"] == sl), 3) for sl in {x["spin"] for x in hb}}
        g = {}
        for x in hb:
            for k, v in x["groups"].items():
                g[k] = g.get(k, 0) + (1 - x["occ"]) * v
        tot = sum(1 - x["occ"] for x in hb)
        out["hole_character"] = {k: round(v / tot, 3) for k, v in sorted(g.items(), key=lambda t: -t[1])} if tot > 0 else {}
    if P is None or not P.get("done"):
        return out
    eP = P["edges"]; svP = vmaj_spin(P)
    VBM_P, CBM_P = eP["top_occ_per_spin"][svP], eP["bot_emp_per_spin"][svP]
    out["host"] = {"VBM": VBM_P, "CBM": CBM_P, "gap": CBM_P - VBM_P,
                   "win_VB": VBM_P - eP["top_occ_per_spin"][1 - svP], "win_CB": eP["bot_emp_per_spin"][1 - svP] - CBM_P}
    rd, rp = ref(r), ref(P)
    if rd is None or rp is None:
        return out
    sh = rd - rp
    out["align_shift_eV"] = sh
    rf = r.get("refs") or {}
    if rf.get("E_ref_K3p_far_mean") is not None and (P.get("refs") or {}).get("E_ref_K3p_far_mean") is not None:
        out["align_shift_K3p_eV"] = rf["E_ref_K3p_far_mean"] - P["refs"]["E_ref_K3p_far_mean"]
    top, bot = e["top_occ_per_spin"], e["bot_emp_per_spin"]
    # aligned per-spin edges relative to the pristine host edges
    out["opp_top_rel_hostVBM"] = top[so] - sh - VBM_P          # must be <= -0.3
    out["opp_bot_rel_hostCBM"] = bot[so] - sh - CBM_P          # must be >= +0.3
    out["same_top_rel_hostVBM"] = top[sv] - sh - VBM_P         # > 0: occupied same-spin states above the host VBM
    out["same_bot_rel_hostCBM"] = bot[sv] - sh - CBM_P         # < 0: empty same-spin states below the host CBM
    out["win_VB_aligned"] = VBM_P - (top[so] - sh)
    out["win_CB_aligned"] = (bot[so] - sh) - CBM_P
    # in-gap / near-edge states from the band table (aligned)
    states = []
    for skey, rows in (r.get("band_table") or {}).items():
        sp = int(skey[1:])
        for row in rows:
            lo, hi = row["Emin"] - sh, row["Emax"] - sh
            ingap = hi > VBM_P + 0.05 and lo < CBM_P - 0.05
            opp_near = sp == so and (hi > VBM_P - EDGE_TOL and lo < CBM_P + EDGE_TOL)
            if ingap or opp_near:
                states.append({"spin": lab(sp, sv, M), "E_rel_hostVBM": [round(lo - VBM_P, 3), round(hi - VBM_P, 3)],
                               "E_rel_hostCBM": [round(lo - CBM_P, 3), round(hi - CBM_P, 3)], "occ": row["occ"],
                               "groups": row["groups"], "region": row["region"], "top_atoms": row["top_atoms"][:3]})
    out["near_gap_states"] = states
    # kill flags
    winP_VB, winP_CB = out["host"]["win_VB"], out["host"]["win_CB"]
    out["kill_opp_near_edge"] = bool(out["opp_top_rel_hostVBM"] > -EDGE_TOL or out["opp_bot_rel_hostCBM"] < EDGE_TOL
                                     or any(x["spin"] != "V-maj" for x in states))
    out["win_ratio_VB"] = out["win_VB_aligned"] / winP_VB if winP_VB > 0 else None
    out["win_ratio_CB"] = out["win_CB_aligned"] / winP_CB if winP_CB > 0 else None
    out["win_ratio_VB_cell"] = e["win_VB"] / winP_VB if (winP_VB > 0 and e["vbm_spin"] == sv) else None
    out["win_ratio_CB_cell"] = e["win_CB"] / winP_CB if (winP_CB > 0 and e["cbm_spin"] == sv) else None
    out["kill_window_cut"] = bool(min(out["win_ratio_VB"], out["win_ratio_CB"]) < WIN_FRAC)
    # formation energies
    EP, Ed = P["scf"]["energy_eV"], s["energy_eV"]
    d = r["defect"]
    muK = ((refs or {}).get("K_bcc") or {}).get("per_atom_eV")
    muW = ((refs or {}).get("H2O") or {}).get("energy_eV")
    if d in ("AS", "FLIP"):
        out["E_form_eV"] = Ed - EP
    elif d.startswith("KVAC") and muK is not None:
        out["E_form_eV"] = Ed - EP + muK
        out["E_form_def"] = "E(KVAC) - E(P) + mu_K(bcc K), neutral, K-rich limit"
    elif d == "VAC" and muK is not None and muW is not None:
        out["E_form_open_eV"] = Ed - EP + 3 * muK - 6 * muW
        out["E_form_def"] = "E(VAC) - E(P) + 3 mu_K(bcc) - 6 mu_H2O(gas); add mu[M(CN)6] (undefined here) for E_f"
    # hole / charge bookkeeping vs pristine (ortho-atomic d occupations, sphere moments)
    hub =(r.get("scf") or {}).get("hub") or {}
    hubP = (P.get("scf") or {}).get("hub") or {}
    sp = (r.get("relaxed") or {}).get("species") or []
    spP = (P.get("relaxed") or {}).get("species") or []
    for el in ("V", M):
        nd = [v[2] for k, v in hub.items() if int(k) - 1 < len(sp) and sp[int(k) - 1] == el]
        ndP = [v[2] for k, v in hubP.items() if int(k) - 1 < len(spP) and spP[int(k) - 1] == el]
        if nd and ndP:
            out[f"dN_d_{el}_total"] = sum(nd) - sum(ndP) * len(nd) / len(ndP)
    return out


def fmt(x, f="%.2f"):
    return "-" if x is None else (f % x)


def pct(x):
    return "-" if x is None else f"{100 * x:.0f}%"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--npz", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()
    peeks = fetch(a.npz) if a.fetch else (json.loads((CACHE / "relax_peek.json").read_text()) if (CACHE / "relax_peek.json").exists() else {})
    refs = None
    try:
        refs = json.loads((CACHE / f"jobs/{JOBS['refs']['job']}/results/refs.json").read_text())
    except Exception:
        pass
    print("== status ==")
    for cid in JOBS:
        st = status(cid); r = load(cid)
        stage = "-"
        if r:
            stage = "done" if r.get("done") else ("projwfc" if r.get("nscf") else ("nscf" if r.get("scf") else ("scf" if (r.get("relax") or {}).get("done") else "relax")))
        pk = peeks.get(cid)
        extra = ""
        if pk and stage in ("relax", "-"):
            extra = (f" relax: {pk['n_scf']} SCF, symops {pk['symops']}, nk {pk['nk']}, M {pk['M']}, |M| {pk['absM']}, F {pk['force']}, "
                     f"dE {fmt(pk['dE_relax_meV'], '%.0f')} meV, iters {pk['scf_iters']}, cpu {fmt(pk['cpu_s'], '%.0f')} s"
                     + (f" ERR {pk['errors']}" if pk["errors"] else ""))
        print(f"  {cid:16s} {st.get('state', '?'):8s} rc={st.get('rc')} h={fmt(st.get('hours'))} stage={stage}{extra}")
    if refs:
        print("== refs ==", {k: (v.get("per_atom_eV") or v.get("energy_eV")) for k, v in refs.items() if isinstance(v, dict)})
    allres = {}
    for m in ("KVCr", "KVMo"):
        Mel = "Cr" if m == "KVCr" else "Mo"
        P = load(f"{m}_P")
        print(f"\n== {m} (PBE+U V 3 / {Mel} {3 if m == 'KVCr' else 2}, 110 Ry, 4 f.u. cells; windows from the 4x4x4 nscf) ==")
        for cid in [c for c in JOBS if c.startswith(m + "_")]:
            o = analyse(cid, P, refs)
            if o is None:
                print(f"  {cid:16s} (not done)"); continue
            allres[cid] = o
            line = (f"  {cid:16s} M={fmt(o['M'], '%+.3f')} (exp {fmt(o['M_expected'], '%+.0f')}) |M|={fmt(o['absM'])} gap={fmt(o['gap'])} "
                    f"VBM:{o['vbm']} CBM:{o['cbm']} win(cell) {fmt(o['win_VB_cell'])}/{fmt(o['win_CB_cell'])}")
            if "align_shift_eV" in o:
                line += (f" | aligned to P (shift {o['align_shift_eV']:+.3f}): opp-spin top {o['opp_top_rel_hostVBM']:+.2f} vs VBM, "
                         f"opp-spin bottom {o['opp_bot_rel_hostCBM']:+.2f} vs CBM; win {o['win_VB_aligned']:.2f}/{o['win_CB_aligned']:.2f} "
                         f"({pct(o['win_ratio_VB'])}/{pct(o['win_ratio_CB'])} of P)")
            if "E_form_eV" in o:
                line += f" | E_f = {o['E_form_eV']:+.3f} eV"
            if "E_form_open_eV" in o:
                line += f" | E(VAC)-E(P)+3muK-6muH2O = {o['E_form_open_eV']:+.3f} eV (+ mu[M(CN)6])"
            print(line)
            if o.get("relax"):
                rl = o["relax"]
                print(f"      relax: conv={rl.get('bfgs_converged')} nSCF={rl.get('n_scf')} max_disp={fmt(rl.get('max_disp_A'), '%.3f')} A"
                      + (f" E(relax end)-E(scf)={o['E_relax_vs_scf_meV']:+.1f} meV" if o.get("E_relax_vs_scf_meV") is not None else ""))
            mm = o.get("metal_moments") or []
            if mm:
                print("      metals (role m_sphere / m_d): " + ", ".join(f"{x['role']}#{x['i']} {fmt(x['m_sphere'], '%+.2f')}/{fmt(x['m_d'], '%+.2f')}" for x in mm))
            if o.get("partial_bands"):
                print(f"      partially occupied bands: holes per spin {o['holes_per_spin']}, hole character {o['hole_character']}")
            for k in ("dN_d_V_total", f"dN_d_{Mel}_total"):
                if k in o:
                    print(f"      {k} vs pristine = {o[k]:+.3f} e/cell")
            for st_ in o.get("near_gap_states", []):
                print(f"      state {st_['spin']:7s} E-VBM_host {st_['E_rel_hostVBM']} E-CBM_host {st_['E_rel_hostCBM']} occ {st_['occ']} "
                      f"region {st_['region']} {st_['groups']} top {st_['top_atoms']}")
            if "kill_opp_near_edge" in o:
                print(f"      KILL opposite-spin within {EDGE_TOL} eV of an edge / in gap: {o['kill_opp_near_edge']}; "
                      f"KILL window cut > 50 %: {o['kill_window_cut']}")
        # K vacancy variants: relative energies
        ks = {c: allres[c]["E_eV"] for c in allres if c.startswith(m + "_KVAC")}
        if len(ks) > 1:
            b = min(ks.values())
            print("  K-vacancy variants (E - min, meV): " + ", ".join(f"{c.split('_', 1)[1]} {1000 * (v - b):.0f} (M {allres[c]['M']:+.2f})" for c, v in ks.items()))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(allres, indent=1))


if __name__ == "__main__":
    main()
