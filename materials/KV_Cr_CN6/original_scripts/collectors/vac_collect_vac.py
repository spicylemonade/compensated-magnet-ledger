"""Collect the [M(CN)6]-vacancy cells of the PBA defect gate (originals and the fast _r reruns).

usage: python tracks/lcm/pba/vac/collect_vac.py [--fetch] [--npz] [--json out.json]

Reuses the analysis of tracks/lcm/pba/defects/collect.py unchanged (alignment on the far-field CN 3sigma level of the
pristine 4 f.u. cell, kill flags: opposite-spin state within 0.3 eV of a host edge, window < 50 % of pristine), with
its own cache in tracks/lcm/pba/vac/results/ (fixed paths only, never lists the volume). Jobs:
  lcm/pba/def/{KVCr,KVMo}_P      pristine 4 f.u. references (done)
  lcm/pba/def/{KVCr,KVMo}_VAC    original relax -> scf -> nscf -> projwfc (slow relax, still running at 09:40 PDT)
  lcm/pba/def/{KVCr,KVMo}_VAC_r  SKIP_RELAX rerun at the last relax positions (vac/manifest.json)
  lcm/pba/def/KVCr_VAC_w1        KVCr_VAC_r + one zeolitic H2O at the vacancy centre (P1, short relax; vac/manifest.json)
  lcm/pba/def/refs               mu_K (bcc K), mu_H2O
Writes results/vac_summary.json (per cell: M, gap, edge spins, windows, aligned windows, opposite-spin edges vs host
edges, near-gap states with character, capped-V vs bulk-V moments, kill flags).
"""
import argparse
import importlib.util
import json
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DEF = HERE.parent / "defects"
_spec = importlib.util.spec_from_file_location("defcol", DEF / "collect.py")
col = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(col)

VMAN = json.loads((HERE / "manifest.json").read_text())
KEEP = ["KVCr_P", "KVMo_P", "KVCr_VAC", "KVMo_VAC", "refs"]
JOBS = {k: col.JOBS[k] for k in KEEP}
for j in VMAN["jobs"]:
    JOBS[j["id"]] = j
col.JOBS = JOBS
col.CACHE = HERE / "results"
col.CACHE.mkdir(exist_ok=True)

_mm = col.metal_moments


def metal_moments(r):
    """role lookup in defects/inputs or vac/inputs (the _r cells)"""
    sp = (r.get("relaxed") or {}).get("species") or []
    sm = (r.get("scf") or {}).get("site_moments") or []
    dm = (r.get("scf") or {}).get("d_moment") or {}
    role = None
    for d in (DEF / "inputs", HERE / "inputs"):
        p = d / f"{r['id']}.json"
        if p.exists():
            role = json.loads(p.read_text())["role"]
            break
    out = []
    for i, s in enumerate(sp):
        if s in ("V", "Cr", "Mo"):
            out.append((i, s, role[i] if role else s, sm[i] if i < len(sm) else None, dm.get(str(i + 1), dm.get(i + 1))))
    return out


col.metal_moments = metal_moments


def capped(cid):
    """indices of the water-capped V (from the region list of the input)"""
    for d in (DEF / "inputs", HERE / "inputs"):
        p = d / f"{cid}.json"
        if p.exists():
            c = json.loads(p.read_text())
            return set(i for i in (c.get("region") or []) if c["species"][i] == "V")
    return set()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    ap.add_argument("--npz", action="store_true")
    ap.add_argument("--json")
    a = ap.parse_args()
    peeks = col.fetch(a.npz) if a.fetch else (json.loads((col.CACHE / "relax_peek.json").read_text())
                                               if (col.CACHE / "relax_peek.json").exists() else {})
    refs = None
    try:
        refs = json.loads((col.CACHE / f"jobs/{JOBS['refs']['job']}/results/refs.json").read_text())
    except Exception:
        pass
    print("== status ==")
    for cid in JOBS:
        st = col.status(cid); r = col.load(cid)
        stage = "-"
        if r:
            stage = "done" if r.get("done") else ("projwfc" if r.get("nscf") else ("nscf" if r.get("scf") else
                                                  ("scf" if (r.get("relax") or {}).get("done") else "relax")))
        pk = peeks.get(cid)
        extra = ""
        if pk and stage in ("relax", "-"):
            extra = f" relax: {pk['n_scf']} SCF, M {pk['M']}, F {pk['force']}, dE {col.fmt(pk['dE_relax_meV'], '%.0f')} meV"
        print(f"  {cid:12s} {st.get('state', '?'):8s} rc={st.get('rc')} h={col.fmt(st.get('hours'))} stage={stage}{extra}")
    summ = {}
    for m in ("KVCr", "KVMo"):
        P = col.load(f"{m}_P")
        for cid in [c for c in JOBS if c.startswith(f"{m}_VAC")]:
            o = col.analyse(cid, P, refs)
            if o is None:
                print(f"\n  {cid}: not done"); continue
            cv = capped(cid)
            mm = o.get("metal_moments") or []
            o["capped_V"] = [x for x in mm if x["i"] in cv]
            o["bulk_V"] = [x for x in mm if x["el"] == "V" and x["i"] not in cv]
            summ[cid] = o
            print(f"\n  {cid}: M = {o['M']:+.3f} (expected {o['M_expected']:+.0f}), |M| {o['absM']:.2f}, gap {o['gap']:.2f} eV, "
                  f"VBM {o['vbm']} CBM {o['cbm']}, cell windows {o['win_VB_cell']:.2f}/{o['win_CB_cell']:.2f}")
            if "align_shift_eV" in o:
                print(f"     aligned (shift {o['align_shift_eV']:+.3f}; K3p {col.fmt(o.get('align_shift_K3p_eV'), '%+.3f')}): "
                      f"opp-spin VB top {o['opp_top_rel_hostVBM']:+.2f} vs host VBM, opp-spin CB bottom {o['opp_bot_rel_hostCBM']:+.2f} "
                      f"vs host CBM; same-spin top {o['same_top_rel_hostVBM']:+.2f}, bottom {o['same_bot_rel_hostCBM']:+.2f}; "
                      f"windows {o['win_VB_aligned']:.2f}/{o['win_CB_aligned']:.2f} ({col.pct(o['win_ratio_VB'])}/{col.pct(o['win_ratio_CB'])})")
            print("     capped V (m_sphere/m_d): " + ", ".join(f"#{x['i']} {col.fmt(x['m_sphere'], '%+.2f')}/{col.fmt(x['m_d'], '%+.2f')}" for x in o["capped_V"])
                  + " | bulk V: " + ", ".join(f"#{x['i']} {col.fmt(x['m_sphere'], '%+.2f')}/{col.fmt(x['m_d'], '%+.2f')}" for x in o["bulk_V"]))
            print("     M sites: " + ", ".join(f"#{x['i']} {col.fmt(x['m_sphere'], '%+.2f')}/{col.fmt(x['m_d'], '%+.2f')}" for x in mm if x["el"] != "V"))
            for s in o.get("near_gap_states", []):
                print(f"     state {s['spin']:7s} E-VBM {s['E_rel_hostVBM']} E-CBM {s['E_rel_hostCBM']} occ {s['occ']} region {s['region']} "
                      f"{s['groups']} top {s['top_atoms']}")
            if "E_form_open_eV" in o:
                print(f"     E(VAC) - E(P) + 3 mu_K - 6 mu_H2O = {o['E_form_open_eV']:+.3f} eV (+ mu[M(CN)6])")
            if "kill_opp_near_edge" in o:
                print(f"     KILL opp-spin near edge: {o['kill_opp_near_edge']}; KILL window cut: {o['kill_window_cut']}")
    (col.CACHE / "vac_summary.json").write_text(json.dumps(summ, indent=1))
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
