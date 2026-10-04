"""Track L PBA T_C model: Heisenberg model of the d3/d3 cyanide double perovskites, calibrated on measured T_C.

Inputs (all read, nothing submitted):
  * collect.py JSON snapshot (E(FM) - E(LCM) per magnetic ion at every PBE+U point; HSE when finished; FM = fixed-M
    d3/d3 reference where available)  -> tc/collect_snapshot.json  (--fetch reruns collect.py --fetch)
  * 2-f.u. supercell energies (jexch_submit.py; jobs lcm/pba/tc/J2_*) -> J1, J2C, J2N, J3, J4C, J4N at base U
  * classical MC tables (tc_mc_runs.py) in tc/mc/*.json: pure-model T_C, vacancy dilution T_C(p), NNN runs
Model: H = sum_<ij> J_ij S_i.S_j (J > 0 AF), S = 3/2 on both sublattices, simple-cubic metal net (rock-salt colouring).
  Broken-symmetry (Ising) mapping of collinear DFT states: E(FM) - E(LCM) = (6 J1 + 8 J3) S^2 per magnetic ion.
Temperature levels (all in K for the DFT J):
  MF   : two-sublattice molecular field, k T = S(S+1)/3 [ (aA+aB)/2 + sqrt(((aA-aB)/2)^2 + b^2) ],
         b = 6 J1 + 8 J3, aX = -(12 J2X + 6 J4X)
  MC   : classical spins of length S (S^2 scale; Schart 2024 convention), Binder crossings; pure J1: 1.4430 J1 S^2
  RPA  : Tyablikov, quantum S = 3/2: k T = S(S+1)/3 / < A/(A^2 - B^2) >_BZ  (A = J_AB(0) - J_AA(0) + J_AA(q),
         B = J_AB(q); intra couplings averaged over the two sublattices); pure J1: 1.3189 J1 S(S+1)
  QMC~ : quantum-corrected classical, T_MC * S(S+1)/S^2 * theta_Q(3/2)/theta_cl, theta from the S = 1/2 sc-AF QMC value
         (k T_N = 0.946 J, Sandvik PRL 80, 5196 (1998)) and the classical limit 1.443, interpolated in 1/S(S+1) or 1/S
Calibration: f = T_exp / T_model(sample composition). Prediction: T = f * T_model(target composition).
usage: python tc_model.py [--fetch]   -> tc/summary.json (+ table on stdout)"""
import argparse
import importlib.util
import json
import math
import pathlib
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
PBA = HERE.parent
ROOT = PBA.parents[2]
MC = HERE / "mc"
KB = 0.0861733  # meV/K
S = 1.5
SS1, S2 = S * (S + 1), S * S
TH_CL = 1.4430                     # classical sc Heisenberg (unit spins)
TH_Q12 = 0.946 / 0.75              # S = 1/2 sc AF QMC, per S(S+1)
WATSON = 1.516386                  # sc Watson integral <1/(1-gamma)>

# ---------------------------------------------------------------- experimental anchors (checked 2026-10-03)
EXP = {
    "CrCr": {"compound": "Cr(III)[Cr(III)(CN)6]0.996(7) (NMF route, Fm-3m, a = 10.42 A)",
             "ref": "Schart, Torres-Cavanillas, Wheeler, ..., Pasta, Inorg. Chem. 2024, doi 10.1021/acs.inorgchem.4c03856 (PMC11615938)",
             "T": 240.0, "T_err": 10.0, "how": "neutron (111) magnetic intensity vanishes at 240 +- 10 K; chi*T minimum 265 K",
             "p": 0.996, "p_err": 0.007, "w": 0.0,
             "note": "NOT the classic Cr(II)[Cr(III)(CN)6]2/3.10/3H2O (Mallah 1993, also 240 K): that one is HS Cr(II) S = 2 "
                     "on N with 1/3 vacancies, not d3/d3. Schart: HSE06 J/kB = -100 K (H = -J sum_<ij>), classical MC "
                     "T_C 323 K (S^2 convention), Theta = -836 K, M_nc 0.038 muB/f.u."},
    "KVCr": {"compound": "K V(II)[Cr(III)(CN)6].2H2O.0.1 KO3SCF3 (sol-gel, crystalline, a = 10.55 A)",
             "ref": "Holmes & Girolami, JACS 121, 5593 (1999); Verdaguer & Girolami review (Magnetism: Molecules to Materials V, 2005) Table 9.8",
             "T": 376.0, "T_err": 11.0, "how": "M(T); 365 K after heat treatment (spread used as error)",
             "p": 0.979, "p_err": 0.021, "w": 0.0,
             "note": "z = 1 by analysis, M_sat 0.7 kG cm3/mol = 0.125 muB/f.u. = |w - 3v| (w = V(III), v = [Cr(CN)6] "
                     "vacancy): bracket v in [0, 0.042] (all-vacancy limit) -> p = 0.979 +- 0.021; hydrated (zeolitic 2 H2O)"},
}
CHECKS = {
    "VCr_Ferlay": {"compound": "V(II)0.42V(III)0.58[Cr(III)(CN)6]0.86.2.8H2O (amorphous)", "ref": "Ferlay et al., Nature 378, 701 (1995)",
                   "T": 315.0, "p": 0.86, "w": 0.58, "S_w": 1.0,
                   "rJ_lit": 31.24 / 75.56,
                   "rJ_note": "J(V(III)-Cr)/J(V(II)-Cr) = -31.24/-75.56 cm-1 (cluster DFT, Table 9.3 of the Verdaguer-Girolami review, ref. 103)"},
    "VCr_23": {"compound": "V(II)[Cr(III)(CN)6]2/3.3.5H2O (crystalline, a = 10.54 A; same paper as KVCr)",
               "ref": "Holmes & Girolami 1999 (review Table 9.8)", "T": 330.0, "p": 2 / 3, "w": 0.0},
    "VCr_Cs": {"compound": "Cs0.82V(II)[Cr(III)(CN)6]0.92-0.94 (crystalline, a = 10.65 A)", "ref": "Holmes & Girolami 1999",
               "T": 337.0, "p": 0.93, "w": 0.0},
    "VMo_Magott": {"compound": "{[K(crypt-222)]0.34 V(II)1.37[Mo(III)(CN)6](BF4)0.08}.xMeCN (amorphous, 27 % [Mo(CN)6] vacancies)",
                   "ref": "Magott et al., Adv. Sci. 2025/26, doi 10.1002/advs.202511285 (PMC12850227)",
                   "T": 417.0, "T_lo": 340.0, "T_hi": 454.0, "p": 1 / 1.37, "w": 0.0,
                   "how": "M(T) measured 2-340 K only (decomposes above 340 K); T_c extrapolated: 417 K (Bloch exponent "
                          "1.927 fit), 454 K (Bloch 3/2); hysteresis up to 340 K -> hard lower bound 340 K"},
}


# ---------------------------------------------------------------- data
def run_collect(fetch):
    snap = HERE / "collect_snapshot.json"
    cmd = [sys.executable, str(PBA / "collect.py"), "--json", str(snap)] + (["--fetch"] if fetch else [])
    subprocess.run(cmd, cwd=str(ROOT), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
    return json.loads(snap.read_text())


def dE_table(snap):
    """meV per magnetic ion; FM = fixed-M (d3/d3) where present"""
    out = {}
    for m in ("KVCr", "KVMo", "CrCr", "KVCr_2H2O"):
        for tag, r in (snap["pbeU"].get(m) or {}).items():
            if not r.get("lcm"):
                continue
            v = r.get("dEfix") if r.get("dEfix") is not None else r.get("dE")
            if v is None:
                continue
            flag = []
            fm = r.get("fm") or {}
            if fm.get("M") is not None and abs(abs(fm["M"]) - 6) > 0.05:
                flag.append(f"free FM charge-transferred (M = {fm['M']:.2f}); fixed-M value used")
            if (r.get("fmfix") or {}).get("gap") is not None and r["fmfix"]["gap"] < 0.05:
                flag.append(f"fixed-M FM gap {r['fmfix']['gap']:.2f} eV (d3/d3 FM unstable to V->M transfer)")
            out.setdefault(m, {})[tag] = {"dE": float(v), "flags": flag}
        h = (snap.get("hse") or {}).get(m)
        if h and h.get("dE") is not None:
            out.setdefault(m, {})["HSE"] = {"dE": float(h["dE"]), "flags": []}
    return out


BASE = {"KVCr": "U3_3", "KVMo": "U3_2", "CrCr": "U3"}
RY = 13.605693122994


def hse_prelim(fetch):
    """Preliminary HSE06 E(FMfix) - E(LCM) for KVCr / KVMo from the RUNNING jobs' pw.x outputs (fixed paths): last
    '!' energy of each, accepted only after >= 2 EXX outer loops with dexx < 1e-4 Ry in BOTH runs. Cached in
    tc/results/hse_prelim.json. Superseded by collect.py's deep_hse.json + FMfix results once those exist."""
    import re
    cache = HERE / "results" / "hse_prelim.json"
    old = json.loads(cache.read_text()) if cache.exists() else {}
    if fetch:
        spec = importlib.util.spec_from_file_location("mrun", ROOT / "infra" / "mrun.py")
        mrun = importlib.util.module_from_spec(spec); spec.loader.exec_module(mrun)
        for m in ("KVCr", "KVMo"):
            rec = {}
            for o, path in (("LCM", f"jobs/lcm/pba/hse/{m}_hse_LCM/hse_lcm.out"),
                            ("FMfix", f"jobs/lcm/pba/hse/{m}_hse_FMfix/{m}_hse_FMfix.out")):
                try:
                    t = mrun._read(path).decode(errors="replace")
                except Exception:
                    continue
                e = re.findall(r"^!+\s+total energy\s+=\s+([-0-9.]+)\s+Ry", t, re.M)
                dx = re.findall(r"est\. exchange err \(dexx\)\s+=\s+([-0-9.Ee+]+)", t)
                tm = re.findall(r"total magnetization\s+=\s+([-0-9.]+)", t)
                rec[o] = {"E_Ry": float(e[-1]) if e else None, "E_prev_Ry": float(e[-2]) if len(e) > 1 else None,
                          "n_bang": len(e), "dexx": float(dx[-1]) if dx else None,
                          "M": float(tm[-1]) if tm else None, "done": "JOB DONE" in t, "path": path}
            if rec:
                old[m] = rec
        old["t"] = __import__("time").strftime("%Y-%m-%d %H:%M")
        cache.parent.mkdir(parents=True, exist_ok=True); cache.write_text(json.dumps(old, indent=1))
    out = {}
    for m in ("KVCr", "KVMo"):
        r = old.get(m) or {}
        L, F = r.get("LCM"), r.get("FMfix")
        if not (L and F) or None in (L["E_Ry"], F["E_Ry"], L["dexx"], F["dexx"]):
            continue
        if L["n_bang"] < 3 or F["n_bang"] < 3 or L["dexx"] > 2e-3 or F["dexx"] > 2e-3:
            continue
        # convergence error: the last EXX outer-loop energy change of each run (changes shrink ~10x per loop)
        err = sum(abs(x["E_Ry"] - x["E_prev_Ry"]) for x in (L, F) if x.get("E_prev_Ry") and not x.get("done")) * RY * 1000 / 2
        out[m] = {"dE": (F["E_Ry"] - L["E_Ry"]) * RY * 1000 / 2, "dE_err_upper": err, "dexx": [L["dexx"], F["dexx"]],
                  "t": old.get("t")}
    return out


def load_j2(fetch):
    """least-squares J's from the 2-f.u. supercells (per compound, base U). Returns {m: {...}} or {}."""
    man = HERE / "manifest_tc.json"
    if not man.exists():
        return {}
    man = json.loads(man.read_text())
    spec = importlib.util.spec_from_file_location("jexch", HERE / "jexch_submit.py")
    jx = importlib.util.module_from_spec(spec); spec.loader.exec_module(jx)
    mrun = None
    out = {}
    for job in man["jobs"]:
        m = job["material"]
        cands = {c["id"]: c for c in json.loads((ROOT / job["input"]).read_text())}
        rows = []
        for path in job["results"]:
            loc = HERE / "results" / path
            if fetch and not loc.exists():
                if mrun is None:
                    s2 = importlib.util.spec_from_file_location("mrun", ROOT / "infra" / "mrun.py")
                    mrun = importlib.util.module_from_spec(s2); s2.loader.exec_module(mrun)
                try:
                    txt = mrun._read(path)
                    d = json.loads(txt)
                    if d.get("converged_scf"):
                        loc.parent.mkdir(parents=True, exist_ok=True); loc.write_bytes(txt)
                except Exception:
                    pass
            if not loc.exists():
                continue
            d = json.loads(loc.read_text())
            cid = pathlib.Path(path).stem
            if not d.get("converged_scf") or d.get("energy_eV") is None:
                continue
            c = cands[cid]
            rows.append({"id": cid, "cell": c["meta"]["cell"], "config": c["meta"]["config"], "E": d["energy_eV"],
                         "M": d.get("total_mag"), "gap": (d.get("glob") or {}).get("gap"), "pc": jx.pair_coeffs(c)})
        if not rows:
            continue
        cells = sorted({r["cell"] for r in rows})
        y = np.array([r["E"] * 1000.0 for r in rows])
        E0 = np.array([[1.0 if r["cell"] == cc else 0.0 for cc in cells] for r in rows])
        # nested models: all six; J4C = J4N tied; no J4 (J4 is constant inside the type-I cell)
        models = [("full", ["J1", "J2C", "J2N", "J3", "J4C", "J4N"], None),
                  ("J4 tied", ["J1", "J2C", "J2N", "J3", "J4"], "tie"),
                  ("no J4", ["J1", "J2C", "J2N", "J3"], None)]
        for mname, keys, mode in models:
            cols = []
            for k in keys:
                if k == "J4":
                    cols.append([(r["pc"]["J4C"] + r["pc"]["J4N"]) * S2 for r in rows])
                else:
                    cols.append([r["pc"][k] * S2 for r in rows])
            A = np.hstack([np.array(cols).T, E0])
            if np.linalg.matrix_rank(A) == A.shape[1] and len(rows) >= A.shape[1]:
                break
        else:
            out[m] = {"n_configs": len(rows), "status": "underdetermined", "configs": [r["id"] for r in rows]}
            continue
        sol, res, rk, sv = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ sol
        dof = max(len(y) - A.shape[1], 0)
        sig = float(np.sqrt((resid ** 2).sum() / dof)) if dof > 0 else None
        cov = (sig ** 2) * np.linalg.inv(A.T @ A) if sig else None
        J = {k: float(sol[i]) for i, k in enumerate(keys)}
        Jerr = {k: float(np.sqrt(cov[i, i])) for i, k in enumerate(keys)} if cov is not None else {}
        if "J4" in J:
            J["J4C"] = J["J4N"] = J.pop("J4"); Jerr["J4C"] = Jerr["J4N"] = Jerr.pop("J4", None)
        # FM - LCM in the supercells (per ion) for the cross-check against the 15-atom value
        chk = {}
        for cc in cells:
            e = {r["config"]: r["E"] for r in rows if r["cell"] == cc}
            if "FM" in e and "LCM" in e:
                chk[cc] = (e["FM"] - e["LCM"]) * 1000 / 4
        out[m] = {"n_configs": len(rows), "status": "ok", "fit_model": mname, "J_meV": J, "J_err_meV": Jerr, "fit_rms_meV_per_cell": sig,
                  "dE_FM_LCM_supercell_meV_per_ion": chk, "configs": {r["id"]: {"M": r["M"], "gap": r["gap"]} for r in rows},
                  "convention": "H = sum_<ij> J S_i.S_j, J > 0 AF, S^2 = 9/4 folded (broken-symmetry mapping)"}
    return out


# ---------------------------------------------------------------- model temperatures
def shell_vecs(v):
    import itertools
    s = set()
    for p in itertools.permutations(v):
        for sg in itertools.product((1, -1), repeat=3):
            s.add(tuple(a * b for a, b in zip(p, sg)))
    return np.array(sorted(s), float)


_Q = None


def _qgrid(n=48):
    global _Q
    if _Q is None or _Q[0] != n:
        g = (np.arange(n) + 0.5) / n * 2 * np.pi - np.pi
        q = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
        _Q = (n, q)
    return _Q[1]


def rpa_factor(J):
    """k T_RPA / S(S+1)  (meV) for couplings J (meV; J1, J3 inter; J2C/J2N/J4C/J4N intra averaged)."""
    q = _qgrid()
    def jq(v, Jv):
        if Jv == 0:
            return 0.0
        return Jv * np.cos(q @ shell_vecs(v).T).sum(1)
    J2 = 0.5 * (J.get("J2C", 0) + J.get("J2N", 0)); J4 = 0.5 * (J.get("J4C", 0) + J.get("J4N", 0))
    JAB = jq((1, 0, 0), J["J1"]) + jq((1, 1, 1), J.get("J3", 0))
    JAA = jq((1, 1, 0), J2) + jq((2, 0, 0), J4)
    JAB0 = 6 * J["J1"] + 8 * J.get("J3", 0); JAA0 = 12 * J2 + 6 * J4
    A = JAB0 - JAA0 + JAA
    den = A * A - JAB * JAB
    den = np.where(np.abs(den) < 1e-12, 1e-12, den)
    raw = 1.0 / (3.0 * np.mean(A / den))
    # the q -> 0 / (pi,pi,pi) singularities converge slowly on a midpoint grid: normalise by the pure-J1 sum on the same
    # grid and the exact pure-J1 value 2 J1 / W (W = sc Watson integral)
    g = np.cos(q @ shell_vecs((1, 0, 0)).T).sum(1) / 6.0
    raw1 = 1.0 / (3.0 * np.mean(1.0 / (6.0 * (1 - g * g))))
    return raw * (2.0 / WATSON) / raw1


def mf_factor(J):
    """k T_MF / S(S+1) (meV), two-sublattice MF"""
    b = 6 * J["J1"] + 8 * J.get("J3", 0)
    aA = -(12 * J.get("J2N", 0) + 6 * J.get("J4N", 0)); aC = -(12 * J.get("J2C", 0) + 6 * J.get("J4C", 0))
    return ((aA + aC) / 2 + math.sqrt(((aA - aC) / 2) ** 2 + b * b)) / 3.0


def quantum_ratio():
    """T_Q / T_classical(S^2) for S = 3/2: S(S+1)/S^2 * theta_Q/theta_cl, two interpolations"""
    x12, x32 = 0.75, SS1
    t1 = TH_CL - (TH_CL - TH_Q12) * x12 / x32          # in 1/S(S+1)
    t2 = TH_CL - (TH_CL - TH_Q12) * 0.5 / 1.5          # in 1/S
    th = 0.5 * (t1 + t2); dth = 0.5 * abs(t1 - t2)
    return SS1 / S2 * th / TH_CL, SS1 / S2 * dth / TH_CL, th


def mc_json(name):
    p = MC / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def dilution_curve():
    pts = []
    for f in sorted(MC.glob("dil_p*.json")):
        d = json.loads(f.read_text())
        if d.get("Tc"):
            pts.append((d["p"], d["Tc"], d["Tc_err"]))
    pts.sort()
    return pts


def D_of_p(p, curve):
    """T_C(p)/T_C(1) from the MC table (linear interpolation); None outside the table"""
    if not curve:
        return None, None
    ps = np.array([c[0] for c in curve]); ts = np.array([c[1] for c in curve]); es = np.array([c[2] for c in curve])
    t1 = ts[np.argmax(ps)]
    if p > ps.max() + 1e-9 or p < ps.min() - 1e-9:
        return None, None
    v = float(np.interp(p, ps, ts)) / t1
    e = float(np.interp(p, ps, es)) / t1
    return v, e


def classical_factor(Jrel, tag):
    """t = k T_MC / (J1 S^2) for relative couplings Jrel (J1 = 1). Uses an explicit MC run if present
    (mc/model_<tag>.json), else the pure value corrected linearly with the nnn_* sensitivity runs."""
    d = mc_json(f"model_{tag}") if tag else None
    if d and d.get("Tc"):
        return d["Tc"], d["Tc_err"], "explicit MC"
    t0 = TH_CL
    sens = {}
    for k, name in (("J2", "nnn_J2p05"), ("J3", "nnn_J3p05"), ("J4", "nnn_J4p05"), ("J2C", "nnn_J2Cp05")):
        dd = mc_json(name)
        if dd and dd.get("Tc"):
            sens[k] = (dd["Tc"] - t0) / 0.05
    t = t0
    j2 = 0.5 * (Jrel.get("J2C", 0) + Jrel.get("J2N", 0)); j4 = 0.5 * (Jrel.get("J4C", 0) + Jrel.get("J4N", 0))
    used = []
    for k, v in (("J2", j2), ("J3", Jrel.get("J3", 0)), ("J4", j4)):
        if v == 0:
            continue
        if k in sens:
            t += sens[k] * v; used.append(k)
        else:
            # MF-like fallback: t scales with the MF factor
            used.append(k + "(MF)")
            base = mf_factor({"J1": 1.0}); t *= mf_factor({**{"J1": 1.0}, k if k != "J2" else "J2C": v,
                                                         **({"J2N": v} if k == "J2" else {})}) / base
    return t, 0.006 * t, "linear NNN sensitivity " + ",".join(used) if used else "pure J1 (1.4430)"


def perc_pc():
    p = MC / "percolation.json"
    if not p.exists():
        return None
    rows = json.loads(p.read_text())
    xs = [r[0] for r in rows]; ys = [r[1][-1] for r in rows]  # largest L
    for i in range(len(xs) - 1):
        if ys[i] < 0.5 <= ys[i + 1]:
            pc = xs[i] + (xs[i + 1] - xs[i]) * (0.5 - ys[i]) / (ys[i + 1] - ys[i])
            return {"p_c": pc, "err": 0.01, "method": "spanning probability 0.5 crossing, union-find, L = 24/48, 6 samples",
                    "note": "N sites linked only through occupied C sites = fcc site percolation with 1st+2nd neighbours (~0.136)"}
    return None


def model_T(dE_meV, Jrel, tag=None):
    """T (K) at the levels MF / MC / RPA / Q for one DFT dE (meV/ion, = (6 J1 + 8 J3) S^2) and relative couplings"""
    eff = 6 * Jrel.get("J1", 1.0) + 8 * Jrel.get("J3", 0.0)
    J1 = dE_meV / (eff * S2)
    J = {k: v * J1 for k, v in Jrel.items()}
    J["J1"] = J1
    t, te, how = classical_factor(Jrel, tag)
    qr, qre, th = quantum_ratio()
    T_mc = t * J1 * S2 / KB
    return {"J1_meV": J1, "J_meV": J, "T_MF": SS1 * mf_factor(J) / KB, "T_MC": T_mc, "T_MC_err": te * J1 * S2 / KB,
            "T_RPA": SS1 * rpa_factor(J) / KB, "T_Q": T_mc * qr, "mc_source": how}


LEVELS = ("T_MC", "T_Q", "T_RPA", "T_MF")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fetch", action="store_true")
    a = ap.parse_args()
    snap = run_collect(a.fetch)
    dE = dE_table(snap)
    hp = hse_prelim(a.fetch)
    for m, v in hp.items():
        if "HSE" not in dE.get(m, {}):
            dE.setdefault(m, {})["HSE"] = {"dE": v["dE"], "flags": [f"PRELIMINARY HSE (running job; dexx {v['dexx'][0]:.0e}/{v['dexx'][1]:.0e} Ry; "
                                                                    f"last-loop change <= {v['dE_err_upper']:.1f} meV/ion; {v['t']})"]}
    j2 = load_j2(a.fetch)
    curve = dilution_curve()
    qr, qre, thq = quantum_ratio()

    # relative couplings per compound (J1 = 1), from the supercells when available
    Jrel = {}
    for m in ("KVCr", "KVMo", "CrCr"):
        r = j2.get(m)
        if r and r.get("status") == "ok":
            J = r["J_meV"]
            Jrel[m] = {k: v / J["J1"] for k, v in J.items()}
        else:
            Jrel[m] = {"J1": 1.0}

    Jrel["KVCr_2H2O"] = Jrel["KVCr"]  # hydrated model cell: same relative couplings as anhydrous KVCr
    # model temperatures for every (compound, U point / HSE)
    T = {}
    for m, pts in dE.items():
        for tag, v in pts.items():
            T.setdefault(m, {})[tag] = {"dE_meV_per_ion": v["dE"], "flags": v["flags"],
                                        **model_T(v["dE"], Jrel[m], tag=m if m != "KVCr_2H2O" else "KVCr")}
            # J1-only variant (NNN ignored) for comparison
            T[m][tag]["T_MC_J1only"] = TH_CL * v["dE"] / 6.0 / KB

    # composition corrections of the calibrant samples (dilution factor from MC)
    cal = {}
    EXPV = dict(EXP)
    if "KVCr_2H2O" in T:
        EXPV["KVCr_2H2O"] = dict(EXP["KVCr"], note="same sample as KVCr; DFT cell = KV[Cr(CN)6].2H2O (P1, 21 atoms, U 3/3, free FM)")
    for m, e in EXPV.items():
        D, De = D_of_p(e["p"], curve)
        if D is None:
            D, De = 1.0, 0.0
        # derivative for the p error
        D_lo, _ = D_of_p(max(e["p"] - e["p_err"], 0.5), curve)
        D_hi, _ = D_of_p(min(e["p"] + e["p_err"], 1.0), curve)
        Dp_err = 0.5 * abs((D_hi or D) - (D_lo or D))
        rows = {}
        for tag, r in T.get(m, {}).items():
            rows[tag] = {lv: e["T"] / (r[lv] * D) for lv in LEVELS}
            rows[tag]["T_MC_J1only"] = e["T"] / (r["T_MC_J1only"] * D)
        cal[m] = {"exp": e, "D_sample": D, "D_err": math.hypot(De, Dp_err), "f": rows}

    # calibration ratio test: T_model(KVCr sample) / T_model(CrCr sample) vs 376/240 (JUDGE: 1.57 +- 0.25)
    ratio_test = []
    for a_, b_ in (("U3_3", "U3"), ("U2_2", "U2"), ("U4_4", "U4"), ("HSE", "HSE")):
        if a_ in T.get("KVCr", {}) and b_ in T.get("CrCr", {}):
            r_dE = T["KVCr"][a_]["dE_meV_per_ion"] / T["CrCr"][b_]["dE_meV_per_ion"]
            r_T = (T["KVCr"][a_]["T_MC"] * cal["KVCr"]["D_sample"]) / (T["CrCr"][b_]["T_MC"] * cal["CrCr"]["D_sample"])
            ratio_test.append({"pair": f"KVCr {a_} / CrCr {b_}", "dE_ratio": r_dE, "T_model_ratio_samples": r_T,
                               "exp": 376 / 240, "dev_rel": r_T / (376 / 240) - 1,
                               "flags": T["KVCr"][a_]["flags"] + T["CrCr"][b_]["flags"]})

    # predictions for stoichiometric KVMo, per calibrant and U pairing
    preds = []
    def add(label, cm, ctag, ttag, kind):
        if ctag not in T.get(cm, {}) or ttag not in T.get("KVMo", {}):
            return
        e = EXPV[cm]
        rec = {"label": label, "calibrant": cm, "cal_point": ctag, "KVMo_point": ttag, "kind": kind,
               "flags": T["KVMo"][ttag]["flags"]}
        for lv in LEVELS + ("T_MC_J1only",):
            f = e["T"] / (T[cm][ctag][lv] * cal[cm]["D_sample"])
            rec[lv] = f * T["KVMo"][ttag][lv]
        rel = math.hypot(e["T_err"] / e["T"], cal[cm]["D_err"])
        rec["cal_rel_err"] = rel
        preds.append(rec)

    add("base, KVCr-calibrated (3/3 vs 3/2)", "KVCr", "U3_3", "U3_2", "primary")
    add("base, CrCr-calibrated (U 3 vs 3/2)", "CrCr", "U3", "U3_2", "primary")
    for uv in (2, 4):
        add(f"U_V {uv}: KVMo({uv},2) vs KVCr({uv},2)", "KVCr", f"U{uv}_2", f"U{uv}_2", "U-spread")
        add(f"U_V {uv}: KVMo({uv},2) vs KVCr({uv},4)", "KVCr", f"U{uv}_4", f"U{uv}_2", "U-spread")
        add(f"U_V {uv}: KVMo({uv},0) vs KVCr({uv},2) [U_Mo = 0, FM transfer-unstable]", "KVCr", f"U{uv}_2", f"U{uv}_0", "flagged")
        add(f"U {uv} diagonal: KVMo({uv},2) vs CrCr({uv})", "CrCr", f"U{uv}", f"U{uv}_2", "U-spread")
    add("base, KVCr.2H2O-calibrated (hydrated cell, U 3/3)", "KVCr_2H2O", "U3_3", "U3_2", "hydrated")
    add("HSE06, KVCr-calibrated", "KVCr", "HSE", "HSE", "HSE")
    add("HSE06, CrCr-calibrated", "CrCr", "HSE", "HSE", "HSE")

    def stats(sel, lv="T_MC"):
        v = np.array([p[lv] for p in preds if p["kind"] in sel])
        return (float(v.mean()), float(v.min()), float(v.max()), len(v)) if len(v) else None

    prim = [p for p in preds if p["kind"] == "primary"]
    central = float(np.mean([p["T_MC"] for p in prim])) if prim else None
    spreadU = stats(("primary", "U-spread"))
    # error budget (relative, 1 sigma-ish): U spread (half range / central), calibrant disagreement, exp T + sample p,
    # NNN/model (spread across MC/Q/RPA/MF after calibration), HSE (if present)
    budget = {}
    if central:
        budget["U (half-range of primary + U-spread pairings)"] = 0.5 * (spreadU[2] - spreadU[1]) / central
        if len(prim) == 2:
            budget["calibrant (KVCr vs CrCr)"] = 0.5 * abs(prim[0]["T_MC"] - prim[1]["T_MC"]) / central
        budget["exp. T_C + sample composition (calibrant)"] = float(np.mean([p["cal_rel_err"] for p in prim]))
        lv = [np.mean([p[l] for p in prim]) for l in LEVELS]
        budget["Heisenberg level (MC/Q/RPA/MF after calibration)"] = 0.5 * (max(lv) - min(lv)) / central
        nn = np.mean([abs(p["T_MC"] - p["T_MC_J1only"]) for p in prim])
        has_j2 = any(Jrel[m] != {"J1": 1.0} for m in Jrel)
        has_j4 = any("J4C" in Jrel[m] for m in Jrel)
        # uncertainty of the beyond-NN correction: 30 % of its size with J4 fitted, 50 % without (J4 not resolved),
        # 5 % absolute if no supercell data at all
        budget["beyond-NN couplings (uncertainty of the J2/J3/J4 correction)"] = \
            float((0.3 if has_j4 else 0.5) * nn / central) if has_j2 else 0.05
        h = [p for p in preds if p["kind"] == "HSE"]
        if h:
            budget["functional (HSE06 vs PBE+U base)"] = abs(np.mean([p["T_MC"] for p in h]) - central) / central
    tot = math.sqrt(sum(v * v for v in budget.values())) if budget else None

    # consistency checks against vacancy-rich V-Cr samples, and the V1.37[Mo(CN)6] parent
    checks = {}
    fV = cal["KVCr"]["f"].get("U3_3", {}).get("T_MC")
    base_VCr = T["KVCr"]["U3_3"]["T_MC"] * fV if fV else None  # = 376 K / D_sample(KVCr): stoichiometric V-Cr in the model
    D_KVMo73 = mc_json("model_KVMo_p0.730")  # dilution with the KVMo NNN set (if run)
    for k, c in CHECKS.items():
        D, De = D_of_p(c["p"], curve)
        rec = {"exp": c, "D_MC": D}
        if k.startswith("VCr") and base_VCr:
            rec["T_stoich_VCr_model"] = base_VCr
            if c.get("w"):
                for nm, key in (("check_Ferlay_lit", "T_pred_rJlit"), ("check_Ferlay_rJ1", "T_pred_rJ1")):
                    mm = mc_json(nm)
                    if mm and mm.get("Tc"):
                        rec[key] = base_VCr * mm["Tc"] / TH_CL
                        rec[key + "_ratio_exp_over_pred"] = c["T"] / rec[key]
            elif D is not None:
                rec["T_pred"] = base_VCr * D
                rec["ratio_exp_over_pred"] = c["T"] / rec["T_pred"]
        if k == "VMo_Magott" and D is not None:
            if D_KVMo73 and D_KVMo73.get("Tc") and mc_json("model_KVMo"):
                rec["D_MC_KVMo_NNN"] = D_KVMo73["Tc"] / mc_json("model_KVMo")["Tc"]  # check only (J1-only curve used)
            rec["T_pred_from_DFTcal_model"] = (central or 0) * D
            rec["ratio_exp_over_pred"] = c["T"] / rec["T_pred_from_DFTcal_model"] if central else None
            emp = 1.0 - (1.0 - 330.0 / 376.0) * (1.0 - c["p"]) / (1.0 / 3.0)  # Holmes-Girolami crystalline V-Cr: T(2/3)/T(1)
            rec["D_empirical_VCr"] = emp
            rec["T_stoich_parent_cal"] = {"MC_dilution": {"417": c["T"] / D, "340": c["T_lo"] / D, "454": c["T_hi"] / D},
                                          "empirical_VCr_dilution": {"417": c["T"] / emp, "340": c["T_lo"] / emp, "454": c["T_hi"] / emp}}
            # implied exchange ratio V-Mo / V-Cr after composition correction (vs the PBE+U dE ratio)
            if base_VCr:
                rec["implied_T_ratio_VMo_over_VCr"] = {"MC_dilution": c["T"] / D / base_VCr, "empirical": c["T"] / emp / base_VCr}
            rec["_D"] = D; rec["_emp"] = emp
        checks[k] = rec

    # B: the same Heisenberg model calibrated directly on the V-Mo parent (417 K at p = 0.73)
    anc = None
    mg = checks.get("VMo_Magott", {})
    if mg.get("_D"):
        c = CHECKS["VMo_Magott"]
        Ds = [mg["_D"], mg["_emp"]]
        vals = [c["T"] / d for d in Ds]
        Bc = float(np.mean(vals))
        sD = 0.5 * abs(vals[0] - vals[1]) / Bc
        sT = (c["T_hi"] - c["T"]) / c["T"]           # Bloch-exponent ambiguity (417 vs 454 K) as 1 sigma
        Be = Bc * math.hypot(sD, sT)
        anc = {"T_K": Bc, "T_err_K": Be, "variants": {"MC dilution": vals[0], "empirical V-Cr dilution": vals[1]},
               "hard_floor_K": c["T_lo"] / max(Ds), "note": "T(1:1) = T_C(parent)/D(0.73); parent T_C from a Bloch extrapolation of "
                                                            "M(T) measured to 340 K only"}
    # C: combination of A (DFT-ratio, V-Cr/Cr-Cr calibrated) and B (V-Mo parent calibrated), Birge-inflated
    comb = None
    if central and anc:
        A_, sA = central, central * tot
        B_, sB = anc["T_K"], anc["T_err_K"]
        w1, w2 = 1 / sA ** 2, 1 / sB ** 2
        mu = (w1 * A_ + w2 * B_) / (w1 + w2)
        se = 1 / math.sqrt(w1 + w2)
        chi2 = w1 * (A_ - mu) ** 2 + w2 * (B_ - mu) ** 2
        comb = {"T_K": mu, "T_err_K": se * max(1.0, math.sqrt(chi2)), "chi2_1dof": chi2, "tension_A_over_B": A_ / B_,
                "range_K": [min(A_ - sA, B_ - sB), max(A_ + sA, B_ + sB)]}
    for k in checks:
        checks[k].pop("_D", None); checks[k].pop("_emp", None)

    # vacancy table for KVMo. A(p) = A * D_MC(p); B(p) = parent T_C rescaled from p = 0.73 with each dilution law (mean);
    # C(p) = Birge-weighted combination of A(p) and B(p)
    def D_emp(p):
        return 1.0 - (1.0 - 330.0 / 376.0) * (1.0 - p) / (1.0 / 3.0)
    vac = []
    p0 = CHECKS["VMo_Magott"]["p"]
    D0, _ = D_of_p(p0, curve)
    for p in (1.0, 0.97, 0.95, 0.9, 0.85, 0.8, p0, 2 / 3, 0.6, 0.5):
        D, De = D_of_p(p, curve)
        if D is None:
            continue
        row = {"p": p, "vacancy_fraction": 1 - p, "D_MC": D, "D_emp_VCr": D_emp(p), "M_uncomp_muB_per_fu": 3 * (1 - p)}
        if central:
            row["T_A_DFTcal"] = central * D
            row["T_A_err"] = central * D * (tot or 0)
        if anc and D0:
            v = [CHECKS["VMo_Magott"]["T"] * D / D0, CHECKS["VMo_Magott"]["T"] * D_emp(p) / D_emp(p0)]
            Bp = float(np.mean(v))
            row["T_B_parentcal"] = Bp
            sT = (CHECKS["VMo_Magott"]["T_hi"] - CHECKS["VMo_Magott"]["T"]) / CHECKS["VMo_Magott"]["T"]
            row["T_B_err"] = Bp * math.hypot(0.5 * abs(v[0] - v[1]) / Bp, sT)
        if central and anc and D0:
            a1, s1 = row["T_A_DFTcal"], row["T_A_err"]; b1, s2 = row["T_B_parentcal"], row["T_B_err"]
            w1, w2 = 1 / s1 ** 2, 1 / s2 ** 2
            mu = (w1 * a1 + w2 * b1) / (w1 + w2)
            chi2 = w1 * (a1 - mu) ** 2 + w2 * (b1 - mu) ** 2
            row["T_C_combined"] = mu
            row["T_C_err"] = max(1.0, math.sqrt(chi2)) / math.sqrt(w1 + w2)
        vac.append(row)

    summary = {
        "what": "Heisenberg T_C model for d3/d3 PBAs (Track L); classical MC (S^2), quantum-corrected, RPA (S = 3/2), MF",
        "date": __import__("time").strftime("%Y-%m-%d %H:%M"),
        "conventions": {"H": "sum_<ij> J_ij S_i.S_j, J > 0 AF, S = 3/2", "mapping": "E(FM) - E(LCM) = (6 J1 + 8 J3) S^2 per ion",
                        "T_MC": "classical spins |S| = 3/2 (S^2 scale), 1.4430 J1 S^2 for J1 only",
                        "T_Q_over_T_MC": [qr, qre], "theta_Q(3/2)": thq,
                        "T_RPA_pure_J1": "1.3189 J1 S(S+1)"},
        "dE_meV_per_ion": dE, "J_supercell": j2, "J_relative_used": Jrel, "model_T_K": T,
        "dilution_curve_MC": [{"p": p, "Tc_over_J1S2": t, "err": e} for p, t, e in curve],
        "calibration": cal, "calibration_ratio_test": ratio_test, "KVMo_predictions": preds,
        "KVMo_stoichiometric": {"A_DFT_ratio_VCr_CrCr_calibrated": {"T_K": central, "T_err_K": (central or 0) * (tot or 0),
                                                                     "U_pairing_stats_mean_min_max_n": spreadU,
                                                                     "error_budget_rel": budget, "total_rel_err": tot},
                                "B_parent_calibrated": anc, "C_combined": comb},
        "checks": checks, "KVMo_vs_vacancy": vac,
        "percolation_threshold_C_sublattice": perc_pc(),
        "finite_T_compensation": {m: (mc_json(f"mt_{m}") or {}).get("rows") for m in ("KVCr", "KVMo")},
        "nnn_sensitivity_MC": {k: (mc_json(k) or {}).get("Tc") for k in ("nnn_J2p05", "nnn_J2m05", "nnn_J3p05", "nnn_J4p05", "nnn_J2Cp05")},
        "explicit_MC_models": {k: {"Tc_over_J1S2": (mc_json(k) or {}).get("Tc"), "err": (mc_json(k) or {}).get("Tc_err"),
                                   "J": (mc_json(k) or {}).get("J")} for k in ("model_KVCr", "model_KVMo", "model_CrCr", "model_KVMo_p0.730")},
    }
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1, default=float))

    # ---- print
    print("== dE(FM-LCM) meV/ion, J1_eff, model T (K) ==")
    for m in T:
        for tag, r in T[m].items():
            print(f"  {m:5s} {tag:6s} dE {r['dE_meV_per_ion']:6.1f}  J1 {r['J1_meV']:6.2f} meV  MC {r['T_MC']:6.0f}  Q {r['T_Q']:6.0f}  "
                  f"RPA {r['T_RPA']:6.0f}  MF {r['T_MF']:6.0f}  {'; '.join(r['flags'])}")
    print("== calibration f = T_exp / T_model(sample) ==")
    for m, c in cal.items():
        for tag, f in c["f"].items():
            print(f"  {m:5s} {tag:6s} D_sample {c['D_sample']:.3f}  f_MC {f['T_MC']:.3f}  f_Q {f['T_Q']:.3f}  f_RPA {f['T_RPA']:.3f}  f_MF {f['T_MF']:.3f}")
    print("== ratio test T_model(KVCr)/T_model(CrCr) vs 376/240 = 1.567 ==")
    for r in ratio_test:
        print(f"  {r['pair']:24s} dE ratio {r['dE_ratio']:.3f}  model-T ratio {r['T_model_ratio_samples']:.3f}  dev {100 * r['dev_rel']:+.1f} %  {'; '.join(r['flags'])}")
    print("== stoichiometric KV[Mo(CN)6] ==")
    for p in preds:
        print(f"  {p['label']:58s} MC {p['T_MC']:5.0f}  Q {p['T_Q']:5.0f}  RPA {p['T_RPA']:5.0f}  MF {p['T_MF']:5.0f}  "
              f"(J1-only {p['T_MC_J1only']:5.0f}) {'; '.join(p['flags'])}")
    if central:
        print(f"  central {central:.0f} K, total rel err {tot:.3f} -> +- {central * tot:.0f} K; budget {json.dumps({k: round(v, 3) for k, v in budget.items()})}")
    print("== checks ==")
    for k, c in checks.items():
        print(" ", k, json.dumps({kk: (round(v, 3) if isinstance(v, float) else v) for kk, v in c.items() if kk != "exp"}, default=str))
    print("== KV[Mo(CN)6] stoichiometric: A (DFT ratio) / B (parent-calibrated) / C (combined) ==")
    if central:
        print(f"  A {central:.0f} +- {central * tot:.0f} K")
    if anc:
        print(f"  B {anc['T_K']:.0f} +- {anc['T_err_K']:.0f} K  variants {json.dumps({k: round(v) for k, v in anc['variants'].items()})}  floor {anc['hard_floor_K']:.0f} K")
    if comb:
        print(f"  C {comb['T_K']:.0f} +- {comb['T_err_K']:.0f} K  (chi2 {comb['chi2_1dof']:.2f}, A/B {comb['tension_A_over_B']:.2f})")
    print("== KVMo vs vacancy fraction ==")
    for r in vac:
        print(f"  p {r['p']:.3f}  D {r['D_MC']:.3f}  A {r.get('T_A_DFTcal', float('nan')):5.0f} +- {r.get('T_A_err', 0):3.0f}  "
              f"B {r.get('T_B_parentcal', float('nan')):5.0f} +- {r.get('T_B_err', 0):3.0f}  C {r.get('T_C_combined', float('nan')):5.0f} +- {r.get('T_C_err', 0):3.0f}  "
              f"M_uncomp {r['M_uncomp_muB_per_fu']:.2f}")


if __name__ == "__main__":
    main()
