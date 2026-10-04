"""Assemble ord2_summary.json from ce2_*.json, relaxlevel_*.json, torder_*.json, sro_*.json and the live job states.
Rerun recipe after more jobs finish:
  python fetch.py YBMFO YBCFO
  python analyze2.py YBMFO YBMFO_P4n_h110.json ../../../infra/jobs_src/lcm_ord2/tasks_ybmfo_v4.json
  python relaxlevel.py YBMFO YBMFO_P4n_h110.json ../../../infra/jobs_src/lcm_ord2/tasks_ybmfo_v4.json [still-relaxing ids]
  python analyze2.py YBMFO ...            (again: picks up relaxlevel_YBMFO.json for the 'best' dataset)
  (same three steps for YBCFO with YBCFO_lay_h110.json / tasks_ybcfo_v4.json)
  bash run_torder.sh ; python sro.py ce2_YBMFO.json E_PM_best_bw_bayes best 16 ; python sro.py ce2_YBCFO.json E_PM_best_bw_bayes best 16
  python make_summary.py"""
import json, os, glob, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REG = os.path.join(HERE, "../../../infra/job_registry.jsonl")
MAIN = "E_PM_best_bw_bayes"
VARIANTS = [("E_PM_best_bw_bayes", "best"), ("E_PM_best_bw", "best"), ("E_PM_best", "ip1+ap+Y"), ("E_PM_low", "best"),
            ("E_PM_all_bayes", "best"), ("E_PM_clean_bayes", "best")]


def jl(p):
    return json.load(open(os.path.join(HERE, p))) if os.path.exists(os.path.join(HERE, p)) else None


def job_states():
    try:
        import modal
        v = modal.Volume.from_name("magdisc-data")
    except Exception:
        return {}
    jobs = []
    for l in open(REG):
        try:
            j = json.loads(l)["job"]
        except Exception:
            continue
        if j.startswith("lcm/ord2/v4/") and j not in jobs:
            jobs.append(j)
    out = {}
    for j in jobs:
        try:
            out[j.split("/")[-1]] = json.loads(b"".join(v.read_file(f"jobs/{j}/status.json"))).get("state")
        except Exception:
            out[j.split("/")[-1]] = "unknown"
    return out


summ = {"generated": time.strftime("%Y-%m-%d %H:%M"), "status": "PROVISIONAL", "compounds": {}}
for tag in ("YBMFO", "YBCFO"):
    ce = jl(f"ce2_{tag}.json")
    if ce is None:
        continue
    arr = ce["arrangements"]
    tab = {}
    for a, e in sorted(arr.items(), key=lambda kv: kv[1]["best"]["E"] if kv[1].get("best") else 1e9):
        b = e.get("best")
        tab[a] = {"E_best_meV_fu": None if not b else round(b["E"], 1), "sigma": None if not b else round(b["sigma"], 1),
                  "source": None if not b else b["src"],
                  "E_low_protocol": round(e["low"]["E"], 1) if e.get("low") else None,
                  "E_all_protocol": round(e["all"]["E"], 1) if e.get("all") else None,
                  "E_clean_protocol": round(e["clean"]["E"], 1) if e.get("clean") else None,
                  "E_relaxlevel": round(e["relaxlevel"]["E"], 1) if e.get("relaxlevel") else None,
                  "n_runs_nonFM": e.get("n"), "ct_runs": [c for c, v in e.get("configs", {}).items() if v["ct_sites"]],
                  "Pi": {k: round(v, 3) for k, v in e["Pi"].items() if k != "ip4"}}
    m = ce["ce"][MAIN + "_best"]
    ent = {"reference": f"{tag}_RS_A (rock-salt, 36-atom 1x1x2 cell) = 0",
           "exchange_fit_clean_noFM": {"J_meV_unit_spins": ce["J_clean_noFM_meV"], "rms_meV_cell": round(ce["J_fit_rms_meV_cell"], 1)},
           "n_runs": ce["n_runs"], "n_charge_transfer_runs": ce["n_ct_runs"],
           "FM_above_PM_line_meV_fu": {"mean": round(float(np.mean(ce["fm_heisenberg_dev_meV_fu"])), 1)} if ce["fm_heisenberg_dev_meV_fu"] else None,
           "PM_ordering_energies_meV_per_fu": tab,
           "CE_main": {"dataset": MAIN, "classes": m["classes"], "prior_scale": m.get("prior_scale"),
                       "ECI_meV_per_bond": {k: round(v, 2) for k, v in m["V_meV_per_bond"].items()},
                       "ECI_boot_std": {k: round(v, 2) for k, v in m["V_boot_std"].items()},
                       "wrms_meV_fu": round(m["wrms"], 1), "LOO_CV_weighted_meV_fu": round(m["cv"], 1), "LOO_CV_unweighted": round(m["cv_unw"], 1),
                       "LOO_resid_meV_fu": m["loo_resid"],
                       "convention": "E/f.u. = c0 + sum_k V_k (bonds_k per f.u.) Pi_k, sigma = +1 Mn/Cu, -1 Fe; V > 0 favours hetero pairs"},
           "CE_variants": {}, "T_order": {}, "SRO": None}
    for key, csel in VARIANTS:
        lst = ce["ce"].get(key)
        if not lst:
            continue
        e = ce["ce"][key + "_best"] if csel == "best" else [x for x in lst if "+".join(x["classes"]) == csel][0]
        tf = jl(f"torder_{tag}_{key}{'_' + csel if csel != 'best' else ''}.json")
        ent["CE_variants"][f"{key}:{csel}"] = {"ECI": {k: round(v, 1) for k, v in e["V_meV_per_bond"].items()}, "CVw": round(e["cv"], 1),
                                                "wrms": round(e["wrms"], 1)}
        if tf:
            ent["T_order"][f"{key}:{csel}"] = {"q_ground": tf["q_ground"], "Tpeak_C_by_L": tf["Tpeak"], "T_order_K_L16": round(tf["T_order_K"]),
                                               "boot_mean": round(tf["boot_T16_mean"]), "boot_std": round(tf["boot_T16_std"]),
                                               "boot_p16_p84": [round(x) for x in tf["boot_T16_p16_p84"]]}
    s = jl(f"sro_{tag}_{MAIN}.json")
    if s:
        ent["SRO"] = {"CE": MAIN, "L": s["L"], "q_ground": s["q_ground"],
                      "by_T": [{"T_K": r["T"], "eta_q0": round(r["eta_q0"], 3), "x_antisite_RS": round(r["x_AS_RS"], 3) if r["x_AS_RS"] is not None else None,
                                "hetero_frac_ip1_ap_Y": [round(r["hetero_frac"][k], 3) for k in ("ip1", "ap", "Y")]} for r in s["T"]]}
    summ["compounds"][tag] = ent
summ["job_states"] = job_states()
# ---- interpretation (2026-10-02 analysis pass; update the numbers/text when the pending jobs land) ----
summ["benchmark_YBCFO_experiment"] = {
    "source": "Morin et al., Nat. Commun. 7, 13758 (2016)",
    "synthesis": "solid state, 1150 C (1423 K) 50 h in O2; cooled at 5/100/300/500 K/h or quenched in liquid N2",
    "observable": "P4mm split-site Cu/Fe occupation difference = polar (LAY-type, q=(0,0,1)) bipyramid-orientation order parameter",
    "eta_LAY_exp": {"quenched": 0.10, "5 K/h": 0.16},
    "LOG_expectation": "LOG.md: BPR_F (all bipyramids Cu-O-Fe, random orientation) is the experimentally proposed state; TRACKL_LOG: YBaCuFeO5 disordered in-plane experimentally",
}
summ["conclusions"] = {
    "YBMFO_T_order_K": {"main_CE": 915, "main_CE_L8_12_16": [915, 965, 915], "cooling_run_L16_jump": "973-1073 K",
                        "bootstrap_p16_p84": [800, 1136], "CE_variant_range_RS_ground_state": [836, 1210],
                        "headline": "T_order = 950 (+250/-150) K  (~680 C; 530-930 C)"},
    "YBMFO_SRO_at_900_1300C": "no LRO (eta_RS <= 0.04, finite-size floor); hetero NN-pair fractions ip1/ap/Y = 0.56-0.59 / 0.60-0.66 / 0.58-0.61 (random 0.50, rock-salt 1.00), i.e. Warren-Cowley alpha_ip1 = -0.12..-0.18",
    "YBMFO_equilibrium_antisites_below_T_order": {"873 K": 0.13, "773 K": 0.064, "673 K": 0.03},
    "YBCFO_T_order_K": {"main_CE": 1118, "bootstrap_p16_p84": [900, 1200], "alt_CE": 1358, "ordered_state": "LAY (polar P4mm bipyramid orientation)",
                        "eta_LAY_equilibrium_main_CE": {"1423 K": "~0 (floor)", "1273 K": 0.06, "1173 K": 0.19, "1073 K": 0.71}},
    "benchmark_verdict": "Protocol T_order(YBCFO) = 1120 K (0.79 x the 1423 K sintering T; 900-1360 K) is consistent with the weak measured LRO (0.10-0.16) only if B-site exchange freezes at >~1150 K during cooling; no gross overestimate, but a 10-25 % overestimate cannot be excluded. Practical implication: ordering that needs equilibration below ~1150 K is not reached by cooling/annealing in these 112 Cu/Fe (and by analogy Mn/Fe) oxides.",
    "YBMFO_verdict": "NOT achievable by conventional annealing: T_order (950 K, 800-1210 K) lies below the 1173-1573 K synthesis/anneal window and below the YBaCuFeO5 benchmark's own T_order (1120 K), where slow cooling yields only eta ~ 0.16. Rock-salt order >~90 % (x_AS < 5 %) needs equilibration at <~800 K (marginal at best, only at the top of the error bar). Consistent with the B-disordered GdBaMnFeO5 / NdBaMnFeO5 literature.",
    "why_lower_than_earlier_estimates": "Earlier lane estimates (Ising 2000 K from G/FM ordering energies; PM-proxy 1300 K) assumed ionic Mn2+/Fe3+ everywhere. In the PM protocol many non-rock-salt arrangements lower their energy by Mn2+ + Fe3+ -> Mn3+ + Fe2+ charge transfer (restoring a 2+/3+ checkerboard); random SQS_D is only 116 +- 15 meV/f.u. above rock salt (vs 230 predicted), an NN antisite pair costs ~0.2 eV in the PM state (0.84 eV in G order), and in-plane stripes (STRp_C, +25) and Y-homo stacking (RSapY, +18 meV/f.u.) lie close to rock salt.",
    "provisional": ["SQS_D enters at relax level (r0, 70 Ry/half k, calibrated +-4 meV/f.u. on the D cell via AS1_D); protocol SCFs running (main + helpers h1/h2)",
                    "SQS_F still relaxing; AS1_D has 1/5 non-FM configs (helpers running); YBCFO BPR_F 1/5 configs, YBCFO SQS_D/CHKmix_A still relaxing (relax-level, unconverged), RSYap_A relax-level only",
                    "LAY_A / LAYY_A protocol SCFs are metastable CT states 115-390 meV/f.u. above the relax-level state; relax-level used",
                    "pair CE quality is poor (weighted LOO-CV 39 meV/f.u. YBMFO, 33 YBCFO); T_order spread across CE variants (836-1210 K) is the dominant uncertainty",
                    "E_PM depends on the exchange model by 10-40 meV/f.u. for near-RS arrangements (species-resolved J preferred: rms 27 vs 75 meV/cell); STRp_C is +25 (species J) but -13 meV/f.u. with species-independent J",
                    "the ionic-only (clean) CE predicts a non-RS ground state q=(0,1,1/2) (in-plane stripes + Y-homo stacking), an extrapolation outside the training set: needs a DFT check (STRp with Y-homo stacking, F-type 72-atom cell)",
                    "no vibrational/magnetic-entropy contributions; PM = Heisenberg intercept (FM excluded: +81 meV/f.u. above the PM line on average)"]}
json.dump(summ, open(os.path.join(HERE, "ord2_summary.json"), "w"), indent=1, default=float)
print("wrote ord2_summary.json;", {k: len(v["PM_ordering_energies_meV_per_fu"]) for k, v in summ["compounds"].items()})
