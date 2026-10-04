"""Write tracks/lcm/pba/vac/summary.json from cached results (run collect_vac.py --fetch and
tracks/lcm/pba/defects/collect.py --fetch first). No volume access here.

Contents: vacancy-cell numbers (M_cell, gap, edge spins, windows, aligned windows, opposite-spin edges, defect states),
the K-vacancy polaron variants (harvest), the DFT lability energies with Boltzmann fractions, and the realistic
defect-concentration scenarios with their M(0) and window consequences (README sections 3-4).
"""
import json
import math
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
DEFC = HERE.parent / "defects" / "results" / "jobs" / "lcm" / "pba" / "def"
KB = 8.617333e-5


def jl(p):
    try:
        return json.loads(pathlib.Path(p).read_text())
    except Exception:
        return None


def r3(x):
    return None if x is None else round(float(x), 3)


def vac_block(o):
    if not o:
        return None
    keep = {k: r3(o.get(k)) for k in ("M", "M_expected", "absM", "gap", "win_VB_cell", "win_CB_cell", "align_shift_eV",
                                       "align_shift_K3p_eV", "opp_top_rel_hostVBM", "opp_bot_rel_hostCBM",
                                       "same_top_rel_hostVBM", "same_bot_rel_hostCBM", "win_VB_aligned", "win_CB_aligned",
                                       "win_ratio_VB", "win_ratio_CB", "E_eV")}
    keep.update(vbm=o.get("vbm"), cbm=o.get("cbm"), kill_opp_near_edge=o.get("kill_opp_near_edge"),
                kill_window_cut=o.get("kill_window_cut"), host=o.get("host"),
                capped_V=[{"i": x["i"], "m_sphere": r3(x["m_sphere"]), "m_d": r3(x["m_d"])} for x in o.get("capped_V", [])],
                bulk_V=[{"i": x["i"], "m_sphere": r3(x["m_sphere"]), "m_d": r3(x["m_d"])} for x in o.get("bulk_V", [])],
                M_sites=[{"i": x["i"], "el": x["el"], "m_sphere": r3(x["m_sphere"]), "m_d": r3(x["m_d"])}
                         for x in o.get("metal_moments", []) if x["el"] != "V"],
                near_gap_states=o.get("near_gap_states"), relax=o.get("relax"))
    return keep


def main():
    vs = jl(HERE / "results" / "vac_summary.json") or {}
    man = jl(HERE / "manifest.json") or {}
    out = {"lane": "Track L (session 66190d) PBA vacancy + building-block lability, tracks/lcm/pba/vac/",
           "protocol": "PBE+U (V 3, Cr 3, Mo 2 eV ortho-atomic), PseudoDojo SR NC 110 Ry, 4 f.u. cells at the stage-1 "
                       "lattice; scf 2x2x2 shifted, nscf 4x4x4 Gamma-centred; alignment on the far-field CN 3sigma level "
                       "of the pristine 4 f.u. cell (defects/collect.py analyse(), unchanged)",
           "jobs": [{k: j.get(k) for k in ("job", "call_id", "t_submit", "fn", "cmd", "positions_from", "positions_n_scf")}
                    for j in man.get("jobs", [])]}
    vac = {}
    for m in ("KVCr", "KVMo"):
        best = f"{m}_VAC" if vs.get(f"{m}_VAC") else f"{m}_VAC_r"
        vac[m] = {"source": best if vs.get(best) else None,
                  "note": ("final relaxed cell" if best.endswith("VAC") else
                           "SKIP_RELAX rerun at the last positions of the still-running relax (residual max force "
                           "~0.004 Ry/bohr on water H and K; energy still descending a few meV/step)"),
                  "VAC": vac_block(vs.get(f"{m}_VAC")), "VAC_r": vac_block(vs.get(f"{m}_VAC_r")),
                  "VAC_w1_zeolitic_water": vac_block(vs.get(f"{m}_VAC_w1"))}
    # void ("cavity") state: lowest empty band with water weight > 0.5 in each spin, aligned to the host VBM (CN 3sigma)
    for m in ("KVCr", "KVMo"):
        for key in ("VAC_r", "VAC", "VAC_w1_zeolitic_water"):
            b = vac[m].get(key)
            cid = {"VAC_r": f"{m}_VAC_r", "VAC": f"{m}_VAC", "VAC_w1_zeolitic_water": f"{m}_VAC_w1"}[key]
            r = jl(HERE / "results" / "jobs" / "lcm" / "pba" / "def" / cid / "results" / f"{cid}.json")
            if not b or not r or b.get("align_shift_eV") is None:
                continue
            sh, vbm, cbm = b["align_shift_eV"], b["host"]["VBM"], b["host"]["CBM"]
            k3 = b.get("align_shift_K3p_eV")
            sv = 1  # V moments negative -> V-majority = spin index 1 (checked by collect_vac)
            vs_ = {}
            for sk, rows in r["band_table"].items():
                w = [x for x in rows if x["occ"] < 0.02 and (x["groups"].get("H_w", 0) + x["groups"].get("O_w", 0)
                                                            + x["groups"].get("H_z", 0) + x["groups"].get("O_z", 0)) > 0.5]
                if w:
                    x = min(w, key=lambda y: y["Emin"])
                    vs_["V-maj" if int(sk[1:]) == sv else "M-maj"] = {
                        "E_rel_hostVBM": [r3(x["Emin"] - sh - vbm), r3(x["Emax"] - sh - vbm)],
                        "E_rel_hostCBM": [r3(x["Emin"] - sh - cbm), r3(x["Emax"] - sh - cbm)],
                        "E_rel_hostCBM_K3p_alignment": None if k3 is None else r3(x["Emin"] - k3 - cbm),
                        "groups": x["groups"], "region": x["region"]}
            b["void_state"] = vs_
    out["vacancy_cells"] = vac
    # K-vacancy polaron variants (harvest from the defects cache)
    K = jl(DEFC / "KVMo_KVAC" / "results" / "KVMo_KVAC.json")
    pm = jl(DEFC / "KVMo_KVAC_polM" / "results" / "KVMo_KVAC_polM.json")
    pv = jl(DEFC / "KVMo_KVAC_polV" / "results" / "KVMo_KVAC_polV.json")
    peek = jl(HERE.parent / "defects" / "results" / "relax_peek.json") or {}
    pol = {}
    if K:
        EK = K["scf"]["energy_eV"]
        if pm and pm.get("scf"):
            pol["E_polM_minus_KVAC_eV"] = r3(pm["scf"]["energy_eV"] - EK)
            if (pm.get("relax") or {}).get("energy_eV") is not None:
                pol["E_polM_relax_end_minus_KVAC_eV"] = r3(pm["relax"]["energy_eV"] - EK)
            pol["polM_state"] = ("done" if pm.get("done") else "nscf/projwfc pending") + \
                "; relax stopped at RELAX_MAXSEC unconverged (total force 0.051 Ry/bohr), so this is an upper bound"
        if pv and pv.get("scf"):
            pol["E_polV_minus_KVAC_eV"] = r3(pv["scf"]["energy_eV"] - EK)
            pol["polV_state"] = "done" if pv.get("done") else "scf done"
        elif peek.get("KVMo_KVAC_polV", {}).get("E_last_eV") is not None:
            pol["E_polV_minus_KVAC_eV_relax_in_progress"] = r3(peek["KVMo_KVAC_polV"]["E_last_eV"] - EK)
            pol["polV_state"] = "relax still running (last BFGS energy)"
    out["K_vacancy_polarons_KVMo"] = pol
    # lability energetics (Gate-1 PBE+U, 4 f.u. cells) -> Boltzmann fractions
    E = {"FLIP_KVCr": 0.177, "FLIP_KVMo": 0.182, "AS_pair_KVCr": 0.927, "AS_pair_KVMo": 0.565,
         "isomer_KCrV_per_fu_base": 0.626, "isomer_KCrV_per_fu_range": [0.442, 0.810]}
    bz = {k: {str(T): float("%.2g" % math.exp(-v / (KB * T))) for T in (300, 340, 373)}
          for k, v in E.items() if isinstance(v, float)}
    out["lability_DFT"] = {"energies_eV": E, "boltzmann_exp(-E/kT)": bz,
                           "reading": "equilibrium populations are negligible for antisites and ~1e-3 per CN for "
                                      "flips; real antisite/flip contents are kinetic (set by building-block lability "
                                      "during assembly), not thermodynamic"}
    out["concentration_scenarios"] = {
        "rule": "M(0) = |3v + u - w| muB per f.u. (lane sign: M = w - 3v - u; v = [M(CN)6] vacancy fraction, "
                "w = V(III) fraction, u = Mo(IV) fraction); antisites and CN flips keep M = 0 (isospin)",
        "KVCr": {"v": [0.0, 0.08], "v_best_estimate": 0.04, "w": [0.0, 0.05], "antisite_pairs_per_fu": "<1e-4",
                 "CN_flip_per_CN": "<=1e-3", "M0_muB_per_fu": [0.0, 0.25], "M0_measured": 0.125,
                 "basis": "Holmes & Girolami 1999 z = 1.00 nominal, M_sat 0.125; best alkali V-Cr samples z = 0.92-0.95; "
                          "inert [Cr(CN)6]3- (k_aq(Cr3+) 2.4e-6 s-1); V(II) prefers the N end (Nelson & Miller 2008), "
                          "Lewis-acid labilisation needs a C-philic acid (Fe(II), Co(II): Avendano 2010)"},
        "KVMo_1to1_target": {"v": [0.02, 0.15], "v_best_estimate": 0.06, "w_plus_u": [0.0, 0.05],
                             "antisite_like_per_fu": "1e-3 to 1e-2 (cyanide-limited, dark, aprotic); more if [Mo(CN)6]3- decomposes",
                             "CN_flip_per_CN": "<=2e-3", "M0_muB_per_fu": [0.05, 0.45],
                             "basis": "no 1:1 phase exists; parent V1.37[Mo(CN)6] has v = 0.27 with K(crypt)+ (too big to template); "
                                      "K/Cs-templated V-Cr reach v = 0.05-0.08; NMF route reached v = 0.004 for Cr-Cr; "
                                      "[Mo(CN)6]3- is labile (4d3, 7-coordination, protic and photo decomposition)"}}
    out["T_C_vacancy_effect"] = "-3 % per 5 % [M(CN)6] vacancies (tc/README.md section 5), p_c(C sublattice) ~ 0.135"
    (HERE / "summary.json").write_text(json.dumps(out, indent=1))
    print("wrote", HERE / "summary.json")


if __name__ == "__main__":
    main()
