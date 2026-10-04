"""Check every number in the ledger against the raw files in this repository.

    python tools/verify.py                 # re-derive each claim from the raw pw.x outputs / data files
    python tools/verify.py --repro         # same checks, but on the independent re-runs in reproduce/results/
    python tools/verify.py --write-ledger  # regenerate ledger/claims.csv and the claim tables in LEDGER.md

Tiers
    A   re-derived here from raw Quantum ESPRESSO output files (seconds, only numpy needed)
    B   read from a recorded analysis file in this repository (the script that produced it is included)
    L2  re-run end to end by an included script (exchange fit + Monte Carlo, cluster-expansion Monte Carlo);
        the check here reads that script's output if you have run it, otherwise the recorded file
    E   experimental fact from the literature (not checkable by computation; cited)
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from qe_parse import final_cell, last_bands, parse_pw, resolve, spin_windows  # noqa: E402

K = "materials/KV_Cr_CN6/runs/"
Y = "materials/YBaMnFeO5/runs/"
KA = K + "A_pbeu_relax_scf_nscf_Ugrid/"

_cache: dict = {}


def _pw(rel):
    if ("pw", rel) not in _cache:
        _cache[("pw", rel)] = parse_pw(resolve(ROOT / rel))
    return _cache[("pw", rel)]


def _win(rel):
    if ("win", rel) not in _cache:
        _cache[("win", rel)] = spin_windows(last_bands(resolve(ROOT / rel)))
    return _cache[("win", rel)]


def _json(rel):
    return json.load(open(ROOT / rel))


# ---- extractors -------------------------------------------------------------------------------------------
def mag(rel):
    return _pw(rel)["total_mag"]


def site(rel, i):
    return _pw(rel)["site_moments"][i]


def win(rel, key):
    return _win(rel)[key]


def _hull(before):
    """YBaMnFeO5 energy above the hull (meV/atom) recomputed from the raw outputs of the hull runs."""
    from crosscheck_raw import HULL_FILE_WRITTEN, hull_distance
    return round(hull_distance(before=HULL_FILE_WRITTEN if before else None)[0], 2)


def dE_ion(ref, other, n):
    """(E_other - E_ref) per magnetic ion, meV."""
    return (_pw(other)["energy_eV"] - _pw(ref)["energy_eV"]) / n * 1000.0


def fcc_a(rel):
    c, _, _ = final_cell(resolve(ROOT / rel))
    return float(np.linalg.norm(c[0]) * np.sqrt(2.0))


def cell_len(rel, i):
    c, _, _ = final_cell(resolve(ROOT / rel))
    return float(np.linalg.norm(c[i]))


def ugrid(key, fn):
    vals = [win(KA + f"KVCr_{u}LCM_nscf.out", key) for u in ("", "U2_2_", "U2_4_", "U4_2_", "U4_4_")]
    return fn(vals)


def jget(rel, *path):
    d = _json(rel)
    for p in path:
        d = d[p]
    return d


def first_existing(*rels):
    for r in rels:
        if (ROOT / r).exists():
            return r
    return rels[-1]


def ugrid_unipolar():
    return all(win(KA + f"KVCr_{u}LCM_nscf.out", "unipolar") for u in ("", "U2_2_", "U2_4_", "U4_2_", "U4_4_"))


def text_has(rel, s):
    return s in (ROOT / rel).read_text()


CE = "materials/YBaMnFeO5/cation_order_cluster_expansion/"
EX = "materials/YBaMnFeO5/exchange_fit_and_Neel_T/"

# ---- the ledger ------------------------------------------------------------------------------------------
# expected = the value stated in the agents' dossiers (or 'new' if first computed for this ledger)
PBEU_K = "QE 7.5, PBE+U (ortho-atomic, U_V = U_Cr = 3 eV), PseudoDojo NC, 110 Ry"
HSE_K = "QE 7.5, HSE06 (no U), 90 Ry, q-mesh 2x2x2, k 4x4x4"
PBEU_Y = "QE 7.5, PBE+U (ortho-atomic, U_Mn = U_Fe = 4 eV), PseudoDojo NC"
HSE_Y = "QE 7.5, HSE06 (no U), 75 Ry on the 110-Ry geometry, q-mesh 2x2x1"

CLAIMS = [
    # ===================================== KV[Cr(CN)6] =====================================
    dict(id="K01", mat="KV[Cr(CN)6]", tier="A", what="Relaxed cubic lattice constant of the ideal anhydrous crystal (exp. 10.55 A, hydrated powder)",
         unit="A", expected=10.680, tol=0.002, method=PBEU_K + ", vc-relax", files=[KA + "KVCr_vcrelax.out"],
         f=lambda: fcc_a(KA + "KVCr_vcrelax.out"), repro=None),
    dict(id="K02", mat="KV[Cr(CN)6]", tier="A", what="Net spin moment of the ideal crystal (zero = fully compensated)",
         unit="Bohr mag/cell", expected=0.0, tol=0.01, method=PBEU_K + ", scf 5x5x5", files=[KA + "KVCr_LCM.out"],
         f=lambda: mag(KA + "KVCr_LCM.out"), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.repro.out", "mag")),
    dict(id="K03", mat="KV[Cr(CN)6]", tier="A", what="V moment (N-bound site), integrated in a small sphere: antiparallel to Cr",
         unit="Bohr mag", expected=-2.02, tol=0.02, method=PBEU_K, files=[KA + "KVCr_LCM.out"],
         f=lambda: site(KA + "KVCr_LCM.out", 1), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.repro.out", ("site", 1))),
    dict(id="K04", mat="KV[Cr(CN)6]", tier="A", what="Cr moment (C-bound site), integrated in a small sphere",
         unit="Bohr mag", expected=2.32, tol=0.02, method=PBEU_K, files=[KA + "KVCr_LCM.out"],
         f=lambda: site(KA + "KVCr_LCM.out", 2), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.repro.out", ("site", 2))),
    dict(id="K05", mat="KV[Cr(CN)6]", tier="A", what="Band gap", unit="eV", expected=1.95, tol=0.01,
         method=PBEU_K + ", nscf 8x8x8", files=[KA + "KVCr_LCM_nscf.out"], f=lambda: win(KA + "KVCr_LCM_nscf.out", "gap"),
         repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out", "gap")),
    dict(id="K06", mat="KV[Cr(CN)6]", tier="A", what="Both band edges in the same spin channel (unipolar)", unit="bool",
         expected=True, tol=None, method=PBEU_K + ", nscf 8x8x8", files=[KA + "KVCr_LCM_nscf.out"],
         f=lambda: win(KA + "KVCr_LCM_nscf.out", "unipolar"), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out", "unipolar")),
    dict(id="K07", mat="KV[Cr(CN)6]", tier="A", what="Hole spin window (top of valence band that is 100 % one spin)", unit="eV",
         expected=2.02, tol=0.01, method=PBEU_K + ", nscf 8x8x8", files=[KA + "KVCr_LCM_nscf.out"],
         f=lambda: win(KA + "KVCr_LCM_nscf.out", "win_VB"), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out", "win_VB")),
    dict(id="K08", mat="KV[Cr(CN)6]", tier="A", what="Electron spin window (bottom of conduction band that is 100 % one spin)", unit="eV",
         expected=1.15, tol=0.01, method=PBEU_K + ", nscf 8x8x8", files=[KA + "KVCr_LCM_nscf.out"],
         f=lambda: win(KA + "KVCr_LCM_nscf.out", "win_CB"), repro=("kvcr_pbeu", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out", "win_CB")),
    dict(id="K09", mat="KV[Cr(CN)6]", tier="A", what="Robust to the Hubbard U choice: smallest hole window over the U grid (U_V, U_Cr in {2,3,4} eV)",
         unit="eV", expected=1.39, tol=0.01, method=PBEU_K.replace("U_V = U_Cr = 3 eV", "5 (U_V, U_Cr) points"), files=[KA + "KVCr_*_LCM_nscf.out"],
         f=lambda: ugrid("win_VB", min), repro=None),
    dict(id="K10", mat="KV[Cr(CN)6]", tier="A", what="Smallest electron window over the U grid", unit="eV", expected=0.87, tol=0.01,
         method="as K09", files=[KA + "KVCr_*_LCM_nscf.out"], f=lambda: ugrid("win_CB", min), repro=None),
    dict(id="K11", mat="KV[Cr(CN)6]", tier="A", what="Unipolar at all five U points", unit="bool", expected=True, tol=None,
         method="as K09", files=[KA + "KVCr_*_LCM_nscf.out"], f=ugrid_unipolar, repro=None),
    dict(id="K12", mat="KV[Cr(CN)6]", tier="A", what="Compensated state is far below the ferromagnet: E(FM) - E(compensated) per magnetic ion",
         unit="meV", expected=152.2, tol=0.2, method=PBEU_K, files=[KA + "KVCr_LCM.out", KA + "KVCr_FM.out"],
         f=lambda: dE_ion(KA + "KVCr_LCM.out", KA + "KVCr_FM.out", 2),
         repro=("kvcr_pbeu", ("A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.repro.out", "A_pbeu_relax_scf_nscf_Ugrid/KVCr_FM.repro.out", 2), "dE")),
    dict(id="K13", mat="KV[Cr(CN)6]", tier="A", what="HSE06: net spin moment", unit="Bohr mag/cell", expected=0.0, tol=0.01, method=HSE_K,
         files=[K + "C_hse06_ideal_LCM/hse_lcm.out"], f=lambda: mag(K + "C_hse06_ideal_LCM/hse_lcm.out"),
         repro=("kvcr_hse", "C_hse06_ideal_LCM/hse_lcm.repro.out", "mag")),
    dict(id="K14", mat="KV[Cr(CN)6]", tier="A", what="HSE06: band gap", unit="eV", expected=2.09, tol=0.01, method=HSE_K,
         files=[K + "C_hse06_ideal_LCM/hse_lcm.out"], f=lambda: win(K + "C_hse06_ideal_LCM/hse_lcm.out", "gap"),
         repro=("kvcr_hse", "C_hse06_ideal_LCM/hse_lcm.repro.out", "gap")),
    dict(id="K15", mat="KV[Cr(CN)6]", tier="A", what="HSE06: unipolar band edges", unit="bool", expected=True, tol=None, method=HSE_K,
         files=[K + "C_hse06_ideal_LCM/hse_lcm.out"], f=lambda: win(K + "C_hse06_ideal_LCM/hse_lcm.out", "unipolar"),
         repro=("kvcr_hse", "C_hse06_ideal_LCM/hse_lcm.repro.out", "unipolar")),
    dict(id="K16", mat="KV[Cr(CN)6]", tier="A", what="HSE06: hole spin window", unit="eV", expected=2.64, tol=0.01, method=HSE_K,
         files=[K + "C_hse06_ideal_LCM/hse_lcm.out"], f=lambda: win(K + "C_hse06_ideal_LCM/hse_lcm.out", "win_VB"),
         repro=("kvcr_hse", "C_hse06_ideal_LCM/hse_lcm.repro.out", "win_VB")),
    dict(id="K17", mat="KV[Cr(CN)6]", tier="A", what="HSE06: electron spin window", unit="eV", expected=1.57, tol=0.01, method=HSE_K,
         files=[K + "C_hse06_ideal_LCM/hse_lcm.out"], f=lambda: win(K + "C_hse06_ideal_LCM/hse_lcm.out", "win_CB"),
         repro=("kvcr_hse", "C_hse06_ideal_LCM/hse_lcm.repro.out", "win_CB")),
    dict(id="K18", mat="KV[Cr(CN)6]", tier="A", what="HSE06: E(FM, moment fixed at 6) - E(compensated) per magnetic ion", unit="meV",
         expected=181.1, tol=0.2, method=HSE_K, files=[K + "C_hse06_ideal_LCM/hse_lcm.out", K + "D_hse06_ideal_FM_fixed_moment/KVCr_hse_FMfix.out"],
         f=lambda: dE_ion(K + "C_hse06_ideal_LCM/hse_lcm.out", K + "D_hse06_ideal_FM_fixed_moment/KVCr_hse_FMfix.out", 2), repro=None),
    dict(id="K19", mat="KV[Cr(CN)6]", tier="A", what="HSE06, dihydrate KV[Cr(CN)6]*2H2O (one water arrangement; exact-exchange loop stopped before full convergence): hole window",
         unit="eV", expected=2.43, tol=0.01, method=HSE_K.replace("q-mesh 2x2x2, k 4x4x4", "q-mesh 1, k-spacing 0.30/A"),
         files=[K + "F_hse06_dihydrate_LCM/hse_lcm.out"], f=lambda: win(K + "F_hse06_dihydrate_LCM/hse_lcm.out", "win_VB"), repro=None),
    dict(id="K20", mat="KV[Cr(CN)6]", tier="A", what="HSE06, dihydrate: electron window (still unipolar)", unit="eV", expected=1.42, tol=0.01,
         method="as K19", files=[K + "F_hse06_dihydrate_LCM/hse_lcm.out"], f=lambda: win(K + "F_hse06_dihydrate_LCM/hse_lcm.out", "win_CB"), repro=None),
    dict(id="K21", mat="KV[Cr(CN)6]", tier="A", what="PBE+U, same dihydrate: hole window collapses (methods disagree about water)", unit="eV",
         expected=0.93, tol=0.01, method=PBEU_K + ", nscf 6x6x6", files=[K + "G_pbeu_dihydrate/KVCr_2H2O_LCM_nscf.out"],
         f=lambda: win(K + "G_pbeu_dihydrate/KVCr_2H2O_LCM_nscf.out", "win_VB"), repro=None),
    dict(id="K22", mat="KV[Cr(CN)6]", tier="A", what="HSE06, water-filled [Cr(CN)6] vacancy (65-atom cell): net moment = -3 (composition, not symmetry, sets the moment)",
         unit="Bohr mag/cell", expected=-3.0, tol=0.01, method=HSE_K.replace("q-mesh 2x2x2, k 4x4x4", "q-mesh 1"),
         files=[K + "H_hse06_vacancy_with_water/hse_lcm.out"], f=lambda: mag(K + "H_hse06_vacancy_with_water/hse_lcm.out"), repro=None),
    dict(id="K23", mat="KV[Cr(CN)6]", tier="A", what="HSE06, water-filled vacancy: electron window (in-cell)", unit="eV", expected=0.74, tol=0.01,
         method="as K22", files=[K + "H_hse06_vacancy_with_water/hse_lcm.out"], f=lambda: win(K + "H_hse06_vacancy_with_water/hse_lcm.out", "win_CB"), repro=None),
    dict(id="K24", mat="KV[Cr(CN)6]", tier="A", what="HSE06, water-filled vacancy: hole window (in-cell)", unit="eV", expected=2.80, tol=0.01,
         method="as K22", files=[K + "H_hse06_vacancy_with_water/hse_lcm.out"], f=lambda: win(K + "H_hse06_vacancy_with_water/hse_lcm.out", "win_VB"), repro=None),
    dict(id="K25", mat="KV[Cr(CN)6]", tier="A", what="PBE+U, V/Cr antisite pair (60-atom cell): moment stays 0 ...", unit="Bohr mag/cell",
         expected=0.0, tol=0.01, method=PBEU_K, files=[K + "K_pbeu_defect_antisite_pair/KVCr_AS_lcm.out"],
         f=lambda: mag(K + "K_pbeu_defect_antisite_pair/KVCr_AS_lcm.out"), repro=None),
    dict(id="K26", mat="KV[Cr(CN)6]", tier="A", what="... but the hole window shrinks to (in-cell; 0.13 eV after alignment to the host)", unit="eV",
         expected=0.23, tol=0.01, method=PBEU_K + ", nscf 4x4x4", files=[K + "K_pbeu_defect_antisite_pair/KVCr_AS_lcm_nscf.out"],
         f=lambda: win(K + "K_pbeu_defect_antisite_pair/KVCr_AS_lcm_nscf.out", "win_VB"), repro=None),
    dict(id="K27", mat="KV[Cr(CN)6]", tier="A", what="Control, same-element Cr[Cr(CN)6] (HSE06): edges in opposite spins (bipolar) ...", unit="bool",
         expected=False, tol=None, method=HSE_K, files=[K + "Q_hse06_control_CrCr_same_element/hse_lcm.out"],
         f=lambda: win(K + "Q_hse06_control_CrCr_same_element/hse_lcm.out", "unipolar"), repro=None),
    dict(id="K28", mat="KV[Cr(CN)6]", tier="A", what="... with a tiny hole window (why two different metals matter)", unit="eV", expected=0.10, tol=0.01,
         method=HSE_K, files=[K + "Q_hse06_control_CrCr_same_element/hse_lcm.out"],
         f=lambda: win(K + "Q_hse06_control_CrCr_same_element/hse_lcm.out", "win_VB"), repro=None),
    dict(id="K29", mat="KV[Cr(CN)6]", tier="B", what="Spin space group of the ideal crystal (findspingroup): compensated ferrimagnet, Zeeman-type (s-wave) splitting",
         unit="label", expected="216.216.1.1.L", tol=None, method="findspingroup on the relaxed P1 cell",
         files=["materials/KV_Cr_CN6/analysis_data/pba_spingroup_check.txt"],
         f=lambda: "216.216.1.1.L" if text_has("materials/KV_Cr_CN6/analysis_data/pba_spingroup_check.txt", '"fsg_ssg": "216.216.1.1.L"') else "missing",
         repro=None),
    dict(id="K30", mat="KV[Cr(CN)6]", tier="E", what="Measured magnetic ordering temperature of KV[Cr(CN)6]*2H2O (Holmes & Girolami, JACS 121, 5593 (1999))",
         unit="K", expected=376, tol=None, method="experiment (SQUID magnetometry, hydrated powder)", files=["literature"], f=None, repro=None),
    dict(id="K31", mat="KV[Cr(CN)6]", tier="E", what="Measured saturation moment at 5 K (expected 0 for perfect stoichiometry)", unit="Bohr mag/f.u.",
         expected=0.125, tol=None, method="experiment (Holmes & Girolami 1999, via Verdaguer & Girolami 2005 review)", files=["literature"], f=None, repro=None),
    # ===================================== YBaMnFeO5 =====================================
    dict(id="Y01", mat="YBaMnFeO5", tier="A", what="Relaxed cell, rock-salt-ordered P4/n: a", unit="A", expected=5.665, tol=0.001,
         method=PBEU_Y + ", 110 Ry, vc-relax, G-type", files=[Y + "A_pbeu_vcrelax_110Ry/vc110.out"],
         f=lambda: cell_len(Y + "A_pbeu_vcrelax_110Ry/vc110.out", 0), repro=None),
    dict(id="Y02", mat="YBaMnFeO5", tier="A", what="Relaxed cell: c", unit="A", expected=7.689, tol=0.001, method="as Y01",
         files=[Y + "A_pbeu_vcrelax_110Ry/vc110.out"], f=lambda: cell_len(Y + "A_pbeu_vcrelax_110Ry/vc110.out", 2), repro=None),
    dict(id="Y03", mat="YBaMnFeO5", tier="A", what="Net spin moment of the G-type compensated state", unit="Bohr mag/cell", expected=0.0, tol=0.01,
         method=PBEU_Y + ", 90 Ry, 36-atom cell", files=[Y + "H_pbeu_stacking_vs_U/U44_G/single.out"],
         f=lambda: mag(Y + "H_pbeu_stacking_vs_U/U44_G/single.out"), repro=("ybmfo_stack_U44", "U44_G/single.repro.out", "mag")),
    dict(id="Y04", mat="YBaMnFeO5", tier="A", what="Mn moment (sphere-integrated); Fe moment is about -3.7, antiparallel", unit="Bohr mag",
         expected=4.00, tol=0.02, method="as Y03", files=[Y + "H_pbeu_stacking_vs_U/U44_G/single.out"],
         f=lambda: site(Y + "H_pbeu_stacking_vs_U/U44_G/single.out", 8), repro=("ybmfo_stack_U44", "U44_G/single.repro.out", ("site", 8))),
    dict(id="Y05", mat="YBaMnFeO5", tier="A", what="Closest competing magnetic order (spin flip across the Y layer) above the ground state, U = 4 eV",
         unit="meV/mag. ion", expected=8.9, tol=0.1, method="as Y03",
         files=[Y + "H_pbeu_stacking_vs_U/U44_G/single.out", Y + "H_pbeu_stacking_vs_U/U44_Yflip/single.out"],
         f=lambda: dE_ion(Y + "H_pbeu_stacking_vs_U/U44_G/single.out", Y + "H_pbeu_stacking_vs_U/U44_Yflip/single.out", 8),
         repro=("ybmfo_stack_U44", ("U44_G/single.repro.out", "U44_Yflip/single.repro.out", 8), "dE")),
    dict(id="Y06", mat="YBaMnFeO5", tier="A", what="Same competitor at U = 0 eV", unit="meV/mag. ion", expected=16.4, tol=0.1, method="as Y03, U = 0",
         files=[Y + "H_pbeu_stacking_vs_U/U00_*"], f=lambda: dE_ion(Y + "H_pbeu_stacking_vs_U/U00_G/single.out", Y + "H_pbeu_stacking_vs_U/U00_Yflip/single.out", 8), repro=None),
    dict(id="Y07", mat="YBaMnFeO5", tier="A", what="Same competitor at U = 6 eV (smallest margin in the U scan)", unit="meV/mag. ion", expected=7.0, tol=0.1,
         method="as Y03, U = 6", files=[Y + "H_pbeu_stacking_vs_U/U66_*"],
         f=lambda: dE_ion(Y + "H_pbeu_stacking_vs_U/U66_G/single.out", Y + "H_pbeu_stacking_vs_U/U66_Yflip/single.out", 8), repro=None),
    dict(id="Y08", mat="YBaMnFeO5", tier="A", what="Ferromagnet above the ground state, U = 4 eV (72-atom cell, 75 Ry)", unit="meV/mag. ion",
         expected=199.8, tol=0.2, method=PBEU_Y + ", 75 Ry, 72-atom cell (new number; dossier range +160 to +330 over U = 0-6 eV)",
         files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G.out", Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_FM.out"],
         f=lambda: dE_ion(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G.out", Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_FM.out", 16), repro=None),
    dict(id="Y09", mat="YBaMnFeO5", tier="A", what="HSE06: competitor (Y-layer flip) above the ground state", unit="meV/mag. ion", expected=8.5, tol=0.1,
         method=HSE_Y, files=[Y + "B_hse06_G_type_LCM/hse_lcm.out", Y + "C_hse06_Yflip_competitor/hse_lcm.out"],
         f=lambda: dE_ion(Y + "B_hse06_G_type_LCM/hse_lcm.out", Y + "C_hse06_Yflip_competitor/hse_lcm.out", 8), repro=None),
    dict(id="Y10", mat="YBaMnFeO5", tier="A", what="HSE06: band gap", unit="eV", expected=2.35, tol=0.01, method=HSE_Y,
         files=[Y + "B_hse06_G_type_LCM/hse_lcm.out"], f=lambda: win(Y + "B_hse06_G_type_LCM/hse_lcm.out", "gap"), repro=None),
    dict(id="Y11", mat="YBaMnFeO5", tier="A", what="HSE06: unipolar band edges", unit="bool", expected=True, tol=None, method=HSE_Y,
         files=[Y + "B_hse06_G_type_LCM/hse_lcm.out"], f=lambda: win(Y + "B_hse06_G_type_LCM/hse_lcm.out", "unipolar"), repro=None),
    dict(id="Y12", mat="YBaMnFeO5", tier="A", what="HSE06: hole spin window", unit="eV", expected=1.00, tol=0.01, method=HSE_Y,
         files=[Y + "B_hse06_G_type_LCM/hse_lcm.out"], f=lambda: win(Y + "B_hse06_G_type_LCM/hse_lcm.out", "win_VB"), repro=None),
    dict(id="Y13", mat="YBaMnFeO5", tier="A", what="HSE06: electron spin window", unit="eV", expected=1.40, tol=0.01, method=HSE_Y,
         files=[Y + "B_hse06_G_type_LCM/hse_lcm.out"], f=lambda: win(Y + "B_hse06_G_type_LCM/hse_lcm.out", "win_CB"), repro=None),
    dict(id="Y14", mat="YBaMnFeO5", tier="A", what="PBE+U band gap (72-atom cell)", unit="eV", expected=1.40, tol=0.01,
         method=PBEU_Y + ", 75 Ry on the 110-Ry geometry, nscf", files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out"],
         f=lambda: win(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out", "gap"), repro=None),
    dict(id="Y15", mat="YBaMnFeO5", tier="A", what="PBE+U windows: hole", unit="eV", expected=0.45, tol=0.01, method="as Y14",
         files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out"],
         f=lambda: win(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out", "win_VB"), repro=None),
    dict(id="Y16", mat="YBaMnFeO5", tier="A", what="PBE+U windows: electron", unit="eV", expected=1.43, tol=0.01, method="as Y14",
         files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out"],
         f=lambda: win(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_G_nscf.out", "win_CB"), repro=None),
    dict(id="Y17", mat="YBaMnFeO5", tier="A", what="Worst-case spin disorder (2 of 16 spins flipped): electron edge keeps its spin, window", unit="eV",
         expected=0.155, tol=0.005, method="as Y14", files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out"],
         f=lambda: win(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out", "win_CB"), repro=None),
    dict(id="Y18", mat="YBaMnFeO5", tier="A", what="... while the hole edge switches to the other spin (holes need good magnetic order)", unit="bool",
         expected=False, tol=None, method="as Y14", files=[Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out"],
         f=lambda: win(Y + "D_pbeu_72atom_G_FM_and_spin_flips/YBMFO110_snap_flip2far_nscf.out", "unipolar"), repro=None),
    dict(id="Y19", mat="YBaMnFeO5", tier="A", what="One neighbouring Mn/Fe swap (12.5 % of sites, 72-atom cell): gap of the other spin channel collapses (dense grid)",
         unit="eV", expected=0.01, tol=0.005, method=PBEU_Y + ", 75 Ry, relaxed defect cell, nscf",
         files=[Y + "F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite_nscf.out"],
         f=lambda: win(Y + "F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite_nscf.out", "gap_dn"), repro=None),
    dict(id="Y20", mat="YBaMnFeO5", tier="A", what="... same, on the coarser self-consistent grid", unit="eV", expected=0.30, tol=0.005, method="as Y19",
         files=[Y + "F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite.out"],
         f=lambda: win(Y + "F_pbeu_antisite_pair_nearest/YBMFO_antisite_swap_NN_Gsite.out", "gap_dn"), repro=None),
    dict(id="Y21", mat="YBaMnFeO5", tier="L2", what="Exchange fit from 25 raw energies: in-plane Mn-O-Fe coupling J_ip", unit="meV",
         expected=-32.6, tol=0.1, method="least squares, 10 pair classes, 22 'clean' configurations (rerun_fit_and_mc.py)",
         files=[EX + "rerun_fit_and_mc.py", Y + "I_pbeu_exchange_fit_configs/", Y + "H_pbeu_stacking_vs_U/U44_*"],
         f=lambda: jget(first_existing(EX + "rerun_result.json", EX + "rerun_result_quick.json"), "J_meV", "ip")
         if (ROOT / EX / "rerun_result.json").exists() or (ROOT / EX / "rerun_result_quick.json").exists()
         else dict(zip(jget(EX + "jfit110/fits.json", "clean_M<=10.5_10", "classes"), jget(EX + "jfit110/fits.json", "clean_M<=10.5_10", "J_meV")))["ip"],
         repro=None),
    dict(id="Y22", mat="YBaMnFeO5", tier="L2", what="Neel temperature, classical Monte Carlo on the fitted model (Binder crossings, L = 6/8/10)", unit="K",
         expected=417, tol=15, method="rerun_fit_and_mc.py (full mode)", files=[EX + "rerun_fit_and_mc.py", EX + "jfit110/mc_clean_Mle10.5_10.json"],
         f=lambda: float(np.mean(jget(EX + "rerun_result.json", "binder_crossings_K"))) if (ROOT / EX / "rerun_result.json").exists()
         else float(np.mean(jget(EX + "jfit110/mc_clean_Mle10.5_10.json", "binder_first_crossings"))), repro=None),
    dict(id="Y23", mat="YBaMnFeO5", tier="B", what="Calibrated Neel temperature: Y22 x (measured/simulated ratio for YFeO3 with the same protocol, 1.17-1.19)",
         unit="K", expected=490, tol=15, method="T_MC x mean f_MC(YFeO3); see Neel_T_calibration/",
         files=["materials/YBaMnFeO5/analysis_data/calibration_summary.json"],
         f=lambda: CLAIMS_BY_ID["Y22"]["f"]() * float(np.mean(jget("materials/YBaMnFeO5/analysis_data/calibration_summary.json",
                                                                    "YBMFO_calibration", "per_calibrant", "YFeO3_best", "f_MC"))),
         repro=None),
    dict(id="Y24", mat="YBaMnFeO5", tier="A", what="Energy above the 0 K convex hull, as the agents computed it: the 28 competing phases that had finished by then (superseded by Y24b)",
         unit="meV/atom", expected=2.6, tol=0.05,
         method=PBEU_Y + ", 110 Ry vc-relax of every phase; lowest-energy mix at fixed composition (tools/crosscheck_raw.py)",
         files=[Y + "K_hull_competing_phases_raw/", "materials/YBaMnFeO5/analysis_data/hull_QE110.json"],
         f=lambda: _hull(before=True), repro=None),
    dict(id="Y24b", mat="YBaMnFeO5", tier="A", what="Energy above the 0 K convex hull with all 33 competing phases that completed (BaFe2O4 and Ba6Y2Fe4O15 finished later and lower the hull; Ba2Fe2O5 never finished)",
         unit="meV/atom", expected=13.7, tol=0.05, method="as Y24",
         files=[Y + "K_hull_competing_phases_raw/"], f=lambda: _hull(before=False), repro=None),
    dict(id="Y25", mat="YBaMnFeO5", tier="L2", what="DECISIVE: Mn/Fe order-disorder temperature from the paramagnetic cluster expansion (C peak, L = 16)",
         unit="K", expected=915, tol=60, method="canonical MC of a fitted cluster expansion (21 arrangements, 96 DFT runs); rerun_torder.sh",
         files=[CE + "ce2_YBMFO.json", CE + "torder.py", CE + "rerun_torder.sh"],
         f=lambda: jget(first_existing(CE + "rerun_torder_YBMFO_E_PM_best_bw_bayes.json", CE + "torder_YBMFO_E_PM_best_bw_bayes.json"), "T_order_K"),
         repro=None),
    dict(id="Y26", mat="YBaMnFeO5", tier="B", what="Headline order-disorder temperature with uncertainty (bootstrap + model variants)", unit="K",
         expected="950 (+250/-150)", tol=None, method="ord2_summary.json", files=[CE + "ord2_summary.json"],
         f=lambda: "950 (+250/-150)" if "950 (+250/-150)" in jget(CE + "ord2_summary.json", "conclusions", "YBMFO_T_order_K", "headline") else "missing",
         repro=None),
    dict(id="Y27", mat="YBaMnFeO5", tier="E", what="Typical synthesis/anneal window for these layered perovskites (900-1300 C); B-site exchange inferred to freeze below ~1150 K (YBaCuFeO5 benchmark: Morin et al., Nat. Commun. 7, 13758 (2016))",
         unit="K", expected="1173-1573", tol=None, method="literature-based estimate (the 1150 K freeze-out is an analogy from a Cu/Fe compound)",
         files=["literature"], f=None, repro=None),
]


CLAIMS_BY_ID = {c["id"]: c for c in CLAIMS}


def _val_repro(c):
    case, spec, kind = c["repro"]
    base = f"reproduce/results/{case}/"
    if kind == "dE":
        a, b, n = spec
        return (parse_pw(resolve(ROOT / (base + b)))["energy_eV"] - parse_pw(resolve(ROOT / (base + a)))["energy_eV"]) / n * 1000.0
    rel = base + spec
    if kind == "mag":
        return parse_pw(resolve(ROOT / rel))["total_mag"]
    if isinstance(kind, tuple) and kind[0] == "site":
        return parse_pw(resolve(ROOT / rel))["site_moments"][kind[1]]
    return spin_windows(last_bands(resolve(ROOT / rel)))[kind]


def _status(val, c):
    if val is None:
        return "n/a"
    exp, tol = c["expected"], c["tol"]
    if isinstance(exp, bool) or isinstance(val, bool):
        return "PASS" if bool(val) == bool(exp) else "FAIL"
    if isinstance(exp, str):
        return "PASS" if str(val) == exp else "FAIL"
    if tol is None:
        return "n/a"
    return "PASS" if abs(float(val) - float(exp)) <= tol else "FAIL"


def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        return f"{v:.3f}" if abs(v) < 100 else f"{v:.1f}"
    return str(v)


def run(repro=False):
    rows = []
    for c in CLAIMS:
        if repro and not c["repro"]:
            continue
        try:
            if repro:
                val = _val_repro(c)
            else:
                val = c["f"]() if c["f"] else None
            st = _status(val, c)
        except FileNotFoundError as e:
            val, st = None, f"MISSING ({pathlib.Path(str(e)).name})"
        except Exception as e:  # pragma: no cover
            val, st = None, f"ERROR {e!r}"
        rows.append((c, val, st))
    return rows


def print_table(rows, title):
    print(f"\n{title}\n" + "=" * len(title))
    print(f"{'id':4s} {'tier':4s} {'status':8s} {'computed':>12s} {'recorded':>14s}  {'unit':16s} claim")
    for c, val, st in rows:
        print(f"{c['id']:4s} {c['tier']:4s} {st:8s} {_fmt(val):>12s} {_fmt(c['expected']):>14s}  {c['unit']:16s} {c['what'][:90]}")
    n = {s: sum(1 for _, _, x in rows if x == s) for s in ("PASS", "FAIL")}
    miss = sum(1 for _, _, x in rows if x.startswith(("MISSING", "ERROR")))
    other = len(rows) - n["PASS"] - n["FAIL"] - miss
    print(f"\n{n['PASS']} pass, {n['FAIL']} fail, {miss} missing input files, {other} not computable here (experiment / recorded-only)")
    return n["FAIL"] == 0 and miss == 0


def write_ledger(rows):
    out = ROOT / "ledger" / "claims.csv"
    out.parent.mkdir(exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "material", "tier", "claim", "recorded_value", "unit", "recomputed_value", "status", "tolerance", "method", "files"])
        for c, val, st in rows:
            w.writerow([c["id"], c["mat"], c["tier"], c["what"], _fmt(c["expected"]), c["unit"], _fmt(val), st,
                        "" if c["tol"] is None else c["tol"], c["method"], "; ".join(c["files"])])
    md = []
    for mat in ("KV[Cr(CN)6]", "YBaMnFeO5"):
        md.append(f"\n#### {mat}\n\n| id | tier | claim | recorded | recomputed | status | raw files |\n|---|---|---|---|---|---|---|")
        for c, val, st in rows:
            if c["mat"] != mat:
                continue
            files = "<br>".join(f"`{f.replace('materials/', '')}`" for f in c["files"])
            md.append(f"| {c['id']} | {c['tier']} | {c['what']} ({c['unit']}) | {_fmt(c['expected'])} | {_fmt(val)} | {st} | {files} |")
    led = ROOT / "LEDGER.md"
    txt = led.read_text() if led.exists() else ""
    a, b = "<!-- CLAIMS-TABLE:BEGIN -->", "<!-- CLAIMS-TABLE:END -->"
    block = a + "\n" + "\n".join(md) + "\n\n" + b
    if a in txt and b in txt:
        txt = txt.split(a)[0] + block + txt.split(b)[1]
    else:
        txt += "\n" + block + "\n"
    led.write_text(txt)
    print(f"wrote {out.relative_to(ROOT)} and the claim tables in LEDGER.md")


if __name__ == "__main__":
    if "--repro" in sys.argv:
        ok = print_table(run(repro=True), "Independent re-runs (reproduce/results) vs recorded values")
    else:
        rows = run()
        ok = print_table(rows, "Ledger claims re-derived from the raw files in this repository")
        if "--write-ledger" in sys.argv:
            write_ledger(rows)
    sys.exit(0 if ok else 1)
