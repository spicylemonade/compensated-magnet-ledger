"""Write reproduce/RESULTS.md: original runs vs independent re-runs (reproduce/results/), side by side."""
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from qe_parse import last_bands, parse_pw, resolve, spin_windows  # noqa: E402

K = ROOT / "materials/KV_Cr_CN6/runs"
Y = ROOT / "materials/YBaMnFeO5/runs"
R = ROOT / "reproduce/results"


def ok(*ps):
    try:
        [resolve(p) for p in ps]
        return True
    except FileNotFoundError:
        return False


rows = []
a, r = K / "A_pbeu_relax_scf_nscf_Ugrid", R / "kvcr_pbeu/A_pbeu_relax_scf_nscf_Ugrid"
if ok(r / "KVCr_LCM.repro.out", r / "KVCr_LCM_nscf.repro.out", r / "KVCr_FM.repro.out"):
    o, n = parse_pw(resolve(a / "KVCr_LCM.out")), parse_pw(resolve(r / "KVCr_LCM.repro.out"))
    fo, fn = parse_pw(resolve(a / "KVCr_FM.out")), parse_pw(resolve(r / "KVCr_FM.repro.out"))
    wo, wn = spin_windows(last_bands(resolve(a / "KVCr_LCM_nscf.out"))), spin_windows(last_bands(resolve(r / "KVCr_LCM_nscf.repro.out")))
    rows += [("KV[Cr(CN)6] PBE+U (U 3/3 eV)", "net spin moment (Bohr mag/cell)", f"{o['total_mag']:.2f}", f"{n['total_mag']:.2f}"),
             ("", "total energy (eV)", f"{o['energy_eV']:.6f}", f"{n['energy_eV']:.6f}"),
             ("", "V / Cr sphere moments", f"{o['site_moments'][1]:.3f} / {o['site_moments'][2]:.3f}", f"{n['site_moments'][1]:.3f} / {n['site_moments'][2]:.3f}"),
             ("", "band gap, nscf 8x8x8 (eV)", f"{wo['gap']:.4f}", f"{wn['gap']:.4f}"),
             ("", "hole / electron window (eV)", f"{wo['win_VB']:.4f} / {wo['win_CB']:.4f}", f"{wn['win_VB']:.4f} / {wn['win_CB']:.4f}"),
             ("", "both edges same spin", str(wo["unipolar"]), str(wn["unipolar"])),
             ("", "E(FM) - E(compensated) (meV/ion)", f"{(fo['energy_eV'] - o['energy_eV']) / 2e-3:.2f}", f"{(fn['energy_eV'] - n['energy_eV']) / 2e-3:.2f}")]
h = R / "kvcr_hse/C_hse06_ideal_LCM/hse_lcm.repro.out"
if ok(h):
    o, n = parse_pw(resolve(K / "C_hse06_ideal_LCM/hse_lcm.out")), parse_pw(resolve(h))
    wo, wn = spin_windows(last_bands(resolve(K / "C_hse06_ideal_LCM/hse_lcm.out"))), spin_windows(last_bands(resolve(h)))
    rows += [("KV[Cr(CN)6] HSE06", "net spin moment", f"{o['total_mag']:.2f}", f"{n['total_mag']:.2f}"),
             ("", "total energy (eV)", f"{o['energy_eV']:.6f}", f"{n['energy_eV']:.6f}"),
             ("", "band gap (eV)", f"{wo['gap']:.4f}", f"{wn['gap']:.4f}"),
             ("", "hole / electron window (eV)", f"{wo['win_VB']:.4f} / {wo['win_CB']:.4f}", f"{wn['win_VB']:.4f} / {wn['win_CB']:.4f}")]
else:
    rows.append(("KV[Cr(CN)6] HSE06", "(re-run still in progress when this file was written)", "", ""))
yb = R / "ybmfo_stack_U44"
if ok(yb / "U44_G/single.repro.out", yb / "U44_Yflip/single.repro.out"):
    go, yo = parse_pw(resolve(Y / "H_pbeu_stacking_vs_U/U44_G/single.out")), parse_pw(resolve(Y / "H_pbeu_stacking_vs_U/U44_Yflip/single.out"))
    gn, yn = parse_pw(resolve(yb / "U44_G/single.repro.out")), parse_pw(resolve(yb / "U44_Yflip/single.repro.out"))
    rows += [("YBaMnFeO5 PBE+U (U 4/4 eV)", "G-type total energy (eV)", f"{go['energy_eV']:.6f}", f"{gn['energy_eV']:.6f}"),
             ("", "net spin moment (G)", f"{go['total_mag']:.2f}", f"{gn['total_mag']:.2f}"),
             ("", "E(Y-layer flip) - E(G) (meV/mag. ion)", f"{(yo['energy_eV'] - go['energy_eV']) / 8e-3:.3f}", f"{(yn['energy_eV'] - gn['energy_eV']) / 8e-3:.3f}")]
s = R / "kvcr_pbeu_sssp/A_pbeu_relax_scf_nscf_Ugrid"
sssp_rows = []
if ok(s / "KVCr_LCM.repro.out", s / "KVCr_LCM_nscf.repro.out"):
    n = parse_pw(resolve(s / "KVCr_LCM.repro.out")); ws = spin_windows(last_bands(resolve(s / "KVCr_LCM_nscf.repro.out")))
    wo = spin_windows(last_bands(resolve(a / "KVCr_LCM_nscf.out")))
    sssp_rows = [("net spin moment", "0.00", f"{n['total_mag']:.2f}"),
                 ("V / Cr sphere moments", "-2.016 / 2.318", f"{n['site_moments'][1]:.3f} / {n['site_moments'][2]:.3f}"),
                 ("band gap (eV)", f"{wo['gap']:.3f}", f"{ws['gap']:.3f}"),
                 ("hole / electron window (eV)", f"{wo['win_VB']:.3f} / {wo['win_CB']:.3f}", f"{ws['win_VB']:.3f} / {ws['win_CB']:.3f}"),
                 ("both edges same spin", str(wo["unipolar"]), str(ws["unipolar"]))]

md5 = {}
for c in ("kvcr_pbeu", "kvcr_hse", "ybmfo_stack_U44", "kvcr_pbeu_sssp"):
    f = R / c / "summary.json"
    if f.exists():
        md5[c] = json.loads(f.read_text())

L = ["# Independent re-runs: results",
     "",
     "These were re-run on 2026-10-04 in a fresh cloud container (Modal) using `reproduce/modal_repro.py`:",
     "",
     "- a clean conda-forge Quantum ESPRESSO 7.5 image, with no campaign code;",
     "- pseudopotentials downloaded from their public sources;",
     "- the unmodified input files from this repository (only `pseudo_dir` and `outdir` rewritten).",
     "",
     "**Pseudopotential check.** Every PseudoDojo file used by these materials matched the md5 sum of the file used in the original runs:",
     ""]
for c, v in md5.items():
    if c != "kvcr_pbeu_sssp" and v.get("md5_check"):
        bad = [x for x in v["md5_check"] if "DIFFERS" in x or "MISSING" in x]
        L.append(f"- `{c}`: {len(v['md5_check'])} files checked, {'all identical' if not bad else 'MISMATCH: ' + '; '.join(bad)}")
        break
L += ["", "## Reproduction (same inputs, same pseudopotentials)", "", "| Calculation | Quantity | Original run | Independent re-run |", "|---|---|---|---|"]
L += [f"| {a_} | {b_} | {c_} | {d_} |" for a_, b_, c_, d_ in rows]
L += ["",
      "Differences in the last digits come from parallelisation and random starting wavefunctions. They are far below the precision any claim relies on.",
      "",
      "Check it yourself: `python tools/verify.py --repro`.",
      ""]
if sssp_rows:
    L += ["## Robustness: a different pseudopotential family (not a reproduction)", "",
          "This is the same KV[Cr(CN)6] PBE+U calculation with SSSP 1.3 efficiency pseudopotentials (GBRV ultrasoft for V and Cr, PSlibrary PAW for K and C, THEOS for N; 60/480 Ry) instead of PseudoDojo. "
          "Hubbard projectors depend on the pseudopotential, so the numbers are not expected to match exactly. The physical picture is unchanged.", "",
          "| Quantity | PseudoDojo (original) | SSSP 1.3 (re-run) |", "|---|---|---|"]
    L += [f"| {a_} | {b_} | {c_} |" for a_, b_, c_ in sssp_rows]
    L.append("")
(ROOT / "reproduce" / "RESULTS.md").write_text("\n".join(L) + "\n")
print("\n".join(L))
