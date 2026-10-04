"""Write ledger/runs.csv: one row per included Quantum ESPRESSO input, with its settings and the final energy,
magnetisation and convergence of its output.

'converged' means: scf runs, the SCF converged; relax / vc-relax runs, the geometry optimisation also converged;
hybrid (HSE06) runs, the outer exact-exchange loop also converged (QE prints a '!!' total energy). It is 'n/a' for
nscf and bands runs, which are not self-consistent, and empty when no output is included.

    python tools/build_runs_index.py
"""
import csv
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
from qe_parse import parse_pw, resolve  # noqa: E402

DESC = json.load(open(ROOT / "ledger" / "run_folder_descriptions.json"))


def _converged(calc, r):
    if not r:
        return ""
    if calc in ("nscf", "bands"):
        return "n/a"
    if calc in ("relax", "vc-relax"):
        return bool(r["converged"] and r["bfgs_converged"])
    return bool(r["converged"])
rows = []
for mat in ("KV_Cr_CN6", "YBaMnFeO5"):
    for d in sorted((ROOT / "materials" / mat / "runs").iterdir()):
        if not d.is_dir():
            continue
        for inp in sorted(d.rglob("*.in")):
            try:
                r = parse_pw(resolve(inp.with_suffix(".out")))
                out_path = str(resolve(inp.with_suffix(".out")).relative_to(ROOT))
            except FileNotFoundError:
                r, out_path = {}, ""
            t = inp.read_text(errors="ignore")
            g = lambda pat: (m.group(1) if (m := re.search(pat, t)) else "")
            kp = t.split("K_POINTS")[1].splitlines()[1].strip() if "K_POINTS" in t else ""
            rows.append({
                "material": mat, "run_folder": d.name, "input": str(inp.relative_to(ROOT)), "output": out_path,
                "calculation": g(r"calculation\s*=\s*'(\S+)'"),
                "functional": "HSE06" if "input_dft = 'hse'" in t else ("PBE+U" if "HUBBARD" in t else ("PBE" if "&SYSTEM" in t else "")),
                "hubbard_U_eV": ";".join(f"{a}={b}" for a, b in re.findall(r"U\s+(\S+)-\S+\s+([\d.]+)", t)),
                "ecutwfc_Ry": g(r"ecutwfc\s*=\s*([\d.]+)"), "nat": g(r"nat\s*=\s*(\d+)"), "kpoints": kp,
                "q_mesh": "x".join(g(rf"nqx{i}\s*=\s*(\d+)") for i in (1, 2, 3)) if "nqx1" in t else "",
                "energy_eV": f"{r['energy_eV']:.6f}" if r.get("energy_eV") else "",
                "total_mag": r.get("total_mag", ""), "abs_mag": r.get("abs_mag", ""),
                "converged": _converged(g(r"calculation\s*=\s*'(\S+)'") or "scf", r),
                "exx_loop_converged": "" if r.get("exx_converged") is None else r.get("exx_converged"),
                "final_marker": r.get("final_marker", "") or "", "description": DESC.get(f"{mat}/{d.name}", ""),
            })
with open(ROOT / "ledger" / "runs.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(len(rows), "inputs indexed;", sum(1 for r in rows if r["exx_loop_converged"] is False), "hybrid runs whose exact-exchange loop did not converge")
