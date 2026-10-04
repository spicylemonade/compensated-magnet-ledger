"""Track L PBA defect Gate 1: chemical-potential references at the defect protocol (PseudoDojo SR NC, 110 Ry, no U).
  K_bcc : bcc K metal, vc-relax (k spacing 0.12 /A, Marzari-Vanderbilt 0.02 Ry), E per atom -> mu_K (K-rich limit)
  H2O   : isolated water, relax in a 14 A cubic box (Gamma), nspin 1 -> mu_H2O (gas, 0 K, no ZPE)
usage: python pba_refs.py   -> results/refs.json"""
import json
import os
import re

import numpy as np
from ase import Atoms
from ase.build import bulk
from qeutil import write_pw, run_pw, parse_pw, save_json, SCRATCH

NC = int(os.environ.get("NCORES", 16))
os.makedirs("results", exist_ok=True)
out = json.load(open("results/refs.json")) if os.path.exists("results/refs.json") else {}

if not out.get("K_bcc", {}).get("energy_eV"):
    at = bulk("K", "bcc", a=5.28)
    write_pw("K_bcc.in", at, calculation="vc-relax", pset="dojo_sr", ecutwfc=110.0, kspacing=0.12, smearing="mv",
             degauss=0.02, prefix="K_bcc", conv_thr=1e-9, cell={"press_conv_thr": 0.2})
    r = run_pw("K_bcc.in", ncores=NC, npool=min(NC, 16), timeout=4 * 3600)
    p = parse_pw(r["out"])
    txt = open(r["out"], errors="replace").read()
    vol = re.findall(r"unit-cell volume\s+=\s+([0-9.]+)", txt)
    out["K_bcc"] = {"energy_eV": p.get("energy_eV"), "per_atom_eV": p.get("energy_eV"), "pressure_kbar": p.get("pressure_kbar"),
                    "bfgs_converged": p.get("bfgs_converged"), "volume_bohr3": float(vol[-1]) if vol else None,
                    "a_A": (2 * float(vol[-1]) * 0.529177210903 ** 3) ** (1 / 3) if vol else None, "sec": r["seconds"], "errors": p.get("errors")}
    save_json(out, "results/refs.json")
    print("K_bcc", out["K_bcc"], flush=True)

if not out.get("H2O", {}).get("energy_eV"):
    L = 14.0
    th = np.radians(104.5 / 2)
    at = Atoms("OH2", positions=[[0, 0, 0], [0.97 * np.sin(th), 0, 0.97 * np.cos(th)], [-0.97 * np.sin(th), 0, 0.97 * np.cos(th)]],
               cell=[L, L, L], pbc=True)
    at.translate([L / 2, L / 2, L / 2])
    write_pw("H2O.in", at, calculation="relax", pset="dojo_sr", ecutwfc=110.0, kpts=(1, 1, 1), smearing="gaussian",
             degauss=0.001, prefix="H2O", conv_thr=1e-9)
    r = run_pw("H2O.in", ncores=NC, npool=1, timeout=4 * 3600)
    p = parse_pw(r["out"])
    out["H2O"] = {"energy_eV": p.get("energy_eV"), "bfgs_converged": p.get("bfgs_converged"), "box_A": L, "sec": r["seconds"],
                  "errors": p.get("errors")}
    save_json(out, "results/refs.json")
    print("H2O", out["H2O"], flush=True)
out["done"] = True
save_json(out, "results/refs.json")
