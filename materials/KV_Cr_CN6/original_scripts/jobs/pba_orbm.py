"""Track L PBA Gate 0, SOC follow-up (best effort): cell ORBITAL magnetization of the LCM state of KV[Mo(CN)6].
Runs after `lcm_deep.py <cfg> soc` (which gives the spin residue |M_spin| for m||z and m||x, and the MAE) in the same job.
  1. noncollinear + SOC PBE+U SCF (dojo_fr, moments along z, same settings as lcm_deep soc) with the scratch kept
  2. nscf with lorbm = .true. (modern-theory / Kubo orbital magnetization, NC pseudopotentials only) on the full
     uniform k grid (nosym, noinv), nbnd ~1.5x the occupied spinor bands
output: orbm.json {scf: {...}, nscf: {rc, sec, orbm_lines, errors}}. Any failure is recorded, never raised: the SOC gate
then falls back to |M_spin| from deep_soc.json plus the ligand-field g-shift bound (kinetic.md 4.2).
2026-10-03 04:55 fix (job lcm/pba/soc/KVMo_soc_r): the first lorbm nscf (in lcm/pba/soc/KVMo_soc) stopped after 1 s with
"Error in routine iosys (1): Berry Phase/electric fields only for insulators!" because it inherited smeared occupations.
QE requires occupations = 'fixed' for lorbm; the nscf now sets occupations 'fixed' (degauss 0). The SCF is unchanged
(smeared, like lcm_deep soc; SOC-z gap 1.37 eV, so the smearing is immaterial)."""
import json, os, re, sys, shutil
import numpy as np
from ase import Atoms
from qeutil import write_pw, run_pw, parse_pw, save_json, SCRATCH, zvalence

NC = int(os.environ.get("NCORES", 32))
cfg = json.load(open(sys.argv[1]))
ECUT = cfg.get("ecut", 90.0); KSP = cfg.get("kspacing", 0.25); MOM = cfg.get("moment_mag", 3.0)
KORB = cfg.get("orbm_kspacing", 0.30)
at = Atoms(cfg["species"], cell=cfg["lattice"], scaled_positions=cfg["frac"], pbc=True)
U0 = {k: tuple(v) for k, v in cfg["hubbard"].items()}
mom = [(s * MOM, 0.0, 0.0) if s != 0 else (0.0, 0.0, 0.0) for s in cfg["spins"]]
out = {"id": cfg["id"]}
pref = "orbm"
try:
    write_pw(f"{pref}.in", at, pset="dojo_fr", ecutwfc=ECUT, kspacing=KSP, moments=mom, smearing="gaussian", degauss=0.005,
             hubbard=U0, prefix=pref, conv_thr=1e-8 * len(at), tstress=False, tprnfor=False, noncolin=True, lspinorb=True,
             electrons={"mixing_beta": 0.2, "electron_maxstep": 400})
    r = run_pw(f"{pref}.in", ncores=NC, npool=max(1, NC // 4), timeout=4 * 3600)
    p = parse_pw(r["out"]); p["sec"] = r["seconds"]
    out["scf"] = {k: p.get(k) for k in ("energy_eV", "converged_scf", "total_mag_vec", "abs_mag", "site_moments", "homo_lumo", "sec", "errors")}
    save_json(out, "orbm.json")
    nel = sum(zvalence(el, "dojo_fr") for el in cfg["species"])
    nbnd = int(np.ceil(1.5 * nel / 2.0)) * 2
    out["nbnd"] = nbnd
    write_pw(f"{pref}_nscf.in", at, calculation="nscf", pset="dojo_fr", ecutwfc=ECUT, kspacing=KORB, moments=mom,
             smearing="gaussian", degauss=0.005, hubbard=U0, prefix=pref, tstress=False, tprnfor=False, noncolin=True,
             lspinorb=True, nbnd=nbnd, control={"lorbm": True},
             system={"nosym": True, "noinv": True, "occupations": "fixed", "degauss": 0.0})
    r2 = run_pw(f"{pref}_nscf.in", ncores=NC, npool=max(1, NC // 4), timeout=4 * 3600)
    txt = open(r2["out"], errors="replace").read()
    lines = txt.splitlines()
    idx = [i for i, ln in enumerate(lines) if re.search(r"orbital magnet|M_LC|M_IC|M_tot|Kubo", ln, re.I)]
    keep = sorted({j for i in idx for j in range(i, min(len(lines), i + 8))})
    out["nscf"] = {"rc": r2["rc"], "sec": r2["seconds"], "job_done": "JOB DONE" in txt,
                   "orbm_lines": [lines[j] for j in keep][:80],
                   "errors": re.findall(r"Error in routine.*\n.*", txt)[:3]}
except Exception as e:  # best effort
    out["exception"] = repr(e)
save_json(out, "orbm.json")
shutil.rmtree(SCRATCH / pref, ignore_errors=True)
print("orbm", json.dumps(out.get("nscf", out.get("exception")))[:2000], flush=True)
