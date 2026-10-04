"""Track L PBA Gate 0: stage-1 chain for one cyanide double perovskite, run inside ONE Modal container.
usage: python pba_chain.py chain.json          (env ECUT/KSCF/KNSCF/NPOOL/NPOOL_NSCF/SCF_MAXSEC/RELAX_MAXSEC as lcm_stage1.py)
chain.json: {"id", "A", "MN", "MC", "base": <lcm_stage1 candidate with relax_config "LCM">,
             "grid": [{"tag": "U2_4", "hubbard": {el: [orb, U]}}], "grid_nscf": bool}
steps (all resumable; lcm_stage1.py skips finished candidates / configs and continues an interrupted vc-relax):
  1. lcm_stage1.py [base]  -> results/<id>.json  (PBE+U vc-relax in LCM, SCF LCM + FM, dense nscf LCM)
  2. exact F-43m / Fm-3m projection of the relaxed cell (pba_geom.symmetrize) -> results/<id>_sym.json
     (if the relax failed, the input geometry is used and flagged "relax_failed")
  3. lcm_stage1.py [grid points] at the symmetrized geometry -> results/<id>_<tag>.json
     (LCM SCF + dense nscf (if grid_nscf) + FM SCF auto-added by lcm_stage1)
lcm_stage1.py is the unchanged pivot/garnet copy (md5 be16bdd5...)."""
import json, os, subprocess, sys, time
import numpy as np
from pba_geom import symmetrize, pba_cell

ch = json.load(open(sys.argv[1]))
cid = ch["id"]
os.makedirs("results", exist_ok=True)


def stage1(cands, tag):
    path = f"chain_{tag}.json"
    json.dump(cands, open(path, "w"), indent=1)
    t0 = time.time()
    r = subprocess.run([sys.executable, "lcm_stage1.py", path], env=os.environ.copy())
    print(f"[chain {cid}] stage1 {tag}: rc={r.returncode} {time.time() - t0:.0f}s", flush=True)
    return r.returncode


# 1. base: relax + LCM/FM + nscf
stage1([ch["base"]], "base")
res = json.load(open(f"results/{cid}.json")) if os.path.exists(f"results/{cid}.json") else {}

# 2. symmetrize the relaxed geometry
b = ch["base"]
sym = {"id": cid, "t": time.time()}
if res.get("relaxed"):
    rl = res["relaxed"]
    cell, info = symmetrize(rl["lattice"], rl["frac"], rl["species"], ch["MN"], ch["MC"], ch.get("A"))
    sym.update(cell=cell, info=info, source="relaxed", relax_parse=res.get("relax_parse"))
else:
    cell, info = symmetrize(b["lattice"], b["frac"], b["species"], ch["MN"], ch["MC"], ch.get("A"))
    sym.update(cell=cell, info=info, source="input", relax_failed=True, relax_err=res.get("relax_err"))
json.dump(sym, open(f"results/{cid}_sym.json", "w"), indent=1)
print(f"[chain {cid}] symmetrized ({sym['source']}): a={info['a']:.4f} xN={info['xN']:.5f} xC={info['xC']:.5f} "
      f"bonds={ {k: round(v, 3) for k, v in info['bonds'].items()} } max_dev={info['max_atom_dev_A']:.2e} A", flush=True)

# 3. U grid at the symmetrized geometry
grid = []
for g in ch.get("grid", []):
    c = {k: v for k, v in b.items() if k not in ("relax_config",)}
    c.update(id=f"{cid}_{g['tag']}", lattice=cell["lattice"], species=cell["species"], frac=cell["frac"],
             hubbard=g["hubbard"], nscf_configs=(["LCM"] if ch.get("grid_nscf", True) else []))
    grid.append(c)
if grid:
    stage1(grid, "grid")
print(f"[chain {cid}] done", flush=True)
