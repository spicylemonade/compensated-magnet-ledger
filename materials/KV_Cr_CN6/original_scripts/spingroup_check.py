"""F5 spin-group classification of a collinear order with three independent tools (amcheck, findspingroup, amscreen),
reusing the shared wrappers in data/magndata/magndata_lib.py (read-only). Writes a P1 MAGNDATA-style mcif.
usage: python tracks/lcm/spingroup_check.py cand.json CONFIGNAME [CONFIG2 ...]"""
import json, sys, os, pathlib, warnings
import numpy as np
warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "data" / "magndata"))
import magndata_lib as ml
from pymatgen.core import Structure, Lattice
TMP = "/private/tmp/claude-501/-Users-geby-chemistry/0128b63c-3ab4-49af-9683-b3bd1a9fe21d/scratchpad"


def write_p1_mcif(path, S, spins, mom=4.0):
    a, b, c = S.lattice.abc; al, be, ga = S.lattice.angles
    L = ["#\\#CIF_2.0", "data_lcmcheck", "_parent_space_group.name_H-M_alt  'P 1'", "_parent_space_group.IT_number 1",
         "loop_", "_parent_propagation_vector.id", "_parent_propagation_vector.kxkykz", "k1 [0 0 0]",
         '_space_group_magn.number_BNS  1.1', '_space_group_magn.name_BNS  "P 1"',
         f"_cell_length_a {a:.6f}", f"_cell_length_b {b:.6f}", f"_cell_length_c {c:.6f}",
         f"_cell_angle_alpha {al:.4f}", f"_cell_angle_beta {be:.4f}", f"_cell_angle_gamma {ga:.4f}",
         "loop_", "_space_group_symop_magn_operation.id", "_space_group_symop_magn_operation.xyz", "1 x,y,z,+1",
         "loop_", "_space_group_symop_magn_centering.id", "_space_group_symop_magn_centering.xyz", "1 x,y,z,+1",
         "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z", "_atom_site_occupancy"]
    labs = []
    for i, site in enumerate(S):
        lab = f"{site.specie.symbol}{i+1}"; labs.append(lab)
        x, y, z = site.frac_coords % 1.0
        L.append(f"{lab} {site.specie.symbol} {x:.6f} {y:.6f} {z:.6f} 1")
    L += ["loop_", "_atom_site_moment.label", "_atom_site_moment.crystalaxis_x", "_atom_site_moment.crystalaxis_y", "_atom_site_moment.crystalaxis_z"]
    # moment along c (crystal axis z) ; magnitude mom
    for lab, s in zip(labs, spins):
        if s != 0:
            L.append(f"{lab} 0 0 {s * mom:.3f}")
    open(path, "w").write("\n".join(L) + "\n")


def check(c, spins, tag, mom=4.0):
    S = Structure(Lattice(c["lattice"]), c["species"], c["frac"])
    p = os.path.join(TMP, f"lcmcheck_{tag}.mcif")
    write_p1_mcif(p, S, spins, mom)
    L = S.lattice.matrix; frac = S.frac_coords % 1.0
    species = [s.specie.symbol for s in S]; numbers = [s.specie.Z for s in S]
    out = {}
    try:
        frac_s, _ = ml.symmetrize_positions(L, frac, numbers)
    except Exception:
        frac_s = frac
    try:
        out["amcheck"] = ml.run_amcheck(L, frac_s, species, numbers, spins)
    except Exception as e:
        out["amcheck"] = f"error:{e}"
    try:
        v, w, info = ml.run_amscreen(L, frac_s, species, spins)
        out["amscreen"] = v + (f"({w})" if w else "")
    except Exception as e:
        out["amscreen"] = f"error:{e}"
    try:
        r = ml.run_findspingroup(p)
        out["findspingroup"] = r["verdict"] + (f"({r['wave']})" if r.get("wave") else "")
        out["fsg_phase"] = r["phase"]; out["fsg_ssg"] = r["ssg_index"]; out["fsg_splitting_wo_soc"] = r["spinsplitting_wo_soc"]
    except Exception as e:
        out["findspingroup"] = f"error:{type(e).__name__}:{str(e)[:80]}"
    return out


if __name__ == "__main__":
    c = json.load(open(sys.argv[1]))
    c = c[0] if isinstance(c, list) else c
    for name in sys.argv[2:]:
        cf = [x for x in c["configs"] if x["name"] == name][0]
        print(c["id"], name, json.dumps(check(c, cf["spins"], tag=f"{c['id']}_{name}")))
