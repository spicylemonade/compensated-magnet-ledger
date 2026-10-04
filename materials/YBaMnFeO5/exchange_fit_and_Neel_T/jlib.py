"""Shared helpers for the direction-resolved exchange model of rock-salt YBaMnFeO5 (Track L)."""
import json
import numpy as np
from pymatgen.core import Structure, Lattice
MAG = ("Mn", "Fe")


def cls(si, sj, d, dz):
    pair = "".join(sorted((si, sj)))
    if pair == "FeMn":
        if dz < 0.6 and d < 4.3: return "ip"
        if dz > 3.8 and d < 4.3: return "ap"
        if 3.3 < dz < 3.8 and d < 3.9: return "Y"
        if d < 7.0: return f"FeMn_{round(d, 1)}"
        return None
    if d < 6.0:
        return f"{pair}_{'ip' if dz < 0.6 else 'oop'}_{round(d, 1)}"
    return None

def features(sc, spins, classes, rmax=7.0):
    f = np.zeros(len(classes))
    idx = [i for i, s in enumerate(sc) if s.specie.symbol in MAG]
    for i in idx:
        for nb in sc.get_neighbors(sc[i], rmax):
            j = nb.index
            if j not in idx: continue
            k = cls(sc[i].specie.symbol, sc[j].specie.symbol, nb.nn_distance, abs(nb.coords[2] - sc[i].coords[2]))
            if k in classes:
                f[classes.index(k)] += 0.5 * spins[i] * spins[j]
    return f

CLASSES = ["ip", "ap", "Y", "MnMn_ip_5.6", "FeFe_ip_5.6", "MnMn_oop_5.4", "FeFe_oop_5.4", "MnMn_oop_5.6", "FeFe_oop_5.6", "FeMn_6.7", "FeMn_6.9"]
