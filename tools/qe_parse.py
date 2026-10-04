"""Minimal, dependency-light parser for Quantum ESPRESSO pw.x (v7.x) text output.

Only numpy is needed. Reads plain or gzipped files. Written for this ledger so that every number can be
re-derived from the raw outputs without the campaign's own (cloud-bound) tooling.

Main entry points
    parse_pw(path)          -> dict(energy_eV, total_mag, abs_mag, site_moments, fermi_eV, converged, ...)
    last_bands(path)        -> dict(k, eig[nk, nspin, nbnd], occ[nk, nspin, nbnd] or None, fermi_eV)
    spin_windows(bands)     -> dict(gap, vbm_spin, cbm_spin, win_VB, win_CB, gap_up, gap_dn, unipolar)
    final_cell(path)        -> (cell[3x3] in Angstrom, species list, fractional positions) from a (vc-)relax

Conventions (identical to the original analysis code)
    spin index 0 = QE "SPIN UP", 1 = "SPIN DOWN"
    a state is occupied if its occupation is > 0.5 (equivalently, below the Fermi level for Gaussian smearing)
    win_VB = VBM(edge spin) - VBM(other spin)   > 0 : the top win_VB of the valence band holds only one spin
    win_CB = CBM(other spin) - CBM(edge spin)   > 0 : the bottom win_CB of the conduction band holds only one spin
    'unipolar' = valence and conduction band edges in the same spin channel
"""
from __future__ import annotations

import gzip
import pathlib
import re

import numpy as np

RY_EV = 13.605693122994
BOHR_A = 0.529177210903
_NUM = re.compile(r"-?\d+\.\d+")


def _lines(path):
    p = pathlib.Path(path)
    if not p.exists() and pathlib.Path(str(p) + ".gz").exists():
        p = pathlib.Path(str(p) + ".gz")
    op = gzip.open if p.suffix == ".gz" else open
    with op(p, "rt", errors="ignore") as fh:
        return fh.read().splitlines()


def resolve(path):
    """Return the existing path, trying a '.gz' suffix if the plain file is absent."""
    p = pathlib.Path(path)
    if p.exists():
        return p
    q = pathlib.Path(str(p) + ".gz")
    if q.exists():
        return q
    raise FileNotFoundError(path)


def parse_pw(path):
    L = _lines(path)
    out = {"path": str(path)}
    energies = [float(m.group(1)) for l in L if (m := re.match(r"^!+\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", l))]
    out["energy_eV"] = energies[-1] * RY_EV if energies else None
    out["final_marker"] = None
    for l in reversed(L):
        if re.match(r"^!+\s+total energy", l):
            out["final_marker"] = l.split()[0]
            break
    tm = [float(m.group(1)) for l in L if (m := re.search(r"total magnetization\s+=\s+(-?\d+\.\d+)", l))]
    am = [float(m.group(1)) for l in L if (m := re.search(r"absolute magnetization\s+=\s+(-?\d+\.\d+)", l))]
    out["total_mag"] = tm[-1] if tm else None
    out["abs_mag"] = am[-1] if am else None
    # last "Magnetic moment per site" block
    idx = [i for i, l in enumerate(L) if "Magnetic moment per site" in l]
    mom = []
    if idx:
        for l in L[idx[-1] + 1:]:
            m = re.search(r"atom\s+\d+\s+\(R=\s*[\d.]+\)\s+charge=\s*(-?[\d.]+)\s+magn=\s*(-?[\d.]+)", l)
            if not m:
                if mom:
                    break
                continue
            mom.append(float(m.group(2)))
    out["site_moments"] = mom
    ef = [float(m.group(1)) for l in L if (m := re.search(r"the Fermi energy is\s+(-?\d+\.\d+)\s+ev", l))]
    ef2 = [(float(m.group(1)), float(m.group(2))) for l in L
           if (m := re.search(r"the spin up/dw Fermi energies are\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+ev", l))]
    out["fermi_eV"] = ef[-1] if ef else (ef2[-1] if ef2 else None)
    out["converged"] = any("convergence has been achieved" in l for l in L) and not any("convergence NOT achieved" in l for l in L)
    out["job_done"] = any("JOB DONE" in l for l in L)
    nat = [int(m.group(1)) for l in L[:400] if (m := re.search(r"number of atoms/cell\s+=\s+(\d+)", l))]
    out["nat"] = nat[0] if nat else None
    v = [l for l in L[:5] if "Program PWSCF" in l]
    out["code"] = v[0].strip() if v else None
    dexx = [float(m.group(1)) for l in L if (m := re.search(r"est\. exchange err \(dexx\)\s+=\s+(-?[\d.Ee+-]+)", l))]
    out["dexx_last_Ry"] = dexx[-1] if dexx else None
    press = [float(m.group(1)) for l in L if (m := re.search(r"P=\s*(-?[\d.]+)", l)) and "total   stress" in l]
    out["pressure_kbar"] = press[-1] if press else None
    out["bfgs_converged"] = any("bfgs converged" in l for l in L) or any("End final coordinates" in l for l in L)
    return out


def _parse_block(L, start, stop):
    """Parse one eigenvalue printout (SPIN UP / SPIN DOWN sections) between line indices."""
    spins = {}
    cur = None
    for i in range(start, stop):
        l = L[i]
        if "SPIN UP" in l:
            cur = 0; spins[0] = []; continue
        if "SPIN DOWN" in l:
            cur = 1; spins[1] = []; continue
        if cur is None:
            continue
        if re.search(r"^\s+k =", l) and "bands (ev)" in l:
            kk = _NUM.findall(l.split("(")[0].split("=")[1])
            spins[cur].append({"k": [float(x) for x in kk[:3]], "e": [], "o": [], "mode": "e"})
            continue
        if "occupation numbers" in l and spins.get(cur):
            spins[cur][-1]["mode"] = "o"; continue
        if spins.get(cur):
            nums = _NUM.findall(l)
            if nums and not re.search(r"[a-zA-Z]", l):
                rec = spins[cur][-1]
                (rec["e"] if rec["mode"] == "e" else rec["o"]).extend(float(x) for x in nums)
    return spins


def last_bands(path):
    """Eigenvalues (and occupations) of the last complete band printout in a pw.x output."""
    L = _lines(path)
    nks = [int(m.group(1)) for l in L if (m := re.search(r"number of k points=\s*(\d+)", l))]
    nk_expected = nks[-1] if nks else None
    ends = [i for i, l in enumerate(L) if "End of self-consistent calculation" in l or "End of band structure calculation" in l]
    for e in reversed(ends):
        stop = len(L)
        fermi = None
        for j in range(e + 1, len(L)):
            m = re.search(r"the Fermi energy is\s+(-?\d+\.\d+)\s+ev", L[j])
            m2 = re.search(r"the spin up/dw Fermi energies are\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+ev", L[j])
            m3 = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", L[j])
            if m or m2 or m3:
                stop = j
                fermi = float(m.group(1)) if m else ((float(m2.group(1)), float(m2.group(2))) if m2 else
                                                     0.5 * (float(m3.group(1)) + float(m3.group(2))))
                break
            if "End of self-consistent calculation" in L[j] or "End of band structure calculation" in L[j]:
                stop = j
                break
        sp = _parse_block(L, e, stop)
        if 0 not in sp or 1 not in sp:
            continue
        nk = len(sp[0])
        if nk == 0 or len(sp[1]) != nk or (nk_expected and nk != nk_expected):
            continue
        nb = min(min(len(r["e"]) for r in sp[s]) for s in (0, 1))
        eig = np.array([[sp[s][ik]["e"][:nb] for s in (0, 1)] for ik in range(nk)])
        have_occ = all(len(r["o"]) >= nb for s in (0, 1) for r in sp[s])
        occ = np.array([[sp[s][ik]["o"][:nb] for s in (0, 1)] for ik in range(nk)]) if have_occ else None
        return {"k": np.array([r["k"] for r in sp[0]]), "eig": eig, "occ": occ, "fermi_eV": fermi, "line": e}
    raise ValueError(f"no complete spin-polarised eigenvalue block found in {path}")


def spin_windows(b):
    E = b["eig"]
    if b["occ"] is not None:
        occd = b["occ"] > 0.5
    else:
        ef = b["fermi_eV"]
        if isinstance(ef, tuple):
            occd = np.stack([E[:, 0, :] < ef[0], E[:, 1, :] < ef[1]], axis=1)
        else:
            occd = E < ef
    vb = np.where(occd, E, -np.inf)
    cb = np.where(~occd, E, np.inf)
    tv = vb.max(axis=(0, 2))
    bc = cb.min(axis=(0, 2))
    sv = int(np.argmax(tv)); sc = int(np.argmin(bc))
    return {"gap": float(bc.min() - tv.max()), "vbm_spin": sv, "cbm_spin": sc,
            "win_VB": float(tv[sv] - tv[1 - sv]), "win_CB": float(bc[1 - sc] - bc[sc]),
            "gap_up": float(bc[0] - tv[0]), "gap_dn": float(bc[1] - tv[1]), "unipolar": sv == sc,
            "n_occ_per_spin": [int(occd[0, s].sum()) for s in (0, 1)], "nk": int(E.shape[0])}


def final_cell(path):
    """Final cell (Angstrom), species and fractional coordinates of a relax / vc-relax output."""
    L = _lines(path)
    alat = None
    for l in L[:400]:
        m = re.search(r"lattice parameter \(alat\)\s+=\s+([\d.]+)\s+a\.u\.", l)
        if m:
            alat = float(m.group(1)) * BOHR_A
            break
    cell = None
    species, frac = [], []
    for i, l in enumerate(L):
        if l.startswith("CELL_PARAMETERS"):
            unit = l
            rows = [[float(x) for x in _NUM.findall(L[i + j])] for j in (1, 2, 3)]
            c = np.array(rows)
            if "alat" in unit:
                m = re.search(r"alat=\s*([\d.]+)", unit)
                c = c * (float(m.group(1)) * BOHR_A if m else alat)
            elif "bohr" in unit:
                c = c * BOHR_A
            cell = c
        if l.startswith("ATOMIC_POSITIONS"):
            species, frac = [], []
            for l2 in L[i + 1:]:
                parts = l2.split()
                if len(parts) < 4 or not re.match(r"^[A-Z][a-z]?\d*$", parts[0]):
                    break
                species.append(re.sub(r"\d+$", "", parts[0]))
                frac.append([float(x) for x in parts[1:4]])
    return cell, species, np.array(frac)
