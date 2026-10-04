"""Independent re-run of selected ledger calculations in a clean cloud container (Modal).

What it does, per case:
  1. builds a fresh container image with Quantum ESPRESSO 7.5 from conda-forge (no campaign code, no campaign files);
  2. downloads the pseudopotentials from their public source (on your machine), ships them to the container
     and checks their md5 sums there against reproduce/pseudo_md5.txt;
  3. runs pw.x on the input files stored in this repository (only pseudo_dir/outdir are rewritten);
  4. writes the outputs to reproduce/results/<case>/ (and keeps a copy on the Modal volume `lcm-ledger-repro`).

Usage (needs a Modal account; `pip install modal`):
    modal run reproduce/modal_repro.py --case kvcr_pbeu
    modal run reproduce/modal_repro.py --case all
Then check the numbers with:
    python tools/verify.py --repro

Without Modal, the same thing by hand: reproduce/get_pseudos.sh + reproduce/run_qe.sh (see reproduce/README.md).
"""
import json
import pathlib
import time

import modal

REPO = pathlib.Path(__file__).resolve().parent.parent
K = "materials/KV_Cr_CN6/runs/"
Y = "materials/YBaMnFeO5/runs/"

# Each case: list of input files (run in order) + pseudopotential family + cores.
CASES = {
    # KV[Cr(CN)6], PBE+U (U_V = U_Cr = 3 eV): LCM scf, LCM nscf (8x8x8, for the spin windows), FM scf (exchange energy)
    "kvcr_pbeu": {"inputs": [K + "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.in",
                             K + "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.in",
                             K + "A_pbeu_relax_scf_nscf_Ugrid/KVCr_FM.in"],
                  "pseudos": "dojo", "cpu": 32, "nk": 8},
    # KV[Cr(CN)6], HSE06 (no U), ideal anhydrous cell
    "kvcr_hse": {"inputs": [K + "C_hse06_ideal_LCM/hse_lcm.in"], "pseudos": "dojo", "cpu": 64, "nk": 8},
    # YBaMnFeO5, PBE+U (U_Mn = U_Fe = 4 eV): G-type LCM vs the closest competitor (Y-layer spin flip)
    "ybmfo_stack_U44": {"inputs": [Y + "H_pbeu_stacking_vs_U/U44_G/single.in",
                                   Y + "H_pbeu_stacking_vs_U/U44_Yflip/single.in"],
                        "pseudos": "dojo", "cpu": 32, "nk": 16},
    # Completion, not reproduction: the KV[Cr(CN)6]*2H2O HSE06 run whose exact-exchange loop was stopped before
    # convergence in the original campaign. Identical input, run until QE's own EXX criterion is met.
    "kvcr_hydrate_hse": {"inputs": [K + "F_hse06_dihydrate_LCM/hse_lcm.in"], "pseudos": "dojo", "cpu": 64, "nk": 8,
                         "timeout": 20 * 3600},
    # Robustness, not reproduction: KV[Cr(CN)6] PBE+U with a different pseudopotential family (SSSP 1.3 efficiency)
    "kvcr_pbeu_sssp": {"inputs": [K + "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.in",
                                  K + "A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.in"],
                       "pseudos": "sssp", "cpu": 32, "nk": 8},
}

DOJO_URL = "https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz"
SSSP_TGZ = "https://archive.materialscloud.org/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_efficiency.tar.gz?download=1"
SSSP_JSON = "https://archive.materialscloud.org/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_efficiency.json?download=1"

app = modal.App("lcm-ledger-repro")
vol = modal.Volume.from_name("lcm-ledger-repro", create_if_missing=True)
image = (
    modal.Image.micromamba(python_version="3.11")
    .apt_install("wget", "ca-certificates")
    .micromamba_install("qe=7.5", "openmpi", channels=["conda-forge"])
    .env({"OMP_NUM_THREADS": "1", "OMPI_ALLOW_RUN_AS_ROOT": "1", "OMPI_ALLOW_RUN_AS_ROOT_CONFIRM": "1",
          "OMPI_MCA_rmaps_base_oversubscribe": "1", "PMIX_MCA_gds": "hash"})
)


def _rewrite(text, pseudo_dir, outdir, sssp=None):
    """Change pseudo_dir/outdir; for the SSSP robustness case also swap pseudopotential files and cutoffs."""
    import re
    out = []
    in_species = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("pseudo_dir"):
            line = f"  pseudo_dir = '{pseudo_dir}'"
        elif s.startswith("outdir"):
            line = f"  outdir = '{outdir}'"
        elif sssp and s.startswith("ecutwfc"):
            line = f"  ecutwfc = {sssp['ecutwfc']:.1f}"
        elif sssp and s.startswith("ecutrho"):
            line = f"  ecutrho = {sssp['ecutrho']:.1f}"
        if s.startswith("ATOMIC_SPECIES"):
            in_species = True
        elif in_species and (s == "" or re.match(r"^[A-Z_]+(\s|$)", s) and not re.match(r"^[A-Z][a-z]?\d*\s+\d", s)):
            in_species = False
        elif in_species and sssp:
            lab, mass, _ = s.split()
            el = re.sub(r"\d+$", "", lab)
            line = f"  {lab} {mass} {sssp['files'][el]}"
        out.append(line)
    return "\n".join(out) + "\n"


@app.function(image=image, cpu=32, memory=65536, volumes={"/out": vol}, timeout=8 * 3600)
def run_case(case: str, files: dict, md5_table: str, cpu: int, nk: int, pseudos: str, archives: dict) -> dict:
    import hashlib, json, os, re, subprocess
    t_start = time.time()
    work = pathlib.Path("/tmp/work"); work.mkdir(parents=True, exist_ok=True)
    pdir = pathlib.Path("/tmp/pseudo"); pdir.mkdir(parents=True, exist_ok=True)
    log = []
    sssp = None
    if pseudos == "dojo":
        (pdir / "d.tgz").write_bytes(archives["dojo"])
        subprocess.run(f"cd {pdir} && tar xzf d.tgz && rm d.tgz", shell=True, check=True)
        for f in list(pdir.rglob("*.upf")):
            if f.parent != pdir:
                f.rename(pdir / f.name)
        for row in md5_table.splitlines():
            if not row.strip() or row.startswith("#"):
                continue
            el, fn, want = row.split()[:3]
            p = pdir / fn
            got = hashlib.md5(p.read_bytes()).hexdigest() if p.exists() else "MISSING"
            log.append(f"{fn:8s} {'ok' if got == want else 'DIFFERS'} {got}")
    else:  # SSSP 1.3 efficiency (robustness check with a different pseudopotential family)
        (pdir / "s.tgz").write_bytes(archives["sssp_tgz"])
        subprocess.run(f"cd {pdir} && tar xzf s.tgz && rm s.tgz", shell=True, check=True)
        for f in list(pdir.rglob("*")):
            if f.is_file() and f.parent != pdir:
                f.rename(pdir / f.name)
        meta = json.loads(archives["sssp_json"].decode())
        els = set()
        for txt in files.values():
            m = re.search(r"ATOMIC_SPECIES\n(.*?)\n(?:CELL_PARAMETERS|ATOMIC_POSITIONS)", txt, re.S)
            for l in m.group(1).splitlines():
                els.add(re.sub(r"\d+$", "", l.split()[0]))
        ecut = max(meta[e]["cutoff_wfc"] for e in els)
        sssp = {"files": {e: meta[e]["filename"] for e in els}, "ecutwfc": float(ecut), "ecutrho": float(8 * ecut)}
        for e in sorted(els):
            fp = pdir / meta[e]["filename"]
            got = hashlib.md5(fp.read_bytes()).hexdigest() if fp.exists() else "MISSING"
            log.append(f"{meta[e]['filename']} {'ok' if got == meta[e]['md5'] else 'DIFFERS'} {got}")
        log.append(f"SSSP efficiency 1.3: {sssp}")
    results = {"case": case, "md5_check": log, "runs": []}
    vdir = pathlib.Path("/out") / case
    vdir.mkdir(parents=True, exist_ok=True)
    for rel, txt in files.items():
        d = work / pathlib.Path(rel).parent
        d.mkdir(parents=True, exist_ok=True)
        prefix = re.search(r"prefix\s*=\s*'([^']+)'", txt).group(1)
        outdir = f"/tmp/scratch/{abs(hash(str(pathlib.Path(rel).parent))) % 10**8}/{prefix}"
        pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)  # pw.x creates only the last directory level
        inp = d / (pathlib.Path(rel).stem + ".repro.in")
        inp.write_text(_rewrite(txt, str(pdir), outdir, sssp))
        out = inp.with_suffix(".out")
        # MPI only: Modal sets OMP_NUM_THREADS to the core count, which would start cpu x cpu threads.
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        cmd = (f"mpirun --allow-run-as-root -np {cpu} -x OMP_NUM_THREADS --bind-to none --mca pml ob1 --mca btl self,sm "
               f"--mca osc ^ucx pw.x -nk {nk} -in {inp} > {out} 2>&1")
        t0 = time.time()
        rc = subprocess.run(cmd, shell=True, cwd=d, env=env).returncode
        rec = {"input": rel, "rc": rc, "seconds": round(time.time() - t0, 1),
               "repro_in": inp.read_text(), "repro_out": out.read_text(errors="ignore")}
        results["runs"].append(rec)
        tgt = vdir / pathlib.Path(rel).parent.name
        tgt.mkdir(parents=True, exist_ok=True)
        (tgt / inp.name).write_text(rec["repro_in"])
        (tgt / out.name).write_text(rec["repro_out"])
        vol.commit()
    results["total_seconds"] = round(time.time() - t_start, 1)
    m = re.search(r"Program PWSCF (v\.[^ ]+)", results["runs"][0]["repro_out"]) if results["runs"] else None
    results["pw_version"] = m.group(1) if m else "unknown"
    (vdir / "summary.json").write_text(json.dumps({k: v for k, v in results.items() if k != "runs"} |
                                                  {"runs": [{k: r[k] for k in ("input", "rc", "seconds")} for r in results["runs"]]}, indent=1))
    vol.commit()
    return results


def _download(url, tries=3):
    """Fetch a public archive on the local machine. Integrity is checked by md5 inside the container, so a
    server with an incomplete TLS certificate chain (pseudo-dojo.org at the time of writing) is tolerated."""
    import shutil, ssl, subprocess, tempfile, urllib.request
    for i in range(tries):
        if shutil.which("curl"):
            with tempfile.NamedTemporaryFile() as tmp:
                if subprocess.run(["curl", "-sSL", "--max-time", "300", url, "-o", tmp.name]).returncode == 0:
                    data = pathlib.Path(tmp.name).read_bytes()
                    if data:
                        return data
        for ctx in (None, ssl._create_unverified_context()):
            try:
                with urllib.request.urlopen(url, timeout=300, context=ctx) as r:
                    return r.read()
            except Exception as e:
                print(f"download failed ({e}); retry {i + 1}/{tries}")
        time.sleep(5 * (i + 1))
    raise RuntimeError(f"could not download {url}")


@app.local_entrypoint()
def main(case: str = "kvcr_pbeu"):
    names = [n for n in CASES if n != "kvcr_hydrate_hse"] if case == "all" else case.split(",")
    md5_table = (REPO / "reproduce" / "pseudo_md5.txt").read_text()
    archives = {}
    if any(CASES[n]["pseudos"] == "dojo" for n in names):
        archives["dojo"] = _download(DOJO_URL)
    if any(CASES[n]["pseudos"] == "sssp" for n in names):
        archives["sssp_tgz"] = _download(SSSP_TGZ)
        archives["sssp_json"] = _download(SSSP_JSON)
    calls = []
    for name in names:
        c = CASES[name]
        files = {rel: (REPO / rel).read_text() for rel in c["inputs"]}
        fn = run_case.with_options(cpu=c["cpu"], memory=2048 * c["cpu"], timeout=c.get("timeout", 8 * 3600))
        need = {"dojo": ["dojo"], "sssp": ["sssp_tgz", "sssp_json"]}[c["pseudos"]]
        calls.append((name, fn.spawn(name, files, md5_table, c["cpu"], c["nk"], c["pseudos"], {k: archives[k] for k in need})))
    for name, call in calls:
        try:
            res = call.get()
        except Exception as e:
            print(f"[{name}] FAILED: {e!r}")
            continue
        dest = REPO / "reproduce" / "results" / name
        dest.mkdir(parents=True, exist_ok=True)
        for r in res["runs"]:
            sub = dest / pathlib.Path(r["input"]).parent.name
            sub.mkdir(parents=True, exist_ok=True)
            stem = pathlib.Path(r["input"]).stem
            (sub / f"{stem}.repro.in").write_text(r["repro_in"])
            (sub / f"{stem}.repro.out").write_text(r["repro_out"])
        summary = {k: v for k, v in res.items() if k != "runs"}
        summary["runs"] = [{k: r[k] for k in ("input", "rc", "seconds")} for r in res["runs"]]
        (dest / "summary.json").write_text(json.dumps(summary, indent=1))
        print(f"[{name}] done in {res['total_seconds']} s; md5: {res['md5_check'][:3]}...")
