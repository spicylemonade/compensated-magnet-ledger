"""Fetch lcm_ord2 results from the Modal volume into tracks/lcm/lcm_ord2/results_<TAG>/ (merging helper-job runs).
usage: python fetch.py TAG [TAG...]     (run from the project root or anywhere)
Job dirs are taken from infra/job_registry.jsonl (lcm/ord2/v4/<TAG>_*), and results are read from the known path
jobs/<job>/results/<ID>.json (no VolumeListFiles calls: listdir is rate-limited)."""
import json, os, sys, pathlib
from concurrent.futures import ThreadPoolExecutor
import modal

HERE = pathlib.Path(__file__).resolve().parent
REG = HERE.parents[2] / "infra" / "job_registry.jsonl"
v = modal.Volume.from_name("magdisc-data")
ROOT = "lcm/ord2/v4"


def rd(p):
    return json.loads(b"".join(v.read_file(p)))


def task_id(job):
    b = job.split("/")[-1]
    return b[:-3] if b.endswith(("_h1", "_h2", "_h3", "_h4")) else b


allj = []
for l in open(REG):
    try:
        j = json.loads(l)["job"]
    except Exception:
        continue
    if j.startswith(ROOT + "/") and j not in allj:
        allj.append(j)

for tag in sys.argv[1:]:
    out = HERE / f"results_{tag}"
    out.mkdir(exist_ok=True)
    jobs = sorted(j for j in allj if j.split("/")[-1].startswith(tag + "_"))

    def get(j):
        tid = task_id(j)
        try:
            return j, rd(f"jobs/{j}/results/{tid}.json")
        except Exception:
            return j, None

    with ThreadPoolExecutor(8) as ex:
        got = list(ex.map(get, jobs))
    merged = {}
    for j, d in got:
        if d is None:
            continue
        tid = d["id"]
        if tid not in merged:
            merged[tid] = d
        else:
            m = merged[tid]
            for k, r in d.get("runs", {}).items():
                if r.get("energy_eV") is not None and (k not in m["runs"] or m["runs"][k].get("energy_eV") is None):
                    m["runs"][k] = r
                elif k not in m["runs"]:
                    m["runs"][k] = r
            if not m.get("relaxed") and d.get("relaxed"):
                m["relaxed"] = d["relaxed"]; m["relax"] = d.get("relax")
        merged[tid].setdefault("jobs", []).append(j)
    for tid, d in merged.items():
        json.dump(d, open(out / f"{tid}.json", "w"))
    print(tag, len(merged), {t: len([r for r in d.get("runs", {}).values() if r.get("energy_eV") is not None]) for t, d in sorted(merged.items())})
