"""Charge-transfer probe: re-run selected configurations of finished arrangements WITHOUT the fixed total moment (FIXM=0),
from the d5 start (Hubbard_occ Fe = 5) and from QE's default start (Fe d6 = Fe2+-like), at the relaxed geometry of the
constrained run. Builds a staging dir with ord2.py, tasks_ctp.json and results/<CTP id>.json (relaxed geometry
pre-filled, so the relax is skipped).
usage: python ctprobe.py STAGEDIR TAG ID1:cfg1,cfg2[:starts] [ID2:...]   (starts default d5+qe)"""
import json, os, shutil, sys
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "../../.."))
stage, tag = sys.argv[1], sys.argv[2]
os.makedirs(stage + "/results", exist_ok=True)
shutil.copy(f"{ROOT}/infra/jobs_src/lcm_ord2/ord2.py", stage)
tasks = {t["id"]: t for t in json.load(open(f"{ROOT}/infra/jobs_src/lcm_ord2/tasks_{tag.lower()}.json"))["tasks"]}
out = []
for spec in sys.argv[3:]:
    parts = spec.split(":")
    aid, cfgs = parts[0], parts[1].split(",")
    starts = parts[2].split("+") if len(parts) > 2 else ["d5", "qe"]
    d = json.load(open(f"{HERE}/results_{tag}/{aid}.json"))
    t = tasks[aid]
    for st in starts:
        nid = f"CTP_{aid}_{st}"
        nt = dict(t); nt["id"] = nid; nt["start"] = st
        nt["configs"] = [c for c in t["configs"] if c["name"] in cfgs]
        out.append(nt)
        json.dump({"id": nid, "relaxed": d["relaxed"], "relax": {"from": aid}, "runs": {}}, open(f"{stage}/results/{nid}.json", "w"))
json.dump({"tasks": out}, open(f"{stage}/tasks_ctp.json", "w"))
print(" ".join(t["id"] for t in out))
