"""Launch helper SCF jobs for a big (72-atom) task once its relax is done in lcm/ord2/v4/<ID>:
   v4/<ID>_h1 runs --configs=FM,G and v4/<ID>_h2 runs --configs=r3,r2 at the same relaxed geometry (staged results JSON).
usage: python helpers.py ID [ID ...]"""
import json, os, shutil, subprocess, sys
import modal
v = modal.Volume.from_name("magdisc-data")
ROOT = "<campaign_root>"
for tid in sys.argv[1:]:
    d = json.loads(b"".join(v.read_file(f"jobs/lcm/ord2/v4/{tid}/results/{tid}.json")))
    if not d.get("relaxed"):
        print(tid, "not relaxed yet"); continue
    tag = tid.split("_")[0].lower()
    for h, cf in (("h1", "FM,G"), ("h2", "r3,r2")):
        st = f"{ROOT}/infra/jobs_src/lcm_ord2_helpers/{tid}_{h}"
        os.makedirs(st + "/results", exist_ok=True)
        shutil.copy(f"{ROOT}/infra/jobs_src/lcm_ord2/ord2.py", st)
        shutil.copy(f"{ROOT}/infra/jobs_src/lcm_ord2/tasks_{tag}_v4.json", st)
        json.dump({"id": tid, "relaxed": d["relaxed"], "relax": d.get("relax"), "runs": {}}, open(f"{st}/results/{tid}.json", "w"))
        cmd = ["python", f"{ROOT}/infra/mrun.py", "submit", "--fn", "run_cpu32", "--job", f"lcm/ord2/v4/{tid}_{h}", "--dir", st,
               "--cmd", f"python ord2.py tasks_{tag}_v4.json {tid} --configs={cf}"]
        print(subprocess.run(cmd, capture_output=True, text=True).stdout.strip()[-120:])
