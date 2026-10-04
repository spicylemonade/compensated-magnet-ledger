"""Compact progress table for lcm_ord2 v2 jobs: relax BFGS steps (number of '!' energies), last force, finished SCF configs."""
import modal, json, re, sys, time
v = modal.Volume.from_name("magdisc-data")
ROOT = "jobs/lcm/ord2/v4"
pref = sys.argv[1] if len(sys.argv) > 1 else ""
rows = []
for e in sorted(v.listdir(ROOT), key=lambda e: e.path):
    j = e.path.split("/")[-1]
    if not j.startswith(pref): continue
    files = {x.path.split("/")[-1]: x for x in v.listdir(e.path)}
    st = {}
    try: st = json.loads(b"".join(v.read_file(e.path + "/status.json")))
    except Exception: pass
    rel = [f for f in files if f.endswith("_relax.out")]
    nst, F, done = 0, "-", False
    if rel:
        t = b"".join(v.read_file(e.path + "/" + rel[0])).decode(errors="replace")
        nst = len(re.findall(r"^!", t, re.M)); F = (re.findall(r"Total force =\s+([0-9.]+)", t) or ["-"])[-1]; done = "End of BFGS" in t
    runs = []
    try:
        rr = [x.path for x in v.listdir(e.path + "/results")]
        if rr:
            d = json.loads(b"".join(v.read_file(rr[0])))
            runs = [k for k, r in d.get("runs", {}).items() if r.get("energy_eV") is not None]
    except Exception: pass
    scfs = [f for f in files if f.endswith(".out") and "_relax" not in f]
    el = (time.time() - st.get("start", time.time())) / 3600
    print("%-18s %-7s %5.2fh relax %2d steps F=%-9s%s scf_done=%s cur=%s" % (j, st.get("state", "?")[:7], el, nst, F, " CONV" if done else "", ",".join(runs), ",".join(f.split("_")[-1][:-4] for f in scfs if f.split("_")[-1][:-4] not in runs)))
