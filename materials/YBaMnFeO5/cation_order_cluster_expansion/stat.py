"""quick progress of lcm_ord2 jobs: state, elapsed, relax BFGS steps so far, finished SCF runs"""
import modal, json, re, sys, time
v = modal.Volume.from_name("magdisc-data")
ROOT = "jobs/lcm/ord2/v4"
pref = sys.argv[1] if len(sys.argv) > 1 else ""
for e in sorted(v.listdir(ROOT), key=lambda e: e.path):
    j = e.path
    if pref and not j.split("/")[-1].startswith(pref): continue
    try: st = json.loads(b"".join(v.read_file(j + "/status.json")))
    except Exception: st = {}
    files = {x.path.split("/")[-1]: x.size for x in v.listdir(j)}
    nst = ""; 
    rel = [f for f in files if f.endswith("_relax.out")]
    if rel:
        t = b"".join(v.read_file(j + "/" + rel[0])).decode(errors="replace")
        nst = "relax: %d scf, F=%s E=%s" % (len(re.findall(r"^!", t, re.M)), (re.findall(r"Total force =\s+([0-9.]+)", t) or ["-"])[-1], (re.findall(r"^!\s+total energy\s+=\s+([-0-9.]+)", t, re.M) or ["-"])[-1])
        if "End of BFGS" in t: nst += " DONE"
    scf = [f for f in files if f.endswith(".out") and "_relax" not in f]
    el = (time.time() - st.get("start", time.time())) / 3600
    print("%-28s %-8s %5.2fh %s | scf outs: %s" % (j.split("/")[-1], st.get("state"), el if st.get("state") == "running" else st.get("hours", 0) or 0, nst, ",".join(sorted(f.split("_")[-1][:-4] for f in scf))))
