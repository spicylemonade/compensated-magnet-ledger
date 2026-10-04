"""SCF-convergence survey of the current relax/scf outputs of lcm_ord2 jobs (prefix filter)."""
import modal, re, sys
v = modal.Volume.from_name("magdisc-data")
ROOT = "jobs/lcm/ord2/v4"
pref = sys.argv[1] if len(sys.argv) > 1 else ""
for e in sorted(v.listdir(ROOT), key=lambda e: e.path):
    j = e.path.split("/")[-1]
    if not j.startswith(pref): continue
    outs = sorted([x.path for x in v.listdir(e.path) if x.path.endswith(".out")], key=lambda p: ("relax" not in p, p))
    if not outs: print(j, "no out"); continue
    msg = []
    for o in outs:
        t = b"".join(v.read_file(o)).decode(errors="replace")
        acc = re.findall(r"estimated scf accuracy\s+<\s+([0-9.Ee+-]+)", t)
        tm = re.findall(r"total cpu time spent up to now is\s+([0-9.]+)", t)
        nn = len(re.findall(r"^!", t, re.M))
        msg.append("%s: it %d !%d acc %s t %s%s" % (o.split("_")[-1][:-4], len(acc), nn, acc[-1] if acc else "-", tm[-1] if tm else "-", " BFGSdone" if "End of BFGS" in t else ""))
    print("%-16s %s" % (j, " | ".join(msg)))
