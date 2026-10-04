"""Seed geometries for the v4 restart: for every task take the best ionic positions reached in lcm/ord2/v3/<id>
(converged relax in results/<id>.json if clean, else the last positions of the relax output with the smallest final
total force). Writes infra/jobs_src/lcm_ord2/tasks_<tag>_v4.json (frac replaced, 'seed' records the source)."""
import json, re, sys
import modal
v = modal.Volume.from_name("magdisc-data")
ROOT = "<campaign_root>/infra/jobs_src/lcm_ord2"


def rd(p):
    return b"".join(v.read_file(p))


for tag in ("ybmfo", "ybcfo"):
    T = json.load(open(f"{ROOT}/tasks_{tag}.json"))
    out = []
    for t in T["tasks"]:
        tid = t["id"]; base = f"jobs/lcm/ord2/v3/{tid}"
        try:
            files = [e.path for e in v.listdir(base)]
        except Exception:
            out.append(t); print(tid, "no dir"); continue
        nat = len(t["species"])
        best = None
        try:
            d = json.loads(rd(f"{base}/results/{tid}.json"))
            r = d.get("relax", {})
            if d.get("relaxed") and r.get("bfgs_converged") and not r.get("errors"):
                best = (0.0, d["relaxed"]["frac"], "results-converged")
        except Exception:
            pass
        if best is None:
            for f in files:
                if "_relax.out" not in f:
                    continue
                txt = rd(f).decode(errors="replace")
                blocks = list(re.finditer(r"ATOMIC_POSITIONS \(crystal\)\n", txt))
                forces = re.findall(r"Total force =\s+([0-9.]+)", txt)
                if not blocks or not forces:
                    continue
                lines = txt[blocks[-1].end():].splitlines()[:nat]
                try:
                    fr = [[float(x) for x in l.split()[1:4]] for l in lines]
                except Exception:
                    continue
                if len(fr) != nat:
                    continue
                F = float(forces[-1])
                if best is None or F < best[0]:
                    best = (F, fr, f.split("/")[-1] + f" ({len(blocks)} steps, F={F})")
        nt = dict(t)
        if best:
            nt["frac"] = best[1]; nt["seed"] = best[2]
        else:
            nt["seed"] = "initial"
        out.append(nt)
        print(tid, nt["seed"])
    json.dump({"tasks": out}, open(f"{ROOT}/tasks_{tag}_v4.json", "w"))
