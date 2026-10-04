"""Build the web version of the blog post.

    python docs/build_page.py            -> docs/blog.html (page fragment, as published) and docs/index.html (full page)

Numbers come from ledger/claims.csv (written by tools/verify.py --write-ledger) and, for the re-run table, from the
outputs in reproduce/results/ (parsed with tools/qe_parse.py). Figures come from docs/figures.py."""
import csv
import html
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "docs"))
import figures  # noqa: E402
from qe_parse import last_bands, parse_pw, resolve, spin_windows  # noqa: E402

REPO_URL = "https://github.com/spicylemonade/compensated-magnet-ledger"
CLAIMS = list(csv.DictReader(open(ROOT / "ledger" / "claims.csv")))
C = {r["id"]: r for r in CLAIMS}


def lid(*ids):
    return " ".join(f'<a class="lid" href="#c-{i}" title="{html.escape(C[i]["claim"])}">{i}</a>' for i in ids)


# ------------------------------------------------------------------ re-run table (independent container)
def repro_rows():
    base = ROOT / "reproduce" / "results"
    K = ROOT / "materials/KV_Cr_CN6/runs"
    Y = ROOT / "materials/YBaMnFeO5/runs"
    rows = []

    def pw(p):
        return parse_pw(resolve(p))

    def win(p):
        return spin_windows(last_bands(resolve(p)))

    def have(*ps):
        try:
            for p in ps:
                resolve(p)
            return True
        except FileNotFoundError:
            return False

    a = K / "A_pbeu_relax_scf_nscf_Ugrid"
    r = base / "kvcr_pbeu/A_pbeu_relax_scf_nscf_Ugrid"
    if have(r / "KVCr_LCM.repro.out", r / "KVCr_LCM_nscf.repro.out"):
        o, n = pw(a / "KVCr_LCM.out"), pw(r / "KVCr_LCM.repro.out")
        wo, wn = win(a / "KVCr_LCM_nscf.out"), win(r / "KVCr_LCM_nscf.repro.out")
        rows += [
            ("KV[Cr(CN)₆], PBE+U", "net spin moment (μB/cell)", f"{o['total_mag']:.2f}", f"{n['total_mag']:.2f}"),
            ("", "total energy (eV)", f"{o['energy_eV']:.5f}", f"{n['energy_eV']:.5f}"),
            ("", "band gap (eV)", f"{wo['gap']:.3f}", f"{wn['gap']:.3f}"),
            ("", "hole / electron window (eV)", f"{wo['win_VB']:.3f} / {wo['win_CB']:.3f}", f"{wn['win_VB']:.3f} / {wn['win_CB']:.3f}"),
        ]
        if have(r / "KVCr_FM.repro.out"):
            fo, fn = pw(a / "KVCr_FM.out"), pw(r / "KVCr_FM.repro.out")
            rows.append(("", "E(FM) − E(compensated) (meV/ion)", f"{(fo['energy_eV'] - o['energy_eV']) / 2 * 1000:.1f}",
                         f"{(fn['energy_eV'] - n['energy_eV']) / 2 * 1000:.1f}"))
    else:
        rows.append(("KV[Cr(CN)₆], PBE+U", "", "", "running"))
    h = base / "kvcr_hse/C_hse06_ideal_LCM/hse_lcm.repro.out"
    if have(h):
        wo, wn = win(K / "C_hse06_ideal_LCM/hse_lcm.out"), win(h)
        rows.append(("KV[Cr(CN)₆], HSE06", "gap; hole / electron window (eV)",
                     f"{wo['gap']:.3f}; {wo['win_VB']:.3f} / {wo['win_CB']:.3f}", f"{wn['gap']:.3f}; {wn['win_VB']:.3f} / {wn['win_CB']:.3f}"))
    else:
        rows.append(("KV[Cr(CN)₆], HSE06", "gap; hole / electron window (eV)", "2.091; 2.636 / 1.567", "running"))
    yb = base / "ybmfo_stack_U44"
    if have(yb / "U44_G/single.repro.out", yb / "U44_Yflip/single.repro.out"):
        go, yo = pw(Y / "H_pbeu_stacking_vs_U/U44_G/single.out"), pw(Y / "H_pbeu_stacking_vs_U/U44_Yflip/single.out")
        gn, yn = pw(yb / "U44_G/single.repro.out"), pw(yb / "U44_Yflip/single.repro.out")
        rows.append(("YBaMnFeO₅, PBE+U", "closest competitor above ground state (meV/ion)",
                     f"{(yo['energy_eV'] - go['energy_eV']) / 8 * 1000:.2f}", f"{(yn['energy_eV'] - gn['energy_eV']) / 8 * 1000:.2f}"))
    else:
        rows.append(("YBaMnFeO₅, PBE+U", "closest competitor above ground state (meV/ion)", "8.88", "running"))
    s = base / "kvcr_pbeu_sssp/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.repro.out"
    if have(s):
        ws = win(s)
        rows.append(("Robustness: other pseudopotentials (SSSP 1.3)", "gap; hole / electron window (eV)",
                     f"{1.953:.3f}; {2.023:.3f} / {1.146:.3f}", f"{ws['gap']:.3f}; {ws['win_VB']:.3f} / {ws['win_CB']:.3f}"))
    out = []
    for a_, b_, c_, d_ in rows:
        out.append(f"<tr><th scope=\"row\">{a_}</th><td>{b_}</td><td class=\"num\">{c_}</td><td class=\"num\">{d_}</td></tr>")
    return "\n".join(out)


# ------------------------------------------------------------------ figure 4: KV[Cr(CN)6] windows under real-world conditions (HTML dot plot)
DOTS = [  # condition, hole (HSE, PBE+U), electron (HSE, PBE+U), ledger ids
    ("Ideal crystal", (2.64, 2.02), (1.57, 1.15), ("K16", "K07"), ("K17", "K08")),
    ("With water (·2H₂O)", (2.43, 0.93), (1.42, 0.92), ("K19", "K21"), ("K20", None)),
    ("Water-filled vacancy", (2.80, 2.04), (0.74, 0.47), ("K24", None), ("K23", None)),
]


def dotplot():
    def track(hse, pbe, ids):
        lo, hi = sorted((hse, pbe))
        t = lambda v: f"{v / 3.0 * 100:.2f}%"
        tip_h = f"HSE06: {hse:.2f} eV" + (f" (ledger {ids[0]})" if ids[0] else "")
        tip_p = f"PBE+U: {pbe:.2f} eV" + (f" (ledger {ids[1]})" if ids[1] else "")
        return (f'<div class="track"><span class="seg" style="left:{t(lo)};width:{(hi - lo) / 3.0 * 100:.2f}%"></span>'
                f'<span class="dot dot-pbe" style="left:{t(pbe)}" title="{tip_p}"><span class="dv">{pbe:.2f}</span></span>'
                f'<span class="dot dot-hse" style="left:{t(hse)}" title="{tip_h}"><span class="dv">{hse:.2f}</span></span></div>')

    panels = []
    for k, name in ((1, "Hole window"), (2, "Electron window")):
        rows = []
        for cond, hw, ew, hid, eid in DOTS:
            v, ids = (hw, hid) if k == 1 else (ew, eid)
            rows.append(f'<div class="drow"><span class="dlabel">{cond}</span>{track(v[0], v[1], ids)}</div>')
        axis = ('<div class="drow daxis"><span class="dlabel"></span><div class="track">'
                + "".join(f'<span class="tick" style="left:{x / 3 * 100:.2f}%">{x if x < 3 else "3 eV"}</span>' for x in (0, 1, 2, 3))
                + '</div></div>')
        panels.append(f'<div class="dpanel"><h4>{name}</h4>{"".join(rows)}{axis}</div>')
    legend = ('<div class="dlegend"><span><i class="dot dot-hse"></i>HSE06 (slower, usually more accurate)</span>'
              '<span><i class="dot dot-pbe"></i>PBE+U (faster)</span></div>')
    return f'<div class="dots" role="img" aria-label="Spin windows of KV[Cr(CN)6] under three conditions, comparing two methods">{legend}{"".join(panels)}</div>'


def ledger_table():
    rows = []
    for r in CLAIMS:
        st = r["status"]
        badge = {"PASS": "ok", "FAIL": "bad"}.get(st, "na")
        files = "<br>".join(f"<code>{html.escape(f.replace('materials/', ''))}</code>" for f in r["files"].split("; "))
        rows.append(f'<tr id="c-{r["id"]}"><td class="mono">{r["id"]}</td><td>{html.escape(r["material"])}</td>'
                    f'<td class="mono">{r["tier"]}</td><td>{html.escape(r["claim"])}</td>'
                    f'<td class="num">{html.escape(r["recorded_value"])} <span class="u">{html.escape(r["unit"])}</span></td>'
                    f'<td class="num">{html.escape(r["recomputed_value"])}</td><td><span class="st st-{badge}">{"pass" if st == "PASS" else ("fail" if st == "FAIL" else "n/a")}</span></td>'
                    f'<td class="files">{files}</td></tr>')
    return "\n".join(rows)


CSS = r"""
/* Layout: one reading column (~66ch); figures may run wider; ledger IDs ride inline as small mono tags that jump to the ledger table */
:root{
  --paper:#f5f6f8; --ink:#161a22; --ink-2:#454c5a; --ink-3:#5f6676; --rule:#d8dce3; --panel:#eceff3; --panel-2:#e3e7ed;
  --spin-up:#b07a00; --spin-dn:#23508e; --spin-up-tint:#f3e6c6; --spin-dn-tint:#d9e3f2;
  --ok:#2f6b3a; --bad:#a3302a;
  --f-display:"Schibsted Grotesk","Helvetica Neue",Arial,sans-serif;
  --f-body:"Source Serif 4","Iowan Old Style",Georgia,serif;
  --f-mono:"IBM Plex Mono",ui-monospace,"SFMono-Regular",Menlo,monospace;
}
@media (prefers-color-scheme: dark){ :root:not([data-theme="light"]){
  --paper:#12151b; --ink:#e8eaef; --ink-2:#bfc4cf; --ink-3:#9aa0ac; --rule:#2b303a; --panel:#1a1e26; --panel-2:#222731;
  --spin-up:#b5872a; --spin-dn:#5b8ad6; --spin-up-tint:#3a2f17; --spin-dn-tint:#1c2a40; --ok:#79b98a; --bad:#e2817b; color-scheme:dark } }
:root[data-theme="dark"]{
  --paper:#12151b; --ink:#e8eaef; --ink-2:#bfc4cf; --ink-3:#9aa0ac; --rule:#2b303a; --panel:#1a1e26; --panel-2:#222731;
  --spin-up:#b5872a; --spin-dn:#5b8ad6; --spin-up-tint:#3a2f17; --spin-dn-tint:#1c2a40; --ok:#79b98a; --bad:#e2817b; color-scheme:dark }
html{-webkit-text-size-adjust:100%}
body{background:var(--paper); color:var(--ink); font-family:var(--f-body); font-size:18px; line-height:1.62; margin:0}
.wrap{max-width:46rem; margin:0 auto; padding-inline:20px; padding-block:40px 80px}
.wide{max-width:min(60rem, 100%); margin-inline:auto}
h1,h2,h3,h4{font-family:var(--f-display); text-wrap:balance; color:var(--ink); line-height:1.15}
h1{font-size:clamp(2.1rem,6vw,3.3rem); font-weight:800; letter-spacing:-0.02em; margin:0}
h2{font-size:1.6rem; font-weight:700; margin:2.6em 0 .5em; letter-spacing:-0.01em}
h3{font-size:1.15rem; font-weight:700; margin:1.8em 0 .4em}
h4{font-size:.95rem; font-weight:700; margin:0 0 .6em; letter-spacing:.01em}
p{margin:0 0 1em}
a{color:inherit; text-decoration-color:var(--spin-dn); text-underline-offset:3px}
a:hover{color:var(--spin-dn)}
a:focus-visible,button:focus-visible,summary:focus-visible{outline:2px solid var(--spin-dn); outline-offset:2px; border-radius:3px}
.eyebrow{font-family:var(--f-mono); font-size:.75rem; letter-spacing:.12em; text-transform:uppercase; color:var(--ink-3); margin:0 0 1rem}
.hero{display:grid; grid-template-columns:1fr auto; gap:24px; align-items:end}
.hero > div{min-width:0}
.dek{font-size:1.2rem; color:var(--ink-2); margin:.8rem 0 0; max-width:36rem}
.eq{display:flex; align-items:center; gap:10px; font-family:var(--f-display); font-weight:700; font-size:1.3rem; color:var(--ink-2); padding-bottom:6px}
.chip{display:inline-flex; align-items:center; gap:6px; padding:6px 12px; border-radius:999px; font-size:1.05rem}
.chip-up{background:var(--spin-up-tint); color:var(--ink); box-shadow:inset 0 0 0 1.5px var(--spin-up)}
.chip-dn{background:var(--spin-dn-tint); color:var(--ink); box-shadow:inset 0 0 0 1.5px var(--spin-dn)}
.chip b{font-size:1.25rem; line-height:1}
.chip-up b{color:var(--spin-up)} .chip-dn b{color:var(--spin-dn)}
.byline{display:flex; flex-wrap:wrap; gap:6px 18px; font-family:var(--f-display); font-size:.88rem; color:var(--ink-3); margin:1.4rem 0 2.2rem; padding-top:1rem; border-top:1px solid var(--rule)}
.byline strong{color:var(--ink-2); font-weight:600}
.lid{font-family:var(--f-mono); font-size:.68em; vertical-align:.12em; padding:1px 5px; border-radius:4px; background:var(--panel-2); color:var(--ink-2); text-decoration:none; white-space:nowrap}
.lid:hover{background:var(--spin-dn-tint); color:var(--ink)}
.up{color:var(--spin-up); font-weight:700} .dn{color:var(--spin-dn); font-weight:700}
ul,ol{padding-left:1.3em; margin:0 0 1.1em} li{margin:.3em 0}
.wish li::marker{color:var(--ink-3)}
figure{margin:2rem 0}
figcaption{font-family:var(--f-display); font-size:.86rem; line-height:1.45; color:var(--ink-3); margin-top:.7rem; max-width:44rem}
.fig-svg{display:block; width:100%; height:auto}
.scroll{overflow-x:auto; -webkit-overflow-scrolling:touch}
.scroll .fig-svg{min-width:560px}
.panels{display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:12px}
.panels .fig-svg{background:var(--panel); border-radius:10px}
@media (max-width:760px){ .panels{grid-template-columns:minmax(0,1fr)} .panels .fig-svg{max-width:340px; margin-inline:auto} }
/* SVG classes (colours from tokens so both themes work) */
.f-text{font:13px var(--f-display); fill:var(--ink)}
.f-title{font:700 16px var(--f-display); fill:var(--ink)}
.f-sub,.f-small{font:12px var(--f-display); fill:var(--ink-3)}
.f-strong{font-weight:700} .f-strong-t{font:700 13px var(--f-display); fill:var(--ink)}
.f-atom{font:600 11px var(--f-display); fill:var(--ink)}
.s-up{stroke:var(--spin-up)} .s-up-fill{fill:var(--spin-up)} .s-dn{stroke:var(--spin-dn)} .s-dn-fill{fill:var(--spin-dn)}
.f-stroke3{stroke-width:3; stroke-linecap:round} .f-stroke2{stroke-width:1.5}
.atom-n{fill:var(--paper); stroke:var(--ink-3); stroke-width:1.2}
.atom-up{fill:var(--spin-up-tint); stroke:var(--spin-up); stroke-width:1.2}
.atom-dn{fill:var(--spin-dn-tint); stroke:var(--spin-dn); stroke-width:1.2}
.vb-up{fill:var(--spin-up); fill-opacity:.38; stroke:var(--spin-up); stroke-width:1.2}
.cb-up{fill:var(--spin-up); fill-opacity:.07; stroke:var(--spin-up); stroke-width:1.2; stroke-dasharray:3 3}
.vb-dn{fill:var(--spin-dn); fill-opacity:.38; stroke:var(--spin-dn); stroke-width:1.2}
.cb-dn{fill:var(--spin-dn); fill-opacity:.07; stroke:var(--spin-dn); stroke-width:1.2; stroke-dasharray:3 3}
.f-axis{stroke:var(--ink-3); stroke-width:1} .f-guide{stroke:var(--ink-3); stroke-width:1; stroke-dasharray:4 4}
.f-ink{stroke:var(--ink)} .f-ink-fill{fill:var(--ink)}
.bar-a{fill:var(--ink); fill-opacity:.78} .bar-b{fill:var(--ink); fill-opacity:.18; stroke:var(--ink); stroke-width:1} .bar-m{fill:var(--ink-3); fill-opacity:.55}
/* dot plot */
.dots{background:var(--panel); border-radius:10px; padding:18px 18px 10px; font-family:var(--f-display)}
.dlegend{display:flex; flex-wrap:wrap; gap:6px 22px; font-size:.84rem; color:var(--ink-2); margin-bottom:14px}
.dlegend span{display:inline-flex; align-items:center; gap:8px}
.dlegend .dot{position:static; transform:none}
.dpanel{margin-bottom:12px}
.drow{display:grid; grid-template-columns:10.5rem minmax(0,1fr); gap:12px; align-items:center; min-height:62px}
.dlabel{font-size:.86rem; color:var(--ink-2)}
.track{position:relative; height:2px; background:var(--rule); margin:0 18px 0 6px}
.seg{position:absolute; top:-1px; height:4px; background:var(--ink-3); opacity:.55; border-radius:2px}
.dot{position:absolute; top:50%; width:12px; height:12px; border-radius:50%; transform:translate(-50%,-50%); display:inline-block; box-sizing:border-box}
.dot-hse{background:var(--ink); box-shadow:0 0 0 2px var(--panel)}
.dot-pbe{background:var(--panel); border:2px solid var(--ink)}
.dv{position:absolute; bottom:15px; left:50%; transform:translateX(-50%); font-size:.74rem; font-family:var(--f-mono); color:var(--ink-2); white-space:nowrap}
.dot-pbe .dv{bottom:auto; top:15px}
.daxis{min-height:26px} .daxis .track{background:transparent; height:0}
.tick{position:absolute; transform:translateX(-50%); font-size:.72rem; font-family:var(--f-mono); color:var(--ink-3); top:-6px}
.unit{position:absolute; right:-18px; top:-6px; font-size:.72rem; font-family:var(--f-mono); color:var(--ink-3)}
@media (max-width:560px){ .drow{grid-template-columns:minmax(0,1fr); gap:2px; padding-bottom:16px} .dlabel{padding-top:4px} }
/* tables */
.tbl{overflow-x:auto; margin:1.2rem 0 1.6rem}
table{border-collapse:collapse; width:100%; font-family:var(--f-display); font-size:.9rem}
th,td{text-align:left; vertical-align:top; padding:9px 10px; border-bottom:1px solid var(--rule)}
thead th{font-weight:600; color:var(--ink-2); font-size:.84rem; border-bottom:1.5px solid var(--ink-3)}
tbody th{font-weight:600; color:var(--ink)}
td.num{font-variant-numeric:tabular-nums; font-family:var(--f-mono); font-size:.84rem; white-space:nowrap}
.score td:nth-child(2), .score td:nth-child(3){min-width:12rem}
.mono{font-family:var(--f-mono); font-size:.8rem}
.u{color:var(--ink-3); font-family:var(--f-display); font-size:.78rem}
.files code{font-size:.72rem; word-break:break-all}
.st{display:inline-block; font-family:var(--f-mono); font-size:.72rem; padding:1px 7px; border-radius:999px}
.st-ok{color:var(--ok); box-shadow:inset 0 0 0 1px var(--ok)} .st-bad{color:var(--bad); box-shadow:inset 0 0 0 1px var(--bad)} .st-na{color:var(--ink-3); box-shadow:inset 0 0 0 1px var(--ink-3)}
/* code */
pre{position:relative; background:var(--panel); border-radius:8px; padding:14px 16px; overflow-x:auto; font-family:var(--f-mono); font-size:.8rem; line-height:1.55; margin:.6rem 0 1.2rem}
pre button{position:absolute; top:8px; right:8px; font-family:var(--f-display); font-size:.74rem; padding:3px 9px; border-radius:6px; border:1px solid var(--rule); background:var(--paper); color:var(--ink-2); cursor:pointer}
pre button:hover{color:var(--ink)}
code{font-family:var(--f-mono); font-size:.86em}
.levels{display:grid; gap:4px}
.level h3{margin-top:1.2em}
details{border-top:1px solid var(--rule); padding-top:.8rem; margin-top:1.4rem}
summary{font-family:var(--f-display); font-weight:700; cursor:pointer; font-size:1.05rem}
.ledger-table{font-size:.82rem} .ledger-table td{padding:7px 8px}
.ledger-table tr:target{background:var(--spin-dn-tint)}
.fine{font-family:var(--f-display); font-size:.85rem; color:var(--ink-3); border-top:1px solid var(--rule); margin-top:3rem; padding-top:1rem}
.fine li{margin:.25em 0}
.note{font-family:var(--f-display); font-size:.9rem; color:var(--ink-2); background:var(--panel); border-radius:8px; padding:12px 14px; margin:1.2rem 0}
@media (max-width:640px){ body{font-size:17px} .hero{grid-template-columns:minmax(0,1fr)} .eq{padding-bottom:0} .wrap{padding-block:28px 64px} }
@media (prefers-reduced-motion:reduce){ *{scroll-behavior:auto} }
"""

JS = r"""
<script>
document.querySelectorAll('pre[data-copy]').forEach(function(pre){
  var b=document.createElement('button'); b.type='button'; b.textContent='Copy';
  b.addEventListener('click',function(){
    var t=pre.querySelector('code').innerText;
    var done=function(){b.textContent='Copied'; setTimeout(function(){b.textContent='Copy'},1500)};
    var fallback=function(){var r=document.createRange(); r.selectNodeContents(pre.querySelector('code')); var s=getSelection(); s.removeAllRanges(); s.addRange(r); b.textContent='Selected';};
    try{ navigator.clipboard.writeText(t).then(done,fallback);}catch(e){fallback();}
  });
  pre.appendChild(b);
});
</script>
"""


def build():
    f1 = "".join(svg for _, svg in figures.fig1_panels())
    body = BODY
    rep = {
        "{{FIG1}}": f1, "{{FIG2}}": figures.fig2(), "{{FIG3}}": figures.fig3(), "{{FIG4}}": dotplot(),
        "{{REPRO}}": repro_rows(), "{{LEDGER}}": ledger_table(),
        "{{REPO_A}}": (f'<a href="{REPO_URL}">repository on GitHub</a>' if REPO_URL else "repository (its link will be added here once it is published)"),
        "{{REPO_B}}": (f'the <a href="{REPO_URL}">repository</a> (LEDGER.md)' if REPO_URL else "the repository's LEDGER.md"),
        "{{NPASS}}": str(sum(1 for r in CLAIMS if r["status"] == "PASS")),
        "{{NRAW}}": str(sum(1 for r in CLAIMS if r["status"] == "PASS" and r["tier"] == "A")),
        "{{NFAIL}}": str(sum(1 for r in CLAIMS if r["status"] == "FAIL")),
    }
    import re
    body = re.sub(r"\{\{L:([A-Za-z0-9,]+)\}\}", lambda m: lid(*m.group(1).split(",")), body)
    for k, v in rep.items():
        body = body.replace(k, v)
    head = ('<title>Zero-Sum Magnets</title>\n<meta name="description" content="Two room-temperature Luttinger-compensated magnet '
            'candidates, YBaMnFeO5 and KV[Cr(CN)6]: what AI agents computed, how sure we are, and a ledger to check the numbers.">\n'
            '<link rel="preconnect" href="https://fonts.googleapis.com"><link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
            '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=Schibsted+Grotesk:wght@400;600;700;800'
            '&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">\n')
    frag = head + "<style>" + CSS + "</style>\n" + body + JS
    (ROOT / "docs" / "blog.html").write_text(frag)
    full = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            + head + "<style>" + CSS + "</style>\n</head>\n<body>\n" + body + JS + "</body>\n</html>\n")
    (ROOT / "docs" / "index.html").write_text(full)
    print("wrote docs/blog.html and docs/index.html", len(frag), "bytes")


BODY = r"""
<main class="wrap">
<header>
  <p class="eyebrow">Research · October 2026</p>
  <div class="hero">
    <div>
      <h1>Two magnets that add up to zero</h1>
      <p class="dek">One of them has been on a shelf since 1999. What a team of AI agents computed, how sure we are, and where every number can be checked.</p>
    </div>
    <div class="eq" aria-label="chromium spin up plus vanadium spin down equals zero">
      <span class="chip chip-up"><b>↑</b>Cr</span><span>+</span><span class="chip chip-dn"><b>↓</b>V</span><span>= 0</span>
    </div>
  </div>
  <p class="byline"><span><strong>Written by</strong> Geby Jaff</span><span><strong>Computed by</strong> Claude Opus 5.5 agents on cloud computers</span><span><strong>Checked by</strong> a ledger: {{NPASS}} of {{NPASS}} computable numbers re-derive</span></p>
</header>

<p>Most of us know two kinds of magnet, even if we don't know their names. The fridge kind is a <strong>ferromagnet</strong>: trillions of tiny atomic magnets all point the same way, so their pull adds up to something you can feel. The other kind, an <strong>antiferromagnet</strong>, is magnetic on the inside, but neighbouring atomic magnets point in opposite directions and cancel exactly. You can't stick it to anything.</p>
<p>For years, people building the next generation of computer memory have wanted something in between. Over the past few days, a team of AI agents and I went looking for materials that might be it. We ended up with two candidates. One we designed from scratch, which our own calculations suggest may be hard to make. The other was first made in a lab in 1999 and turns out to have this property, something we found no earlier paper pointing out. This post explains what we found, how sure we are, and where all the raw data is.</p>

<h2>A 90-second magnet primer</h2>
<p>Electrons have a property called <em>spin</em> that makes each one a tiny magnet pointing "up" or "down". Spin can carry information. That is the idea behind <strong>spintronics</strong>, which already gives us hard-drive read heads and a kind of memory chip called MRAM.</p>
<ul>
  <li><strong>Ferromagnets</strong> work well for spintronics because the electrons that carry current are <em>sorted by spin</em>: at the energies that matter, there are more up-electrons than down-electrons, or only one kind. But they produce stray magnetic fields that disturb their neighbours, and they are comparatively slow to switch.</li>
  <li><strong>Antiferromagnets</strong> have no stray field and can switch roughly a thousand times faster. But in an ordinary antiferromagnet every up-site has an identical down-site, so the electrons are <em>not</em> sorted by spin. That makes them hard to use for spintronics.</li>
</ul>
<p>There is a third option that physicists have only recently named. Imagine an antiferromagnet in which the up-pointing atoms and the down-pointing atoms are <strong>not equivalent</strong> (for example two different elements, or the same element in two different kinds of site), each carrying exactly the same amount of magnetism. The totals still cancel, so ideally there is no stray field. But because the atoms are not equivalent, the electrons can tell them apart, and they end up <strong>sorted by spin, as in a ferromagnet</strong>.</p>

<figure class="wide">
  <div class="panels">{{FIG1}}</div>
  <figcaption>Three kinds of magnet. Arrows are atomic magnets (<span class="up">↑</span> up, <span class="dn">↓</span> down). The small blocks underneath are the electron energy bands for each spin. In the antiferromagnet the two spins' bands are identical. In the other two they sit at different energies, so the electrons at the band edges have one spin.</figcaption>
</figure>

<p>These are called <strong>Luttinger-compensated magnets</strong>. The name refers to a theorem, Luttinger's, which guarantees that in an insulating material of this kind the cancellation is exact rather than approximate: each spin direction holds a whole number of electrons. (Strictly, that holds for the spin of a perfect crystal near absolute zero; smaller effects such as spin–orbit coupling, and heat, can leave a slight imbalance.) The name comes from a 2022 editorial by the physicist Igor Mazin, and very few real examples are known. The only one confirmed by neutron experiments to be an insulator orders at −225 °C.</p>
<p>The wish list for a useful one is short to write down and hard to satisfy:</p>
<ul class="wish">
  <li><strong>Zero net spin</strong>, fixed by the chemistry rather than by luck.</li>
  <li><strong>Stays magnetic above room temperature.</strong></li>
  <li><strong>A semiconductor</strong>: it has a band gap, like silicon, so you can control how many charge carriers it has.</li>
  <li><strong>Spin-sorted carriers</strong> at both edges of the gap, over an energy range that is large compared with the thermal jiggling at room temperature (0.026 eV).</li>
</ul>
<p>That energy range is the <strong>spin window</strong>: the slice of energy at the edge of the gap where every available electron state has the same spin. A window of 1 eV is about 40 times the room-temperature jiggle, so in a perfect crystal the carriers would stay sorted.</p>

<figure class="wide">
  <div class="scroll">{{FIG2}}</div>
  <figcaption>What a spin window looks like, drawn to scale for ideal KV[Cr(CN)₆] (HSE06 calculation; ledger {{L:K14,K16,K17}}). Blocks stand for filled and empty bands, with their internal structure left out. Both band edges belong to the spin-down channel, and the room-temperature thermal energy is the thin sliver at bottom right.</figcaption>
</figure>

<h2>How the search worked, and who did it</h2>
<p>I set up several Claude agents to run in parallel, each with its own lane of a broader search for unusual magnets. This post comes from the lane that hunted for Luttinger-compensated semiconductors. The agents used <strong>density functional theory (DFT)</strong>, the standard way to compute how electrons arrange themselves in a crystal, through a free program called Quantum ESPRESSO on rented cloud computers. That lane submitted about 750 computing jobs over three days, at two levels of theory:</p>
<ul>
  <li><strong>PBE+U</strong>: a fast approximation with a tunable parameter "U". The agents always checked several values.</li>
  <li><strong>HSE06</strong>: a slower, usually more accurate method, used as the tie-breaker.</li>
</ul>
<p>The agents were set up to argue with themselves. Before running a decisive calculation, they wrote down what result would kill the idea. Separate "referee" agents then tried to tear each claim apart, and several claims were retracted along the way. That record of corrections is part of the repository.</p>

<h2>YBaMnFeO₅: a blueprint that may be hard to build</h2>
<p><strong>The idea.</strong> Take a well-known family of layered oxide crystals and put manganese (Mn) and iron (Fe) on the magnetic sites in a perfect 3D checkerboard, so that every Mn is surrounded by Fe and vice versa. Here Mn²⁺ and Fe³⁺ both have five unpaired electrons, so their atomic magnets are the same size, and the checkerboard makes them point opposite ways. Yttrium, barium and oxygen fill in the rest, all of them cheap and earth-abundant.</p>
<p><strong>On paper, it is close to ideal.</strong> The calculations predict:</p>
<ul>
  <li>zero net spin in the perfect crystal {{L:Y03}};</li>
  <li>a band gap of about <strong>2.35 eV</strong> {{L:Y10}};</li>
  <li>both edges of the gap carrying the same spin {{L:Y11}}, with windows of <strong>1.0 eV and 1.4 eV</strong> {{L:Y12,Y13}};</li>
  <li>magnetic order up to roughly <strong>420 K</strong> (about 145 °C) in the raw simulation, or about <strong>490 K</strong> after calibrating against a known relative {{L:Y22,Y23}};</li>
  <li>a position close to, but not on, the edge of thermodynamic stability: about 14 meV per atom above the most stable mix of competing compounds that were computed, which is typical of compounds that have been made {{L:Y24b}}. (An earlier count, before every competing compound had finished, gave 2.6 meV per atom {{L:Y24}}.)</li>
</ul>
<p><strong>The catch is the checkerboard.</strong> Mn and Fe sit next to each other in the periodic table, are nearly the same size, and differ by one unit of charge. That gives them little reason to keep to their own squares. When the agents simulated how the atoms arrange themselves at different temperatures, the checkerboard melted into a random mix at around <strong>950 K</strong> {{L:Y25,Y26}}. To make this kind of oxide you heat it to roughly 900–1300 °C, and much below that the metal atoms are effectively frozen in place. So by the time the material is cool enough for the checkerboard to be favoured, the atoms may no longer be able to move into it, and standard synthesis would likely give a scrambled crystal.</p>

<figure class="wide">
  <div class="scroll">{{FIG3}}</div>
  <figcaption>Why the useful form of YBaMnFeO₅ may be hard to reach. The whisker on the top bar is the uncertainty of the computed melting point (800–1200 K). The lower bar starts where the metal atoms become mobile, a threshold inferred from a related compound. Ledger {{L:Y25,Y26,Y27}}.</figcaption>
</figure>

<p>This isn't just theory. Every chemically similar compound whose atomic arrangement has been checked came out scrambled, including versions with gadolinium or neodymium in place of yttrium and one with cobalt in place of iron. And scrambling ruins the effect: in the simulation, swapping a single neighbouring Mn/Fe pair pushed states of the opposite spin into the gap and shrank it to almost nothing, 0.01 eV {{L:Y19,Y20}}.</p>
<p><strong>Verdict.</strong> It is a beautiful blueprint and a useful lesson, and the agents' own review downgraded it to a design study, because there is no known way yet to make the ordered crystal. The lesson is that the difference between the two magnetic sublattices has to be <em>enforced by strong chemistry</em>, not left to delicate atomic ordering. That lesson led straight to the second candidate.</p>

<h2>KV[Cr(CN)₆]: hiding in plain sight since 1999</h2>
<p><strong>What it is.</strong> KV[Cr(CN)₆] belongs to the same family as <strong>Prussian blue</strong>, the 300-year-old pigment. Picture a cubic scaffold of cyanide groups, each one carbon atom joined to one nitrogen atom. Chromium atoms hold the carbon ends, vanadium atoms hold the nitrogen ends, and potassium ions sit in the holes. In 1999 the chemists Stephen Holmes and Gregory Girolami made it and found that their sample stays magnetic up to <strong>376 K</strong> (103 °C), or 365 K after it had been heated {{L:K30}}, unusually high for a magnet assembled from molecular building blocks. They designed it so that the vanadium and chromium magnets, three unpaired electrons each, would cancel. They measured almost zero: about 2 % of what you'd get if all the spins lined up {{L:K31}}. Small leftovers like this are common in real samples, where a few building blocks are missing. Nobody has measured its band gap or spin sorting.</p>
<p><strong>What we added.</strong> The agents computed its electronic structure. In an ideal crystal it is a Luttinger-compensated magnet with exactly the spin-sorted structure described above:</p>
<ul>
  <li><strong>zero net spin</strong> in the perfect crystal {{L:K02,K13}}, with a symmetry analysis classifying it as this type of magnet {{L:K29}};</li>
  <li>a <strong>band gap of about 2.1 eV</strong> {{L:K14}};</li>
  <li><strong>both band edges carrying the same spin</strong> {{L:K15}}, with windows of <strong>2.6 eV</strong> for holes and <strong>1.6 eV</strong> for electrons {{L:K16,K17}}, 60 to 100 times the room-temperature jiggle. The cheaper method agrees on the picture, with somewhat smaller windows at every value of U tried {{L:K07,K08,K09,K10,K11}};</li>
  <li>the cyanide bridge <strong>fixing which metal sits where</strong>, because chromium strongly prefers the carbon end and vanadium the nitrogen end. That is exactly the chemical enforcement YBaMnFeO₅ lacked.</li>
</ul>
<p>The fact that vanadium and chromium are different elements is what makes the windows large. When the agents ran the same structure with chromium on both sites, Cr[Cr(CN)₆], the windows shrank to 0.1–0.4 eV and the two band edges took opposite spins {{L:K27,K28}}.</p>
<p>Its zero net magnetism is not new: Holmes and Girolami designed it that way, and a 2008 calculation (Kabalan and colleagues) found it to be an insulator with zero net moment and a gap of about 1 eV. What a literature search (October 2026) did not find is any earlier work resolving which spin sits at its band edges, measuring its spin windows, or describing any Prussian-blue-type compound as a Luttinger-compensated magnet. The closest earlier work, a 2024 study of Cr[Cr(CN)₆], noticed unequal spin-up and spin-down densities but went no further.</p>
<p>It also fills a gap others have pointed out. A 2025 paper that predicted two cyanide Luttinger-compensated semiconductors, Mn(CN)₂ and Co(CN)₂, found that they lose their magnetic order at 210 K and 75 K, and named room-temperature order as the open goal. KV[Cr(CN)₆] has already been made, and its 1999 sample stayed ordered up to 376 K.</p>

<h3>Limits of the result</h3>
<ol>
  <li><strong>The calculations are for a perfect crystal.</strong> The real 1999 material is a powder with water in its holes, and there is a single published report of it.</li>
  <li><strong>The two methods disagree about water.</strong> With water in the model, HSE06 says the effect survives (windows about 2.4 and 1.4 eV), while PBE+U says the hole window shrinks by more than half (chart below). HSE06 is the more reliable of the two here, because PBE+U is known to misplace water's energy levels. The HSE06 water calculation was stopped just short of full convergence.</li>
  <li><strong>Missing building blocks.</strong> Real samples often lack some [Cr(CN)₆] units, and each missing unit adds magnetism, so "exactly zero" depends on getting the composition right {{L:K22}}.</li>
  <li><strong>"Semiconductor" is on paper.</strong> Nobody has measured this compound's band gap, conductivity or spin polarisation. Its electrons move in narrow bands, so charge carriers will be sluggish. Think of a material that holds spin-sorted charges, not a fast transistor material.</li>
  <li><strong>Room temperature is close to its 376 K limit</strong>, so its magnetic order is only about 60 % complete there, which would dilute the effect.</li>
  <li><strong>The general physics is known.</strong> Magnets of this kind have long been understood to have spin-split electrons. What is new is pointing at this specific, real, above-room-temperature compound and putting numbers on it.</li>
</ol>

<figure>
  {{FIG4}}
  <figcaption>Spin windows of KV[Cr(CN)₆] in three situations, from the two methods. Where the two dots sit far apart, the methods disagree. The vacancy values are measured inside the defect model; aligned to the perfect crystal's bands, the PBE+U electron window there drops to 0.07 eV. Hover a dot for its ledger ID.</figcaption>
</figure>


<h2>Scorecard</h2>
<div class="tbl">
<table class="score">
  <thead><tr><th scope="col"></th><th scope="col">YBaMnFeO₅</th><th scope="col">KV[Cr(CN)₆]</th></tr></thead>
  <tbody>
    <tr><th scope="row">Origin</th><td>designed in this project</td><td>made by Holmes &amp; Girolami, 1999</td></tr>
    <tr><th scope="row">Net spin (ideal crystal, 0 K)</th><td>zero</td><td>zero</td></tr>
    <tr><th scope="row">Band gap (HSE06)</th><td class="num">2.35 eV</td><td class="num">2.09 eV</td></tr>
    <tr><th scope="row">Spin windows, holes / electrons</th><td class="num">1.0 / 1.4 eV</td><td class="num">2.6 / 1.6 eV</td></tr>
    <tr><th scope="row">Magnetic up to</th><td>≈ 420 K raw, ≈ 490 K calibrated (predicted)</td><td><strong>376 K (measured on the hydrated powder)</strong></td></tr>
    <tr><th scope="row">Can it be made?</th><td>may be hard: the atoms tend to scramble the checkerboard</td><td>yes, once, as a hydrated powder</td></tr>
    <tr><th scope="row">Biggest open question</th><td>is there any route to the ordered crystal?</td><td>does the spin sorting survive in real, wet, imperfect samples?</td></tr>
    <tr><th scope="row">Ever measured spin-sorted?</th><td>no</td><td>no</td></tr>
  </tbody>
</table>
</div>

<h2>How to check every number</h2>
<p>Everything is packaged as a <strong>computational ledger</strong> in a {{REPO_A}}. It holds the exact input files for nearly 900 calculations with their raw, unedited outputs, the relaxed crystal structures, the analysis scripts, and a list of every claim with the files it came from. The list is reproduced at the bottom of this page. It also records the mistakes the agents caught and the caveats they flagged. There are three levels of checking:</p>
<div class="levels">
  <div class="level">
    <h3>1. Check the arithmetic: seconds, on a laptop</h3>
    <p>Recompute {{NRAW}} numbers directly from the raw outputs and check the rest against the included scripts' results and recorded analysis files. Today it reports {{NPASS}} pass and {{NFAIL}} fail.</p>
<pre data-copy><code>pip install numpy
python tools/verify.py</code></pre>
  </div>
  <div class="level">
    <h3>2. Re-run the models: minutes, on a laptop</h3>
    <p>Redo the magnetic-ordering simulation (recorded 417 K; the re-run gives 414 K) and the checkerboard-melting simulation that rules out YBaMnFeO₅ (it reproduces the recorded 915–965 K exactly).</p>
<pre data-copy><code>pip install numpy pymatgen numba
python materials/YBaMnFeO5/exchange_fit_and_Neel_T/rerun_fit_and_mc.py
bash   materials/YBaMnFeO5/cation_order_cluster_expansion/rerun_torder.sh 3</code></pre>
  </div>
  <div class="level">
    <h3>3. Re-run the quantum calculations: hours, on a cluster or in the cloud</h3>
    <p>Download the same pseudopotentials, check them against the original checksums, and run the original input files with Quantum ESPRESSO 7.5.</p>
<pre data-copy><code>conda install -c conda-forge qe=7.5 openmpi
bash reproduce/get_pseudos.sh pseudo
NP=16 NK=8 PSEUDO=$PWD/pseudo bash reproduce/run_qe.sh \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM.in \
  materials/KV_Cr_CN6/runs/A_pbeu_relax_scf_nscf_Ugrid/KVCr_LCM_nscf.in</code></pre>
  </div>
</div>
<p>We ran level 3 ourselves in a fresh cloud machine: a clean install, pseudopotentials downloaded from the public source, and the inputs from the repository.</p>
<div class="tbl">
<table>
  <thead><tr><th scope="col">Calculation</th><th scope="col">Quantity</th><th scope="col">Original</th><th scope="col">Fresh re-run</th></tr></thead>
  <tbody>
{{REPRO}}
  </tbody>
</table>
</div>
<p class="note">The last row is a robustness check, not a reproduction. It uses a different family of pseudopotentials, so small differences are expected. The picture is unchanged: zero net spin, both band edges in the same spin channel, and windows within 0.2 eV.</p>
<p>If you find an error, open an issue on the repository.</p>

<h2>What would settle it</h2>
<p>The most decisive test is for KV[Cr(CN)₆], because it already exists.</p>
<ul>
  <li><strong>Remake it</strong> and measure its composition and magnetism carefully. It is the cheapest step, and it confirms the starting point.</li>
  <li><strong>Element-specific X-ray magnetic measurements</strong> (XMCD at the vanadium and chromium edges) and magneto-optics on samples with near-zero magnetisation. These would show the two sublattices cancelling while the electronic signal remains.</li>
  <li><strong>Spin-resolved photoemission</strong>, which knocks electrons out with light and measures their spin. The prediction is that the top ~2 eV of filled states are essentially all one spin, flipping when the magnetic order is reversed.</li>
</ul>

<h2>What I took away</h2>
<p>AI agents are good at breadth, and at being systematically sceptical when they are set up to be. In three days they ran hundreds of calculations, checked the literature, and set aside their own favourite idea, and only then did they find the better candidate in a 1999 paper. The next step belongs to the lab. Nothing here has been measured yet, and the open question, whether real KV[Cr(CN)₆] keeps its spin-sorted edges, needs an experiment. The ledger gives anyone running those experiments the exact numbers to test against.</p>

<details id="ledger">
  <summary>The ledger: every claim, its value, and the raw files it comes from</summary>
  <p class="note">Tier A: re-derived from raw Quantum ESPRESSO outputs. L2: re-run end to end by an included script. B: read from a recorded analysis file. E: experimental or literature value. "Recorded" is what the agents wrote down. "Recomputed" is what <code>tools/verify.py</code> gets from the files.</p>
  <div class="tbl">
  <table class="ledger-table">
    <thead><tr><th scope="col">ID</th><th scope="col">Material</th><th scope="col">Tier</th><th scope="col">Claim</th><th scope="col">Recorded</th><th scope="col">Recomputed</th><th scope="col">Status</th><th scope="col">Raw files</th></tr></thead>
    <tbody>
{{LEDGER}}
    </tbody>
  </table>
  </div>
</details>

<footer class="fine">
  <ul>
    <li>Computations: Quantum ESPRESSO 7.5 with PseudoDojo pseudopotentials on Modal cloud CPUs, set up, run and analysed by Claude agents (Anthropic), 1–4 October 2026.</li>
    <li>Experimental facts about KV[Cr(CN)₆]: S. M. Holmes and G. S. Girolami, <em>J. Am. Chem. Soc.</em> 121, 5593 (1999).</li>
    <li>The term "Luttinger-compensated": I. Mazin, <em>Phys. Rev. X</em> editorial (2022). Closest prior computational work on this family: Schart et al., <em>Inorg. Chem.</em> 63, 22856 (2024).</li>
    <li>Full list of caveats, files and scripts: {{REPO_B}}.</li>
  </ul>
</footer>
</main>
"""

if __name__ == "__main__":
    build()
