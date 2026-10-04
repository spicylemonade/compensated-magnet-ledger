"""SVG figures for the blog post. Each figure is written two ways:
  - inline (class-styled, colours come from the page's CSS tokens) for docs/blog.html;
  - standalone files in docs/img/ with an embedded light palette, for BLOG.md on GitHub.
All quantitative marks are drawn to scale from ledger values (IDs noted next to each number)."""

# ---------------------------------------------------------------- helpers
def _t(x, y, s, cls="f-text", anchor="start", extra=""):
    return f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}"{extra}>{s}</text>'


def _arrow(x, up, cls):
    if up:
        return (f'<line x1="{x}" y1="106" x2="{x}" y2="74" class="{cls} f-stroke3"/>'
                f'<polygon points="{x-6},77 {x+6},77 {x},64" class="{cls}-fill"/>')
    return (f'<line x1="{x}" y1="64" x2="{x}" y2="96" class="{cls} f-stroke3"/>'
            f'<polygon points="{x-6},93 {x+6},93 {x},106" class="{cls}-fill"/>')


# ---------------------------------------------------------------- figure 1: three kinds of magnet
PANELS = {
    "fm": dict(title="Ferromagnet", sub="the fridge-magnet kind", atoms=[("Fe", 1, "n")] * 6,
               net="large", sorted_="yes", bands="split"),
    "afm": dict(title="Antiferromagnet", sub="one element, cancels out", atoms=[("Mn", 1, "n"), ("Mn", 0, "n")] * 3,
                net="zero", sorted_="no", bands="same"),
    "lcm": dict(title="Luttinger-compensated", sub="two elements, cancels out", atoms=[("Cr", 1, "up"), ("V", 0, "dn")] * 3,
                net="zero", sorted_="yes", bands="split"),
}


def panel(kind):
    p = PANELS[kind]
    out = [_t(150, 24, p["title"], "f-title", "middle"), _t(150, 44, p["sub"], "f-sub", "middle")]
    for i, (el, up, tint) in enumerate(p["atoms"]):
        x = 50 + i * 40
        out.append(_arrow(x, up, "s-up" if up else "s-dn"))
        out.append(f'<circle cx="{x}" cy="128" r="14" class="atom-{tint}"/>')
        out.append(_t(x, 132, el, "f-atom", "middle"))
    out.append(_t(150, 176, f'Net magnetism: <tspan class="f-strong">{p["net"]}</tspan>', "f-text", "middle"))
    # band strips: two spin columns; bottom = filled (valence), top = empty (conduction)
    if p["bands"] == "same":
        geo = {"up": (262, 240), "dn": (262, 240)}
    else:
        geo = {"up": (272, 222), "dn": (254, 238)}
    for col, x in (("up", 100), ("dn", 166)):
        vb_top, cb_bot = geo[col]
        out.append(f'<rect x="{x}" y="{vb_top}" width="34" height="{296 - vb_top}" class="vb-{col}"/>')
        out.append(f'<rect x="{x}" y="206" width="34" height="{cb_bot - 206}" class="cb-{col}"/>')
        out.append(_t(x + 17, 312, "spin ↑" if col == "up" else "spin ↓", "f-small", "middle"))
    out.append(_t(92, 226, "empty", "f-small", "end"))
    out.append(_t(92, 288, "filled", "f-small", "end"))
    out.append(_t(150, 334, f'Electrons sorted by spin: <tspan class="f-strong">{p["sorted_"]}</tspan>', "f-text", "middle"))
    return "\n".join(out)


def fig1_panels():
    """three separate inline SVGs (they wrap to one column on phones)"""
    labels = {"fm": "A ferromagnet: identical atomic magnets all point up, so the net magnetism is large and the spin-up and spin-down electron bands sit at different energies.",
              "afm": "An antiferromagnet: identical atoms alternate up and down, so the net is zero and the two spin bands are identical.",
              "lcm": "A Luttinger-compensated magnet: two different elements alternate up and down with equal magnets, so the net is zero but the two spin bands still sit at different energies."}
    return [(k, f'<svg viewBox="0 0 300 344" role="img" aria-label="{labels[k]}" class="fig-svg">{panel(k)}</svg>') for k in ("fm", "afm", "lcm")]


# ---------------------------------------------------------------- figure 2: spin window, to scale (KV[Cr(CN)6], HSE06)
def y_of(E):
    return 216 - 40 * E


def fig2():
    gap, wvb, wcb = 2.09, 2.64, 1.57          # ledger K14, K16, K17
    o = []
    o.append(_t(240, 26, "spin ↑ (Cr sites)", "f-text", "middle"))
    o.append(_t(360, 26, "spin ↓ (V sites)", "f-text", "middle"))
    # axis
    o.append('<line x1="150" y1="40" x2="150" y2="352" class="f-axis"/>')
    for E in range(-3, 5):
        y = y_of(E)
        o.append(f'<line x1="146" y1="{y}" x2="150" y2="{y}" class="f-axis"/>')
        o.append(_t(140, y + 4, ("−" if E < 0 else "") + str(abs(E)), "f-small", "end"))
    o.append(_t(104, 196, "energy (eV)", "f-small", "middle", ' transform="rotate(-90 104 196)"'))
    # bands
    o.append(f'<rect x="320" y="{y_of(0)}" width="80" height="{352 - y_of(0):.1f}" class="vb-dn"/>')
    o.append(f'<rect x="320" y="40" width="80" height="{y_of(gap) - 40:.1f}" class="cb-dn"/>')
    o.append(f'<rect x="200" y="{y_of(-wvb):.1f}" width="80" height="{352 - y_of(-wvb):.1f}" class="vb-up"/>')
    o.append(f'<rect x="200" y="40" width="80" height="{y_of(gap + wcb) - 40:.1f}" class="cb-up"/>')
    o.append(_t(240, 54, "empty", "f-small", "middle"))
    o.append(_t(240, 344, "filled", "f-small", "middle"))
    o.append(_t(360, 92, "empty", "f-small", "middle"))
    o.append(_t(360, 300, "filled", "f-small", "middle"))
    # guides at the two band edges
    for E in (0, gap):
        o.append(f'<line x1="196" y1="{y_of(E):.1f}" x2="404" y2="{y_of(E):.1f}" class="f-guide"/>')
    # brackets
    def bracket(x, y1, y2):
        return (f'<line x1="{x}" y1="{y1:.1f}" x2="{x}" y2="{y2:.1f}" class="f-ink f-stroke2"/>'
                f'<line x1="{x-5}" y1="{y1:.1f}" x2="{x+5}" y2="{y1:.1f}" class="f-ink f-stroke2"/>'
                f'<line x1="{x-5}" y1="{y2:.1f}" x2="{x+5}" y2="{y2:.1f}" class="f-ink f-stroke2"/>')
    o.append(bracket(418, y_of(0), y_of(-wvb)))
    o.append(_t(430, (y_of(0) + y_of(-wvb)) / 2 - 4, "hole window", "f-strong-t"))
    o.append(_t(430, (y_of(0) + y_of(-wvb)) / 2 + 13, "2.64 eV, every state spin ↓", "f-small"))
    o.append(bracket(418, y_of(gap), y_of(gap + wcb)))
    o.append(_t(430, (y_of(gap) + y_of(gap + wcb)) / 2 - 4, "electron window", "f-strong-t"))
    o.append(_t(430, (y_of(gap) + y_of(gap + wcb)) / 2 + 13, "1.57 eV, every state spin ↓", "f-small"))
    o.append(bracket(306, y_of(gap), y_of(0)))
    o.append(_t(298, (y_of(gap) + y_of(0)) / 2 - 4, "band gap", "f-strong-t", "end"))
    o.append(_t(298, (y_of(gap) + y_of(0)) / 2 + 13, "2.09 eV", "f-small", "end"))
    # scale: 1 eV vs room-temperature thermal energy
    o.append(_t(470, 320, "for scale", "f-small"))
    o.append('<line x1="476" y1="330" x2="476" y2="370" class="f-ink f-stroke2"/>')
    o.append(_t(484, 354, "1 eV", "f-small"))
    o.append('<rect x="524" y="369" width="12" height="1.2" class="f-ink-fill"/>')
    o.append(_t(542, 364, "0.026 eV", "f-small"))
    o.append(_t(542, 380, "(room temp.)", "f-small"))
    label = ("Band edges of ideal KV[Cr(CN)6] from the HSE06 calculation, drawn to scale: in the spin-down channel the gap is 2.09 eV; "
             "the top 2.64 eV of filled states and the bottom 1.57 eV of empty states contain only spin-down states.")
    return f'<svg viewBox="0 0 640 392" role="img" aria-label="{label}" class="fig-svg">' + "\n".join(o) + "</svg>"


# ---------------------------------------------------------------- figure 3: why the ordered form of YBaMnFeO5 may be hard to make
def x_of(T):
    return 40 + T * 640 / 1700


def fig3():
    o = []
    # room temperature
    o.append(f'<line x1="{x_of(300):.1f}" y1="40" x2="{x_of(300):.1f}" y2="196" class="f-guide"/>')
    o.append(_t(x_of(300) - 4, 36, "room temperature", "f-small", "end"))
    # row A: checkerboard stable (central 950 K, whisker 800-1200 K)
    o.append(_t(44, 54, "Mn/Fe checkerboard is the stable arrangement (computed)", "f-text"))
    o.append(f'<rect x="40" y="62" width="{x_of(950) - 40:.1f}" height="18" rx="3" class="bar-a"/>')
    o.append(f'<line x1="{x_of(800):.1f}" y1="71" x2="{x_of(1200):.1f}" y2="71" class="f-ink f-stroke2"/>')
    for T in (800, 1200):
        o.append(f'<line x1="{x_of(T):.1f}" y1="65" x2="{x_of(T):.1f}" y2="77" class="f-ink f-stroke2"/>')
    o.append(_t(x_of(950), 96, "melts ≈ 950 K", "f-small", "middle"))
    # row B: atoms mobile (freeze-out ~1150 K; synthesis 1173-1573 K)
    o.append(_t(x_of(1150), 120, "metal atoms can rearrange", "f-text"))
    o.append(f'<rect x="{x_of(1150):.1f}" y="128" width="{x_of(1700) - x_of(1150):.1f}" height="18" rx="3" class="bar-b"/>')
    o.append(f'<line x1="{x_of(1173):.1f}" y1="156" x2="{x_of(1573):.1f}" y2="156" class="f-ink f-stroke2"/>')
    for T in (1173, 1573):
        o.append(f'<line x1="{x_of(T):.1f}" y1="151" x2="{x_of(T):.1f}" y2="161" class="f-ink f-stroke2"/>')
    o.append(_t((x_of(1173) + x_of(1573)) / 2, 176, "typical synthesis 1173–1573 K", "f-small", "middle"))
    # the gap between them
    o.append(f'<line x1="{x_of(950):.1f}" y1="137" x2="{x_of(1150) - 4:.1f}" y2="137" class="f-ink f-stroke2"/>')
    o.append(f'<polygon points="{x_of(950):.1f},137 {x_of(950) + 7:.1f},133 {x_of(950) + 7:.1f},141" class="f-ink-fill"/>')
    o.append(f'<polygon points="{x_of(1150) - 2:.1f},137 {x_of(1150) - 9:.1f},133 {x_of(1150) - 9:.1f},141" class="f-ink-fill"/>')
    o.append(_t((x_of(950) + x_of(1150)) / 2, 128, "≈ 200 K", "f-small", "middle"))
    # magnetic order range
    o.append(f'<rect x="{x_of(417):.1f}" y="191" width="{x_of(490) - x_of(417):.1f}" height="10" rx="2" class="bar-m"/>')
    o.append(_t(x_of(453), 238, "magnetic order lost at ≈ 420–490 K (predicted)", "f-small", "middle"))
    # axis
    o.append('<line x1="40" y1="196" x2="680" y2="196" class="f-axis"/>')
    for T in (0, 300, 600, 900, 1200, 1500):
        o.append(f'<line x1="{x_of(T):.1f}" y1="196" x2="{x_of(T):.1f}" y2="201" class="f-axis"/>')
        o.append(_t(x_of(T), 216, f"{T}", "f-small", "middle"))
    o.append(_t(680, 216, "K", "f-small", "end"))
    label = ("Temperature scale for YBaMnFeO5: the Mn/Fe checkerboard is stable only below about 950 K, but the metal atoms can only "
             "rearrange above about 1150 K, and syntheses run at 1173 to 1573 K, so the ordered form may be hard to reach by heating and cooling.")
    return f'<svg viewBox="0 0 720 248" role="img" aria-label="{label}" class="fig-svg">' + "\n".join(o) + "</svg>"


# ---------------------------------------------------------------- standalone files for GitHub
FILE_STYLE = """<style>
.f-text{font:13px/1.2 'Söhne',-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;fill:#191917}
.f-title{font:700 16px 'Söhne',-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;fill:#191917}
.f-sub,.f-small{font:12px 'Söhne',-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;fill:#56524e}
.f-strong{font-weight:700}.f-strong-t{font:700 13px 'Söhne',-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;fill:#191917}
.f-atom{font:600 11px 'Söhne',-apple-system,'Segoe UI',Helvetica,Arial,sans-serif;fill:#191917}
.s-up{stroke:#9c9691}.s-up-fill{fill:#9c9691}.s-dn{stroke:#386b46}.s-dn-fill{fill:#386b46}
.f-stroke3{stroke-width:3;stroke-linecap:round}.f-stroke2{stroke-width:1.5}
.atom-n{fill:#f3f2f2;stroke:#9c9691;stroke-width:1.2}.atom-up{fill:#f3f2f2;stroke:#9c9691;stroke-width:1.2}.atom-dn{fill:#d0e7d6;stroke:#386b46;stroke-width:1.2}
.vb-up{fill:#9c9691;fill-opacity:.35;stroke:#9c9691;stroke-width:1.2}.cb-up{fill:#9c9691;fill-opacity:.07;stroke:#9c9691;stroke-width:1.2;stroke-dasharray:3 3}
.vb-dn{fill:#386b46;fill-opacity:.35;stroke:#386b46;stroke-width:1.2}.cb-dn{fill:#386b46;fill-opacity:.07;stroke:#386b46;stroke-width:1.2;stroke-dasharray:3 3}
.f-axis{stroke:#9c9691;stroke-width:1}.f-guide{stroke:#9c9691;stroke-width:1;stroke-dasharray:4 4}
.f-ink{stroke:#191917}.f-ink-fill{fill:#191917}
.bar-a{fill:#386b46}.bar-b{fill:#d0e7d6;stroke:#488959;stroke-width:1}.bar-m{fill:#9c9691;fill-opacity:.7}
</style>"""


def standalone(svg, w, h):
    inner = svg.split(">", 1)[1].rsplit("</svg>", 1)[0]
    vb = svg.split('viewBox="')[1].split('"')[0]
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="{vb}" width="{w}" height="{h}">{FILE_STYLE}'
            f'<rect x="0" y="0" width="100%" height="100%" rx="10" fill="#ffffff"/>{inner}</svg>\n')


if __name__ == "__main__":
    import pathlib
    d = pathlib.Path(__file__).resolve().parent / "img"
    d.mkdir(exist_ok=True)
    p = fig1_panels()
    inner = "".join(f'<g transform="translate({i * 310},0)">{panel(k)}</g>' for i, (k, _) in enumerate(p))
    combo = f'<svg viewBox="0 0 920 344" role="img" aria-label="Three kinds of magnet">{inner}</svg>'
    (d / "three_magnets.svg").write_text(standalone(combo, 920, 344))
    (d / "spin_window.svg").write_text(standalone(fig2(), 640, 392))
    (d / "ybamnfeo5_temperatures.svg").write_text(standalone(fig3(), 720, 248))
    print("wrote", [f.name for f in d.iterdir()])
