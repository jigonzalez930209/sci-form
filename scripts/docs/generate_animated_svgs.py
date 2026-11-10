#!/usr/bin/env python3
"""Generate the animated, theme-aware SVG diagrams used by the documentation.

The diagrams deliberately avoid a filled background: they inherit the page
surface. Colours are CSS custom properties defined on ``.sf-d`` so the same
file renders correctly in light and dark mode. All motion is CSS/SMIL based
and is disabled by ``prefers-reduced-motion``.

Run from the repository root:

    python3 scripts/docs/generate_animated_svgs.py

The script rewrites ``docs/public/svg/*.svg`` in place. It is deterministic,
so the generated files are reproducible.
"""

from __future__ import annotations

import math
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs/public/svg"

# ---------------------------------------------------------------------------
# Shared palette / animation CSS
# ---------------------------------------------------------------------------

CSS = """
.sf-d{
  --sf-ink:#0f172a; --sf-muted:#5b6b82; --sf-line:#9aa9bd;
  --sf-card:rgba(255,255,255,.62); --sf-card2:rgba(255,255,255,.34);
  --sf-accent:#2563eb; --sf-accent2:#0d9488; --sf-warm:#d97706;
  --sf-violet:#7c3aed; --sf-rose:#e11d48;
  --sf-soft:rgba(37,99,235,.12); --sf-glow:rgba(37,99,235,.20);
  font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;
}
html.dark .sf-d{
  --sf-ink:#e6edf7; --sf-muted:#9fb0c7; --sf-line:#5b6b82;
  --sf-card:rgba(15,23,42,.52); --sf-card2:rgba(15,23,42,.28);
  --sf-accent:#60a5fa; --sf-accent2:#2dd4bf; --sf-warm:#fbbf24;
  --sf-violet:#a78bfa; --sf-rose:#fb7185;
  --sf-soft:rgba(96,165,250,.14); --sf-glow:rgba(96,165,250,.20);
}
@keyframes sf-dash{to{stroke-dashoffset:-28}}
@keyframes sf-pulse{0%,100%{opacity:.55}50%{opacity:1}}
@keyframes sf-halo{0%,100%{transform:scale(.92);opacity:.14}50%{transform:scale(1.08);opacity:.34}}
@keyframes sf-float{0%,100%{transform:translateY(0)}50%{transform:translateY(-5px)}}
@keyframes sf-draw{to{stroke-dashoffset:0}}
@keyframes sf-spin{to{transform:rotate(360deg)}}
@keyframes sf-blink{0%,100%{opacity:.25}50%{opacity:1}}
@keyframes sf-rise{0%{transform:translateY(6px);opacity:0}100%{transform:translateY(0);opacity:1}}
@keyframes sf-scan{0%{transform:translateX(0);opacity:0}12%{opacity:1}88%{opacity:1}100%{transform:translateX(var(--sf-span,260px));opacity:0}}
.sf-flow{stroke-dasharray:4 6;animation:sf-dash 1.3s linear infinite}
.sf-pulse{transform-box:fill-box;transform-origin:center;animation:sf-pulse 2.8s ease-in-out infinite}
.sf-halo{transform-box:fill-box;transform-origin:center;animation:sf-halo 3.4s ease-in-out infinite}
.sf-float{transform-box:fill-box;transform-origin:center;animation:sf-float 4.5s ease-in-out infinite}
.sf-spin{transform-box:fill-box;transform-origin:center;animation:sf-spin 12s linear infinite}
.sf-draw{stroke-dasharray:400;stroke-dashoffset:400;animation:sf-draw 2s ease forwards}
.sf-blink{animation:sf-blink 2.4s ease-in-out infinite}
.sf-rise{animation:sf-rise .9s ease both}
.sf-d1{animation-delay:.25s}.sf-d2{animation-delay:.5s}.sf-d3{animation-delay:.75s}
.sf-d4{animation-delay:1s}.sf-d5{animation-delay:1.25s}.sf-d6{animation-delay:1.5s}
@media (prefers-reduced-motion:reduce){
  .sf-flow,.sf-pulse,.sf-halo,.sf-float,.sf-spin,.sf-draw,.sf-blink,.sf-rise{animation:none!important}
  .sf-draw{stroke-dashoffset:0!important}
  .sf-particle{display:none}
}
"""


def _attrs(**kwargs) -> str:
    out = []
    for key, value in kwargs.items():
        if value is None:
            continue
        key = key.rstrip("_").replace("_", "-")
        out.append(f'{key}="{value}"')
    return " ".join(out)


def esc(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


def text(
    x: float,
    y: float,
    value: str,
    *,
    size: int = 12,
    fill: str = "var(--sf-ink)",
    anchor: str = "middle",
    weight: int = 600,
    mono: bool = False,
    cls: str = "",
    opacity: float | None = None,
    rotate: float | None = None,
) -> str:
    family = (
        'ui-monospace,SFMono-Regular,Menlo,monospace'
        if mono
        else 'ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif'
    )
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
        f'font-weight="{weight}" font-family=\'{family}\' fill="{fill}"'
        f'{f" opacity=\"{opacity}\"" if opacity is not None else ""}'
        f'{f" class=\"{cls}\"" if cls else ""}{transform}>{esc(value)}</text>'
    )


def rect(
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    rx: float = 12,
    fill: str = "var(--sf-card)",
    stroke: str = "var(--sf-line)",
    sw: float = 1.4,
    cls: str = "",
    opacity: float | None = None,
) -> str:
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{sw}"'
        f'{f" class=\"{cls}\"" if cls else ""}'
        f'{f" opacity=\"{opacity}\"" if opacity is not None else ""}/>'
    )


def circle(
    cx: float,
    cy: float,
    r: float,
    *,
    fill: str = "var(--sf-accent)",
    stroke: str | None = None,
    sw: float = 1.4,
    cls: str = "",
    opacity: float | None = None,
) -> str:
    return (
        f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{fill}"'
        f'{f" stroke=\"{stroke}\" stroke-width=\"{sw}\"" if stroke else ""}'
        f'{f" class=\"{cls}\"" if cls else ""}'
        f'{f" opacity=\"{opacity}\"" if opacity is not None else ""}/>'
    )


def path(
    d: str,
    *,
    fill: str = "none",
    stroke: str = "var(--sf-line)",
    sw: float = 1.6,
    cls: str = "",
    marker: bool = False,
    opacity: float | None = None,
    dash: str | None = None,
) -> str:
    marker_attr = f' marker-end="url(#{MARKER})"' if marker else ""
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return (
        f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" '
        f'stroke-linecap="round" stroke-linejoin="round"'
        f'{marker_attr}{dash_attr}'
        f'{f" class=\"{cls}\"" if cls else ""}'
        f'{f" opacity=\"{opacity}\"" if opacity is not None else ""}/>'
    )


def particle(d: str, *, dur: float = 2.8, delay: float = 0, r: float = 2.6, fill: str = "var(--sf-accent)") -> str:
    return (
        f'<circle r="{r}" fill="{fill}" class="sf-particle" opacity="0">'
        f'<animateMotion dur="{dur}s" begin="{delay}s" repeatCount="indefinite" path="{d}"/>'
        f'<animate attributeName="opacity" values="0;1;1;0" dur="{dur}s" '
        f'begin="{delay}s" repeatCount="indefinite"/></circle>'
    )


def connector(
    d: str,
    *,
    color: str = "var(--sf-line)",
    sw: float = 1.6,
    flow: bool = True,
    move: bool = True,
    dur: float = 2.8,
    delay: float = 0,
    marker: bool = True,
) -> str:
    cls = "sf-flow" if flow else ""
    out = path(d, stroke=color, sw=sw, cls=cls, marker=marker)
    if move:
        out += particle(d, dur=dur, delay=delay, fill=color)
    return out


def arrow_between(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    *,
    color: str = "var(--sf-line)",
    curve: float = 0,
    flow: bool = True,
    move: bool = True,
    delay: float = 0,
    dur: float = 2.8,
    marker: bool = True,
) -> str:
    if curve == 0:
        d = f"M{x1},{y1} L{x2},{y2}"
    else:
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy) or 1
        nx, ny = -dy / length, dx / length
        cx, cy = mx + nx * curve, my + ny * curve
        d = f"M{x1},{y1} Q{cx},{cy} {x2},{y2}"
    return connector(d, color=color, flow=flow, move=move, delay=delay, dur=dur, marker=marker)


def node(
    x: float,
    y: float,
    w: float,
    h: float,
    label: str,
    *,
    accent: str = "var(--sf-accent)",
    sub: str | None = None,
    icon: str = "",
    cls: str = "sf-rise",
    delay: int = 0,
) -> str:
    out = rect(x, y, w, h, fill="var(--sf-card)", stroke="var(--sf-line)")
    out += rect(x + 1.2, y + 1.2, w - 2.4, 3, rx=1.5, fill=accent, stroke="none", opacity=0.9)
    if icon:
        out += icon
    label_y = y + h / 2 + (5 if sub else 8)
    out += text(x + w / 2, label_y, label, size=12, fill="var(--sf-ink)", weight=650, cls=cls)
    if sub:
        out += text(x + w / 2, y + h / 2 + 22, sub, size=9.5, fill="var(--sf-muted)", weight=450, cls=cls)
    if delay:
        out = out.replace('class="sf-rise"', f'class="sf-rise sf-d{delay}"')
    return out


def pill(x: float, y: float, w: float, label: str, *, color: str = "var(--sf-accent)", size: int = 10, h: float = 22) -> str:
    return (
        rect(x, y, w, h, rx=h / 2, fill="var(--sf-soft)", stroke=color, sw=1)
        + text(x + w / 2, y + h / 2 + size * 0.35, label, size=size, fill=color, weight=650)
    )


def molecule_water(cx: float, cy: float, scale: float = 1.0, *, cls: str = "sf-float") -> str:
    o = (cx, cy)
    h1 = (cx - 28 * scale, cy + 24 * scale)
    h2 = (cx + 28 * scale, cy + 24 * scale)
    out = f'<g class="{cls}">'
    out += path(f"M{o[0]},{o[1]} L{h1[0]},{h1[1]}", stroke="var(--sf-line)", sw=2.4)
    out += path(f"M{o[0]},{o[1]} L{h2[0]},{h2[1]}", stroke="var(--sf-line)", sw=2.4)
    out += circle(*o, 15 * scale, fill="var(--sf-rose)", opacity=0.18)
    out += circle(*o, 9 * scale, fill="var(--sf-rose)")
    out += circle(*h1, 6 * scale, fill="var(--sf-accent)")
    out += circle(*h2, 6 * scale, fill="var(--sf-accent)")
    out += text(o[0], o[1] + 3.5, "O", size=int(10 * scale), fill="#fff", weight=700)
    out += text(h1[0], h1[1] + 3, "H", size=int(8 * scale), fill="#fff", weight=700)
    out += text(h2[0], h2[1] + 3, "H", size=int(8 * scale), fill="#fff", weight=700)
    out += "</g>"
    return out


def benzene(cx: float, cy: float, r: float = 46, *, rotate: float = 0, cls: str = "sf-spin") -> str:
    pts = []
    for i in range(6):
        angle = math.radians(rotate + 60 * i - 90)
        pts.append((cx + r * math.cos(angle), cy + r * math.sin(angle)))
    d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
    out = f'<g class="{cls}">'
    out += path(d, stroke="var(--sf-line)", sw=2.2)
    out += circle(cx, cy, r * 0.56, fill="none", stroke="var(--sf-accent)", sw=1.6, opacity=0.75)
    for x, y in pts:
        out += circle(x, y, 4.5, fill="var(--sf-ink)", opacity=0.85)
    out += "</g>"
    return out


def axes(x: float, y: float, w: float, h: float, *, xlabel: str = "", ylabel: str = "") -> str:
    out = path(f"M{x},{y} L{x},{y+h} M{x},{y+h} L{x+w},{y+h}", stroke="var(--sf-line)", sw=1.2)
    out += path(f"M{x},{y+h} l4,-4 M{x},{y+h} l-4,-4", stroke="var(--sf-line)", sw=1.1)
    out += path(f"M{x},{y} l-4,4 M{x},{y} l4,4", stroke="var(--sf-line)", sw=1.1)
    if xlabel:
        out += text(x + w / 2, y + h + 22, xlabel, size=10, fill="var(--sf-muted)", weight=500)
    if ylabel:
        out += text(x - 18, y + h / 2, ylabel, size=10, fill="var(--sf-muted)", weight=500, rotate=-90)
    return out


def polyline(points: list[tuple[float, float]], *, color: str = "var(--sf-accent)", sw: float = 2.4, cls: str = "sf-draw") -> str:
    d = "M" + " L".join(f"{x:.2f},{y:.2f}" for x, y in points)
    return path(d, stroke=color, sw=sw, cls=cls)


def curve_path(points: list[tuple[float, float]], *, color: str = "var(--sf-accent)", sw: float = 2.4, cls: str = "sf-draw") -> str:
    if not points:
        return ""
    d = f"M{points[0][0]:.2f},{points[0][1]:.2f}"
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        mx, my = (x0 + x1) / 2, (y0 + y1) / 2
        d += f" Q{mx:.2f},{my:.2f} {x1:.2f},{y1:.2f}"
    return path(d, stroke=color, sw=sw, cls=cls)


def sparkle(cx: float, cy: float, r: float = 3, color: str = "var(--sf-accent2)") -> str:
    return path(
        f"M{cx},{cy-r} L{cx+r*0.35},{cy-r*0.35} L{cx+r},{cy} "
        f"L{cx+r*0.35},{cy+r*0.35} L{cx},{cy+r} "
        f"L{cx-r*0.35},{cy+r*0.35} L{cx-r},{cy} "
        f"L{cx-r*0.35},{cy-r*0.35} Z",
        fill=color,
        stroke="none",
        cls="sf-pulse",
    )


def begin(slug: str) -> None:
    global MARKER
    MARKER = f"sf-{slug}-arrow"


def svg(slug: str, w: int, h: int, body: str, *, label: str = "") -> str:
    # Marker ids are prefixed per diagram to avoid collisions when several
    # inline SVGs share the same document.
    global MARKER
    MARKER = f"sf-{slug}-arrow"
    return f'''<svg class="sf-d" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" role="img" aria-label="{esc(label)}">
<style>{CSS}</style>
<defs>
  <marker id="{MARKER}" markerWidth="9" markerHeight="7" refX="8" refY="3.5" orient="auto">
    <path d="M0,0 L0,7 L9,3.5 Z" fill="var(--sf-line)"/>
  </marker>
  <radialGradient id="sf-{slug}-glow" cx="50%" cy="45%" r="55%">
    <stop offset="0%" stop-color="var(--sf-accent)" stop-opacity="0.16"/>
    <stop offset="100%" stop-color="var(--sf-accent)" stop-opacity="0"/>
  </radialGradient>
  <linearGradient id="sf-{slug}-line" x1="0" y1="0" x2="1" y2="0">
    <stop offset="0%" stop-color="var(--sf-accent)"/>
    <stop offset="100%" stop-color="var(--sf-accent2)"/>
  </linearGradient>
</defs>
<rect x="0" y="0" width="{w}" height="{h}" fill="url(#sf-{slug}-glow)" stroke="none"/>
{body}
</svg>
'''


def write(slug: str, w: int, h: int, body: str, *, label: str = "") -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{slug}.svg").write_text(svg(slug, w, h, body, label=label), encoding="utf-8")


MARKER = "sf-arrow"


# ---------------------------------------------------------------------------
# Experimental modules
# ---------------------------------------------------------------------------

def diagram_experimental_gsm() -> None:
    begin("experimental-gsm")
    b = []

    def energy_y(x: float) -> float:
        # Two wells at x=150/650 and a saddle at x=390.
        return (
            95
            + 118 * math.exp(-((x - 150) ** 2) / 5200)
            + 118 * math.exp(-((x - 650) ** 2) / 5200)
            - 42 * math.exp(-((x - 390) ** 2) / 5200)
        )

    pts = [(x, energy_y(x)) for x in range(55, 735, 4)]
    b.append(curve_path(pts, color="var(--sf-line)", sw=2.2, cls="sf-draw"))
    b.append(text(130, 245, "reactant", size=11, fill="var(--sf-muted)"))
    b.append(text(390, 28, "transition state", size=11, fill="var(--sf-accent)", weight=700))
    b.append(text(650, 245, "product", size=11, fill="var(--sf-muted)"))
    # Nodes on the path, growing from both ends.
    node_xs = [145, 230, 310, 390, 470, 550, 645]
    for i, x in enumerate(node_xs):
        y = energy_y(x)
        b.append(circle(x, y, 6.5, fill="var(--sf-card)", stroke="var(--sf-accent)", sw=2, cls=f"sf-rise sf-d{min(i+1,6)}"))
    path_d = "M" + " L".join(f"{x},{energy_y(x):.1f}" for x in node_xs)
    b.append(path(path_d, stroke="var(--sf-accent)", sw=2.2, cls="sf-flow"))
    b.append(particle(path_d, dur=3.4))
    # Two growth arrows from the ends.
    b.append(arrow_between(120, 200, 250, 170, color="var(--sf-accent2)", curve=-18, delay=0.2))
    b.append(arrow_between(660, 200, 530, 170, color="var(--sf-accent2)", curve=18, delay=0.6))
    b.append(pill(342, 252, 96, "grow", color="var(--sf-accent2)"))
    write("experimental-gsm", 780, 280, "".join(b), label="A string of states grows from both ends of a reaction path towards the saddle point.")


def diagram_experimental_randnla() -> None:
    begin("experimental-randnla")
    b = []
    # Full matrix.
    x0, y0, cell = 70, 70, 18
    b.append(rect(x0 - 8, y0 - 8, 7 * cell + 16, 7 * cell + 16, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for r in range(7):
        for c in range(7):
            hot = c < 2 and (r + c) % 3 != 0
            b.append(rect(x0 + c * cell, y0 + r * cell, cell - 2, cell - 2, rx=2,
                          fill="var(--sf-accent)" if hot else "var(--sf-soft)", stroke="none",
                          opacity=0.85 if hot else 0.5, cls=f"sf-rise sf-d{min((r+c)%6+1,6)}"))
    b.append(text(x0 + 3.5 * cell, y0 + 7 * cell + 28, "H, S  ·  N × N", size=10, fill="var(--sf-muted)"))
    # Sketch.
    x1 = 340
    b.append(rect(x1 - 8, y0 - 8, 3 * cell + 16, 7 * cell + 16, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for r in range(7):
        for c in range(3):
            b.append(rect(x1 + c * cell, y0 + r * cell, cell - 2, cell - 2, rx=2,
                          fill="var(--sf-accent2)" if c == 1 else "var(--sf-soft)", stroke="none", opacity=0.8,
                          cls=f"sf-rise sf-d{min((r+c)%6+1,6)}"))
    b.append(text(x1 + 1.5 * cell, y0 + 7 * cell + 28, "sketch  ·  N × k", size=10, fill="var(--sf-muted)"))
    # Small matrix.
    x2 = 590
    b.append(rect(x2 - 8, y0 + 18, 4 * cell + 16, 4 * cell + 16, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for r in range(4):
        for c in range(4):
            b.append(rect(x2 + c * cell, y0 + 26 + r * cell, cell - 2, cell - 2, rx=2,
                          fill="var(--sf-accent2)" if r == c else "var(--sf-soft)", stroke="none", opacity=0.85))
    b.append(text(x2 + 2 * cell, y0 + 4 * cell + 48, "k × k", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(x0 + 7 * cell + 12, y0 + 3.5 * cell, x1 - 12, y0 + 3.5 * cell, color="var(--sf-accent)", delay=0.2))
    b.append(arrow_between(x1 + 3 * cell + 12, y0 + 3.5 * cell, x2 - 12, y0 + 3.5 * cell, color="var(--sf-accent2)", delay=0.6))
    b.append(pill(248, 210, 74, "project", color="var(--sf-accent)"))
    b.append(pill(490, 210, 74, "solve", color="var(--sf-accent2)"))
    write("experimental-randnla", 780, 280, "".join(b), label="A large Hamiltonian is projected onto a random sketch before solving a small eigenproblem.")


def diagram_experimental_kpm() -> None:
    begin("experimental-kpm")
    b = []
    # Spectrum sticks.
    base_y = 205
    for i, (x, h) in enumerate([(86, 52), (108, 86), (132, 120), (158, 74), (184, 105), (210, 48), (238, 92)]):
        b.append(path(f"M{x},{base_y} L{x},{base_y-h}", stroke="var(--sf-accent)", sw=3, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(path(f"M70,{base_y} L260,{base_y}", stroke="var(--sf-line)", sw=1.2))
    b.append(text(160, 236, "spectrum", size=10, fill="var(--sf-muted)"))
    # Chebyshev waves.
    waves = []
    for i in range(80):
        x = 300 + i * 2.6
        y = 150 + 38 * math.sin(i / 4.0) * math.exp(-i / 90.0)
        waves.append((x, y))
    b.append(curve_path(waves, color="var(--sf-accent2)", sw=2.1, cls="sf-draw"))
    b.append(text(410, 236, "Chebyshev expansion", size=10, fill="var(--sf-muted)"))
    # DOS curve.
    dos = []
    for i in range(90):
        x = 565 + i * 2.2
        y = 205 - 70 * math.exp(-((i - 28) ** 2) / 90) - 95 * math.exp(-((i - 55) ** 2) / 70) - 50 * math.exp(-((i - 78) ** 2) / 55)
        dos.append((x, y))
    b.append(path(f"M565,205 L760,205", stroke="var(--sf-line)", sw=1.1))
    b.append(curve_path(dos, color="var(--sf-warm)", sw=2.6, cls="sf-draw"))
    b.append(text(662, 236, "density of states", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(266, 150, 292, 150, color="var(--sf-accent)", delay=0.2))
    b.append(arrow_between(520, 150, 556, 150, color="var(--sf-accent2)", delay=0.6))
    b.append(pill(330, 36, 76, "moments", color="var(--sf-accent2)"))
    b.append(pill(602, 36, 76, "DOS", color="var(--sf-warm)"))
    write("experimental-kpm", 780, 280, "".join(b), label="A discrete spectrum is expanded in Chebyshev polynomials and turned into a smooth density of states.")


def diagram_experimental_mbh() -> None:
    begin("experimental-mbh")
    b = []
    # Molecule with rigid blocks.
    b.append(benzene(165, 135, 42, rotate=0, cls="sf-spin"))
    b.append(circle(165, 135, 58, fill="none", stroke="var(--sf-accent)", sw=1.2, opacity=0.5, cls="sf-halo"))
    for i, (x, y) in enumerate([(120, 100), (210, 100), (120, 170), (210, 170)]):
        b.append(circle(x, y, 7, fill="var(--sf-accent2)", opacity=0.8, cls=f"sf-pulse sf-d{min(i+1,6)}"))
    b.append(text(165, 220, "rigid blocks", size=10, fill="var(--sf-muted)"))
    # Projection matrix.
    x0, y0, cell = 330, 82, 17
    b.append(rect(x0 - 8, y0 - 8, 7 * cell + 16, 7 * cell + 16, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for r in range(7):
        for c in range(7):
            diag = r == c
            b.append(rect(x0 + c * cell, y0 + r * cell, cell - 2, cell - 2, rx=2,
                          fill="var(--sf-accent)" if diag else "var(--sf-soft)", stroke="none", opacity=0.85 if diag else 0.4))
    b.append(text(x0 + 3.5 * cell, y0 + 7 * cell + 28, "projection", size=10, fill="var(--sf-muted)"))
    # Reduced modes.
    for i, h in enumerate([44, 68, 92, 58, 78, 36]):
        x = 570 + i * 24
        b.append(path(f"M{x},205 L{x},{205-h}", stroke="var(--sf-accent2)", sw=5, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(path("M555,205 L715,205", stroke="var(--sf-line)", sw=1.1))
    b.append(text(635, 232, "reduced modes", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(240, 140, 315, 140, color="var(--sf-accent)", delay=0.2))
    b.append(arrow_between(470, 140, 548, 140, color="var(--sf-accent2)", delay=0.6))
    b.append(pill(360, 36, 80, "3N → 6B", color="var(--sf-accent)"))
    write("experimental-mbh", 780, 280, "".join(b), label="Rigid molecular blocks reduce the full Cartesian Hessian to a small set of vibrational modes.")


def diagram_experimental_cga() -> None:
    begin("experimental-cga")
    b = []
    # Point and conformal embedding.
    b.append(circle(150, 160, 9, fill="var(--sf-accent)", cls="sf-pulse"))
    b.append(text(150, 208, "point", size=10, fill="var(--sf-muted)"))
    # Basis blades.
    for i, (dx, dy, color) in enumerate([(34, -20, "var(--sf-accent)"), (-26, -30, "var(--sf-accent2)"), (8, 34, "var(--sf-violet)")]):
        b.append(path(f"M150,160 L{150+dx},{160+dy}", stroke=color, sw=2, cls="sf-flow"))
        b.append(text(150 + dx * 1.25, 160 + dy * 1.25, ["e₁", "e₂", "e₃"][i], size=9, fill=color))
    b.append(path("M150,108 A52,52 0 1 1 149.9,108", stroke="var(--sf-line)", sw=1, dash="4 6", cls="sf-spin"))
    # Motor rotation.
    b.append(arrow_between(215, 160, 330, 160, color="var(--sf-accent)", delay=0.2))
    b.append(text(272, 145, "motor", size=10, fill="var(--sf-accent)", weight=700))
    # Rotated object (a small molecule) with an arc.
    b.append(benzene(490, 150, 38, rotate=18, cls="sf-spin"))
    b.append(path("M420,110 A90,90 0 0 1 560,110", stroke="var(--sf-accent2)", sw=2, cls="sf-flow"))
    b.append(path("M560,110 l-10,-2 l6,9", stroke="var(--sf-accent2)", sw=2))
    b.append(text(490, 222, "rotated subtree", size=10, fill="var(--sf-muted)"))
    b.append(pill(410, 40, 100, "G(4,1) motor", color="var(--sf-accent2)"))
    write("experimental-cga", 780, 280, "".join(b), label="A conformal motor rotates a molecular subtree around a bond axis.")


def diagram_experimental_eeq() -> None:
    begin("experimental-eeq")
    b = []
    # Atom cloud with coordination shells.
    for i, (x, y, r) in enumerate([(150, 150, 14), (220, 110, 11), (225, 195, 11), (300, 150, 9)]):
        b.append(circle(x, y, r + 8, fill="none", stroke="var(--sf-accent)", sw=1.1, opacity=0.35, cls=f"sf-halo sf-d{min(i+1,6)}"))
        b.append(circle(x, y, r, fill="var(--sf-accent2)" if i % 2 else "var(--sf-accent)", opacity=0.85))
        b.append(text(x, y + r + 20, ["+", "−", "+", "−"][i], size=12, fill="var(--sf-muted)", weight=700))
    b.append(path("M150,150 L220,110 M150,150 L225,195 M220,110 L300,150 M225,195 L300,150", stroke="var(--sf-line)", sw=1.2, cls="sf-flow"))
    b.append(text(225, 245, "coordination", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(340, 150, 430, 150, color="var(--sf-accent)", delay=0.2))
    # Electronegativity curve.
    b.append(axes(455, 78, 250, 130, xlabel="CN", ylabel="χ"))
    pts = [(455 + i * 3.1, 175 - 46 * (1 - math.exp(-i / 18))) for i in range(81)]
    b.append(curve_path(pts, color="var(--sf-warm)", sw=2.5))
    b.append(circle(455 + 52 * 3.1, 175 - 46 * (1 - math.exp(-52 / 18)), 5, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(text(580, 238, "CN-dependent hardness", size=10, fill="var(--sf-muted)"))
    b.append(pill(700, 46, 54, "q", color="var(--sf-accent2)"))
    write("experimental-eeq", 780, 280, "".join(b), label="Coordination numbers shift atomic electronegativity and produce geometry-dependent partial charges.")


def diagram_experimental_alpb() -> None:
    begin("experimental-alpb")
    b = []
    # Solute cavity with ripples.
    for r, op in [(88, 0.18), (68, 0.26), (48, 0.34)]:
        b.append(circle(235, 150, r, fill="none", stroke="var(--sf-accent2)", sw=1.2, opacity=op, cls="sf-halo"))
    b.append(molecule_water(235, 145, 1.0))
    b.append(text(235, 242, "solute cavity", size=10, fill="var(--sf-muted)"))
    # Born radii rays.
    for angle in (20, 90, 160, 250, 320):
        rad = math.radians(angle)
        x1, y1 = 235 + 46 * math.cos(rad), 150 + 46 * math.sin(rad)
        x2, y2 = 235 + 76 * math.cos(rad), 150 + 76 * math.sin(rad)
        b.append(path(f"M{x1:.1f},{y1:.1f} L{x2:.1f},{y2:.1f}", stroke="var(--sf-accent)", sw=1.5, cls="sf-flow"))
    b.append(arrow_between(340, 150, 430, 150, color="var(--sf-accent)", delay=0.2))
    # GB kernel curve.
    b.append(axes(455, 78, 240, 130, xlabel="r", ylabel="f"))
    pts = [(455 + i * 3.0, 195 - 95 / (1 + math.exp(-(i - 35) / 5))) for i in range(80)]
    b.append(curve_path(pts, color="var(--sf-accent)", sw=2.5))
    b.append(text(575, 238, "generalized Born kernel", size=10, fill="var(--sf-muted)"))
    b.append(pill(520, 36, 100, "solvent", color="var(--sf-accent2)"))
    write("experimental-alpb", 780, 280, "".join(b), label="A solute cavity with Born radii is coupled to a continuum solvent through the generalized Born kernel.")


def diagram_experimental_d4() -> None:
    begin("experimental-d4")
    b = []
    # Two fragments with dispersion waves.
    b.append(benzene(190, 160, 40, rotate=8, cls="sf-float"))
    b.append(benzene(560, 140, 34, rotate=-14, cls="sf-float sf-d2"))
    b.append(text(190, 228, "fragment A", size=10, fill="var(--sf-muted)"))
    b.append(text(560, 208, "fragment B", size=10, fill="var(--sf-muted)"))
    for i in range(5):
        y = 110 + i * 18
        b.append(path(f"M255,{y} C330,{y-18} 420,{y+18} 500,{y}", stroke="var(--sf-accent)" if i % 2 == 0 else "var(--sf-accent2)", sw=1.6, opacity=0.65, cls=f"sf-flow sf-d{min(i+1,6)}"))
    b.append(text(378, 244, "dynamic C6 / C8", size=10, fill="var(--sf-muted)"))
    b.append(pill(338, 42, 82, "dispersion", color="var(--sf-warm)"))
    write("experimental-d4", 780, 280, "".join(b), label="Two molecular fragments exchange dynamic dispersion interactions.")


def diagram_experimental_sdr() -> None:
    begin("experimental-sdr")
    b = []
    # Distance constraints.
    for i, (x, y) in enumerate([(150, 95), (110, 180), (195, 180)]):
        b.append(circle(x, y, 8, fill="var(--sf-accent)", opacity=0.85, cls=f"sf-pulse sf-d{min(i+1,6)}"))
    b.append(path("M150,95 L110,180 L195,180 Z", stroke="var(--sf-line)", sw=1.4))
    b.append(text(150, 228, "distance bounds", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(240, 150, 330, 150, color="var(--sf-accent)", delay=0.2))
    # PSD cone.
    b.append(path("M400,210 L470,90 L540,210 Z", fill="var(--sf-soft)", stroke="var(--sf-line)", sw=1.4))
    b.append(path("M400,210 L470,90 L540,210", stroke="var(--sf-accent)", sw=1.8, cls="sf-draw"))
    b.append(circle(470, 150, 7, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(circle(470, 150, 18, fill="none", stroke="var(--sf-accent2)", sw=1, opacity=0.6, cls="sf-halo"))
    b.append(text(470, 238, "positive semidefinite cone", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(560, 150, 640, 150, color="var(--sf-accent2)", delay=0.6))
    # Embedded triangle.
    b.append(path("M675,180 L710,110 L745,180 Z", stroke="var(--sf-accent2)", sw=2.2, cls="sf-draw"))
    for x, y in [(675, 180), (710, 110), (745, 180)]:
        b.append(circle(x, y, 5.5, fill="var(--sf-accent2)"))
    b.append(text(710, 218, "embedding", size=10, fill="var(--sf-muted)"))
    b.append(pill(438, 36, 62, "project", color="var(--sf-accent)"))
    write("experimental-sdr", 780, 280, "".join(b), label="Distance constraints are projected onto the positive-semidefinite cone to obtain a stable embedding.")


def diagram_experimental_cpm() -> None:
    begin("experimental-cpm")
    b = []
    # Electrode slab.
    b.append(rect(80, 70, 180, 120, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for i in range(4):
        b.append(rect(92, 82 + i * 26, 156, 14, rx=4, fill="var(--sf-accent)", opacity=0.25 + i * 0.1))
    b.append(text(170, 215, "electrode", size=10, fill="var(--sf-muted)"))
    # Potential arrow.
    b.append(arrow_between(280, 130, 360, 130, color="var(--sf-warm)", delay=0.2))
    b.append(text(320, 112, "V", size=12, fill="var(--sf-warm)", weight=700))
    # Ions.
    for i, (x, y, sign) in enumerate([(410, 90, "+"), (450, 150, "−"), (410, 200, "+"), (490, 105, "−"), (490, 190, "+")]):
        b.append(circle(x, y, 17, fill="var(--sf-soft)", stroke="var(--sf-line)", sw=1.2, cls=f"sf-float sf-d{min(i+1,6)}"))
        b.append(text(x, y + 5, sign, size=15, fill="var(--sf-accent2)" if sign == "+" else "var(--sf-rose)", weight=700))
    b.append(text(450, 238, "electrolyte", size=10, fill="var(--sf-muted)"))
    # Charge response bars.
    for i, h in enumerate([40, 66, 88, 56]):
        x = 580 + i * 28
        b.append(rect(x, 170 - h, 16, h, rx=4, fill="var(--sf-accent2)" if i % 2 else "var(--sf-accent)", opacity=0.8, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(path("M565,170 L685,170", stroke="var(--sf-line)", sw=1.1))
    b.append(text(625, 196, "induced charge", size=10, fill="var(--sf-muted)"))
    write("experimental-cpm", 780, 280, "".join(b), label="An applied electrode potential induces a charge response in the electrolyte.")


def diagram_experimental_riemannian() -> None:
    begin("experimental-riemannian")
    b = []
    # PSD cone with geodesic.
    b.append(path("M120,220 L250,70 L380,220 Z", fill="var(--sf-soft)", stroke="var(--sf-line)", sw=1.5))
    b.append(path("M120,220 L250,70 L380,220", stroke="var(--sf-accent)", sw=1.8, cls="sf-draw"))
    b.append(text(250, 248, "positive semidefinite manifold", size=10, fill="var(--sf-muted)"))
    # Geodesic path.
    geo = [(150, 200), (210, 140), (250, 116), (300, 142), (345, 192)]
    b.append(curve_path(geo, color="var(--sf-accent2)", sw=2.2))
    b.append(circle(150, 200, 6, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(150, 226, "start", size=9, fill="var(--sf-muted)"))
    # Moving point.
    b.append('<circle r="7" fill="var(--sf-warm)"><animateMotion dur="4s" repeatCount="indefinite" path="M150,200 Q210,140 250,116 Q300,142 345,192"/></circle>')
    b.append(circle(345, 192, 6, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(355, 218, "minimum", size=9, fill="var(--sf-muted)"))
    # Objective contours.
    for r in (26, 44, 62):
        b.append(circle(250, 130, r, fill="none", stroke="var(--sf-violet)", sw=1.1, opacity=0.35 - r * 0.002))
    b.append(text(250, 40, "Riemannian L-BFGS", size=11, fill="var(--sf-violet)", weight=700))
    b.append(pill(430, 70, 90, "retract", color="var(--sf-accent2)"))
    b.append(arrow_between(400, 110, 455, 110, color="var(--sf-accent2)", delay=0.4))
    write("experimental-riemannian", 780, 280, "".join(b), label="An optimizer follows a geodesic on the positive-semidefinite manifold towards a constrained minimum.")


def diagram_validation_pipeline() -> None:
    begin("validation-pipeline")
    b = []
    checks = [
        ("energy", "E/N < 0.05", "var(--sf-accent)"),
        ("tetrahedral", "volume > 0.50", "var(--sf-accent2)"),
        ("chirality", "sign match", "var(--sf-violet)"),
        ("planarity", "SP2 check", "var(--sf-warm)"),
        ("double bond", "not linear", "var(--sf-rose)"),
    ]
    x = 55
    for i, (name, rule, color) in enumerate(checks):
        b.append(rect(x, 78, 124, 102, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(circle(x + 62, 116, 18, fill="var(--sf-soft)", stroke=color, sw=1.6, cls="sf-pulse"))
        b.append(text(x + 62, 121, "✓", size=15, fill=color, weight=700))
        b.append(text(x + 62, 152, name, size=10.5, fill="var(--sf-ink)", weight=650))
        b.append(text(x + 62, 168, rule, size=8.5, fill="var(--sf-muted)", weight=450))
        if i < len(checks) - 1:
            b.append(arrow_between(x + 128, 129, x + 150, 129, color=color, delay=0.2 + i * 0.15))
        x += 150
    b.append(text(390, 40, "accept or retry", size=11, fill="var(--sf-muted)", weight=600))
    b.append(pill(330, 210, 120, "next conformer", color="var(--sf-accent2)"))
    write("validation-pipeline", 790, 260, "".join(b), label="Five geometric checks accept a conformer or send it back for another embedding.")


def diagram_validation_tetrahedral() -> None:
    begin("validation-tetrahedral")
    b = []
    # Valid tetrahedron.
    b.append(path("M100,190 L190,190 L145,115 Z", stroke="var(--sf-line)", sw=1.5))
    b.append(path("M145,115 L155,70 L100,190", stroke="var(--sf-line)", sw=1.5))
    b.append(path("M155,70 L190,190", stroke="var(--sf-line)", sw=1.5, dash="4 4"))
    for x, y in [(100, 190), (190, 190), (145, 115), (155, 70)]:
        b.append(circle(x, y, 7, fill="var(--sf-accent2)", opacity=0.85))
    b.append(circle(145, 142, 8, fill="var(--sf-accent)"))
    b.append(text(145, 230, "volume large", size=10, fill="var(--sf-muted)"))
    b.append(text(145, 48, "accept", size=10, fill="var(--sf-accent2)", weight=700))
    # Invalid collapsed tetrahedron.
    b.append(path("M470,190 L560,190 L515,178 Z", stroke="var(--sf-line)", sw=1.5))
    b.append(path("M515,178 L525,174 L470,190", stroke="var(--sf-line)", sw=1.5))
    for x, y in [(470, 190), (560, 190), (515, 178), (525, 174)]:
        b.append(circle(x, y, 7, fill="var(--sf-rose)", opacity=0.75))
    b.append(circle(515, 183, 8, fill="var(--sf-warm)"))
    b.append(text(515, 230, "volume near zero", size=10, fill="var(--sf-muted)"))
    b.append(text(515, 48, "reject", size=10, fill="var(--sf-rose)", weight=700))
    b.append(arrow_between(250, 130, 400, 130, color="var(--sf-line)", delay=0.3))
    b.append(text(325, 115, "collapse", size=10, fill="var(--sf-muted)"))
    write("validation-tetrahedral", 650, 260, "".join(b), label="A tetrahedral centre with large volume is accepted; a nearly coplanar centre is rejected.")


def diagram_validation_double_bond() -> None:
    begin("validation-double-bond")
    b = []
    # Valid bent geometry.
    b.append(path("M100,160 L220,160", stroke="var(--sf-ink)", sw=3))
    b.append(path("M220,160 L280,110", stroke="var(--sf-ink)", sw=2.4))
    b.append(path("M220,160 L280,210", stroke="var(--sf-ink)", sw=2.4))
    b.append(circle(220, 160, 6, fill="var(--sf-accent)"))
    b.append(text(190, 232, "bent", size=10, fill="var(--sf-muted)"))
    b.append(text(190, 58, "accept", size=10, fill="var(--sf-accent2)", weight=700))
    # Invalid linear geometry.
    b.append(path("M470,160 L590,160", stroke="var(--sf-ink)", sw=3))
    b.append(path("M530,160 L650,160", stroke="var(--sf-ink)", sw=2.4))
    b.append(path("M530,160 L590,160", stroke="var(--sf-rose)", sw=2.4, dash="4 4"))
    b.append(circle(530, 160, 6, fill="var(--sf-rose)"))
    b.append(text(530, 232, "collinear", size=10, fill="var(--sf-muted)"))
    b.append(text(530, 58, "reject", size=10, fill="var(--sf-rose)", weight=700))
    b.append(arrow_between(315, 150, 405, 150, color="var(--sf-line)", delay=0.3))
    write("validation-double-bond", 780, 260, "".join(b), label="Double-bond substituents must stay bent; a collinear arrangement is rejected.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    # Experimental modules
    diagram_experimental_gsm()
    diagram_experimental_randnla()
    diagram_experimental_kpm()
    diagram_experimental_mbh()
    diagram_experimental_cga()
    diagram_experimental_eeq()
    diagram_experimental_alpb()
    diagram_experimental_d4()
    diagram_experimental_sdr()
    diagram_experimental_cpm()
    diagram_experimental_riemannian()
    diagram_validation_pipeline()
    diagram_validation_tetrahedral()
    diagram_validation_double_bond()
    # Conformer pipeline and distance geometry
    diagram_pipeline_overview()
    diagram_eht_pipeline()
    diagram_distance_geometry_pipeline()
    diagram_distance_geometry_matrix()
    diagram_bounds_matrix_overview()
    diagram_bounds_matrix_triangle()
    diagram_bounds_matrix_1_4()
    diagram_embedding_eigendecomp()
    diagram_etkdg_pipeline()
    diagram_etkdg_torsion()
    # Force fields, parsing, charges and spectra
    diagram_force_fields_overview()
    diagram_force_fields_torsion()
    diagram_uff_energy_terms()
    diagram_mmff94_energy_terms()
    diagram_smiles_pipeline()
    diagram_smiles_phenol()
    diagram_smarts_pattern()
    diagram_esp_grid()
    diagram_population_mulliken()
    diagram_population_lowdin()
    diagram_dipole_vector()
    diagram_dos_smearing()
    diagram_kabsch_alignment()
    diagram_materials_unit_cell()
    diagram_materials_assembly()
    diagram_transport_arrow()
    diagram_transport_chunked()
    diagram_uvvis_stda_pipeline()
    diagram_spectroscopy_overview()
    diagram_broadening_functions()
    diagram_ir_hessian_pipeline()
    diagram_nmr_pipeline()
    diagram_nmr_karplus()

    # Backend map, SCF convergence, gradients and benchmarks
    diagram_energy_surface_map()
    diagram_pm3_scf_diis()
    diagram_gradient_central_difference()
    diagram_ks_gradient_units()
    diagram_benchmarks_rmsd()
    diagram_benchmarks_timing()
    diagram_benchmarks_success()
    diagram_ani_architecture()
    diagram_hf3c_architecture()



# ---------------------------------------------------------------------------
# Conformer pipeline and distance geometry
# ---------------------------------------------------------------------------

def matrix_panel(
    x: float,
    y: float,
    rows: int,
    cols: int,
    cell: float,
    *,
    seed: int = 7,
    color: str = "var(--sf-accent)",
    diag: bool = False,
    label: str = "",
    cls: str = "",
) -> str:
    out = rect(x - 7, y - 7, cols * cell + 14, rows * cell + 14, rx=10, fill="var(--sf-card2)", stroke="var(--sf-line)")
    for r in range(rows):
        for c in range(cols):
            hot = (r == c) if diag else ((r * 13 + c * 7 + seed) % 11 < 3)
            out += rect(
                x + c * cell,
                y + r * cell,
                cell - 1.8,
                cell - 1.8,
                rx=1.8,
                fill=color if hot else "var(--sf-soft)",
                stroke="none",
                opacity=0.82 if hot else 0.34,
                cls=cls,
            )
    if label:
        out += text(x + cols * cell / 2, y + rows * cell + 24, label, size=10, fill="var(--sf-muted)")
    return out


def orbital(cx: float, cy: float, r: float = 26, *, color: str = "var(--sf-accent)", rotate: float = 0) -> str:
    out = f'<g transform="rotate({rotate} {cx} {cy})" class="sf-float">'
    out += path(
        f"M{cx},{cy} C{cx-r*0.9},{cy-r*0.9} {cx-r*0.9},{cy+r*0.1} {cx},{cy} "
        f"C{cx+r*0.9},{cy-r*0.1} {cx+r*0.9},{cy+r*0.9} {cx},{cy} Z",
        fill=color,
        stroke="none",
        opacity=0.55,
    )
    out += circle(cx, cy, 4.5, fill="var(--sf-ink)")
    out += "</g>"
    return out


def diagram_pipeline_overview() -> None:
    begin("pipeline-overview")
    b = []
    steps = [
        ("SMILES", "parse", "var(--sf-accent)"),
        ("graph", "atoms + bonds", "var(--sf-accent2)"),
        ("bounds", "lower / upper", "var(--sf-violet)"),
        ("embed", "sample distances", "var(--sf-warm)"),
        ("refine", "force field", "var(--sf-rose)"),
        ("validate", "geometry checks", "var(--sf-accent)"),
    ]
    positions = [(70, 95), (310, 95), (550, 95), (550, 300), (310, 300), (70, 300)]
    for i, ((title, sub, color), (x, y)) in enumerate(zip(steps, positions)):
        b.append(rect(x, y, 180, 110, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, y + 1.5, 177, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 90, y + 48, title, size=15, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 90, y + 70, sub, size=10, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 90, y + 92, 5, fill=color, cls="sf-pulse"))
    # Winding connectors.
    b.append(arrow_between(255, 150, 300, 150, color="var(--sf-accent)", delay=0.2))
    b.append(arrow_between(495, 150, 540, 150, color="var(--sf-accent2)", delay=0.4))
    b.append(arrow_between(640, 210, 640, 290, color="var(--sf-violet)", delay=0.6))
    b.append(arrow_between(540, 355, 495, 355, color="var(--sf-rose)", delay=0.8))
    b.append(arrow_between(300, 355, 255, 355, color="var(--sf-accent)", delay=1.0))
    # A small 3D-conformer glyph in the empty band between the two rows.
    b.append(path("M355,248 L385,228 L420,244 L392,264 Z", stroke="var(--sf-accent2)", sw=1.6, cls="sf-draw"))
    for x, y in [(355, 248), (385, 228), (420, 244), (392, 264)]:
        b.append(circle(x, y, 4.5, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(388, 286, "3D", size=10, fill="var(--sf-muted)"))
    write("pipeline-overview", 800, 470, "".join(b), label="Six-stage conformer pipeline: parse SMILES, build distance bounds, embed, refine, validate and return a 3D conformer.")


def diagram_eht_pipeline() -> None:
    begin("eht-pipeline")
    b = []
    cards = [
        ("basis", "STO-3G", "var(--sf-accent)"),
        ("S and H", "overlap + Hamiltonian", "var(--sf-accent2)"),
        ("orthonormalize", "S⁻¹ᐟ²", "var(--sf-violet)"),
        ("diagonalize", "H′ C = C E", "var(--sf-warm)"),
        ("orbitals", "HOMO / LUMO", "var(--sf-rose)"),
    ]
    x = 35
    for i, (title, sub, color) in enumerate(cards):
        b.append(rect(x, 72, 132, 112, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 73.5, 129, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 66, 114, title, size=11.5, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 66, 132, sub, size=8.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 66, 158, 16, fill="var(--sf-soft)", stroke=color, sw=1.5, cls="sf-pulse"))
        b.append(text(x + 66, 163, str(i + 1), size=12, fill=color, weight=700))
        if i < len(cards) - 1:
            b.append(arrow_between(x + 136, 128, x + 154, 128, color=color, delay=0.2 + i * 0.15))
        x += 154
    b.append(pill(330, 212, 116, "one-shot eigensolve", color="var(--sf-accent2)"))
    write("eht-pipeline", 800, 260, "".join(b), label="Extended Hückel theory: build overlap and Hamiltonian, orthonormalize, diagonalize and read the orbitals.")


def diagram_distance_geometry_pipeline() -> None:
    begin("distance-geometry-pipeline")
    b = []
    cards = [
        ("bounds", "lower / upper", "var(--sf-accent)"),
        ("sample", "random distances", "var(--sf-accent2)"),
        ("embed", "Gram + eigen", "var(--sf-violet)"),
        ("refine", "bounds force field", "var(--sf-warm)"),
    ]
    x = 55
    for i, (title, sub, color) in enumerate(cards):
        b.append(rect(x, 78, 150, 112, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 79.5, 147, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 75, 120, title, size=13, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 75, 140, sub, size=9.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 75, 166, 10, fill="var(--sf-soft)", stroke=color, sw=1.4, cls="sf-pulse"))
        if i < len(cards) - 1:
            b.append(arrow_between(x + 154, 134, x + 176, 134, color=color, delay=0.2 + i * 0.2))
        x += 176
    # Matrix-to-coordinate glyph.
    b.append(matrix_panel(600, 206, 4, 4, 14, seed=3, color="var(--sf-violet)", label=""))
    b.append(arrow_between(665, 260, 710, 260, color="var(--sf-violet)", delay=0.4))
    b.append(path("M725,280 L755,220 L785,280 Z", stroke="var(--sf-accent2)", sw=2, cls="sf-draw"))
    b.append(text(755, 300, "3D", size=10, fill="var(--sf-muted)"))
    write("distance-geometry-pipeline", 800, 320, "".join(b), label="Distance geometry pipeline: bounds, random distance sampling, embedding and force-field refinement.")


def diagram_distance_geometry_matrix() -> None:
    begin("distance-geometry-matrix")
    b = []
    # Distance matrix.
    b.append(matrix_panel(45, 70, 6, 6, 21, seed=1, color="var(--sf-accent)", label="distances"))
    b.append(arrow_between(190, 145, 250, 145, color="var(--sf-accent)", delay=0.2))
    # Gram matrix.
    b.append(matrix_panel(270, 70, 6, 6, 21, seed=5, color="var(--sf-accent2)", diag=True, label="Gram"))
    b.append(arrow_between(415, 145, 475, 145, color="var(--sf-accent2)", delay=0.4))
    # Spectrum.
    b.append(axes(505, 80, 115, 115, xlabel="λ"))
    for i, h in enumerate([70, 55, 42, 30, 20, 12]):
        b.append(path(f"M{520+i*15},195 L{520+i*15},{195-h}", stroke="var(--sf-violet)", sw=4, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(text(562, 232, "eigenvalues", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(632, 145, 675, 145, color="var(--sf-violet)", delay=0.6))
    # 3D coordinates.
    b.append(path("M700,190 L735,95 L775,180 Z", stroke="var(--sf-accent2)", sw=2, cls="sf-draw"))
    for x, y in [(700, 190), (735, 95), (775, 180)]:
        b.append(circle(x, y, 6, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(738, 222, "coordinates", size=10, fill="var(--sf-muted)"))
    b.append(pill(248, 36, 72, "square", color="var(--sf-accent)"))
    write("distance-geometry-matrix", 800, 280, "".join(b), label="A distance matrix becomes a Gram matrix; its leading eigenvalues give 3D coordinates.")


def diagram_bounds_matrix_overview() -> None:
    begin("bounds-matrix-overview")
    b = []
    x0, y0, cell, n = 180, 70, 34, 6
    b.append(rect(x0 - 12, y0 - 12, n * cell + 24, n * cell + 24, rx=14, fill="var(--sf-card2)", stroke="var(--sf-line)"))
    for r in range(n):
        for c in range(n):
            lower = r > c
            upper = r < c
            fill = "var(--sf-accent)" if lower else ("var(--sf-accent2)" if upper else "var(--sf-soft)")
            b.append(rect(x0 + c * cell, y0 + r * cell, cell - 3, cell - 3, rx=4, fill=fill, stroke="none", opacity=0.55 if lower else (0.75 if upper else 0.3)))
    b.append(text(x0 + n * cell / 2, y0 + n * cell + 34, "atom pair (i, j)", size=10, fill="var(--sf-muted)"))
    b.append(text(x0 - 26, y0 + 3 * cell, "i", size=11, fill="var(--sf-muted)"))
    b.append(text(x0 + 3 * cell, y0 - 22, "j", size=11, fill="var(--sf-muted)"))
    # Highlight a pair.
    r, c = 4, 1
    b.append(rect(x0 + c * cell - 2, y0 + r * cell - 2, cell + 1, cell + 1, rx=5, fill="none", stroke="var(--sf-warm)", sw=2.4, cls="sf-pulse"))
    b.append(arrow_between(430, 130, 520, 130, color="var(--sf-warm)", delay=0.2))
    # Bounds interval.
    b.append(text(650, 72, "distance interval", size=10.5, fill="var(--sf-muted)"))
    b.append(path("M560,120 L760,120", stroke="var(--sf-line)", sw=2))
    b.append(path("M600,100 L600,140 M740,100 L740,140", stroke="var(--sf-accent)", sw=2))
    b.append(path("M600,120 L740,120", stroke="var(--sf-accent)", sw=7, opacity=0.35))
    b.append(text(600, 92, "lower", size=9, fill="var(--sf-accent)"))
    b.append(text(740, 92, "upper", size=9, fill="var(--sf-accent)"))
    b.append(text(670, 160, "sampled d", size=9.5, fill="var(--sf-warm)"))
    b.append(circle(670, 120, 6, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(pill(540, 190, 112, "triangle smooth", color="var(--sf-accent2)"))
    write("bounds-matrix-overview", 800, 325, "".join(b), label="The bounds matrix stores a lower and an upper distance for every atom pair.")


def diagram_bounds_matrix_triangle() -> None:
    begin("bounds-matrix-triangle")
    b = []
    a, bpt, c = (120, 220), (320, 90), (520, 220)
    b.append(path(f"M{a[0]},{a[1]} L{bpt[0]},{bpt[1]} L{c[0]},{c[1]} Z", stroke="var(--sf-line)", sw=2))
    for p, name in [(a, "A"), (bpt, "B"), (c, "C")]:
        b.append(circle(*p, 12, fill="var(--sf-accent)", opacity=0.85))
        b.append(text(p[0], p[1] + 4, name, size=10, fill="#fff", weight=700))
    b.append(text(220, 145, "d₁", size=11, fill="var(--sf-accent2)", weight=700))
    b.append(text(420, 145, "d₂", size=11, fill="var(--sf-accent2)", weight=700))
    b.append(text(320, 250, "d₁₃", size=11, fill="var(--sf-warm)", weight=700))
    b.append(path("M140,205 C230,160 410,160 500,205", stroke="var(--sf-warm)", sw=2.2, cls="sf-flow"))
    b.append(pill(265, 34, 110, "triangle inequality", color="var(--sf-accent)"))
    b.append(text(320, 282, "smooth impossible bounds", size=10, fill="var(--sf-muted)"))
    write("bounds-matrix-triangle", 700, 300, "".join(b), label="Triangle smoothing propagates a distance bound through a third atom.")


def diagram_bounds_matrix_1_4() -> None:
    begin("bounds-matrix-1-4")
    b = []
    # Cis conformer.
    b.append(path("M90,190 L140,120 L210,150 L260,80", stroke="var(--sf-ink)", sw=3))
    b.append(path("M90,190 L140,120 L210,150 L260,80", stroke="var(--sf-line)", sw=1, dash="3 5"))
    for x, y in [(90, 190), (140, 120), (210, 150), (260, 80)]:
        b.append(circle(x, y, 7, fill="var(--sf-accent)"))
    b.append(text(175, 235, "cis", size=11, fill="var(--sf-muted)"))
    b.append(text(175, 54, "short 1–4 distance", size=9.5, fill="var(--sf-warm)"))
    # Trans conformer.
    b.append(path("M440,190 L490,120 L560,150 L610,80", stroke="var(--sf-ink)", sw=3))
    for x, y in [(440, 190), (490, 120), (560, 150), (610, 80)]:
        b.append(circle(x, y, 7, fill="var(--sf-accent2)"))
    b.append(text(525, 235, "trans", size=11, fill="var(--sf-muted)"))
    b.append(text(525, 54, "long 1–4 distance", size=9.5, fill="var(--sf-accent2)"))
    b.append(arrow_between(310, 145, 395, 145, color="var(--sf-violet)", delay=0.3))
    b.append(text(352, 128, "rotate", size=10, fill="var(--sf-violet)"))
    b.append(pill(312, 260, 112, "dihedral → bounds", color="var(--sf-violet)"))
    write("bounds-matrix-1-4", 780, 300, "".join(b), label="The 1–4 distance depends on the dihedral angle, so bounds are dihedral-dependent.")


def diagram_embedding_eigendecomp() -> None:
    begin("embedding-eigendecomp")
    b = []
    b.append(matrix_panel(40, 70, 6, 6, 19, seed=8, color="var(--sf-accent)", label="sampled d"))
    b.append(arrow_between(170, 140, 225, 140, color="var(--sf-accent)", delay=0.2))
    b.append(matrix_panel(245, 70, 6, 6, 19, seed=4, color="var(--sf-accent2)", diag=True, label="Gram"))
    b.append(arrow_between(375, 140, 430, 140, color="var(--sf-accent2)", delay=0.4))
    b.append(axes(455, 80, 120, 110, xlabel="λ"))
    for i, h in enumerate([72, 54, 40, 27, 17, 10]):
        b.append(path(f"M{470+i*16},190 L{470+i*16},{190-h}", stroke="var(--sf-violet)", sw=4, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(text(515, 228, "eigenvalues", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(585, 140, 630, 140, color="var(--sf-violet)", delay=0.6))
    b.append(path("M650,190 L690,80 L735,185 Z", stroke="var(--sf-accent2)", sw=2.2, cls="sf-draw"))
    for x, y in [(650, 190), (690, 80), (735, 185)]:
        b.append(circle(x, y, 6, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(692, 222, "3D structure", size=10, fill="var(--sf-muted)"))
    b.append(pill(200, 32, 110, "random sample", color="var(--sf-accent)"))
    write("embedding-eigendecomp", 800, 280, "".join(b), label="Sampled distances form a Gram matrix whose leading eigenvectors place the atoms in 3D.")


def diagram_etkdg_pipeline() -> None:
    begin("etkdg-pipeline")
    b = []
    # Left: distance geometry.
    b.append(rect(40, 80, 210, 160, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(145, 108, "distance geometry", size=12, fill="var(--sf-ink)", weight=700))
    b.append(matrix_panel(82, 126, 4, 4, 16, seed=2, color="var(--sf-accent)"))
    b.append(arrow_between(165, 190, 220, 190, color="var(--sf-accent)", delay=0.2))
    # Middle: torsion histogram.
    b.append(rect(245, 80, 250, 160, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(370, 108, "torsion knowledge", size=12, fill="var(--sf-ink)", weight=700))
    b.append(axes(280, 130, 180, 78))
    for i, h in enumerate([24, 40, 58, 34, 50, 22, 36]):
        b.append(rect(295 + i * 22, 208 - h, 12, h, rx=3, fill="var(--sf-accent2)", opacity=0.65, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(arrow_between(500, 160, 550, 160, color="var(--sf-accent2)", delay=0.4))
    # Right: refined conformer.
    b.append(rect(555, 80, 205, 160, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(657, 108, "force-field refinement", size=11.5, fill="var(--sf-ink)", weight=700))
    b.append(benzene(657, 165, 38, rotate=12, cls="sf-spin"))
    b.append(pill(320, 258, 130, "ETKDG", color="var(--sf-accent)"))
    write("etkdg-pipeline", 800, 300, "".join(b), label="Distance geometry plus CSD torsion preferences plus force-field refinement produces the ETKDG conformer.")


def diagram_etkdg_torsion() -> None:
    begin("etkdg-torsion")
    b = []
    # Bond with rotating substituent.
    b.append(path("M120,160 L220,160", stroke="var(--sf-ink)", sw=3))
    b.append(path("M220,160 L300,100", stroke="var(--sf-accent)", sw=2.6, cls="sf-float"))
    b.append(path("M220,160 L300,220", stroke="var(--sf-accent2)", sw=2.6, cls="sf-float sf-d2"))
    b.append(circle(220, 160, 7, fill="var(--sf-accent)"))
    b.append(path("M240,120 A55,55 0 0 1 300,110", stroke="var(--sf-violet)", sw=1.8, dash="4 5", cls="sf-flow"))
    b.append(text(270, 98, "φ", size=11, fill="var(--sf-violet)", weight=700))
    # Potential curve.
    b.append(axes(390, 75, 280, 140, xlabel="φ", ylabel="E"))
    pts = []
    for i in range(120):
        x = 390 + i * 2.3
        phi = i / 119 * 2 * math.pi
        y = 190 - (38 + 30 * math.cos(3 * phi) + 8 * math.cos(phi) + 6 * math.cos(2 * phi))
        pts.append((x, y))
    b.append(curve_path(pts, color="var(--sf-warm)", sw=2.5))
    # Mark the three minima.
    for phi in (math.pi / 3, math.pi, 5 * math.pi / 3):
        i = int(phi / (2 * math.pi) * 119)
        b.append(circle(390 + i * 2.3, 190 - (38 + 30 * math.cos(3 * phi) + 8 * math.cos(phi) + 6 * math.cos(2 * phi)), 5, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(text(510, 240, "anti", size=10, fill="var(--sf-accent2)"))
    b.append(text(610, 240, "gauche", size=10, fill="var(--sf-accent)"))
    b.append(pill(500, 36, 110, "CSD preference", color="var(--sf-accent2)"))
    write("etkdg-torsion", 700, 280, "".join(b), label="Torsion preferences bias the dihedral angle towards the low-energy anti and gauche wells.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Force fields, parsing, charges and spectra
# ---------------------------------------------------------------------------

def diagram_force_fields_overview() -> None:
    begin("force-fields-overview")
    b = []
    b.append(axes(80, 50, 460, 170, xlabel="distance", ylabel="E"))
    # Allowed window around the minimum.
    b.append(rect(258, 50, 104, 170, rx=8, fill="var(--sf-accent2)", stroke="none", opacity=0.10))
    b.append(text(310, 42, "allowed", size=9.5, fill="var(--sf-accent2)", weight=650))
    # Morse-like potential: repulsive wall, minimum, zero asymptote.
    pts = []
    for i in range(120):
        x = 90 + i * 3.6
        r = 0.8 + i * 0.02
        v = (1 - math.exp(-1.6 * (r - 2.0))) ** 2 - 1
        y = 130 - 52 * v
        pts.append((x, max(55, min(215, y))))
    b.append(curve_path(pts, color="var(--sf-warm)", sw=2.8))
    # Oscillating particle in the well.
    b.append('<circle r="7" fill="var(--sf-accent)"><animateMotion dur="3.6s" repeatCount="indefinite" path="M286,178 Q306,166 326,178 Q306,190 286,178"/></circle>')
    b.append(text(150, 235, "too close", size=9.5, fill="var(--sf-rose)"))
    b.append(text(500, 235, "too far", size=9.5, fill="var(--sf-accent2)"))
    write("force-fields-overview", 600, 280, "".join(b), label="A force field balances repulsion at short range against attraction at long range.")


def diagram_force_fields_torsion() -> None:
    begin("force-fields-torsion")
    b = []
    # Rotating bond glyph.
    b.append(path("M70,160 L160,160", stroke="var(--sf-ink)", sw=3))
    b.append(path("M160,160 L245,105", stroke="var(--sf-accent)", sw=2.6, cls="sf-float"))
    b.append(path("M160,160 L245,215", stroke="var(--sf-accent2)", sw=2.6, cls="sf-float sf-d2"))
    b.append(circle(160, 160, 7, fill="var(--sf-accent)"))
    b.append(path("M180,115 A60,60 0 0 1 250,100", stroke="var(--sf-violet)", sw=1.8, dash="4 5", cls="sf-flow"))
    b.append(text(215, 95, "φ", size=11, fill="var(--sf-violet)", weight=700))
    # Curves.
    b.append(axes(330, 55, 320, 165, xlabel="φ", ylabel="E"))
    v3 = []
    v2 = []
    for i in range(120):
        x = 330 + i * 2.55
        phi = i / 119 * 2 * math.pi
        v3.append((x, 190 - (30 + 28 * math.cos(3 * phi))))
        v2.append((x, 190 - (22 + 22 * math.cos(2 * phi))))
    b.append(curve_path(v3, color="var(--sf-accent)", sw=2.4))
    b.append(curve_path(v2, color="var(--sf-warm)", sw=2.2))
    b.append(text(600, 235, "V3 (sp3)", size=9.5, fill="var(--sf-accent)"))
    b.append(text(455, 235, "V2 (amide)", size=9.5, fill="var(--sf-warm)"))
    b.append(pill(420, 26, 90, "torsion terms", color="var(--sf-violet)"))
    write("force-fields-torsion", 700, 280, "".join(b), label="Torsion potentials with two-fold and three-fold periodicities.")


def diagram_uff_energy_terms() -> None:
    begin("uff-energy-terms")
    b = []
    cards = [
        ("stretch", "spring"),
        ("bend", "angle"),
        ("torsion", "dihedral"),
        ("van der Waals", "contact"),
    ]
    x = 45
    for i, (title, sub) in enumerate(cards):
        b.append(rect(x, 70, 170, 170, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(x + 85, 98, title, size=11, fill="var(--sf-ink)", weight=700))
        if i == 0:
            b.append(path(f"M{x+40},155 Q{x+85},115 {x+130},155", stroke="var(--sf-accent)", sw=2.4, cls="sf-flow"))
            b.append(circle(x + 40, 155, 7, fill="var(--sf-accent)"))
            b.append(circle(x + 130, 155, 7, fill="var(--sf-accent)"))
        elif i == 1:
            b.append(path(f"M{x+40},180 L{x+85},120 L{x+130},180", stroke="var(--sf-accent2)", sw=2.6))
            b.append(circle(x + 85, 120, 8, fill="var(--sf-accent2)", cls="sf-pulse"))
            b.append(path(f"M{x+65},180 A32,32 0 0 1 {x+105},180", stroke="var(--sf-warm)", sw=1.6, dash="3 4"))
        elif i == 2:
            b.append(path(f"M{x+40},180 L{x+75},130 L{x+110},175 L{x+140},125", stroke="var(--sf-violet)", sw=2.4))
            b.append(path(f"M{x+75},130 A28,28 0 0 1 {x+110},128", stroke="var(--sf-violet)", sw=1.6, dash="3 4", cls="sf-flow"))
        else:
            b.append(circle(x + 60, 165, 12, fill="var(--sf-rose)", opacity=0.7))
            b.append(circle(x + 115, 165, 12, fill="var(--sf-accent)", opacity=0.7))
            for r in (24, 34, 44):
                b.append(circle(x + 88, 165, r, fill="none", stroke="var(--sf-line)", sw=1, opacity=0.4))
        b.append(text(x + 85, 218, sub, size=9.5, fill="var(--sf-muted)"))
        x += 185
    write("uff-energy-terms", 800, 270, "".join(b), label="UFF energy terms: bond stretch, angle bend, torsion and van der Waals contact.")


def diagram_mmff94_energy_terms() -> None:
    begin("mmff94-energy-terms")
    b = []
    cards = [
        ("bond", "quartic"),
        ("angle", "cubic"),
        ("stretch-bend", "coupling"),
        ("out-of-plane", "improper"),
        ("nonbonded", "vdW + charge"),
    ]
    x = 30
    for i, (title, sub) in enumerate(cards):
        b.append(rect(x, 64, 148, 166, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(x + 74, 90, title, size=10.5, fill="var(--sf-ink)", weight=700))
        cy = 145
        if i == 0:
            b.append(path(f"M{x+30},{cy+20} L{x+118},{cy+20}", stroke="var(--sf-line)", sw=1.2))
            b.append(path(f"M{x+42},{cy+55} Q{x+74},{cy-30} {x+106},{cy+55}", stroke="var(--sf-accent)", sw=2.4))
            b.append(circle(x + 74, cy - 2, 8, fill="var(--sf-accent)", cls="sf-pulse"))
        elif i == 1:
            b.append(path(f"M{x+30},{cy+42} L{x+74},{cy-28} L{x+118},{cy+42}", stroke="var(--sf-accent2)", sw=2.4))
            b.append(path(f"M{x+54},{cy+42} A28,28 0 0 1 {x+94},{cy+42}", stroke="var(--sf-warm)", sw=1.6, dash="3 4"))
        elif i == 2:
            b.append(path(f"M{x+30},{cy+40} L{x+74},{cy-18} L{x+118},{cy+40}", stroke="var(--sf-violet)", sw=2.4))
            b.append(path(f"M{x+52},{cy+10} L{x+96},{cy+10}", stroke="var(--sf-violet)", sw=1.6, dash="3 4", cls="sf-flow"))
        elif i == 3:
            b.append(path(f"M{x+74},{cy-35} L{x+122},{cy+28} L{x+26},{cy+28} Z", fill="var(--sf-soft)", stroke="var(--sf-line)", sw=1.4))
            b.append(circle(x + 74, cy + 8, 8, fill="var(--sf-warm)", cls="sf-pulse"))
            b.append(path(f"M{x+74},{cy-35} L{x+74},{cy+50}", stroke="var(--sf-warm)", sw=1.4, dash="3 4"))
        else:
            b.append(circle(x + 45, cy + 16, 11, fill="var(--sf-rose)", opacity=0.7))
            b.append(circle(x + 105, cy + 16, 11, fill="var(--sf-accent)", opacity=0.7))
            for r in (22, 32):
                b.append(circle(x + 75, cy + 16, r, fill="none", stroke="var(--sf-line)", sw=1, opacity=0.4))
            b.append(text(x + 75, cy + 56, "+ / −", size=9, fill="var(--sf-muted)"))
        b.append(text(x + 74, 210, sub, size=9, fill="var(--sf-muted)"))
        x += 158
    write("mmff94-energy-terms", 820, 260, "".join(b), label="MMFF94 energy terms: bond, angle, stretch-bend, out-of-plane and nonbonded interactions.")


def diagram_smiles_pipeline() -> None:
    begin("smiles-pipeline")
    b = []
    # Token chips.
    tokens = ["c1", "c", "c", "c", "c", "c1", "O"]
    x = 45
    for i, tok in enumerate(tokens):
        w = 42
        b.append(rect(x, 105, w, 34, rx=8, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(x + w / 2, 127, tok, size=11, fill="var(--sf-accent)", weight=700, mono=True))
        x += 52
        if i < len(tokens) - 1:
            b.append(path(f"M{x-8},122 L{x+2},122", stroke="var(--sf-line)", sw=1.2))
    b.append(text(230, 168, "token stream", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(410, 122, 480, 122, color="var(--sf-accent)", delay=0.2))
    # Graph: benzene ring plus the O substituent on the rightmost carbon.
    b.append(benzene(600, 135, 58, rotate=0, cls="sf-spin"))
    b.append(path("M658,135 L684,135", stroke="var(--sf-line)", sw=2))
    b.append(circle(697, 135, 10, fill="var(--sf-rose)", opacity=0.9, cls="sf-pulse"))
    b.append(text(697, 139, "O", size=10, fill="#fff", weight=700))
    b.append(text(620, 225, "molecular graph", size=10, fill="var(--sf-muted)"))
    b.append(pill(260, 205, 112, "parse + rings", color="var(--sf-accent2)"))
    write("smiles-pipeline", 800, 270, "".join(b), label="A SMILES token stream becomes a molecular graph with atoms, bonds and rings.")


def diagram_smiles_phenol() -> None:
    begin("smiles-phenol")
    b = []
    # Benzene ring with an explicit C-O-H substituent on the rightmost carbon.
    b.append(benzene(300, 155, 88, rotate=0, cls="sf-spin"))
    b.append(path("M388,155 L452,155", stroke="var(--sf-line)", sw=2))
    b.append(circle(467, 155, 15, fill="var(--sf-rose)", opacity=0.9, cls="sf-pulse"))
    b.append(text(467, 160, "O", size=12, fill="#fff", weight=700))
    b.append(path("M477,147 L492,112", stroke="var(--sf-line)", sw=1.8))
    b.append(text(500, 105, "H", size=12, fill="var(--sf-accent)", weight=700))
    b.append(text(300, 285, "phenol", size=13, fill="var(--sf-ink)", weight=700))
    b.append(pill(220, 38, 100, "aromatic ring", color="var(--sf-accent)"))
    b.append(pill(430, 205, 84, "hydroxyl", color="var(--sf-rose)"))
    write("smiles-phenol", 700, 310, "".join(b), label="Phenol: a benzene ring with a hydroxyl substituent.")


def diagram_smarts_pattern() -> None:
    begin("smarts-pattern")
    b = []
    # Query card.
    b.append(rect(40, 70, 220, 170, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(150, 98, "query pattern", size=11, fill="var(--sf-ink)", weight=700))
    b.append(path("M80,170 L120,120 L165,155 L210,105", stroke="var(--sf-violet)", sw=2.2))
    for x, y in [(80, 170), (120, 120), (165, 155), (210, 105)]:
        b.append(circle(x, y, 7, fill="var(--sf-violet)", opacity=0.8))
    b.append(path("M120,120 L165,155", stroke="var(--sf-warm)", sw=5, opacity=0.5, cls="sf-pulse"))
    b.append(text(150, 218, "SMARTS", size=9.5, fill="var(--sf-muted)", mono=True))
    b.append(arrow_between(270, 155, 340, 155, color="var(--sf-violet)", delay=0.2))
    # Molecule scan.
    b.append(rect(360, 70, 400, 170, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(benzene(500, 150, 44, rotate=10, cls="sf-spin"))
    b.append(path("M544,150 L610,150 L650,110", stroke="var(--sf-ink)", sw=2.4))
    b.append(circle(650, 110, 8, fill="var(--sf-rose)"))
    b.append(path("M544,150 L610,150", stroke="var(--sf-warm)", sw=6, opacity=0.55, cls="sf-pulse"))
    b.append(rect(375, 80, 5, 150, rx=2.5, fill="var(--sf-accent2)", stroke="none", opacity=0.8, cls="sf-blink"))
    b.append(text(560, 218, "match", size=10, fill="var(--sf-warm)", weight=700))
    b.append(pill(430, 36, 110, "substructure scan", color="var(--sf-accent2)"))
    write("smarts-pattern", 800, 270, "".join(b), label="A SMARTS pattern scans a molecule and highlights the matching substructure.")


def diagram_esp_grid() -> None:
    begin("esp-grid")
    b = []
    # 3D-ish grid.
    for i in range(9):
        x = 80 + i * 30
        b.append(path(f"M{x},80 L{x},230", stroke="var(--sf-line)", sw=0.8, opacity=0.35))
    for j in range(6):
        y = 80 + j * 30
        b.append(path(f"M80,{y} L320,{y}", stroke="var(--sf-line)", sw=0.8, opacity=0.35))
    # Potential lobes.
    b.append(circle(200, 155, 58, fill="var(--sf-rose)", opacity=0.16, cls="sf-halo"))
    b.append(circle(200, 155, 34, fill="var(--sf-accent)", opacity=0.18, cls="sf-halo sf-d2"))
    b.append(molecule_water(200, 150, 0.9))
    b.append(text(200, 258, "solute in a voxel grid", size=10, fill="var(--sf-muted)"))
    # Scan plane.
    b.append(rect(80, 80, 6, 150, rx=3, fill="var(--sf-warm)", stroke="none", opacity=0.8, cls="sf-blink"))
    b.append(arrow_between(345, 155, 420, 155, color="var(--sf-accent)", delay=0.2))
    # Color legend.
    b.append(text(500, 80, "ESP", size=11, fill="var(--sf-ink)", weight=700))
    b.append(path("M460,110 L620,110", stroke="var(--sf-line)", sw=10, opacity=0.6))
    b.append(path("M460,110 L540,110", stroke="var(--sf-rose)", sw=10, opacity=0.7))
    b.append(path("M540,110 L620,110", stroke="var(--sf-accent)", sw=10, opacity=0.7))
    b.append(text(460, 136, "−", size=14, fill="var(--sf-rose)", weight=700))
    b.append(text(620, 136, "+", size=14, fill="var(--sf-accent)", weight=700))
    b.append(text(540, 180, "electrostatic potential", size=10, fill="var(--sf-muted)"))
    b.append(pill(470, 210, 120, "grid + isosurface", color="var(--sf-accent2)"))
    write("esp-grid", 780, 290, "".join(b), label="The electrostatic potential is sampled on a voxel grid around the solute.")


def diagram_population_mulliken() -> None:
    begin("population-mulliken")
    b = []
    b.append(matrix_panel(45, 70, 5, 5, 22, seed=2, color="var(--sf-accent)", label="density P"))
    b.append(text(100, 210, "P = C · n · Cᵀ", size=9, fill="var(--sf-muted)", mono=True))
    b.append(arrow_between(190, 135, 250, 135, color="var(--sf-accent)", delay=0.2))
    b.append(matrix_panel(275, 70, 5, 5, 22, seed=6, color="var(--sf-accent2)", label="P · S"))
    b.append(arrow_between(420, 135, 480, 135, color="var(--sf-accent2)", delay=0.4))
    # Atom charges.
    for i, (x, y, q, color) in enumerate([(550, 120, "+0.3", "var(--sf-rose)"), (620, 160, "−0.2", "var(--sf-accent)"), (690, 120, "+0.1", "var(--sf-rose)")]):
        b.append(circle(x, y, 22, fill="var(--sf-soft)", stroke=color, sw=1.6, cls=f"sf-pulse sf-d{min(i+1,6)}"))
        b.append(text(x, y + 4, q, size=10, fill=color, weight=700))
    b.append(text(620, 215, "atomic charges", size=10, fill="var(--sf-muted)"))
    b.append(pill(470, 36, 120, "population analysis", color="var(--sf-accent2)"))
    write("population-mulliken", 760, 250, "".join(b), label="Mulliken analysis contracts the density matrix with the overlap matrix into atomic charges.")


def diagram_population_lowdin() -> None:
    begin("population-lowdin")
    b = []
    b.append(matrix_panel(35, 70, 5, 5, 20, seed=4, color="var(--sf-accent)", label="overlap S"))
    b.append(arrow_between(165, 135, 215, 135, color="var(--sf-accent)", delay=0.2))
    b.append(rect(235, 85, 110, 100, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(290, 125, "S⁻¹ᐟ²", size=18, fill="var(--sf-accent2)", weight=700))
    b.append(text(290, 150, "orthogonalizer", size=8.5, fill="var(--sf-muted)"))
    b.append(arrow_between(355, 135, 405, 135, color="var(--sf-accent2)", delay=0.4))
    b.append(matrix_panel(430, 70, 5, 5, 20, seed=8, color="var(--sf-violet)", diag=True, label="P̃"))
    b.append(arrow_between(560, 135, 610, 135, color="var(--sf-violet)", delay=0.6))
    for i, (x, y, q, color) in enumerate([(650, 110, "+0.2", "var(--sf-rose)"), (700, 160, "−0.1", "var(--sf-accent)")]):
        b.append(circle(x, y, 19, fill="var(--sf-soft)", stroke=color, sw=1.5, cls=f"sf-pulse sf-d{min(i+1,6)}"))
        b.append(text(x, y + 4, q, size=9, fill=color, weight=700))
    b.append(text(675, 210, "charges", size=10, fill="var(--sf-muted)"))
    b.append(pill(250, 220, 120, "symmetric orthogonalization", color="var(--sf-accent2)"))
    write("population-lowdin", 760, 250, "".join(b), label="Löwdin analysis orthogonalizes the basis before contracting the density matrix.")


def diagram_dipole_vector() -> None:
    begin("dipole-vector")
    b = []
    b.append(molecule_water(300, 160, 1.5))
    b.append(text(240, 208, "δ−", size=13, fill="var(--sf-rose)", weight=700))
    b.append(text(220, 125, "δ+", size=12, fill="var(--sf-accent)", weight=700))
    b.append(text(380, 125, "δ+", size=12, fill="var(--sf-accent)", weight=700))
    # Dipole arrow.
    b.append(path("M300,160 L470,160", stroke="var(--sf-warm)", sw=4, marker=True, cls="sf-draw"))
    b.append(text(480, 155, "μ", size=16, fill="var(--sf-warm)", weight=700))
    b.append(path("M300,160 A60,60 0 0 1 360,105", stroke="var(--sf-warm)", sw=1.5, dash="4 5", cls="sf-flow"))
    b.append(pill(430, 210, 110, "dipole moment", color="var(--sf-warm)"))
    write("dipole-vector", 700, 260, "".join(b), label="Partial charges on water create a molecular dipole moment.")


def diagram_dos_smearing() -> None:
    begin("dos-smearing")
    b = []
    # Sticks.
    for i, (x, h) in enumerate([(70, 50), (95, 80), (120, 110), (150, 70), (180, 95), (210, 45)]):
        b.append(path(f"M{x},210 L{x},{210-h}", stroke="var(--sf-accent)", sw=3, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(path("M55,210 L230,210", stroke="var(--sf-line)", sw=1.1))
    b.append(text(142, 236, "eigenvalues", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(245, 150, 305, 150, color="var(--sf-accent)", delay=0.2))
    # Gaussian curves.
    for k, (center, height, color) in enumerate([(420, 62, "var(--sf-accent2)"), (470, 88, "var(--sf-warm)"), (520, 55, "var(--sf-violet)")]):
        pts = []
        for i in range(70):
            x = 350 + i * 4
            y = 210 - height * math.exp(-((x - center) ** 2) / (2 * 22 ** 2))
            pts.append((x, y))
        b.append(curve_path(pts, color=color, sw=1.8))
    b.append(text(470, 236, "broadened peaks", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(560, 150, 620, 150, color="var(--sf-accent2)", delay=0.4))
    # Total DOS.
    pts = []
    for i in range(80):
        x = 640 + i * 2
        y = 210 - 40 * math.exp(-((x - 690) ** 2) / 900) - 55 * math.exp(-((x - 725) ** 2) / 700) - 30 * math.exp(-((x - 755) ** 2) / 600)
        pts.append((x, y))
    b.append(path("M630,210 L795,210", stroke="var(--sf-line)", sw=1.1))
    b.append(curve_path(pts, color="var(--sf-accent)", sw=2.6))
    b.append(text(712, 236, "total DOS", size=10, fill="var(--sf-muted)"))
    b.append(pill(380, 36, 100, "Gaussian smearing", color="var(--sf-accent2)"))
    write("dos-smearing", 800, 260, "".join(b), label="Discrete eigenvalues are broadened into a continuous density of states.")


def diagram_kabsch_alignment() -> None:
    begin("kabsch-alignment")
    b = []
    # Before.
    pts_a = [(120, 120), (185, 95), (210, 165), (145, 195)]
    pts_b = [(145, 105), (220, 125), (185, 205), (115, 175)]
    for x, y in pts_a:
        b.append(circle(x, y, 9, fill="var(--sf-accent)", opacity=0.75))
    for x, y in pts_b:
        b.append(circle(x, y, 9, fill="var(--sf-warm)", opacity=0.75))
    b.append(text(165, 232, "before", size=10, fill="var(--sf-muted)"))
    # Rotation arrow.
    b.append(arrow_between(260, 150, 365, 150, color="var(--sf-violet)", curve=-35, delay=0.2))
    b.append(text(310, 112, "SVD", size=10, fill="var(--sf-violet)", weight=700))
    # After (aligned).
    b.append(path("M430,120 L495,95 L520,165 L455,195 Z", stroke="var(--sf-accent2)", sw=1.6, dash="4 5"))
    for x, y in [(430, 120), (495, 95), (520, 165), (455, 195)]:
        b.append(circle(x, y, 9, fill="var(--sf-accent2)", opacity=0.8, cls="sf-pulse"))
    b.append(text(475, 232, "after", size=10, fill="var(--sf-muted)"))
    b.append(text(640, 100, "RMSD", size=11, fill="var(--sf-ink)", weight=700))
    b.append(text(640, 126, "0.12 Å", size=18, fill="var(--sf-accent2)", weight=700))
    b.append(pill(550, 180, 120, "optimal rotation", color="var(--sf-accent2)"))
    write("kabsch-alignment", 760, 270, "".join(b), label="Kabsch alignment finds the optimal rotation between two point clouds.")


def diagram_materials_unit_cell() -> None:
    begin("materials-unit-cell")
    b = []
    # 3D cell.
    b.append(path("M130,200 L260,200 L330,120 L200,120 Z", fill="var(--sf-soft)", stroke="var(--sf-line)", sw=1.6))
    b.append(path("M130,200 L130,95 L200,40 L200,120", stroke="var(--sf-line)", sw=1.6))
    b.append(path("M200,40 L330,40 L330,120", stroke="var(--sf-line)", sw=1.6))
    b.append(path("M260,200 L260,95 L330,40", stroke="var(--sf-line)", sw=1.6))
    # Vectors.
    b.append(arrow_between(200, 160, 290, 160, color="var(--sf-accent)", delay=0.2))
    b.append(text(250, 148, "a", size=11, fill="var(--sf-accent)", weight=700))
    b.append(arrow_between(200, 160, 130, 95, color="var(--sf-accent2)", delay=0.4))
    b.append(text(150, 115, "b", size=11, fill="var(--sf-accent2)", weight=700))
    b.append(arrow_between(200, 160, 200, 55, color="var(--sf-violet)", delay=0.6))
    b.append(text(212, 95, "c", size=11, fill="var(--sf-violet)", weight=700))
    b.append(text(230, 240, "unit cell", size=10, fill="var(--sf-muted)"))
    # Angle arcs.
    b.append(path("M240,160 A40,40 0 0 0 220,140", stroke="var(--sf-warm)", sw=1.5, dash="3 4"))
    b.append(text(255, 130, "γ", size=10, fill="var(--sf-warm)"))
    b.append(pill(430, 60, 90, "a, b, c", color="var(--sf-accent)"))
    b.append(pill(430, 110, 90, "α, β, γ", color="var(--sf-warm)"))
    write("materials-unit-cell", 700, 270, "".join(b), label="A crystallographic unit cell defined by three lattice vectors and three angles.")


def diagram_materials_assembly() -> None:
    begin("materials-assembly")
    b = []
    steps = [
        ("topology", "nodes + edges"),
        ("placement", "secondary units"),
        ("framework", "periodic cell"),
    ]
    x = 55
    for i, (title, sub) in enumerate(steps):
        b.append(rect(x, 70, 200, 180, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(x + 100, 98, title, size=12, fill="var(--sf-ink)", weight=700))
        b.append(text(x + 100, 118, sub, size=9, fill="var(--sf-muted)"))
        if i == 0:
            nodes = [(x + 60, 175), (x + 140, 155), (x + 115, 215), (x + 175, 205)]
            for a, bb in [(0, 1), (0, 2), (1, 3), (2, 3)]:
                b.append(path(f"M{nodes[a][0]},{nodes[a][1]} L{nodes[bb][0]},{nodes[bb][1]}", stroke="var(--sf-line)", sw=1.5, cls="sf-flow"))
            for px, py in nodes:
                b.append(circle(px, py, 7, fill="var(--sf-accent)", cls="sf-pulse"))
        elif i == 1:
            for j, (px, py) in enumerate([(x + 55, 180), (x + 145, 180)]):
                b.append(rect(px - 25, py - 25, 50, 50, rx=8, fill="var(--sf-soft)", stroke="var(--sf-accent2)", sw=1.4, cls=f"sf-rise sf-d{min(j+1,6)}"))
                b.append(circle(px, py, 6, fill="var(--sf-accent2)", cls="sf-pulse"))
            b.append(path(f"M{x+80},180 L{x+120},180", stroke="var(--sf-accent2)", sw=1.8, cls="sf-flow"))
        else:
            b.append(path(f"M{x+45},200 L{x+100},200 L{x+155},200", stroke="var(--sf-accent)", sw=1.4))
            b.append(path(f"M{x+45},200 L{x+75},150 L{x+130},150 L{x+155},200", stroke="var(--sf-accent)", sw=1.4))
            b.append(path(f"M{x+45},200 L{x+75},250 L{x+130},250 L{x+155},200", stroke="var(--sf-accent)", sw=1.4))
            for px, py in [(x + 45, 200), (x + 75, 150), (x + 130, 150), (x + 155, 200), (x + 75, 250), (x + 130, 250)]:
                b.append(circle(px, py, 5, fill="var(--sf-accent2)", opacity=0.8))
        if i < 2:
            b.append(arrow_between(x + 205, 160, x + 240, 160, color="var(--sf-accent2)", delay=0.2 + i * 0.2))
        x += 250
    b.append(pill(300, 270, 160, "periodic framework", color="var(--sf-accent)"))
    write("materials-assembly", 800, 310, "".join(b), label="Topology nodes are decorated with secondary building units to assemble a periodic framework.")


def diagram_transport_arrow() -> None:
    begin("transport-arrow")
    b = []
    # Row-oriented.
    b.append(text(170, 50, "row-oriented", size=11, fill="var(--sf-muted)", weight=650))
    for i in range(3):
        b.append(rect(60, 70 + i * 52, 220, 40, rx=8, fill="var(--sf-card)", stroke="var(--sf-line)"))
        for j, color in enumerate(["var(--sf-accent)", "var(--sf-accent2)", "var(--sf-violet)", "var(--sf-warm)"]):
            b.append(rect(72 + j * 50, 80 + i * 52, 38, 20, rx=4, fill=color, stroke="none", opacity=0.5 + 0.1 * j))
    b.append(arrow_between(300, 145, 380, 145, color="var(--sf-accent)", delay=0.2))
    b.append(text(340, 128, "transpose", size=9.5, fill="var(--sf-muted)"))
    # Column-oriented.
    b.append(text(570, 50, "column-oriented", size=11, fill="var(--sf-muted)", weight=650))
    for j, color in enumerate(["var(--sf-accent)", "var(--sf-accent2)", "var(--sf-violet)", "var(--sf-warm)"]):
        b.append(rect(430 + j * 55, 70, 44, 156, rx=8, fill="var(--sf-card)", stroke="var(--sf-line)"))
        for i in range(3):
            b.append(rect(439 + j * 55, 80 + i * 48, 26, 30, rx=4, fill=color, stroke="none", opacity=0.6, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(pill(460, 250, 150, "contiguous buffers", color="var(--sf-accent2)"))
    write("transport-arrow", 700, 280, "".join(b), label="Row-oriented records are transposed into contiguous columnar buffers.")


def diagram_transport_chunked() -> None:
    begin("transport-chunked")
    b = []
    # Main array.
    b.append(rect(45, 85, 210, 80, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
    for i in range(10):
        b.append(rect(55 + i * 19, 100, 15, 50, rx=3, fill="var(--sf-accent)", opacity=0.55, cls=f"sf-rise sf-d{min(i%6+1,6)}"))
    b.append(text(150, 190, "large result", size=10, fill="var(--sf-muted)"))
    # Chunks.
    for i in range(3):
        y = 60 + i * 66
        b.append(arrow_between(265, 125, 330, y + 24, color="var(--sf-accent)", curve=(-12 if i == 0 else 12 if i == 2 else 0), delay=0.2 + i * 0.2))
        b.append(rect(345, y, 100, 48, rx=8, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(395, y + 28, f"chunk {i}", size=9.5, fill="var(--sf-ink)", weight=650))
    # Workers.
    for i in range(3):
        b.append(rect(520, 55 + i * 66, 120, 52, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(circle(550, 81 + i * 66, 9, fill="var(--sf-accent2)", cls="sf-pulse"))
        b.append(text(590, 86 + i * 66, f"worker {i}", size=9, fill="var(--sf-muted)"))
        b.append(arrow_between(450, 84 + i * 66, 515, 81 + i * 66, color="var(--sf-accent2)", delay=0.4 + i * 0.2))
    b.append(pill(270, 250, 150, "streamed chunks", color="var(--sf-accent2)"))
    write("transport-chunked", 700, 280, "".join(b), label="A large result is streamed in chunks to independent workers.")


def diagram_uvvis_stda_pipeline() -> None:
    begin("uvvis-stda-pipeline")
    b = []
    steps = [
        ("geometry", "3D structure", "var(--sf-accent)"),
        ("orbitals", "EHT / xTB", "var(--sf-accent2)"),
        ("transitions", "CIS-like", "var(--sf-violet)"),
        ("oscillator", "intensity", "var(--sf-warm)"),
        ("spectrum", "ε(λ)", "var(--sf-rose)"),
    ]
    x = 30
    for i, (title, sub, color) in enumerate(steps):
        b.append(rect(x, 70, 132, 120, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 71.5, 129, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 66, 110, title, size=11, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 66, 128, sub, size=8.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 66, 158, 14, fill="var(--sf-soft)", stroke=color, sw=1.5, cls="sf-pulse"))
        if i < len(steps) - 1:
            b.append(arrow_between(x + 136, 130, x + 154, 130, color=color, delay=0.2 + i * 0.15))
        x += 154
    # Spectrum glyph.
    b.append(axes(650, 145, 120, 45))
    pts = []
    for i in range(60):
        xx = 650 + i * 2
        yy = 190 - 34 * math.exp(-((xx - 690) ** 2) / 160) - 24 * math.exp(-((xx - 735) ** 2) / 200)
        pts.append((xx, yy))
    b.append(curve_path(pts, color="var(--sf-rose)", sw=2))
    b.append(pill(320, 220, 120, "sTDA UV-Vis", color="var(--sf-rose)"))
    write("uvvis-stda-pipeline", 800, 250, "".join(b), label="UV-Vis sTDA pipeline from geometry to molecular orbitals, transitions, oscillator strengths and spectrum.")


def diagram_spectroscopy_overview() -> None:
    begin("spectroscopy-overview")
    b = []
    b.append(rect(55, 120, 150, 100, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(130, 160, "3D geometry", size=12, fill="var(--sf-ink)", weight=700))
    b.append(text(130, 180, "elements + coords", size=9, fill="var(--sf-muted)"))
    branches = [
        (320, 65, "UV-Vis", "electronic transitions", "var(--sf-violet)"),
        (320, 155, "IR", "vibrational modes", "var(--sf-warm)"),
        (320, 245, "NMR", "chemical shifts", "var(--sf-accent2)"),
    ]
    for i, (x, y, title, sub, color) in enumerate(branches):
        b.append(arrow_between(210, 170, x - 20, y + 35, color=color, curve=(-18 if i == 0 else 18 if i == 2 else 0), delay=0.2 + i * 0.15))
        b.append(rect(x, y, 200, 70, rx=12, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, y + 1.5, 197, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 100, y + 30, title, size=12, fill="var(--sf-ink)", weight=700))
        b.append(text(x + 100, y + 50, sub, size=9, fill="var(--sf-muted)"))
        b.append(circle(x + 180, y + 35, 6, fill=color, cls="sf-pulse"))
    # Spectrum glyphs.
    b.append(axes(580, 55, 150, 50))
    for i, h in enumerate([26, 38, 20]):
        b.append(path(f"M{610+i*30},105 L{610+i*30},{105-h}", stroke="var(--sf-violet)", sw=3))
    b.append(axes(580, 145, 150, 50))
    for i, h in enumerate([18, 30, 24, 12]):
        b.append(path(f"M{600+i*24},195 L{600+i*24},{195-h}", stroke="var(--sf-warm)", sw=3))
    b.append(axes(580, 235, 150, 50))
    for i, h in enumerate([14, 26, 20, 32]):
        b.append(path(f"M{600+i*24},285 L{600+i*24},{285-h}", stroke="var(--sf-accent2)", sw=3))
    write("spectroscopy-overview", 780, 330, "".join(b), label="One 3D geometry feeds three spectroscopy tracks: UV-Vis, infrared and nuclear magnetic resonance.")


def diagram_broadening_functions() -> None:
    begin("broadening-functions")
    b = []
    b.append(axes(60, 55, 580, 150, xlabel="Δ", ylabel="intensity"))
    # Gaussian.
    g = []
    for i in range(120):
        x = 60 + i * 4.8
        y = 190 - 120 * math.exp(-((x - 350) ** 2) / (2 * 32 ** 2))
        g.append((x, y))
    b.append(curve_path(g, color="var(--sf-accent)", sw=2.6))
    # Lorentzian.
    l = []
    for i in range(120):
        x = 60 + i * 4.8
        y = 190 - 120 / (1 + ((x - 350) / 28) ** 2)
        l.append((x, y))
    b.append(curve_path(l, color="var(--sf-warm)", sw=2.6))
    b.append(text(350, 235, "Gaussian", size=10, fill="var(--sf-accent)"))
    b.append(text(470, 235, "Lorentzian", size=10, fill="var(--sf-warm)"))
    b.append(pill(250, 26, 100, "line shape", color="var(--sf-violet)"))
    write("broadening-functions", 700, 260, "".join(b), label="Gaussian and Lorentzian line shapes differ mainly in their wings.")


def diagram_ir_hessian_pipeline() -> None:
    begin("ir-hessian-pipeline")
    b = []
    steps = [
        ("geometry", "N atoms", "var(--sf-accent)"),
        ("displacements", "6N energies", "var(--sf-accent2)"),
        ("Hessian", "3N × 3N", "var(--sf-violet)"),
        ("modes", "mass-weighted", "var(--sf-warm)"),
        ("spectrum", "IR intensity", "var(--sf-rose)"),
    ]
    x = 30
    for i, (title, sub, color) in enumerate(steps):
        b.append(rect(x, 70, 132, 120, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 71.5, 129, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 66, 112, title, size=11, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 66, 130, sub, size=8.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 66, 160, 14, fill="var(--sf-soft)", stroke=color, sw=1.5, cls="sf-pulse"))
        if i < len(steps) - 1:
            b.append(arrow_between(x + 136, 130, x + 154, 130, color=color, delay=0.2 + i * 0.15))
        x += 154
    # Mode glyph.
    for i, h in enumerate([20, 34, 46, 30]):
        b.append(path(f"M{660+i*22},190 L{660+i*22},{190-h}", stroke="var(--sf-rose)", sw=3, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(pill(330, 220, 140, "numerical Hessian", color="var(--sf-accent2)"))
    write("ir-hessian-pipeline", 800, 250, "".join(b), label="Infrared spectrum from a numerical Hessian: displacements, Hessian, normal modes and intensities.")


def diagram_nmr_pipeline() -> None:
    begin("nmr-pipeline")
    b = []
    steps = [
        ("SMILES", "topology", "var(--sf-accent)"),
        ("HOSE", "environment", "var(--sf-accent2)"),
        ("shifts", "ppm", "var(--sf-violet)"),
        ("coupling", "J(φ)", "var(--sf-warm)"),
        ("spectrum", "peaks", "var(--sf-rose)"),
    ]
    x = 35
    for i, (title, sub, color) in enumerate(steps):
        b.append(rect(x, 70, 132, 112, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 71.5, 129, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 66, 112, title, size=11, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 66, 130, sub, size=8.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 66, 155, 13, fill="var(--sf-soft)", stroke=color, sw=1.5, cls="sf-pulse"))
        if i < len(steps) - 1:
            b.append(arrow_between(x + 136, 126, x + 154, 126, color=color, delay=0.2 + i * 0.15))
        x += 154
    # Spectrum.
    b.append(axes(660, 150, 100, 38))
    for i, h in enumerate([26, 18, 34, 12]):
        b.append(path(f"M{675+i*22},188 L{675+i*22},{188-h}", stroke="var(--sf-rose)", sw=2.5))
    b.append(pill(320, 210, 140, "empirical + Karplus", color="var(--sf-accent2)"))
    write("nmr-pipeline", 800, 240, "".join(b), label="NMR pipeline from SMILES to HOSE environments, chemical shifts, couplings and a broadened spectrum.")


def diagram_nmr_karplus() -> None:
    begin("nmr-karplus")
    b = []
    # Rotating dihedral.
    b.append(path("M70,160 L160,160", stroke="var(--sf-ink)", sw=3))
    b.append(path("M160,160 L235,115", stroke="var(--sf-accent)", sw=2.4, cls="sf-float"))
    b.append(path("M160,160 L235,205", stroke="var(--sf-accent2)", sw=2.4, cls="sf-float sf-d2"))
    b.append(circle(160, 160, 7, fill="var(--sf-accent)"))
    b.append(path("M180,120 A55,55 0 0 1 240,108", stroke="var(--sf-violet)", sw=1.6, dash="4 5", cls="sf-flow"))
    b.append(text(210, 104, "φ", size=10, fill="var(--sf-violet)", weight=700))
    # Karplus curve.
    b.append(axes(330, 55, 320, 160, xlabel="φ", ylabel="³J"))
    pts = []
    for i in range(120):
        x = 330 + i * 2.6
        phi = i / 119 * 2 * math.pi
        y = 185 - (7 + 7 * math.cos(phi) - 5 * math.cos(2 * phi) + 4 * math.cos(3 * phi))
        pts.append((x, y))
    b.append(curve_path(pts, color="var(--sf-accent)", sw=2.5))
    b.append(circle(330 + 2.6 * 50, 185 - (7 + 7 * math.cos(math.pi / 3) - 5 * math.cos(math.pi / 3) + 4 * math.cos(math.pi)), 5, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(text(390, 235, "0°", size=9, fill="var(--sf-muted)"))
    b.append(text(490, 235, "90°", size=9, fill="var(--sf-muted)"))
    b.append(text(590, 235, "180°", size=9, fill="var(--sf-muted)"))
    b.append(pill(430, 26, 100, "Karplus", color="var(--sf-accent2)"))
    write("nmr-karplus", 700, 260, "".join(b), label="The Karplus relation maps a dihedral angle to a vicinal proton-proton coupling constant.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# Backend map, SCF convergence, gradients and benchmarks
# ---------------------------------------------------------------------------

def diagram_energy_surface_map() -> None:
    begin("energy-surface-map")
    b = []
    rows = [
        ("pm3", "MOPAC-parity PM3", "stationary density", "var(--sf-accent)"),
        ("xtb / gfn2", "GFN2-xTB", "warm-started FD", "var(--sf-accent2)"),
        ("gfn1", "GFN1-xTB", "finite differences", "var(--sf-violet)"),
        ("gfn0", "GFN0-xTB", "finite differences", "var(--sf-warm)"),
    ]
    y = 55
    for i, (name, surface, grad, color) in enumerate(rows):
        b.append(rect(35, y, 118, 48, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(94, y + 30, name, size=10.5, fill="var(--sf-ink)", weight=700, mono=True))
        b.append(arrow_between(158, y + 24, 225, y + 24, color=color, delay=0.15 * i))
        b.append(rect(230, y, 205, 48, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(231.5, y + 1.5, 202, 3, rx=1.5, fill=color, stroke="none"))
        b.append(text(332, y + 30, surface, size=10.5, fill="var(--sf-ink)", weight=650))
        b.append(arrow_between(440, y + 24, 510, y + 24, color=color, delay=0.3 + 0.15 * i))
        b.append(rect(515, y, 245, 48, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(text(637, y + 30, grad, size=10, fill="var(--sf-muted)", weight=600))
        y += 62
    b.append(text(94, 34, "public name", size=9.5, fill="var(--sf-muted)", weight=650))
    b.append(text(332, 34, "energy surface", size=9.5, fill="var(--sf-muted)", weight=650))
    b.append(text(637, 34, "gradient path", size=9.5, fill="var(--sf-muted)", weight=650))
    b.append(pill(300, 312, 190, "name == energy == gradient", color="var(--sf-accent2)"))
    write("energy-surface-map", 800, 350, "".join(b), label="Every public method name maps to exactly one energy surface and one gradient path.")


def diagram_pm3_scf_diis() -> None:
    begin("pm3-scf-diis")
    b = []
    b.append(axes(70, 55, 360, 170, xlabel="iteration", ylabel="log residual"))
    # Both curves start at a large residual and fall; the damped loop stalls.
    damp = []
    for i in range(60):
        x = 70 + i * 6
        y = 65 + 68 * (1 - math.exp(-i / 10)) + 4 * math.sin(i / 3) * 0.18
        damp.append((x, y))
    b.append(curve_path(damp, color="var(--sf-rose)", sw=2.2))
    diis = []
    for i in range(60):
        x = 70 + i * 6
        y = 65 + 132 * (1 - math.exp(-i / 6))
        diis.append((x, min(205, y)))
    b.append(curve_path(diis, color="var(--sf-accent2)", sw=2.6))
    b.append(text(285, 150, "damped fixed point", size=9.5, fill="var(--sf-rose)"))
    b.append(text(245, 95, "Pulay DIIS", size=10, fill="var(--sf-accent2)", weight=700))
    # Iteration counters.
    b.append(text(575, 105, "200", size=26, fill="var(--sf-rose)", weight=700))
    b.append(text(575, 126, "budget exhausted", size=9, fill="var(--sf-muted)"))
    b.append(text(575, 178, "27", size=26, fill="var(--sf-accent2)", weight=700))
    b.append(text(575, 199, "converged", size=9, fill="var(--sf-muted)"))
    b.append(pill(520, 232, 160, "same fixed point", color="var(--sf-accent)"))
    write("pm3-scf-diis", 720, 280, "".join(b), label="Pulay DIIS converges the PM3 self-consistent field in far fewer iterations than damped fixed-point mixing.")


def diagram_gradient_central_difference() -> None:
    begin("gradient-central-difference")
    b = []
    # Central atom and displacements.
    b.append(circle(250, 155, 16, fill="var(--sf-rose)", opacity=0.8))
    b.append(text(250, 160, "A", size=11, fill="#fff", weight=700))
    b.append(circle(150, 155, 13, fill="var(--sf-accent)", opacity=0.7, cls="sf-float"))
    b.append(text(150, 160, "−h", size=9, fill="#fff", weight=700))
    b.append(circle(350, 155, 13, fill="var(--sf-accent2)", opacity=0.7, cls="sf-float sf-d2"))
    b.append(text(350, 160, "+h", size=9, fill="#fff", weight=700))
    b.append(path("M170,155 L233,155", stroke="var(--sf-line)", sw=1.6, marker=True))
    b.append(path("M267,155 L330,155", stroke="var(--sf-line)", sw=1.6, marker=True))
    b.append(text(250, 200, "displace one coordinate", size=9.5, fill="var(--sf-muted)"))
    # Energy boxes.
    b.append(rect(60, 230, 180, 52, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(150, 262, "E(R − h)", size=11, fill="var(--sf-accent)", weight=700))
    b.append(rect(260, 230, 180, 52, rx=10, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(350, 262, "E(R + h)", size=11, fill="var(--sf-accent2)", weight=700))
    b.append(arrow_between(150, 210, 150, 226, color="var(--sf-accent)", delay=0.2))
    b.append(arrow_between(350, 210, 350, 226, color="var(--sf-accent2)", delay=0.4))
    # Slope.
    b.append(arrow_between(450, 256, 545, 256, color="var(--sf-warm)", delay=0.6))
    b.append(rect(555, 220, 170, 72, rx=12, fill="var(--sf-card)", stroke="var(--sf-warm)"))
    b.append(text(640, 252, "∂E/∂x", size=16, fill="var(--sf-warm)", weight=700))
    b.append(text(640, 272, "(E₊ − E₋) / 2h", size=9, fill="var(--sf-muted)", mono=True))
    b.append(pill(300, 36, 190, "two energy evaluations", color="var(--sf-accent)"))
    write("gradient-central-difference", 760, 310, "".join(b), label="A central difference displaces one coordinate by plus and minus h and divides the energy difference by 2h.")


def diagram_ks_gradient_units() -> None:
    begin("ks-gradient-units")
    b = []
    b.append(rect(50, 90, 200, 86, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(150, 124, "dE/dR", size=16, fill="var(--sf-accent)", weight=700))
    b.append(text(150, 148, "Hartree / Bohr", size=10, fill="var(--sf-muted)"))
    b.append(arrow_between(255, 133, 335, 133, color="var(--sf-warm)", delay=0.2))
    b.append(text(295, 116, "× 0.5292", size=10, fill="var(--sf-warm)", weight=700))
    b.append(rect(340, 90, 200, 86, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(440, 124, "dE/dR", size=16, fill="var(--sf-accent2)", weight=700))
    b.append(text(440, 148, "Hartree / Å", size=10, fill="var(--sf-muted)"))
    b.append(text(150, 210, "Å per Bohr", size=9.5, fill="var(--sf-muted)"))
    b.append(text(440, 210, "same physical force", size=9.5, fill="var(--sf-muted)"))
    b.append(path("M300,230 L380,230", stroke="var(--sf-line)", sw=2))
    b.append(circle(300, 230, 6, fill="var(--sf-rose)"))
    b.append(circle(380, 230, 6, fill="var(--sf-accent2)"))
    b.append(text(340, 252, "unit conversion, not physics", size=9.5, fill="var(--sf-muted)"))
    b.append(pill(220, 32, 180, "the 0.5292 factor", color="var(--sf-warm)"))
    write("ks-gradient-units", 620, 280, "".join(b), label="A KS-DFT gradient in Hartree per Bohr equals the energy derivative per angstroms multiplied by 0.5292.")


def diagram_benchmarks_rmsd() -> None:
    begin("benchmarks-rmsd")
    b = []
    b.append(text(20, 30, "embed success by category", size=11, fill="var(--sf-ink)", weight=700, anchor="start"))
    cats = [("simple", 0.99), ("alcohol", 0.97), ("amine", 0.98), ("aromatic", 0.95), ("halogen", 0.93)]
    x0, y = 50, 70
    for i, (name, value) in enumerate(cats):
        yy = y + i * 42
        b.append(text(x0 + 66, yy + 18, name, size=9.5, fill="var(--sf-muted)", weight=550))
        b.append(rect(x0 + 90, yy, 420, 22, rx=11, fill="var(--sf-soft)", stroke="none"))
        b.append(rect(x0 + 90, yy, 420 * value, 22, rx=11, fill="var(--sf-accent)", stroke="none", opacity=0.75, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x0 + 525, yy + 16, f"{value*100:.0f}%", size=10, fill="var(--sf-ink)", weight=700))
    b.append(text(300, 290, "success rate", size=9.5, fill="var(--sf-muted)"))
    write("benchmarks-rmsd", 620, 310, "".join(b), label="Conformer embedding succeeds for almost every molecule in each chemical category.")


def diagram_benchmarks_timing() -> None:
    begin("benchmarks-timing")
    b = []
    b.append(axes(60, 55, 500, 170, xlabel="atoms", ylabel="time (ms)"))
    pts = [(70 + i * 38, 200 - (12 + 1.55 * i ** 1.72)) for i in range(12)]
    pts = [(x, max(50, y)) for x, y in pts]
    b.append(curve_path(pts, color="var(--sf-accent)", sw=2.6))
    for x, y in pts[::2]:
        b.append(circle(x, y, 4, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(310, 275, "Floyd–Warshall + BFGS dominate", size=9.5, fill="var(--sf-muted)"))
    b.append(pill(430, 28, 120, "scaling", color="var(--sf-accent)"))
    write("benchmarks-timing", 620, 300, "".join(b), label="Runtime grows with molecule size; triangle smoothing and BFGS dominate.")


def diagram_benchmarks_success() -> None:
    begin("benchmarks-success")
    b = []
    b.append(axes(60, 55, 500, 170, xlabel="RMSD (Å)", ylabel="count"))
    heights = [18, 52, 96, 74, 44, 22, 10]
    for i, h in enumerate(heights):
        x = 85 + i * 62
        b.append(rect(x, 225 - h, 44, h, rx=5, fill="var(--sf-accent)" if i % 2 == 0 else "var(--sf-accent2)", opacity=0.7, cls=f"sf-rise sf-d{min(i+1,6)}"))
    for i, label in enumerate(["0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6"]):
        b.append(text(107 + i * 62, 244, label, size=8.5, fill="var(--sf-muted)"))
    b.append(text(310, 278, "heavy-atom pairwise distribution", size=9.5, fill="var(--sf-muted)"))
    write("benchmarks-success", 620, 300, "".join(b), label="Most conformers have a small heavy-atom RMSD; the distribution has a short right tail.")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------




# ---------------------------------------------------------------------------
# ANI and HF-3c architecture diagrams
# ---------------------------------------------------------------------------

def diagram_ani_architecture() -> None:
    begin("ani-architecture")
    b = []
    steps = [
        ("geometry", "atoms + cell", "var(--sf-accent)"),
        ("neighbors", "cell list", "var(--sf-accent2)"),
        ("AEV", "radial + angular", "var(--sf-violet)"),
        ("network", "per element", "var(--sf-warm)"),
        ("forces", "analytic", "var(--sf-rose)"),
    ]
    x = 30
    for i, (title, sub, color) in enumerate(steps):
        b.append(rect(x, 70, 140, 120, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)"))
        b.append(rect(x + 1.5, 71.5, 137, 4, rx=2, fill=color, stroke="none"))
        b.append(text(x + 70, 112, title, size=11.5, fill="var(--sf-ink)", weight=700, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(x + 70, 130, sub, size=8.5, fill="var(--sf-muted)", weight=450, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(circle(x + 70, 158, 14, fill="var(--sf-soft)", stroke=color, sw=1.5, cls="sf-pulse"))
        if i < len(steps) - 1:
            b.append(arrow_between(x + 144, 130, x + 162, 130, color=color, delay=0.2 + i * 0.15))
        x += 162
    # Small network glyph in the network card.
    for layer, count in enumerate([3, 4, 2]):
        for j in range(count):
            cx = 568 + layer * 24
            cy = 170 + (j - (count - 1) / 2) * 13
            b.append(circle(cx, cy, 3, fill="var(--sf-warm)", opacity=0.8))
    b.append(pill(330, 220, 170, "translation + rotation invariant", color="var(--sf-accent2)"))
    write("ani-architecture", 850, 250, "".join(b), label="ANI inference: geometry, neighbor list, atomic environment vectors, per-element network and analytic forces.")


def diagram_hf3c_architecture() -> None:
    begin("hf3c-architecture")
    b = []
    # Integral and SCF stack.
    b.append(rect(45, 70, 190, 150, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(140, 98, "minimal basis", size=11, fill="var(--sf-ink)", weight=700))
    b.append(text(140, 116, "STO-3G integrals", size=9, fill="var(--sf-muted)"))
    for i in range(3):
        b.append(rect(75, 136 + i * 24, 130, 14, rx=4, fill="var(--sf-accent)", opacity=0.25 + i * 0.1, cls=f"sf-rise sf-d{min(i+1,6)}"))
    b.append(arrow_between(240, 145, 300, 145, color="var(--sf-accent)", delay=0.2))
    b.append(rect(315, 70, 190, 150, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(410, 98, "HF + DIIS", size=11, fill="var(--sf-ink)", weight=700))
    b.append(text(410, 116, "self-consistent field", size=9, fill="var(--sf-muted)"))
    b.append(path("M350,170 Q380,140 410,170 Q440,200 470,170", stroke="var(--sf-accent2)", sw=2, cls="sf-flow"))
    b.append(circle(410, 170, 5, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(arrow_between(505, 145, 565, 145, color="var(--sf-accent2)", delay=0.4))
    # Correction stack.
    b.append(rect(580, 55, 220, 180, rx=16, fill="var(--sf-card)", stroke="var(--sf-line)"))
    b.append(text(690, 82, "empirical corrections", size=11, fill="var(--sf-ink)", weight=700))
    for i, (label, color) in enumerate([("D3(BJ)", "var(--sf-warm)"), ("gCP", "var(--sf-violet)"), ("SRB", "var(--sf-rose)")]):
        b.append(rect(610, 100 + i * 36, 160, 26, rx=8, fill="var(--sf-soft)", stroke=color, sw=1.1, cls=f"sf-rise sf-d{min(i+1,6)}"))
        b.append(text(690, 118 + i * 36, label, size=9.5, fill=color, weight=650))
    b.append(pill(300, 250, 190, "HF-3c energy", color="var(--sf-accent2)"))
    write("hf3c-architecture", 850, 285, "".join(b), label="HF-3c: minimal-basis integrals, self-consistent field with DIIS, and D3, gCP and SRB corrections.")


if __name__ == "__main__":
    main()
