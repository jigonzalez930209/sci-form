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

MB_A = [-200, -100, -170, 15]
MB_a = [-1, -1, -6.5, 0.7]
MB_b = [0, 0, 11, 0.6]
MB_c = [-10, -10, -6.5, 0.7]
MB_x0 = [1, 0, -0.5, -1]
MB_y0 = [0, 0.5, 1.5, 1]


def muller_brown(x: float, y: float) -> float:
    return sum(A * math.exp(a * (x - x0) ** 2 + b_ * (x - x0) * (y - y0) + c * (y - y0) ** 2) for A, a, b_, c, x0, y0 in zip(MB_A, MB_a, MB_b, MB_c, MB_x0, MB_y0))


def mb_grad(x: float, y: float, h: float = 1e-5) -> tuple[float, float]:
    return ((muller_brown(x + h, y) - muller_brown(x - h, y)) / (2 * h), (muller_brown(x, y + h) - muller_brown(x, y - h)) / (2 * h))


def mb_string(n: int = 15, iters: int = 4000) -> list:
    """Simplified string method on Müller-Brown (equal-arc reparametrization)."""
    A, B = (-0.558, 1.442), (0.623, 0.028)
    pts = [(A[0] + (B[0] - A[0]) * i / (n - 1), A[1] + (B[1] - A[1]) * i / (n - 1)) for i in range(n)]
    for _ in range(iters):
        new = [pts[0]]
        for i in range(1, n - 1):
            gx, gy = mb_grad(*pts[i])
            new.append((pts[i][0] - 1e-4 * gx, pts[i][1] - 1e-4 * gy))
        new.append(pts[-1])
        seg = [math.dist(new[i], new[i + 1]) for i in range(n - 1)]
        tot = sum(seg); cum = [0.0]
        for d in seg:
            cum.append(cum[-1] + d)
        out = []
        for i in range(n):
            t = tot * i / (n - 1)
            j = max(k for k in range(n - 1) if cum[k] <= t + 1e-12) if i < n - 1 else n - 2
            f = (t - cum[j]) / (seg[j] or 1)
            out.append((new[j][0] + f * (new[j + 1][0] - new[j][0]), new[j][1] + f * (new[j + 1][1] - new[j][1])))
        pts = out
    return pts


def diagram_experimental_gsm() -> None:
    slug = "experimental-gsm"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Growing string on the Müller–Brown surface (the F82 validation surface)", size=15, weight=700))
    x0_, y0_, w_, h_ = 40, 56, 520, 320
    xr, yr = (-1.5, 1.1), (-0.4, 2.0)
    sx = lambda x: x0_ + (x - xr[0]) / (xr[1] - xr[0]) * w_
    sy = lambda y: y0_ + h_ - (y - yr[0]) / (yr[1] - yr[0]) * h_
    nx, ny = 52, 32
    cw, ch = w_ / nx, h_ / ny
    for j in range(ny):
        for i in range(nx):
            X = xr[0] + (i + 0.5) / nx * (xr[1] - xr[0]); Y = yr[1] - (j + 0.5) / ny * (yr[1] - yr[0])
            v = max(min(muller_brown(X, Y), 100), -150)
            t = (v + 150) / 250
            col = f"rgb({int(40 + 200 * t)},{int(90 + 120 * (1 - abs(t - 0.5) * 2))},{int(235 - 190 * t)})"
            b.append(f'<rect x="{x0_ + i * cw:.1f}" y="{y0_ + j * ch:.1f}" width="{cw + 0.4:.1f}" height="{ch + 0.4:.1f}" fill="{col}" opacity="0.55"/>')
    b.append(rect(x0_, y0_, w_, h_, rx=4, fill="none", sw=1.2))
    pts = mb_string()
    es = [muller_brown(*p_) for p_ in pts]
    ts = max(range(len(pts)), key=lambda i: es[i])
    half = len(pts) // 2
    left, right = pts[: half + 1], pts[half:]
    b.append(path("M" + " L".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in left), stroke="var(--sf-accent2)", sw=2.6))
    b.append(path("M" + " L".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in right), stroke="var(--sf-violet)", sw=2.6))
    for i, (x, y) in enumerate(pts):
        b.append(circle(sx(x), sy(y), 4.5, fill="var(--sf-accent2)" if i <= half else "var(--sf-violet)", stroke="#fff", sw=1))
    b.append(circle(sx(pts[ts][0]), sy(pts[ts][1]), 8, fill="none", stroke="var(--sf-rose)", sw=2.4, cls="sf-pulse"))
    b.append(text(sx(pts[ts][0]) + 12, sy(pts[ts][1]) - 8, "highest node → TS", size=10, fill="var(--sf-rose)", weight=700, anchor="start"))
    for (x, y), lbl in [((-0.558, 1.442), "A −146.7"), ((-0.050, 0.467), "C −80.8"), ((0.623, 0.028), "B −108.2")]:
        b.append(text(sx(x), sy(y) + (22 if lbl[0] != "A" else -14), lbl, size=9.5, fill="var(--sf-ink)", weight=700))
    b.append(text(x0_ + w_ / 2, y0_ + h_ + 18, "x", size=9, fill="var(--sf-muted)"))
    b.append(text(x0_ - 12, y0_ + h_ / 2, "y", size=9, fill="var(--sf-muted)"))
    # right panel: algorithm + energy profile
    x = 584
    b.append(rect(x, 56, 292, 320, rx=14, sw=1.3))
    b.append(text(x + 14, 78, "growth from both ends", size=11, weight=700, anchor="start"))
    b.append(text(x + 22, 98, "τ̂ = (xᵢ₊₁ − xᵢ₋₁)/|…|", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 22, 116, "g⊥ = g − (g·τ̂) τ̂ → relax nodes", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 22, 134, "add node Δs along the path", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 22, 152, "join when tips < Δs/2", size=9.5, weight=500, anchor="start"))
    px, py, pw, ph = x + 34, 180, 236, 130
    emin, emax = -150, max(es) + 10
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    prof = [(px + i / (len(es) - 1) * pw, py + ph - (e - emin) / (emax - emin) * ph) for i, e in enumerate(es)]
    b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in prof), stroke="var(--sf-ink)", sw=2))
    for i, (a, c) in enumerate(prof):
        b.append(circle(a, c, 3.2, fill="var(--sf-accent2)" if i <= half else "var(--sf-violet)"))
    b.append(circle(prof[ts][0], prof[ts][1], 6, fill="none", stroke="var(--sf-rose)", sw=2))
    b.append(text(px + pw / 2, py + ph + 14, "node index along the string", size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(x + 146, 350, f"E_TS = {es[ts]:.1f} · barrier from A = {es[ts] - es[0]:.1f}", size=10, fill="var(--sf-rose)", weight=700))
    b.append(text(x + 146, 366, "validated: barrier 1.6e-3, TS location 0.028", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 400, "".join(b), label="Growing string on the Müller-Brown surface: nodes grown from minima A and B meet, the highest node marks the transition-state region and the energy profile along the string gives the barrier.")


def diagram_experimental_randnla() -> None:
    import numpy as np
    slug = "experimental-randnla"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Randomized Nyström: N × N → sketch N × k → k × k eigenproblem", size=15, weight=700))
    # shapes to scale (N = 600, k = 160)
    N, k = 600, 160
    sc = 0.24
    def block(x, y, w, h, col, lbl, sub):
        b.append(rect(x, y, w * sc, h * sc, rx=4, fill=col, stroke=col, sw=1.2, opacity=0.25))
        b.append(text(x + w * sc / 2, y + h * sc / 2 + 4, lbl, size=11, weight=800))
        b.append(text(x + w * sc / 2, y + h * sc + 16, sub, size=9, fill="var(--sf-muted)", weight=500))
    y = 64
    block(30, y, N, N, "var(--sf-accent)", "S", "N × N (600 × 600)")
    b.append(text(186, y + 76, "·", size=22, weight=800))
    block(200, y, k, N, "var(--sf-warm)", "Ω", "N × k")
    b.append(text(250, y + 76, "=", size=18, weight=800))
    block(266, y, k, N, "var(--sf-accent2)", "Y", "SΩ")
    b.append(arrow_between(330, y + 72, 380, y + 72, color="var(--sf-accent2)"))
    block(390, y, k, N, "var(--sf-violet)", "Q", "Y = QR")
    b.append(arrow_between(434, y + 72, 484, y + 72, color="var(--sf-violet)"))
    block(494, y + 54, k, k, "var(--sf-rose)", "B", "B = QᵀSQ, k × k")
    b.append(text(700, y + 70, "S ≈ Q B Qᵀ", size=13, weight=800, anchor="start"))
    b.append(text(700, y + 92, "k = 50 + ⌈log₂ N⌉ by default", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(700, y + 110, "here k = 160 = N/4 + 10", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(700, y + 136, "validated N = 600, k = 160:", size=9.5, weight=600, anchor="start"))
    b.append(text(700, y + 154, "error 9.6e-4, ~80 vs ~143 ms", size=9.5, fill="var(--sf-accent2)", weight=700, anchor="start"))
    # synthetic spectrum test computed here
    rng = np.random.default_rng(7)
    n2 = 300
    Qm, _ = np.linalg.qr(rng.standard_normal((n2, n2)))
    lam = np.exp(-np.arange(n2) / 18.0) * 10
    S = (Qm * lam) @ Qm.T
    exact = np.sort(np.linalg.eigvalsh(S))[::-1]
    def nystrom(kk, q):
        Om = rng.standard_normal((n2, kk))
        Y = S @ Om
        for _ in range(q):
            Y = S @ (S @ Y) if False else S @ Y
        Qy, _ = np.linalg.qr(Y)
        Bm = Qy.T @ S @ Qy
        return np.sort(np.linalg.eigvalsh(Bm))[::-1]
    a0 = nystrom(40, 0); a2 = nystrom(40, 2)
    x = 24; yb = 262
    b.append(rect(x, yb, 852, 200, rx=14, sw=1.3))
    b.append(text(x + 14, yb + 22, "synthetic check computed for this figure: N = 300, eigenvalues 10·e^(−i/18), k = 40", size=10.5, weight=700, anchor="start"))
    px, py, pw, ph = x + 60, yb + 38, 520, 130
    sxi = lambda i: px + i / 45 * pw
    syl = lambda v: py + ph - (math.log10(max(v, 1e-6)) + 0.3) / 1.7 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    for e_ in (1, 0):
        b.append(text(px - 6, syl(10 ** e_) + 3, f"{10 ** e_:g}", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    for i in range(0, 46, 10):
        b.append(text(sxi(i), py + ph + 13, str(i), size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(path("M" + " L".join(f"{sxi(i):.1f},{syl(v):.1f}" for i, v in enumerate(exact[:45])), stroke="var(--sf-ink)", sw=2.4))
    for arr, col in [(a0, "var(--sf-warm)"), (a2, "var(--sf-accent2)")]:
        for i, v in enumerate(arr[:40]):
            b.append(circle(sxi(i), syl(v), 2.8, fill=col))
    err0 = abs(a0[:4] - exact[:4]).max() / exact[0]; err2 = abs(a2[:4] - exact[:4]).max() / exact[0]
    lx = px + pw + 30
    b.append(path(f"M{lx},{py+10} l20,0", stroke="var(--sf-ink)", sw=2.4)); b.append(text(lx + 26, py + 14, "exact eigenvalues", size=9.5, weight=500, anchor="start"))
    b.append(circle(lx + 10, py + 32, 3, fill="var(--sf-warm)")); b.append(text(lx + 26, py + 36, f"q = 0: top-4 rel. err {err0:.1e}", size=9.5, weight=500, anchor="start"))
    b.append(circle(lx + 10, py + 54, 3, fill="var(--sf-accent2)")); b.append(text(lx + 26, py + 58, f"q = 2: top-4 rel. err {err2:.1e}", size=9.5, weight=500, anchor="start"))
    b.append(text(lx, py + 92, "leading eigenpairs are recovered;", size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(lx, py + 108, "power iterations q sharpen the tail", size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(lx, py + 128, "error > 0.1 % → exact fallback", size=9.5, fill="var(--sf-rose)", weight=600, anchor="start"))
    b.append(text(px + pw / 2, py + ph + 28, "eigenvalue index (y: log scale)", size=8.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 478, "".join(b), label="Randomized Nyström drawn to scale for N 600 and k 160, and a synthetic check where a rank-40 sketch with two power iterations recovers the leading eigenvalues; errors above 0.1 percent trigger the exact fallback.")


def diagram_experimental_kpm() -> None:
    import numpy as np
    slug = "experimental-kpm"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "KPM density of states for a 200-site tight-binding ring (computed for this figure)", size=15, weight=700))
    n = 200
    H = np.zeros((n, n))
    for i in range(n):
        H[i, (i + 1) % n] = H[(i + 1) % n, i] = -1.0
    a_, b0 = 2.05, 0.0
    Ht = (H - b0 * np.eye(n)) / a_
    def moments(M):
        mu = np.zeros(M)
        T0 = np.eye(n); T1 = Ht.copy()
        mu[0] = np.trace(T0) / n; mu[1] = np.trace(T1) / n
        for k in range(2, M):
            T0, T1 = T1, 2 * Ht @ T1 - T0
            mu[k] = np.trace(T1) / n
        return mu
    def jackson(M):
        k = np.arange(M)
        return ((M - k + 1) * np.cos(np.pi * k / (M + 1)) + np.sin(np.pi * k / (M + 1)) / np.tan(np.pi / (M + 1))) / (M + 1)
    def dos(mu, g, xs):
        out = []
        for x in xs:
            Tk = np.cos(np.arange(len(mu)) * np.arccos(x))
            out.append((g[0] * mu[0] + 2 * np.sum(g[1:] * mu[1:] * Tk[1:])) / (np.pi * np.sqrt(1 - x * x)))
        return np.array(out) / a_
    xs = np.linspace(-0.97, 0.97, 300)
    Es = xs * a_
    mu100 = moments(100)
    d20 = dos(mu100[:20], jackson(20), xs); d100 = dos(mu100, jackson(100), xs); dnoj = dos(mu100[:20], np.ones(20), xs)
    ev = np.linalg.eigvalsh(H)
    # left: Chebyshev polynomials
    b.append(rect(24, 50, 260, 330, rx=14, sw=1.3))
    b.append(text(38, 72, "T_{k+1} = 2H̃T_k − T_{k−1}", size=10.5, weight=700, anchor="start"))
    px, py, pw, ph = 48, 92, 220, 150
    for k, col in [(1, "var(--sf-accent)"), (2, "var(--sf-accent2)"), (3, "var(--sf-violet)"), (5, "var(--sf-warm)")]:
        pts = [(px + (x + 1) / 2 * pw, py + ph / 2 - math.cos(k * math.acos(x)) * ph / 2 * 0.95) for x in [-1 + 2 * i / 100 for i in range(101)]]
        b.append(path("M" + " L".join(f"{u:.1f},{v:.1f}" for u, v in pts), stroke=col, sw=1.8))
        b.append(text(px + pw + 2, pts[-1][1] + 3, f"T{k}", size=8.5, fill=col, weight=700, anchor="start"))
    b.append(path(f"M{px},{py+ph/2} L{px+pw},{py+ph/2}", sw=0.8, dash="2 3"))
    b.append(text(154, 266, "H̃ = (H − b)/a ∈ [−1, 1]", size=9.5, weight=600))
    b.append(text(154, 290, "μₖ = Tr[T_k(H̃)]/N", size=10, weight=700))
    b.append(text(154, 310, "1 sparse mat-vec per moment", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(154, 330, "stochastic trace: n_v ±1 probes", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(154, 360, "here: exact trace, a = 2.05", size=9, fill="var(--sf-muted)", weight=450))
    # right: DOS
    x = 304
    b.append(rect(x, 50, 572, 330, rx=14, sw=1.3))
    px, py, pw, ph = x + 50, 80, 480, 230
    ymax = 0.6
    sx = lambda e: px + (e + 2.1) / 4.2 * pw
    sy = lambda v: py + ph - min(v, ymax) / ymax * ph
    hist, edges = np.histogram(ev, bins=42, range=(-2.1, 2.1))
    bw = edges[1] - edges[0]
    for hv, e0 in zip(hist, edges[:-1]):
        hh = hv / n / bw
        b.append(rect(sx(e0) + 0.5, sy(hh), sx(e0 + bw) - sx(e0) - 1, py + ph - sy(hh), rx=1, fill="var(--sf-line)", stroke="none", opacity=0.45))
    for arr, col, sw_, dash in [(dnoj, "var(--sf-rose)", 1.2, "3 3"), (d20, "var(--sf-warm)", 1.8, None), (d100, "var(--sf-accent)", 2.4, None)]:
        b.append(path("M" + " L".join(f"{sx(E):.1f},{sy(max(v, 0)):.1f}" for E, v in zip(Es, arr)), stroke=col, sw=sw_, dash=dash, cls="sf-draw"))
    b.append(path(f"M{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    for e in (-2, -1, 0, 1, 2):
        b.append(text(sx(e), py + ph + 13, str(e), size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(px + pw, py + ph + 26, "E (hopping units)", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    lg = [("exact eigenvalue histogram", "var(--sf-line)", None), ("M = 20, no kernel (Gibbs ringing)", "var(--sf-rose)", "3 3"), ("M = 20, Jackson", "var(--sf-warm)", None), ("M = 100, Jackson", "var(--sf-accent)", None)]
    for i, (t_, col, dash) in enumerate(lg):
        yy = py + 104 + i * 18
        b.append(path(f"M{sx(-0.85)},{yy} l18,0", stroke=col, sw=2.4 if i else 6, dash=dash))
        b.append(text(sx(-0.85) + 24, yy + 3, t_, size=9, weight=500, anchor="start"))
    b.append(text(x + 286, 364, "van Hove peaks at band edges ±2 sharpen as M grows; validated moments 3.5e-15, DOS L2 1.7e-14", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="Kernel polynomial method on a 200-site ring: Chebyshev polynomials, and the density of states from 20 and 100 moments with the Jackson kernel compared with the exact eigenvalue histogram and an undamped expansion that rings.")


WATER_EEQ = [-0.547, 0.276, 0.271]
WATER_EEQ_CN = [1.476, 0.757, 0.719]


def benzene_atoms(cx: float, cy: float, r: float) -> tuple[list, list]:
    C = [(cx + r * math.cos(math.radians(90 + 60 * i)), cy - r * math.sin(math.radians(90 + 60 * i))) for i in range(6)]
    H = [(cx + 1.62 * r * math.cos(math.radians(90 + 60 * i)), cy - 1.62 * r * math.sin(math.radians(90 + 60 * i))) for i in range(6)]
    return C, H


def diagram_experimental_mbh() -> None:
    slug = "experimental-mbh"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Mobile Block Hessian on benzene: 36 Cartesian DOF → 24", size=15, weight=700))
    # left: benzene with rigid ring block and flexible H
    b.append(rect(24, 50, 300, 330, rx=14, sw=1.3))
    cx, cy = 174, 200
    C, H = benzene_atoms(cx, cy, 52)
    b.append(circle(cx, cy, 68, fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1.4, opacity=0.13))
    for i in range(6):
        b.append(path(f"M{C[i][0]:.1f},{C[i][1]:.1f} L{C[(i+1)%6][0]:.1f},{C[(i+1)%6][1]:.1f}", stroke="var(--sf-ink)", sw=2.4))
        b.append(path(f"M{C[i][0]:.1f},{C[i][1]:.1f} L{H[i][0]:.1f},{H[i][1]:.1f}", sw=1.6))
    for (hx, hy) in H:
        for dx, dy, col in [(14, 0, "var(--sf-warm)"), (0, -14, "var(--sf-warm)"), (-9, 9, "var(--sf-warm)")]:
            b.append(path(f"M{hx:.1f},{hy:.1f} l{dx},{dy}", stroke=col, sw=1.3))
        b.append(atom_ball(slug, hx, hy, "H", 6))
    for (x_, y_) in C:
        b.append(atom_ball(slug, x_, y_, "C", 9, label=False))
    b.append(path(f"M{cx-30},{cy} A30,30 0 1 1 {cx},{cy+30}", stroke="var(--sf-accent)", sw=1.8, cls="sf-flow"))
    b.append(text(174, 74, "ring = one rigid block", size=11, fill="var(--sf-accent)", weight=700))
    b.append(text(174, 330, "block: 3 translations + 3 rotations", size=9.5, weight=600))
    b.append(text(174, 348, "each H: x, y, z (flexible)", size=9.5, fill="var(--sf-warm)", weight=600))
    # middle: L matrix drawn to scale 36 x 24
    x = 344
    b.append(rect(x, 50, 260, 330, rx=14, sw=1.3))
    b.append(text(x + 130, 74, "projection L (36 × 24)", size=11, weight=700))
    cs = 6.2
    mx, my = x + 56, 92
    for r_ in range(36):
        for c_ in range(24):
            atom = r_ // 3
            if atom < 6:
                hot = c_ < 6
                col = "var(--sf-accent)"
            else:
                hot = c_ == 6 + (atom - 6) * 3 + (r_ % 3)
                col = "var(--sf-warm)"
            if hot:
                b.append(f'<rect x="{mx + c_ * cs:.1f}" y="{my + r_ * cs * 0.62:.1f}" width="{cs - 0.8:.1f}" height="{cs * 0.62 - 0.6:.1f}" fill="{col}" opacity="0.8"/>')
    b.append(rect(mx - 1, my - 1, 24 * cs + 2, 36 * cs * 0.62 + 2, rx=2, fill="none", sw=0.8))
    b.append(text(mx - 6, my + 18 * cs * 0.62 * 0.33, "C", size=9, fill="var(--sf-accent)", weight=700, anchor="end"))
    b.append(text(mx - 6, my + 18 * cs * 0.62 * 1.4, "H", size=9, fill="var(--sf-warm)", weight=700, anchor="end"))
    b.append(text(mx + 3 * cs, my - 6, "block 6", size=8.5, fill="var(--sf-accent)", weight=700))
    b.append(text(mx + 15 * cs, my - 6, "flexible 18", size=8.5, fill="var(--sf-warm)", weight=700))
    b.append(text(x + 130, 254, "translation + rotation columns", size=9.5, weight=500))
    b.append(text(x + 130, 272, "about the block COM", size=9.5, weight=500))
    b.append(text(x + 130, 300, "H_MBH = LᵀHL,  M_MBH = LᵀML", size=10, weight=700))
    b.append(text(x + 130, 322, "H_MBH c = λ M_MBH c", size=10, weight=700))
    b.append(text(x + 130, 346, "ν = √λ · 108.59 / 2πc", size=9.5, fill="var(--sf-muted)", weight=500))
    # right: DOF bars
    x = 624
    b.append(rect(x, 50, 252, 330, rx=14, sw=1.3))
    b.append(text(x + 126, 74, "n_DOF = 6·n_blocks + 3·n_flex", size=10.5, weight=700))
    for i, (lbl, v, col) in enumerate([("full Hessian", 36, "var(--sf-line)"), ("MBH", 24, "var(--sf-accent2)")]):
        bx = x + 50 + i * 90
        h = v / 36 * 180
        b.append(rect(bx, 290 - h, 60, h, rx=6, fill=col, stroke="none", opacity=0.8))
        b.append(text(bx + 30, 284 - h, str(v), size=14, weight=800))
        b.append(text(bx + 30, 306, lbl, size=9.5, weight=600))
    b.append(text(x + 126, 334, "= 6·1 + 3·6 = 24  → matrix 44 % of the cells", size=9.5, weight=500))
    b.append(text(x + 126, 356, "validated: zero modes ≤ 5 cm⁻¹, stretch 1.3e-5", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="Mobile Block Hessian on benzene: the six-carbon ring is one rigid block with six collective coordinates and the six hydrogens keep three each, so the 36 by 24 projection reduces the Hessian from 36 to 24 degrees of freedom.")


def diagram_experimental_cga() -> None:
    slug = "experimental-cga"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Conformal geometric algebra G(4,1): points, motors and the sandwich product", size=15, weight=700))
    # left: grade structure
    b.append(rect(24, 50, 270, 330, rx=14, sw=1.3))
    b.append(text(159, 74, "32 components by grade", size=11, weight=700))
    grades = [(0, 1, "scalar"), (1, 5, "vectors"), (2, 10, "bivectors"), (3, 10, "trivectors"), (4, 5, "4-vectors"), (5, 1, "pseudoscalar")]
    for i, (g, n, nm) in enumerate(grades):
        y = 96 + i * 40
        b.append(text(44, y + 15, f"{g}", size=11, fill="var(--sf-violet)", weight=800))
        for k in range(n):
            b.append(rect(60 + k * 14, y + 4, 11, 16, rx=2, fill="var(--sf-violet)" if g in (0, 2, 4) else "var(--sf-line)", stroke="none", opacity=0.75))
        b.append(text(278, y + 15, f"{nm} · {n}", size=9, weight=500, anchor="end"))
    b.append(text(159, 352, "motors live in the even grades (0, 2, 4)", size=9.5, fill="var(--sf-violet)", weight=600))
    b.append(text(159, 368, "e₁² = e₂² = e₃² = e₊² = 1, e₋² = −1", size=9, fill="var(--sf-muted)", weight=450))
    # middle: point embedding with a real point
    x = 314
    p_ = (1.20, 0.50, 0.00)
    sq = sum(v * v for v in p_)
    b.append(rect(x, 50, 260, 330, rx=14, sw=1.3))
    b.append(text(x + 130, 74, "lift a 3D point", size=11, weight=700))
    b.append(text(x + 130, 98, "P = p + ½|p|² e∞ + eₒ", size=11, weight=700))
    b.append(text(x + 130, 120, f"p = ({p_[0]:.2f}, {p_[1]:.2f}, {p_[2]:.2f}) Å", size=10, weight=500))
    comps = [("e₁", p_[0]), ("e₂", p_[1]), ("e₃", p_[2]), ("e∞", 0.5 * sq), ("eₒ", 1.0)]
    for i, (nm, v) in enumerate(comps):
        y = 140 + i * 30
        b.append(rect(x + 40, y, 180, 24, rx=6, fill="var(--sf-card2)", sw=0.8))
        b.append(text(x + 56, y + 16, nm, size=10.5, fill="var(--sf-violet)", weight=700, anchor="start"))
        b.append(text(x + 208, y + 16, f"{v:+.4f}", size=10.5, weight=600, anchor="end", mono=True))
    b.append(text(x + 130, 312, "P·P = 0 (null vector) — distances", size=9.5, weight=500))
    b.append(text(x + 130, 330, "become inner products: P·Q = −½|p − q|²", size=9.5, weight=500))
    # right: motor rotating a substituent about a bond (Rodrigues-equivalent numbers)
    x = 594
    b.append(rect(x, 50, 282, 330, rx=14, sw=1.3))
    b.append(text(x + 141, 74, "motor M = T R about the C–C axis", size=11, weight=700))
    cx, cy = x + 141, 160
    # side view of the C–C axis
    ax0, ax1 = cx - 110, cx - 10
    b.append(path(f"M{ax0},{cy} L{ax1},{cy}", stroke="var(--sf-ink)", sw=3))
    b.append(path(f"M{ax0-12},{cy} L{ax1+18},{cy}", stroke="var(--sf-violet)", sw=1, dash="2 3"))
    b.append(atom_ball(slug, ax0, cy, "C", 11)); b.append(atom_ball(slug, ax1, cy, "C", 11))
    b.append(text((ax0 + ax1) / 2, cy + 28, "axis L̂", size=9, fill="var(--sf-violet)", weight=600))
    # end-on view along the axis: H before (dashed) and after a 60° motor
    ex, ey, R = cx + 76, cy, 40
    b.append(circle(ex, ey, R, fill="none", stroke="var(--sf-line)", sw=1, opacity=0.7))
    b.append(circle(ex, ey, 4, fill="var(--sf-ink)"))
    for ang, solid in [(90, False), (30, True)]:
        hx, hy = ex + R * math.cos(math.radians(ang)), ey - R * math.sin(math.radians(ang))
        b.append(path(f"M{ex},{ey} L{hx:.1f},{hy:.1f}", stroke="var(--sf-violet)" if solid else "var(--sf-line)", sw=2, dash=None if solid else "4 3"))
        b.append(atom_ball(slug, hx, hy, "H", 6) if solid else circle(hx, hy, 6, fill="none", stroke="var(--sf-line)", sw=1.2))
    b.append(path(f"M{ex + 26 * math.cos(math.radians(90)):.1f},{ey - 26:.1f} A26,26 0 0 1 {ex + 26 * math.cos(math.radians(30)):.1f},{ey - 26 * math.sin(math.radians(30)):.1f}", stroke="var(--sf-violet)", sw=1.6, cls="sf-flow"))
    b.append(text(ex + 26, ey - 30, "60°", size=10, fill="var(--sf-violet)", weight=700, anchor="start"))
    b.append(text(ex, ey + R + 16, "view along the axis", size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(cx, cy + 74, "X′ = M X M̃  (sandwich)", size=10.5, weight=700))
    b.append(text(cx, cy + 92, "R = cos(θ/2) − sin(θ/2) L̂", size=9.5, weight=500))
    b.append(text(cx, cy + 108, "T = 1 − ½ t e∞", size=9.5, weight=500))
    b.append(text(cx, cy + 132, "vs Rodrigues + clifford Cl(4,1):", size=9.5, weight=600))
    b.append(text(cx, cy + 150, "motors 3.4e-15, dihedral 6.6e-15", size=9.5, fill="var(--sf-accent2)", weight=700))
    b.append(text(cx, cy + 166, "16 rigid motions", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="Conformal geometric algebra: 32 multivector components by grade, the lift of the point (1.2, 0.5, 0) into a null vector, and a motor rotating a hydrogen 60 degrees about a C–C axis through the sandwich product.")


def diagram_experimental_eeq() -> None:
    slug = "experimental-eeq"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "EEQ charges for water — real sci-form output (sci-form eeq)", size=15, weight=700))
    # left: CN Fermi function
    energy_panel(b, 24, 50, 300, 330, "coordination number term", "1 / (1 + exp(−16 (r_cov/r − 1)))",
                 [(lambda r: 1 / (1 + math.exp(-16 * (1 / r - 1))), "var(--sf-accent)", 2.4, None)], (0.5, 2.0), (0, 1.05),
                 xlabel="r / r_cov", marks=[(1.0, "r = r_cov → ½", "var(--sf-accent)")],
                 note="smooth, differentiable bond count")
    # middle: linear system
    x = 344
    b.append(rect(x, 50, 260, 330, rx=14, sw=1.3))
    b.append(text(x + 130, 74, "(N + 1) × (N + 1) system", size=11, weight=700))
    b.append(text(x + 130, 96, "μ = χᵢ + ηᵢqᵢ + Σⱼ γ(rᵢⱼ) qⱼ", size=10, weight=600))
    cs = 34; mx, my = x + 34, 116
    labs = ["O", "H", "H", "λ"]
    for i in range(4):
        for j in range(4):
            if i == 3 and j == 3:
                fill, op, t_ = "var(--sf-card2)", None, "0"
            elif i == 3 or j == 3:
                fill, op, t_ = "var(--sf-warm)", 0.35, "1"
            elif i == j:
                fill, op, t_ = "var(--sf-violet)", 0.6, "η"
            else:
                fill, op, t_ = "var(--sf-accent)", 0.3, "γ"
            b.append(rect(mx + j * cs, my + i * cs, cs - 3, cs - 3, rx=4, fill=fill, stroke="none" if op else "var(--sf-line)", sw=0.6, opacity=op))
            b.append(text(mx + j * cs + cs / 2 - 1.5, my + i * cs + cs / 2 + 3, t_, size=10, weight=700))
        b.append(text(mx - 10, my + i * cs + cs / 2 + 3, labs[i], size=9, fill="var(--sf-muted)"))
    b.append(text(mx + 4 * cs + 14, my + 2 * cs, "·", size=18, weight=800))
    vec = ["q_O", "q_H", "q_H", "μ"]
    for i, v in enumerate(vec):
        b.append(rect(mx + 4 * cs + 26, my + i * cs, 44, cs - 3, rx=4, fill="var(--sf-card2)", sw=0.6))
        b.append(text(mx + 4 * cs + 48, my + i * cs + cs / 2 + 3, v, size=9.5, weight=700))
    b.append(text(x + 130, 274, "γ = erf(√2 r / σᵢⱼ) / r", size=10, weight=600))
    b.append(text(x + 130, 296, "last row: Σ qᵢ = Q_total", size=9.5, fill="var(--sf-warm)", weight=600))
    b.append(text(x + 130, 318, "χᵢ shifted by CNᵢ (geometry)", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 130, 348, "one dense solve, no SCF", size=9.5, fill="var(--sf-muted)", weight=450))
    # right: water with real charges and CN
    x = 624
    b.append(rect(x, 50, 252, 330, rx=14, sw=1.3))
    b.append(text(x + 126, 74, "water (embed seed 42)", size=11, weight=700))
    b.append(water_2d(slug, x + 126, 172, 70, charges=WATER_EEQ, r=16))
    for i, (lb, q, cn) in enumerate(zip(["O", "H₁", "H₂"], WATER_EEQ, WATER_EEQ_CN)):
        y = 260 + i * 20
        b.append(text(x + 30, y, lb, size=10, weight=700, anchor="start"))
        b.append(text(x + 130, y, f"q = {q:+.3f}", size=10, weight=500, anchor="end", mono=True))
        b.append(text(x + 230, y, f"CN = {cn:.3f}", size=10, weight=500, anchor="end", mono=True))
    b.append(text(x + 126, 338, "Σq = 1e−16 · vs multicharge 0.6.0:", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 126, 354, "charges agree to 6.9e−12 e", size=9, fill="var(--sf-accent2)", weight=700))
    write(slug, 900, 396, "".join(b), label="EEQ for water: the Fermi-function coordination number, the bordered linear system with the charge-neutrality row, and sci-form's charges O −0.547 and H +0.276 and +0.271 with coordination numbers 1.48, 0.76 and 0.72.")


def diagram_experimental_alpb() -> None:
    slug = "experimental-alpb"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Generalized Born + ALPB correction (formulas as implemented)", size=15, weight=700))
    Ri, Rj = 1.6, 2.2
    fgb = lambda r: math.sqrt(r * r + Ri * Rj * math.exp(-r * r / (4 * Ri * Rj)))
    energy_panel(b, 24, 50, 300, 330, "GB kernel f_GB (R_i = 1.6, R_j = 2.2 Å)", "√(r² + RᵢRⱼ exp(−r²/4RᵢRⱼ))",
                 [(fgb, "var(--sf-accent)", 2.4, None), (lambda r: r, "var(--sf-line)", 1.4, "4 3")], (0, 8), (0, 8.5), xlabel="r (Å)",
                 marks=[], note="r = 0: √(RᵢRⱼ) = 1.88 Å (self term) · large r: → r")
    A = 0.571412
    x = 344
    b.append(rect(x, 50, 260, 330, rx=14, sw=1.3))
    b.append(text(x + 130, 74, "ALPB scaling of E_GB", size=11, weight=700))
    b.append(text(x + 130, 96, "E_ALPB = E_GB (ε − 1)/(ε + A·x)", size=10, weight=600))
    b.append(text(x + 130, 114, "A = 0.571412 · x = E_GB / E_Coulomb", size=9, fill="var(--sf-muted)", weight=500))
    px, py, pw, ph = x + 40, 134, 200, 150
    sx = lambda e: px + math.log10(e) / math.log10(100) * pw
    sy = lambda v: py + ph - v / 1.0 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    for xx, col in [(0.0, "var(--sf-line)"), (0.5, "var(--sf-accent2)"), (1.0, "var(--sf-violet)")]:
        pts = [(sx(1 + i * 0.99), sy((1 + i * 0.99 - 1) / (1 + i * 0.99 + A * xx) / 1.0)) for i in range(100)]
        b.append(path("M" + " L".join(f"{u:.1f},{v:.1f}" for u, v in pts), stroke=col, sw=2))
    for e in (1, 10, 78.5):
        b.append(text(sx(e), py + ph + 13, f"{e:g}", size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(path(f"M{sx(78.5):.1f},{py} L{sx(78.5):.1f},{py+ph}", stroke="var(--sf-accent)", sw=1, dash="3 3"))
    b.append(text(sx(78.5) - 4, py + 12, "water", size=8.5, fill="var(--sf-accent)", weight=700, anchor="end"))
    b.append(text(px + pw, py + ph + 26, "solvent ε", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(x + 130, 336, "curves: x = 0, 0.5, 1", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 130, 354, "factor → 1 for a perfect conductor", size=9, fill="var(--sf-muted)", weight=450))
    x = 624
    b.append(rect(x, 50, 252, 330, rx=14, sw=1.3))
    b.append(text(x + 126, 74, "Born radii (HCT + OBC)", size=11, weight=700))
    b.append(text(x + 126, 96, "1/Rᵢ = 1/ρᵢ − tanh(αΨ − βΨ² + γΨ³)/ρᵢ", size=9, weight=600))
    b.append(text(x + 126, 114, "α 1.0 · β 0.8 · γ 4.85", size=9, fill="var(--sf-muted)", weight=500))
    cx, cy = x + 126, 196
    for rr, op in [(66, 0.10), (46, 0.16)]:
        b.append(circle(cx, cy, rr, fill="var(--sf-accent)", opacity=op))
    b.append(water_2d(slug, cx, cy, 40, r=11))
    b.append(text(cx + 50, cy - 50, "buried → larger R", size=8.5, fill="var(--sf-accent)", weight=600, anchor="start"))
    b.append(text(x + 126, 290, "ΔG_np = γ_SA · SASA", size=10, weight=700))
    b.append(text(x + 126, 308, "water SASA = 123.4 Å² (sci-form)", size=9.5, weight=500))
    b.append(text(x + 126, 334, "vs tblite 0.7.0 ALPB:", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 126, 350, "electrostatic within 3.14 kcal/mol (tol. 5)", size=9, fill="var(--sf-warm)", weight=600))
    write(slug, 900, 396, "".join(b), label="Generalized Born with ALPB: the GB kernel interpolating from the Born self term to the Coulomb distance, the ALPB scaling factor versus solvent dielectric, and OBC Born radii with the non-polar SASA term.")


def tet_min_volume(tilt_deg: float) -> float:
    """Min |triple product| of unit neighbor→center vectors for a C3v centre.

    tilt = angle of the three base bonds below the plane through the centre
    (19.47° is the ideal tetrahedron, 0° a flat centre); the 4th bond is axial.
    """
    t = math.radians(tilt_deg)
    vs = [(math.cos(t) * math.cos(math.radians(a)), math.cos(t) * math.sin(math.radians(a)), -math.sin(t)) for a in (0, 120, 240)] + [(0.0, 0.0, 1.0)]
    def trip(a, b_, c):
        return a[0] * (b_[1] * c[2] - b_[2] * c[1]) - a[1] * (b_[0] * c[2] - b_[2] * c[0]) + a[2] * (b_[0] * c[1] - b_[1] * c[0])
    import itertools
    return min(abs(trip(*cmb)) for cmb in itertools.combinations(vs, 3))


def diagram_experimental_d4() -> None:
    slug = "experimental-d4"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "D4 dispersion: BJ-damped two-body term, CN-weighted C₆ and the three-body ATM sign", size=15, weight=700))
    a1, a2, s8 = 0.3957, 4.4628, 0.7875  # D4(PBE), a2 in bohr
    C6, ratio = 18.0, 20.0  # illustrative C–C values (Ha·bohr⁶, bohr²)
    C8 = C6 * ratio
    R0 = a1 * math.sqrt(ratio) + a2
    def e2(r):
        return -(C6 / (r ** 6 + R0 ** 6) + s8 * C8 / (r ** 8 + R0 ** 8)) * 627.5
    def raw(r):
        return -(C6 / r ** 6) * 627.5
    energy_panel(b, 24, 50, 330, 330, "C–C pair, D4(PBE) BJ damping", "−s₆C₆/(r⁶ + R₀⁶) − s₈C₈/(r⁸ + R₀⁸)",
                 [(raw, "var(--sf-line)", 1.4, "4 3"), (e2, "var(--sf-warm)", 2.4, None)], (2.0, 14.0), (-0.9, 0.05), xlabel="r (bohr)",
                 marks=[(R0, f"R₀ = {R0:.2f}", "var(--sf-warm)")], note="dashed: undamped −C₆/r⁶ diverges; BJ stays finite")
    # CN weights
    x = 374
    b.append(rect(x, 50, 250, 330, rx=14, sw=1.3))
    b.append(text(x + 125, 74, "C₆ follows the coordination number", size=11, weight=700))
    b.append(text(x + 125, 94, "C₆ᴬᴮ = Σ W_a W_b C₆^ref", size=10, weight=600))
    px, py, pw, ph = x + 30, 112, 200, 150
    refs = [0.0, 0.99, 1.99, 2.99, 3.98]
    cols = ["var(--sf-line)", "var(--sf-accent)", "var(--sf-accent2)", "var(--sf-violet)", "var(--sf-warm)"]
    sx = lambda c: px + c / 4.5 * pw
    sy = lambda v: py + ph - v * ph
    def weights(cn):
        g = [math.exp(-6 * (cn - r_) ** 2) for r_ in refs]
        z = sum(g) or 1
        return [v / z for v in g]
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    for k in range(5):
        pts = [(sx(i * 4.5 / 120), sy(weights(i * 4.5 / 120)[k])) for i in range(121)]
        b.append(path("M" + " L".join(f"{u:.1f},{v:.1f}" for u, v in pts), stroke=cols[k], sw=1.8))
        b.append(text(sx(refs[k]), py - 4, f"{refs[k]:.0f}", size=8.5, fill=cols[k], weight=700))
    for c in (0, 1, 2, 3, 4):
        b.append(text(sx(c), py + ph + 13, str(c), size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(px + pw, py + ph + 26, "CN of carbon", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(x + 125, 322, "Gaussian weights over reference", size=9.5, weight=500))
    b.append(text(x + 125, 340, "CNs (sp → sp³ carbon)", size=9.5, weight=500))
    b.append(text(x + 125, 364, "same Fermi CN as EEQ", size=9, fill="var(--sf-muted)", weight=450))
    # ATM sign
    x = 644
    b.append(rect(x, 50, 232, 330, rx=14, sw=1.3))
    b.append(text(x + 116, 74, "three-body ATM", size=11, weight=700))
    b.append(text(x + 116, 94, "∝ (3 cosθ_A cosθ_B cosθ_C + 1)/(r³r³r³)", size=8.5, weight=600))
    for k, (pts, lbl, val, col) in enumerate([
        ([(x + 56, 180), (x + 106, 180), (x + 81, 137)], "equilateral", 3 * 0.125 + 1, "var(--sf-rose)"),
        ([(x + 136, 170), (x + 166, 170), (x + 196, 170)], "linear", 3 * 1 * 1 * -1 + 1, "var(--sf-accent2)"),
    ]):
        b.append(path("M" + " L".join(f"{u},{v}" for u, v in pts) + " Z", stroke=col, sw=1.8, fill=col, opacity=0.9))
        for (u, v) in pts:
            b.append(atom_ball(slug, u, v, "C", 7, label=False))
        cxk = sum(u for u, _ in pts) / 3
        b.append(text(cxk, 210, lbl, size=9.5, weight=600))
        b.append(text(cxk, 228, f"factor {val:+.3f}", size=9.5, fill=col, weight=700))
    b.append(text(x + 116, 262, "equilateral → repulsive (+)", size=9.5, fill="var(--sf-rose)", weight=600))
    b.append(text(x + 116, 280, "linear → attractive (−)", size=9.5, fill="var(--sf-accent2)", weight=600))
    b.append(text(x + 116, 316, "vs dftd4 4.3.0: D4(PBE)", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 116, 332, "6.2e−10 Ha · D3(HF-3c) 6.8e−9", size=9, fill="var(--sf-accent2)", weight=700))
    b.append(text(x + 116, 362, "C₆, C₈ in the left plot illustrative", size=8.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="D4 dispersion: the Becke-Johnson damped C6 and C8 pair energy with PBE parameters compared with the divergent undamped term, Gaussian coordination-number weights for carbon references, and the sign of the three-body ATM term for equilateral and linear triples.")


def diagram_experimental_sdr() -> None:
    import numpy as np
    slug = "experimental-sdr"
    begin(slug)
    b = [element_defs(slug)]
    rng = np.random.default_rng(5)
    P = rng.uniform(-2, 2, (10, 3))
    n = len(P)
    D2 = ((P[:, None] - P[None]) ** 2).sum(-1)
    known = [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < 0.7 or j == i + 1]
    mask = np.zeros((n, n), bool)
    for i, j in known:
        mask[i, j] = mask[j, i] = True
    np.fill_diagonal(mask, True)
    Dk = np.where(mask, D2, D2[mask].mean())
    Jc = np.eye(n) - np.ones((n, n)) / n
    hist, spectra = [], []
    Gold = None
    for it in range(3000):
        G = -0.5 * Jc @ Dk @ Jc
        w, V = np.linalg.eigh(G)
        if it == 0:
            spectra.append(np.sort(w)[::-1])
        w3 = np.where(np.arange(n) >= n - 3, np.clip(w, 0, None), 0.0)
        X = (V * w3) @ V.T
        dg = np.diag(X)
        Dn = dg[:, None] + dg[None] - 2 * X
        Dk = np.where(mask, D2, Dn)
        if Gold is not None:
            hist.append(np.linalg.norm(X - Gold) / np.linalg.norm(Gold))
            if hist[-1] < 1e-6:
                break
        Gold = X
    spectra.append(np.sort(np.linalg.eigvalsh(X))[::-1])
    w, V = np.linalg.eigh(X)
    C = V[:, -3:] * np.sqrt(np.clip(w[-3:], 0, None))
    Dr = ((C[:, None] - C[None]) ** 2).sum(-1)
    err_known = max(abs(math.sqrt(Dr[i, j]) - math.sqrt(D2[i, j])) for i, j in known)
    err_all = max(abs(math.sqrt(max(Dr[i, j], 0)) - math.sqrt(D2[i, j])) for i in range(n) for j in range(i + 1, n))
    b.append(text(450, 30, f"SDR embedding by alternating projections — run for this figure: {n} points, {len(known)} of {n*(n-1)//2} distances", size=14, weight=700))
    # left: known-distance graph
    b.append(rect(24, 50, 250, 330, rx=14, sw=1.3))
    b.append(text(149, 74, f"{len(known)} of {n*(n-1)//2} distances known", size=11, weight=700))
    pts = [(149 + p_[0] * 46, 206 - p_[1] * 46) for p_ in P]
    for i, j in known:
        b.append(path(f"M{pts[i][0]:.1f},{pts[i][1]:.1f} L{pts[j][0]:.1f},{pts[j][1]:.1f}", stroke="var(--sf-accent)", sw=1.2, opacity=0.6))
    for (u, v) in pts:
        b.append(circle(u, v, 6, fill="var(--sf-accent)"))
    b.append(text(149, 350, "dᵢⱼ² = Xᵢᵢ + Xⱼⱼ − 2Xᵢⱼ", size=10, weight=600))
    b.append(text(149, 368, "X = Gram matrix (unknown)", size=9, fill="var(--sf-muted)", weight=450))
    # middle: convergence
    x = 294
    b.append(rect(x, 50, 300, 330, rx=14, sw=1.3))
    b.append(text(x + 150, 74, "‖X_{k+1} − X_k‖_F / ‖X_k‖_F", size=10.5, weight=700))
    px, py, pw, ph = x + 50, 92, 230, 200
    m = len(hist)
    sx = lambda k: px + k / max(m - 1, 1) * pw
    sy = lambda v: py + (0 - math.log10(max(v, 1e-9))) / 8 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    b.append(path("M" + " L".join(f"{sx(k):.1f},{sy(v):.1f}" for k, v in enumerate(hist)), stroke="var(--sf-accent2)", sw=2.2, cls="sf-draw"))
    b.append(path(f"M{px},{sy(1e-6):.1f} L{px+pw},{sy(1e-6):.1f}", stroke="var(--sf-rose)", sw=1.1, dash="4 4"))
    b.append(text(px + pw - 2, sy(1e-6) - 5, "stop 1e−6", size=9, fill="var(--sf-rose)", weight=700, anchor="end"))
    for e_ in (0, -2, -4, -6, -8):
        b.append(text(px - 6, sy(10 ** e_) + 3, f"1e{e_}" if e_ else "1", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(px + pw, py + ph + 14, f"{m} iterations", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(x + 150, 330, "P_Ω: restore known dᵢⱼ² ↔ P₊: keep top-3 λ ≥ 0", size=9.5, weight=600))
    b.append(text(x + 150, 350, f"max error: known {err_known:.1e} Å · unknown {err_all:.1e} Å", size=9.5, fill="var(--sf-accent2)", weight=700))
    # right: spectrum before/after
    x = 614
    b.append(rect(x, 50, 262, 330, rx=14, sw=1.3))
    b.append(text(x + 131, 74, "eigenvalues of X", size=11, weight=700))
    for k, (sp, lbl, col) in enumerate([(spectra[0], "Gram of the initial guess", "var(--sf-rose)"), (spectra[-1], "converged", "var(--sf-accent2)")]):
        bx0, by0 = x + 24, 110 + k * 130
        mx = max(abs(spectra[0]).max(), spectra[-1].max())
        b.append(text(bx0, by0 - 6, lbl, size=9.5, fill=col, weight=700, anchor="start"))
        b.append(path(f"M{bx0},{by0+50} L{bx0+214},{by0+50}", sw=0.8))
        for i, v in enumerate(sp):
            h = v / mx * 48
            b.append(rect(bx0 + 6 + i * 26, by0 + 50 - max(h, 0), 18, max(abs(h), 1.2), rx=2, fill=col if v >= 0 else "var(--sf-rose)", stroke="none", opacity=0.8))
    b.append(text(x + 131, 362, "rank 3 remains → C = V₃ Σ₃^½ (unknown dᵢⱼ recovered)", size=9, weight=600))
    write(slug, 900, 396, "".join(b), label="Semidefinite relaxation run on ten random points with most distances known: alternating distance and rank-3 PSD projections converge, negative eigenvalues vanish and the recovered coordinates also reproduce the unknown distances.")


def diagram_experimental_cpm() -> None:
    import numpy as np
    slug = "experimental-cpm"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Constant-potential charges: solve M q = μ − χ (toy two-plate cell computed for this figure)", size=14, weight=700))
    # two 3x3 plates 6 Å apart
    pos = [(x, y, z) for z in (0.0, 6.0) for x in (0, 2.46, 4.92) for y in (0, 2.46, 4.92)]
    n = len(pos)
    eta, chi = 10.0, 0.0
    M = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            M[i, j] = eta if i == j else 14.3996 / math.dist(pos[i], pos[j])
    def solve(V):
        mu = np.array([V / 2 if p_[2] == 0 else -V / 2 for p_ in pos])
        return np.linalg.solve(M, mu - chi)
    q2 = solve(2.0)
    Vs = np.linspace(-2, 2, 9)
    Q = [solve(V)[:9].sum() for V in Vs]
    Ccap = (Q[-1] - Q[0]) / (Vs[-1] - Vs[0])
    J = M - np.diag(np.diag(M))
    rho = max(abs(np.linalg.eigvals(J / eta)))
    # left: plates with charges
    b.append(rect(24, 50, 300, 330, rx=14, sw=1.3))
    b.append(text(174, 74, "plates at +V/2 and −V/2 (V = 2 V)", size=11, weight=700))
    for k, (z, yb) in enumerate([(0.0, 120), (6.0, 270)]):
        for i in range(3):
            for j in range(3):
                idx = k * 9 + i * 3 + j
                qv = q2[idx]
                u = 86 + i * 60 + j * 22
                v = yb + j * 10
                col = "var(--sf-accent)" if qv > 0 else "var(--sf-rose)"
                b.append(circle(u, v, 6 + abs(qv) * 40, fill=col, opacity=0.8, cls="sf-pulse"))
        b.append(text(300, yb + 14, "+V/2" if k == 0 else "−V/2", size=10, fill="var(--sf-accent)" if k == 0 else "var(--sf-rose)", weight=700, anchor="end"))
    b.append(path("M174,170 L174,250", stroke="var(--sf-line)", sw=1.2, dash="3 3"))
    b.append(text(184, 214, "6 Å", size=9, fill="var(--sf-muted)", anchor="start"))
    b.append(text(174, 330, f"Q_plate = {q2[:9].sum():+.3f} e  /  {q2[9:].sum():+.3f} e", size=10, weight=700))
    b.append(text(174, 350, "edge atoms carry more charge", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(174, 366, "η = 10, χ = 0, ε = 1 (toy values)", size=9, fill="var(--sf-muted)", weight=450))
    # middle: Q(V)
    energy_panel(b, 344, 50, 270, 330, "charge response Q(V)", f"C = ∂Q/∂V = {Ccap:.3f} e/V",
                 [(lambda V: float(np.interp(V, Vs, Q)), "var(--sf-accent2)", 2.4, None)], (-2, 2), (min(Q) * 1.1, max(Q) * 1.1), xlabel="V (volts)",
                 note="linear: one LU solve per potential")
    # right: why direct solve
    x = 634
    b.append(rect(x, 50, 242, 330, rx=14, sw=1.3))
    b.append(text(x + 121, 74, "why LU, not fixed point", size=11, weight=700))
    b.append(text(x + 121, 98, "Mᵢᵢ = ηᵢ, Mᵢⱼ = 14.3996/(ε rᵢⱼ)", size=9.5, weight=600))
    b.append(text(x + 121, 130, "fixed point q ← (q + q_new)/2", size=9.5, weight=500))
    b.append(text(x + 121, 148, "converges only if ρ(J/η) < 1", size=9.5, weight=500))
    b.append(text(x + 121, 178, f"this toy: ρ(J/η) = {rho:.2f}", size=11, fill="var(--sf-rose)" if rho > 1 else "var(--sf-accent2)", weight=800))
    b.append(text(x + 121, 200, "(diverges)" if rho > 1 else "(converges)", size=9.5, fill="var(--sf-rose)" if rho > 1 else "var(--sf-accent2)", weight=600))
    b.append(text(x + 121, 232, "graphene bilayer: charges → 1e157", size=9.5, fill="var(--sf-rose)", weight=600))
    b.append(text(x + 121, 250, "before 2026-10-13 → LU since", size=9.5, weight=500))
    b.append(text(x + 121, 286, "vs LAMMPS electrode/conp:", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 121, 302, "plate charge scale 1.76×,", size=9, fill="var(--sf-warm)", weight=600))
    b.append(text(x + 121, 318, "per-atom cosine 0.44", size=9, fill="var(--sf-warm)", weight=600))
    b.append(text(x + 121, 348, "(experimental: not yet parity)", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="Constant-potential method on a toy two-plate cell: the linear system gives opposite plate charges with larger edge charges, a linear charge-voltage response whose slope is the capacitance, and the spectral-radius reason the damped fixed point was replaced by a direct LU solve.")


def diagram_experimental_riemannian() -> None:
    slug = "experimental-riemannian"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Optimizing on the PSD set: the trace-2 slice of 2 × 2 symmetric matrices is a disk", size=15, weight=700))
    # X = [[1+u, v],[v, 1-u]] is PSD  ⇔  u² + v² ≤ 1
    cx, cy, R = 250, 220, 130
    su = lambda u: cx + u * R
    sv = lambda v: cy - v * R
    b.append(rect(24, 50, 452, 340, rx=14, sw=1.3))
    b.append(circle(cx, cy, R, fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1.6, opacity=0.14))
    b.append(circle(cx, cy, R, fill="none", stroke="var(--sf-accent)", sw=1.6))
    b.append(text(cx, cy + R + 22, "u² + v² ≤ 1  ⇔  X ⪰ 0 · boundary = rank 1", size=9.5, fill="var(--sf-accent)", weight=600))
    D = (1.45, 0.55)
    b.append(circle(su(D[0]), sv(D[1]), 6, fill="var(--sf-rose)"))
    b.append(text(su(D[0]) - 4, sv(D[1]) - 10, "target D (indefinite)", size=9.5, fill="var(--sf-rose)", weight=700, anchor="end"))
    # projected gradient iterates on f = ||X − D||² from a start inside
    u, v = -0.5, -0.35
    pts = [(u, v)]
    raw = []
    for _ in range(12):
        gu, gv = 2 * (u - D[0]), 2 * (v - D[1])
        nu, nv = u - 0.12 * gu, v - 0.12 * gv
        raw.append((nu, nv))
        r_ = math.hypot(nu, nv)
        if r_ > 1:
            nu, nv = nu / r_, nv / r_
        u, v = nu, nv
        pts.append((u, v))
    for (a_, b_), (ru, rv) in zip(pts[:-1], raw):
        if math.hypot(ru, rv) > 1:
            b.append(path(f"M{su(a_):.1f},{sv(b_):.1f} L{su(ru):.1f},{sv(rv):.1f}", stroke="var(--sf-warm)", sw=1.2, dash="3 3"))
            b.append(circle(su(ru), sv(rv), 3, fill="var(--sf-warm)"))
    b.append(path("M" + " L".join(f"{su(a_):.1f},{sv(b_):.1f}" for a_, b_ in pts), stroke="var(--sf-violet)", sw=2.2))
    for a_, b_ in pts:
        b.append(circle(su(a_), sv(b_), 4, fill="var(--sf-violet)"))
    opt = (D[0] / math.hypot(*D), D[1] / math.hypot(*D))
    b.append(circle(su(opt[0]), sv(opt[1]), 7, fill="none", stroke="var(--sf-accent2)", sw=2.4, cls="sf-pulse"))
    b.append(text(su(opt[0]) + 10, sv(opt[1]) + 18, "optimum = P₊(D)", size=9.5, fill="var(--sf-accent2)", weight=700, anchor="start"))
    b.append(text(su(pts[0][0]) - 6, sv(pts[0][1]) + 16, "start", size=9.5, fill="var(--sf-violet)", weight=700, anchor="end"))
    b.append(text(60, 80, "dashed orange: raw step leaves the set", size=9, fill="var(--sf-warm)", weight=600, anchor="start"))
    b.append(text(60, 96, "violet: retracted iterates stay feasible", size=9, fill="var(--sf-violet)", weight=600, anchor="start"))
    x = 496
    b.append(rect(x, 50, 380, 340, rx=14, sw=1.3))
    b.append(text(x + 16, 76, "on the real problem (metric matrix G)", size=11, weight=700, anchor="start"))
    rows = [
        ("manifold", "{X ⪰ 0, rank X = r}"),
        ("Riemannian gradient", "grad f = P_T(∇f)"),
        ("Stiefel tangent proj.", "ξ − X·sym(Xᵀξ)"),
        ("retraction", "P₊(X + ξ + ½ ξ X⁻¹ ξ)"),
        ("P₊", "zero the negative eigenvalues"),
        ("L-BFGS", "m = 5 pairs, vector transport"),
        ("line search", "Armijo along the geodesic"),
    ]
    for i, (k_, v_) in enumerate(rows):
        y = 104 + i * 32
        b.append(rect(x + 14, y - 16, 352, 26, rx=7, fill="var(--sf-card2)", sw=0.6))
        b.append(text(x + 24, y + 1, k_, size=10, fill="var(--sf-violet)", weight=700, anchor="start"))
        b.append(text(x + 356, y + 1, v_, size=10, weight=500, anchor="end"))
    b.append(text(x + 190, 344, "vs SciPy L-BFGS-B: projection 8.2e−12,", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 190, 362, "objective/gradient 1.2e−16", size=9.5, fill="var(--sf-accent2)", weight=700))
    write(slug, 900, 406, "".join(b), label="Optimization on the positive-semidefinite set shown exactly as the unit disk of trace-2, 2 by 2 matrices: raw gradient steps toward an indefinite target leave the disk, retraction brings them back and the iterates converge to the PSD projection of the target.")


def diagram_validation_pipeline() -> None:
    slug = "validation-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "Five geometric checks after minimization — any failure sends the attempt back to embedding", size=15, weight=700))
    W, H, Y = 196, 236, 54
    cards = [("1", "Energy / atom", "var(--sf-accent)"), ("2", "Tetrahedral", "var(--sf-accent2)"), ("3", "Chiral sign", "var(--sf-violet)"), ("4", "Planarity", "var(--sf-warm)"), ("5", "Double bond", "var(--sf-rose)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    info = [
        ("E_total / N < 0.05", "bounds-FF energy after BFGS", "catches: bad embedding"),
        ("min |Vₖ| > 0.50", "(0.25 in small rings)", "C/N, degree 4, ≥ 2 rings"),
        ("sign V = sign V_target", "|V| > 0.2 |V_target|", "@ → V > 0 · @@ → V < 0"),
        ("E_oop < 0.7 · N_improper", "sp: angle > 175°", "catches: puckered sp²"),
        ("cos θ + 1 > 10⁻³", "θ: substituent vs C=C", "catches: linear substituent"),
    ]
    for i, (a_, b2, c_) in enumerate(info):
        x = xs[i]
        b.append(rect(x + 12, Y + 52, W - 24, 34, rx=8, fill="var(--sf-soft)", stroke=cards[i][2], sw=1.1))
        b.append(text(x + W / 2, Y + 74, a_, size=10.5, weight=700))
        b.append(text(x + W / 2, Y + 108, b2, size=9.5, weight=500))
        b.append(text(x + W / 2, Y + 128, c_, size=9.5, fill="var(--sf-muted)", weight=450))
        b.append(circle(x + 66, Y + 176, 13, fill="var(--sf-accent2)"))
        b.append(path(f"M{x+60},{Y+176} l4,4 l8,-8", stroke="#fff", sw=2))
        b.append(text(x + 66, Y + 204, "pass", size=9, fill="var(--sf-accent2)", weight=700))
        b.append(circle(x + 130, Y + 176, 13, fill="var(--sf-rose)"))
        b.append(path(f"M{x+125},{Y+171} l10,10 M{x+135},{Y+171} l-10,10", stroke="#fff", sw=2))
        b.append(text(x + 130, Y + 204, "retry", size=9, fill="var(--sf-rose)", weight=700))
        b.append(path(f"M{x+130},{Y+189} L{x+130},{Y+H+18}", stroke="var(--sf-rose)", sw=1.2, dash="3 3"))
    y = Y + H + 18
    b.append(retry_path(slug, f"M{xs[4]+130},{y} L40,{y} L40,{Y + H / 2} L22,{Y + H / 2}".replace("L22", "L23")))
    b.append(text(xs[2] + 98, y + 22, "retry budget 10N attempts (30 atoms → 300) · most molecules pass in 1–5, strained / many stereocentres 10–50", size=10, fill="var(--sf-rose)", weight=600))
    write(slug, 1100, y + 40, "".join(b), label="Conformer validation: energy per atom below 0.05, tetrahedral triple products above 0.50, chiral volume sign with 20 percent tolerance, out-of-plane energy below 0.7 per improper and double-bond substituent check, each failure retrying within a 10N budget.")


def diagram_validation_tetrahedral() -> None:
    slug = "validation-tetrahedral"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Tetrahedral check: smallest triple product of unit vectors vs. how flat the centre is", size=15, weight=700))
    ideal = tet_min_volume(19.47)
    # left: two centres drawn
    for k, (tilt, title, col) in enumerate([(19.47, "ideal sp³ (tilt 19.5°)", "var(--sf-accent2)"), (3.0, "flattened (tilt 3°)", "var(--sf-rose)")]):
        x = 24 + k * 220; cx, cy = x + 100, 196
        b.append(rect(x, 50, 200, 300, rx=14, sw=1.3))
        b.append(text(x + 100, 74, title, size=11, fill=col, weight=700))
        t = math.radians(tilt)
        L = 62
        base = [(math.cos(t) * math.cos(math.radians(a)), math.cos(t) * math.sin(math.radians(a)), -math.sin(t)) for a in (200, 320, 80)]
        proj = lambda v: (cx + L * (v[0] * 0.95 + v[1] * 0.25), cy - L * (v[2] * 1.0 + v[1] * 0.35))
        top = proj((0, 0, 1))
        pts = [proj(v) for v in base]
        b.append(path(f"M{pts[0][0]:.1f},{pts[0][1]:.1f} L{pts[1][0]:.1f},{pts[1][1]:.1f} L{pts[2][0]:.1f},{pts[2][1]:.1f} Z", fill=col, stroke=col, sw=1, opacity=0.15))
        for q in pts + [top]:
            b.append(path(f"M{cx},{cy} L{q[0]:.1f},{q[1]:.1f}", stroke="var(--sf-ink)", sw=2.2))
        for q in pts + [top]:
            b.append(atom_ball(slug, q[0], q[1], "C", 8, label=False))
        b.append(atom_ball(slug, cx, cy, "C", 11))
        v = tet_min_volume(tilt)
        b.append(text(x + 100, 300, f"min |Vₖ| = {v:.3f}", size=12, fill=col, weight=800))
        b.append(text(x + 100, 322, "accept" if v > 0.5 else "reject (< 0.50)", size=10, fill=col, weight=700))
    # right: curve
    x = 464
    b.append(rect(x, 50, 412, 300, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "min |V̂_AC · (V̂_BC × V̂_DC)| over the 4 triples", size=10.5, weight=700, anchor="start"))
    px, py, pw, ph = x + 50, 90, 330, 200
    sx = lambda tdeg: px + tdeg / 25 * pw
    sy = lambda v: py + ph - v / 0.85 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    b.append(rect(px, sy(0.5), pw, py + ph - sy(0.5), rx=0, fill="var(--sf-rose)", stroke="none", opacity=0.07))
    b.append(path("M" + " L".join(f"{sx(i/4):.1f},{sy(tet_min_volume(i/4)):.1f}" for i in range(0, 101)), stroke="var(--sf-ink)", sw=2.4))
    for th, lbl, col in [(0.5, "0.50 normal", "var(--sf-rose)"), (0.25, "0.25 small rings", "var(--sf-warm)")]:
        b.append(path(f"M{px},{sy(th):.1f} L{px+pw},{sy(th):.1f}", stroke=col, sw=1.2, dash="5 4"))
        b.append(text(px + pw - 4, sy(th) - 5, lbl, size=9, fill=col, weight=700, anchor="end"))
    b.append(circle(sx(19.47), sy(ideal), 5, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(sx(19.47) - 6, sy(ideal) - 10, f"ideal {ideal:.3f}", size=9.5, fill="var(--sf-accent2)", weight=700, anchor="end"))
    b.append(circle(sx(3), sy(tet_min_volume(3)), 5, fill="var(--sf-rose)"))
    for tv in (0, 5, 10, 15, 20, 25):
        b.append(text(sx(tv), py + ph + 13, f"{tv}°", size=8.5, fill="var(--sf-muted)", weight=500))
    for vv in (0, 0.25, 0.5, 0.75):
        b.append(text(px - 6, sy(vv) + 3, f"{vv:g}", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(px + pw, py + ph + 26, "tilt of the three base bonds below the centre", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(text(x + 206, 340, "applies to C/N with degree 4 in ≥ 2 SSSR rings, not in a 3-ring · centre must also lie inside (tol 0.30)", size=8.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 366, "".join(b), label=f"Tetrahedral check: the smallest triple product of unit bond vectors is {ideal:.3f} for an ideal sp3 centre and falls below the 0.50 threshold as the centre flattens; a 3 degree tilt is rejected.")


def diagram_validation_double_bond() -> None:
    slug = "validation-double-bond"
    begin(slug)
    b = [element_defs(slug)]
    crit = math.degrees(math.acos(-1 + 1e-3))
    b.append(text(450, 30, f"Double-bond check: reject when a substituent is collinear with C=C (θ > {crit:.1f}°)", size=15, weight=700))
    for k, (th, title, col) in enumerate([(122.0, "normal sp² geometry", "var(--sf-accent2)"), (179.2, "collinear substituent", "var(--sf-rose)")]):
        x = 24 + k * 220; cy = 190
        b.append(rect(x, 50, 200, 300, rx=14, sw=1.3))
        b.append(text(x + 100, 74, title, size=11, fill=col, weight=700))
        c1, c2 = (x + 70, cy), (x + 140, cy)
        b.append(path(f"M{c1[0]},{c1[1]-3} L{c2[0]},{c2[1]-3} M{c1[0]},{c1[1]+3} L{c2[0]},{c2[1]+3}", stroke="var(--sf-ink)", sw=2))
        a = math.radians(th)
        # substituent on C1, angle θ measured between C1→sub and C1→C2 bond vector
        sub = (c1[0] + 52 * math.cos(math.pi - a + math.pi) * -1, c1[1] - 52 * math.sin(a))
        sub = (c1[0] + 52 * math.cos(a), c1[1] - 52 * math.sin(a))
        b.append(path(f"M{c1[0]},{c1[1]} L{sub[0]:.1f},{sub[1]:.1f}", stroke=col, sw=2.4))
        b.append(path(f"M{c1[0]+18},{c1[1]} A18,18 0 0 0 {c1[0] + 18 * math.cos(a):.1f},{c1[1] - 18 * math.sin(a):.1f}", stroke=col, sw=1.4))
        b.append(path(f"M{c2[0]},{c2[1]} l26,-40 M{c2[0]},{c2[1]} l26,40", stroke="var(--sf-ink)", sw=2))
        for q in (c1, c2):
            b.append(atom_ball(slug, q[0], q[1], "C", 11))
        b.append(atom_ball(slug, sub[0], sub[1], "C", 8, label=False, cls="sf-pulse"))
        b.append(text(c1[0] + 4, c1[1] - 26, "θ", size=11, fill=col, weight=700))
        v = math.cos(a) + 1
        b.append(text(x + 100, 290, f"θ = {th:.1f}°", size=11, weight=700))
        b.append(text(x + 100, 312, f"cos θ + 1 = {v:.2e}" if v < 0.01 else f"cos θ + 1 = {v:.3f}", size=11, fill=col, weight=800))
        b.append(text(x + 100, 332, "accept" if v > 1e-3 else "reject (< 10⁻³)", size=10, fill=col, weight=700))
    x = 464
    b.append(rect(x, 50, 412, 300, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "cos θ + 1 on a log scale", size=11, weight=700, anchor="start"))
    px, py, pw, ph = x + 56, 90, 320, 200
    sx = lambda tdeg: px + (tdeg - 110) / 70 * pw
    sy = lambda v: py + (0.5 - math.log10(max(v, 1e-6))) / 6.5 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    b.append(rect(sx(crit), py, px + pw - sx(crit), ph, rx=0, fill="var(--sf-rose)", stroke="none", opacity=0.1))
    b.append(path("M" + " L".join(f"{sx(110 + i * 0.35):.1f},{sy(math.cos(math.radians(110 + i * 0.35)) + 1):.1f}" for i in range(200)), stroke="var(--sf-ink)", sw=2.4))
    b.append(path(f"M{px},{sy(1e-3):.1f} L{px+pw},{sy(1e-3):.1f}", stroke="var(--sf-rose)", sw=1.2, dash="5 4"))
    b.append(text(px + 6, sy(1e-3) - 6, "threshold 10⁻³", size=9, fill="var(--sf-rose)", weight=700, anchor="start"))
    b.append(text(sx(crit) - 4, py + 14, f"θ > {crit:.2f}° rejected", size=9, fill="var(--sf-rose)", weight=700, anchor="end"))
    for tv in (110, 120, 140, 160, 180):
        b.append(text(sx(tv), py + ph + 13, f"{tv}°", size=8.5, fill="var(--sf-muted)", weight=500))
    for e_ in (0, -2, -4, -6):
        b.append(text(px - 6, sy(10 ** e_) + 3, f"1e{e_}" if e_ else "1", size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
    b.append(circle(sx(122), sy(math.cos(math.radians(122)) + 1), 5, fill="var(--sf-accent2)"))
    b.append(circle(sx(179.2), sy(math.cos(math.radians(179.2)) + 1), 5, fill="var(--sf-rose)"))
    b.append(text(x + 206, 340, "θ: angle between substituent bond and the double-bond vector, both ends checked", size=8.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 366, "".join(b), label=f"Double-bond check: a 122 degree substituent gives cos θ + 1 = 0.47 and passes, a 179.2 degree one gives about 1e-4 and is rejected; the 1e-3 threshold corresponds to θ above {crit:.1f} degrees.")


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


# ---------------------------------------------------------------------------
# Shared visual vocabulary for the redesigned diagrams. Every diagram should
# build on these so cards, atoms and loops look the same across the docs.
# ---------------------------------------------------------------------------
ELEMENT_FILL = {
    "C": ("#94a3b8", "#1e293b"),
    "O": ("#fda4af", "#be123c"),
    "N": ("#93c5fd", "#1d4ed8"),
    "H": ("#ffffff", "#cbd5e1"),
    "S": ("#fde68a", "#b45309"),
}


def element_defs(slug: str) -> str:
    """Shaded-sphere gradients for atoms plus a rose marker for retry loops."""
    out = "<defs>"
    for el, (hi, lo) in ELEMENT_FILL.items():
        out += (
            f'<radialGradient id="sf-{slug}-{el}" cx="35%" cy="35%" r="70%">'
            f'<stop offset="0%" stop-color="{hi}"/><stop offset="100%" stop-color="{lo}"/></radialGradient>'
        )
    out += (
        f'<marker id="sf-{slug}-arrow-rose" markerWidth="9" markerHeight="7" refX="8" refY="3.5" orient="auto">'
        '<path d="M0,0 L0,7 L9,3.5 Z" fill="var(--sf-rose)"/></marker>'
    )
    return out + "</defs>"


def atom_ball(slug: str, cx: float, cy: float, el: str, r: float = 11, *, label: bool = True, cls: str = "") -> str:
    out = circle(cx, cy, r, fill=f"url(#sf-{slug}-{el})", stroke="var(--sf-line)" if el == "H" else None, sw=0.6, cls=cls)
    if label and el != "H":
        out += text(cx, cy + r * 0.36, el, size=max(8, int(r * 0.95)), fill="#fff", weight=700)
    return out


def step_card(x: float, y: float, w: float, h: float, n: int, title: str, color: str) -> str:
    out = rect(x, y, w, h, rx=14, fill="var(--sf-card)", stroke="var(--sf-line)", sw=1.3)
    out += rect(x + 1.5, y + 1.5, w - 3, 4, rx=2, fill=color, stroke="none")
    out += circle(x + 20, y + 24, 10, fill=color)
    out += text(x + 20, y + 28, str(n), size=11, fill="#fff", weight=700)
    out += text(x + 38, y + 28, title, size=12, fill="var(--sf-ink)", weight=700, anchor="start")
    return out


def retry_path(slug: str, d: str) -> str:
    return (
        f'<path d="{d}" fill="none" stroke="var(--sf-rose)" stroke-width="1.6" stroke-dasharray="5 5" '
        f'stroke-linecap="round" stroke-linejoin="round" marker-end="url(#sf-{slug}-arrow-rose)" class="sf-flow"/>'
    )


def stage_row(b: list, x0: float, y: float, w: float, h: float, gap: float, cards: list) -> list[float]:
    """Row of tagged stage cards joined by flowing arrows. Returns card x positions."""
    xs = [x0 + i * (w + gap) for i in range(len(cards))]
    for i, (tag, title, color) in enumerate(cards):
        x = xs[i]
        b.append(rect(x, y, w, h, rx=14, sw=1.3))
        b.append(rect(x + 1.5, y + 1.5, w - 3, 4, rx=2, fill=color, stroke="none"))
        tw = 10 + 6.5 * len(tag)
        b.append(rect(x + 10, y + 14, tw, 20, rx=6, fill=color, stroke="none"))
        b.append(text(x + 10 + tw / 2, y + 28, tag, size=10.5, fill="#fff", weight=700))
        b.append(text(x + 18 + tw, y + 28, title, size=11, weight=700, anchor="start"))
        if i < len(cards) - 1:
            b.append(arrow_between(x + w + 2, y + h / 2, x + w + gap - 3, y + h / 2, color=color, delay=0.2 * i))
    return xs


CAPTION = dict(size=9.5, fill="var(--sf-muted)", weight=450)


# Real sci-form output for water (embed("O", seed 42) → EHT). Used by every
# electronic-structure diagram so the numbers agree across pages.
WATER_XYZ = [(-0.0057, 0.3852, 0.0), (-0.7961, -0.1947, 0.0), (0.8018, -0.1905, 0.0)]  # O, H, H (Å)
WATER_EHT = [-32.71, -15.90, -15.26, -14.80, -8.96, -8.26]  # eV, 8 electrons → HOMO = 4th
WATER_MULLIKEN = [-1.168, 0.584, 0.584]
WATER_LOWDIN = [-1.096, 0.548, 0.548]
WATER_DIPOLE_D = (0.045, -3.242, 0.0)


def water_2d(slug: str, cx: float, cy: float, k: float, *, charges: list | None = None, r: float = 13) -> str:
    out = ""
    pts = [(cx + x * k, cy - y * k) for x, y, _ in WATER_XYZ]
    for hx, hy in pts[1:]:
        out += path(f"M{pts[0][0]:.1f},{pts[0][1]:.1f} L{hx:.1f},{hy:.1f}", stroke="var(--sf-ink)", sw=3)
    out += atom_ball(slug, *pts[0], "O", r)
    for hx, hy in pts[1:]:
        out += atom_ball(slug, hx, hy, "H", r * 0.7)
    if charges:
        for (x, y), q in zip(pts, charges):
            col = "var(--sf-rose)" if q < 0 else "var(--sf-accent)"
            out += text(x, y + (-r - 8 if q < 0 else r + 16), f"{q:+.3f}", size=10, fill=col, weight=700)
    return out


def diagram_pipeline_overview() -> None:
    slug = "pipeline-overview"
    begin(slug)
    b = [element_defs(slug)]
    cols, rows, W, H = [40, 330, 620], [78, 288, 498], 260, 150
    phases = [
        ("Phase 1 · Topology → Bounds", "var(--sf-accent)", ["Parse SMILES", "Molecular graph", "Bounds matrix + smoothing"]),
        ("Phase 2 · Embedding → Optimization", "var(--sf-violet)", ["Sample distances", "Metric matrix → eigenvectors", "BFGS bounds minimization"]),
        ("Phase 3 · Refinement → Output", "var(--sf-accent2)", ["ETKDG force field", "Validation checks", "3D conformer"]),
    ]
    muted = dict(size=10, fill="var(--sf-muted)", weight=450)
    b.append(text(460, 32, "ETKDGv2 conformer pipeline — ethanol (CCO) as the running example", size=15, weight=700))
    for r, (label, color, titles) in enumerate(phases):
        y = rows[r]
        b.append(text(40, y - 14, label, size=12, fill=color, weight=700, anchor="start"))
        for c in range(3):
            b.append(step_card(cols[c], y, W, H, r * 3 + c + 1, titles[c], color))
        for c in range(2):
            b.append(arrow_between(cols[c] + W + 4, y + H / 2, cols[c + 1] - 6, y + H / 2, color=color, delay=0.2 * (r * 2 + c)))
    for r in range(2):
        y1, y2 = rows[r] + H, rows[r + 1]
        b.append(connector(f"M750,{y1+4} L750,{y1+22} L170,{y1+22} L170,{y2-26}", delay=0.5 + r))

    # 1 · SMILES tokens
    x, y = cols[0], rows[0]
    for i, (tok, col) in enumerate([("C", "var(--sf-ink)"), ("C", "var(--sf-ink)"), ("O", "var(--sf-rose)")]):
        b.append(rect(x + 70 + i * 42, y + 50, 34, 34, rx=6, fill="var(--sf-card2)"))
        b.append(text(x + 87 + i * 42, y + 73, tok, size=18, fill=col, weight=700, mono=True))
    b.append(text(x + 130, y + 110, '"CCO" → three heavy-atom tokens', **muted))
    b.append(text(x + 130, y + 128, "implicit H filled from valence", **muted))

    # 2 · Molecular graph with explicit hydrogens
    x, y = cols[1], rows[0]
    A = {"C1": (x + 80, y + 92), "C2": (x + 130, y + 70), "O": (x + 180, y + 92)}
    for hx, hy, p in [(x + 52, y + 78, "C1"), (x + 64, y + 120, "C1"), (x + 94, y + 124, "C1"),
                      (x + 116, y + 46, "C2"), (x + 148, y + 48, "C2"), (x + 208, y + 78, "O")]:
        b.append(path(f"M{A[p][0]},{A[p][1]} L{hx},{hy}", sw=1.3))
        b.append(atom_ball(slug, hx, hy, "H", 5))
    b.append(path(f"M{A['C1'][0]},{A['C1'][1]} L{A['C2'][0]},{A['C2'][1]} L{A['O'][0]},{A['O'][1]}", stroke="var(--sf-ink)", sw=2.4))
    for k, (ax, ay) in A.items():
        b.append(atom_ball(slug, ax, ay, k[0]))
    b.append(text(x + 130, y + 140, "9 atoms · 8 bonds · all sp³", **muted))

    # 3 · Bounds matrix: upper triangle = u, lower = l, smoothed by Floyd–Warshall
    x, y = cols[2], rows[0]
    mx, my, cs = x + 28, y + 46, 22
    labs = ["C", "C", "O", "H"]
    for i in range(4):
        b.append(text(mx - 9, my + i * cs + 15, labs[i], size=9, fill="var(--sf-muted)"))
        b.append(text(mx + i * cs + 10, my - 4, labs[i], size=9, fill="var(--sf-muted)"))
        for j in range(4):
            fill = "var(--sf-card2)" if i == j else ("var(--sf-violet)" if j > i else "var(--sf-accent)")
            op = 0.0 if i == j else 0.25 + 0.6 * abs(i - j) / 3
            b.append(rect(mx + j * cs, my + i * cs, cs - 2, cs - 2, rx=3, fill=fill, sw=0.6, opacity=round(op, 2) if op else None))
    b.append(rect(mx - 2, my - 2, 4 * cs + 2, 4 * cs + 2, rx=4, fill="none", stroke="var(--sf-warm)", sw=1.6, cls="sf-blink"))
    lx = x + 132
    for k, (s_, c_) in enumerate([("upper uᵢⱼ", "var(--sf-violet)"), ("lower lᵢⱼ", "var(--sf-accent)")]):
        b.append(rect(lx, y + 52 + k * 18, 10, 10, rx=2, fill=c_, stroke="none", opacity=0.7))
        b.append(text(lx + 16, y + 61 + k * 18, s_, size=10, weight=500, anchor="start"))
    b.append(text(lx, y + 102, "1-2 bonds, 1-3 angles,", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(lx, y + 116, "1-4 torsions, VDW", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(lx, y + 136, "Floyd–Warshall smoothing", size=9.5, fill="var(--sf-warm)", anchor="start"))

    # 4 · Random distances inside [l, u]
    x, y = cols[0], rows[1]
    for k, (lo, hi, d) in enumerate([(0.30, 0.42, 0.36), (0.45, 0.60, 0.55), (0.40, 0.88, 0.62)]):
        yy, x0, L = y + 60 + k * 26, x + 40, 190
        b.append(path(f"M{x0},{yy} L{x0+L},{yy}", sw=1))
        b.append(rect(x0 + lo * L, yy - 5, (hi - lo) * L, 10, rx=5, fill="var(--sf-violet)", stroke="none", opacity=0.28))
        b.append(circle(x0 + d * L, yy, 4.5, fill="var(--sf-warm)", cls=f"sf-pulse sf-d{k+1}"))
        b.append(text(x0 - 8, yy + 4, ["d₁₂", "d₁₃", "d₁₄"][k], size=9, fill="var(--sf-muted)", anchor="end"))
    b.append(text(x + 130, y + 140, "dᵢⱼ ~ U[lᵢⱼ, uᵢⱼ]   (MinstdRand)", **muted))

    # 5 · Metric matrix and eigenvalue spectrum
    x, y = cols[1], rows[1]
    b.append(text(x + 72, y + 60, "T = ½(D₀ᵢ + D₀ⱼ − dᵢⱼ²)", size=10, weight=600))
    b.append(text(x + 72, y + 84, "power iteration:", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 72, y + 98, "keep top 3 (or 4D)", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 72, y + 114, "xᵢₖ = √λₖ · vᵢₖ", size=10, fill="var(--sf-accent2)", weight=600))
    base = y + 110
    for k, v in enumerate([1.0, 0.72, 0.48, 0.12, 0.05, -0.03]):
        bx, h = x + 150 + k * 16, abs(v) * 50
        if v >= 0:
            b.append(rect(bx, base - h, 11, h, rx=2, fill="var(--sf-accent2)" if k < 3 else "var(--sf-line)", stroke="none", cls="sf-pulse" if k < 3 else ""))
        else:
            b.append(rect(bx, base, 11, h + 2, rx=2, fill="var(--sf-rose)", stroke="none"))
    b.append(path(f"M{x+146},{base} L{x+248},{base}", sw=1))
    b.append(text(x + 130, y + 140, "kept λ ≤ 0 → retry", size=10, fill="var(--sf-rose)"))

    # 6 · BFGS on the bounds-violation energy
    x, y = cols[2], rows[1]
    x0, y0 = x + 30, y + 120
    b.append(path(f"M{x0},{y0-80} L{x0},{y0} L{x0+200},{y0}", sw=1.1))
    b.append(polyline([(x0 + i * 10, y0 - 8 - 70 * 0.88 ** (i * 1.6) + (3 if i % 3 == 1 else 0)) for i in range(21)], color="var(--sf-warm)", sw=2, cls=""))
    b.append(path(f"M{x0},{y0-14} L{x0+200},{y0-14}", stroke="var(--sf-rose)", sw=1.1, dash="4 4"))
    b.append(text(x0 + 198, y0 - 19, "0.05 per atom", size=9, fill="var(--sf-rose)", anchor="end"))
    b.append(text(x0 + 196, y0 - 64, "E_bounds", size=9, fill="var(--sf-muted)", anchor="end"))
    b.append(text(x + 130, y + 140, "minimize squared bound violations", **muted))

    # 7 · Torsion potential from CSD preferences
    x, y = cols[0], rows[2]
    x0, y0 = x + 30, y + 86
    pts = [(x0 + i * 2, y0 + 24 - 3.2 * 7 * (1 + math.cos(3 * math.radians(i * 3.6)))) for i in range(101)]  # 7(1 + cos 3φ)
    b.append(polyline(pts, color="var(--sf-accent2)", sw=2, cls=""))
    b.append(path(f"M{x0},{y0+30} L{x0+200},{y0+30}", sw=1))
    for a_, lbl in [(0, "0°"), (100, "180°"), (200, "360°")]:
        b.append(text(x0 + a_, y0 + 41, lbl, size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(x + 130, y + 58, "V(φ) = Σₖ Vₖ (1 + sₖ cos kφ)", size=10, weight=600))
    b.append(text(x + 130, y + 142, "837 CSD SMARTS + UFF inversion", size=10, fill="var(--sf-muted)", weight=450))

    # 8 · Geometry checks
    x, y = cols[1], rows[2]
    for k, s_ in enumerate(["tetrahedral volume", "chiral volume sign", "sp² planarity", "double-bond geometry"]):
        yy = y + 58 + k * 20
        b.append(circle(x + 40, yy - 4, 7, fill="var(--sf-accent2)", cls=f"sf-pulse sf-d{k+1}"))
        b.append(path(f"M{x+36.5},{yy-4} l2.5,2.5 l4.5,-5", stroke="#fff"))
        b.append(text(x + 56, yy, s_, size=10.5, weight=500, anchor="start"))
    b.append(text(x + 130, y + 140, "any failure → retry from step 4", size=10, fill="var(--sf-rose)"))

    # 9 · Final 3D conformer
    x, y = cols[2], rows[2]
    cx, cy = x + 130, y + 90
    heavy = [(-48, 10, "C"), (0, -12, "C"), (46, 12, "O")]
    b.append('<g class="sf-float">')
    for hx, hy, p in [(-72, -6, 0), (-62, 36, 0), (-34, 32, 0), (-8, -40, 1), (18, -32, 1), (72, 0, 2)]:
        b.append(path(f"M{cx+heavy[p][0]},{cy+heavy[p][1]} L{cx+hx},{cy+hy}", sw=1.6))
        b.append(atom_ball(slug, cx + hx, cy + hy, "H", 5.5))
    b.append(path(f"M{cx-48},{cy+10} L{cx},{cy-12} L{cx+46},{cy+12}", stroke="var(--sf-ink)", sw=3))
    for ax, ay, el in heavy:
        b.append(atom_ball(slug, cx + ax, cy + ay, el))
    b.append("</g>")
    b.append(text(x + 130, y + 140, "x, y, z for all 9 atoms", **muted))

    # Retry loop back into embedding
    y8 = rows[2] + H
    b.append(retry_path(slug, f"M460,{y8+4} L460,{y8+30} L18,{y8+30} L18,{rows[1]+H/2} L34,{rows[1]+H/2}"))
    for c in (1, 2):
        b.append(path(f"M{cols[c]+130},{rows[1]+H+2} L{cols[c]+130},{rows[1]+H+14}", stroke="var(--sf-rose)", sw=1.2, dash="2 3"))
    b.append(text(472, y8 + 24, "retry loop · up to 10N attempts", size=10.5, fill="var(--sf-rose)", anchor="start"))
    b.append(text(472, y8 + 38, "after N/4 consecutive failures → random box placement in [−5, 5]³", size=9.5, fill="var(--sf-rose)", weight=500, anchor="start"))
    write(slug, 920, 730, "".join(b), label="ETKDGv2 nine-step pipeline in three phases, illustrated with ethanol: parse SMILES, molecular graph, bounds matrix with triangle smoothing, distance sampling, metric-matrix eigendecomposition, BFGS bounds minimization, ETKDG torsion force field, validation checks and 3D output, with the retry loop back to embedding.")


def diagram_eht_pipeline() -> None:
    slug = "eht-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    W, H, Y, gap, x0 = 196, 236, 64, 18, 24
    cards = [
        ("B1", "Atomic basis", "var(--sf-accent)"),
        ("B2", "Overlap S + Hamiltonian H", "var(--sf-accent2)"),
        ("B3", "Solve HC = SCE", "var(--sf-violet)"),
        ("B4", "MO on 3D grid", "var(--sf-warm)"),
        ("B5", "Marching Cubes", "var(--sf-rose)"),
    ]
    b.append(text(550, 32, "Extended Hückel pipeline — water (6 valence AOs, 8 electrons)", size=15, weight=700))
    xs = stage_row(b, x0, Y, W, H, gap, cards)
    cap = CAPTION

    # B1 · STO (exact) vs its 3-Gaussian fit
    x = xs[0]; ox, oy, pw, ph = x + 22, Y + 150, 160, 92
    b.append(path(f"M{ox},{oy-ph} L{ox},{oy} L{ox+pw},{oy}", sw=1.1))
    # Real STO-3G coefficients for a ζ = 1 1s Slater orbital (r in bohr, 0..4).
    al, dk = (2.227660, 0.405771, 0.109818), (0.154329, 0.535328, 0.444635)
    def g3(r):
        return sum(d * (2 * a_ / math.pi) ** 0.75 * math.exp(-a_ * r * r) for a_, d in zip(al, dk))
    s0 = 1 / math.sqrt(math.pi)
    sto = [(ox + i * pw / 60, oy - ph * 0.9 * math.exp(-i / 15)) for i in range(61)]
    gfit = [(ox + i * pw / 60, oy - ph * 0.9 * g3(i / 15) / s0) for i in range(61)]
    b.append(polyline(sto, color="var(--sf-accent)", sw=2.2, cls=""))
    b.append(polyline(gfit, color="var(--sf-warm)", sw=1.6, cls="sf-flow"))
    b.append(text(ox + pw, oy + 13, "r", size=9, fill="var(--sf-muted)", anchor="end"))
    b.append(path(f"M{x+96},{Y+64} l14,0", stroke="var(--sf-accent)", sw=2.2))
    b.append(text(x + 114, Y + 67, "Slater e^−ζr", size=9, anchor="start", weight=500))
    b.append(path(f"M{x+96},{Y+80} l14,0", stroke="var(--sf-warm)", sw=1.6, dash="4 6"))
    b.append(text(x + 114, Y + 83, "Σ₃ dₖ e^−αₖr²", size=9, anchor="start", weight=500))
    b.append(text(x + W / 2, Y + 176, "STO-3G per valence AO", size=10, weight=600))
    b.append(text(x + W / 2, Y + 194, "O: 2s 2pₓ 2p_y 2p_z · H: 1s ×2", **cap))
    b.append(text(x + W / 2, Y + 212, "VSIP Hᵢᵢ: O2p −14.8 · H1s −13.6 eV", **cap))

    # B2 · 6×6 matrix with VSIP diagonal and WH off-diagonals
    x = xs[1]; mx, my, cs = x + 26, Y + 50, 19
    labs = ["s", "x", "y", "z", "H", "H"]
    vals = [[1, 0, 0, 0, .52, .52], [0, 1, 0, 0, .38, -.38], [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, .30, .30], [.52, .38, 0, .30, 1, .18], [.52, -.38, 0, .30, .18, 1]]
    for i in range(6):
        b.append(text(mx - 8, my + i * cs + 13, labs[i], size=8.5, fill="var(--sf-muted)"))
        b.append(text(mx + i * cs + 9, my - 4, labs[i], size=8.5, fill="var(--sf-muted)"))
        for j in range(6):
            v = vals[i][j]
            if i == j:
                fill, op = "var(--sf-violet)", 0.85
            elif v == 0:
                fill, op = "var(--sf-card2)", None
            else:
                fill, op = "var(--sf-accent2)" if v > 0 else "var(--sf-rose)", round(0.25 + abs(v), 2)
            b.append(rect(mx + j * cs, my + i * cs, cs - 2, cs - 2, rx=2.5, fill=fill, sw=0.5, opacity=op))
    b.append(rect(mx - 2, my - 2, 6 * cs + 2, 6 * cs + 2, rx=4, fill="none", stroke="var(--sf-accent2)", sw=1.2, cls="sf-blink"))
    lx = x + 150
    for k, (lbl, c_) in enumerate([("Hᵢᵢ", "var(--sf-violet)"), ("Sᵢⱼ>0", "var(--sf-accent2)"), ("Sᵢⱼ<0", "var(--sf-rose)")]):
        b.append(rect(lx, Y + 56 + k * 18, 9, 9, rx=2, fill=c_, stroke="none", opacity=0.8))
        b.append(text(lx + 13, Y + 64 + k * 18, lbl, size=8.5, anchor="start", weight=500))
    b.append(text(x + W / 2, Y + 186, "Hᵢⱼ = ½ K Sᵢⱼ (Hᵢᵢ + Hⱼⱼ)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 204, "Wolfsberg–Helmholtz, K = 1.75", **cap))
    b.append(text(x + W / 2, Y + 220, "zero blocks = symmetry-forbidden", **cap))

    # B3 · Löwdin + MO level diagram with occupation
    x = xs[2]
    b.append(text(x + W / 2, Y + 54, "H′ = S⁻¹ᐟ² H S⁻¹ᐟ²,  C = S⁻¹ᐟ² C′", size=9.5, weight=600))
    levels = [(e, 2 if i < 4 else 0) for i, e in enumerate(WATER_EHT)]
    slots = [Y + 176, Y + 150, Y + 132, Y + 114, Y + 84, Y + 68]  # schematic spacing, energies labelled
    for k, (e, occ) in enumerate(levels):
        yy = slots[k]
        homo, lumo = k == 3, k == 4
        col = "var(--sf-violet)" if occ else "var(--sf-line)"
        b.append(path(f"M{x+60},{yy} L{x+110},{yy}", stroke=col, sw=2.2))
        if occ:
            b.append(path(f"M{x+78},{yy+5} l0,-11 l-2.5,3.5 M{x+92},{yy-6} l0,11 l2.5,-3.5", stroke="var(--sf-ink)", sw=1.2))
        b.append(text(x + 54, yy + 3, f"{e:.1f}", size=8, fill="var(--sf-muted)", anchor="end", weight=500))
        if homo or lumo:
            b.append(text(x + 116, yy + 3, "HOMO" if homo else "LUMO", size=9, fill="var(--sf-violet)" if homo else "var(--sf-rose)", anchor="start", weight=700))
    yh, yl = slots[3], slots[4]
    b.append(path(f"M{x+160},{yh} L{x+160},{yl}", stroke="var(--sf-warm)", sw=1.4, dash="3 3", cls="sf-flow"))
    b.append(text(x + 166, (yh + yl) / 2 + 3, "gap", size=9, fill="var(--sf-warm)", anchor="start"))
    b.append(text(x + W / 2, Y + 204, "Löwdin orthogonalization", size=10, weight=600))
    b.append(text(x + W / 2, Y + 220, "eV · 8 e⁻ fill 4 lowest MOs", **cap))

    # B4 · voxel grid sampling ψ(r)
    x = xs[3]; gx, gy, n, c = x + 30, Y + 50, 8, 17
    for i in range(n):
        for j in range(n):
            dx, dy = i - 3.5, j - 3.5
            v = (math.exp(-((dx + 1.6) ** 2 + dy ** 2) / 3) - math.exp(-((dx - 1.6) ** 2 + dy ** 2) / 3))
            fill = "var(--sf-accent)" if v > 0 else "var(--sf-rose)"
            op = round(min(0.9, abs(v) * 1.1), 2)
            b.append(rect(gx + i * c, gy + j * c, c - 1.5, c - 1.5, rx=2, fill=fill if op > 0.06 else "var(--sf-card2)", stroke="none", opacity=op if op > 0.06 else None))
    b.append(rect(gx - 2, gy - 2, n * c + 2, n * c + 2, rx=4, fill="none", sw=1))
    b.append(rect(gx + 2 * c, gy, c - 1.5, n * c - 1.5, rx=2, fill="none", stroke="var(--sf-warm)", sw=1.4, cls="sf-pulse"))
    b.append(text(x + W / 2, Y + 204, "ψₘ(r) = Σᵢ Cᵢₘ χᵢ(r)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 220, "spacing 0.15 Å · Rayon chunks", **cap))

    # B5 · isosurface lobes around water
    x = xs[4]; cx, cy = x + W / 2, Y + 112
    b.append('<g class="sf-float">')
    b.append(path(f"M{cx-6},{cy-4} C{cx-70},{cy-60} {cx-86},{cy+10} {cx-30},{cy+24} C{cx-14},{cy+28} {cx-4},{cy+12} {cx-6},{cy-4} Z", fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1.2, opacity=0.38))
    b.append(path(f"M{cx+6},{cy-4} C{cx+70},{cy-60} {cx+86},{cy+10} {cx+30},{cy+24} C{cx+14},{cy+28} {cx+4},{cy+12} {cx+6},{cy-4} Z", fill="var(--sf-rose)", stroke="var(--sf-rose)", sw=1.2, opacity=0.38))
    for (hx, hy) in [(cx - 34, cy + 22), (cx + 34, cy + 22)]:
        b.append(path(f"M{cx},{cy} L{hx},{hy}", stroke="var(--sf-ink)", sw=2.4))
        b.append(atom_ball(slug, hx, hy, "H", 6))
    b.append(atom_ball(slug, cx, cy, "O", 10))
    b.append("</g>")
    b.append(text(x + 44, Y + 66, "+ψ", size=10, fill="var(--sf-accent)", weight=700))
    b.append(text(x + W - 44, Y + 66, "−ψ", size=10, fill="var(--sf-rose)", weight=700))
    b.append(text(x + W / 2, Y + 186, "|ψ| = isovalue (e.g. 0.05)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 204, "256-case voxel table →", **cap))
    b.append(text(x + W / 2, Y + 220, "triangle mesh, two phases", **cap))
    write(slug, 1100, 324, "".join(b), label="Extended Hückel pipeline for water: STO-3G basis, overlap and Wolfsberg-Helmholtz Hamiltonian, Löwdin-orthogonalized eigenproblem with HOMO and LUMO, orbital sampled on a 3D grid, and Marching Cubes isosurface with two phases.")


# Butane carbon skeleton (C1–C4, gauche, φ = 65°) used as the worked example
# in the distance-geometry diagrams. All numbers below come from this geometry.
BUTANE_PAIRS = [  # label, lower, upper, picked distance (Å) — UFF sp³ C–C r = 1.514 Å
    ("1-2 bond", 1.50, 1.52, 1.51),
    ("1-3 angle", 2.43, 2.51, 2.47),
    ("1-4 torsion", 2.52, 3.81, 2.95),
]
BUTANE_D = [[0, 1.51, 2.47, 2.95], [1.51, 0, 1.51, 2.47], [2.47, 1.51, 0, 1.51], [2.95, 2.47, 1.51, 0]]
BUTANE_T = [[2.54, 0.59, -1.32, -1.82], [0.59, 0.94, -0.21, -1.32], [-1.32, -0.21, 0.94, 0.59], [-1.82, -1.32, 0.59, 2.54]]
BUTANE_EIG = [5.249, 1.449, 0.257, 0.000]
BUTANE_X = [(-1.47, -0.60, -0.15), (-0.68, 0.60, 0.32), (0.68, 0.60, -0.32), (1.47, -0.60, 0.15)]


def butane_glyph(slug: str, cx: float, cy: float, scale: float = 34) -> str:
    pts = [(cx + x * scale, cy - y * scale, z) for x, y, z in BUTANE_X]
    out = path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y, _ in pts), stroke="var(--sf-ink)", sw=2.6)
    for i, (x, y, z) in enumerate(pts):
        out += atom_ball(slug, x, y, "C", 9 + 3 * z)
        out += text(x, y + 22 if i in (0, 3) else y - 15, f"C{i+1}", size=9, fill="var(--sf-muted)")
    return out


def interval_bars(b: list, x: float, y: float, w: float, *, picked: bool, lo_ax: float = 1.2, hi_ax: float = 4.2) -> None:
    sx = lambda v: x + (v - lo_ax) / (hi_ax - lo_ax) * w
    for k, (lbl, lo, hi, d) in enumerate(BUTANE_PAIRS):
        yy = y + k * 30
        b.append(text(x, yy - 7, lbl, size=8.5, fill="var(--sf-muted)", anchor="start", weight=500))
        b.append(path(f"M{x},{yy} L{x+w},{yy}", sw=0.8))
        b.append(rect(sx(lo), yy - 5, max(sx(hi) - sx(lo), 3), 10, rx=5, fill="var(--sf-violet)", stroke="none", opacity=0.35))
        if picked:
            b.append(circle(sx(d), yy, 4.5, fill="var(--sf-warm)", cls=f"sf-pulse sf-d{k+1}"))
            b.append(text(x + w, yy - 7, f"d = {d:.2f}", size=8.5, fill="var(--sf-warm)", anchor="end"))
        else:
            b.append(text(x + w, yy - 7, f"{lo:.2f}–{hi:.2f}", size=8.5, fill="var(--sf-violet)", anchor="end"))
    yb = y + 3 * 30 - 6
    for v in (1.5, 2.5, 3.5):
        b.append(text(sx(v), yb, f"{v} Å", size=8, fill="var(--sf-muted)", weight=500))


def diagram_distance_geometry_pipeline() -> None:
    slug = "distance-geometry-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    W, H, Y = 196, 236, 64
    cards = [
        ("1", "Bounds lᵢⱼ ≤ dᵢⱼ ≤ uᵢⱼ", "var(--sf-accent)"),
        ("2", "Triangle smoothing", "var(--sf-accent2)"),
        ("3", "Pick distances", "var(--sf-violet)"),
        ("4", "Embed", "var(--sf-warm)"),
        ("5", "Refine", "var(--sf-rose)"),
    ]
    b.append(text(550, 32, "Distance geometry — butane carbon skeleton (C1–C2–C3–C4)", size=15, weight=700))
    xs = stage_row(b, 24, Y, W, H, 18, cards)

    # 1 · bounds per topological distance
    x = xs[0]
    interval_bars(b, x + 16, Y + 70, 164, picked=False)
    b.append(text(x + W / 2, Y + 186, "from topology only", size=10, weight=600))
    b.append(text(x + W / 2, Y + 204, "1-4 spans cis (2.52) → trans (3.81)", **CAPTION))
    b.append(text(x + W / 2, Y + 220, "≥ 1-5: VDW lower bounds", **CAPTION))

    # 2 · triangle inequality on bounds
    x = xs[1]; P = {"i": (x + 40, Y + 140), "j": (x + 156, Y + 140), "k": (x + 98, Y + 64)}
    for a_, c_, col in [("i", "k", "var(--sf-line)"), ("k", "j", "var(--sf-line)"), ("i", "j", "var(--sf-accent2)")]:
        b.append(path(f"M{P[a_][0]},{P[a_][1]} L{P[c_][0]},{P[c_][1]}", stroke=col, sw=2 if a_ == "i" and c_ == "j" else 1.6, cls="sf-flow" if c_ == "j" and a_ == "i" else ""))
    for n_, (px, py) in P.items():
        b.append(circle(px, py, 9, fill="var(--sf-card2)", stroke="var(--sf-ink)", sw=1.2))
        b.append(text(px, py + 4, n_, size=10, weight=700))
    b.append(text(x + 52, Y + 98, "uᵢₖ", size=9, fill="var(--sf-muted)"))
    b.append(text(x + 146, Y + 98, "uₖⱼ", size=9, fill="var(--sf-muted)"))
    b.append(text(x + 98, Y + 158, "uᵢⱼ", size=9, fill="var(--sf-accent2)"))
    b.append(text(x + W / 2, Y + 184, "uᵢⱼ ← min(uᵢⱼ, uᵢₖ + uₖⱼ)", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 200, "lᵢⱼ ← max(lᵢⱼ, lᵢₖ − uₖⱼ, lₖⱼ − uᵢₖ)", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 220, "Floyd–Warshall · lᵢⱼ > uᵢⱼ ⇒ infeasible", **CAPTION))

    # 3 · random distance inside each interval
    x = xs[2]
    interval_bars(b, x + 16, Y + 70, 164, picked=True)
    b.append(text(x + W / 2, Y + 186, "dᵢⱼ = lᵢⱼ + r (uᵢⱼ − lᵢⱼ)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 204, "r ~ U(0,1), MinstdRand", **CAPTION))
    b.append(text(x + W / 2, Y + 220, "sₙ₊₁ = 48271 sₙ mod (2³¹−1)", **CAPTION))

    # 4 · embed (metric matrix → eigenvectors → coordinates)
    x = xs[3]
    b.append('<g class="sf-float">' + butane_glyph(slug, x + W / 2, Y + 110, 40) + "</g>")
    b.append(text(x + W / 2, Y + 186, "T = VΛVᵀ,  xᵢₖ = √λₖ vᵢₖ", size=10, weight=600))
    b.append(text(x + W / 2, Y + 204, "λ = 5.25, 1.45, 0.26, 0.00", **CAPTION))
    b.append(text(x + W / 2, Y + 220, "4D when chiral centers exist", **CAPTION))

    # 5 · refine with the bounds force field
    x = xs[4]; ox, oy = x + 26, Y + 150
    b.append(path(f"M{ox},{oy-84} L{ox},{oy} L{ox+150},{oy}", sw=1.1))
    b.append(polyline([(ox + i * 7.5, oy - 6 - 76 * math.exp(-i / 4.2)) for i in range(21)], color="var(--sf-rose)", sw=2, cls=""))
    b.append(text(ox + 6, oy - 72, "E_bounds", size=9, fill="var(--sf-muted)", anchor="start"))
    b.append(text(ox + 150, oy + 13, "BFGS steps", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    b.append(text(x + W / 2, Y + 186, "Σ (d² / u² − 1)² + chiral + 4D", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 204, "phase 1: w_chiral 1.0 · w_4D 0.1", **CAPTION))
    b.append(text(x + W / 2, Y + 220, "phase 2: w_chiral 0.2 · w_4D 1.0", **CAPTION))
    write(slug, 1100, 324, "".join(b), label="Distance geometry on butane: topological bounds, triangle-inequality smoothing, random distance picking, metric-matrix embedding and bounds force-field refinement.")


def number_matrix(b: list, x: float, y: float, M: list, cell: float, *, diverging: bool, color: str, vmax: float, labels: list[str]) -> None:
    n = len(M)
    for i in range(n):
        b.append(text(x - 10, y + i * cell + cell / 2 + 3, labels[i], size=8.5, fill="var(--sf-muted)"))
        b.append(text(x + i * cell + cell / 2, y - 6, labels[i], size=8.5, fill="var(--sf-muted)"))
        for j in range(n):
            v = M[i][j]
            fill = color if (v >= 0 or not diverging) else "var(--sf-rose)"
            op = round(0.12 + 0.6 * abs(v) / vmax, 2)
            b.append(rect(x + j * cell, y + i * cell, cell - 2, cell - 2, rx=3, fill=fill, stroke="none", opacity=op))
            b.append(text(x + j * cell + cell / 2 - 1, y + i * cell + cell / 2 + 3, f"{v:.2f}", size=8.5, weight=600, mono=True))


def diagram_distance_geometry_matrix() -> None:
    slug = "distance-geometry-matrix"
    begin(slug)
    b = [element_defs(slug)]
    W, H, Y = 250, 250, 64
    cards = [
        ("D", "Squared distances dᵢⱼ²", "var(--sf-accent)"),
        ("T", "Metric (Gram) matrix", "var(--sf-accent2)"),
        ("λ", "Power iteration", "var(--sf-violet)"),
        ("X", "Coordinates", "var(--sf-warm)"),
    ]
    b.append(text(550, 32, "From distances to coordinates — butane carbons, real numbers", size=15, weight=700))
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    labs = ["C1", "C2", "C3", "C4"]
    D2 = [[round(v * v, 2) for v in row] for row in BUTANE_D]
    number_matrix(b, xs[0] + 46, Y + 58, D2, 40, diverging=False, color="var(--sf-accent)", vmax=9, labels=labs)
    b.append(text(xs[0] + W / 2, Y + 236, "d₁₄ = 2.95 Å (gauche pick, φ = 65°)", **CAPTION))
    number_matrix(b, xs[1] + 46, Y + 58, BUTANE_T, 40, diverging=True, color="var(--sf-accent2)", vmax=2.6, labels=labs)
    b.append(text(xs[1] + W / 2, Y + 230, "Tᵢⱼ = ½(D₀ᵢ + D₀ⱼ − Dᵢⱼ) = xᵢ·xⱼ", size=9.5, weight=600))
    b.append(text(xs[1] + W / 2, Y + 245, "negative = atoms on opposite sides", **CAPTION))
    # eigenvalues
    x = xs[2]; base = Y + 170
    for k, v in enumerate(BUTANE_EIG):
        bx = x + 40 + k * 46; h = v / BUTANE_EIG[0] * 100
        keep = k < 3
        b.append(rect(bx, base - max(h, 1.5), 30, max(h, 1.5), rx=3, fill="var(--sf-violet)" if keep else "var(--sf-line)", stroke="none", cls="sf-pulse" if keep else ""))
        b.append(text(bx + 15, base - h - 6, f"{v:.2f}", size=9, weight=650, fill="var(--sf-violet)" if keep else "var(--sf-muted)"))
        b.append(text(bx + 15, base + 14, f"λ{k+1}", size=9, fill="var(--sf-muted)"))
    b.append(path(f"M{x+34},{base} L{x+222},{base}", sw=1))
    b.append(text(x + W / 2, Y + 206, "v ← Tv / ‖Tv‖,  λ = vᵀTv", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 222, "deflate T ← T − λvvᵀ, repeat ×3", **CAPTION))
    b.append(text(x + W / 2, Y + 242, "λ₄ = 0 → exactly 3D-embeddable", size=9.5, fill="var(--sf-accent2)", weight=600))
    # coordinates
    x = xs[3]
    b.append('<g class="sf-float">' + butane_glyph(slug, x + W / 2, Y + 110, 50) + "</g>")
    for i, (cx_, cy_, cz_) in enumerate(BUTANE_X):
        b.append(text(x + 24, Y + 182 + i * 15, f"C{i+1}  {cx_:+.2f} {cy_:+.2f} {cz_:+.2f}", size=9, weight=500, anchor="start", mono=True))
    b.append(text(x + W - 14, Y + 182, "xᵢₖ =", size=9, fill="var(--sf-muted)", anchor="end"))
    b.append(text(x + W - 14, Y + 196, "√λₖ vᵢₖ", size=9, fill="var(--sf-muted)", anchor="end"))
    write(slug, 1100, 336, "".join(b), label="Worked butane example: squared distance matrix, metric matrix T, eigenvalues 5.25, 1.45, 0.26 and 0 from power iteration, and recovered 3D coordinates.")


def planar_14(r: float, theta_deg: float, phi_deg: float) -> float:
    """1-4 distance for equal bonds r, equal angles θ and dihedral φ (exact geometry)."""
    th, ph = math.radians(theta_deg), math.radians(phi_deg)
    p0 = (r * math.cos(th), r * math.sin(th), 0.0)
    p3 = (r - r * math.cos(th), r * math.sin(th) * math.cos(ph), r * math.sin(th) * math.sin(ph))
    return math.dist(p0, p3)


TOPO_COLOR = {1: "var(--sf-accent)", 2: "var(--sf-accent2)", 3: "var(--sf-violet)", 4: "var(--sf-warm)"}


def diagram_bounds_matrix_overview() -> None:
    slug = "bounds-matrix-overview"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 32, "Bounds matrix for pentane carbons (C1…C5), in Å", size=15, weight=700))
    U = {1: 1.52, 2: 2.51, 3: 3.81, 4: 5.02}
    L = {1: 1.50, 2: 2.43, 3: 2.52, 4: 2.38}
    x0, y0, cs, n = 70, 90, 54, 5
    b.append(rect(x0 - 14, y0 - 30, n * cs + 28, n * cs + 44, rx=14, fill="var(--sf-card)", sw=1.3))
    for i in range(n):
        b.append(text(x0 - 6, y0 + i * cs + cs / 2 + 3, f"C{i+1}", size=9, fill="var(--sf-muted)", anchor="end"))
        b.append(text(x0 + i * cs + cs / 2, y0 - 10, f"C{i+1}", size=9, fill="var(--sf-muted)"))
        for j in range(n):
            x, y = x0 + j * cs, y0 + i * cs
            if i == j:
                b.append(rect(x, y, cs - 3, cs - 3, rx=5, fill="var(--sf-card2)", sw=0.6))
                b.append(text(x + cs / 2 - 1, y + cs / 2 + 3, "0", size=10, fill="var(--sf-muted)", mono=True))
                continue
            t_ = abs(i - j)
            upper = j > i
            b.append(rect(x, y, cs - 3, cs - 3, rx=5, fill=TOPO_COLOR[t_], stroke="none", opacity=0.55 if upper else 0.25))
            b.append(text(x + cs / 2 - 1, y + cs / 2 + 4, f"{(U if upper else L)[t_]:.2f}", size=10.5, weight=650, mono=True))
    # triangle labels
    b.append(text(x0 + n * cs - 4, y0 + n * cs + 6, "", size=1))
    b.append(text(x0 + n * cs * 0.78, y0 + n * cs * 0.2 - 2, "", size=1))
    b.append(text(x0 + n * cs / 2, y0 + n * cs + 6, "upper triangle uᵢⱼ (i < j) · lower triangle lᵢⱼ · diagonal 0", size=9.5, fill="var(--sf-muted)", weight=500))
    # highlighted 1-4 pair (both triangles)
    for (i, j) in [(0, 3), (3, 0)]:
        b.append(rect(x0 + j * cs - 2, y0 + i * cs - 2, cs + 1, cs + 1, rx=6, fill="none", stroke="var(--sf-ink)", sw=2, cls="sf-pulse"))
    # legend by source
    lx, ly = 410, 84
    rows = [
        (1, "1-2 bond", "UFF rᵢⱼ ± 0.01", "1.514 ± 0.01"),
        (2, "1-3 angle", "law of cosines ± 0.04", "θ = 109.47°"),
        (3, "1-4 torsion", "[d_cis, d_trans]", "2.52 … 3.81"),
        (4, "1-5 VDW", "0.70 · (r_vdw + r_vdw)", "0.70 · 3.40 = 2.38"),
    ]
    b.append(text(lx, ly - 12, "where each entry comes from", size=11, weight=700, anchor="start"))
    for k, (t_, name, rule, ex) in enumerate(rows):
        y = ly + k * 46
        b.append(rect(lx, y, 400, 38, rx=10, fill="var(--sf-card)", sw=1))
        b.append(rect(lx + 10, y + 11, 16, 16, rx=4, fill=TOPO_COLOR[t_], stroke="none", opacity=0.75))
        b.append(text(lx + 36, y + 24, name, size=11, weight=700, anchor="start"))
        b.append(text(lx + 128, y + 24, rule, size=10, fill="var(--sf-muted)", weight=500, anchor="start"))
        b.append(text(lx + 390, y + 24, ex, size=10, fill=TOPO_COLOR[t_], anchor="end", mono=True))
    # smoothing callout for 1-5 upper bound
    y = ly + 4 * 46 + 6
    b.append(rect(lx, y, 400, 54, rx=10, fill="var(--sf-soft)", stroke="var(--sf-warm)", sw=1.2))
    b.append(text(lx + 14, y + 22, "1-5 upper bound before smoothing: 10⁶ Å", size=10.5, weight=600, anchor="start"))
    b.append(text(lx + 14, y + 40, "after Floyd–Warshall: min(u₁₃ + u₃₅) = 2.51 + 2.51 = 5.02", size=10, fill="var(--sf-warm)", anchor="start"))
    # 1-4 interval
    y = ly + 5 * 46 + 30
    b.append(text(lx, y, "highlighted pair C1–C4:", size=10, weight=600, anchor="start"))
    sx = lambda v: lx + 150 + (v - 2.2) / (4.0 - 2.2) * 240
    b.append(path(f"M{sx(2.2)},{y-4} L{sx(4.0)},{y-4}", sw=1))
    b.append(rect(sx(2.52), y - 10, sx(3.81) - sx(2.52), 12, rx=6, fill="var(--sf-violet)", stroke="none", opacity=0.35))
    b.append(text(sx(2.52), y + 14, "cis 2.52", size=9, fill="var(--sf-violet)"))
    b.append(text(sx(3.81), y + 14, "trans 3.81", size=9, fill="var(--sf-violet)"))
    write(slug, 900, 400, "".join(b), label="Pentane carbon bounds matrix with real values: upper bounds above the diagonal, lower bounds below, colored by source (1-2 bond, 1-3 angle, 1-4 torsion, 1-5 van der Waals), and the 1-5 upper bound tightened by triangle smoothing.")


def diagram_bounds_matrix_triangle() -> None:
    slug = "bounds-matrix-triangle"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 32, "1-3 bounds from the law of cosines — ethanol C–C–O", size=15, weight=700))
    # main triangle drawn to scale (60 px per Å)
    k = 120
    C2 = (230, 230)
    th = math.radians(109.47)
    C1 = (C2[0] - 1.514 * k, C2[1])
    O = (C2[0] - 1.407 * k * math.cos(th), C2[1] - 1.407 * k * math.sin(th))
    b.append(path(f"M{C1[0]:.1f},{C1[1]:.1f} L{C2[0]},{C2[1]} L{O[0]:.1f},{O[1]:.1f}", stroke="var(--sf-ink)", sw=3))
    b.append(path(f"M{C1[0]:.1f},{C1[1]:.1f} L{O[0]:.1f},{O[1]:.1f}", stroke="var(--sf-warm)", sw=2, dash="6 5", cls="sf-flow"))
    # angle arc
    r_ = 30
    a1, a2 = math.pi, math.pi - th
    b.append(path(f"M{C2[0]-r_},{C2[1]} A{r_},{r_} 0 0 1 {C2[0]+r_*math.cos(a2):.1f},{C2[1]-r_*math.sin(a2):.1f}", stroke="var(--sf-violet)", sw=1.6))
    b.append(text(C2[0] - 6, C2[1] - 38, "θ = 109.47°", size=10.5, fill="var(--sf-violet)", anchor="end"))
    b.append(text((C1[0] + C2[0]) / 2, C2[1] + 24, "d₁ = 1.51 Å (C–C)", size=10, fill="var(--sf-accent)"))
    b.append(text((C2[0] + O[0]) / 2 + 12, (C2[1] + O[1]) / 2 + 4, "d₂ = 1.41 Å (C–O)", size=10, fill="var(--sf-accent)", anchor="start"))
    b.append(text((C1[0] + O[0]) / 2 - 10, (C1[1] + O[1]) / 2 - 8, "d₁₃ = 2.39 Å", size=11, fill="var(--sf-warm)", anchor="end", weight=700))
    b.append(atom_ball(slug, *C1, "C", 12))
    b.append(atom_ball(slug, *C2, "C", 12))
    b.append(atom_ball(slug, O[0], O[1], "O", 12))
    # formula card
    x = 470
    b.append(rect(x, 64, 400, 92, rx=12, sw=1.2))
    b.append(text(x + 200, 92, "d₁₃ = √(d₁² + d₂² − 2 d₁ d₂ cos θ)", size=12.5, weight=700))
    b.append(text(x + 200, 116, "= √(1.514² + 1.407² − 2·1.514·1.407·cos 109.47°)", size=10, fill="var(--sf-muted)", weight=500, mono=True))
    b.append(text(x + 200, 140, "bounds [2.35, 2.43] Å  (± 0.04)", size=11, fill="var(--sf-warm)", weight=700))
    # hybridization strip
    b.append(text(x, 184, "θ by hybridization of the middle atom (C–C–C, UFF radii):", size=10.5, weight=600, anchor="start"))
    hyb = [("sp³", 109.47, 1.514, "2.47"), ("sp²", 120.0, 1.464, "2.54"), ("sp", 180.0, 1.412, "2.82")]
    for i, (name, ang, r, d) in enumerate(hyb):
        cx, cy = x + 66 + i * 134, 248
        L = 34
        t = math.radians(ang)
        pA = (cx - L / 2 - L * 0.0, cy)
        mid = (cx, cy)
        pA = (cx - L, cy) if ang == 180 else (cx - L * math.sin(t / 2), cy + L * math.cos(t / 2) * 0.0 + L * math.cos(t / 2))
        pB = (cx + L, cy) if ang == 180 else (cx + L * math.sin(t / 2), cy + L * math.cos(t / 2))
        b.append(rect(cx - 60, cy - 36, 120, 106, rx=10, fill="var(--sf-card)", sw=1))
        b.append(path(f"M{pA[0]:.1f},{pA[1]:.1f} L{mid[0]},{mid[1]} L{pB[0]:.1f},{pB[1]:.1f}", stroke="var(--sf-ink)", sw=2.2))
        b.append(path(f"M{pA[0]:.1f},{pA[1]:.1f} L{pB[0]:.1f},{pB[1]:.1f}", stroke="var(--sf-warm)", sw=1.2, dash="3 3"))
        for q in (pA, mid, pB):
            b.append(atom_ball(slug, q[0], q[1], "C", 6, label=False))
        b.append(text(cx, cy - 18, f"{name} · {ang:g}°", size=10, weight=700))
        b.append(text(cx, cy + 56, f"d₁₃ = {d} Å", size=10, fill="var(--sf-warm)", mono=True))
    b.append(text(450, 346, "multiple paths → union: lowest lower bound, highest upper bound · ±0.04 doubled per large ring sp² atom", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 360, "".join(b), label="Law of cosines for 1-3 bounds in ethanol: C–C 1.51 Å, C–O 1.41 Å and 109.47° give 2.39 Å with ±0.04 Å tolerance; reference 1-3 distances for sp3, sp2 and sp middle atoms.")


def diagram_bounds_matrix_1_4() -> None:
    slug = "bounds-matrix-1-4"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 32, "1-4 bounds come from the cis and trans extremes — butane C1…C4", size=15, weight=700))
    r, th = 1.514, 109.47
    k = 48
    def chain(ox, oy, trans):
        t = math.radians(th)
        p1 = (0.0, 0.0); p2 = (r, 0.0)
        p0 = (r * math.cos(t), r * math.sin(t))
        p3 = (r - r * math.cos(t), (-1 if trans else 1) * r * math.sin(t))
        pts = [p0, p1, p2, p3]
        return [(ox + x * k, oy - y * k) for x, y in pts]
    for idx, (trans, label, d, col) in enumerate([(False, "cis · φ = 0°", 2.52, "var(--sf-rose)"), (True, "trans · φ = 180°", 3.81, "var(--sf-accent2)")]):
        bx = 40 + idx * 230
        b.append(rect(bx, 60, 210, 230, rx=14, sw=1.3))
        oy = 175 if trans else 215
        pts = chain(bx + 68, oy, trans)
        b.append(path(f"M{pts[0][0]:.1f},{pts[0][1]:.1f} L{pts[3][0]:.1f},{pts[3][1]:.1f}", stroke=col, sw=2, dash="5 4", cls="sf-flow"))
        b.append(path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts), stroke="var(--sf-ink)", sw=3))
        for i, (x, y) in enumerate(pts):
            b.append(atom_ball(slug, x, y, "C", 10))
            b.append(text(x, y + (24 if y >= oy else -16), f"C{i+1}", size=9, fill="var(--sf-muted)"))
        b.append(text(bx + 105, 84, label, size=11.5, weight=700))
        b.append(text(bx + 105, 270, f"d₁₄ = {d:.2f} Å", size=12, fill=col, weight=700))
    # d14(φ) curve
    px, py, pw, ph = 520, 80, 340, 170
    lo, hi = 2.3, 4.0
    sy = lambda v: py + ph - (v - lo) / (hi - lo) * ph
    sx = lambda a: px + a / 360 * pw
    b.append(rect(sx(0), sy(3.81), pw, sy(2.52) - sy(3.81), rx=0, fill="var(--sf-violet)", stroke="none", opacity=0.14))
    b.append(axes(px, py, pw, ph, xlabel="dihedral φ (°)", ylabel="d₁₄ (Å)"))
    b.append(curve_path([(sx(a), sy(planar_14(r, th, a))) for a in range(0, 361, 10)], color="var(--sf-violet)", sw=2.2, cls=""))
    for v, lbl, col in [(2.52, "l₁₄ = d_cis", "var(--sf-rose)"), (3.81, "u₁₄ = d_trans", "var(--sf-accent2)")]:
        b.append(path(f"M{px},{sy(v):.1f} L{px+pw},{sy(v):.1f}", stroke=col, sw=1.2, dash="4 4"))
        b.append(text(px + pw - 4, sy(v) - 5, f"{lbl} = {v}", size=9.5, fill=col, anchor="end"))
    for a in (0, 60, 180, 300, 360):
        b.append(text(sx(a), py + ph + 13, f"{a}", size=8.5, fill="var(--sf-muted)", weight=500))
    g = planar_14(r, th, 65)
    b.append(circle(sx(65), sy(g), 5, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(text(sx(65) + 8, sy(g) + 14, f"gauche 65° → {g:.2f}", size=9, fill="var(--sf-warm)", anchor="start"))
    b.append(text(690, 296, "default: [min, max] · E/Z forced · ring sp²–sp² planar → [d_cis, d_cis]", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 312, "".join(b), label="Butane cis (2.52 Å) and trans (3.81 Å) conformations set the 1-4 lower and upper bounds; every dihedral gives a distance inside that band.")


def power_iteration_trace(T: list, iters: int = 12, seed_vec: tuple = (0.6, -0.2, 0.5, 0.6)) -> tuple[list, list]:
    """Real power iteration with deflation on a small matrix. Returns (eigs, rayleigh traces)."""
    n = len(T)
    M = [row[:] for row in T]
    eigs, traces = [], []
    for k in range(3):
        v = list(seed_vec[k:] + seed_vec[:k])
        nv = math.sqrt(sum(x * x for x in v)); v = [x / nv for x in v]
        tr = []
        for _ in range(iters):
            w = [sum(M[i][j] * v[j] for j in range(n)) for i in range(n)]
            tr.append(sum(v[i] * w[i] for i in range(n)))
            nw = math.sqrt(sum(x * x for x in w)); v = [x / nw for x in w]
        lam = sum(v[i] * sum(M[i][j] * v[j] for j in range(n)) for i in range(n))
        eigs.append(lam); traces.append(tr)
        M = [[M[i][j] - lam * v[i] * v[j] for j in range(n)] for i in range(n)]
    return eigs, traces


def diagram_embedding_eigendecomp() -> None:
    slug = "embedding-eigendecomp"
    begin(slug)
    b = [element_defs(slug)]
    eigs, traces = power_iteration_trace(BUTANE_T)
    b.append(text(450, 32, "Power iteration with deflation on butane's metric matrix T", size=15, weight=700))
    # convergence plot
    px, py, pw, ph = 70, 76, 470, 220
    b.append(rect(px - 46, py - 20, pw + 70, ph + 64, rx=14, sw=1.3))
    b.append(axes(px, py, pw, ph, xlabel="iteration", ylabel="Rayleigh quotient vᵀTv"))
    vmax = 6.0
    sy = lambda v: py + ph - v / vmax * ph
    sx = lambda i: px + 10 + i * (pw - 20) / 11
    cols = ["var(--sf-violet)", "var(--sf-accent2)", "var(--sf-warm)"]
    for k, tr in enumerate(traces):
        b.append(path(f"M{px},{sy(eigs[k]):.1f} L{px+pw},{sy(eigs[k]):.1f}", stroke=cols[k], sw=1, dash="3 4", opacity=0.7))
        b.append(polyline([(sx(i), sy(v)) for i, v in enumerate(tr)], color=cols[k], sw=2, cls=""))
        for i, v in enumerate(tr):
            b.append(circle(sx(i), sy(v), 2.6, fill=cols[k]))
        b.append(text(px + pw - 4, sy(eigs[k]) - 6, f"λ{k+1} = {eigs[k]:.3f}", size=10, fill=cols[k], anchor="end"))
    for v in (0, 2, 4, 6):
        b.append(text(px - 8, sy(v) + 3, str(v), size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    # algorithm card
    x = 600
    b.append(rect(x, 56, 276, 150, rx=14, sw=1.3))
    b.append(text(x + 16, 80, "for k = 1, 2, 3:", size=10.5, weight=700, anchor="start", mono=True))
    lines = ["v ← random unit (Mt19937)", "repeat ≤ 200:", "  w = T v,  λ = v·w", "  v = w / ‖w‖  (tol 10⁻⁶)", "T ← T − λ v vᵀ   (deflate)"]
    for i, ln in enumerate(lines):
        b.append(text(x + 24, 102 + i * 20, ln, size=10, weight=500, anchor="start", mono=True, fill="var(--sf-ink)" if i != 4 else "var(--sf-rose)"))
    # result card
    b.append(rect(x, 218, 276, 98, rx=14, sw=1.3))
    b.append(text(x + 16, 240, "xᵢₖ = √λₖ · vᵢₖ   (check λₖ > 10⁻³)", size=10, weight=650, anchor="start"))
    b.append('<g class="sf-float">' + butane_glyph(slug, x + 138, 282, 26) + "</g>")
    b.append(text(450, 366, "Each deflation removes the found direction, so the next run converges to the next-largest eigenvalue; λ₄ = 0 is never needed.", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 380, "".join(b), label="Real power-iteration traces on butane's metric matrix: three deflated runs converge to 5.25, 1.45 and 0.26, which scale the eigenvectors into 3D coordinates.")


def torsion_v(phi_deg: float) -> float:
    """CSD pattern [!#1][CX4]!@;-[CX4][!#1]: V = 7 (1 + cos 3φ)."""
    return 7.0 * (1 + math.cos(3 * math.radians(phi_deg)))


def newman(cx: float, cy: float, phi: float, r: float = 18) -> str:
    """Newman projection: front C (dot) with 3 bonds, back C (circle) rotated by φ; methyls highlighted."""
    out = circle(cx, cy, r, fill="var(--sf-card)", stroke="var(--sf-ink)", sw=1.3)
    for k in range(3):
        a = math.radians(-90 + k * 120)
        col, sw_ = ("var(--sf-violet)", 2.4) if k == 0 else ("var(--sf-line)", 1.4)
        out += path(f"M{cx},{cy} L{cx + 1.6*r*math.cos(a):.1f},{cy + 1.6*r*math.sin(a):.1f}", stroke=col, sw=sw_)
        a2 = math.radians(-90 + phi + k * 120)
        col2, sw2 = ("var(--sf-violet)", 2.4) if k == 0 else ("var(--sf-line)", 1.4)
        out += path(f"M{cx + r*math.cos(a2):.1f},{cy + r*math.sin(a2):.1f} L{cx + 1.6*r*math.cos(a2):.1f},{cy + 1.6*r*math.sin(a2):.1f}", stroke=col2, sw=sw2)
    return out + circle(cx, cy, 2.6, fill="var(--sf-ink)")


def diagram_etkdg_torsion() -> None:
    slug = "etkdg-torsion"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Butane C1–C2–C3–C4 torsion: V(φ) = 7 (1 + cos 3φ)", size=15, weight=700))
    b.append(text(450, 50, "matched pattern [!#1:1][CX4:2]!@;-[CX4:3][!#1:4]  ·  s₃ = +1, V₃ = 7.0, other Vₖ = 0", size=10, fill="var(--sf-muted)", weight=500, mono=True))
    px, py, pw, ph = 80, 140, 740, 170
    sx = lambda a: px + (a + 180) / 360 * pw
    sy = lambda v: py + ph - v / 15 * ph
    b.append(axes(px, py, pw, ph, xlabel="dihedral φ (°)", ylabel="V (kcal/mol)"))
    curve = [(sx(a), sy(torsion_v(a))) for a in range(-180, 181, 5)]
    b.append(curve_path(curve, color="var(--sf-accent2)", sw=2.4, cls=""))
    for a in range(-180, 181, 60):
        b.append(text(sx(a), py + ph + 13, f"{a}", size=8.5, fill="var(--sf-muted)", weight=500))
    for v in (0, 7, 14):
        b.append(text(px - 8, sy(v) + 3, str(v), size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    # Newman projections above key angles
    for a, name, col in [(-180, "anti", "var(--sf-accent2)"), (-60, "gauche⁻", "var(--sf-accent2)"), (0, "syn (eclipsed)", "var(--sf-rose)"), (60, "gauche⁺", "var(--sf-accent2)"), (120, "eclipsed", "var(--sf-rose)"), (180, "anti", "var(--sf-accent2)")]:
        cx = min(max(sx(a), px + 30), px + pw - 30)
        b.append(newman(cx, 92, a))
        b.append(text(cx, 128 if a not in (0, 120) else 128, name, size=9, fill=col))
    # DG starting point rolling into nearest well
    start = 95
    d = f"M{sx(start):.1f},{sy(torsion_v(start)):.1f} " + " ".join(f"L{sx(a):.1f},{sy(torsion_v(a)):.1f}" for a in range(start, 59, -5))
    b.append(circle(sx(start), sy(torsion_v(start)), 5, fill="var(--sf-warm)"))
    b.append(text(sx(start) + 8, sy(torsion_v(start)) - 6, f"DG start φ = {start}°", size=9.5, fill="var(--sf-warm)", anchor="start"))
    b.append(particle(d, dur=2.4, r=4, fill="var(--sf-warm)"))
    b.append(path(d, stroke="var(--sf-warm)", sw=1.4, dash="3 4", cls="sf-flow"))
    b.append(text(450, 352, "BFGS slides each dihedral into its nearest well. With V₃ alone, anti and gauche wells are equally deep (0 kcal/mol);", size=9.5, fill="var(--sf-muted)", weight=500))
    b.append(text(450, 367, "eclipsed geometries (0°, ±120°) cost 14 kcal/mol and are effectively forbidden.", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 380, "".join(b), label="Butane torsion energy 7(1 + cos 3φ) from the matched CSD SMARTS pattern, with Newman projections at anti, gauche and eclipsed angles and a distance-geometry start relaxing into the nearest well.")


def dihedral_hist(b: list, x: float, y: float, w: float, h: float, counts: list, color: str) -> None:
    b.append(path(f"M{x},{y+h} L{x+w},{y+h}", sw=1))
    bw = w / len(counts)
    mx = max(counts)
    for i, c in enumerate(counts):
        hh = c / mx * h
        b.append(rect(x + i * bw + 1, y + h - hh, bw - 2, hh, rx=1.5, fill=color, stroke="none", opacity=0.75))
    for a, lbl in [(0, "−180"), (0.5, "0"), (1, "180")]:
        b.append(text(x + a * w, y + h + 12, lbl, size=8, fill="var(--sf-muted)", weight=500))


def diagram_etkdg_pipeline() -> None:
    slug = "etkdg-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "What ETKDG adds on top of plain distance geometry", size=15, weight=700))
    # left: plain DG histogram (near uniform)
    b.append(rect(24, 54, 220, 250, rx=14, sw=1.3))
    b.append(text(134, 78, "after DG embedding", size=11.5, weight=700))
    b.append(text(134, 94, "sp³–sp³ dihedrals are ~random", size=9.5, fill="var(--sf-muted)", weight=450))
    uni = [5, 6, 5, 7, 6, 5, 6, 7, 5, 6, 6, 5, 7, 6, 5, 6, 6, 5]
    dihedral_hist(b, 44, 120, 180, 120, uni, "var(--sf-line)")
    b.append(text(134, 280, "φ (°)", size=9, fill="var(--sf-muted)"))
    # middle: force field terms
    x = 290
    b.append(rect(x, 54, 330, 250, rx=14, sw=1.3))
    b.append(text(x + 165, 78, "ETKDG 3D force field → BFGS (≤ 300 iter)", size=11.5, weight=700))
    terms = [
        ("torsion", "Σ Vₖ (1 + sₖ cos kφ)", "837 CSD SMARTS, else basic rules", "var(--sf-accent2)"),
        ("inversion", "K (1 − sin Y) ×3 perms", "planar sp² centers", "var(--sf-violet)"),
        ("distance", "k/2 · max(0, |d − d₀| − ε)²", "1-2 k = 100, ε = 0.01 Å", "var(--sf-accent)"),
        ("linear angle", "k (θ − 180°)²", "sp atoms", "var(--sf-warm)"),
    ]
    for i, (name, f_, note, col) in enumerate(terms):
        y = 94 + i * 50
        b.append(rect(x + 14, y, 302, 42, rx=9, fill="var(--sf-card2)", stroke=col, sw=1.1))
        b.append(text(x + 26, y + 18, name, size=10.5, fill=col, weight=700, anchor="start"))
        b.append(text(x + 120, y + 18, f_, size=10, weight=600, anchor="start"))
        b.append(text(x + 120, y + 33, note, size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(arrow_between(248, 180, 286, 180, color="var(--sf-accent2)"))
    b.append(arrow_between(624, 180, 652, 180, color="var(--sf-accent2)", delay=0.4))
    # right: ETKDG histogram peaked at ±60, 180
    b.append(rect(656, 54, 220, 250, rx=14, sw=1.3))
    b.append(text(766, 78, "after ETKDG refinement", size=11.5, weight=700))
    b.append(text(766, 94, "dihedrals sit in CSD wells", size=9.5, fill="var(--sf-muted)", weight=450))
    peaked = [10, 3, 0, 1, 4, 9, 4, 1, 0, 0, 1, 4, 9, 4, 1, 0, 3, 10]
    dihedral_hist(b, 676, 120, 180, 120, peaked, "var(--sf-accent2)")
    for a, lbl in [(-60, "g⁻"), (60, "g⁺"), (180, "anti")]:
        b.append(text(676 + (a + 180) / 360 * 180, 114, lbl, size=8.5, fill="var(--sf-accent2)"))
    b.append(text(766, 280, "φ (°)", size=9, fill="var(--sf-muted)"))
    b.append(text(450, 326, "Illustrative histograms for sp³–sp³ bonds · skip refinement when initial E < 10⁻⁵ · result goes to validation", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 340, "".join(b), label="Plain distance geometry leaves sp3–sp3 dihedrals near random; the ETKDG force field of CSD torsions, inversions, distance and linear-angle constraints, minimized by BFGS, moves them into the experimental wells.")


def energy_panel(
    b: list, x: float, y: float, w: float, h: float, title: str, formula: str,
    curves: list, xr: tuple, yr: tuple, *, xlabel: str = "", note: str = "", marks: list | None = None,
) -> None:
    """Shared plot card: title, formula, axes and one or more (fn, color, width, dash) curves."""
    b.append(rect(x, y, w, h, rx=14, sw=1.3))
    b.append(text(x + 14, y + 22, title, size=11.5, weight=700, anchor="start"))
    b.append(text(x + 14, y + 40, formula, size=9.5, fill="var(--sf-ink)", weight=500, anchor="start"))
    px, py, pw, ph = x + 30, y + 54, w - 46, h - 96
    sx = lambda v: px + (v - xr[0]) / (xr[1] - xr[0]) * pw
    sy = lambda v: py + ph - (min(max(v, yr[0]), yr[1]) - yr[0]) / (yr[1] - yr[0]) * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    if yr[0] < 0 < yr[1]:
        b.append(path(f"M{px},{sy(0):.1f} L{px+pw},{sy(0):.1f}", sw=0.8, dash="2 3"))
    for fn, color, sw_, dash in curves:
        pts = []
        for i in range(161):
            v = xr[0] + (xr[1] - xr[0]) * i / 160
            e = fn(v)
            if e is None:
                continue
            pts.append((sx(v), sy(e)))
        b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in pts), stroke=color, sw=sw_, dash=dash, cls="sf-draw"))
    for (mv, lbl, col) in marks or []:
        b.append(path(f"M{sx(mv):.1f},{py} L{sx(mv):.1f},{py+ph}", stroke=col, sw=1, dash="3 3"))
        b.append(text(sx(mv), py - 3, lbl, size=8.5, fill=col))
    if xlabel:
        b.append(text(px + pw, py + ph + 13, xlabel, size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    if note:
        b.append(text(x + w / 2, y + h - 10, note, size=9, fill="var(--sf-muted)", weight=450))


def diagram_force_fields_overview() -> None:
    slug = "force-fields-overview"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Bounds-violation force field: E = E_dist + w_chiral·E_chiral + w_4D·E_4D", size=15, weight=700))
    l, u = 2.52, 3.81  # butane 1-4 bounds

    def e_dist(d):
        if d > u:
            return (d * d / (u * u) - 1) ** 2
        if d < l:
            return (2 * l * l / (l * l + d * d) - 1) ** 2
        return 0.0
    energy_panel(b, 24, 50, 380, 270, "E_dist for one pair (butane C1–C4)", "upper: (d²/u² − 1)²   lower: (2l²/(l² + d²) − 1)²",
                 [(e_dist, "var(--sf-accent)", 2.4, None)], (1.0, 5.4), (0, 1.1), xlabel="d (Å)",
                 marks=[(l, "l = 2.52", "var(--sf-rose)"), (u, "u = 3.81", "var(--sf-accent2)")],
                 note="flat bottom: zero force anywhere inside [l, u] · pairs with u − l > 5 Å skipped")
    # chiral volume
    x = 424
    b.append(rect(x, 50, 220, 270, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "E_chiral", size=11.5, weight=700, anchor="start"))
    b.append(text(x + 14, 90, "V = v₁ · (v₂ × v₃)", size=9.5, weight=500, anchor="start"))
    cx, cy = x + 110, 170
    nb = [(cx - 50, cy + 40, "var(--sf-accent)"), (cx + 52, cy + 34, "var(--sf-accent2)"), (cx + 6, cy - 58, "var(--sf-violet)"), (cx - 24, cy + 6, "var(--sf-warm)")]
    b.append(path(f"M{nb[0][0]},{nb[0][1]} L{nb[1][0]},{nb[1][1]} L{nb[2][0]},{nb[2][1]} Z", fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1, opacity=0.15))
    for i, (nx, ny, col) in enumerate(nb[:3]):
        b.append(path(f"M{nb[3][0]},{nb[3][1]} L{nx},{ny}", stroke=col, sw=2, marker=True))
        b.append(text(nx + (10 if i else -10), ny + 4, f"v{i+1}", size=9.5, fill=col))
    b.append(circle(nb[3][0], nb[3][1], 5, fill="var(--sf-ink)"))
    b.append(text(x + 110, 252, "penalize V outside [V_lower, V_upper]", size=9.5, weight=600))
    b.append(text(x + 110, 270, "sign of V encodes @ vs @@", size=9, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 110, 300, "(V − V_bound)²", size=10, fill="var(--sf-violet)", weight=700))
    # 4D + phases
    x = 664
    b.append(rect(x, 50, 212, 270, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "E_4D = Σ xᵢ₄²", size=11.5, weight=700, anchor="start"))
    b.append(text(x + 14, 90, "pull the 4th coordinate to 0", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    for i, (h4, col) in enumerate([(0.9, "var(--sf-line)"), (0.6, "var(--sf-line)"), (0.25, "var(--sf-warm)"), (0.03, "var(--sf-accent2)")]):
        bx = x + 30 + i * 42
        b.append(rect(bx, 200 - h4 * 90, 26, max(h4 * 90, 2), rx=3, fill=col, stroke="none", opacity=0.8))
    b.append(path(f"M{x+24},200 L{x+196},200", sw=1))
    b.append(text(x + 106, 216, "|x₄| during minimization →", size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(x + 106, 246, "phase 1: w_chiral 1.0 · w_4D 0.1", size=9.5, weight=600))
    b.append(text(x + 106, 264, "phase 2: w_chiral 0.2 · w_4D 1.0", size=9.5, weight=600))
    b.append(text(x + 106, 290, "then keep first 3 columns", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 336, "".join(b), label="Bounds-violation force field: flat-bottom distance penalty for the butane 1-4 pair (zero between 2.52 and 3.81 Å), signed chiral volume penalty, and the 4D penalty collapsed in two phases.")


TORSION_EXAMPLE = [(1, 1.0, 4.0), (2, 1.0, 3.1), (3, 1.0, 3.9), (4, -1.0, -0.8), (6, 1.0, 0.7)]


def diagram_force_fields_torsion() -> None:
    slug = "force-fields-torsion"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "M6 torsion: V(φ) = Σₖ Vₖ (1 + sₖ cos kφ) — one real CSD pattern", size=15, weight=700))
    b.append(text(450, 50, "[C:1][CX4H2:2]!@;-[NX3:3][C:4]   e.g. the C–C–N–C torsion of triethylamine", size=10, fill="var(--sf-muted)", weight=500, mono=True))
    cols = {1: "var(--sf-accent)", 2: "var(--sf-violet)", 3: "var(--sf-warm)", 4: "var(--sf-rose)", 6: "var(--sf-accent2)"}
    comps = [((lambda ph, k=k, sk=sk, V=V: V * (1 + sk * math.cos(k * math.radians(ph)))), cols[k], 1.3, "4 3") for k, sk, V in TORSION_EXAMPLE]
    total = lambda ph: sum(V * (1 + sk * math.cos(k * math.radians(ph))) for k, sk, V in TORSION_EXAMPLE)
    energy_panel(b, 24, 66, 560, 290, "components (dashed) and total (solid)", "V₁ = 4.0  V₂ = 3.1  V₃ = 3.9  V₄ = −0.8 (s₄ = −1)  V₆ = 0.7",
                 comps + [(total, "var(--sf-ink)", 2.6, None)], (-180, 180), (-2, 24), xlabel="φ (°)",
                 marks=[(-75, "−75°", "var(--sf-accent2)"), (75, "75°", "var(--sf-accent2)"), (-162, "−162°", "var(--sf-muted)"), (162, "162°", "var(--sf-muted)")],
                 note="minima at ±75° (6.9) and ±162° (7.3); syn 0° peaks at 23.4 — anti is split by a 0.3 bump")
    lx = 60
    for i, (k, sk, V) in enumerate(TORSION_EXAMPLE):
        b.append(path(f"M{lx + i*100},{378} l16,0", stroke=cols[k], sw=2, dash="4 3"))
        b.append(text(lx + i * 100 + 20, 381, f"k = {k}", size=9.5, fill=cols[k], anchor="start"))
    # Chebyshev card
    x = 604
    b.append(rect(x, 66, 272, 290, rx=14, sw=1.3))
    b.append(text(x + 14, 88, "one cos, five recurrences", size=11.5, weight=700, anchor="start"))
    b.append(text(x + 14, 106, "c = cos φ  (from the dot product)", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    rows = ["cos 2φ = 2c² − 1", "cos 3φ = 4c³ − 3c", "cos 4φ = 8c⁴ − 8c² + 1", "cos 5φ = 16c⁵ − 20c³ + 5c", "cos 6φ = 32c⁶ − 48c⁴ + 18c² − 1"]
    for i, r_ in enumerate(rows):
        b.append(rect(x + 14, 120 + i * 34, 244, 26, rx=7, fill="var(--sf-card2)", sw=0.8))
        b.append(text(x + 24, 137 + i * 34, r_, size=10, weight=500, anchor="start", mono=True))
    b.append(text(x + 136, 300, "no atan2, no extra trig calls", size=9.5, fill="var(--sf-accent2)", weight=600))
    b.append(text(x + 136, 318, "per torsion per energy evaluation", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 396, "".join(b), label="Real CSD torsion pattern with V1 to V6 Fourier components and their sum, plus the Chebyshev recurrences that compute cos kφ from a single cosine.")


def diagram_uff_energy_terms() -> None:
    slug = "uff-energy-terms"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "UFF energy terms with real parameters (sp³ carbon C_3 unless noted)", size=15, weight=700))
    r0, kb = 1.514, 699.7
    t0 = math.radians(109.47)
    C2 = 1 / (4 * math.sin(t0) ** 2); C1 = -4 * C2 * math.cos(t0); C0 = C2 * (2 * math.cos(t0) ** 2 + 1)
    W, H = 280, 200
    panels = [
        ("stretch", "½ k (r − r₀)²   k = 664.12 Z*ᵢZ*ⱼ / r₀³ = 700", [(lambda r: 0.5 * kb * (r - r0) ** 2, "var(--sf-accent)", 2.2, None)], (1.2, 1.85), (0, 40), "r (Å)", [(r0, "1.514", "var(--sf-accent)")], "C_3–C_3 harmonic well"),
        ("bend", "k (C₀ + C₁ cos θ + C₂ cos 2θ)", [(lambda a: 100 * (C0 + C1 * math.cos(math.radians(a)) + C2 * math.cos(2 * math.radians(a))), "var(--sf-accent2)", 2.2, None)], (60, 180), (0, 70), "θ (°)", [(109.47, "109.47", "var(--sf-accent2)")], "k = 100 for display · minimum at θ₀"),
        ("torsion", "½ V (1 − cos nφ₀ cos nφ)   n = 3, V = 2.119", [(lambda a: 0.5 * 2.119 * (1 + math.cos(3 * math.radians(a))), "var(--sf-violet)", 2.2, None)], (-180, 180), (0, 2.5), "φ (°)", [(60, "60", "var(--sf-violet)"), (180, "", "var(--sf-violet)")], "staggered minima at ±60°, 180°"),
        ("van der Waals", "ε[(x/r)¹² − 2(x/r)⁶]   x = 3.851, ε = 0.105", [(lambda r: 0.105 * ((3.851 / r) ** 12 - 2 * (3.851 / r) ** 6), "var(--sf-warm)", 2.2, None)], (3.3, 7.0), (-0.13, 0.15), "r (Å)", [(3.851, "x = 3.85", "var(--sf-warm)")], "well depth ε at r = x"),
        ("inversion (sp² C)", "k (1 − cos ω)   k = 6, planar ω₀ = 0", [(lambda w: 6 * (1 - math.cos(math.radians(w))), "var(--sf-rose)", 2.2, None)], (-40, 40), (0, 1.6), "ω (°)", [(0, "0", "var(--sf-rose)")], "planar sp² center has zero energy"),
        ("electrostatic", "332.0637 qᵢqⱼ / (ε r)   q = +0.3 / −0.3", [(lambda r: 332.0637 * -0.09 / r, "var(--sf-accent)", 2.2, None), (lambda r: 332.0637 * 0.09 / r, "var(--sf-line)", 1.4, "4 3")], (1.5, 8), (-22, 22), "r (Å)", [], "solid: opposite charges · dashed: like charges"),
    ]
    for i, (t_, f_, cv, xr, yr, xl, mk, nt) in enumerate(panels):
        x = 24 + (i % 3) * (W + 16); y = 50 + (i // 3) * (H + 14)
        energy_panel(b, x, y, W, H, t_, f_, cv, xr, yr, xlabel=xl, marks=mk, note=nt)
    b.append(text(450, 486, "E_total = E_stretch + E_bend + E_torsion + E_vdW + E_oop + E_elec  ·  minimized with L-BFGS using analytic gradients", size=10, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 500, "".join(b), label="The six UFF energy terms drawn from their formulas with sp3 carbon parameters: harmonic stretch at 1.514 Å, cosine-Fourier bend at 109.47°, threefold torsion, Lennard-Jones 12-6, inversion and Coulomb.")


def diagram_mmff94_energy_terms() -> None:
    slug = "mmff94-energy-terms"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "MMFF94 terms vs. their simple counterparts (C–C–C–C, MMFF type 1)", size=15, weight=700))
    kb, r0, cs = 4.258, 1.508, -2.0
    ka, th0, cb = 0.851, 109.608, -0.006981317
    V1, V2, V3 = 0.103, 0.681, 0.332
    Rs, eps = 3.9, 0.07
    W, H = 280, 200
    def es(r):
        d = r - r0
        return 143.9325 * 0.5 * kb * d * d * (1 + cs * d + 7 / 12 * cs * cs * d * d)
    def eb(t):
        d = t - th0
        return 0.043844 * 0.5 * ka * d * d * (1 + cb * d)
    def et(ph):
        p_ = math.radians(ph)
        return 0.5 * (V1 * (1 + math.cos(p_)) + V2 * (1 - math.cos(2 * p_)) + V3 * (1 + math.cos(3 * p_)))
    def e147(r):
        return eps * (1.07 * Rs / (r + 0.07 * Rs)) ** 7 * (1.12 * Rs ** 7 / (r ** 7 + 0.12 * Rs ** 7) - 2)
    def lj(r):
        return eps * ((Rs / r) ** 12 - 2 * (Rs / r) ** 6)
    panels = [
        ("bond stretch (quartic)", "143.93·½k Δr²[1 + c_s Δr + 7/12 c_s² Δr²]", [(lambda r: 143.9325 * 0.5 * kb * (r - r0) ** 2, "var(--sf-line)", 1.4, "4 3"), (es, "var(--sf-accent)", 2.4, None)], (1.28, 1.80), (0, 45), "r (Å)", [(r0, "r₀ 1.508", "var(--sf-accent)")], "dashed: harmonic · stretching is softer"),
        ("angle bend (cubic)", "0.043844·½k Δθ²[1 + c_b Δθ]", [(lambda t: 0.043844 * 0.5 * ka * (t - th0) ** 2, "var(--sf-line)", 1.4, "4 3"), (eb, "var(--sf-accent2)", 2.4, None)], (72, 160), (0, 40), "θ (°)", [(th0, "θ₀ 109.6", "var(--sf-accent2)")], "dashed: harmonic · softer toward linear"),
        ("torsion (3-term)", "½[V₁(1+cos φ) + V₂(1−cos 2φ) + V₃(1+cos 3φ)]", [(et, "var(--sf-violet)", 2.4, None)], (-180, 180), (0, 1.3), "φ (°)", [], "V₁ 0.103 · V₂ 0.681 · V₃ 0.332 kcal/mol"),
        ("vdW buffered 14-7", "ε(1.07R*/(R+0.07R*))⁷ (1.12R*⁷/(R⁷+0.12R*⁷) − 2)", [(lj, "var(--sf-line)", 1.4, "4 3"), (e147, "var(--sf-warm)", 2.4, None)], (3.25, 7.0), (-0.09, 0.12), "R (Å)", [(Rs, "R* 3.9", "var(--sf-warm)")], "dashed: LJ 12-6 · softer wall"),
        ("out-of-plane", "½·143.93·(π/180)²·k_oop·χ²", [(lambda c: 0.5 * 143.9325 * (math.pi / 180) ** 2 * 0.04 * c * c, "var(--sf-rose)", 2.4, None)], (-40, 40), (0, 4.5), "χ (°)", [(0, "planar", "var(--sf-rose)")], "sp² centers only · skipped if untyped"),
        ("buffered Coulomb", "332.0716 qᵢqⱼ / (R + 0.05)", [(lambda r: 332.0716 * -0.09 / r, "var(--sf-line)", 1.4, "4 3"), (lambda r: 332.0716 * -0.09 / (r + 0.05), "var(--sf-accent)", 2.4, None)], (0.3, 6), (-110, 0), "R (Å)", [], "dashed: bare 1/R · 1-4 pairs scaled 0.75"),
    ]
    for i, (t_, f_, cv, xr, yr, xl, mk, nt) in enumerate(panels):
        x = 24 + (i % 3) * (W + 16); y = 50 + (i // 3) * (H + 14)
        energy_panel(b, x, y, W, H, t_, f_, cv, xr, yr, xlabel=xl, marks=mk, note=nt)
    b.append(text(450, 486, "+ stretch–bend 2.5121(k_ijk Δr_ij + k_kji Δr_jk)Δθ couples bonds to the angle between them · examples use illustrative ε, R*, q, k_oop", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 500, "".join(b), label="MMFF94 terms drawn from their formulas next to simpler counterparts: quartic stretch versus harmonic, cubic bend, three-term torsion, buffered 14-7 versus Lennard-Jones, out-of-plane and buffered Coulomb.")


def hex_points(cx: float, cy: float, r: float) -> list[tuple[float, float]]:
    """Flat-sided hexagon vertices, index 0 at the right (where a substituent attaches)."""
    return [(cx + r * math.cos(math.radians(60 * i)), cy - r * math.sin(math.radians(60 * i))) for i in range(6)]


def phenol_graph(slug: str, cx: float, cy: float, r: float, *, hydrogens: bool, indices: bool, kekule: bool = True) -> str:
    """Phenol drawn as a static graph: SMILES order c1(0) c(1)… c1(5) O(6), O on atom 5."""
    pts = hex_points(cx, cy, r)
    # SMILES atom order walks the ring; atom 5 closes ring bond 1 back to atom 0 and carries O.
    order = [1, 2, 3, 4, 5, 0]  # hexagon slot for SMILES atom i
    P = [pts[order[i]] for i in range(6)]
    out = ""
    for i in range(6):
        a_, b_ = P[i], P[(i + 1) % 6]
        closing = i == 5
        out += path(f"M{a_[0]:.1f},{a_[1]:.1f} L{b_[0]:.1f},{b_[1]:.1f}", stroke="var(--sf-warm)" if closing else "var(--sf-ink)", sw=2.4, dash="5 4" if closing else None)
        if kekule and i % 2 == 0:
            mx, my = (a_[0] + b_[0]) / 2, (a_[1] + b_[1]) / 2
            ix, iy = cx + (mx - cx) * 0.78, cy + (my - cy) * 0.78
            dx, dy = (b_[0] - a_[0]) * 0.32, (b_[1] - a_[1]) * 0.32
            out += path(f"M{ix-dx:.1f},{iy-dy:.1f} L{ix+dx:.1f},{iy+dy:.1f}", stroke="var(--sf-ink)", sw=1.6)
    ox, oy = P[5][0] + r * 0.95, P[5][1]
    out += path(f"M{P[5][0]:.1f},{P[5][1]:.1f} L{ox:.1f},{oy:.1f}", stroke="var(--sf-ink)", sw=2.4)
    if hydrogens:
        for i in range(5):
            vx, vy = P[i][0] - cx, P[i][1] - cy
            hx, hy = P[i][0] + vx * 0.55, P[i][1] + vy * 0.55
            out += path(f"M{P[i][0]:.1f},{P[i][1]:.1f} L{hx:.1f},{hy:.1f}", sw=1.4)
            out += atom_ball(slug, hx, hy, "H", 6)
        hx, hy = ox + r * 0.38, oy - r * 0.42
        out += path(f"M{ox:.1f},{oy:.1f} L{hx:.1f},{hy:.1f}", sw=1.4)
        out += atom_ball(slug, hx, hy, "H", 6)
    for i in range(6):
        out += atom_ball(slug, P[i][0], P[i][1], "C", 11)
        if indices:
            vx, vy = P[i][0] - cx, P[i][1] - cy
            out += text(P[i][0] - vx * 0.36, P[i][1] - vy * 0.36 + 3, str(i), size=9, fill="var(--sf-accent)", weight=700)
    out += atom_ball(slug, ox, oy, "O", 12, cls="sf-pulse")
    if indices:
        out += text(ox, oy + 26, "6", size=9, fill="var(--sf-accent)", weight=700)
    return out


def diagram_smiles_pipeline() -> None:
    slug = "smiles-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    W, H, Y = 254, 290, 60
    cards = [("1", "Tokenize", "var(--sf-accent)"), ("2", "Build graph", "var(--sf-accent2)"), ("3", "Perceive", "var(--sf-violet)"), ("4", "Molecule", "var(--sf-warm)")]
    b.append(text(550, 30, "From the string  c1ccccc1O  to a molecular graph", size=15, weight=700))
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 · tokens with their roles
    x = xs[0]
    toks = [("c", "aromatic C", "var(--sf-ink)"), ("1", "ring open", "var(--sf-warm)"), ("c", "", "var(--sf-ink)"), ("c", "", "var(--sf-ink)"), ("c", "", "var(--sf-ink)"),
            ("c", "", "var(--sf-ink)"), ("c", "", "var(--sf-ink)"), ("1", "ring close", "var(--sf-warm)"), ("O", "aliphatic O", "var(--sf-rose)")]
    for i, (t_, role, col) in enumerate(toks):
        tx = x + 16 + i * 25
        b.append(rect(tx, Y + 52, 21, 28, rx=5, fill="var(--sf-card2)", stroke=col if col != "var(--sf-ink)" else "var(--sf-line)", sw=1))
        b.append(text(tx + 10.5, Y + 71, t_, size=13, fill=col, weight=700, mono=True))
    atom_idx = [0, None, 1, 2, 3, 4, 5, None, 6]
    for i, ai in enumerate(atom_idx):
        if ai is not None:
            b.append(text(x + 26.5 + i * 25, Y + 94, str(ai), size=8.5, fill="var(--sf-accent)", weight=700))
    rows = [("c", "aromatic atom → aromatic flag", "var(--sf-ink)"), ("1", "ring-bond digit: opens, then closes", "var(--sf-warm)"), ("O", "organic subset, implicit H", "var(--sf-rose)")]
    for k, (t_, d_, col) in enumerate(rows):
        b.append(text(x + 22, Y + 130 + k * 22, t_, size=11, fill=col, weight=700, mono=True))
        b.append(text(x + 38, Y + 130 + k * 22, d_, size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(x + W / 2, Y + 230, "blue numbers = atom index", size=9.5, fill="var(--sf-accent)", weight=600))
    b.append(text(x + W / 2, Y + 250, "9 tokens → 7 heavy atoms", **CAPTION))
    # 2 · graph with ring closure
    x = xs[1]
    b.append(phenol_graph(slug, x + 110, Y + 140, 46, hydrogens=False, indices=True, kekule=False))
    b.append(text(x + W / 2, Y + 230, "ring digit 1 pairs atom 5 → atom 0", size=9.5, fill="var(--sf-warm)", weight=600))
    b.append(text(x + W / 2, Y + 250, "7 atoms · 7 bonds (6 aromatic, 1 single)", **CAPTION))
    # 3 · perception: hybridization, implicit H, SSSR
    x = xs[2]
    items = [
        ("hybridization", "aromatic C → SP2 · O (2 single) → SP3"),
        ("implicit H", "nH = v − Σ order − |q|"),
        ("", "C0–C4: 1 H each · C5: 0 · O: 1"),
        ("SSSR", "μ = |E| − |V| + 1 = 7 − 7 + 1 = 1"),
        ("", "one 6-membered ring {0…5}"),
    ]
    yy = Y + 60
    for head, body in items:
        if head:
            b.append(text(x + 16, yy, head, size=10.5, fill="var(--sf-violet)", weight=700, anchor="start"))
            yy += 18
        b.append(text(x + 16, yy, body, size=9.5, weight=500, anchor="start", mono=bool(head == "" or "=" in body)))
        yy += 26
    b.append(text(x + W / 2, Y + 250, "all needed by bounds + ETKDG", **CAPTION))
    # 4 · final molecule
    x = xs[3]
    b.append(phenol_graph(slug, x + 112, Y + 140, 40, hydrogens=True, indices=False))
    b.append(text(x + W / 2, Y + 230, "13 atoms · 13 bonds · 1 ring", size=10, weight=600))
    b.append(text(x + W / 2, Y + 250, "Molecule { graph, ring_info }", size=9.5, fill="var(--sf-muted)", weight=500, mono=True))
    write(slug, 1100, 370, "".join(b), label="Parsing c1ccccc1O: nine tokens with atom indices and ring-bond digits, ring closure from atom 5 to atom 0, hybridization, implicit hydrogens and SSSR perception, and the final 13-atom phenol molecule.")


def diagram_smiles_phenol() -> None:
    slug = "smiles-phenol"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Phenol, c1ccccc1O, as the parser returns it", size=15, weight=700))
    b.append(rect(24, 50, 440, 320, rx=14, sw=1.3))
    b.append(phenol_graph(slug, 200, 210, 72, hydrogens=True, indices=True))
    b.append(text(244, 360, "orange dashed = ring-closure bond (digit 1)", size=9.5, fill="var(--sf-warm)", weight=600))
    x = 484
    b.append(rect(x, 50, 392, 320, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "per-atom properties", size=11.5, weight=700, anchor="start"))
    hdr = ["idx", "elem", "arom", "hyb", "nH", "ring"]
    cx_ = [x + 20, x + 64, x + 120, x + 176, x + 232, x + 276]
    for c_, h_ in zip(cx_, hdr):
        b.append(text(c_, 98, h_, size=9.5, fill="var(--sf-muted)", weight=700, anchor="start"))
    rows = [(str(i), "C", "yes", "SP2", "1" if i < 5 else "0", "6") for i in range(6)] + [("6", "O", "no", "SP3", "1", "—")]
    for r_, row in enumerate(rows):
        y = 118 + r_ * 22
        if r_ == 5 or r_ == 6:
            b.append(rect(x + 12, y - 15, 368, 21, rx=5, fill="var(--sf-soft)", stroke="none"))
        for c_, v in zip(cx_, row):
            b.append(text(c_, y, v, size=10, weight=500, anchor="start", mono=True, fill="var(--sf-rose)" if v == "O" else "var(--sf-ink)"))
    b.append(text(x + 16, 290, "totals: 6 C + 1 O + 6 H = 13 atoms", size=10.5, weight=700, anchor="start"))
    b.append(text(x + 16, 310, "bonds: 6 ring + C–O + 6 X–H = 13", size=10, weight=500, anchor="start"))
    b.append(text(x + 16, 330, "SSSR: 1 ring of size 6", size=10, weight=500, anchor="start"))
    b.append(text(x + 16, 350, "C5 has no H: its 3rd neighbor is O", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 386, "".join(b), label="Phenol from c1ccccc1O with atom indices, the ring-closure bond, explicit hydrogens and a per-atom table: six aromatic SP2 carbons, one SP3 oxygen, 13 atoms and one six-membered ring.")


def diagram_smarts_pattern() -> None:
    slug = "smarts-pattern"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Torsion matching on butane's central bond — first match wins", size=15, weight=700))
    # target molecule with mapped atoms
    b.append(rect(24, 50, 300, 320, rx=14, sw=1.3))
    b.append(text(174, 74, "target: butane, bond C2–C3", size=11.5, weight=700))
    pts = [(70, 230), (130, 160), (210, 220), (270, 150)]
    b.append(path("M" + " L".join(f"{x},{y}" for x, y in pts), stroke="var(--sf-ink)", sw=3))
    b.append(path(f"M{pts[1][0]},{pts[1][1]} L{pts[2][0]},{pts[2][1]}", stroke="var(--sf-warm)", sw=9, opacity=0.35, cls="sf-pulse"))
    for i, (x, y) in enumerate(pts):
        b.append(atom_ball(slug, x, y, "C", 13))
        b.append(rect(x - 13, y + (18 if i % 2 == 0 else -40), 26, 20, rx=6, fill="var(--sf-violet)", stroke="none"))
        b.append(text(x, y + (32 if i % 2 == 0 else -26), f":{i+1}", size=10, fill="#fff", weight=700))
    b.append(text(174, 300, "rotatable: single, not in ring (!@)", size=9.5, fill="var(--sf-muted)", weight=500))
    b.append(text(174, 318, "both centers CX4 (4 connections)", size=9.5, fill="var(--sf-muted)", weight=500))
    b.append(text(174, 346, "→ coefficients for this torsion", size=10, fill="var(--sf-warm)", weight=650))
    # priority list of real patterns
    x = 344
    b.append(rect(x, 50, 532, 320, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "CSD table, scanned in priority order", size=11.5, weight=700, anchor="start"))
    tries = [
        ("[O:1]=[C:2]!@;-[O:3]~[CH0:4]", "no C=O in butane", False),
        ("[CX3:1]=[CX3:2]!@;-[CH2:3][C:4]", ":2 must be CX3 (sp²); C2 is CX4", False),
        ("[cH0:1][c:2]!@;-[CX4H2:3][CX4:4]", ":1, :2 must be aromatic", False),
        ("[!#1:1][CX4:2]!@;-[CX4:3][!#1:4]", "all four atoms map ✓", True),
    ]
    for i, (pat, why, ok) in enumerate(tries):
        y = 92 + i * 56
        col = "var(--sf-accent2)" if ok else "var(--sf-rose)"
        b.append(rect(x + 14, y, 504, 46, rx=9, fill="var(--sf-soft)" if ok else "var(--sf-card2)", stroke=col, sw=1.4 if ok else 0.8, cls="sf-pulse" if ok else ""))
        b.append(circle(x + 34, y + 23, 9, fill=col))
        b.append(path(f"M{x+30},{y+23} l3,3 l6,-6" if ok else f"M{x+30},{y+19} l8,8 M{x+38},{y+19} l-8,8", stroke="#fff", sw=1.8))
        b.append(text(x + 52, y + 20, pat, size=10.5, weight=600, anchor="start", mono=True))
        b.append(text(x + 52, y + 37, why, size=9.5, fill=col, weight=500, anchor="start"))
    b.append(text(x + 266, 330, "V₃ = 7.0, s₃ = +1  →  V(φ) = 7 (1 + cos 3φ)", size=11, fill="var(--sf-accent2)", weight=700))
    b.append(text(x + 266, 352, "no match in all 837 → basic-knowledge rule (SP3–SP3: V₃ = 7.0)", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 386, "".join(b), label="SMARTS torsion matching on butane: three real CSD patterns fail for stated reasons and [!#1][CX4]!@;-[CX4][!#1] matches first, giving V3 = 7.0.")


def diagram_esp_grid() -> None:
    slug = "esp-grid"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "ESP of water on the z = 0 slice of its grid (Mulliken charges, Φ = Σ q_A / |r − R_A|)", size=14, weight=700))
    pad, h = 3.0, 0.5
    xs_ = [p_[0] for p_ in WATER_XYZ]; ys_ = [p_[1] for p_ in WATER_XYZ]
    x0, y0 = min(xs_) - pad, min(ys_) - pad
    nx = int((max(xs_) - min(xs_) + 2 * pad) / h) + 1
    ny = int((max(ys_) - min(ys_) + 2 * pad) / h) + 1
    nz = int((0 + 2 * pad) / h) + 1
    cell = 21
    gx, gy = 40, 60
    R = 0.25
    for j in range(ny):
        for i in range(nx):
            X, Y = x0 + i * h, y0 + (ny - 1 - j) * h
            phi = 0.0
            for (ax, ay, _), q in zip(WATER_XYZ, WATER_MULLIKEN):
                d = max(math.hypot(X - ax, Y - ay), 0.35)
                phi += q / d
            phi *= 0.529177  # e/Å → a.u. (e/bohr)
            t = max(-1.0, min(1.0, phi / R))
            if t < 0:
                col = f"rgb(255,{int((1 + t) * 255)},{int((1 + t) * 255)})"
            else:
                col = f"rgb({int((1 - t) * 255)},{int((1 - t) * 255)},255)"
            b.append(f'<rect x="{gx + i * cell}" y="{gy + j * cell}" width="{cell - 1.2}" height="{cell - 1.2}" rx="2" fill="{col}"/>')
    k = cell / h
    ox, oy = gx + (0 - x0) * k + cell / 2 - 0.6, gy + (ny - 1 - (0 - y0) / h) * cell + cell / 2 - 0.6
    b.append(water_2d(slug, ox, oy, k, r=10))
    b.append(rect(gx - 2, gy - 2, nx * cell + 2, ny * cell + 2, rx=4, fill="none", sw=1.2, cls="sf-pulse"))
    b.append(text(gx + nx * cell / 2, gy + ny * cell + 18, f"{nx} × {ny} points shown (spacing 0.5 Å, padding 3.0 Å)", size=9.5, fill="var(--sf-muted)", weight=500))
    # right column
    x = 40 + nx * cell + 30
    b.append(rect(x, 56, 900 - x - 24, 300, rx=14, sw=1.3))
    b.append(text(x + 16, 80, "grid size", size=11.5, weight=700, anchor="start"))
    b.append(text(x + 16, 100, "Nᵢ = ⌊(maxᵢ − minᵢ + 2·pad) / h⌋ + 1", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 16, 120, f"water: {nx} × {ny} × {nz} = {nx*ny*nz} points", size=10, fill="var(--sf-accent)", weight=700, anchor="start"))
    b.append(text(x + 16, 152, "color map, R = 0.25 a.u.", size=11.5, weight=700, anchor="start"))
    for i in range(41):
        t = -1 + i / 20
        col = f"rgb(255,{int((1 + t) * 255)},{int((1 + t) * 255)})" if t < 0 else f"rgb({int((1 - t) * 255)},{int((1 - t) * 255)},255)"
        b.append(f'<rect x="{x + 16 + i * 5}" y="164" width="5.2" height="16" fill="{col}"/>')
    b.append(text(x + 16, 194, "−R", size=9, fill="var(--sf-rose)", anchor="start"))
    b.append(text(x + 118, 194, "0", size=9, fill="var(--sf-muted)"))
    b.append(text(x + 221, 194, "+R", size=9, fill="var(--sf-accent)", anchor="end"))
    b.append(text(x + 16, 222, "red: electron-rich (O lone pairs)", size=10, fill="var(--sf-rose)", weight=600, anchor="start"))
    b.append(text(x + 16, 240, "blue: electron-poor (O–H hydrogens)", size=10, fill="var(--sf-accent)", weight=600, anchor="start"))
    b.append(text(x + 16, 272, "charges: O −1.168, H +0.584", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(x + 16, 290, "points within 0.35 Å clamped for display", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(x + 16, 320, "parallel: Rayon chunks over points", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    b.append(text(x + 16, 338, "export: Gaussian .cube (bohr)", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    write(slug, 900, 400, "".join(b), label=f"Electrostatic potential of water computed from sci-form Mulliken charges on the z = 0 slice of a {nx} by {ny} by {nz} grid, colored red for negative and blue for positive potential.")


def population_bars(b: list, x: float, y: float, w: float, pops: list, Z: list, color: str) -> None:
    labels = ["O", "H", "H"]
    for i, (pop, z, lb) in enumerate(zip(pops, Z, labels)):
        yy = y + i * 40
        b.append(text(x, yy + 14, lb, size=11, weight=700, fill="var(--sf-rose)" if lb == "O" else "var(--sf-ink)"))
        sc = w / 8
        b.append(rect(x + 16, yy, z * sc, 20, rx=4, fill="none", stroke="var(--sf-line)", sw=1, opacity=None))
        b.append(rect(x + 16, yy, pop * sc, 20, rx=4, fill=color, stroke="none", opacity=0.6, cls="sf-rise"))
        b.append(text(x + 22 + max(pop, z) * sc, yy + 14, f"{pop:.3f} e⁻  (Z = {z})", size=9.5, weight=600, anchor="start"))


def diagram_population_mulliken() -> None:
    slug = "population-mulliken"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Mulliken populations of water (sci-form EHT)", size=15, weight=700))
    # left: matrices
    b.append(rect(24, 50, 300, 300, rx=14, sw=1.3))
    b.append(text(174, 74, "q_A = Z_A − Σ_{μ∈A} (PS)_μμ", size=11.5, weight=700))
    labs = ["s", "x", "y", "z", "H", "H"]
    owner = [0, 0, 0, 0, 1, 2]
    cs = 30; mx, my = 66, 100
    for i in range(6):
        b.append(text(mx - 10, my + i * cs + 19, labs[i], size=9, fill="var(--sf-muted)"))
        b.append(text(mx + i * cs + 14, my - 6, labs[i], size=9, fill="var(--sf-muted)"))
        for j in range(6):
            same = owner[i] == owner[j]
            diag = i == j
            fill = "var(--sf-violet)" if diag else ("var(--sf-accent)" if same else "var(--sf-warm)")
            op = 0.8 if diag else (0.25 if same else 0.45)
            b.append(rect(mx + j * cs, my + i * cs, cs - 2, cs - 2, rx=3, fill=fill, stroke="none", opacity=op))
    for k, (lbl, col) in enumerate([("(PS)μμ — summed per atom", "var(--sf-violet)"), ("O–H overlap pop. split 50 / 50", "var(--sf-warm)")]):
        b.append(rect(40, 290 + k * 22, 10, 10, rx=2, fill=col, stroke="none", opacity=0.8))
        b.append(text(56, 299 + k * 22, lbl, size=9.5, weight=500, anchor="start"))
    # middle: populations
    x = 344
    b.append(rect(x, 50, 300, 300, rx=14, sw=1.3))
    b.append(text(x + 150, 74, "gross populations (valence e⁻)", size=11.5, weight=700))
    pops = [6 - WATER_MULLIKEN[0], 1 - WATER_MULLIKEN[1], 1 - WATER_MULLIKEN[2]]
    population_bars(b, x + 22, 100, 160, pops, [6, 1, 1], "var(--sf-violet)")
    b.append(text(x + 150, 250, "Tr(PS) = 7.168 + 0.416 + 0.416 = 8", size=10, weight=600))
    b.append(text(x + 150, 270, "= N_electrons → Σ q_A = 0", size=10, fill="var(--sf-accent2)", weight=600))
    b.append(text(x + 150, 300, "outline = Z (valence electrons)", size=9.5, fill="var(--sf-muted)", weight=450))
    # right: charges on molecule
    x = 664
    b.append(rect(x, 50, 212, 300, rx=14, sw=1.3))
    b.append(text(x + 106, 74, "Mulliken charges", size=11.5, weight=700))
    b.append(water_2d(slug, x + 106, 190, 70, charges=WATER_MULLIKEN, r=16))
    b.append(text(x + 106, 300, "EHT exaggerates polarity", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(x + 106, 316, "(qualitative trend: O ≪ H)", size=9.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 366, "".join(b), label="Mulliken analysis of water from sci-form EHT: the PS diagonal summed per atom gives 7.168 electrons on O and 0.416 on each H, so charges are −1.168 and +0.584 and sum to zero.")


def diagram_population_lowdin() -> None:
    slug = "population-lowdin"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Löwdin populations: orthogonalize first, then read the diagonal", size=15, weight=700))
    W, H, Y = 196, 236, 56
    cards = [("1", "Diagonalize S", "var(--sf-accent)"), ("2", "S¹ᐟ² = UΛ¹ᐟ²Uᵀ", "var(--sf-accent2)"), ("3", "P̃ = S¹ᐟ² P S¹ᐟ²", "var(--sf-violet)"), ("4", "Charges", "var(--sf-warm)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 eigenvalues of S
    x = xs[0]
    ev = [1.86, 1.27, 1.0, 1.0, 0.61, 0.26]
    for k, v in enumerate(ev):
        bx = x + 24 + k * 26
        b.append(rect(bx, Y + 160 - v * 50, 18, v * 50, rx=3, fill="var(--sf-accent)", stroke="none", opacity=0.75))
    b.append(path(f"M{x+20},{Y+160} L{x+180},{Y+160}", sw=1))
    b.append(text(x + W / 2, Y + 184, "S = UΛUᵀ, all λ > 0", size=10, weight=600))
    b.append(text(x + W / 2, Y + 202, "illustrative overlap spectrum", **CAPTION))
    # 2 orthogonal basis
    x = xs[1]; cx, cy = x + 98, Y + 120
    b.append(path(f"M{cx},{cy} L{cx+60},{cy-12}", stroke="var(--sf-line)", sw=1.6, dash="4 3"))
    b.append(path(f"M{cx},{cy} L{cx+30},{cy-54}", stroke="var(--sf-line)", sw=1.6, dash="4 3"))
    b.append(path(f"M{cx},{cy} L{cx+62},{cy}", stroke="var(--sf-accent2)", sw=2.2, marker=True))
    b.append(path(f"M{cx},{cy} L{cx},{cy-62}", stroke="var(--sf-accent2)", sw=2.2, marker=True))
    b.append(text(cx + 66, cy - 16, "χ₁", size=9.5, fill="var(--sf-muted)", anchor="start"))
    b.append(text(cx + 34, cy - 58, "χ₂", size=9.5, fill="var(--sf-muted)", anchor="start"))
    b.append(text(cx + 70, cy + 4, "λ₁", size=9.5, fill="var(--sf-accent2)", anchor="start"))
    b.append(text(cx - 6, cy - 66, "λ₂", size=9.5, fill="var(--sf-accent2)", anchor="end"))
    b.append(text(x + W / 2, Y + 184, "symmetric: λ → √λ per eigenvector", size=10, weight=600))
    b.append(text(x + W / 2, Y + 202, "orthonormal AOs closest to originals", **CAPTION))
    # 3 diagonal only
    x = xs[2]; mx, my, cs = x + 40, Y + 52, 20
    for i in range(6):
        for j in range(6):
            b.append(rect(mx + j * cs, my + i * cs, cs - 2, cs - 2, rx=2.5, fill="var(--sf-violet)" if i == j else "var(--sf-card2)", stroke="none" if i == j else "var(--sf-line)", sw=0.5, opacity=0.85 if i == j else None))
    b.append(text(x + W / 2, Y + 184, "0 ≤ P̃_μμ ≤ 2 always", size=10, weight=600))
    b.append(text(x + W / 2, Y + 202, "no 50/50 overlap split", **CAPTION))
    # 4 charges compared
    x = xs[3]
    b.append(water_2d(slug, x + 98, Y + 112, 58, charges=WATER_LOWDIN, r=14))
    b.append(text(x + W / 2, Y + 184, "Löwdin O −1.096 · H +0.548", size=10, fill="var(--sf-violet)", weight=700))
    b.append(text(x + W / 2, Y + 202, "Mulliken O −1.168 · H +0.584", **CAPTION))
    b.append(text(x + W / 2, Y + 220, "Σ q = 0 in both schemes", **CAPTION))
    write(slug, 900, 308, "".join(b), label="Löwdin analysis: diagonalize the overlap, build the symmetric orthogonalizer, read the diagonal of the transformed density; water gives O −1.096 and H +0.548 versus Mulliken −1.168 and +0.584.")


def diagram_dipole_vector() -> None:
    slug = "dipole-vector"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "μ = Σ q_A R_A for water (sci-form Mulliken charges, Å)", size=15, weight=700))
    b.append(rect(24, 50, 420, 320, rx=14, sw=1.3))
    cx, cy, k = 234, 170, 120
    b.append(path(f"M{cx-180},{cy} L{cx+180},{cy}", sw=0.8, dash="2 4"))
    b.append(path(f"M{cx},{cy-110} L{cx},{cy+170}", sw=0.8, dash="2 4"))
    b.append(text(cx + 176, cy - 6, "x", size=9, fill="var(--sf-muted)"))
    b.append(text(cx + 8, cy - 100, "y", size=9, fill="var(--sf-muted)", anchor="start"))
    b.append(water_2d(slug, cx, cy, k, charges=WATER_MULLIKEN, r=16))
    # dipole arrow from O toward H side (negative → positive, page convention)
    mx, my = WATER_DIPOLE_D[0], WATER_DIPOLE_D[1]
    ox, oy = cx + WATER_XYZ[0][0] * k, cy - WATER_XYZ[0][1] * k
    L = 38
    x1, y1, x2, y2 = ox + 34, oy + 8, ox + 34 + mx * L, oy + 8 - my * L
    b.append(path(f"M{x1:.1f},{y1:.1f} L{x2:.1f},{y2 - 12:.1f}", stroke="var(--sf-warm)", sw=3.5, cls="sf-pulse"))
    b.append(path(f"M{x2-8:.1f},{y2-14:.1f} L{x2:.1f},{y2:.1f} L{x2+8:.1f},{y2-14:.1f} Z", fill="var(--sf-warm)", stroke="var(--sf-warm)", sw=1))
    b.append(text(ox + 44, oy + 60, "μ", size=13, fill="var(--sf-warm)", weight=800, anchor="start"))
    b.append(text(234, 352, "arrow: from negative (O) toward positive (H) charge", size=9.5, fill="var(--sf-muted)", weight=500))
    # table
    x = 464
    b.append(rect(x, 50, 412, 320, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "component sums (e·Å)", size=11.5, weight=700, anchor="start"))
    hdr = ["atom", "q", "x", "y", "q·y"]
    cols = [x + 20, x + 84, x + 156, x + 228, x + 310]
    for c_, h_ in zip(cols, hdr):
        b.append(text(c_, 100, h_, size=9.5, fill="var(--sf-muted)", weight=700, anchor="start"))
    tot = 0.0
    for i, (lb, q, (X, Y, _)) in enumerate(zip(["O", "H", "H"], WATER_MULLIKEN, WATER_XYZ)):
        y = 124 + i * 24
        tot += q * Y
        for c_, v in zip(cols, [lb, f"{q:+.3f}", f"{X:+.3f}", f"{Y:+.3f}", f"{q*Y:+.3f}"]):
            b.append(text(c_, y, v, size=10, weight=500, anchor="start", mono=True))
    b.append(path(f"M{x+16},{192} L{x+396},{192}", sw=1))
    b.append(text(x + 20, 212, "μ_y = Σ q·y", size=10, weight=700, anchor="start"))
    b.append(text(cols[4], 212, f"{tot:+.3f}", size=10, weight=700, anchor="start", mono=True))
    b.append(text(x + 20, 244, f"× 4.80321 → μ_y = {tot*4.80321:+.2f} D", size=11, fill="var(--sf-warm)", weight=700, anchor="start"))
    b.append(text(x + 20, 266, "μ_x = +0.04 D, μ_z = 0", size=10, weight=500, anchor="start"))
    b.append(text(x + 20, 292, "|μ| = 3.24 D", size=13, fill="var(--sf-warm)", weight=800, anchor="start"))
    b.append(text(x + 20, 316, "experiment 1.85 D: EHT charges are too large,", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 20, 332, "so use EHT dipoles for trends, not values", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 386, "".join(b), label="Water dipole from sci-form Mulliken charges and coordinates: the q·y sum gives −0.675 e·Å, which is 3.24 D, larger than the experimental 1.85 D because EHT charges are exaggerated.")


def diagram_dos_smearing() -> None:
    slug = "dos-smearing"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Water DOS: six EHT eigenvalues broadened by Gaussians", size=15, weight=700))
    px, py, pw, ph = 70, 70, 560, 250
    lo, hi = -36.0, -4.0
    sx = lambda e: px + (e - lo) / (hi - lo) * pw
    def dos(E, sg):
        return sum(math.exp(-((E - e) ** 2) / (2 * sg * sg)) / (sg * math.sqrt(2 * math.pi)) for e in WATER_EHT)
    ymax = dos(-15.3, 0.3) * 1.1
    sy = lambda v: py + ph - v / ymax * ph
    b.append(rect(px - 46, py - 20, pw + 70, ph + 66, rx=14, sw=1.3))
    b.append(rect(sx(lo), py, sx(-14.8 + 0.01) - sx(lo), ph, rx=0, fill="var(--sf-violet)", stroke="none", opacity=0.07))
    b.append(axes(px, py, pw, ph, xlabel="E (eV)", ylabel="D(E) (states/eV)"))
    for e in WATER_EHT:
        occ = e <= -14.8
        b.append(path(f"M{sx(e):.1f},{py+ph} L{sx(e):.1f},{py+ph-26}", stroke="var(--sf-violet)" if occ else "var(--sf-rose)", sw=2.2))
    for sg, col, sw_, dash in [(0.1, "var(--sf-line)", 1.2, "3 3"), (0.3, "var(--sf-accent)", 2.4, None), (0.5, "var(--sf-accent2)", 1.6, "6 4")]:
        pts = [(sx(lo + (hi - lo) * i / 400), sy(min(dos(lo + (hi - lo) * i / 400, sg), ymax))) for i in range(401)]
        b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in pts), stroke=col, sw=sw_, dash=dash, cls="sf-draw"))
    for e in (-35, -30, -25, -20, -15, -10, -5):
        b.append(text(sx(e), py + ph + 13, str(e), size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(sx(-14.8), py - 6, "HOMO −14.80", size=9, fill="var(--sf-violet)", anchor="end"))
    b.append(text(sx(-8.96), py - 6, "LUMO −8.96", size=9, fill="var(--sf-rose)", anchor="start"))
    # legend / formula
    x = 676
    b.append(rect(x, 50, 200, 316, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "D(E) = Σᵢ G(E − εᵢ, σ)", size=10.5, weight=700, anchor="start"))
    for k, (lbl, col, dash) in enumerate([("σ = 0.1 eV (clipped)", "var(--sf-line)", "3 3"), ("σ = 0.3 eV (default)", "var(--sf-accent)", None), ("σ = 0.5 eV", "var(--sf-accent2)", "6 4")]):
        y = 102 + k * 22
        b.append(path(f"M{x+16},{y} l22,0", stroke=col, sw=2, dash=dash))
        b.append(text(x + 44, y + 3, lbl, size=9.5, weight=500, anchor="start"))
    b.append(path(f"M{x+27},{178} l0,-12", stroke="var(--sf-violet)", sw=2.2))
    b.append(text(x + 44, 176, "occupied εᵢ", size=9.5, weight=500, anchor="start"))
    b.append(path(f"M{x+27},{200} l0,-12", stroke="var(--sf-rose)", sw=2.2))
    b.append(text(x + 44, 198, "virtual εᵢ", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 16, 234, "∫ D(E) dE = 6 orbitals", size=10, weight=600, anchor="start"))
    b.append(text(x + 16, 254, "grid: min εᵢ − 5σ … max + 5σ", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 16, 272, "PDOS: weights w_Ai = Σ|C_μi|²", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 16, 290, "Σ_A D_A(E) = D(E)", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 16, 318, "−15.90 / −15.26 pair: shoulder", size=9.5, fill="var(--sf-accent)", weight=600, anchor="start"))
    b.append(text(x + 16, 334, "at 0.3 eV, single peak at 0.5 eV", size=9.5, fill="var(--sf-accent)", weight=600, anchor="start"))
    write(slug, 900, 380, "".join(b), label="Water density of states from its six sci-form EHT eigenvalues, broadened with Gaussian widths 0.1, 0.3 and 0.5 eV; nearby levels merge as sigma grows.")


def kabsch_example():
    """Real Kabsch run: butane target Q, reference P = rotated, shifted, perturbed copy."""
    import numpy as np
    Q = np.array(BUTANE_X)
    th = math.radians(40)
    Rz = np.array([[math.cos(th), -math.sin(th), 0], [math.sin(th), math.cos(th), 0], [0, 0, 1]])
    noise = np.array([[0.05, -0.04, 0.02], [-0.03, 0.06, -0.01], [0.02, -0.05, 0.04], [-0.04, 0.03, -0.05]])
    P = Q @ Rz.T + np.array([2.2, 0.9, 0.0]) + noise
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    H = Pc.T @ Qc
    U, S_, Vt = np.linalg.svd(H)
    d = float(np.sign(np.linalg.det(Vt.T @ U.T)))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    aligned = Pc @ R.T
    rmsd = float(np.sqrt(((aligned - Qc) ** 2).sum(1).mean()))
    return P, Q, Pc, Qc, H, R, d, aligned, rmsd


def chain_2d(pts: list, color: str, *, label: str = "", dash: str | None = None, r: float = 7) -> str:
    out = path("M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts), stroke=color, sw=2.4, dash=dash)
    for i, (x, y) in enumerate(pts):
        out += circle(x, y, r, fill=color, opacity=0.9)
        if label:
            out += text(x, y + 3, str(i + 1), size=8, fill="#fff", weight=700)
    return out


def diagram_kabsch_alignment() -> None:
    slug = "kabsch-alignment"
    begin(slug)
    b = [element_defs(slug)]
    P, Q, Pc, Qc, H, R, d, aligned, rmsd = kabsch_example()
    b.append(text(550, 30, "Kabsch on butane: reference P (rotated 40°, shifted, ±0.06 Å noise) onto target Q", size=15, weight=700))
    W, Hc, Y = 196, 250, 54
    cards = [("1", "Inputs", "var(--sf-accent)"), ("2", "Center", "var(--sf-accent2)"), ("3", "H = P̄ᵀQ̄", "var(--sf-violet)"), ("4", "R from SVD", "var(--sf-warm)"), ("5", "Align + RMSD", "var(--sf-rose)")]
    xs = stage_row(b, 24, Y, W, Hc, 18, cards)
    proj = lambda M, cx, cy, k: [(cx + v[0] * k, cy - v[1] * k) for v in M]
    # 1 inputs
    x = xs[0]
    b.append(chain_2d(proj(Q, x + 70, Y + 130, 26), "var(--sf-accent2)", label="1"))
    b.append(chain_2d(proj(P, x + 70, Y + 130, 26), "var(--sf-rose)", label="1", dash="4 3"))
    b.append(text(x + W / 2, Y + 214, "Q target · P reference", size=10, weight=600))
    b.append(text(x + W / 2, Y + 232, "same atom order, N = 4", **CAPTION))
    # 2 centered
    x = xs[1]
    b.append(chain_2d(proj(Qc, x + 98, Y + 130, 30), "var(--sf-accent2)", label="1"))
    b.append(chain_2d(proj(Pc, x + 98, Y + 130, 30), "var(--sf-rose)", label="1", dash="4 3"))
    b.append(circle(x + 98, Y + 130, 3, fill="var(--sf-ink)"))
    b.append(text(x + W / 2, Y + 214, "subtract centroids", size=10, weight=600))
    b.append(text(x + W / 2, Y + 232, "rotation left: 40° about z", **CAPTION))
    # 3 H matrix numbers
    x = xs[2]
    for i in range(3):
        for j in range(3):
            v = H[i][j]
            b.append(rect(x + 26 + j * 50, Y + 70 + i * 34, 46, 30, rx=4, fill="var(--sf-violet)", stroke="none", opacity=round(0.1 + min(abs(v), 6) / 8, 2)))
            b.append(text(x + 49 + j * 50, Y + 89 + i * 34, f"{v:+.2f}", size=9.5, weight=600, mono=True))
    b.append(text(x + W / 2, Y + 214, "3 × 3 cross-covariance", size=10, weight=600))
    b.append(text(x + W / 2, Y + 232, "H = UΣVᵀ", **CAPTION))
    # 4 R matrix
    x = xs[3]
    for i in range(3):
        for j in range(3):
            v = R[i][j]
            b.append(rect(x + 26 + j * 50, Y + 70 + i * 34, 46, 30, rx=4, fill="var(--sf-warm)", stroke="none", opacity=round(0.1 + abs(v) * 0.6, 2)))
            b.append(text(x + 49 + j * 50, Y + 89 + i * 34, f"{v:+.2f}", size=9.5, weight=600, mono=True))
    ang = math.degrees(math.atan2(R[1][0], R[0][0]))
    b.append(text(x + W / 2, Y + 194, f"d = sign det(VUᵀ) = {d:+.0f}", size=10, weight=600))
    b.append(text(x + W / 2, Y + 214, "R = V·diag(1, 1, d)·Uᵀ", size=10, weight=600))
    b.append(text(x + W / 2, Y + 232, f"recovers {ang:+.1f}° about z", **CAPTION))
    # 5 aligned
    x = xs[4]
    b.append(chain_2d(proj(Qc, x + 98, Y + 130, 30), "var(--sf-accent2)", label="1", r=8))
    b.append(chain_2d(proj(aligned, x + 98, Y + 130, 30), "var(--sf-rose)", dash="4 3", r=4))
    b.append(text(x + W / 2, Y + 194, f"RMSD = {rmsd:.3f} Å", size=12, fill="var(--sf-rose)", weight=800))
    b.append(text(x + W / 2, Y + 214, "√(1/N Σ |R P̄ᵢ − Q̄ᵢ|²)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 232, "only the noise remains", **CAPTION))
    write(slug, 1100, 324, "".join(b), label=f"Real Kabsch alignment of butane: centroids removed, cross-covariance H and rotation R from SVD recover the 40 degree rotation, leaving RMSD {rmsd:.3f} angstrom from the added noise.")


CELL = dict(a=5.0, b=6.0, c=7.0, alpha=80.0, beta=95.0, gamma=105.0)


def lattice_matrix(a, b_, c, al, be, ga):
    al, be, ga = map(math.radians, (al, be, ga))
    A = (a, 0.0, 0.0)
    B = (b_ * math.cos(ga), b_ * math.sin(ga), 0.0)
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0))
    return [A, B, (cx, cy, cz)]


def diagram_materials_unit_cell() -> None:
    slug = "materials-unit-cell"
    begin(slug)
    b = [element_defs(slug)]
    M = lattice_matrix(*CELL.values())
    V = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    b.append(text(450, 30, "Triclinic cell a = 5, b = 6, c = 7 Å · α = 80°, β = 95°, γ = 105°", size=15, weight=700))
    # oblique projection
    k = 24
    pr = lambda v: (150 + v[0] * k + v[2] * k * -0.45 * 0 + v[1] * k * 0.55 * 0 + v[1] * k * 0.0, 0)
    def P(v):
        x, y, z = v
        return (160 + k * (x + 0.55 * y * 0.0) + k * 0.42 * y, 300 - k * 0.32 * y - k * z)
    O = (0, 0, 0)
    a_, b_, c_ = M
    add = lambda *vs: tuple(sum(c) for c in zip(*vs))
    corners = {"O": O, "a": a_, "b": b_, "c": c_, "ab": add(a_, b_), "ac": add(a_, c_), "bc": add(b_, c_), "abc": add(a_, b_, c_)}
    edges = [("O", "a"), ("O", "b"), ("O", "c"), ("a", "ab"), ("a", "ac"), ("b", "ab"), ("b", "bc"), ("c", "ac"), ("c", "bc"), ("ab", "abc"), ("ac", "abc"), ("bc", "abc")]
    b.append(rect(24, 50, 420, 330, rx=14, sw=1.3))
    for u, v in edges:
        if u == "O":
            continue
        p1, p2 = P(corners[u]), P(corners[v])
        b.append(path(f"M{p1[0]:.1f},{p1[1]:.1f} L{p2[0]:.1f},{p2[1]:.1f}", sw=1.3, dash="4 3" if "b" in (u + v) and u in ("b",) else None))
    for name, col in [("a", "var(--sf-accent)"), ("b", "var(--sf-accent2)"), ("c", "var(--sf-violet)")]:
        p2 = P(corners[name]); p0 = P(O)
        b.append(path(f"M{p0[0]:.1f},{p0[1]:.1f} L{p2[0]:.1f},{p2[1]:.1f}", stroke=col, sw=2.4))
        b.append(circle(p2[0], p2[1], 4, fill=col))
        b.append(text(p2[0] + 10, p2[1] + 4, f"{name} = {CELL[name]:.0f} Å", size=10.5, fill=col, weight=700, anchor="start"))
    f = (0.5, 0.5, 0.5)
    r = tuple(sum(M[i][j] * f[i] for i in range(3)) for j in range(3))
    pc = P(r)
    b.append(circle(pc[0], pc[1], 6, fill="var(--sf-warm)", cls="sf-pulse"))
    b.append(text(pc[0] + 10, pc[1] - 8, "f = (½, ½, ½)", size=9.5, fill="var(--sf-warm)", anchor="start"))
    b.append(text(234, 366, "α: b∧c · β: a∧c · γ: a∧b", size=9.5, fill="var(--sf-muted)", weight=500))
    # right column: matrix + transforms
    x = 464
    b.append(rect(x, 50, 412, 330, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "lattice matrix M (rows a, b, c), Å", size=11.5, weight=700, anchor="start"))
    for i, row in enumerate(M):
        for j, v in enumerate(row):
            b.append(rect(x + 60 + j * 86, 88 + i * 30, 80, 26, rx=5, fill="var(--sf-card2)", sw=0.6))
            b.append(text(x + 100 + j * 86, 105 + i * 30, f"{v:+.3f}", size=10, weight=600, mono=True))
        b.append(text(x + 40, 105 + i * 30, "abc"[i], size=10.5, fill=["var(--sf-accent)", "var(--sf-accent2)", "var(--sf-violet)"][i], weight=700))
    b.append(text(x + 16, 200, f"V = det M = {V:.2f} Å³", size=11, fill="var(--sf-accent)", weight=700, anchor="start"))
    b.append(text(x + 16, 228, f"r = Mᵀ f = ({r[0]:.2f}, {r[1]:.2f}, {r[2]:.2f}) Å", size=10.5, fill="var(--sf-warm)", weight=600, anchor="start"))
    b.append(text(x + 16, 248, "f = (M⁻¹)ᵀ r  (inverse)", size=10, weight=500, anchor="start"))
    b.append(text(x + 16, 280, "wrap: f ← f − ⌊f⌋  e.g. (1.2, −0.3, 0.5) → (0.2, 0.7, 0.5)", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 16, 300, "minimum image: Δf ← Δf − round(Δf)", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 16, 320, "supercell n_a×n_b×n_c: f ← (f₀ + (i,j,k)) / n", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 16, 352, "a along x, b in the xy-plane, c fills z", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 396, "".join(b), label=f"A triclinic cell with a 5, b 6, c 7 angstrom and angles 80, 95, 105 degrees, its lattice matrix, volume {V:.1f} cubic angstrom and the fractional-to-Cartesian mapping of the cell center.")


def diagram_materials_assembly() -> None:
    slug = "materials-assembly"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "Assembling an sql framework: topology → metal nodes → oriented linkers → crystal", size=15, weight=700))
    W, H, Y = 254, 290, 54
    cards = [("1", "Topology sql", "var(--sf-accent)"), ("2", "Place metal nodes", "var(--sf-violet)"), ("3", "Orient linkers", "var(--sf-warm)"), ("4", "CrystalStructure", "var(--sf-accent2)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    def grid(x, y, s_, fn_node, fn_edge, n=2):
        out = ""
        for i in range(n + 1):
            for j in range(n + 1):
                cx, cy = x + i * s_, y + j * s_
                if i < n:
                    out += fn_edge(cx, cy, cx + s_, cy)
                if j < n:
                    out += fn_edge(cx, cy, cx, cy + s_)
        for i in range(n + 1):
            for j in range(n + 1):
                out += fn_node(x + i * s_, y + j * s_)
        return out
    # 1 topology
    x = xs[0]
    b.append(grid(x + 47, Y + 64, 80, lambda cx, cy: circle(cx, cy, 6, fill="var(--sf-accent)"), lambda x1, y1, x2, y2: path(f"M{x1},{y1} L{x2},{y2}", sw=1.6, dash="4 3")))
    b.append(rect(x + 47, Y + 64, 80, 80, rx=0, fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1.4, opacity=0.12))
    b.append(text(x + W / 2, Y + 252, "1 node · 2 edges per cell · CN 4", size=10, weight=600))
    b.append(text(x + W / 2, Y + 270, "shaded = one unit cell", **CAPTION))
    # 2 nodes: square-planar metal with 4 connection points
    x = xs[1]
    def metal(cx, cy):
        o = ""
        for dx, dy in [(14, 0), (-14, 0), (0, 14), (0, -14)]:
            o += path(f"M{cx},{cy} L{cx+dx},{cy+dy}", stroke="var(--sf-violet)", sw=2)
            o += circle(cx + dx, cy + dy, 3, fill="var(--sf-violet)", opacity=0.6)
        return o + circle(cx, cy, 8, fill="var(--sf-violet)") + text(cx, cy + 3, "M", size=7.5, fill="#fff", weight=700)
    b.append(grid(x + 47, Y + 64, 80, metal, lambda x1, y1, x2, y2: path(f"M{x1},{y1} L{x2},{y2}", sw=1, dash="2 5")))
    b.append(text(x + W / 2, Y + 252, "SBU at each vertex", size=10, weight=600))
    b.append(text(x + W / 2, Y + 270, "square-planar connection points", **CAPTION))
    # 3 linkers rotated (Rodrigues) along edges
    x = xs[2]
    def linker(x1, y1, x2, y2):
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        horiz = y1 == y2
        o = path(f"M{x1},{y1} L{x2},{y2}", stroke="var(--sf-ink)", sw=2)
        if horiz:
            o += f'<polygon points="{mx-12},{my} {mx-6},{my-10} {mx+6},{my-10} {mx+12},{my} {mx+6},{my+10} {mx-6},{my+10}" fill="var(--sf-card)" stroke="var(--sf-warm)" stroke-width="1.6"/>'
        else:
            o += f'<polygon points="{mx},{my-12} {mx+10},{my-6} {mx+10},{my+6} {mx},{my+12} {mx-10},{my+6} {mx-10},{my-6}" fill="var(--sf-card)" stroke="var(--sf-warm)" stroke-width="1.6"/>'
        return o
    b.append(grid(x + 47, Y + 64, 80, lambda cx, cy: circle(cx, cy, 8, fill="var(--sf-violet)"), linker))
    b.append(text(x + W / 2, Y + 252, "Rodrigues rotation: linker axis ∥ edge", size=10, weight=600))
    b.append(text(x + W / 2, Y + 270, "edge types may use different linkers", **CAPTION))
    # 4 crystal with cell box and fractional coords
    x = xs[3]
    b.append(grid(x + 47, Y + 64, 80, lambda cx, cy: circle(cx, cy, 8, fill="var(--sf-violet)"), linker))
    b.append(rect(x + 47, Y + 64, 80, 80, rx=0, fill="none", stroke="var(--sf-accent2)", sw=2.4, cls="sf-pulse"))
    b.append(text(x + 87, Y + 58, "a", size=10, fill="var(--sf-accent2)", weight=700))
    b.append(text(x + 40, Y + 108, "b", size=10, fill="var(--sf-accent2)", weight=700))
    b.append(text(x + W / 2, Y + 252, "all atoms in fractional coordinates", size=10, weight=600))
    b.append(text(x + W / 2, Y + 270, "+ lattice → periodic 2×2 view", **CAPTION))
    write(slug, 1100, 364, "".join(b), label="Framework assembly on the sql net: one node and two edges per cell, square-planar metal nodes at each vertex, linkers rotated onto each edge, and the periodic crystal with its unit cell.")


def diagram_transport_arrow() -> None:
    slug = "transport-arrow"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "pack_conformers: water rows → contiguous typed columns", size=15, weight=700))
    rows = [(x, y, z, -0.0, el) for (x, y, z), el in zip(WATER_XYZ, [8, 1, 1])]
    # left: JSON rows
    b.append(rect(24, 50, 330, 300, rx=14, sw=1.3))
    b.append(text(40, 74, "row-oriented (JSON)", size=11.5, weight=700, anchor="start"))
    for i, (x, y, z, e, el) in enumerate(rows):
        yy = 100 + i * 66
        b.append(rect(40, yy, 298, 54, rx=8, fill="var(--sf-card2)", sw=0.8))
        b.append(text(52, yy + 22, f'{{"x":{x:.3f},"y":{y:.3f},', size=9.5, weight=500, anchor="start", mono=True))
        b.append(text(52, yy + 40, f' "z":{z:.1f},"energy":…,"element":{el}}}', size=9.5, weight=500, anchor="start", mono=True))
    b.append(text(189, 316, "parse every field, objects scattered in heap", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(189, 334, "GPU upload needs copy + repack", size=9.5, fill="var(--sf-rose)", weight=600))
    b.append(arrow_between(360, 200, 396, 200, color="var(--sf-accent)"))
    # right: columns
    x = 404
    b.append(rect(x, 50, 472, 300, rx=14, sw=1.3))
    b.append(text(x + 16, 74, "RecordBatch, num_rows = 3", size=11.5, weight=700, anchor="start"))
    cols = [("x", "Float64", [r[0] for r in rows], "var(--sf-accent)"), ("y", "Float64", [r[1] for r in rows], "var(--sf-accent2)"),
            ("z", "Float64", [r[2] for r in rows], "var(--sf-violet)"), ("energy", "Float64", ["E", "E", "E"], "var(--sf-warm)"), ("element", "Uint8", [8, 1, 1], "var(--sf-rose)")]
    for i, (name, typ, vals, col) in enumerate(cols):
        yy = 92 + i * 44
        b.append(text(x + 16, yy + 20, name, size=10.5, fill=col, weight=700, anchor="start"))
        b.append(text(x + 16, yy + 34, typ, size=8.5, fill="var(--sf-muted)", weight=500, anchor="start"))
        bw = 84 if typ == "Float64" else 40
        for j, v in enumerate(vals):
            b.append(rect(x + 92 + j * (bw + 2), yy + 6, bw, 28, rx=3, fill=col, stroke="none", opacity=0.25))
            b.append(text(x + 92 + j * (bw + 2) + bw / 2, yy + 24, f"{v:+.3f}" if isinstance(v, float) else str(v), size=9.5, weight=600, mono=True))
        nbytes = 8 * 3 if typ == "Float64" else 3
        b.append(text(x + 456, yy + 24, f"{nbytes} B", size=9, fill="var(--sf-muted)", anchor="end", mono=True))
    b.append(text(x + 236, 328, "each column = one contiguous buffer → Float64Array / Uint8Array, zero-copy", size=9.5, fill="var(--sf-accent2)", weight=600))
    write(slug, 900, 366, "".join(b), label="Water conformer rows as JSON objects versus the RecordBatch produced by pack_conformers: five contiguous typed columns x, y, z, energy (Float64) and element (Uint8).")


def diagram_transport_chunked() -> None:
    slug = "transport-chunked"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Chunked streaming and worker split/merge", size=15, weight=700))
    # chunking of water ESP grid: 2912 values, chunk 1000
    b.append(rect(24, 50, 852, 130, rx=14, sw=1.3))
    b.append(text(40, 74, "chunk_esp_grid(water grid, 1000): N_chunks = ⌈2912 / 1000⌉ = 3", size=11.5, weight=700, anchor="start"))
    total, cs_ = 2912, 1000
    x0, w = 40, 820
    for k, (n, col) in enumerate([(1000, "var(--sf-accent)"), (1000, "var(--sf-accent2)"), (912, "var(--sf-violet)")]):
        xs_ = x0 + k * cs_ / total * w
        ww = n / total * w
        b.append(rect(xs_ + 1, 92, ww - 2, 34, rx=6, fill=col, stroke="none", opacity=0.3, cls=f"sf-pulse sf-d{k+1}"))
        b.append(text(xs_ + ww / 2, 114, f"index {k} · total 3 · count {n}", size=9.5, weight=600, mono=True))
    b.append(text(450, 150, "lazy ChunkedIterator: one slice alive at a time → GPU or Web Worker", size=9.5, fill="var(--sf-muted)", weight=500))
    b.append(text(450, 168, "coordinates chunk by atoms: 1000 atoms = 3000 f64 values", size=9.5, fill="var(--sf-muted)", weight=500))
    # split / merge
    b.append(rect(24, 194, 852, 230, rx=14, sw=1.3))
    b.append(text(40, 218, "split_batch: 200 SMILES, max_workers = 8 → n_workers = min(8, ⌈200/64⌉) = 4", size=11.5, weight=700, anchor="start"))
    wc = ["var(--sf-accent)", "var(--sf-accent2)", "var(--sf-violet)", "var(--sf-warm)"]
    for i in range(16):
        b.append(rect(40 + i * 24, 238, 20, 20, rx=4, fill=wc[i % 4], stroke="none", opacity=0.75))
        b.append(text(50 + i * 24, 252, str(i), size=8, fill="#fff", weight=700))
    b.append(text(434, 252, "… item i → worker i mod 4 (round-robin)", size=9.5, fill="var(--sf-muted)", weight=500, anchor="start"))
    for k in range(4):
        cx = 120 + k * 160
        b.append(rect(cx - 60, 284, 120, 64, rx=10, fill="var(--sf-card2)", stroke=wc[k], sw=1.4))
        b.append(text(cx, 306, f"worker {k}", size=10.5, fill=wc[k], weight=700))
        b.append(text(cx, 324, f"items {k}, {k+4}, {k+8}, …", size=9, weight=500, mono=True))
        b.append(text(cx, 340, "50 items", size=9, fill="var(--sf-muted)", weight=450))
        b.append(connector(f"M{cx},348 L{cx},{372 + k * 8} L736,{372 + k * 8}", color=wc[k], delay=0.3 * k))
    b.append(rect(740, 284, 120, 120, rx=10, fill="var(--sf-soft)", stroke="var(--sf-accent2)", sw=1.4))
    b.append(text(800, 316, "merge_results", size=10.5, fill="var(--sf-accent2)", weight=700))
    b.append(text(800, 336, "sort by task_id", size=9.5, weight=500))
    b.append(text(800, 354, "→ original order", size=9.5, weight=500))
    b.append(text(800, 384, "0, 1, 2, 3, …", size=9.5, fill="var(--sf-muted)", weight=500, mono=True))
    write(slug, 900, 440, "".join(b), label="The water ESP grid of 2912 values streamed as three chunks of 1000, 1000 and 912, and 200 SMILES split round-robin over four workers then merged back by task id.")


# Real sci-form outputs used by the spectroscopy diagrams.
# Benzaldehyde (embed seed 42) → compute_stda_uvvis(σ = 0.3 eV): (energy eV, f).
BENZALDEHYDE_STDA = [
    (5.612, 0.2288), (5.29, 0.2235), (6.431, 0.1438), (6.267, 0.132), (3.813, 0.1229), (6.593, 0.1121), (7.026, 0.0978),
    (6.93, 0.0945), (6.561, 0.09), (7.548, 0.0887), (6.235, 0.0778), (6.367, 0.074), (4.039, 0.0688), (6.397, 0.0688),
    (6.139, 0.0664), (0.69, 0.0552), (7.128, 0.0522), (6.757, 0.0509), (1.165, 0.0497), (6.365, 0.0486), (7.356, 0.0482),
    (7.621, 0.0474), (6.789, 0.0467), (1.248, 0.0414), (1.292, 0.0405), (7.099, 0.0376), (7.16, 0.0369), (6.463, 0.0366),
    (6.952, 0.0363), (7.156, 0.0353), (7.229, 0.03), (7.089, 0.0285), (4.43, 0.0214), (7.919, 0.0209), (7.992, 0.0205),
    (7.003, 0.0174), (5.516, 0.0171), (3.522, 0.0146), (5.18, 0.0113), (5.908, 0.0102), (5.319, 0.0099), (5.387, 0.0087),
    (1.375, 0.0084), (0.607, 0.0084), (7.506, 0.0075), (7.538, 0.0071), (7.324, 0.0066), (3.909, 0.0048),
]
# Water (embed seed 42) → compute_vibrational_analysis(method = "xtb"): (cm⁻¹, km/mol).
WATER_IR_XTB = [(1586.9, 90.5), (3204.5, 49.83), (3352.2, 56.54)]
WATER_ZPVE_XTB = 0.505
# Ethanol → predict_nmr_shifts / predict_nmr_couplings (topological), 400 MHz ¹H.
ETHANOL_H = [("CH₃", 0.90, 3, "t", [1, 2, 1]), ("CH₂", 2.01, 2, "quint", [1, 4, 6, 4, 1]), ("OH", 2.50, 1, "t", [1, 2, 1])]
KARPLUS = (7.76, 1.10, 1.40)


def gauss(x: float, s_: float) -> float:
    return math.exp(-x * x / (2 * s_ * s_)) / (s_ * math.sqrt(2 * math.pi))


def lorentz(x: float, g: float) -> float:
    return (g / math.pi) / (x * x + g * g)


def spectrum_plot(b: list, x: float, y: float, w: float, h: float, sticks: list, xr: tuple, fn, *, color: str, stick_color: str, xlabel: str, ticks: list, reverse: bool = False, ymax: float | None = None) -> tuple:
    sx = (lambda v: x + (xr[1] - v) / (xr[1] - xr[0]) * w) if reverse else (lambda v: x + (v - xr[0]) / (xr[1] - xr[0]) * w)
    vals = [fn(xr[0] + (xr[1] - xr[0]) * i / 400) for i in range(401)]
    top = ymax or max(vals) * 1.08
    sy = lambda v: y + h - v / top * h
    smax = max(si for _, si in sticks)
    for pos, inten in sticks:
        if xr[0] <= pos <= xr[1]:
            b.append(path(f"M{sx(pos):.1f},{y+h} L{sx(pos):.1f},{y + h - inten / smax * h * 0.92:.1f}", stroke=stick_color, sw=1.3, opacity=0.7))
    b.append(path("M" + " L".join(f"{sx(xr[0] + (xr[1] - xr[0]) * i / 400):.1f},{sy(v):.1f}" for i, v in enumerate(vals)), stroke=color, sw=2.2))
    b.append(path(f"M{x},{y+h} L{x+w},{y+h}", sw=1))
    for t in ticks:
        b.append(text(sx(t), y + h + 13, f"{t:g}", size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(x + w, y + h + 26, xlabel, size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    return sx, sy


def diagram_uvvis_stda_pipeline() -> None:
    slug = "uvvis-stda-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "sTDA UV-Vis of benzaldehyde (c1ccc(cc1)C=O), real sci-form output", size=15, weight=700))
    W, H, Y = 254, 290, 54
    cards = [("1", "EHT orbitals", "var(--sf-accent)"), ("2", "Enumerate i → a", "var(--sf-violet)"), ("3", "Oscillator strengths", "var(--sf-warm)"), ("4", "Broadened ε(λ)", "var(--sf-accent2)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 levels around the gap
    x = xs[0]
    lv = [-13.4, -12.48, -12.01, -11.93, -11.32, -10.76, -10.64]
    slots = [Y + 214, Y + 184, Y + 160, Y + 136, Y + 104, Y + 80, Y + 62]  # schematic spacing
    for i, e in enumerate(lv):
        yy = slots[i]
        occ = i < 4
        col = "var(--sf-violet)" if occ else "var(--sf-rose)"
        b.append(path(f"M{x+70},{yy:.1f} L{x+150},{yy:.1f}", stroke=col, sw=2.2))
        b.append(text(x + 64, yy + 3, f"{e:.2f}", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    b.append(text(x + 156, slots[3] + 3, "HOMO 19", size=9, fill="var(--sf-violet)", anchor="start"))
    b.append(text(x + 156, slots[4] + 3, "LUMO 20", size=9, fill="var(--sf-rose)", anchor="start"))
    b.append(text(x + W / 2, Y + 250, "40 valence e⁻ · gap 0.61 eV", size=10, weight=600))
    b.append(text(x + W / 2, Y + 268, "εᵢ in eV, frontier window", **CAPTION))
    # 2 transitions
    x = xs[1]
    b.append(text(x + W / 2, Y + 64, "ΔE_ia = ε_a − ε_i", size=11, weight=700))
    b.append(text(x + W / 2, Y + 88, "λ_ia = 1239.84 / ΔE_ia  nm", size=10, weight=500))
    for k, (fr, to, e) in enumerate([(16, 24, 3.813), (16, 25, 5.29), (8, 24, 5.612)]):
        yy = Y + 120 + k * 34
        b.append(rect(x + 20, yy, 214, 28, rx=7, fill="var(--sf-card2)", sw=0.8))
        b.append(text(x + 32, yy + 18, f"{fr} → {to}", size=10, weight=700, anchor="start", mono=True))
        b.append(text(x + 222, yy + 18, f"{e:.2f} eV = {1239.84 / e:.0f} nm", size=10, weight=500, anchor="end"))
    b.append(text(x + W / 2, Y + 250, "50 transitions kept (1–8 eV window)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 268, "every occupied → virtual pair", **CAPTION))
    # 3 stick spectrum
    x = xs[2]
    sx = lambda e: x + 20 + (e - 0.5) / 7.7 * 214
    for e, f in BENZALDEHYDE_STDA:
        b.append(path(f"M{sx(e):.1f},{Y+214} L{sx(e):.1f},{Y + 214 - f / 0.2288 * 140:.1f}", stroke="var(--sf-warm)", sw=1.6))
    b.append(path(f"M{x+20},{Y+214} L{x+234},{Y+214}", sw=1))
    for t in (1, 3, 5, 7):
        b.append(text(sx(t), Y + 227, str(t), size=8.5, fill="var(--sf-muted)", weight=500))
    b.append(text(sx(5.612), Y + 66, "f = 0.229", size=9, fill="var(--sf-warm)"))
    b.append(text(x + W / 2, Y + 250, "f_ia ∝ ΔE_ia |μ_ia|²", size=10, weight=600))
    b.append(text(x + W / 2, Y + 268, "one-center AO transition dipoles", **CAPTION))
    # 4 broadened in wavelength
    x = xs[3]
    def eps_nm(lam):
        E = 1239.84 / lam
        return sum(f * gauss(E - e, 0.3) for e, f in BENZALDEHYDE_STDA)
    spectrum_plot(b, x + 20, Y + 64, 214, 150, [(1239.84 / e, f) for e, f in BENZALDEHYDE_STDA], (160, 420), eps_nm, color="var(--sf-accent2)", stick_color="var(--sf-warm)", xlabel="λ (nm)", ticks=[200, 300, 400])
    b.append(text(x + W / 2, Y + 250, "Gaussian σ = 0.3 eV on the energy axis", size=10, weight=600))
    b.append(text(x + W / 2, Y + 268, "maxima ≈ 191, 225, 319 nm (overlapping lines)", **CAPTION))
    write(slug, 1100, 364, "".join(b), label="Real sTDA output for benzaldehyde: EHT frontier orbitals with a 0.61 eV gap, enumerated occupied-to-virtual transitions, oscillator-strength sticks and the Gaussian-broadened absorption with maxima near 191, 225 and 319 nm.")


def diagram_spectroscopy_overview() -> None:
    slug = "spectroscopy-overview"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "One embedded geometry, three spectroscopy tracks (real sci-form outputs)", size=15, weight=700))
    b.append(rect(24, 60, 160, 330, rx=14, sw=1.3))
    b.append(text(104, 84, "embed(smiles, seed)", size=10, weight=700, mono=True))
    b.append(text(104, 102, "elements + 3D coords", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append('<g class="sf-float">' + water_2d(slug, 104, 230, 46, r=12) + "</g>")
    tracks = [
        ("UV-Vis · sTDA", "benzaldehyde · EHT MOs → i→a → f", "var(--sf-accent)", 60),
        ("IR · numerical Hessian", "water · GFN2-xTB · 3N − 6 = 3 modes", "var(--sf-warm)", 172),
        ("¹H NMR · HOSE + Karplus", "ethanol · 400 MHz · topological J", "var(--sf-accent2)", 284),
    ]
    for name, sub, col, y in tracks:
        b.append(connector(f"M188,225 C230,225 230,{y+53} 268,{y+53}", color=col, delay=y / 200))
        b.append(rect(272, y, 604, 106, rx=14, sw=1.3))
        b.append(rect(273.5, y + 1.5, 4, 103, rx=2, fill=col, stroke="none"))
        b.append(text(290, y + 24, name, size=11.5, fill=col, weight=700, anchor="start"))
        b.append(text(290, y + 42, sub, size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    # UV
    def eps_e(E):
        return sum(f * gauss(E - e, 0.3) for e, f in BENZALDEHYDE_STDA)
    spectrum_plot(b, 520, 70, 340, 66, BENZALDEHYDE_STDA, (0.5, 8.2), eps_e, color="var(--sf-accent)", stick_color="var(--sf-line)", xlabel="E (eV)", ticks=[1, 3, 5, 7])
    # IR (wavenumber decreasing to the right, as plotted by convention)
    def ir(nu):
        return sum(i * lorentz(nu - f, 20) for f, i in WATER_IR_XTB)
    spectrum_plot(b, 520, 182, 340, 66, WATER_IR_XTB, (1000, 4000), ir, color="var(--sf-warm)", stick_color="var(--sf-line)", xlabel="ν̃ (cm⁻¹)", ticks=[1000, 2000, 3000, 4000], reverse=True)
    b.append(text(300, 242, "1587 bend · 3205 / 3352 stretch", size=9.5, weight=600, anchor="start"))
    b.append(text(300, 260, "ZPVE 0.505 eV", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(300, 130, "bands 5.61 · 5.29 · 3.81 eV", size=9.5, weight=600, anchor="start"))
    b.append(text(300, 148, "(221 · 234 · 325 nm)", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    # NMR
    lines = []
    for _, d0, n, _, w8 in ETHANOL_H:
        J = 7.0 / 400
        m = len(w8)
        tot = sum(w8)
        for i, wi in enumerate(w8):
            lines.append((d0 + (i - (m - 1) / 2) * J, n * wi / tot))
    def nmr(d):
        return sum(wi * lorentz(d - di, 0.004) for di, wi in lines)
    spectrum_plot(b, 520, 294, 340, 66, [(d, w8) for d, w8 in lines], (0.0, 3.5), nmr, color="var(--sf-accent2)", stick_color="var(--sf-line)", xlabel="δ (ppm)", ticks=[0, 1, 2, 3], reverse=True)
    b.append(text(300, 354, "CH₃ 0.90 t · CH₂ 2.01 quint · OH 2.50 t", size=9.5, weight=600, anchor="start"))
    b.append(text(300, 372, "³J = 7.0 Hz", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 406, "".join(b), label="One embedded geometry feeding three tracks with real sci-form outputs: sTDA UV-Vis of benzaldehyde, xTB infrared spectrum of water and the predicted 1H NMR spectrum of ethanol.")


def diagram_broadening_functions() -> None:
    slug = "broadening-functions"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Gaussian vs. Lorentzian line shape, both with width parameter 0.3 eV", size=15, weight=700))
    energy_panel(b, 24, 50, 520, 320, "line shapes centred at x = 0 (area = 1)", "G = exp(−x²/2σ²)/(σ√2π)     L = (γ/π)/(x² + γ²)",
                 [(lambda x: gauss(x, 0.3), "var(--sf-accent)", 2.4, None), (lambda x: lorentz(x, 0.3), "var(--sf-rose)", 2.4, None)],
                 (-1.6, 1.6), (0, 1.45), xlabel="x (eV)", marks=[(-0.353, "", "var(--sf-accent)"), (0.353, "FWHM G 0.71", "var(--sf-accent)"), (0.3, "", "var(--sf-rose)"), (-0.3, "FWHM L 0.60", "var(--sf-rose)")],
                 note="peak heights: G 1.33 · L 1.06 eV⁻¹")
    # log tail panel
    x = 564
    b.append(rect(x, 50, 312, 320, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "tails on a log scale", size=11.5, weight=700, anchor="start"))
    px, py, pw, ph = x + 40, 92, 250, 200
    sx = lambda v: px + v / 2.0 * pw
    sy = lambda v: py + ph - (math.log10(max(v, 1e-6)) + 6) / 6.3 * ph
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    for e_ in (0, -2, -4, -6):
        b.append(text(px - 6, sy(10 ** e_) + 3, f"10^{e_}" if e_ else "1", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    for fn, col in [(lambda v: gauss(v, 0.3), "var(--sf-accent)"), (lambda v: lorentz(v, 0.3), "var(--sf-rose)")]:
        pts = [(sx(i / 100), sy(fn(i / 100))) for i in range(201) if fn(i / 100) > 1e-6]
        b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in pts), stroke=col, sw=2.2))
    for t in (0, 0.5, 1, 1.5, 2):
        b.append(text(sx(t), py + ph + 13, f"{t:g}", size=8.5, fill="var(--sf-muted)", weight=500))
    g1, l1 = gauss(1.0, 0.3), lorentz(1.0, 0.3)
    b.append(circle(sx(1), sy(g1), 4, fill="var(--sf-accent)"))
    b.append(circle(sx(1), sy(l1), 4, fill="var(--sf-rose)"))
    b.append(text(sx(1) + 8, sy(l1) - 6, f"L(1 eV) = {l1:.3f}", size=9, fill="var(--sf-rose)", anchor="start"))
    b.append(text(sx(1) + 8, sy(g1) + 12, f"G(1 eV) = {g1:.4f}", size=9, fill="var(--sf-accent)", anchor="start"))
    b.append(text(x + 156, 334, f"at 1 eV the Lorentzian is {l1/g1:.0f}× larger:", size=9.5, weight=600))
    b.append(text(x + 156, 352, "distant bands leak into each other", size=9.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 386, "".join(b), label="Gaussian and Lorentzian line shapes with width 0.3 eV: the Gaussian is taller with FWHM 0.71 eV, the Lorentzian has FWHM 0.60 eV and, at 1 eV from the centre, a tail about twenty times larger.")


def diagram_ir_hessian_pipeline() -> None:
    slug = "ir-hessian-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "IR from a numerical Hessian — water with GFN2-xTB (real sci-form output)", size=15, weight=700))
    W, H, Y = 196, 250, 54
    cards = [("1", "Displace ±δ", "var(--sf-accent)"), ("2", "Hessian 9 × 9", "var(--sf-accent2)"), ("3", "Mass-weight + diag", "var(--sf-violet)"), ("4", "Normal modes", "var(--sf-warm)"), ("5", "Spectrum", "var(--sf-rose)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 displacements
    x = xs[0]
    b.append(water_2d(slug, x + 98, Y + 118, 60, r=12))
    hx, hy = x + 98 + WATER_XYZ[1][0] * 60, Y + 118 - WATER_XYZ[1][1] * 60
    for dx in (-14, 14):
        b.append(circle(hx + dx, hy, 6, fill="none", stroke="var(--sf-accent)", sw=1.3))
    b.append(path(f"M{hx-14},{hy-16} L{hx+14},{hy-16}", stroke="var(--sf-accent)", sw=1.4))
    b.append(text(hx, hy - 22, "±δ", size=9, fill="var(--sf-accent)"))
    b.append(text(x + W / 2, Y + 212, "E(q ± δᵢ ± δⱼ), δ = 0.005 Å", size=10, weight=600))
    b.append(text(x + W / 2, Y + 230, "3N = 9 Cartesian coordinates", **CAPTION))
    # 2 hessian blocks
    x = xs[1]; mx, my, cs = x + 26, Y + 56, 16
    lab = ["Ox", "Oy", "Oz", "Hx", "Hy", "Hz", "Hx", "Hy", "Hz"]
    for i in range(9):
        for j in range(9):
            same = i // 3 == j // 3
            zrow = (i % 3 == 2) or (j % 3 == 2)
            op = 0.08 if zrow and not (i == j) else (0.75 if i == j else (0.45 if same else 0.28))
            b.append(rect(mx + j * cs, my + i * cs, cs - 1.5, cs - 1.5, rx=2, fill="var(--sf-accent2)", stroke="none", opacity=op))
    b.append(text(x + W / 2, Y + 212, "Hᵢⱼ = ∂²E / ∂qᵢ∂qⱼ (central diff.)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 230, "z rows ≈ 0: planar molecule", **CAPTION))
    # 3 mass weighting
    x = xs[2]
    b.append(text(x + W / 2, Y + 70, "H̃ᵢⱼ = Hᵢⱼ / √(mᵢ mⱼ)", size=11, weight=700))
    b.append(text(x + W / 2, Y + 96, "m_O = 15.999 · m_H = 1.008", size=9.5, weight=500))
    b.append(text(x + W / 2, Y + 126, "H̃ L = Λ L", size=11, weight=700))
    b.append(text(x + W / 2, Y + 152, "ν̃ₖ = √(F λₖ) / 2πc", size=10.5, weight=600))
    b.append(text(x + W / 2, Y + 174, "6 modes |ν̃| < 10 cm⁻¹ dropped", size=9.5, fill="var(--sf-rose)", weight=600))
    b.append(text(x + W / 2, Y + 212, "3 translations + 3 rotations", size=10, weight=600))
    b.append(text(x + W / 2, Y + 230, "negative λ → imaginary ν̃", **CAPTION))
    # 4 modes drawn with arrows
    x = xs[3]
    modes = [("1587", "bend", [(0, 0.25), (0.5, -0.4), (-0.5, -0.4)]), ("3205", "sym. stretch", [(0, 0.15), (-0.6, -0.45), (0.6, -0.45)]), ("3352", "asym. stretch", [(0.2, 0), (-0.7, -0.5), (-0.7, 0.4)])]
    for k, (f, name, disp) in enumerate(modes):
        cx, cy = x + 40 + k * 60, Y + 110
        pts = [(cx + px * 26, cy - py * 26) for px, py, _ in WATER_XYZ]
        for hx_, hy_ in pts[1:]:
            b.append(path(f"M{pts[0][0]:.1f},{pts[0][1]:.1f} L{hx_:.1f},{hy_:.1f}", stroke="var(--sf-ink)", sw=1.6))
        for (ax, ay), (dx, dy) in zip(pts, disp):
            b.append(path(f"M{ax:.1f},{ay:.1f} l{dx*22:.1f},{-dy*22:.1f}", stroke="var(--sf-warm)", sw=1.8))
        b.append(atom_ball(slug, *pts[0], "O", 6, label=False))
        for hx_, hy_ in pts[1:]:
            b.append(atom_ball(slug, hx_, hy_, "H", 4))
        b.append(text(cx, cy + 46, f"{f}", size=9.5, fill="var(--sf-warm)", weight=700))
        b.append(text(cx, cy + 60, name, size=8, fill="var(--sf-muted)", weight=500))
    b.append(text(x + W / 2, Y + 212, "3N − 6 = 3 real modes (cm⁻¹)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 230, f"ZPVE = ½ Σ hcν̃ = {WATER_ZPVE_XTB} eV", **CAPTION))
    # 5 spectrum
    x = xs[4]
    def ir(nu):
        return sum(i * lorentz(nu - f, 20) for f, i in WATER_IR_XTB)
    spectrum_plot(b, x + 16, Y + 60, 166, 120, WATER_IR_XTB, (1000, 4000), ir, color="var(--sf-rose)", stick_color="var(--sf-warm)", xlabel="ν̃ (cm⁻¹)", ticks=[1000, 2500, 4000], reverse=True)
    b.append(text(x + W / 2, Y + 212, "Iₖ ∝ |∂μ/∂Qₖ|² · Lorentz γ = 20", size=10, weight=600))
    b.append(text(x + W / 2, Y + 230, "90.5 · 49.8 · 56.5 km/mol", **CAPTION))
    write(slug, 1100, 324, "".join(b), label="Numerical-Hessian IR for water with xTB: finite displacements, the 9 by 9 Hessian, mass weighting and diagonalization with six near-zero modes removed, three normal modes at 1587, 3205 and 3352 per centimetre, and the Lorentzian spectrum.")


def diagram_nmr_pipeline() -> None:
    slug = "nmr-pipeline"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "¹H NMR of ethanol (CCO) — real sci-form prediction at 400 MHz", size=15, weight=700))
    W, H, Y = 254, 300, 54
    cards = [("1", "HOSE environments", "var(--sf-accent)"), ("2", "Shifts δ", "var(--sf-violet)"), ("3", "Couplings J", "var(--sf-warm)"), ("4", "Multiplets + Lorentz", "var(--sf-accent2)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 HOSE spheres around CH2 carbon
    x = xs[0]; cx, cy = x + 110, Y + 130
    for r_, col in [(70, "var(--sf-accent)"), (40, "var(--sf-violet)")]:
        b.append(circle(cx, cy, r_, fill="none", stroke=col, sw=1.2, opacity=0.5))
    C1, C2, O = (cx - 44, cy + 20), (cx, cy), (cx + 44, cy + 20)
    b.append(path(f"M{C1[0]},{C1[1]} L{C2[0]},{C2[1]} L{O[0]},{O[1]}", stroke="var(--sf-ink)", sw=2.4))
    b.append(atom_ball(slug, *C1, "C", 10)); b.append(atom_ball(slug, *C2, "C", 10)); b.append(atom_ball(slug, *O, "O", 10))
    b.append(text(cx, cy - 76, "sphere 2", size=8.5, fill="var(--sf-accent)"))
    b.append(text(cx, cy - 46, "sphere 1", size=8.5, fill="var(--sf-violet)"))
    b.append(text(x + W / 2, Y + 236, "BFS to radius 4 from each atom", size=10, weight=600))
    b.append(text(x + W / 2, Y + 254, "CH₂ carbon: neighbors C and O", **CAPTION))
    # 2 shifts table
    x = xs[1]
    rows = [("H (CH₃) ×3", "0.90", "sp3_methyl"), ("H (CH₂) ×2", "2.01", "sp3_methylene"), ("H (OH)", "2.50", "alcohol_O-H"), ("C1 (CH₃)", "20.0", "methyl_C"), ("C2 (CH₂)", "36.9", "methylene_C")]
    for k, (a_, d_, env) in enumerate(rows):
        yy = Y + 64 + k * 34
        b.append(rect(x + 14, yy - 16, 226, 28, rx=7, fill="var(--sf-card2)", sw=0.8))
        b.append(text(x + 24, yy + 2, a_, size=10, weight=600, anchor="start"))
        b.append(text(x + 150, yy + 2, d_, size=10.5, fill="var(--sf-violet)", weight=700, anchor="end", mono=True))
        b.append(text(x + 158, yy + 2, env, size=8, fill="var(--sf-muted)", weight=500, anchor="start", mono=True))
    b.append(text(x + W / 2, Y + 236, "base + inductive + group corrections", size=10, weight=600))
    b.append(text(x + W / 2, Y + 254, "ppm · ¹H and ¹³C", **CAPTION))
    # 3 couplings
    x = xs[2]
    items = [("²J geminal", "12.0 Hz", "same carbon (equivalent: no split)"), ("³J vicinal", "7.0 Hz", "CH₃–CH₂ and CH₂–OH"), ("⁴J long-range", "1.0 Hz", "CH₃···OH (not resolved)")]
    for k, (t_, v, d_) in enumerate(items):
        yy = Y + 64 + k * 52
        b.append(text(x + 18, yy, t_, size=10.5, fill="var(--sf-warm)", weight=700, anchor="start"))
        b.append(text(x + 236, yy, v, size=10.5, weight=700, anchor="end", mono=True))
        b.append(text(x + 18, yy + 18, d_, size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + W / 2, Y + 236, "no 3D coords → free-rotation ³J", size=10, weight=600))
    b.append(text(x + W / 2, Y + 254, "7 Hz = 0.0175 ppm at 400 MHz", **CAPTION))
    # 4 spectrum with insets
    x = xs[3]
    lines = []
    J = 7.0 / 400
    for _, d0, n, _, w8 in ETHANOL_H:
        m = len(w8); tot = sum(w8)
        for i, wi in enumerate(w8):
            lines.append((d0 + (i - (m - 1) / 2) * J, n * wi / tot))
    def nmr(d):
        return sum(wi * lorentz(d - di, 0.004) for di, wi in lines)
    sx, _ = spectrum_plot(b, x + 18, Y + 70, 218, 100, lines, (0.4, 2.9), nmr, color="var(--sf-accent2)", stick_color="var(--sf-line)", xlabel="δ (ppm)", ticks=[0.5, 1.0, 1.5, 2.0, 2.5], reverse=True)
    for name, d0, n, mult, w8 in ETHANOL_H:
        b.append(text(sx(d0), Y + 62, f"{name} {mult}", size=8.5, fill="var(--sf-accent2)", weight=700))
    b.append(text(x + W / 2, Y + 214, "areas 3 : 2 : 1", size=10, weight=600))
    b.append(text(x + W / 2, Y + 236, "Pascal weights 1:2:1 · 1:4:6:4:1", size=10, weight=600))
    b.append(text(x + W / 2, Y + 254, "γ = 0.004 ppm here (default 0.01)", **CAPTION))
    write(slug, 1100, 374, "".join(b), label="Ethanol 1H NMR from sci-form: HOSE spheres, predicted shifts 0.90, 2.01 and 2.50 ppm, couplings of 12, 7 and 1 Hz, and the Lorentzian spectrum with a triplet, quintet and triplet in a 3:2:1 ratio.")


def diagram_nmr_karplus() -> None:
    slug = "nmr-karplus"
    begin(slug)
    b = [element_defs(slug)]
    A, B, C = KARPLUS
    J = lambda ph: A * math.cos(math.radians(ph)) ** 2 - B * math.cos(math.radians(ph)) + C
    b.append(text(450, 30, "³J(H,H) = 7.76 cos²φ − 1.10 cos φ + 1.40 Hz  (Altona–Sundaralingam)", size=15, weight=700))
    px, py, pw, ph = 80, 150, 740, 190
    sx = lambda a: px + a / 360 * pw
    sy = lambda v: py + ph - v / 11 * ph
    b.append(axes(px, py, pw, ph, xlabel="H–C–C–H dihedral φ (°)", ylabel="³J (Hz)"))
    b.append(curve_path([(sx(a), sy(J(a))) for a in range(0, 361, 5)], color="var(--sf-accent2)", sw=2.6, cls=""))
    avg = A / 2 + C
    b.append(path(f"M{px},{sy(avg):.1f} L{px+pw},{sy(avg):.1f}", stroke="var(--sf-warm)", sw=1.4, dash="5 4"))
    b.append(text(px + pw - 4, sy(avg) - 6, f"free-rotation average A/2 + C = {avg:.2f} Hz", size=9.5, fill="var(--sf-warm)", anchor="end"))
    for v in (0, 2, 4, 6, 8, 10):
        b.append(text(px - 8, sy(v) + 3, str(v), size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    for a in range(0, 361, 60):
        b.append(text(sx(a), py + ph + 13, str(a), size=8.5, fill="var(--sf-muted)", weight=500))
    for a, lbl in [(0, "eclipsed"), (60, "gauche"), (90, "orthogonal"), (180, "anti")]:
        b.append(circle(sx(a), sy(J(a)), 5, fill="var(--sf-accent2)", cls="sf-pulse"))
        b.append(text(sx(a) + (10 if a < 180 else 0), sy(J(a)) - 10, f"{J(a):.2f} Hz", size=9.5, fill="var(--sf-accent2)", weight=700, anchor="start" if a < 180 else "middle"))
        b.append(newman(sx(a) + (26 if a == 0 else (-14 if a == 60 else (14 if a == 90 else 0))), 86, a, r=15))
        b.append(text(sx(a) + (26 if a == 0 else (-14 if a == 60 else (14 if a == 90 else 0))), 126, f"{lbl} {a}°", size=9, fill="var(--sf-muted)"))
    b.append(text(450, 386, "with 3D coordinates the actual dihedral is used; topology only → the 5.28 Hz average", size=9.5, fill="var(--sf-muted)", weight=500))
    write(slug, 900, 400, "".join(b), label="Karplus curve 7.76 cos² φ − 1.10 cos φ + 1.40: 8.06 Hz eclipsed, 2.79 Hz gauche, 1.40 Hz at 90 degrees and 10.26 Hz anti, with the 5.28 Hz free-rotation average and Newman projections.")


def diagram_energy_surface_map() -> None:
    slug = "energy-surface-map"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Each public name → exactly one energy surface → the gradient of that energy", size=15, weight=700))
    for x, h_ in [(40, "public name"), (330, "energy surface"), (640, "gradient used by NEB / MD")]:
        b.append(text(x, 62, h_, size=10.5, fill="var(--sf-muted)", weight=700, anchor="start"))
    rows = [
        (["uff", "mmff94"], "classical force field", "analytic", "var(--sf-accent2)"),
        (["pm3"], "MOPAC-parity PM3 (Z ≤ 17)", "stationary density: rebuild H, G(P), E_nuc at R ± h", "var(--sf-accent)"),
        (["xtb", "gfn2", "gfn0*"], "GFN2-xTB (solve_gfn2)", "warm-started central differences", "var(--sf-violet)"),
        (["gfn1"], "GFN1-xTB (solve_gfn1)", "central differences", "var(--sf-warm)"),
        (["hf3c"], "HF-3c", "central differences", "var(--sf-rose)"),
    ]
    for i, (names, surf, grad, col) in enumerate(rows):
        y = 74 + i * 58
        nx = 40
        for nm in names:
            w = 12 + 8 * len(nm)
            b.append(rect(nx, y + 8, w, 30, rx=8, fill="var(--sf-card)", stroke=col, sw=1.2))
            b.append(text(nx + w / 2, y + 28, nm, size=10.5, weight=700, mono=True))
            nx += w + 8
        b.append(arrow_between(nx + 4, y + 23, 324, y + 23, color=col, delay=0.15 * i))
        b.append(rect(330, y + 4, 280, 38, rx=10, sw=1.2))
        b.append(rect(331.5, y + 5.5, 4, 35, rx=2, fill=col, stroke="none"))
        b.append(text(346, y + 28, surf, size=10.5, weight=700, anchor="start"))
        b.append(arrow_between(614, y + 23, 636, y + 23, color=col, delay=0.15 * i + 0.1))
        b.append(rect(640, y + 4, 236, 38, rx=10, fill="var(--sf-card2)", sw=1))
        b.append(text(652, y + 28, grad, size=9 if len(grad) > 34 else 10, weight=500, anchor="start"))
    y = 74 + 5 * 58 + 6
    b.append(rect(40, y, 836, 66, rx=12, fill="var(--sf-soft)", stroke="var(--sf-warm)", sw=1.1))
    b.append(text(56, y + 22, "* \"gfn0\" is only a legacy alias of the GFN2 surface. The faithful GFN0 is xtb::solve_gfn0 (Rust only, closed shell, no gradient).", size=9.5, weight=600, anchor="start"))
    b.append(text(56, y + 42, "PM3 for Z > 17 falls back to the legacy NDDO path (recorded in notes) with full central differences.", size=9.5, weight=500, anchor="start"))
    b.append(text(56, y + 58, "Absolute energies are not comparable across backends (ethanol: −7146 kcal/mol xtb vs −14379 pm3).", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 446, "".join(b), label="Public method names mapped to their energy surfaces and gradients: uff and mmff94 analytic, pm3 stationary-density, xtb, gfn2 and the legacy gfn0 alias on GFN2 with warm-started differences, gfn1 and hf3c central differences.")


def diagram_pm3_scf_diis() -> None:
    slug = "pm3-scf-diis"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "PM3 SCF on the H₂O₂ gate geometry: Pulay DIIS (27 it.) vs damped fixed point (200 it.)", size=14, weight=700))
    px, py, pw, ph = 80, 70, 520, 260
    sx = lambda i: px + i / 200 * pw
    sy = lambda lg: py + (1 - lg) / 12 * ph  # log10 from +1 to −11
    b.append(rect(px - 56, py - 20, pw + 80, ph + 64, rx=14, sw=1.3))
    b.append(axes(px, py, pw, ph, xlabel="SCF iteration", ylabel="log₁₀ residual"))
    for lg in (0, -3, -6, -9):
        b.append(path(f"M{px},{sy(lg):.1f} L{px+pw},{sy(lg):.1f}", sw=0.6, dash="2 4"))
        b.append(text(px - 8, sy(lg) + 3, f"1e{lg}" if lg else "1", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    for it in (0, 50, 100, 150, 200):
        b.append(text(sx(it), py + ph + 13, str(it), size=8.5, fill="var(--sf-muted)", weight=500))
    # schematic shapes anchored to the documented counts and thresholds
    diis = [(i, 0.3 - i * 0.36 - (0.25 if i % 4 == 1 else 0)) for i in range(28)]
    diis[-1] = (27, -9.2)
    damp = [(i, 0.3 - 6.1 * (1 - math.exp(-i / 22)) + 0.12 * math.sin(i / 3)) for i in range(201)]
    b.append(path("M" + " L".join(f"{sx(i):.1f},{sy(v):.1f}" for i, v in damp), stroke="var(--sf-rose)", sw=2.2))
    b.append(path("M" + " L".join(f"{sx(i):.1f},{sy(v):.1f}" for i, v in diis), stroke="var(--sf-accent2)", sw=2.4))
    b.append(path(f"M{px},{sy(-9):.1f} L{px+pw},{sy(-9):.1f}", stroke="var(--sf-accent2)", sw=1.2, dash="5 4"))
    b.append(text(px + pw - 4, sy(-9) - 5, "DIIS stop: ‖FP − PF‖∞ < 1e−9", size=9, fill="var(--sf-accent2)", anchor="end"))
    b.append(circle(sx(27), sy(-9.2), 5, fill="var(--sf-accent2)", cls="sf-pulse"))
    b.append(text(sx(27) + 8, sy(-9.2) + 16, "27", size=10, fill="var(--sf-accent2)", weight=800, anchor="start"))
    b.append(circle(sx(200), sy(damp[-1][1]), 5, fill="var(--sf-rose)", cls="sf-pulse"))
    b.append(text(sx(200) - 6, sy(damp[-1][1]) - 10, "200 (budget)", size=10, fill="var(--sf-rose)", weight=800, anchor="end"))
    b.append(text(sx(110), sy(-6.2) - 26, "density still moving while ΔE stalls", size=9.5, fill="var(--sf-rose)", weight=600))
    b.append(text(px + pw / 2, py + ph + 36, "curve shapes schematic · iteration counts and thresholds from the H₂O₂ gate", size=9, fill="var(--sf-muted)", weight=450))
    x = 640
    b.append(rect(x, 50, 236, 324, rx=14, sw=1.3))
    b.append(text(x + 14, 74, "1 · Pulay DIIS", size=11.5, fill="var(--sf-accent2)", weight=700, anchor="start"))
    for i, ln in enumerate(["e = FP − PF  (S = 1)", "history: last 8 Fock", "solve (m+1)×(m+1)", "≤ 120 iterations"]):
        b.append(text(x + 22, 96 + i * 18, ln, size=9.5, weight=500, anchor="start", mono=i == 0))
    b.append(text(x + 14, 186, "2 · damped fixed point", size=11.5, fill="var(--sf-rose)", weight=700, anchor="start"))
    for i, ln in enumerate(["α ∈ [0.02, 0.9] adaptive", "|ΔE| < 1e−9 eV and", "max|ΔP| < 1e−7", "≤ 200 iterations"]):
        b.append(text(x + 22, 208 + i * 18, ln, size=9.5, weight=500, anchor="start"))
    b.append(text(x + 14, 300, "same energy to 1e−7 eV", size=10, weight=700, anchor="start"))
    b.append(text(x + 14, 320, "FD gradient scatter:", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 14, 338, "2e−3 → 6e−6 eV/Å with DIIS", size=9.5, fill="var(--sf-accent2)", weight=600, anchor="start"))
    write(slug, 900, 390, "".join(b), label="PM3 SCF on H2O2: Pulay DIIS reaches the 1e-9 commutator threshold in 27 iterations while damped fixed-point mixing uses all 200; both give the same energy, but DIIS cuts finite-difference gradient scatter from 2e-3 to 6e-6 eV per angstrom.")


def diagram_gradient_central_difference() -> None:
    slug = "gradient-central-difference"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Central difference: g_k ≈ [E(R + h e_k) − E(R − h e_k)] / 2h", size=15, weight=700))
    # left: energy curve with secant vs tangent (exaggerated h)
    E = lambda x: 0.35 * (x - 0.2) ** 2 + 0.08 * (x - 0.2) ** 3
    dE = lambda x: 0.7 * (x - 0.2) + 0.24 * (x - 0.2) ** 2
    px, py, pw, ph = 60, 70, 380, 240
    sx = lambda v: px + (v + 1.6) / 3.6 * pw
    sy = lambda v: py + ph - (v + 0.05) / 1.5 * ph
    b.append(rect(24, 50, 440, 320, rx=14, sw=1.3))
    b.append(path(f"M{px},{py} L{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    b.append(path("M" + " L".join(f"{sx(-1.6 + 3.35 * i / 120):.1f},{sy(E(-1.6 + 3.35 * i / 120)):.1f}" for i in range(121)), stroke="var(--sf-ink)", sw=2.2, cls="sf-draw"))
    x0, h = 1.1, 0.6
    for xv, col, lbl in [(x0 - h, "var(--sf-accent)", "R − h"), (x0, "var(--sf-ink)", "R"), (x0 + h, "var(--sf-accent2)", "R + h")]:
        b.append(path(f"M{sx(xv):.1f},{py+ph} L{sx(xv):.1f},{sy(E(xv)):.1f}", sw=1, dash="3 3"))
        b.append(circle(sx(xv), sy(E(xv)), 5, fill=col))
        b.append(text(sx(xv), py + ph + 14, lbl, size=9.5, fill=col, weight=700))
    m = (E(x0 + h) - E(x0 - h)) / (2 * h)
    tx = lambda xv, slope: (sx(xv), sy(E(x0) + slope * (xv - x0)))
    a1, a2 = tx(x0 - 0.55, m), tx(x0 + 0.45, m)
    b.append(path(f"M{a1[0]:.1f},{a1[1]:.1f} L{a2[0]:.1f},{a2[1]:.1f}", stroke="var(--sf-warm)", sw=2, dash="6 4"))
    t1, t2 = tx(x0 - 0.55, dE(x0)), tx(x0 + 0.45, dE(x0))
    b.append(path(f"M{t1[0]:.1f},{t1[1]:.1f} L{t2[0]:.1f},{t2[1]:.1f}", stroke="var(--sf-violet)", sw=1.6))
    b.append(text(px + 12, py + 18, "— exact tangent dE/dx", size=9.5, fill="var(--sf-violet)", weight=600, anchor="start"))
    b.append(text(px + 12, py + 36, "- - secant slope (E₊ − E₋)/2h", size=9.5, fill="var(--sf-warm)", weight=600, anchor="start"))
    b.append(text(244, 356, "h exaggerated for display · error O(h² E‴)", size=9.5, fill="var(--sf-muted)", weight=500))
    # right: error vs h (truncation vs noise)
    x = 484
    b.append(rect(x, 50, 392, 320, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "choosing h: truncation vs. SCF noise", size=11.5, weight=700, anchor="start"))
    qx, qy, qw, qh = x + 52, 94, 310, 200
    lx = lambda lh: qx + (lh + 6) / 4 * qw  # log10 h from −6 to −2
    ly = lambda le: qy + (-3 - le) / 7 * qh  # log10 err from −3 to −10
    b.append(path(f"M{qx},{qy} L{qx},{qy+qh} L{qx+qw},{qy+qh}", sw=1))
    trunc = lambda lh: math.log10(20 * (10 ** lh) ** 2)
    noise = lambda lh: math.log10(1e-12 / (2 * 10 ** lh))
    tot = lambda lh: math.log10(20 * (10 ** lh) ** 2 + 1e-12 / (2 * 10 ** lh))
    for fn, col, sw_, dash in [(trunc, "var(--sf-accent)", 1.4, "4 3"), (noise, "var(--sf-rose)", 1.4, "4 3"), (tot, "var(--sf-ink)", 2.4, None)]:
        pts = [(lx(-6 + 4 * i / 80), ly(fn(-6 + 4 * i / 80))) for i in range(81) if -10 <= fn(-6 + 4 * i / 80) <= -3]
        b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in pts), stroke=col, sw=sw_, dash=dash))
    lh0 = math.log10(5e-4)
    b.append(path(f"M{lx(lh0):.1f},{qy} L{lx(lh0):.1f},{qy+qh}", stroke="var(--sf-accent2)", sw=1.4, dash="3 3"))
    b.append(text(lx(lh0) - 4, qy + 12, "PM3/xTB h = 5·10⁻⁴ Å", size=9, fill="var(--sf-accent2)", weight=700, anchor="end"))
    b.append(circle(lx(lh0), ly(trunc(lh0)), 4.5, fill="var(--sf-accent2)"))
    b.append(text(lx(lh0) + 8, ly(trunc(lh0)) + 14, "≈ 5·10⁻⁶ eV/Å", size=8.5, fill="var(--sf-accent2)", weight=600, anchor="start"))
    for lh in (-6, -5, -4, -3, -2):
        b.append(text(lx(lh), qy + qh + 13, f"1e{lh}", size=8.5, fill="var(--sf-muted)", weight=500))
    for le in (-3, -5, -7, -9):
        b.append(text(qx - 6, ly(le) + 3, f"1e{le}", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    b.append(text(qx + qw, qy + qh + 26, "h (Å)", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    b.append(text(qx + 6, ly(-6.4), "noise ε/2h", size=9, fill="var(--sf-rose)", weight=600, anchor="start"))
    b.append(text(lx(-4.6), ly(-7.4), "truncation ∝ h²", size=9, fill="var(--sf-accent)", weight=600, anchor="start"))
    b.append(text(x + 196, 346, "ε illustrative · DIIS keeps SCF noise below truncation", size=9, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 386, "".join(b), label="Central difference geometry: the secant through E(R−h) and E(R+h) approximates the tangent with O(h²) error, and the error versus h trades truncation against SCF noise; at sci-form's h = 5e-4 angstrom truncation dominates at about 5e-6 eV per angstrom.")


def diagram_ks_gradient_units() -> None:
    slug = "ks-gradient-units"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "The 0.5292 \"discrepancy\" was a unit: 1 Bohr = 0.529177 Å", size=15, weight=700))
    # ruler
    b.append(rect(24, 50, 852, 110, rx=14, sw=1.3))
    k = 600
    b.append(rect(60, 78, k * 0.529177, 22, rx=4, fill="var(--sf-accent)", stroke="none", opacity=0.6, cls="sf-rise"))
    b.append(text(60 + k * 0.529177 / 2, 93, "1 Bohr", size=10, fill="#fff", weight=700))
    b.append(rect(60, 108, k, 22, rx=4, fill="var(--sf-accent2)", stroke="none", opacity=0.6))
    b.append(text(60 + k / 2, 123, "1 Å", size=10, fill="#fff", weight=700))
    b.append(text(60 + k + 14, 98, "same energy change over", size=9.5, weight=500, anchor="start"))
    b.append(text(60 + k + 14, 114, "a shorter step → per-Bohr", size=9.5, weight=500, anchor="start"))
    b.append(text(60 + k + 14, 130, "gradient is 0.5292× smaller", size=9.5, weight=600, anchor="start"))
    # equation row
    cards = [("analytic ks_nuclear_gradient", "0.227 Ha/Bohr", "var(--sf-accent)"), ("FD of solve_ks_dft energy per Å", "0.429 Ha/Å", "var(--sf-accent2)")]
    for i, (t_, v, col) in enumerate(cards):
        x = 24 + i * 440
        b.append(rect(x, 176, 412, 112, rx=14, sw=1.3))
        b.append(rect(x + 1.5, 177.5, 409, 4, rx=2, fill=col, stroke="none"))
        b.append(text(x + 206, 206, t_, size=10.5, weight=700))
        b.append(text(x + 206, 244, v, size=20, fill=col, weight=800, mono=True))
        b.append(text(x + 206, 272, "largest H₂O₂ component", size=9.5, fill="var(--sf-muted)", weight=450))
    b.append(text(450, 236, "=", size=20, weight=800))
    b.append(rect(24, 302, 852, 74, rx=12, fill="var(--sf-soft)", stroke="var(--sf-accent2)", sw=1.1))
    b.append(text(450, 326, "dE/dR [Ha/Bohr] = dE/dR [Ha/Å] × 0.529177    ·    0.429 × 0.5292 = 0.227", size=11.5, weight=700, mono=False))
    b.append(text(450, 348, "F80 gate after converting the reference: agreement 3·10⁻⁷ Ha/Bohr (tolerance 2.3·10⁻⁶)", size=10, fill="var(--sf-accent2)", weight=600))
    b.append(text(450, 366, "no missing XC response term — the gradient was right, the comparison mixed units", size=9.5, fill="var(--sf-muted)", weight=450))
    write(slug, 900, 390, "".join(b), label="One Bohr is 0.529177 angstrom, so a gradient per Bohr equals the gradient per angstrom times 0.5292: on H2O2 the largest component is 0.227 Hartree per Bohr, or 0.429 Hartree per angstrom, and the corrected F80 gate agrees to 3e-7.")


DIVERSE_CATEGORIES = [  # tests/fixtures/diverse_molecules.smi (131 molecules, 27 categories)
    ("aromatic", 15), ("simple", 14), ("halogen", 11), ("drug", 10), ("multifunc", 8), ("sulfur", 7), ("n-heterocycle", 7),
    ("amine", 6), ("o-heterocycle", 5), ("carbonyl", 5), ("stereo", 4), ("ether", 4), ("alcohol", 4), ("polycyclic", 3),
    ("phosphorus", 3), ("ester", 3), ("bridged", 3), ("amide", 3), ("acid", 3), ("strained", 2), ("silicon", 2),
    ("macrocycle", 2), ("fused", 2), ("boron", 2), ("urea", 1), ("selenium", 1), ("aminoacid", 1),
]
# Measured with the Python wheel, embed("C"*n, seed 0..4), median time_ms, 16-thread machine.
ALKANE_TIMING = [(17, 1.37), (32, 7.1), (47, 25.03), (62, 53.77), (92, 130.88), (122, 251.99), (152, 434.26), (182, 637.36), (242, 1188.53), (302, 1771.11)]


def diagram_benchmarks_success() -> None:
    slug = "benchmarks-success"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Diverse benchmark (Tier 1): 131 molecules, 27 categories — embed success 130 / 131", size=15, weight=700))
    cols, per = 3, 9
    for i, (cat, n) in enumerate(DIVERSE_CATEGORIES):
        c_, r_ = i // per, i % per
        x, y = 24 + c_ * 290, 56 + r_ * 34
        ok = n - (1 if cat == "strained" else 0)
        b.append(text(x + 110, y + 18, cat, size=10, weight=600, anchor="end"))
        for k in range(n):
            fill = "var(--sf-accent2)" if k < ok else "var(--sf-rose)"
            b.append(rect(x + 118 + k * 10.5, y + 7, 8.5, 14, rx=2, fill=fill, stroke="none", opacity=0.85, cls=f"sf-rise sf-d{min(k % 6 + 1, 6)}"))
        b.append(text(x + 118 + n * 10.5 + 6, y + 18, f"{ok}/{n}", size=9, fill="var(--sf-rose)" if ok < n else "var(--sf-muted)", weight=700 if ok < n else 500, anchor="start"))
    y = 56 + 9 * 34 + 8
    b.append(rect(24, y, 852, 54, rx=12, fill="var(--sf-soft)", stroke="var(--sf-rose)", sw=1.1))
    b.append(rect(40, y + 12, 8.5, 14, rx=2, fill="var(--sf-rose)", stroke="none", cls="sf-pulse"))
    b.append(text(56, y + 23, "the single failure: cubane C12C3C4C1C5C3C2C45 (strained) — 90° ring angles, no valid conformer after 160 iterations", size=10, weight=600, anchor="start"))
    b.append(text(56, y + 42, "each square = one molecule · parse 131/131 · geometry quality 130/131 · recomputed in CI by test_diverse_molecules", size=9.5, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, y + 70, "".join(b), label="Per-category embed results for the 131-molecule diverse benchmark: every molecule in 26 categories embeds, and the only failure is cubane in the strained category.")


def diagram_benchmarks_rmsd() -> None:
    slug = "benchmarks-rmsd"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "sci-form vs RDKit heavy-atom pairwise RMSD (Tier 3: one recorded local run)", size=15, weight=700))
    bins = [("< 0.1 Å", 82.8, "var(--sf-accent2)"), ("0.1–0.3", 11.6, "var(--sf-accent)"), ("0.3–0.5", 4.0, "var(--sf-warm)"), ("0.5–1.0", 1.6, "var(--sf-rose)"), ("≥ 1.0", 0.0, "var(--sf-line)")]
    px, py, pw, ph = 70, 70, 470, 240
    b.append(rect(24, 50, 540, 320, rx=14, sw=1.3))
    b.append(text(40, 70, "share of molecules per RMSD bin", size=11, weight=700, anchor="start"))
    bw = pw / len(bins)
    for i, (lbl, v, col) in enumerate(bins):
        h = v / 90 * (ph - 20)
        x = px + i * bw + 10
        b.append(rect(x, py + ph - h, bw - 20, max(h, 1.5), rx=4, fill=col, stroke="none", opacity=0.8, cls=f"sf-rise sf-d{i+1}"))
        b.append(text(x + (bw - 20) / 2, py + ph - h - 6, f"{v:.1f}%", size=10.5, weight=700, fill=col))
        b.append(text(x + (bw - 20) / 2, py + ph + 16, lbl, size=9.5, fill="var(--sf-muted)", weight=500))
    b.append(path(f"M{px},{py+ph} L{px+pw},{py+ph}", sw=1))
    b.append(text(294, 356, "bins derived from the cumulative thresholds below", size=9.5, fill="var(--sf-muted)", weight=450))
    x = 584
    b.append(rect(x, 50, 292, 320, rx=14, sw=1.3))
    b.append(text(x + 14, 72, "cumulative", size=11, weight=700, anchor="start"))
    for i, (t_, v) in enumerate([("< 0.1 Å", 82.8), ("< 0.3 Å", 94.4), ("< 0.5 Å", 98.4), ("< 1.0 Å", 100.0)]):
        y = 90 + i * 30
        b.append(text(x + 14, y + 14, t_, size=10, weight=600, anchor="start"))
        b.append(rect(x + 80, y + 3, 150, 14, rx=4, fill="var(--sf-card2)", sw=0.6))
        b.append(rect(x + 80, y + 3, 150 * v / 100, 14, rx=4, fill="var(--sf-accent2)", stroke="none", opacity=0.75))
        b.append(text(x + 278, y + 14, f"{v:.1f}%", size=10, weight=700, anchor="end"))
    b.append(text(x + 14, 230, "mean 0.064 Å · median 0.011 Å", size=10.5, weight=700, anchor="start"))
    b.append(text(x + 14, 256, "hardest categories (mean):", size=10, weight=600, anchor="start"))
    for i, (c_, v) in enumerate([("silicon", 0.543), ("selenium", 0.507), ("strained", 0.182), ("polycyclic", 0.112)]):
        b.append(text(x + 24, 276 + i * 18, f"{c_}", size=9.5, weight=500, anchor="start"))
        b.append(text(x + 278, 276 + i * 18, f"{v:.3f} Å", size=9.5, weight=600, anchor="end", mono=True))
    write(slug, 900, 386, "".join(b), label="RMSD distribution versus RDKit from one recorded run: 82.8 percent below 0.1 angstrom, 94.4 below 0.3, 98.4 below 0.5 and all below 1.0, mean 0.064 and median 0.011 angstrom.")


def diagram_benchmarks_timing() -> None:
    slug = "benchmarks-timing"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(450, 30, "Embed time vs. size — n-alkanes CₙH₂ₙ₊₂, measured (median of 5 seeds)", size=15, weight=700))
    px, py, pw, ph = 90, 70, 480, 250
    lx = lambda n: px + (math.log10(n) - 1) / (2.6 - 1) * pw
    ly = lambda t: py + ph - (math.log10(t) + 0.3) / (3.6 + 0.3) * ph
    b.append(rect(24, 50, 580, 324, rx=14, sw=1.3))
    b.append(axes(px, py, pw, ph, xlabel="atoms (log)"))
    b.append(text(px + 6, py - 6, "ms per conformer (log)", size=9, fill="var(--sf-muted)", weight=500, anchor="start"))
    for n in (10, 30, 100, 300):
        b.append(text(lx(n), py + ph + 13, str(n), size=8.5, fill="var(--sf-muted)", weight=500))
    for t in (1, 10, 100, 1000):
        b.append(text(px - 8, ly(t) + 3, str(t), size=8.5, fill="var(--sf-muted)", weight=500, anchor="end"))
        b.append(path(f"M{px},{ly(t):.1f} L{px+pw},{ly(t):.1f}", sw=0.5, dash="2 4"))
    # fitted slope on the upper half
    import statistics as st
    pts = ALKANE_TIMING[3:]
    xs_ = [math.log10(n) for n, _ in pts]; ys_ = [math.log10(t) for _, t in pts]
    mx_, my_ = st.mean(xs_), st.mean(ys_)
    slope = sum((a - mx_) * (c - my_) for a, c in zip(xs_, ys_)) / sum((a - mx_) ** 2 for a in xs_)
    for ref, lbl, col in [(2, "N²", "var(--sf-line)"), (3, "N³", "var(--sf-rose)")]:
        t0 = 53.77
        n_end = min(302, 62 * (3500 / t0) ** (1 / ref))
        b.append(path(f"M{lx(62):.1f},{ly(t0):.1f} L{lx(n_end):.1f},{ly(t0 * (n_end / 62) ** ref):.1f}", stroke=col, sw=1.2, dash="4 4"))
        b.append(text(lx(n_end) + 4, ly(t0 * (n_end / 62) ** ref) + 3, lbl, size=9, fill=col, weight=700, anchor="start"))
    b.append(polyline([(lx(n), ly(t)) for n, t in ALKANE_TIMING], color="var(--sf-accent)", sw=2.4, cls="sf-draw"))
    for n, t in ALKANE_TIMING:
        b.append(circle(lx(n), ly(t), 3.6, fill="var(--sf-accent)"))
    b.append(text(lx(302) - 4, ly(1771) + 20, "302 atoms: 1.77 s", size=9.5, fill="var(--sf-accent)", weight=700, anchor="end"))
    b.append(text(lx(17) + 6, ly(1.37) - 8, "17 atoms: 1.4 ms", size=9.5, fill="var(--sf-accent)", weight=700, anchor="start"))
    x = 624
    b.append(rect(x, 50, 252, 324, rx=14, sw=1.3))
    b.append(text(x + 14, 74, "what dominates", size=11.5, weight=700, anchor="start"))
    b.append(text(x + 14, 98, f"fitted slope (62–302 atoms): {slope:.2f}", size=10.5, fill="var(--sf-accent)", weight=700, anchor="start"))
    for i, ln in enumerate(["Floyd–Warshall smoothing: O(N³)", "BFGS: O(N²) per iteration", "eigenpairs by power iteration", "measured exponent below 3:", "N³ term not yet dominant here"]):
        b.append(text(x + 14, 126 + i * 20, ln, size=9.5, weight=500, anchor="start", fill="var(--sf-ink)" if i < 3 else "var(--sf-muted)"))
    b.append(text(x + 14, 248, "recorded throughputs", size=10.5, weight=700, anchor="start"))
    b.append(text(x + 14, 270, "diverse 131 (small): 138.9 mol/s", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 14, 290, "ChemBL 10K (≤ 100 heavy): 2.1 mol/s", size=9.5, weight=500, anchor="start"))
    b.append(text(x + 14, 320, "machine-dependent; alkane data", size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    b.append(text(x + 14, 336, "from the Python wheel, 16 threads", size=9, fill="var(--sf-muted)", weight=450, anchor="start"))
    write(slug, 900, 390, "".join(b), label=f"Measured embed time for n-alkanes from 17 to 302 atoms on a log-log plot, rising from 1.4 milliseconds to 1.77 seconds with a fitted exponent of {slope:.2f}, between the N squared and N cubed guides.")


def diagram_ani_architecture() -> None:
    slug = "ani-architecture"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "ANI inference in sci-form: per-atom AEV (4032 values) → element network → E = Σ Eᵢ", size=15, weight=700))
    W, H, Y = 196, 270, 54
    cards = [("1", "Cell-list neighbors", "var(--sf-accent)"), ("2", "Radial AEV", "var(--sf-accent2)"), ("3", "Angular AEV", "var(--sf-violet)"), ("4", "Element network", "var(--sf-warm)"), ("5", "Energy + forces", "var(--sf-rose)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 cutoff sphere
    x = xs[0]; cx, cy = x + 98, Y + 124
    b.append(circle(cx, cy, 62, fill="var(--sf-accent)", stroke="var(--sf-accent)", sw=1.2, opacity=0.1))
    b.append(circle(cx, cy, 62, fill="none", stroke="var(--sf-accent)", sw=1.2))
    for ang, rr, el, inside in [(20, 34, "H", True), (140, 40, "C", True), (250, 50, "O", True), (320, 82, "C", False), (80, 76, "H", False)]:
        px_, py_ = cx + rr * math.cos(math.radians(ang)), cy - rr * math.sin(math.radians(ang))
        if inside:
            b.append(path(f"M{cx},{cy} L{px_:.1f},{py_:.1f}", stroke="var(--sf-accent)", sw=1, dash="3 3"))
        b.append(atom_ball(slug, px_, py_, el, 7 if el != "H" else 5, label=False))
    b.append(atom_ball(slug, cx, cy, "C", 10))
    b.append(text(cx - 40, cy + 78, "R_c = 5.2 Å", size=9.5, fill="var(--sf-accent)", weight=700))
    b.append(text(x + W / 2, Y + 226, "src/ani/neighbor.rs", size=9.5, weight=600, mono=True))
    b.append(text(x + W / 2, Y + 246, "atoms outside R_c ignored", **CAPTION))
    # 2 radial
    x = xs[1]
    for k in range(8):
        rs = 0.9 + k * 0.55
        pts = [(x + 20 + r * 30, Y + 180 - 110 * math.exp(-16 * (r - rs) ** 2) * 0.5 * (math.cos(math.pi * r / 5.2) + 1)) for r in [i * 0.05 for i in range(105)]]
        b.append(path("M" + " L".join(f"{a:.1f},{c:.1f}" for a, c in pts), stroke="var(--sf-accent2)", sw=1.4, opacity=0.8))
    b.append(path(f"M{x+20},{Y+180} L{x+178},{Y+180}", sw=1))
    b.append(text(x + 178, Y + 194, "r (Å)", size=8.5, fill="var(--sf-muted)", anchor="end", weight=500))
    b.append(text(x + W / 2, Y + 226, "7 species × (8 η × 8 Rₛ) = 448", size=10, weight=600))
    b.append(text(x + W / 2, Y + 246, "Gaussians × cosine cutoff f_c(r)", **CAPTION))
    # 3 angular
    x = xs[2]; cx, cy = x + 98, Y + 140
    b.append(path(f"M{cx-50},{cy-50} L{cx},{cy} L{cx+56},{cy-40}", stroke="var(--sf-ink)", sw=2))
    b.append(path(f"M{cx-20},{cy-20} A28,28 0 0 1 {cx+23},{cy-16}", stroke="var(--sf-violet)", sw=1.8))
    b.append(text(cx + 2, cy - 34, "θ", size=12, fill="var(--sf-violet)", weight=700))
    b.append(atom_ball(slug, cx, cy, "C", 10)); b.append(atom_ball(slug, cx - 50, cy - 50, "O", 8, label=False)); b.append(atom_ball(slug, cx + 56, cy - 40, "H", 6))
    b.append(text(x + W / 2, Y + 206, "28 species pairs (lo ≤ hi)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 226, "× 4 η × 4 Rₛ × 1 ζ × 8 θₛ = 3584", size=10, weight=600))
    b.append(text(x + W / 2, Y + 246, "448 + 3584 = 4032 per atom", size=9.5, fill="var(--sf-violet)", weight=700))
    # 4 network per element
    x = xs[3]
    layers = [5, 4, 3, 1]
    for li, n in enumerate(layers):
        for j in range(n):
            nx = x + 40 + li * 40; ny = Y + 130 + (j - (n - 1) / 2) * 24
            if li < len(layers) - 1:
                for j2 in range(layers[li + 1]):
                    ny2 = Y + 130 + (j2 - (layers[li + 1] - 1) / 2) * 24
                    b.append(path(f"M{nx},{ny:.1f} L{nx+40},{ny2:.1f}", sw=0.5, opacity=0.6))
    for li, n in enumerate(layers):
        for j in range(n):
            b.append(circle(x + 40 + li * 40, Y + 130 + (j - (n - 1) / 2) * 24, 6, fill="var(--sf-warm)" if li < 3 else "var(--sf-rose)"))
    b.append(text(x + 40, Y + 76, "AEV", size=8.5, fill="var(--sf-muted)"))
    b.append(text(x + 160, Y + 116, "Eᵢ", size=10, fill="var(--sf-rose)", weight=700))
    b.append(text(x + W / 2, Y + 206, "one net per element (H, C, N, O…)", size=10, weight=600))
    b.append(text(x + W / 2, Y + 226, "CELU / GELU activations", **CAPTION))
    b.append(text(x + W / 2, Y + 246, "in-tree weights: random test model", size=9.5, fill="var(--sf-rose)", weight=600))
    # 5 energy + forces
    x = xs[4]
    b.append(text(x + W / 2, Y + 80, "E = Σᵢ Eᵢ(AEVᵢ)", size=12, weight=700))
    b.append(text(x + W / 2, Y + 112, "Fᵢ = −∂E/∂Rᵢ", size=12, weight=700))
    b.append(text(x + W / 2, Y + 138, "backward through the net", size=9.5, weight=500))
    b.append(text(x + W / 2, Y + 156, "and through the AEV (gradients.rs)", size=9.5, weight=500))
    b.append(rect(x + 14, Y + 180, W - 28, 74, rx=9, fill="var(--sf-soft)", stroke="var(--sf-rose)", sw=1))
    b.append(text(x + W / 2, Y + 200, "ANI-2x .netFormat weights", size=9.5, fill="var(--sf-rose)", weight=700))
    b.append(text(x + W / 2, Y + 218, "cannot be loaded; outputs of", size=9, weight=500))
    b.append(text(x + W / 2, Y + 234, "make_test_model are not", size=9, weight=500))
    b.append(text(x + W / 2, Y + 250, "physically meaningful", size=9, weight=500))
    write(slug, 1100, 344, "".join(b), label="ANI inference: a 5.2 angstrom cell-list cutoff, 448 radial and 3584 angular AEV entries per atom, per-element networks summed to the energy, analytic forces by back-propagation, and the caveat that only a random test model ships in-tree.")


def diagram_hf3c_architecture() -> None:
    slug = "hf3c-architecture"
    begin(slug)
    b = [element_defs(slug)]
    b.append(text(550, 30, "HF-3c for water — real compute_hf3c output (Hartree)", size=15, weight=700))
    W, H, Y = 254, 290, 54
    cards = [("1", "Integrals (Obara–Saika)", "var(--sf-accent)"), ("2", "Roothaan–Hall SCF", "var(--sf-accent2)"), ("3", "Three corrections", "var(--sf-warm)"), ("4", "HF-3c energy", "var(--sf-violet)")]
    xs = stage_row(b, 24, Y, W, H, 18, cards)
    # 1 integrals
    x = xs[0]
    rows = [("S", "overlap"), ("T", "kinetic"), ("V", "nuclear attraction"), ("(μν|λσ)", "ERI, 8-fold symmetry")]
    for k, (sym, d_) in enumerate(rows):
        yy = Y + 56 + k * 34
        b.append(rect(x + 14, yy, 226, 28, rx=7, fill="var(--sf-card2)", sw=0.8))
        b.append(text(x + 26, yy + 18, sym, size=10.5, fill="var(--sf-accent)", weight=700, anchor="start", mono=True))
        b.append(text(x + 228, yy + 18, d_, size=9.5, weight=500, anchor="end"))
    b.append(text(x + W / 2, Y + 214, "7 AOs: O 1s 2s 2p×3 · H 1s ×2", size=10, weight=600))
    b.append(text(x + W / 2, Y + 234, "L ≤ 1 (s, p shells)", **CAPTION))
    b.append(text(x + W / 2, Y + 252, "unique ERIs: 406 of 7⁴ = 2401", **CAPTION))
    # 2 SCF
    x = xs[1]
    b.append(text(x + W / 2, Y + 66, "F′C′ = C′E  (Löwdin)", size=10.5, weight=700))
    b.append(text(x + W / 2, Y + 88, "DIIS extrapolation", size=10, weight=500))
    b.append(text(x + W / 2, Y + 110, "E_elec = ½ Tr[P(H + F)]", size=10.5, weight=700))
    occ = [-20.24, -1.25, -0.61, -0.44, -0.39]
    vir = [0.57, 0.72]
    slot = {-20.24: Y + 204, -1.25: Y + 182, -0.61: Y + 170, -0.44: Y + 158, -0.39: Y + 146, 0.57: Y + 130, 0.72: Y + 120}
    sl = lambda e: slot[e]
    for e in occ + vir:
        yy = sl(e)
        b.append(path(f"M{x+70},{yy:.1f} L{x+150},{yy:.1f}", stroke="var(--sf-accent2)" if e < 0 else "var(--sf-rose)", sw=2))
    b.append(text(x + 156, sl(-0.39) + 3, "HOMO −0.387", size=8.5, fill="var(--sf-accent2)", anchor="start"))
    b.append(text(x + 156, sl(0.57) + 6, "LUMO +0.570", size=8.5, fill="var(--sf-rose)", anchor="start"))
    b.append(text(x + 156, sl(-20.24) + 3, "O 1s −20.24", size=8.5, fill="var(--sf-muted)", anchor="start"))
    b.append(text(x + W / 2, Y + 234, "converged in 9 iterations", size=10, weight=600))
    b.append(text(x + W / 2, Y + 252, "10 electrons · gap 0.956 Ha", **CAPTION))
    # 3 corrections bars
    x = xs[2]
    corr = [("D3(BJ)", -0.00263, "var(--sf-warm)"), ("gCP", 0.00183, "var(--sf-violet)"), ("SRB", 0.17233, "var(--sf-rose)")]
    zero = x + 70
    for k, (nm, v, col) in enumerate(corr):
        yy = Y + 64 + k * 46
        w = v / 0.18 * 150
        b.append(text(x + 14, yy + 16, nm, size=10.5, fill=col, weight=700, anchor="start"))
        b.append(rect(zero + min(w, 0), yy + 4, max(abs(w), 2), 18, rx=3, fill=col, stroke="none", opacity=0.75))
        b.append(text(x + 240, yy + 36, f"{v:+.5f}", size=9.5, weight=600, anchor="end", mono=True))
    b.append(path(f"M{zero},{Y+60} L{zero},{Y+200}", sw=1, dash="2 3"))
    b.append(text(x + W / 2, Y + 214, "dispersion · basis-set superposition", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 232, "· short-range basis", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 252, "D3(BJ) matches dftd4 to 6.8·10⁻⁹ Ha", **CAPTION))
    # 4 total
    x = xs[3]
    terms = [("E_HF", "−74.96293"), ("+ D3", "−0.00263"), ("+ gCP", "+0.00183"), ("+ SRB", "+0.17233")]
    for k, (t_, v) in enumerate(terms):
        yy = Y + 68 + k * 24
        b.append(text(x + 40, yy, t_, size=10.5, weight=600, anchor="start"))
        b.append(text(x + 214, yy, v, size=10.5, weight=600, anchor="end", mono=True))
    b.append(path(f"M{x+36},{Y+150} L{x+218},{Y+150}", sw=1.2))
    b.append(text(x + 40, Y + 172, "E_HF-3c", size=11.5, fill="var(--sf-violet)", weight=800, anchor="start"))
    b.append(text(x + 214, Y + 172, "−74.79141", size=11.5, fill="var(--sf-violet)", weight=800, anchor="end", mono=True))
    b.append(text(x + W / 2, Y + 214, "HF/STO-3G layer: PySCF gate 1.2·10⁻⁸ Ha", size=9.5, weight=600))
    b.append(text(x + W / 2, Y + 234, "no external oracle yet for the full", **CAPTION))
    b.append(text(x + W / 2, Y + 252, "HF-3c energy (unit tests only)", **CAPTION))
    write(slug, 1100, 364, "".join(b), label="HF-3c for water from compute_hf3c: seven-function integrals, a nine-iteration DIIS SCF giving E_HF of −74.96293 Hartree, D3 −0.00263, gCP +0.00183 and SRB +0.17233 corrections, and a total of −74.79141 Hartree.")



if __name__ == "__main__":
    main()
