#!/usr/bin/env python3
"""Check the documentation diagrams for the invariants the theme relies on.

- every ``<SvgDiagram src="/svg/...">`` reference resolves;
- each SVG has the shared ``sf-d`` theme class and the dark-mode override;
- no SVG paints a hard-coded white canvas;
- each SVG contains at least one visible animation;
- no rectangle, circle or text extends outside its viewBox.

Exit code is non-zero on the first failed invariant.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DOCS = ROOT / "docs"
SVG_DIR = DOCS / "public" / "svg"


def fail(message: str) -> None:
    print(f"FAIL: {message}")
    sys.exit(1)


def main() -> None:
    referenced: set[str] = set()
    for page in DOCS.rglob("*.md"):
        if ".vitepress" in page.parts:
            continue
        for match in re.finditer(r'<SvgDiagram\s+src="/svg/([^"]+)"', page.read_text(encoding="utf-8")):
            referenced.add(match.group(1))

    missing = sorted(name for name in referenced if not (SVG_DIR / name).exists())
    if missing:
        fail(f"missing SVG files: {missing}")

    files = sorted(SVG_DIR.glob("*.svg"))
    if not files:
        fail("no SVG files found")

    for path in files:
        text = path.read_text(encoding="utf-8")
        name = path.name
        if '<svg class="sf-d"' not in text:
            fail(f"{name}: missing the sf-d theme class")
        if "html.dark .sf-d" not in text:
            fail(f"{name}: missing the dark-mode palette")
        if re.search(r'<rect[^>]*fill="white"', text) or re.search(r'<rect[^>]*fill="#fff', text):
            fail(f"{name}: hard-coded white rectangle")
        body = text.split("</defs>", 1)[-1]
        if not re.search(r"sf-(flow|pulse|float|spin|draw|halo|blink|rise)|animateMotion", body):
            fail(f"{name}: no visible animation")

        viewbox = re.search(r'viewBox="0 0 (\d+) (\d+)"', text)
        if not viewbox:
            fail(f"{name}: missing viewBox")
        width, height = map(int, viewbox.groups())
        max_x = max_y = 0.0
        for match in re.finditer(
            r'<rect[^>]*x="([\d.]+)"[^>]*y="([\d.]+)"[^>]*width="([\d.]+)"[^>]*height="([\d.]+)"',
            text,
        ):
            x, y, w, h = map(float, match.groups())
            max_x, max_y = max(max_x, x + w), max(max_y, y + h)
        for match in re.finditer(
            r'<circle[^>]*cx="([\d.]+)"[^>]*cy="([\d.]+)"[^>]*r="([\d.]+)"', text
        ):
            x, y, r = map(float, match.groups())
            max_x, max_y = max(max_x, x + r), max(max_y, y + r)
        for match in re.finditer(r'<text[^>]*x="([\d.]+)"[^>]*y="([\d.]+)"', text):
            x, y = map(float, match.groups())
            max_x, max_y = max(max_x, x), max(max_y, y)
        if max_x > width + 1 or max_y > height + 1:
            fail(f"{name}: content ({max_x:.1f}, {max_y:.1f}) exceeds viewBox ({width}, {height})")

    print(f"OK: {len(referenced)} references, {len(files)} animated theme-aware SVGs, no white canvas, no overflow")


if __name__ == "__main__":
    main()
