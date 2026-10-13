#!/usr/bin/env python3
"""Generate `tests/fixtures/validation/hla_couplings.json` from the reference
Haasnoot–de Leeuw–Altona (HLA) web calculator.

Reference implementation: R. Stenutz, "Generalized 3JHH calculation acc.
Haasnoot et al." (https://www.stenutz.eu/conf/haasnoot.php), which implements
C.A.G. Haasnoot, F.A.A.M. de Leeuw, C. Altona, *Tetrahedron* 36 (1980)
2783–2792.

The calculator takes the four substituent slots of the H-C-C-H fragment
(S1/S2/S3 on the near carbon, S4 on the far one, per its trans/gauche schemes)
and the H,H torsion angle, and returns 3J(H,H) in Hz. This script probes a
grid of substituents × angles and stores the results.

Usage:
    python scripts/validation/oracle_hla.py [--out FILE]
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "tests/fixtures/validation/hla_couplings.json"
URL = "https://www.stenutz.eu/conf/haasnoot.php"

ANGLES = [0, 30, 60, 90, 120, 150, 180]
# Single substituent in each slot (the others H) plus the all-H baseline.
PATTERNS = [
    ("all_H", ("H", "H", "H", "H")),
    ("s1_C", ("C", "H", "H", "H")),
    ("s2_C", ("H", "C", "H", "H")),
    ("s3_C", ("H", "H", "C", "H")),
    ("s4_C", ("H", "H", "H", "C")),
    ("s1_O", ("O", "H", "H", "H")),
    ("s1_N", ("N", "H", "H", "H")),
    ("s1_F", ("F", "H", "H", "H")),
    ("s1_Cl", ("Cl", "H", "H", "H")),
    ("s1_Br", ("Br", "H", "H", "H")),
    ("s1_S", ("S", "H", "H", "H")),
    ("s1_Si", ("Si", "H", "H", "H")),
]

COUPLING_RE = re.compile(r"coupling for\s+([-\d.]+)\s*(?:&#176;|°)?\s*is\s*([-\d.]+)\s*Hz")


def query(pattern: tuple[str, str, str, str], angle: float) -> float:
    data = urllib.parse.urlencode(
        {
            "action": "Calculate J",
            "s1": pattern[0],
            "s2": pattern[1],
            "s3": pattern[2],
            "s4": pattern[3],
            "angle": f"{angle:g}",
        }
    ).encode()
    req = urllib.request.Request(URL, data=data, method="POST")
    with urllib.request.urlopen(req, timeout=30) as resp:
        html = resp.read().decode("utf-8", "replace")
    m = COUPLING_RE.search(html)
    if not m:
        raise RuntimeError(f"no coupling in response for {pattern} at {angle}")
    return float(m.group(2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    cases = []
    for name, pattern in PATTERNS:
        for angle in ANGLES:
            j = query(pattern, angle)
            cases.append(
                {
                    "name": name,
                    "substituents": list(pattern),
                    "angle_deg": angle,
                    "j_hz": j,
                }
            )
            print(f"  {name:8s} {angle:5.0f}° → {j:7.3f} Hz")
            time.sleep(0.4)

    fixture = {
        "source": "R. Stenutz HLA calculator (stenutz.eu/conf/haasnoot.php), "
        "implementing Haasnoot, de Leeuw, Altona, Tetrahedron 36 (1980) 2783",
        "units": {"angle": "degrees", "j": "Hz"},
        "note": "Substituent slots S1/S2/S3 (near carbon) and S4 (far carbon) "
        "per the calculator's trans/gauche schemes; the others are H.",
        "cases": cases,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(fixture, indent=1) + "\n")
    print(f"wrote {args.out.relative_to(REPO)} ({len(cases)} cases)")


if __name__ == "__main__":
    main()
