#!/usr/bin/env python3
"""pymatgen oracle for materials cross-validation (fichas F64, F65, F66, F69).

- F69: simulated powder XRD patterns (Cu-Kα) for reference structures.
- F65/F66: CIF texts and PBC distances / minimum-image reference values.

Generates:
  - tests/fixtures/validation/pymatgen_xrd.json
  - tests/fixtures/validation/pymatgen_cif.json
  - tests/fixtures/validation/pymatgen_pbc.json

Usage:
    python scripts/validation/oracle_pymatgen.py [xrd|cif|pbc|all]
"""

from __future__ import annotations

import sys
import time

import numpy as np

from common import ANG2BOHR, Report, emit_fixture_report

try:
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.core import Lattice, Structure
    from pymatgen.io.cif import CifParser, CifWriter
except ImportError as exc:  # pragma: no cover
    print(f"pymatgen is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

from importlib.metadata import version as _pkg_version

REFERENCE = f"pymatgen {_pkg_version('pymatgen')}"

# ---------------------------------------------------------------------------
# Structures (lattice parameters in Angstrom / degrees)
# ---------------------------------------------------------------------------

STRUCTURES = [
    {
        "name": "NaCl",
        "lattice": dict(a=5.6402, b=5.6402, c=5.6402, alpha=90.0, beta=90.0, gamma=90.0),
        "sites": [("Na", [0.0, 0.0, 0.0]), ("Cl", [0.5, 0.5, 0.5])],
    },
    {
        "name": "Si",
        "lattice": dict(a=5.431, b=5.431, c=5.431, alpha=90.0, beta=90.0, gamma=90.0),
        "sites": [
            ("Si", [0.0, 0.0, 0.0]),
            ("Si", [0.25, 0.25, 0.25]),
            ("Si", [0.5, 0.5, 0.0]),
            ("Si", [0.25, 0.75, 0.75]),
            ("Si", [0.75, 0.25, 0.75]),
            ("Si", [0.75, 0.75, 0.25]),
            ("Si", [0.0, 0.5, 0.5]),
            ("Si", [0.5, 0.0, 0.5]),
        ],
    },
    {
        "name": "TiO2 rutile",
        "lattice": dict(a=4.5937, b=4.5937, c=2.9587, alpha=90.0, beta=90.0, gamma=90.0),
        "sites": [
            ("Ti", [0.0, 0.0, 0.0]),
            ("Ti", [0.5, 0.5, 0.5]),
            ("O", [0.3047, 0.3047, 0.0]),
            ("O", [0.6953, 0.6953, 0.0]),
            ("O", [0.1953, 0.8047, 0.5]),
            ("O", [0.8047, 0.1953, 0.5]),
        ],
    },
    {
        "name": "Cu",
        "lattice": dict(a=3.6149, b=3.6149, c=3.6149, alpha=90.0, beta=90.0, gamma=90.0),
        "sites": [("Cu", [0.0, 0.0, 0.0])],
    },
]

Z2SYM = {6: "C", 8: "O", 11: "Na", 14: "Si", 22: "Ti", 29: "Cu", 17: "Cl", 1: "H", 7: "N"}


def build_structure(entry) -> Structure:
    lat = Lattice.from_parameters(
        entry["lattice"]["a"],
        entry["lattice"]["b"],
        entry["lattice"]["c"],
        entry["lattice"]["alpha"],
        entry["lattice"]["beta"],
        entry["lattice"]["gamma"],
    )
    return Structure(lat, [s[0] for s in entry["sites"]], [s[1] for s in entry["sites"]])


def run_xrd() -> Report:
    calc = XRDCalculator(wavelength="CuKa")
    cases = []
    for entry in STRUCTURES:
        structure = build_structure(entry)
        pattern = calc.get_pattern(structure, two_theta_range=(0.0, 90.0))
        cases.append(
            {
                "name": entry["name"],
                "lattice": entry["lattice"],
                "elements": [entry["sites"][i][0] for i in range(len(entry["sites"]))],
                "frac_coords": [entry["sites"][i][1] for i in range(len(entry["sites"]))],
                "two_theta": [float(x) for x in pattern.x],
                "intensity": [float(y) for y in pattern.y],
            }
        )
    payload = {
        "reference": REFERENCE,
        "wavelength": 1.54184,  # pymatgen XRDCalculator CuKa weighted average
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F69",
        observable="powder_xrd_vs_pymatgen",
        reference_tool=REFERENCE,
        units="2θ degrees / relative intensity",
        tolerance=0.01,
        level="L2",
        n_cases=len(cases),
        metric="max_abs_two_theta",
        value=float("nan"),
        notes=(
            "Peak positions tolerance 0.01°; intensity correlation ≥ 0.95 in "
            "the plan. Measured by tests/cross_validation.rs::f69_xrd."
        ),
    )
    emit_fixture_report("pymatgen_xrd.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F65 — CIF corpus (structures written by pymatgen, parsed by both sides)
# ---------------------------------------------------------------------------


def run_cif() -> Report:
    cases = []
    for entry in STRUCTURES:
        structure = build_structure(entry)
        cif_text = str(CifWriter(structure))
        parser = CifParser.from_str(cif_text)
        parsed = parser.parse_structures(primitive=False)[0]
        cases.append(
            {
                "name": entry["name"],
                "cif": cif_text,
                "lattice_params": [
                    parsed.lattice.a,
                    parsed.lattice.b,
                    parsed.lattice.c,
                    parsed.lattice.alpha,
                    parsed.lattice.beta,
                    parsed.lattice.gamma,
                ],
                "n_sites": len(parsed),
                "elements": sorted({site.species_string for site in parsed}),
            }
        )
    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F65",
        observable="cif_roundtrip_vs_pymatgen",
        reference_tool=REFERENCE,
        units="Å / degrees",
        tolerance=1e-6,
        level="L2",
        n_cases=len(cases),
        metric="max_abs_lattice_param",
        value=float("nan"),
        notes=(
            "Cell parameters ≤ 1e-6 Å/deg; round-trip stable. Measured by "
            "tests/cross_validation.rs::f65_cif."
        ),
    )
    emit_fixture_report("pymatgen_cif.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F66 — PBC minimum-image distances vs ASE
# ---------------------------------------------------------------------------


def run_pbc() -> Report:
    from ase import Atoms
    from ase.geometry import find_mic

    cases = []
    for entry in STRUCTURES:
        structure = build_structure(entry)
        lat = np.array(structure.lattice.matrix)
        frac = np.array([site.frac_coords for site in structure])
        cart = frac @ lat
        symbols = [site.species_string for site in structure]
        atoms = Atoms(symbols=symbols, positions=cart, cell=lat, pbc=True)

        pairs = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                d, _ = find_mic(cart[j] - cart[i], lat, pbc=True)
                pairs.append({"i": i, "j": j, "mic_distance": float(np.linalg.norm(d))})

        cases.append(
            {
                "name": entry["name"],
                "lattice": entry["lattice"],
                "elements": symbols,
                "frac_coords": frac.tolist(),
                "pairs": pairs,
            }
        )
    payload = {
        "reference": f"ASE {__import__('ase').__version__} find_mic",
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F66",
        observable="minimum_image_distance_vs_ase",
        reference_tool=f"ASE {__import__('ase').__version__}",
        units="Å",
        tolerance=1e-10,
        level="L0/L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Minimum-image distances ≤ 1e-10 Å; supercell/wrap invariants "
            "checked in the Rust test. Measured by "
            "tests/cross_validation.rs::f66_pbc."
        ),
    )
    emit_fixture_report("pymatgen_pbc.json", payload, report)
    return report


def main() -> int:
    parser = __import__("argparse").ArgumentParser(description=__doc__)
    parser.add_argument(
        "mode", nargs="?", default="all", choices=["xrd", "cif", "pbc", "all"]
    )
    args = parser.parse_args()
    if args.mode in ("xrd", "all"):
        run_xrd()
    if args.mode in ("cif", "all"):
        run_cif()
    if args.mode in ("pbc", "all"):
        run_pbc()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
