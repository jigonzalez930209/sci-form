#!/usr/bin/env python3
"""dftd4 oracle for sci-form D4 cross-validation (ficha F46).

Reuses the deterministic 3D corpus from `oracle_tblite.py` and evaluates the
standalone DFT-D4 dispersion energies with the reference dftd4 library
(PBE functional parameterization).

Generates `tests/fixtures/validation/dftd4_dispersion.json` with per molecule:
  - elements, geometry (Angstrom)
  - pbe: total dispersion energy (Hartree)

Usage:
    python scripts/validation/oracle_dftd4.py
"""

from __future__ import annotations

import sys
import time

import numpy as np

from common import ANG2BOHR, Report, emit_fixture_report

try:
    import dftd4
    from dftd4.interface import DampingParam, DispersionModel
except ImportError as exc:  # pragma: no cover
    print(f"dftd4 is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

# Reuse the fixed corpus and geometry builder from the tblite oracle.
from oracle_tblite import SMILES, embed  # noqa: E402

DFTD4_VERSION = getattr(dftd4, "__version__", "unknown")
REFERENCE = f"dftd4 {DFTD4_VERSION}"

METHODS = ["pbe"]

# Heavy-element cases (Z > 36) with explicit geometries: the SMILES corpus from
# the tblite oracle only covers light elements, so the D4 coverage extension
# (validation F46) is gated with these fixed structures.
HEAVY_CASES = [
    ("PbI2", [82, 53, 53], [[0.0, 0.0, 0.0], [2.8, 0.9, 0.0], [-2.8, 0.9, 0.0]]),
    ("Au2", [79, 79], [[0.0, 0.0, 0.0], [0.0, 0.0, 2.47]]),
    ("HgCl2", [80, 17, 17], [[0.0, 0.0, 0.0], [0.0, 0.0, 2.25], [0.0, 0.0, -2.25]]),
    ("LaH3", [57, 1, 1, 1], [[0.0, 0.0, 0.0], [1.9, 0.0, 0.0], [-0.95, 1.65, 0.0], [-0.95, -1.65, 0.0]]),
    ("Rn2", [86, 86], [[0.0, 0.0, 0.0], [0.0, 0.0, 4.0]]),
    ("UO2", [92, 8, 8], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.8], [0.0, 0.0, -1.8]]),
    ("LrH", [103, 1], [[0.0, 0.0, 0.0], [0.0, 0.0, 1.9]]),
]


def dispersion(method: str, elements, coords) -> float:
    numbers = np.array(elements, dtype=np.int32)
    pos_bohr = np.array(coords).reshape(-1, 3) * ANG2BOHR
    model = DispersionModel(numbers, pos_bohr)
    param = DampingParam(method=method)
    res = model.get_dispersion(param, grad=False)
    return float(res.get("energy"))


def dispersion_d3(elements, coords) -> float:
    """s-dftd3 reference with the HF-3c D3-BJ parameters."""
    from dftd3.interface import DispersionModel, RationalDampingParam

    numbers = np.array(elements, dtype=np.int32)
    # s-dftd3 takes Bohr (like dftd4); our corpus geometries are Angstrom.
    pos_bohr = np.array(coords).reshape(-1, 3) * ANG2BOHR
    model = DispersionModel(numbers, pos_bohr)
    param = RationalDampingParam(s6=1.0, s8=0.8777, a1=0.4171, a2=2.9149)
    res = model.get_dispersion(param, grad=False)
    return float(res.get("energy"))


def run() -> Report:
    cases = []
    corpus = list(SMILES) + [(name, None) for name, _, _ in HEAVY_CASES]
    for name, smiles in corpus:
        if smiles is None:
            _, elements, coords3 = next(c for c in HEAVY_CASES if c[0] == name)
            coords = [v for xyz in coords3 for v in xyz]
        else:
            elements, coords = embed(smiles)
        entry = {"name": name, "elements": elements, "coords": coords}
        for method in METHODS:
            entry[method] = {"e_tot_hartree": dispersion(method, elements, coords)}
        try:
            entry["d3_hf3c"] = {"e_tot_hartree": dispersion_d3(elements, coords)}
        except ImportError:
            pass
        cases.append(entry)

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "methods": METHODS,
        "cases": cases,
    }
    report = Report(
        ficha="F46",
        observable="d4d3_dispersion_fixture_provenance",
        reference_tool=REFERENCE,
        units="hartree",
        tolerance=1e-3 / 627.5094740631,  # 1e-3 kcal/mol in Ha
        level="L4",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Provenance only. Both standalone paths are documented gaps "
            "(element-level C6 without coordination/charge dependence), "
            "tracked by tests/cross_validation.rs::f46_d3_d4_known_gaps_tracked."
        ),
    )
    emit_fixture_report("dftd4_dispersion.json", payload, report)
    return report


def main() -> int:
    run()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
