#!/usr/bin/env python3
"""Fukui/EAS ranking oracle for sci-form reactivity validation (ficha F47).

Builds a series of monosubstituted benzenes + pyridine deterministically
(RDKit ETKDG + MMFF94) and stores, per molecule, the ring-carbon classes
(ortho/meta/para or alpha/beta/gamma) together with the experimentally known
electrophilic-substitution site set.

Ground truth: textbook EAS orientation rules (activating groups ortho/para,
deactivating meta-directors, pyridine beta), confirmed with DFT HOMO Mulliken
populations (RHF/6-31G*) for the borderline cases (nitrobenzene,
benzaldehyde, pyridine) during development.

Generates `tests/fixtures/validation/fukui_eas.json`.

Usage:
    python scripts/validation/oracle_fukui.py
"""

from __future__ import annotations

import json
import time

from common import Report, emit_fixture_report

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError as exc:  # pragma: no cover
    raise SystemExit(f"RDKit is required for this oracle: {exc}")

REFERENCE = "EAS textbook orientation rules (DFT-confirmed)"

# molecule -> (SMILES, expected EAS site class among ring CH)
SERIES = {
    "toluene": ("Cc1ccccc1", "ortho_para"),
    "aniline": ("Nc1ccccc1", "ortho_para"),
    "phenol": ("Oc1ccccc1", "ortho_para"),
    "anisole": ("COc1ccccc1", "ortho_para"),
    "chlorobenzene": ("Clc1ccccc1", "ortho_para"),
    "nitrobenzene": ("O=[N+]([O-])c1ccccc1", "meta"),
    "benzaldehyde": ("O=Cc1ccccc1", "meta"),
    "pyridine": ("c1ccncc1", "beta"),
}


def classify(mol: Chem.Mol, name: str):
    """Return (classes dict atom_idx -> class, candidate ring-C indices)."""
    ring = [a.GetIdx() for a in mol.GetAtoms() if a.IsInRingSize(6)]
    assert len(ring) == 6, f"{name}: no 6-ring found"
    dist = Chem.GetDistanceMatrix(mol)
    if name == "pyridine":
        ref = next(i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
        names = {1: "alpha", 2: "beta", 3: "gamma"}
        classes = {}
        for i in ring:
            if i == ref:
                classes[i] = "n"
            else:
                classes[i] = names[int(dist[ref, i])]
        candidates = [i for i in ring if i != ref]
    else:
        # ipso = ring atom bound to a non-ring heavy atom
        ipso = next(
            i
            for i in ring
            if any(
                n.GetAtomicNum() > 1 and n.GetIdx() not in ring
                for n in mol.GetAtomWithIdx(i).GetNeighbors()
            )
        )
        names = {0: "ipso", 1: "ortho", 2: "meta", 3: "para"}
        classes = {i: names[int(dist[ipso, i])] for i in ring}
        candidates = [i for i in ring if i != ipso]
    return classes, candidates


def run() -> Report:
    cases = []
    for name, (smi, expected) in SERIES.items():
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF47
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"embedding failed for {name}")
        AllChem.MMFFOptimizeMolecule(mol)
        conf = mol.GetConformer()
        elements, coords = [], []
        for atom in mol.GetAtoms():
            elements.append(atom.GetAtomicNum())
            p = conf.GetAtomPosition(atom.GetIdx())
            coords += [p.x, p.y, p.z]
        classes, candidates = classify(mol, name)
        cases.append(
            {
                "name": name,
                "smiles": smi,
                "elements": elements,
                "coords": coords,
                "ring_classes": {str(k): v for k, v in classes.items()},
                "candidates": candidates,
                "expected": expected,
            }
        )

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F47",
        observable="eas_site_ranking",
        reference_tool=REFERENCE,
        units="rank",
        tolerance=0.8,
        level="L4",
        n_cases=len(cases),
        metric="fraction_correct",
        value=float("nan"),
        notes=(
            "Top-1 f- among ring CH must fall in the expected EAS set in "
            ">= 80% of the series (>= 7/8). Measured by "
            "tests/cross_validation.rs::f47_fukui."
        ),
    )
    emit_fixture_report("fukui_eas.json", payload, report)
    return report


def main() -> int:
    run()
    print("fixtures -> tests/fixtures/validation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
