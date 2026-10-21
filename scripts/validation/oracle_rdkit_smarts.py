#!/usr/bin/env python3
"""RDKit oracle for sci-form's SMARTS substructure-matching gate (F75).

Compares `sci_form::smarts::substruct_match` against RDKit's
`Mol.GetSubstructMatches` for the 837 CSD torsion patterns shipped in
`tests/fixtures/smarts_patterns.txt` on a frozen molecule corpus.

The Rust side (`tests/rdkit_parity.rs`) recomputes the sci-form matches from
the same fixture and asserts the per-pattern agreement rate, so CI does not
need RDKit.

Usage:
    python scripts/validation/oracle_rdkit_smarts.py [all]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from common import REPO, Report, emit_fixture_report

try:
    from rdkit import Chem
except ImportError as exc:  # pragma: no cover
    print(f"RDKit is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

RDKIT_VERSION = Chem.rdBase.rdkitVersion
REFERENCE = f"RDKit {RDKIT_VERSION}"

PATTERN_FILE = REPO / "tests" / "fixtures" / "smarts_patterns.txt"

#: Frozen molecule corpus.  Chosen so that the CSD patterns have something to
#: match: aromatic rings (with and without heteroatoms), carbonyls, amides,
#: esters, sulfoxides and a sugar.
MOLECULE_CORPUS = [
    "c1ccccc1",
    "Cc1ccccc1",
    "c1ccncc1",
    "c1cc[nH]c1",
    "c1ccsc1",
    "CC(=O)NC",
    "COC(C)=O",
    "CS(C)=O",
    "CC(C)(C)c1ccccc1",
    "OCC1OC(O)C(O)C(O)C1O",
]

#: Minimum fraction of (pattern, molecule) pairs whose canonical match sets
#: must agree exactly.  Measured value is recorded in the fixture notes.
MIN_AGREEMENT = 0.90


def read_patterns() -> list[str]:
    patterns: list[str] = []
    for line in PATTERN_FILE.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        patterns.append(line.split("\t")[0].strip())
    return patterns


def canonical_matches(matches) -> list[list[int]]:
    """Canonical form: each match reduced to a sorted tuple of atom indices."""
    out = {tuple(sorted(m)) for m in matches}
    return [list(t) for t in sorted(out)]


def run_smarts() -> Report:
    patterns = read_patterns()
    molecules = []
    for smi in MOLECULE_CORPUS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        molecules.append((smi, Chem.AddHs(mol)))

    cases = []
    per_pattern = []
    for pidx, smarts in enumerate(patterns):
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            cases.append({"pattern_index": pidx, "unsupported": True, "matches": []})
            per_pattern.append({"pattern_index": pidx, "n_rdkit_matched": 0})
            continue
        p_matched = 0
        for midx, (smi, mol) in enumerate(molecules):
            raw = mol.GetSubstructMatches(patt, uniquify=True, useChirality=False, maxMatches=100000)
            canon = canonical_matches(raw)
            cases.append(
                {
                    "pattern_index": pidx,
                    "molecule_index": midx,
                    "smiles": smi,
                    "n_rdkit_raw": len(raw),
                    "matches": canon,
                }
            )
            if canon:
                p_matched += 1
        per_pattern.append({"pattern_index": pidx, "n_rdkit_matched": p_matched})

    payload = {
        "reference": f"{REFERENCE} Mol.GetSubstructMatches(uniquify=True, useChirality=False)",
        "rdkit_api": "Chem.Mol.GetSubstructMatches(patt, uniquify=True, useChirality=False)",
        "sci_form_api": "sci_form::smarts::substruct_match(mol, &parse_smarts(pattern)?)",
        "canonical_form": "each match reduced to the sorted tuple of matched atom indices",
        "min_agreement": MIN_AGREEMENT,
        "generated": time.strftime("%Y-%m-%d"),
        "patterns": patterns,
        "molecules": [smi for smi, _ in molecules],
        "cases": cases,
        "per_pattern": per_pattern,
    }
    report = Report(
        ficha="F75",
        observable="smarts_match_agreement",
        reference_tool=f"{REFERENCE} GetSubstructMatches",
        units="fraction",
        tolerance=MIN_AGREEMENT,
        level="L2",
        n_cases=len(cases),
        metric="agreement_geq",
        value=float("nan"),
        notes=(
            "837 CSD torsion patterns x 10 molecules; the canonical match sets "
            "must agree exactly for at least 90% of the (pattern, molecule) pairs"
        ),
    )
    emit_fixture_report("smarts_matches_rdkit.json", payload, report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", nargs="?", default="all", choices=["all", "smarts"])
    args = parser.parse_args()
    run_smarts()
    print(f"\nfixtures -> {REPO / 'tests' / 'fixtures' / 'validation'}")
    print(f"summary  -> {REPO / 'reports' / 'validation' / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
