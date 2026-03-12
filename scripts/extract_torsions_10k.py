#!/usr/bin/env python3
"""Extract experimental torsion parameters from RDKit for all 10k SMILES.

Reads from scripts/10k_smiles.smi and generates CSD-derived M6 Fourier
torsion preferences for each molecule. Merges with existing torsion_params.json.
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)


def main():
    # Load existing torsion params
    existing_file = 'tests/fixtures/torsion_params.json'
    try:
        with open(existing_file) as f:
            results = json.load(f)
        print(f"Loaded {len(results)} existing torsion entries", file=sys.stderr)
    except FileNotFoundError:
        results = {}

    # Load 10k SMILES
    with open('scripts/10k_smiles.smi') as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    print(f"Processing {len(smiles_list)} SMILES...", file=sys.stderr)

    n_new = 0
    n_skip = 0
    n_fail = 0

    for i, smiles in enumerate(smiles_list):
        if smiles in results:
            n_skip += 1
            continue

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                n_fail += 1
                continue
            mol = Chem.AddHs(mol)
            torsions_raw = rdDistGeom.GetExperimentalTorsions(mol)
            torsion_list = []
            for t in torsions_raw:
                torsion_list.append({
                    'atoms': list(t['atomIndices']),
                    'v': list(t['V']),
                    'signs': list(t['signs']),
                })
            if torsion_list:
                results[smiles] = torsion_list
                n_new += 1
        except Exception as e:
            n_fail += 1
            if n_fail <= 10:
                print(f"Error on mol {i} ({smiles}): {e}", file=sys.stderr)

        if (i + 1) % 2000 == 0:
            print(f"  {i + 1}/{len(smiles_list)}... ({n_new} new, {n_skip} existing, {n_fail} fail)", file=sys.stderr)

    with open(existing_file, 'w') as f:
        json.dump(results, f)

    print(f"\nDone: {n_new} new + {n_skip} existing = {len(results)} total torsion entries", file=sys.stderr)
    print(f"Failures: {n_fail}", file=sys.stderr)


if __name__ == '__main__':
    main()
