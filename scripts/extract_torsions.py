#!/usr/bin/env python3
"""Extract experimental torsion parameters from RDKit for test molecules.

Uses RDKit's GetExperimentalTorsions API to get the CSD-derived torsion
preferences (M6 Fourier coefficients) for each rotatable bond.
Outputs a JSON file that the Rust test can consume.

Atom ordering matches reference_coords.json because both use MolFromSmiles+AddHs.
"""

import json
import sys

from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)


def main():
    ref_file = 'tests/fixtures/reference_coords.json'
    with open(ref_file) as f:
        molecules = json.load(f)

    results = {}
    n_with_torsions = 0

    for i, mol_data in enumerate(molecules):
        smiles = mol_data['smiles']
        if smiles in results:
            continue  # skip duplicate SMILES
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
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
                n_with_torsions += 1
        except Exception as e:
            print(f"Error on mol {i} ({smiles}): {e}", file=sys.stderr)

        if (i + 1) % 1000 == 0:
            print(f"Processed {i + 1}/{len(molecules)}...", file=sys.stderr)

    output_file = 'tests/fixtures/torsion_params.json'
    with open(output_file, 'w') as f:
        json.dump(results, f)

    print(f"Extracted torsions for {n_with_torsions}/{len(molecules)} molecules -> {output_file}", file=sys.stderr)


if __name__ == '__main__':
    main()


if __name__ == '__main__':
    main()
