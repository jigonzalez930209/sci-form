#!/usr/bin/env python3
"""Generate RDKit reference coordinates + CSD torsion params for GDB20 50k molecules.

Uses multiprocessing to parallelize across all CPU cores.
Output: tests/fixtures/gdb20_reference.json.gz (coords + torsions)
"""

import sys
import time
from multiprocessing import Pool, cpu_count
from fixture_io import dump_json_gz

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit import RDLogger

RDLogger.logger().setLevel(RDLogger.ERROR)


def process_smiles(smiles):
    """Process a single SMILES: embed + extract CSD torsions."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = Chem.AddHs(mol)

        # Embed with seed=42 (pure ETKDG, no MMFF)
        res = AllChem.EmbedMolecule(mol, randomSeed=42)
        if res != 0:
            return None

        conf = mol.GetConformer()

        # Extract atoms
        atoms = []
        for i, atom in enumerate(mol.GetAtoms()):
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "element": atom.GetAtomicNum(),
                "x": round(pos.x, 6),
                "y": round(pos.y, 6),
                "z": round(pos.z, 6),
                "formal_charge": atom.GetFormalCharge(),
                "hybridization": str(atom.GetHybridization()),
            })

        # Extract bonds
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                "start": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "order": str(bond.GetBondType()),
            })

        # Extract CSD experimental torsions
        torsions = []
        try:
            torsions_raw = rdDistGeom.GetExperimentalTorsions(mol)
            for t in torsions_raw:
                torsions.append({
                    "atoms": list(t["atomIndices"]),
                    "v": [round(x, 8) for x in t["V"]],
                    "signs": list(t["signs"]),
                })
        except Exception:
            pass

        return {
            "smiles": smiles,
            "atoms": atoms,
            "bonds": bonds,
            "torsions": torsions,
        }
    except Exception:
        return None


def main():
    input_file = "GDB20.50000.smi"
    output_file = "tests/fixtures/gdb20_reference.json.gz"

    with open(input_file) as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    total = len(smiles_list)
    ncpu = cpu_count()
    print(f"Processing {total} SMILES with {ncpu} cores...", file=sys.stderr)

    start = time.time()

    with Pool(ncpu) as pool:
        results = []
        for i, result in enumerate(pool.imap(process_smiles, smiles_list, chunksize=100)):
            if result is not None:
                results.append(result)
            if (i + 1) % 5000 == 0:
                elapsed = time.time() - start
                rate = (i + 1) / elapsed
                print(
                    f"  {i+1}/{total} ({len(results)} ok, {i+1-len(results)} fail) "
                    f"{rate:.0f} mol/s, ETA {(total-i-1)/rate:.0f}s",
                    file=sys.stderr,
                )

    elapsed = time.time() - start
    print(
        f"\nDone: {len(results)}/{total} molecules in {elapsed:.1f}s "
        f"({total/elapsed:.0f} mol/s)",
        file=sys.stderr,
    )
    print(f"Embed failures: {total - len(results)}", file=sys.stderr)

    dump_json_gz(output_file, results)
    print(f"Saved to {output_file} ({len(results)} molecules)", file=sys.stderr)


if __name__ == "__main__":
    main()
