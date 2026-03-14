#!/usr/bin/env python3
"""Log which attempt index RDKit uses for a specific molecule.

RDKit's EmbedMolecule internally retries with different random distance samples.
This script exposes that information by reimplementing the retry loop with logging.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit import RDLogger
import sys

RDLogger.logger().setLevel(RDLogger.ERROR)

# Molecules to test (worst failures)
SMILES_LIST = [
    "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O",
    "CCNCCC1(C(C)(N)CC#N)CCC(O)C1(C)CC",
    "C#CCNCC1(C)C2CCC1(CCO)C1(CCCO1)C2",
    "CCC1(CCC(C)(C)O)CCC(C)Nc2cncn21",
]

def test_molecule(smiles):
    """Use RDKit's standard embedding to see which attempt succeeds."""
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    
    # Use EmbedMolecule with seed=42 and capture verbose output
    ps = rdDistGeom.EmbedParameters()
    ps.randomSeed = 42
    ps.useExpTorsionAnglePrefs = True
    ps.useBasicKnowledge = True
    ps.enforceChirality = True
    ps.useRandomCoords = False
    ps.numThreads = 1

    # RDKit's default maxIterations = 10 * numAtoms
    n = mol.GetNumAtoms()
    ps.maxIterations = 10 * n
    
    # We can only get the final result, not per-attempt logging
    # But we can try to match by using different seeds and seeing when embedding succeeds
    confId = AllChem.EmbedMolecule(mol, ps)
    print(f"SMILES: {smiles}")
    print(f"  Atoms: {n}, confId: {confId}")
    
    if confId >= 0:
        conf = mol.GetConformer(confId)
        # Print first few atom positions
        for i in range(min(3, n)):
            pos = conf.GetAtomPosition(i)
            print(f"  Atom {i}: ({pos.x:.6f}, {pos.y:.6f}, {pos.z:.6f})")

if __name__ == "__main__":
    for smi in SMILES_LIST:
        test_molecule(smi)
        print()
