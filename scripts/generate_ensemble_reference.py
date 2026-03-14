#!/usr/bin/env python3
"""
Generate RDKit multi-seed ensemble reference data.
For each molecule, generates conformers with N seeds and saves all coordinates.
This allows fair comparison: our conformer vs closest RDKit conformer (min-RMSD).

Output format (JSON):
[{
  "smiles": "...",
  "atoms": [{element, hybridization, formal_charge}],
  "bonds": [{start, end, order}],
  "torsions": [{atoms, v, signs}],
  "conformers": {
    "42": [[x,y,z], [x,y,z], ...],
    "0": [[x,y,z], ...],
    ...
  }
}]
"""
import json, sys, os
from multiprocessing import Pool
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

SEEDS = [42, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

def hyb_str(h):
    from rdkit.Chem import rdchem
    m = {
        rdchem.HybridizationType.SP: "SP",
        rdchem.HybridizationType.SP2: "SP2",
        rdchem.HybridizationType.SP3: "SP3",
        rdchem.HybridizationType.SP3D: "SP3D",
        rdchem.HybridizationType.SP3D2: "SP3D2",
    }
    return m.get(h, "Unknown")

def bond_order_str(bt):
    from rdkit.Chem import rdchem
    m = {
        rdchem.BondType.SINGLE: "SINGLE",
        rdchem.BondType.DOUBLE: "DOUBLE",
        rdchem.BondType.TRIPLE: "TRIPLE",
        rdchem.BondType.AROMATIC: "AROMATIC",
    }
    return m.get(bt, "SINGLE")

def process_mol(args):
    smi, idx = args
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    
    # Get topology info
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append({
            "element": atom.GetAtomicNum(),
            "hybridization": hyb_str(atom.GetHybridization()),
            "formal_charge": atom.GetFormalCharge(),
        })
    
    bonds = []
    for bond in mol.GetBonds():
        bonds.append({
            "start": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "order": bond_order_str(bond.GetBondType()),
        })
    
    # Get CSD torsion info
    torsions = []
    try:
        expTorsionAngles = rdDistGeom.GetExperimentalTorsions(mol)
        for torsion in expTorsionAngles:
            atom_indices = list(torsion[0])
            signs_and_v = torsion[1]
            signs = [int(signs_and_v[i][0]) for i in range(6)]
            v = [float(signs_and_v[i][1]) for i in range(6)]
            torsions.append({
                "atoms": atom_indices,
                "signs": signs,
                "v": v,
            })
    except Exception:
        pass
    
    # Generate conformers for all seeds
    conformers = {}
    for seed in SEEDS:
        try:
            mol2 = Chem.RWMol(mol)
            res = AllChem.EmbedMolecule(mol2, randomSeed=seed)
            if res == 0:
                conf = mol2.GetConformer()
                coords = []
                for i in range(mol2.GetNumAtoms()):
                    pos = conf.GetAtomPosition(i)
                    coords.append([round(pos.x, 6), round(pos.y, 6), round(pos.z, 6)])
                conformers[str(seed)] = coords
        except Exception:
            pass
    
    if not conformers:
        return None
    
    return {
        "smiles": smi,
        "atoms": atoms,
        "bonds": bonds,
        "torsions": torsions,
        "conformers": conformers,
    }

def main():
    # Load GDB-20 reference
    src = sys.argv[1] if len(sys.argv) > 1 else "tests/fixtures/gdb20_reference.json"
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
    
    print(f"Loading molecules from {src}...")
    with open(src) as f:
        ref_mols = json.load(f)
    
    smiles_list = [(m["smiles"], i) for i, m in enumerate(ref_mols[:limit])]
    print(f"Processing {len(smiles_list)} molecules with {len(SEEDS)} seeds each...")
    
    with Pool() as pool:
        results = pool.map(process_mol, smiles_list)
    
    results = [r for r in results if r is not None]
    
    # Stats
    total_confs = sum(len(r["conformers"]) for r in results)
    avg_seeds = total_confs / len(results) if results else 0
    print(f"Success: {len(results)}/{len(smiles_list)} molecules, "
          f"avg {avg_seeds:.1f} conformers per molecule")
    
    out_path = "tests/fixtures/gdb20_ensemble.json"
    with open(out_path, 'w') as f:
        json.dump(results, f)
    
    sz = os.path.getsize(out_path) / 1e6
    print(f"Saved to {out_path} ({sz:.1f} MB)")

if __name__ == '__main__':
    main()
