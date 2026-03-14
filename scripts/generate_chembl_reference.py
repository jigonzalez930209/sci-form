#!/usr/bin/env python3
"""Generate RDKit reference conformers for the ChemBL 10K largest molecules.

Produces a JSON file with:
  - Atom info (element, hybridization, formal_charge, x/y/z)
  - Bond info (start, end, order)
  - CSD torsion contributions (atom indices, V values, signs)
  - Ring info for each molecule

This is the oracle data for benchmarking our Rust implementation.
"""
import json
import sys
import time
import os

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
except ImportError:
    print("RDKit not available. Install with: pip install rdkit", file=sys.stderr)
    sys.exit(1)

def get_hybridization_str(hyb):
    """Convert RDKit hybridization to string."""
    hyb_map = {
        Chem.rdchem.HybridizationType.SP: "SP",
        Chem.rdchem.HybridizationType.SP2: "SP2",
        Chem.rdchem.HybridizationType.SP3: "SP3",
        Chem.rdchem.HybridizationType.SP3D: "SP3D",
        Chem.rdchem.HybridizationType.SP3D2: "SP3D2",
    }
    return hyb_map.get(hyb, "UNKNOWN")

def get_bond_order_str(bt):
    """Convert RDKit bond type to string."""
    bt_map = {
        Chem.rdchem.BondType.SINGLE: "SINGLE",
        Chem.rdchem.BondType.DOUBLE: "DOUBLE",
        Chem.rdchem.BondType.TRIPLE: "TRIPLE",
        Chem.rdchem.BondType.AROMATIC: "AROMATIC",
    }
    return bt_map.get(bt, "SINGLE")

def get_csd_torsions(mol, conf_id=0):
    """Extract CSD torsion contributions from RDKit."""
    try:
        torsions_raw = rdDistGeom.GetExperimentalTorsions(mol)
        torsions = []
        for t in torsions_raw:
            torsions.append({
                "atoms": list(t["atomIndices"]),
                "v": [round(x, 8) for x in t["V"]],
                "signs": list(t["signs"]),
            })
        return torsions
    except:
        return []

def process_molecule(smiles, seed=42, timeout_sec=120):
    """Process a single molecule and return reference data."""
    import signal
    
    def handler(signum, frame):
        raise TimeoutError("Embedding timed out")
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        
        # Generate 3D conformer using ETKDG
        params = rdDistGeom.ETKDGv3()
        params.randomSeed = seed
        params.useSmallRingTorsions = True
        params.useMacrocycleTorsions = True
        params.useBasicKnowledge = True
        params.enforceChirality = True
        params.numThreads = 1
        params.maxIterations = 0  # auto
        
        # Set alarm for timeout
        old_handler = signal.signal(signal.SIGALRM, handler)
        signal.alarm(timeout_sec)
        
        try:
            conf_id = AllChem.EmbedMolecule(mol, params)
        finally:
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
        if conf_id < 0:
            return None
        
        conf = mol.GetConformer(conf_id)
        
        # Extract atom info
        atoms = []
        for i in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            atoms.append({
                "element": atom.GetAtomicNum(),
                "x": round(pos.x, 6),
                "y": round(pos.y, 6),
                "z": round(pos.z, 6),
                "formal_charge": atom.GetFormalCharge(),
                "hybridization": get_hybridization_str(atom.GetHybridization()),
            })
        
        # Extract bond info
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                "start": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "order": get_bond_order_str(bond.GetBondType()),
            })
        
        # Extract CSD torsion contributions
        torsions = get_csd_torsions(mol, conf_id)
        
        return {
            "smiles": smiles,
            "atoms": atoms,
            "bonds": bonds,
            "torsions": torsions,
        }
    except Exception as e:
        return None

def main():
    input_file = "data/chembl_10k_largest.smi"
    output_file = "tests/fixtures/chembl_10k_reference.json"
    
    limit = int(os.environ.get("CHEMBL_LIMIT", "10000"))
    
    print(f"Reading {input_file}...", file=sys.stderr)
    molecules = []
    with open(input_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                molecules.append((parts[0], parts[1]))
    
    molecules = molecules[:limit]
    print(f"Processing {len(molecules)} molecules...", file=sys.stderr)
    
    results = []
    failed = 0
    start = time.time()
    
    for idx, (smiles, chembl_id) in enumerate(molecules):
        if idx % 100 == 0 and idx > 0:
            elapsed = time.time() - start
            rate = idx / elapsed
            eta = (len(molecules) - idx) / rate if rate > 0 else 0
            print(f"  [{idx}/{len(molecules)}] {len(results)} ok, {failed} fail, "
                  f"{rate:.1f} mol/s, ETA {eta:.0f}s", file=sys.stderr)
        
        result = process_molecule(smiles)
        if result is not None:
            result["chembl_id"] = chembl_id
            results.append(result)
        else:
            failed += 1
    
    elapsed = time.time() - start
    print(f"\nDone: {len(results)} ok, {failed} failed in {elapsed:.1f}s "
          f"({len(results)/elapsed:.1f} mol/s)", file=sys.stderr)
    
    # Write output
    with open(output_file, 'w') as f:
        json.dump(results, f)
    
    print(f"Written to {output_file} ({os.path.getsize(output_file) / 1024 / 1024:.1f} MB)",
          file=sys.stderr)

if __name__ == "__main__":
    main()
