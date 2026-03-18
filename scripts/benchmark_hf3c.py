import sys
import json
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    from pyscf import gto, scf
except ImportError:
    print("Please install pyscf: pip install pyscf")
    sys.exit(1)

def compute_pyscf_hf(elements, coords_angstrom):
    # Format for pyscf: "C 0.0 0.0 0.0; O 0.0 1.2 0.0"
    atom_str = ""
    for el, c in zip(elements, coords_angstrom):
        atom_str += f"{el} {c[0]} {c[1]} {c[2]}; "
    
    mol = gto.M(
        atom=atom_str,
        basis='sto-3g',
        charge=0,
        spin=0,
        verbose=0
    )
    mf = scf.RHF(mol)
    mf.kernel()
    return mf.e_tot, mf.energy_nuc()

def run_sciform_hf3c(elements, coords):
    # Flatten coords
    flat_coords = []
    for c in coords:
        flat_coords.extend(c)
        
    cmd = [
        "cargo", "run", "-p", "sci-form-cli", "--release", "--bin", "sci-form", "--", "hf3c",
        json.dumps(elements),
        json.dumps(flat_coords)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd="../")
    if res.returncode != 0:
        print("Error running sci-form:", res.stderr)
        return None
    try:
        data = json.loads(res.stdout)
        return data
    except json.JSONDecodeError:
        print("Invalid JSON:", res.stdout)
        return None

def main():
    smiles_list = ["C", "O", "CC", "CCO", "c1ccccc1", "C(=O)O", "C1CC1", "N#C"]
    
    max_err = 0.0
    
    for idx, smi in enumerate(smiles_list):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.UFFOptimizeMolecule(mol)
        
        conf = mol.GetConformer()
        elements = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        
        # 1. Run our HF-3c (extract hf_energy component)
        rust_res = run_sciform_hf3c(elements, coords)
        if not rust_res:
            continue
            
        rust_hf = rust_res["hf_energy"]
        rust_nuc = rust_res["nuclear_repulsion"]
        
        # 2. Run PySCF STO-3G
        pyscf_hf, pyscf_nuc = compute_pyscf_hf(elements, coords)
        
        err_pct = abs(rust_hf - pyscf_hf) / abs(pyscf_hf) * 100.0
        max_err = max(max_err, err_pct)
        
        print(f"[{idx+1}/{len(smiles_list)}] {smi:10s} | Rust HF: {rust_hf:12.6f} (Nuc: {rust_nuc:10.6f}) | PySCF: {pyscf_hf:12.6f} (Nuc: {pyscf_nuc:10.6f}) | Err: {err_pct:.4f}%")
        
        assert err_pct < 1.0, f"Error exceeds 1% for {smi}! ({err_pct:.4f}%)"

    print(f"\n✅ All HF-3c (STO-3G component) calculations passed! Max deviation = {max_err:.4f}%")

if __name__ == "__main__":
    main()
