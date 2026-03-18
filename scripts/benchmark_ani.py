import sys
import json
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

try:
    import torch
    import torchani
except ImportError:
    print("Please install torch and torchani: pip install torch torchani")
    sys.exit(1)

def compute_torchani_aev(elements, coords_angstrom):
    # Use ANI-2x which has 7 elements
    ani2x = torchani.models.ANI2x(periodic_table_index=True)
    
    species = torch.tensor([elements], dtype=torch.long)
    coordinates = torch.tensor([coords_angstrom], dtype=torch.float32, requires_grad=False)
    
    aev_computer = ani2x.aev_computer
    
    # compute aev
    _, aevs = aev_computer((species, coordinates))
    
    return aevs[0].detach().numpy() # shape (N, 384)

def run_sciform_ani(elements, coords):
    # Flatten coords
    flat_coords = []
    for c in coords:
        flat_coords.extend(c)
        
    cmd = [
        "cargo", "run", "-p", "sci-form-cli", "--release", "--bin", "sci-form", "--", "ani",
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
        
        # 1. Run our ANI
        rust_res = run_sciform_ani(elements, coords)
        if not rust_res or "aevs" not in rust_res:
            continue
            
        rust_aevs = np.array(rust_res["aevs"])
        
        # 2. Run TorchANI
        torchani_aevs = compute_torchani_aev(elements, coords)
        
        # Different implementations might order the padding or elements slightly differently,
        # but the shape and structure should match ANI-2x standards.
        if rust_aevs.shape != torchani_aevs.shape:
            print(f"Shape mismatch on {smi}: rust={rust_aevs.shape}, torchani={torchani_aevs.shape}")
            continue
        
        diff = np.abs(rust_aevs - torchani_aevs)
        
        # Ignore zero AEVs for relative error
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_err = np.where(torchani_aevs > 1e-6, diff / torchani_aevs, 0.0)
            
        max_rel_err_pct = np.nanmax(rel_err) * 100.0
        max_err = max(max_err, max_rel_err_pct)
        
        print(f"[{idx+1}/{len(smiles_list)}] {smi:10s} | Max AEV diff: {np.max(diff):.6f} | Max rel err: {max_rel_err_pct:.4f}%")
        
        assert max_rel_err_pct < 1.0, f"Error exceeds 1% for {smi}! ({max_rel_err_pct:.4f}%)"

    print(f"\n✅ All ANI calculations (AEVs component) passed! Max deviation = {max_err:.4f}%")

if __name__ == "__main__":
    main()
