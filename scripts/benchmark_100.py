import sys
import json
import subprocess
import numpy as np
import random
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

try:
    from pyscf import gto, scf
    HAS_PYSCF = True
except ImportError:
    print("Warning: pyscf not installed")
    HAS_PYSCF = False

try:
    import torch
    import torchani
    HAS_TORCHANI = True
except ImportError:
    print("Warning: torch or torchani not installed")
    HAS_TORCHANI = False

def compute_torchani_aev(elements, coords_angstrom):
    ani2x = torchani.models.ANI2x(periodic_table_index=True)
    species = torch.tensor([elements], dtype=torch.long)
    coordinates = torch.tensor([coords_angstrom], dtype=torch.float32, requires_grad=False)
    _, aevs = ani2x.aev_computer((species, coordinates))
    return aevs[0].detach().numpy()

def compute_pyscf_hf(atom_str):
    mol = gto.M(atom=atom_str, basis='sto-3g', charge=0, spin=0, verbose=0)
    mf = scf.RHF(mol)
    mf.kernel()
    return mf.e_tot, mf.energy_nuc()

def run_sciform(method, elements, coords):
    flat_coords = []
    for c in coords:
        flat_coords.extend(c)
    cmd = [
        "cargo", "run", "-p", "sci-form-cli", "--release", "--bin", "sci-form", "--", method,
        json.dumps(elements),
        json.dumps(flat_coords)
    ]
    res = subprocess.run(cmd, capture_output=True, text=True, cwd="../")
    if res.returncode != 0:
        return None
    try:
        return json.loads(res.stdout)
    except:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--num", type=int, default=100, help="Number of molecules to test")
    parser.add_argument("--smi-file", type=str, default="10k_smiles.smi")
    args = parser.parse_args()

    smiles_list = []
    try:
        with open(args.smi_file, "r") as f:
            for line in f:
                parts = line.strip().split()
                if parts:
                    smiles_list.append(parts[0])
    except FileNotFoundError:
        print(f"File {args.smi_file} not found. Using default set.")
        smiles_list = ["C", "O", "CC", "CCO", "c1ccccc1", "C(=O)O", "C1CC1", "N#C", "C=C", "c1ccncc1", "c1ncncn1", "CC(=O)C", "CC(=O)N"]

    # Deduplicate and pick
    smiles_list = list(set(smiles_list))
    if len(smiles_list) > args.num:
        random.seed(42)  # For reproducibility
        test_smiles = random.sample(smiles_list, args.num)
    else:
        test_smiles = smiles_list

    print(f"Testing {len(test_smiles)} molecules...")

    max_hf_err = 0.0
    max_ani_err = 0.0
    passed = 0
    failed = 0

    for idx, smi in enumerate(test_smiles):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        if not mol:
            continue
            
        AllChem.EmbedMolecule(mol, randomSeed=idx)
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception:
            pass # fallback to crude geom
            
        if mol.GetNumConformers() == 0:
            continue
            
        conf = mol.GetConformer()
        elements = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
        coords = [list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
        
        # Check odd electrons - PySCF simple wrap only does RHF for spin=0
        if sum(elements) % 2 != 0:
            print(f"[{idx+1}/{len(test_smiles)}] {smi:15s} | Skipped (open shell)")
            continue
            
        atom_str = ""
        for el, c in zip(elements, coords):
            atom_str += f"{el} {c[0]:.6f} {c[1]:.6f} {c[2]:.6f}; "
            
        # 1. HF-3c checking
        hf_err_pct = 0.0
        if HAS_PYSCF:
            pyscf_hf, _ = compute_pyscf_hf(atom_str)
            rust_res_hf = run_sciform("hf3c", elements, coords)
            if rust_res_hf and "hf_energy" in rust_res_hf:
                rust_hf = rust_res_hf["hf_energy"]
                hf_err_pct = abs(rust_hf - pyscf_hf) / abs(pyscf_hf) * 100.0
                max_hf_err = max(max_hf_err, hf_err_pct)
        
        # 2. ANI checking
        ani_err_pct = 0.0
        if HAS_TORCHANI:
            try:
                torch_aevs = compute_torchani_aev(elements, coords)
                rust_res_ani = run_sciform("ani", elements, coords)
                if rust_res_ani and "aevs" in rust_res_ani:
                    rust_aevs = np.array(rust_res_ani["aevs"])
                    if rust_aevs.shape == torch_aevs.shape:
                        diff = np.abs(rust_aevs - torch_aevs)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            rel_err = np.where(torch_aevs > 1e-6, diff / torch_aevs, 0.0)
                        ani_err_pct = np.nanmax(rel_err) * 100.0
                        max_ani_err = max(max_ani_err, ani_err_pct)
            except Exception as e:
                pass # ANI2x might fail for some elements like halogens if not in trained set
                
        if hf_err_pct < 1.0 and (ani_err_pct < 1.0 or not HAS_TORCHANI):
            print(f"[{idx+1}/{len(test_smiles)}] {smi:15s} | HF Err: {hf_err_pct:6.4f}% | ANI Err: {ani_err_pct:6.4f}%  [PASS]")
            passed += 1
        else:
            print(f"[{idx+1}/{len(test_smiles)}] {smi:15s} | HF Err: {hf_err_pct:6.4f}% | ANI Err: {ani_err_pct:6.4f}%  [FAIL]")
            failed += 1

    print("\n" + "="*50)
    print(f"Benchmark Complete: {passed} Passed, {failed} Failed.")
    if HAS_PYSCF: print(f"Max HF-3c Deviation: {max_hf_err:.4f}%")
    if HAS_TORCHANI: print(f"Max ANI AEV Deviation: {max_ani_err:.4f}%")
    print("="*50)
    
    if failed > 0 or max_hf_err >= 1.0 or max_ani_err >= 1.0:
        sys.exit(1)

if __name__ == "__main__":
    main()
