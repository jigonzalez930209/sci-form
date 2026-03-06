import rdkit.Chem as Chem
from rdkit.Chem import AllChem
import json
import numpy as np

def test_rdkit_variance():
    # Load the same 100 molecules
    with open("tests/fixtures/reference_coords.json", "r") as f:
        data = json.load(f)
    
    rmsd_list = []
    
    for item in data[:100]:
        m = Chem.MolFromSmiles(item['smiles'])
        m = Chem.AddHs(m)
        
        # We need to ensure stereochem and conformers are handled exactly
        # Generate two conformers
        res = AllChem.EmbedMultipleConfs(m, numConfs=2, params=AllChem.ETKDGv3())
        
        if len(res) == 2:
            # Calculate RMSD without alignment? ETKDG generates them aligned?
            # RDKit GetBestRMS align them
            rmsd = AllChem.GetBestRMS(m, m, res[0], res[1])
            rmsd_list.append(rmsd)
            
    avg_rmsd = sum(rmsd_list) / len(rmsd_list)
    print(f"Average internal RDKit ETKDGv3 RMSD variance: {avg_rmsd:.3f} A")

if __name__ == "__main__":
    test_rdkit_variance()
