import urllib.request
import os

def download_chembl_smiles():
    print("Downloading ChEMBL dataset sample... This might take a moment.")
    # A known good simple SMILES dataset for benchmarking from RDKit regression tests
    url = "https://raw.githubusercontent.com/rdkit/rdkit/master/Projects/DbCLI/testData/pubchem.200.txt"
    # Actually, let's use a bigger guaranteed RDKit dataset if available, but pubchem.200 only has 200.
    # Let's try to get one of the larger ML datasets like ESOL or FreeSolv, they have ~1-4k.
    # Instead, let's generate exactly 10,000 random-ish SMILES by repeating and mutating, 
    # or just use a known robust URL like the deepchem datasets.
    
    url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/qm9.csv"
    output_file = "scripts/10k_smiles.smi"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    try:
        print(f"Fetching from {url}...")
        response = urllib.request.urlopen(url)
        data = response.read().decode('utf-8')
        
        lines = data.split('\n')
        
        with open(output_file, 'w') as f:
            count = 0
            # QM9 has header: "mol_id,smiles,A,B,C,mu,alpha,homo,lumo,gap,r2,zpve,u0,u298,h298,g298,cv"
            for line in lines[1:]: # Skip header
                line = line.strip()
                if not line:
                    continue
                parts = line.split(',')
                if len(parts) > 1:
                    smiles = parts[1].strip('"')
                    f.write(smiles + '\n')
                    count += 1
                if count >= 10000:
                    break
                    
        print(f"Successfully extracted {count} SMILES to {output_file}")
    except Exception as e:
        print(f"Error downloading dataset: {e}")
        
if __name__ == "__main__":
    download_chembl_smiles()
