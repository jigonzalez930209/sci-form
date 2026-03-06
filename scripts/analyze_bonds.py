import json

def analyze():
    with open("tests/fixtures/rdkit_10k_reference.json", "r") as f:
        data = json.load(f)

    bond_lengths = {}
    
    for mol in data:
        atoms = mol["atoms"]
        bonds = mol["bonds"]
        for b in bonds:
            a1 = atoms[b["start"]]
            a2 = atoms[b["end"]]
            dx = a1["x"] - a2["x"]
            dy = a1["y"] - a2["y"]
            dz = a1["z"] - a2["z"]
            dist = (dx*dx + dy*dy + dz*dz)**0.5
            
            e1 = a1["element"]
            e2 = a2["element"]
            if e1 > e2:
                e1, e2 = e2, e1
                
            key = f"{e1}-{e2}-{b['order']}"
            if key not in bond_lengths:
                bond_lengths[key] = []
            bond_lengths[key].append(dist)

    print("Average 1-2 Distances:")
    for k, v in bond_lengths.items():
        if len(v) > 10:
            avg = sum(v) / len(v)
            print(f'"{k}" => {avg:.3f}, // count: {len(v)}')

if __name__ == "__main__":
    analyze()
