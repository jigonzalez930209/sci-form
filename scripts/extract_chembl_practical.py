"""Extract 10K diverse ChemBL molecules in a practical heavy-atom range (30-100).
These are genuinely large drug-like molecules where distance geometry is feasible.
"""
import gzip
import random
import sys

chembl_file = "data/chembl_36_chemreps.txt.gz"
output_file = "data/chembl_10k_practical.smi"

print(f"Reading {chembl_file}...")
molecules = []
with gzip.open(chembl_file, 'rt') as f:
    header = f.readline()
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        smiles = parts[1]
        chembl_id = parts[0]
        # Count heavy atoms from SMILES (approximate)
        heavy = 0
        i = 0
        in_bracket = False
        while i < len(smiles):
            c = smiles[i]
            if c == '[':
                in_bracket = True
                i += 1
                continue
            if c == ']':
                in_bracket = False
                i += 1
                continue
            if in_bracket:
                if c.isalpha() and c != 'H' and c != 'h':
                    heavy += 1
                    # Skip rest of element symbol
                    while i + 1 < len(smiles) and smiles[i+1].islower() and smiles[i+1] != 'h':
                        i += 1
                i += 1
                continue
            if c.isalpha():
                if c not in ('h',):  # lowercase h is hydrogen in aromatic
                    heavy += 1
                    # Two-letter elements
                    if i + 1 < len(smiles) and smiles[i+1].islower() and smiles[i+1] not in ('c','n','o','s','p'):
                        if c == 'C' and smiles[i+1] == 'l':
                            i += 1  # Cl
                        elif c == 'B' and smiles[i+1] == 'r':
                            i += 1  # Br
            i += 1

        if 30 <= heavy <= 100:
            molecules.append((smiles, chembl_id, heavy))

print(f"Found {len(molecules)} molecules with 30-100 heavy atoms")

# Stratified sampling: equal from each bucket
buckets = {}
for m in molecules:
    bucket = (m[2] // 10) * 10  # 30, 40, 50, 60, 70, 80, 90, 100
    buckets.setdefault(bucket, []).append(m)

target = 10000
per_bucket = target // len(buckets)
selected = []
random.seed(42)

for bucket in sorted(buckets.keys()):
    pool = buckets[bucket]
    n = min(per_bucket, len(pool))
    selected.extend(random.sample(pool, n))
    print(f"  Bucket {bucket}-{bucket+9}: {len(pool)} available, selected {n}")

# Fill remaining from largest buckets
remaining = target - len(selected)
if remaining > 0:
    selected_ids = {m[1] for m in selected}
    all_remaining = []
    for bucket in sorted(buckets.keys(), reverse=True):
        for m in buckets[bucket]:
            if m[1] not in selected_ids:
                all_remaining.append(m)
    extra = random.sample(all_remaining, min(remaining, len(all_remaining)))
    selected.extend(extra)

# Sort by heavy atom count descending
selected.sort(key=lambda m: -m[2])

with open(output_file, 'w') as f:
    for smiles, chembl_id, heavy in selected:
        f.write(f"{smiles}\t{chembl_id}\t{heavy}\n")

print(f"\nWrote {len(selected)} molecules to {output_file}")
print(f"Heavy atom range: {selected[-1][2]}-{selected[0][2]}")
