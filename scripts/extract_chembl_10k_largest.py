#!/usr/bin/env python3
"""Extract the 10,000 largest molecules from ChemBL 36 by heavy atom count."""
import gzip
import sys
import re

def count_heavy_atoms_approx(smiles):
    """Fast approximate heavy atom count from SMILES string."""
    # Remove everything in brackets except the element letter
    # Count: uppercase letters = atoms, but subtract H counts
    count = 0
    i = 0
    s = smiles
    while i < len(s):
        c = s[i]
        if c == '[':
            # Bracketed atom - count as 1 heavy atom (unless it's [H] or [2H] etc.)
            j = s.find(']', i)
            if j == -1:
                break
            bracket = s[i+1:j]
            # Skip if it's just H with optional isotope/charge
            if re.match(r'^\d*H[+-]?\d*$', bracket):
                pass  # Don't count explicit H
            else:
                count += 1
            i = j + 1
        elif c.isupper():
            # Organic subset atom
            if c == 'H':
                i += 1
            elif c == 'B' and i+1 < len(s) and s[i+1] == 'r':
                count += 1
                i += 2
            elif c == 'C' and i+1 < len(s) and s[i+1] == 'l':
                count += 1
                i += 2
            else:
                count += 1
                # Handle lowercase second letter (e.g., not Br/Cl already handled)
                if i+1 < len(s) and s[i+1].islower():
                    i += 2
                else:
                    i += 1
        else:
            i += 1
    return count

# Read all SMILES with atom counts
print("Reading ChemBL 36...", file=sys.stderr)
molecules = []
with gzip.open('data/chembl_36_chemreps.txt.gz', 'rt') as f:
    header = f.readline()  # skip header
    for line_num, line in enumerate(f):
        if line_num % 500000 == 0 and line_num > 0:
            print(f"  processed {line_num} lines...", file=sys.stderr)
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        smiles = parts[1]
        # Skip very long SMILES (likely peptides/polymers) and very short ones
        if len(smiles) > 500 or len(smiles) < 5:
            continue
        # Skip SMILES with disconnected fragments (salts etc.)
        if '.' in smiles:
            continue
        natoms = count_heavy_atoms_approx(smiles)
        molecules.append((natoms, smiles, parts[0]))

print(f"Total valid molecules: {len(molecules)}", file=sys.stderr)

# Sort by heavy atom count descending, take top 10K
molecules.sort(key=lambda x: -x[0])
top10k = molecules[:10000]

print(f"Selected 10K largest: {top10k[0][0]} to {top10k[-1][0]} heavy atoms", file=sys.stderr)

# Write output
with open('data/chembl_10k_largest.smi', 'w') as f:
    for natoms, smiles, chembl_id in top10k:
        f.write(f"{smiles}\t{chembl_id}\t{natoms}\n")

print(f"Written to data/chembl_10k_largest.smi", file=sys.stderr)

# Also show distribution
from collections import Counter
sizes = Counter()
for natoms, _, _ in top10k:
    bucket = (natoms // 10) * 10
    sizes[bucket] += 1
print("\nSize distribution:", file=sys.stderr)
for bucket in sorted(sizes.keys()):
    print(f"  {bucket}-{bucket+9}: {sizes[bucket]}", file=sys.stderr)
