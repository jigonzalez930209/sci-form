#!/usr/bin/env python3
"""Compare 3D FF constraint counts between our code and RDKit.

For a specific molecule, dumps the number of each constraint type
from RDKit's construct3DForceField, so we can verify our Rust code
produces the same set.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
from rdkit import RDLogger
from rdkit.Chem import rdMolDescriptors
import rdkit.Geometry.rdGeometry as Geometry
from rdkit.ForceField import rdForceField

RDLogger.logger().setLevel(RDLogger.ERROR)

# Test molecule - worst failing one
SMILES = "C#CC(C)(COCCC1CCCNCC1)OCC(C)=O"

def main():
    mol = Chem.MolFromSmiles(SMILES)
    mol = Chem.AddHs(mol)
    n = mol.GetNumAtoms()
    
    print(f"SMILES: {SMILES}")
    print(f"Atoms: {n}")
    print(f"Bonds: {mol.GetNumBonds()}")
    
    # Count structural features
    n_aromatic = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    n_sp2 = sum(1 for a in mol.GetAtoms() if a.GetHybridization() == Chem.HybridizationType.SP2)
    n_sp = sum(1 for a in mol.GetAtoms() if a.GetHybridization() == Chem.HybridizationType.SP)
    n_sp3 = sum(1 for a in mol.GetAtoms() if a.GetHybridization() == Chem.HybridizationType.SP3)
    ring_info = mol.GetRingInfo()
    
    print(f"SP: {n_sp}, SP2: {n_sp2}, SP3: {n_sp3}, Aromatic: {n_aromatic}")
    print(f"Rings: {ring_info.NumRings()}, Ring sizes: {[len(r) for r in ring_info.AtomRings()]}")
    
    # Get experimental torsions
    torsions = rdDistGeom.GetExperimentalTorsions(mol)
    print(f"\nCSD Torsions: {len(torsions)}")
    for i, t in enumerate(torsions):
        print(f"  T{i}: atoms={list(t['atomIndices'])}, bond={t['bondIndex']}")
    
    # Count 1-3 angles
    n_angles = 0
    n_sp_angles = 0
    n_improper_angles = 0
    n_normal_angles = 0
    for atom in mol.GetAtoms():
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for a_idx in range(len(neighbors)):
            for b_idx in range(a_idx + 1, len(neighbors)):
                n_angles += 1
                center = atom.GetIdx()
                
                # Check for SP/linear
                bond_a = mol.GetBondBetweenAtoms(center, neighbors[a_idx])
                bond_b = mol.GetBondBetweenAtoms(center, neighbors[b_idx])
                is_triple_a = bond_a.GetBondType() == Chem.BondType.TRIPLE
                is_triple_b = bond_b.GetBondType() == Chem.BondType.TRIPLE
                is_double_a = bond_a.GetBondType() == Chem.BondType.DOUBLE
                is_double_b = bond_b.GetBondType() == Chem.BondType.DOUBLE
                
                is_sp = is_triple_a or is_triple_b or (is_double_a and is_double_b and len(neighbors) == 2)
                
                # Check for improper
                elem = atom.GetAtomicNum()
                hyb = atom.GetHybridization()
                is_improper = (elem in [6, 7, 8] and hyb == Chem.HybridizationType.SP2 and len(neighbors) == 3)
                
                if is_sp:
                    n_sp_angles += 1
                elif is_improper:
                    n_improper_angles += 1
                else:
                    n_normal_angles += 1
    
    print(f"\n1-3 Angles: {n_angles} total")
    print(f"  SP (linear): {n_sp_angles}")
    print(f"  Improper: {n_improper_angles}")
    print(f"  Normal: {n_normal_angles}")
    
    # Count atom pairs marked
    torsion_pairs = set()
    for t in torsions:
        atoms = list(t['atomIndices'])
        i, l = min(atoms[0], atoms[3]), max(atoms[0], atoms[3])
        torsion_pairs.add((i, l))
    
    bond_pairs = set()
    for bond in mol.GetBonds():
        i, j = min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()), max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        bond_pairs.add((i, j))
    
    angle_pairs = set()
    for atom in mol.GetAtoms():
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        for a_idx in range(len(neighbors)):
            for b_idx in range(a_idx + 1, len(neighbors)):
                i, k = min(neighbors[a_idx], neighbors[b_idx]), max(neighbors[a_idx], neighbors[b_idx])
                angle_pairs.add((i, k))
    
    all_marked = torsion_pairs | bond_pairs | angle_pairs
    n_long_range = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) not in all_marked:
                n_long_range += 1
    
    print(f"\nAtom pair marking:")
    print(f"  Torsion 1-4 pairs: {len(torsion_pairs)}")
    print(f"  Bond 1-2 pairs: {len(bond_pairs)}")
    print(f"  Angle 1-3 pairs: {len(angle_pairs)}")
    print(f"  Total marked: {len(all_marked)}")
    print(f"  Long-range pairs: {n_long_range}")
    print(f"  Total possible: {n*(n-1)//2}")
    
    # Count improper/inversion atoms
    n_improper_atoms = 0
    for atom in mol.GetAtoms():
        elem = atom.GetAtomicNum()
        hyb = atom.GetHybridization()
        neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
        if elem in [6, 7, 8] and hyb == Chem.HybridizationType.SP2 and len(neighbors) == 3:
            n_improper_atoms += 1
    print(f"\nImproper (inversion) centers: {n_improper_atoms}")
    print(f"Inversion contributions: {n_improper_atoms * 3}")
    
    # Also count ring flatness torsions that would be added by basic knowledge
    ring_torsions = 0
    done_bonds = set()
    for ring in ring_info.AtomRings():
        rSize = len(ring)
        if rSize < 4 or rSize > 6:
            continue
        for i in range(rSize):
            a1 = ring[i]
            a2 = ring[(i+1) % rSize]
            a3 = ring[(i+2) % rSize]
            a4 = ring[(i+3) % rSize]
            bid = mol.GetBondBetweenAtoms(a2, a3).GetIdx()
            # Check doneBonds (from CSD torsions)
            bond_key = (min(a2, a3), max(a2, a3))
            # Check if CSD covers this bond
            csd_covers = False
            for t in torsions:
                atoms = list(t['atomIndices'])
                tc_bond = (min(atoms[1], atoms[2]), max(atoms[1], atoms[2]))
                if tc_bond == bond_key:
                    csd_covers = True
                    break
            
            if not csd_covers and bid not in done_bonds:
                # Check all 4 atoms SP2
                all_sp2 = all(mol.GetAtomWithIdx(a).GetHybridization() == Chem.HybridizationType.SP2 for a in [a1, a2, a3, a4])
                if all_sp2:
                    ring_torsions += 1
                    done_bonds.add(bid)
                    print(f"  Ring torsion: {a1}-{a2}-{a3}-{a4} (ring size {rSize})")
    
    print(f"\nRing flatness torsions (basic knowledge): {ring_torsions}")
    print(f"Total torsion contributions: {len(torsions) + ring_torsions}")


if __name__ == "__main__":
    main()
