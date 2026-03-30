#!/usr/bin/env python3
"""
Comprehensive SMIRKS reaction comparison: sci-form vs RDKit vs OpenBabel

This script tests a comprehensive set of chemical reactions using SMIRKS patterns
and compares the results from sci-form with those from RDKit and OpenBabel.

The goal is to:
1. Identify any differences in SMIRKS interpretation
2. Document accuracy and compatibility
3. Provide a benchmark for reaction transform quality
"""

import json
import sys
import os
from typing import List, Dict, Optional, Tuple

# Common organic reaction SMIRKS patterns
REACTION_LIBRARY = [
    {
        "name": "Carboxylic acid deprotonation",
        "smirks": "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
        "test_molecules": ["CC(=O)O", "C(=O)O", "CC(C)C(=O)O"],
        "category": "acid-base",
    },
    {
        "name": "Alcohol deprotonation",
        "smirks": "[C:1][OH:2]>>[C:1][O-:2]",
        "test_molecules": ["CCO", "CC(C)O", "c1ccccc1O"],
        "category": "acid-base",
    },
    {
        "name": "Amine protonation",
        "smirks": "[N:1]>>[N:1+]",
        "test_molecules": ["CCN", "c1cccnc1", "NCCN"],
        "category": "acid-base",
    },
    {
        "name": "Ketone to alcohol reduction",
        "smirks": "[C:1]=[O:2]>>[C:1][OH:2]",
        "test_molecules": ["CC(=O)C", "c1ccc(C(=O)C)cc1"],
        "category": "reduction",
    },
    {
        "name": "Aldehyde to alcohol reduction",
        "smirks": "[C:1][C:2]([H:3])=[O:4]>>[C:1][C:2]([H:3])[OH:4]",
        "test_molecules": ["CC=O", "c1ccccc1C=O"],
        "category": "reduction",
    },
    {
        "name": "Alcohol to ketone oxidation",
        "smirks": "[C:1][C:2]([OH:3])[H:4]>>[C:1][C:2](=[O:3])",
        "test_molecules": ["CC(O)C", "c1ccccc1C(O)C"],
        "category": "oxidation",
    },
    {
        "name": "Aromatic halogenation",
        "smirks": "[c:1][H:2]>>[c:1][Cl:2]",
        "test_molecules": ["c1ccccc1", "c1ccncc1"],
        "category": "substitution",
    },
    {
        "name": "Aromatic nitration",
        "smirks": "[c:1][H:2]>>[c:1][N+:2](=[O:3])[O-:4]",
        "test_molecules": ["c1ccccc1", "c1ccc(C)cc1"],
        "category": "substitution",
    },
    {
        "name": "Ester hydrolysis",
        "smirks": "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH:3]",
        "test_molecules": ["CC(=O)OCC", "c1ccccc1C(=O)OC"],
        "category": "hydrolysis",
    },
    {
        "name": "Amide hydrolysis",
        "smirks": "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[OH:3]",
        "test_molecules": ["CC(=O)NC", "c1ccccc1C(=O)NC"],
        "category": "hydrolysis",
    },
]


def test_sciform_reaction(smirks: str, smiles: str) -> Dict:
    """Test a SMIRKS reaction using sci-form."""
    try:
        import sci_form
        result = sci_form.apply_smirks(smirks, smiles)
        return {
            "success": result.success,
            "n_transforms": result.n_transforms,
            "products": result.products if result.success else [],
            "error": None if result.success else result.messages,
        }
    except ImportError:
        return {
            "success": None,
            "error": "sci_form not installed",
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }


def test_rdkit_reaction(smirks: str, smiles: str) -> Dict:
    """Test a SMIRKS reaction using RDKit."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        
        rxn = AllChem.ReactionFromSmarts(smirks)
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            return {"success": False, "error": "Invalid SMILES"}
        
        products = rxn.RunReactants((mol,))
        
        if len(products) > 0:
            # Convert products to SMILES
            product_smiles = []
            for product_set in products:
                for p in product_set:
                    try:
                        smi = Chem.MolToSmiles(p)
                        product_smiles.append(smi)
                    except:
                        pass
            
            return {
                "success": True,
                "n_products": len(products),
                "products": product_smiles,
                "error": None,
            }
        else:
            return {
                "success": False,
                "n_products": 0,
                "error": "No products",
            }
    except ImportError:
        return {
            "success": None,
            "error": "RDKit not installed",
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }


def test_openbabel_reaction(smirks: str, smiles: str) -> Dict:
    """Test reaction using OpenBabel (limited SMIRKS support)."""
    try:
        from openbabel import openbabel as ob
        
        # OpenBabel doesn't have native SMIRKS support like RDKit
        # We can still parse the molecule though
        conv = ob.OBConversion()
        conv.SetInAndOutFormats("smi", "smi")
        
        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return {"success": False, "error": "Invalid SMILES"}
        
        return {
            "success": None,
            "error": "OpenBabel doesn't support SMIRKS transforms directly",
        }
    except ImportError:
        return {
            "success": None,
            "error": "OpenBabel not installed",
        }
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
        }


def compare_reactions():
    """Compare all reactions across implementations."""
    print("=" * 80)
    print("SMIRKS Reaction Comparison: sci-form vs RDKit vs OpenBabel")
    print("=" * 80)
    
    results = []
    
    for reaction in REACTION_LIBRARY:
        print(f"\n{reaction['name']} ({reaction['category']})")
        print(f"  SMIRKS: {reaction['smirks']}")
        print("-" * 80)
        
        for smiles in reaction['test_molecules']:
            print(f"  Test: {smiles}")
            
            # Test with each implementation
            sf_result = test_sciform_reaction(reaction['smirks'], smiles)
            rdkit_result = test_rdkit_reaction(reaction['smirks'], smiles)
            ob_result = test_openbabel_reaction(reaction['smirks'], smiles)
            
            # Display results
            if sf_result['success'] is not None:
                print(f"    sci-form: {'✓' if sf_result['success'] else '✗'} "
                      f"(transforms: {sf_result.get('n_transforms', 0)})")
            else:
                print(f"    sci-form: SKIP ({sf_result.get('error', 'unknown')})")
            
            if rdkit_result['success'] is not None:
                print(f"    RDKit:    {'✓' if rdkit_result['success'] else '✗'} "
                      f"(products: {rdkit_result.get('n_products', 0)})")
            else:
                print(f"    RDKit:    SKIP ({rdkit_result.get('error', 'unknown')})")
            
            # Compare agreement
            if (sf_result['success'] is not None and 
                rdkit_result['success'] is not None):
                if sf_result['success'] == rdkit_result['success']:
                    print("    Agreement: ✓ Both agree")
                else:
                    print("    Agreement: ✗ DIFFERENT RESULTS!")
            
            results.append({
                'reaction': reaction['name'],
                'smirks': reaction['smirks'],
                'smiles': smiles,
                'sci_form': sf_result,
                'rdkit': rdkit_result,
                'openbabel': ob_result,
            })
    
    return results


def generate_report(results: List[Dict]):
    """Generate a summary report."""
    print("\n" + "=" * 80)
    print("SUMMARY REPORT")
    print("=" * 80)
    
    sf_available = any(r['sci_form']['success'] is not None for r in results)
    rdkit_available = any(r['rdkit']['success'] is not None for r in results)
    
    if not sf_available:
        print("⚠ sci-form not available - Python bindings need to be built")
        return
    
    if not rdkit_available:
        print("⚠ RDKit not available - comparison limited")
        print("\nTo install RDKit: conda install -c conda-forge rdkit")
        return
    
    # Count agreements and disagreements
    comparable = [r for r in results 
                  if r['sci_form']['success'] is not None 
                  and r['rdkit']['success'] is not None]
    
    if not comparable:
        print("No comparable results found")
        return
    
    agreements = sum(1 for r in comparable 
                     if r['sci_form']['success'] == r['rdkit']['success'])
    disagreements = len(comparable) - agreements
    
    print(f"\nTotal test cases: {len(comparable)}")
    print(f"Agreements: {agreements} ({100*agreements/len(comparable):.1f}%)")
    print(f"Disagreements: {disagreements} ({100*disagreements/len(comparable):.1f}%)")
    
    if disagreements > 0:
        print("\nDisagreements:")
        for r in comparable:
            if r['sci_form']['success'] != r['rdkit']['success']:
                print(f"  - {r['reaction']}: {r['smiles']}")
                print(f"    sci-form: {r['sci_form']['success']}")
                print(f"    RDKit: {r['rdkit']['success']}")
    
    # Save detailed results
    output_file = "smirks_comparison_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nDetailed results saved to: {output_file}")


def main():
    """Main entry point."""
    results = compare_reactions()
    generate_report(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
