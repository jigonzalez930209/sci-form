#!/usr/bin/env python3
"""
Test SMIRKS/reaction functionality and compare with RDKit.

This test suite validates the sci-form SMIRKS implementation by:
1. Testing basic SMIRKS parsing and application
2. Comparing reaction transforms with RDKit's AllChem.ReactionFromSmarts
3. Documenting any differences in behavior
"""

import sys
import os

# Add project root to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))


def test_smirks_parse():
    """Test basic SMIRKS parsing."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed (Python bindings not built)")
        return False
    
    # Basic acid deprotonation
    transform = sci_form.parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]")
    assert transform is not None
    assert len(transform.reactant_smarts) == 1
    assert len(transform.product_smarts) == 1
    print("  ✓ Basic SMIRKS parsing works")
    return True


def test_smirks_apply():
    """Test applying SMIRKS to molecules."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed")
        return False
    
    # Apply acid deprotonation to acetic acid
    result = sci_form.apply_smirks(
        "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
        "CC(=O)O"
    )
    
    assert result is not None
    assert result.success
    assert result.n_transforms == 1
    print("  ✓ SMIRKS application works")
    return True


def test_smirks_no_match():
    """Test SMIRKS that doesn't match."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed")
        return False
    
    # Try to apply nitrogen pattern to ethanol (no N present)
    result = sci_form.apply_smirks("[N:1]>>[N:1]", "CCO")
    
    assert result is not None
    assert not result.success
    assert result.n_transforms == 0
    print("  ✓ Non-matching SMIRKS handled correctly")
    return True


def test_common_reactions():
    """Test common organic reactions."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed")
        return False
    
    reactions = [
        # (name, smirks, reactant_smiles)
        ("Deprotonation", "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]", "CC(=O)O"),
        ("Alcohol deprotonation", "[C:1][OH:2]>>[C:1][O-:2]", "CCO"),
        ("Aromatic halogenation", "[c:1][H:2]>>[c:1][Cl:2]", "c1ccccc1"),
    ]
    
    for name, smirks, smiles in reactions:
        result = sci_form.apply_smirks(smirks, smiles)
        if result.success:
            print(f"  ✓ {name}: matched successfully")
        else:
            print(f"  ⚠ {name}: {result.messages}")
    
    return True


def compare_with_rdkit():
    """Compare sci-form SMIRKS with RDKit reactions."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed")
        return False
    
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        print("SKIP: RDKit not installed (comparison not possible)")
        return True  # Not a failure, just can't compare
    
    print("\n  Comparing with RDKit:")
    
    test_cases = [
        {
            "name": "Carboxylic acid deprotonation",
            "smirks": "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
            "smiles": "CC(=O)O",
        },
        {
            "name": "Alcohol deprotonation",
            "smirks": "[C:1][OH:2]>>[C:1][O-:2]",
            "smiles": "CCO",
        },
    ]
    
    for test in test_cases:
        # sci-form result
        sf_result = sci_form.apply_smirks(test["smirks"], test["smiles"])
        
        # RDKit result
        try:
            rxn = AllChem.ReactionFromSmarts(test["smirks"])
            mol = Chem.MolFromSmiles(test["smiles"])
            products = rxn.RunReactants((mol,))
            rdkit_matched = len(products) > 0
        except Exception as e:
            print(f"    ⚠ RDKit error for {test['name']}: {e}")
            rdkit_matched = False
        
        # Compare
        if sf_result.success == rdkit_matched:
            print(f"    ✓ {test['name']}: both agree (match={sf_result.success})")
        else:
            print(f"    ⚠ {test['name']}: sci-form={sf_result.success}, RDKit={rdkit_matched}")
    
    return True


def compare_with_openbabel():
    """Compare sci-form SMIRKS with OpenBabel if available."""
    try:
        import sci_form
    except ImportError:
        print("SKIP: sci_form not installed")
        return False
    
    try:
        from openbabel import openbabel as ob
    except ImportError:
        print("SKIP: OpenBabel not installed (comparison not possible)")
        return True  # Not a failure, just can't compare
    
    print("\n  OpenBabel comparison:")
    print("    Note: OpenBabel doesn't have direct SMIRKS support like RDKit")
    print("    (OpenBabel uses reaction SMARTS differently)")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("SMIRKS/Reaction Tests")
    print("=" * 60)
    
    tests = [
        ("Parse SMIRKS", test_smirks_parse),
        ("Apply SMIRKS", test_smirks_apply),
        ("Non-matching SMIRKS", test_smirks_no_match),
        ("Common reactions", test_common_reactions),
        ("RDKit comparison", compare_with_rdkit),
        ("OpenBabel comparison", compare_with_openbabel),
    ]
    
    passed = 0
    failed = 0
    skipped = 0
    
    for name, test_fn in tests:
        print(f"\n[TEST] {name}")
        try:
            result = test_fn()
            if result:
                passed += 1
            elif result is False:
                failed += 1
            else:
                skipped += 1
        except Exception as e:
            failed += 1
            print(f"  ✗ FAILED: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "=" * 60)
    print(f"Results: {passed} passed, {failed} failed, {skipped} skipped")
    print("=" * 60)
    
    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
