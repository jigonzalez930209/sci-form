#!/usr/bin/env python3
"""Cross-language integration test for sci-form.

Tests that Python bindings produce correct results and match
the expected behavior across all API functions.
"""
import sys
import json
import time


def test_version():
    import sci_form
    v = sci_form.version()
    assert v.startswith("sci-form"), f"Bad version: {v}"
    print(f"  version: {v}")


def test_embed_single():
    import sci_form

    # Ethanol
    r = sci_form.embed("CCO", seed=42)
    assert r.is_ok(), f"Embed failed: {r.error}"
    assert r.num_atoms == 9, f"Expected 9 atoms, got {r.num_atoms}"
    assert len(r.coords) == 27, f"Expected 27 coords, got {len(r.coords)}"
    assert len(r.elements) == 9
    assert r.elements[0] == 6  # Carbon
    assert r.elements[2] == 8  # Oxygen
    assert len(r.bonds) > 0

    # Check positions helper
    pos = r.get_positions()
    assert len(pos) == 9
    assert len(pos[0]) == 3

    # Reproducibility: same seed = same result
    r2 = sci_form.embed("CCO", seed=42)
    assert r.coords == r2.coords, "Same seed should give same result"

    print(f"  embed CCO: {r.num_atoms} atoms, {len(r.bonds)} bonds, {r.time_ms:.1f}ms")


def test_embed_molecules():
    import sci_form

    molecules = {
        "c1ccccc1": 12,       # benzene
        "CC(=O)O": 8,          # acetic acid
        "CC(=O)Oc1ccccc1C(=O)O": 21,  # aspirin
        "C": 5,                 # methane
        "N": 4,                 # ammonia
        "O": 3,                 # water
        "C#N": 3,               # hydrogen cyanide
        "C=C": 6,               # ethylene
        "c1ccncc1": 11,         # pyridine
    }

    for smi, expected_atoms in molecules.items():
        r = sci_form.embed(smi, seed=42)
        assert r.is_ok(), f"Failed for {smi}: {r.error}"
        assert r.num_atoms == expected_atoms, f"{smi}: expected {expected_atoms} atoms, got {r.num_atoms}"

    print(f"  embed {len(molecules)} molecules: all OK")


def test_batch():
    import sci_form

    smiles = ["CCO", "c1ccccc1", "CC(=O)O", "C#N", "c1ccncc1"]
    results = sci_form.embed_batch(smiles, seed=42, num_threads=2)

    assert len(results) == len(smiles), f"Expected {len(smiles)} results, got {len(results)}"
    for r in results:
        assert r.is_ok(), f"Batch failed for {r.smiles}: {r.error}"

    # Verify batch results match single results
    for smi in smiles:
        single = sci_form.embed(smi, seed=42)
        batch_r = next(r for r in results if r.smiles == smi)
        assert single.coords == batch_r.coords, f"Batch/single mismatch for {smi}"

    print(f"  batch {len(smiles)} molecules: all match single results")


def test_parse():
    import sci_form

    info = sci_form.parse("CCO")
    assert info["num_atoms"] == 9
    assert info["num_bonds"] == 8
    assert len(info["atoms"]) == 9

    # Check atom info
    carbons = [a for a in info["atoms"] if a["element"] == 6]
    oxygens = [a for a in info["atoms"] if a["element"] == 8]
    hydrogens = [a for a in info["atoms"] if a["element"] == 1]
    assert len(carbons) == 2
    assert len(oxygens) == 1
    assert len(hydrogens) == 6

    print(f"  parse CCO: {info['num_atoms']} atoms, {info['num_bonds']} bonds")


def test_error_handling():
    import sci_form

    # Invalid SMILES
    r = sci_form.embed("INVALID_SMILES_XYZ", seed=42)
    assert not r.is_ok(), "Should have failed for invalid SMILES"
    assert r.error is not None

    # Parse error
    try:
        sci_form.parse("INVALID_SMILES_XYZ")
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

    print("  error handling: OK")


def test_coordinate_sanity():
    import sci_form
    import math

    r = sci_form.embed("CCO", seed=42)
    pos = r.get_positions()

    # Bond lengths should be reasonable (0.8-2.0 Å for typical bonds)
    for i, (_, _, order) in enumerate(r.bonds):
        a, b = r.bonds[i][0], r.bonds[i][1]
        dx = pos[a][0] - pos[b][0]
        dy = pos[a][1] - pos[b][1]
        dz = pos[a][2] - pos[b][2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        assert 0.5 < dist < 3.0, f"Bond {a}-{b} length {dist:.2f} out of range"

    print(f"  coordinate sanity: all bond lengths in [0.5, 3.0] Å")


def test_performance():
    import sci_form

    # Single molecule should be < 50ms
    smiles = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O", "c1ccncc1"]
    for smi in smiles:
        r = sci_form.embed(smi, seed=42)
        assert r.time_ms < 50, f"{smi} took {r.time_ms:.1f}ms (> 50ms target)"

    # Batch performance
    batch_smiles = smiles * 20  # 100 molecules
    start = time.time()
    results = sci_form.embed_batch(batch_smiles, seed=42)
    elapsed_ms = (time.time() - start) * 1000
    avg_ms = elapsed_ms / len(batch_smiles)

    ok = sum(1 for r in results if r.is_ok())
    print(f"  performance: {len(batch_smiles)} mols in {elapsed_ms:.0f}ms ({avg_ms:.1f}ms/mol), {ok}/{len(batch_smiles)} OK")


def test_data_roundtrip():
    """Test that data can be serialized/deserialized for interop."""
    import sci_form

    r = sci_form.embed("CCO", seed=42)

    # Simulate sending to another system via JSON
    data = {
        "smiles": r.smiles,
        "num_atoms": r.num_atoms,
        "coords": r.coords,
        "elements": list(r.elements),
        "bonds": [(a, b, o) for a, b, o in r.bonds],
        "error": r.error,
        "time_ms": r.time_ms,
    }
    json_str = json.dumps(data)
    roundtripped = json.loads(json_str)

    assert roundtripped["num_atoms"] == r.num_atoms
    assert roundtripped["coords"] == list(r.coords)
    assert roundtripped["elements"] == list(r.elements)

    print("  data roundtrip: JSON serialize/deserialize OK")


def main():
    tests = [
        ("Version", test_version),
        ("Single embed", test_embed_single),
        ("Multiple molecules", test_embed_molecules),
        ("Batch parallel", test_batch),
        ("Parse", test_parse),
        ("Error handling", test_error_handling),
        ("Coordinate sanity", test_coordinate_sanity),
        ("Performance", test_performance),
        ("Data roundtrip", test_data_roundtrip),
    ]

    print("=" * 60)
    print("sci-form Python Integration Tests")
    print("=" * 60)

    passed = 0
    failed = 0
    for name, test_fn in tests:
        try:
            test_fn()
            passed += 1
            print(f"  [PASS] {name}")
        except Exception as e:
            failed += 1
            print(f"  [FAIL] {name}: {e}")

    print("=" * 60)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 60)
    sys.exit(1 if failed > 0 else 0)


if __name__ == "__main__":
    main()
