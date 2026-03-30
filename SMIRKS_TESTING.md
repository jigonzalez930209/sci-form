# SMIRKS Reaction Testing Guide

This document describes how to test and validate the SMIRKS reaction transform implementation in sci-form, including comparisons with RDKit and OpenBabel.

## Overview

SMIRKS (SMILES Reaction Specification) is an extension of SMARTS that describes chemical reactions using atom-mapped reactant>>product patterns. The sci-form library implements SMIRKS parsing and pattern matching for chemical reaction transforms.

## Running Tests

### Rust Unit Tests

The core SMIRKS implementation has 23 comprehensive unit tests:

```bash
# Run all SMIRKS unit tests
cargo test --lib smirks

# Run with verbose output
cargo test --lib smirks -- --nocapture
```

### Rust Integration Tests

Integration tests cover common reaction patterns and edge cases:

```bash
# Run SMIRKS integration tests
cargo test --test test_smirks_reactions

# Run specific test
cargo test --test test_smirks_reactions test_acid_base_reactions
```

### Python Integration Tests

Python tests validate the Python bindings and can compare with RDKit:

```bash
# Build Python bindings first (requires maturin)
cd crates/python
maturin develop

# Run Python integration tests
cd ../..
python tests/integration/test_smirks_reactions.py
```

### RDKit/OpenBabel Comparison

Comprehensive comparison with RDKit and OpenBabel:

```bash
# Requires RDKit to be installed
# conda install -c conda-forge rdkit

python scripts/compare_smirks_reactions.py
```

This script tests 10 common reaction types across multiple molecules and generates a JSON report.

## Test Categories

### 1. Acid-Base Reactions

Tests deprotonation and protonation reactions:
- Carboxylic acid deprotonation: `[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]`
- Alcohol deprotonation: `[C:1][OH:2]>>[C:1][O-:2]`
- Amine protonation: `[N:1]>>[N:1+]`

### 2. Oxidation-Reduction

Tests oxidation and reduction transformations:
- Ketone reduction: `[C:1]=[O:2]>>[C:1][OH:2]`
- Aldehyde reduction: `[C:1][C:2]([H:3])=[O:4]>>[C:1][C:2]([H:3])[OH:4]`
- Alcohol oxidation: `[C:1][C:2]([OH:3])[H:4]>>[C:1][C:2](=[O:3])`

### 3. Substitution Reactions

Tests nucleophilic and electrophilic substitutions:
- Aromatic halogenation: `[c:1][H:2]>>[c:1][Cl:2]`
- Aromatic nitration: `[c:1][H:2]>>[c:1][N+:2](=[O:3])[O-:4]`

### 4. Hydrolysis

Tests bond cleavage reactions:
- Ester hydrolysis: `[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH:3]`
- Amide hydrolysis: `[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[OH:3]`

### 5. Edge Cases

Tests error handling and special cases:
- Invalid SMIRKS patterns
- Multi-component reactions (not yet supported)
- Atom map consistency
- Stereochemistry preservation
- Multiple reaction sites

## Test Results

### Current Status

**Rust Tests:**
- Unit tests: 23/23 passing ✓
- Integration tests: 15/15 passing ✓

**Python Tests:**
Requires Python bindings to be built with maturin.

### Comparison with RDKit

The `compare_smirks_reactions.py` script compares sci-form with RDKit on common reactions:

Expected Results:
- **Pattern parsing:** sci-form should parse all valid SMIRKS patterns
- **Pattern matching:** Agreement with RDKit on most organic reactions
- **Product generation:** Currently returns SMARTS patterns (full product generation TBD)

### Known Limitations

1. **Multi-component reactions:** Not yet supported for `apply_smirks`
   - Patterns with multiple reactants parse correctly
   - Application requires single-molecule input (largest fragment)

2. **Product generation:** Returns matched product SMARTS patterns
   - Full molecular product generation is future work
   - Atom mappings are correctly identified

3. **Stereochemistry:** Preserved in matching but not yet in products
   - Stereochemical centers detected correctly
   - Product stereochemistry generation is future work

## Writing New Tests

### Rust Unit Tests

Add tests to `src/smirks.rs`:

```rust
#[test]
fn test_my_reaction() {
    let result = apply_smirks(
        "[C:1]=[O:2]>>[C:1][OH:2]",
        "CC(=O)C"
    ).unwrap();
    assert!(result.success);
}
```

### Rust Integration Tests

Add tests to `tests/test_smirks_reactions.rs`:

```rust
#[test]
fn test_my_reaction_class() {
    // Test multiple related reactions
    let patterns = vec![
        "[C:1]=[O:2]>>[C:1][OH:2]",
        "[C:1]#[N:2]>>[C:1][NH2:2]",
    ];
    
    for pattern in patterns {
        let result = parse_smirks(pattern);
        assert!(result.is_ok());
    }
}
```

### Python Tests

Add tests to `tests/integration/test_smirks_reactions.py`:

```python
def test_my_reaction():
    import sci_form
    result = sci_form.apply_smirks(
        "[C:1]=[O:2]>>[C:1][OH:2]",
        "CC(=O)C"
    )
    assert result.success
```

### Comparison Tests

Add reactions to `scripts/compare_smirks_reactions.py`:

```python
REACTION_LIBRARY.append({
    "name": "My Reaction",
    "smirks": "[C:1]=[O:2]>>[C:1][OH:2]",
    "test_molecules": ["CC(=O)C", "c1ccc(C(=O)C)cc1"],
    "category": "reduction",
})
```

## Performance Benchmarks

Expected performance:
- **Parsing:** < 1ms for typical patterns
- **Matching:** < 10ms for small molecules (< 50 atoms)
- **Batch processing:** Patterns should be parsed once and reused

## Validation Checklist

Before committing changes to SMIRKS:

- [ ] All Rust unit tests pass (`cargo test --lib smirks`)
- [ ] All integration tests pass (`cargo test --test test_smirks_reactions`)
- [ ] Python bindings compile (`cargo build -p sci-form-python`)
- [ ] Python tests pass (if bindings built)
- [ ] Code passes clippy (`cargo clippy -- -D warnings`)
- [ ] Code is formatted (`cargo fmt --check`)
- [ ] Documentation is updated
- [ ] Comparison with RDKit shows expected agreement

## Continuous Integration

The CI pipeline should run:
1. Rust unit tests
2. Rust integration tests
3. Python binding compilation
4. Python integration tests (if RDKit available)
5. Comparison report generation (optional, for analysis)

## Troubleshooting

### "sci_form not installed" error

Build Python bindings:
```bash
cd crates/python
pip install maturin
maturin develop
```

### "RDKit not installed" error

Install RDKit:
```bash
conda install -c conda-forge rdkit
# or
pip install rdkit-pypi
```

### Test failures

1. Check that dependencies are installed
2. Verify Python version (3.9+)
3. Check Rust version (1.77+)
4. Run with `--nocapture` for detailed output
5. Check that the SMIRKS pattern is valid

## Future Improvements

1. **Full product generation:** Generate complete molecular products
2. **Multi-component support:** Handle reactions with multiple reactants
3. **Stereochemistry:** Full stereochemical product generation
4. **Reaction enumeration:** Generate all possible products
5. **Performance:** Optimize pattern matching for large molecules
6. **CLI tool:** Command-line interface for reaction testing

## References

- SMIRKS specification: [Daylight Theory Manual](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)
- RDKit reactions: [RDKit Documentation](https://www.rdkit.org/docs/RDKit_Book.html#chemical-reactions)
- OpenBabel: [OpenBabel Documentation](http://openbabel.org/docs/current/)
