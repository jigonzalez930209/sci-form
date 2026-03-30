# SMIRKS Reaction Examples

This document provides practical examples of using SMIRKS reaction transforms in sci-form.

## Quick Examples

### Rust

```rust
use sci_form::smirks::{parse_smirks, apply_smirks};

// Parse a SMIRKS pattern
let transform = parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]").unwrap();
println!("Reactants: {:?}", transform.reactant_smarts);
println!("Products: {:?}", transform.product_smarts);

// Apply to a molecule
let result = apply_smirks(
    "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
    "CC(=O)O"  // acetic acid
).unwrap();

if result.success {
    println!("Transform applied! {} atoms mapped", result.atom_mapping.len());
} else {
    println!("No match: {:?}", result.messages);
}
```

### Python

```python
import sci_form

# Parse SMIRKS
transform = sci_form.parse_smirks("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]")
print(f"Reactants: {transform.reactant_smarts}")
print(f"Products: {transform.product_smarts}")

# Apply to molecule
result = sci_form.apply_smirks(
    "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
    "CC(=O)O"  # acetic acid
)

if result.success:
    print(f"Success! {result.n_transforms} transforms applied")
    print(f"Products: {result.products}")
else:
    print(f"No match: {result.messages}")
```

## Common Reaction Patterns

### 1. Acid-Base Chemistry

#### Carboxylic Acid Deprotonation
```rust
let result = apply_smirks(
    "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",
    "CC(=O)O"
).unwrap();
```
Transforms acetic acid to acetate ion.

#### Alcohol Deprotonation
```rust
let result = apply_smirks(
    "[C:1][OH:2]>>[C:1][O-:2]",
    "CCO"
).unwrap();
```
Transforms ethanol to ethoxide.

#### Amine Protonation
```rust
let result = apply_smirks(
    "[N:1]>>[N+:1]",
    "CCN"
).unwrap();
```
Transforms amine to ammonium.

### 2. Oxidation-Reduction

#### Ketone Reduction
```rust
let result = apply_smirks(
    "[C:1]=[O:2]>>[C:1][OH:2]",
    "CC(=O)C"  // acetone
).unwrap();
```
Reduces ketone to secondary alcohol.

#### Alcohol Oxidation
```rust
let result = apply_smirks(
    "[C:1][C:2]([OH:3])[H:4]>>[C:1][C:2](=[O:3])",
    "CC(O)C"  // isopropanol
).unwrap();
```
Oxidizes secondary alcohol to ketone.

#### Aldehyde Reduction
```rust
let result = apply_smirks(
    "[C:1][C:2]([H:3])=[O:4]>>[C:1][C:2]([H:3])[OH:4]",
    "CC=O"  // acetaldehyde
).unwrap();
```
Reduces aldehyde to primary alcohol.

### 3. Substitution Reactions

#### Aromatic Halogenation
```rust
let result = apply_smirks(
    "[c:1][H:2]>>[c:1][Cl:2]",
    "c1ccccc1"  // benzene
).unwrap();
```
Substitutes aromatic hydrogen with chlorine.

#### Aromatic Nitration
```rust
let result = apply_smirks(
    "[c:1][H:2]>>[c:1][N+:2](=[O:3])[O-:4]",
    "c1ccccc1"  // benzene
).unwrap();
```
Adds nitro group to aromatic ring.

### 4. Hydrolysis

#### Ester Hydrolysis
```rust
let result = apply_smirks(
    "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH:3]",
    "CC(=O)OCC"  // ethyl acetate
).unwrap();
```
Cleaves ester to carboxylic acid.

#### Amide Hydrolysis
```rust
let result = apply_smirks(
    "[C:1](=[O:2])[N:3]>>[C:1](=[O:2])[OH:3]",
    "CC(=O)NC"  // N-methylacetamide
).unwrap();
```
Cleaves amide to carboxylic acid.

## Advanced Patterns

### Multi-Component Reactions

```rust
// Esterification (not yet supported for apply_smirks)
let transform = parse_smirks(
    "[C:1](=[O:2])[OH:3].[C:4][OH:5]>>[C:1](=[O:2])[O:5][C:4]"
).unwrap();

// Currently only parses, cannot apply to multi-component reactants
// This is a future enhancement
```

### Complex Atom Mapping

```rust
// SN2 substitution with detailed mapping
let transform = parse_smirks(
    "[C:1]([H:2])([H:3])([H:4])[Br:5]>>[C:1]([H:2])([H:3])([H:4])[OH:5]"
).unwrap();

// Maps all atoms including hydrogens
let result = apply_smirks(
    "[C:1][Br:2]>>[C:1][OH:2]",
    "CCBr"
).unwrap();
```

### Functional Group Selectivity

```rust
// Selective transformation of primary alcohols
let result = apply_smirks(
    "[C:1][CH2:2][OH:3]>>[C:1][CH2:2][O-:3]",
    "CCCO"  // propanol (primary alcohol)
).unwrap();

// Won't match secondary or tertiary alcohols
```

## Testing and Validation

### Testing Against Known Reactions

```rust
#[test]
fn test_williamson_ether_synthesis() {
    // Williamson ether synthesis pattern
    let smirks = "[C:1][OH:2].[C:3][Br:4]>>[C:1][O:2][C:3]";
    let transform = parse_smirks(smirks).unwrap();
    
    assert_eq!(transform.reactant_smarts.len(), 2);
    assert_eq!(transform.product_smarts.len(), 1);
}
```

### Comparing with RDKit

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import sci_form

smirks = "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]"
smiles = "CC(=O)O"

# sci-form
sf_result = sci_form.apply_smirks(smirks, smiles)

# RDKit
rxn = AllChem.ReactionFromSmarts(smirks)
mol = Chem.MolFromSmiles(smiles)
rdkit_products = rxn.RunReactants((mol,))

print(f"sci-form success: {sf_result.success}")
print(f"RDKit products: {len(rdkit_products)}")
```

## Performance Tips

### Reusing Parsed Transforms

```rust
// Parse once
let transform = parse_smirks("[C:1][OH:2]>>[C:1][O-:2]").unwrap();

// Apply to multiple molecules
let molecules = vec!["CCO", "CC(C)O", "c1ccccc1O"];

for smiles in molecules {
    // Note: Currently need to re-apply with SMIRKS string
    // Future: apply_transform(transform, smiles)
    let result = apply_smirks(
        "[C:1][OH:2]>>[C:1][O-:2]",
        smiles
    ).unwrap();
    
    if result.success {
        println!("{} -> transformed", smiles);
    }
}
```

### Batch Processing

```python
import sci_form

reactions = [
    ("[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]", "CC(=O)O"),
    ("[C:1][OH:2]>>[C:1][O-:2]", "CCO"),
    ("[N:1]>>[N+:1]", "CCN"),
]

results = []
for smirks, smiles in reactions:
    result = sci_form.apply_smirks(smirks, smiles)
    results.append({
        'smirks': smirks,
        'smiles': smiles,
        'success': result.success,
    })

# Process results
successful = [r for r in results if r['success']]
print(f"{len(successful)}/{len(results)} reactions successful")
```

## Common Pitfalls

### 1. Hydrogen Mapping

```rust
// Explicit hydrogen mapping
let smirks = "[C:1][H:2]>>[C:1][Cl:2]";  // May not match SMILES without explicit H

// Better: Use implicit hydrogen matching
let smirks = "[C:1]>>[C:1][Cl]";  // More flexible
```

### 2. Aromatic vs Aliphatic

```rust
// Aromatic carbon (lowercase 'c')
let aromatic = "[c:1][H:2]>>[c:1][Cl:2]";

// Aliphatic carbon (uppercase 'C')
let aliphatic = "[C:1][H:2]>>[C:1][Cl:2]";

// They match different molecule types!
```

### 3. Multi-Component Reactions

```rust
// Currently not supported for apply_smirks
let multi = "[A:1].[B:2]>>[A:1][B:2]";

// Workaround: Use single-component patterns only
let single = "[A:1]>>[A:1][B]";
```

## Future Enhancements

The following features are planned:

1. **Full Product Generation**: Generate complete molecular products, not just patterns
2. **Multi-Component Support**: Handle reactions with multiple reactants
3. **Stereochemistry**: Preserve and transform stereochemical information
4. **Reaction Enumeration**: Generate all possible products from a pattern
5. **Transform Reuse**: Apply parsed transforms without re-parsing

## References

- [SMIRKS Specification](https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html)
- [RDKit Chemical Reactions](https://www.rdkit.org/docs/RDKit_Book.html#chemical-reactions)
- [sci-form SMIRKS Testing Guide](SMIRKS_TESTING.md)
