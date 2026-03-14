# SMARTS Pattern Matching

SMARTS (SMiles ARbitrary Target Specification) is a pattern language for describing molecular substructures. sci-form includes a complete SMARTS parser and matcher used to apply CSD torsion patterns to molecular torsion bonds.

## Role in the Pipeline

<img src="/svg/smarts-pattern.svg" alt="smarts-pattern" class="svg-diagram" />

For each rotatable bond in the molecule, the torsion matcher:

1. Iterates through 837 SMARTS patterns in priority order
2. Attempts to match each pattern against the molecular substructure around the bond
3. Returns the Fourier coefficients ($V_1..V_6$, sign) from the **first matching pattern**
4. If no CSD pattern matches, falls back to basic knowledge rules

## SMARTS Language

### Atom Primitives

| Syntax | Description | Example |
|--------|-------------|---------|
| `C`, `N`, `O` | Aliphatic element | `C` matches aliphatic carbon |
| `c`, `n`, `o` | Aromatic element | `c` matches aromatic carbon |
| `#6` | Atomic number | `#6` matches any carbon |
| `[*]` | Any atom | wildcard |
| `A` | Any aliphatic | |
| `a` | Any aromatic | |

### Atom Properties (inside brackets)

| Syntax | Description | Example |
|--------|-------------|---------|
| `H3` | Total hydrogen count | `[CH3]` = methyl carbon |
| `D4` | Total degree (heavy + H) | `[D4]` = 4 connections |
| `d3` | Heavy degree | `[d3]` = 3 heavy neighbors |
| `R` | In any ring | `[R]` |
| `r6` | In ring of size 6 | `[r6]` |
| `x2` | Ring bond count | `[x2]` = 2 ring bonds |
| `+`, `-` | Formal charge | `[N+]`, `[O-]` |
| `^2` | Hybridization (SP2) | `[C^2]` = SP2 carbon |

### Logical Operators

| Syntax | Meaning | Example |
|--------|---------|---------|
| `,` | OR | `[C,N]` = C or N |
| `&` | AND (high priority) | `[C&R]` = C in ring |
| `;` | AND (low priority) | `[C;H2]` = C with 2H |
| `!` | NOT | `[!C]` = not carbon |

### Bond Primitives

| Syntax | Description |
|--------|-------------|
| `-` | Single bond |
| `=` | Double bond |
| `#` | Triple bond |
| `:` | Aromatic bond |
| `~` | Any bond |
| `@` | Ring bond |
| `!@` | Not a ring bond |
| (implicit) | Single or aromatic |

### Recursive SMARTS

The `$(...)` syntax embeds a full SMARTS as an atom property:

```
[$(C(=O)N)]    // C that matches C(=O)N — i.e., amide carbon
```

This is used extensively in the CSD pattern library for complex environment matching.

## Matching Algorithm

### Substructure Search

The matcher uses **backtracking** with pruning:

```
function match(pattern, molecule):
    for each atom a in molecule:
        if atom_matches(pattern.root, a):
            mapping = {pattern.root → a}
            if extend_match(pattern, molecule, mapping):
                return mapping
    return None

function extend_match(pattern, molecule, mapping):
    next_pattern_atom = first unmapped pattern atom connected to mapped atoms
    if next_pattern_atom is None:
        return true  // complete match
    
    for each neighbor n of mapping[connected_pattern_atom]:
        if n not in mapping.values():
            if atom_matches(next_pattern_atom, n) and
               bond_matches(pattern_bond, molecule_bond):
                mapping[next_pattern_atom] = n
                if extend_match(pattern, molecule, mapping):
                    return true
                delete mapping[next_pattern_atom]
    return false
```

### Ring Info Optimization

Ring information is computed **once** per molecule and reused across all 837 pattern matches. This is critical for performance since ring queries (`R`, `r6`, `x2`) are common in CSD patterns.

The SSSR (Smallest Set of Smallest Rings) provides:
- `in_ring(atom)` — whether an atom participates in any ring
- `ring_sizes(atom)` — set of ring sizes containing this atom
- `ring_bond_count(atom)` — number of ring bonds

### Thread-Local Pattern Cache

Parsed SMARTS patterns are cached in a **thread-local** store:

```rust
thread_local! {
    static PARSED_PATTERNS: RefCell<Option<Vec<ParsedPattern>>> = RefCell::new(None);
}
```

On first use per thread, all 837 SMARTS strings are parsed into pattern graphs. Subsequent calls reuse the cached patterns. This amortizes the $O(837 \times P)$ parsing cost across molecules.

## Torsion Pattern Format

Each CSD torsion pattern in `torsion_data.rs` encodes:

```rust
TorsionPattern {
    smarts: "[C:1]-[C:2](-[OH])-[C:3]-[C:4]",
    atom_map: [1, 2, 3, 4],  // which atoms define the torsion
    v: [0.0, 3.5, 0.0, 1.2, 0.0, 0.3],  // V₁..V₆
    signs: [1, -1, 1, -1, 1, -1],        // sign for each term
}
```

The **atom map** (`:1`, `:2`, `:3`, `:4`) marks the four atoms that define the torsion angle. The matched atoms are ordered according to this map to correctly orient the torsion energy.

## Pattern Statistics

Of the 837 total patterns:

| Property | Value |
|----------|-------|
| Total patterns | 837 |
| v2 (general) | 365 |
| Macrocycle | 472 |
| Use recursive SMARTS | ~40% |
| Average pattern length | ~35 characters |
| Longest pattern | ~120 characters |
| Typical match time per bond | ~50 μs |

## Bond Filtering

Not all bonds are eligible for CSD torsion matching. The torsion matcher excludes:

- **Non-rotatable bonds** (double, triple, aromatic)
- **Bridged ring bonds** (bonds common to 3+ rings)
- **Fused ring bonds** (bonds in 4+ SSSR rings)
- **Terminal bonds** (hydrogen bonds, degree-1 atoms)

These bonds either have fixed geometry (double bonds) or are constrained by the ring system, so torsion preferences don't apply.
