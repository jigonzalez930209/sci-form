# SMILES Parsing

The first step of the pipeline converts a SMILES (Simplified Molecular Input Line Entry System) string into a molecular graph suitable for 3D coordinate generation.

## From String to Graph

<SvgDiagram src="/svg/smiles-pipeline.svg" alt="smiles-pipeline" />

## The `Molecule` Structure

The parser produces a `Molecule` containing:

```rust
pub struct Molecule {
    pub graph: Graph<Atom, Bond, Undirected>,
    pub ring_info: Vec<Vec<usize>>,  // SSSR ring membership
}
```

### Atom Properties

Each atom stores the properties needed for 3D generation:

| Property | Type | Description |
|----------|------|-------------|
| `element` | `u8` | Atomic number (1=H, 6=C, 7=N, 8=O, ...) |
| `hybridization` | `Hybridization` | SP, SP2, SP3, or Unspecified |
| `formal_charge` | `i8` | Formal charge (−2 to +2 typically) |
| `aromatic` | `bool` | Whether the atom is aromatic |
| `num_implicit_h` | `u8` | Implicit hydrogens added |
| `chiral` | `Option<Chirality>` | `@` (CCW) or `@@` (CW) |
| `in_ring` | `bool` | Whether the atom is in any ring |

### Bond Properties

| Property | Type | Values |
|----------|------|--------|
| `order` | `BondOrder` | Single, Double, Triple, Aromatic |
| `stereo` | `BondStereo` | None, E (trans), Z (cis) |

## Hybridization Assignment

After building the graph, hybridization is assigned based on local environment:

$$
\text{hybridization} = \begin{cases}
\text{SP3} & \text{if all bonds are single} \\
\text{SP2} & \text{if any bond is double or atom is aromatic} \\
\text{SP} & \text{if any bond is triple or two double bonds}
\end{cases}
$$

Special cases:
- **Aromatic atoms** → SP2
- **N with 3 single bonds** → SP2 if in an aromatic ring, SP3 otherwise
- **O with 2 single bonds** → SP3 (e.g., ether)
- **Terminal atoms** (degree 1) → inherit from neighbor

## Implicit Hydrogen Addition

SMILES encodes most hydrogens implicitly. The parser adds them explicitly because 3D coordinate generation requires all atoms:

$$
n_H = v - \sum_{bonds} \text{order} - |q|
$$

where $v$ is the normal valence, and $q$ is the formal charge.

Standard valences used:

| Element | Valences |
|---------|----------|
| C | 4 |
| N | 3, 5 |
| O | 2 |
| S | 2, 4, 6 |
| P | 3, 5 |
| B | 3 |
| F, Cl, Br, I | 1 |
| Si | 4 |
| Se | 2 |

For aromatic atoms, the aromatic bond contributes 1.5 to the bond order sum, but in practice the parser assigns 1 per aromatic bond and adjusts the hydrogen count to maintain proper valence.

## Ring Detection (SSSR)

The **Smallest Set of Smallest Rings** (SSSR) is computed after graph construction. This is essential for:

- Bounds matrix: ring bonds constrain torsion angles
- SMARTS matching: ring-membership queries (`R`, `r`)
- ETKDG: ring torsion patterns differ from chain torsions
- Force field: ring planarity enforcement

The SSSR is found using a modified graph traversal:

1. Compute the cycle rank: $\mu = |E| - |V| + 1$
2. Find all shortest-path back edges
3. Extract $\mu$ independent rings

## Example: Phenol

For the SMILES `c1ccccc1O`:

<SvgDiagram src="/svg/smiles-phenol.svg" alt="smiles-phenol" />

Result: 13 atoms (6C + 1O + 6H), 13 bonds, 1 SSSR ring of size 6, all ring carbons SP2, oxygen SP3.
