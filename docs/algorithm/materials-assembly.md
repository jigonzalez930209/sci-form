# Materials Assembly

sci-form supports assembly of **Metal-Organic Frameworks (MOFs)** and other crystalline structures from **Secondary Building Units (SBUs)** and topological blueprints. This module handles periodic unit cells, coordination geometries, and framework construction.

> The reference, measured deviation, confidence region and limitations for this
> method are recorded in [Scientific validation](/validation/scientific-validation).

## Unit Cell and Lattice Vectors

<SvgDiagram src="/svg/materials-unit-cell.svg" alt="A triclinic cell with its lattice matrix, volume and the fractional-to-Cartesian mapping." />

### Lattice Parameters

A crystallographic unit cell is defined by six parameters:

| Parameter | Symbol | Description |
|-----------|--------|-------------|
| $a, b, c$ | Lengths | Edge lengths in Ångströms |
| $\alpha$ | Angle | Angle between $\vec{b}$ and $\vec{c}$ |
| $\beta$ | Angle | Angle between $\vec{a}$ and $\vec{c}$ |
| $\gamma$ | Angle | Angle between $\vec{a}$ and $\vec{b}$ |

### Lattice Matrix

The lattice vectors are stored as rows of a $3 \times 3$ matrix $M$:

$$
M = \begin{pmatrix} \vec{a} \\ \vec{b} \\ \vec{c} \end{pmatrix}
$$

For the general triclinic case:

$$
\vec{a} = (a, 0, 0)
$$
$$
\vec{b} = (b\cos\gamma, \; b\sin\gamma, \; 0)
$$
$$
\vec{c} = (c\cos\beta, \; c\frac{\cos\alpha - \cos\beta\cos\gamma}{\sin\gamma}, \; c \cdot v)
$$

where $v$ ensures the correct volume.

### Coordinate Transforms

**Fractional → Cartesian**:
$$
\vec{r}_{\text{cart}} = M^T \cdot \vec{f}
$$

**Cartesian → Fractional**:
$$
\vec{f} = (M^{-1})^T \cdot \vec{r}_{\text{cart}}
$$

### Volume

$$
V = \vec{a} \cdot (\vec{b} \times \vec{c}) = \det(M)
$$

### Periodic Boundary Conditions

**Wrapping** fractional coordinates into $[0, 1)$:

$$
\vec{f}_{\text{wrapped}} = \vec{f} - \lfloor \vec{f} \rfloor
$$

**Minimum image distance** under periodicity:

$$
\Delta \vec{f} = \Delta \vec{f} - \text{round}(\Delta \vec{f})
$$
$$
d = |M^T \cdot \Delta \vec{f}|
$$

### Supercell Generation

A $(n_a \times n_b \times n_c)$ supercell replicates atoms at:

$$
\vec{f}_{ijk} = \frac{\vec{f}_0 + (i, j, k)}{(n_a, n_b, n_c)}, \quad i \in [0, n_a), \; j \in [0, n_b), \; k \in [0, n_c)
$$

Lattice vectors scale: $\vec{a}' = n_a \vec{a}$, etc.

## Framework Assembly

<SvgDiagram src="/svg/materials-assembly.svg" alt="Assembling an sql framework: topology, metal nodes, Rodrigues-oriented linkers and the periodic crystal." />

### Assembly Pipeline

1. **Select topology**: Define node positions and edge connectivity (pcu, dia, sql, or an RCSR net from `topology_data`)
2. **Place metal nodes**: Position SBU at each topology vertex with appropriate coordination geometry
3. **Orient linkers**: Use Rodrigues rotation to align linker SBUs along each edge direction
4. **Select linker per edge type**: Inequivalent edges can use different linker SBUs
5. **Build crystal**: Combine all atoms with fractional coordinates into a `CrystalStructure`

### Topologies

Built-in homogeneous nets (`Topology::pcu/dia/sql/...`) and data-driven RCSR nets share the same assembler. `Topology::from_name` accepts both, and `Topology::from_net` converts a `topology_data::NetData` entry into the symmetry-expanded graph used by assembly.

| Topology | Coordination | Nodes | Edges | Description |
|----------|-------------|-------|-------|-------------|
| **pcu** | 6 | 1 | 3 | Primitive cubic — octahedral nodes |
| **dia** | 4 | 2 | 4 | Diamond — tetrahedral nodes |
| **sql** | 4 | 1 | 2 | Square lattice — square planar nodes |
| **bcu** | 8 | 2 | 8 | Body-centred cubic |
| **fcu** | 12 | 4 | 24 | Face-centred cubic |
| **nbo** | 4 | 2 | 4 | NbO net — square-planar nodes |
| **kgm** | 4 | 3 | 6 | Kagomé lattice |
| **hxl** | 6 | 1 | 3 | Hexagonal lattice |
| **pillared_sql** | 6 | 1 | 3 | Pillared square lattice |
| **acs** | 6 | 2 | 6 | RCSR net (`P6₃/mmc`) |
| **afw** | 4 | 9 | 18 | RCSR net with two edge types |
| **sod** | 4 | 12 | 24 | Sodalite net |
| **pts** | 4 | 4 | 8 | PtS net (data-driven expansion) |
| **xah** | 5 / 4 | 4 | 9 | RCSR net with two edge types |
| **tfy** | 4 / 3 | 184 | 288 | Large RCSR net with two edge types |

Nodes/edges are the symmetry-expanded counts in one unit cell.

### SBU Types

**Metal Node**: A central metal atom with connection points arranged in the specified coordination geometry.

```rust
let node = Sbu::metal_node(30, 2.0, CoordinationGeometry::Tetrahedral);
```

The first argument is the atomic number (`30` = Zn).

**Linear Linker**: A bridging unit connecting two metal nodes.

```rust
let linker = Sbu::linear_linker(&[6, 6, 6], 1.4, "BDC");
```

Elements are atomic numbers (`6` = C).

### Rodrigues Rotation

Linkers are oriented from the default $\hat{x}$ direction to the edge direction using the **Rodrigues rotation formula**:

$$
R(\vec{v}) = \vec{v} \cos\theta + (\hat{k} \times \vec{v}) \sin\theta + \hat{k}(\hat{k} \cdot \vec{v})(1 - \cos\theta)
$$

where $\hat{k}$ is the rotation axis (normalized cross product of $\hat{x}$ and the target direction) and $\theta$ is the angle between them.

## API

### Rust

```rust
use sci_form::materials::{CoordinationGeometry, Sbu, Topology, UnitCell};

let cell = UnitCell::cubic(10.0);
let node = Sbu::metal_node(30, 2.0, CoordinationGeometry::Octahedral);
let linker = Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
let topology = Topology::pcu();

let crystal = sci_form::assemble_framework(&node, &linker, &topology, &cell);
println!("{} atoms", crystal.num_atoms());
```

For nets with inequivalent edges, pass one linker per edge type:

```rust
let topology = Topology::from_name("xah")?;
let linkers = [
    Sbu::linear_linker(&[6], 1.4, "short"),
    Sbu::linear_linker(&[7], 1.2, "long"),
];
let crystal =
    sci_form::assemble_framework_with_linkers(&node, &linkers, &topology, &cell);
```

Edge `i` uses `linkers[topology.edge_types[i]]`. If the linker list is shorter than the type index, the last linker is reused. Topologies without `edge_types` (older serialized data) default to type `0`.

### CLI

```bash
# Unit cell info
sci-form cell --a 10 --b 10 --c 10 --alpha 90 --beta 90 --gamma 90

# Framework assembly
sci-form assemble --topology pcu --a 10.0 --metal 30 --geometry octahedral
```

### Python

```python
import sci_form

cell = sci_form.create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
print(cell.volume)

# topology, metal (atomic number), geometry, lattice_a, supercell
crystal = sci_form.assemble_framework("pcu", 30, "octahedral", 10.0, 1)
print(crystal.num_atoms)
```

Python and WASM currently expose the single-linker wrapper; the multi-linker assembler is available from Rust.

## Validation

- **Cubic cell**: Volume $= a^3$, all angles 90°
- **Coordinate roundtrip**: frac → cart → frac recovers original
- **PBC wrapping**: All fractional coords in $[0, 1)$ after wrapping
- **Minimum image**: Distance is always ≤ half the cell diagonal
- **Supercell scaling**: $(2 \times 2 \times 2)$ supercell has 8× the atoms
- **Extended nets (F67)**: node counts, coordination and edge-type counts match pormake's `.cgd` data
- **Multi-edge assembly**: unit tests in `materials::assembly` check that `xah` uses the linker selected by each `edge_types` entry
