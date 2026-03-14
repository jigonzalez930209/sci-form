# Materials Assembly

sci-form supports assembly of **Metal-Organic Frameworks (MOFs)** and other crystalline structures from **Secondary Building Units (SBUs)** and topological blueprints. This module handles periodic unit cells, coordination geometries, and framework construction.

## Unit Cell and Lattice Vectors

<img src="/svg/materials-unit-cell.svg" alt="Crystallographic unit cell" class="svg-diagram" />

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

<img src="/svg/materials-assembly.svg" alt="Framework assembly pipeline" class="svg-diagram" />

### Assembly Pipeline

1. **Select topology**: Define node positions and edge connectivity (pcu, dia, sql, etc.)
2. **Place metal nodes**: Position SBU at each topology vertex with appropriate coordination geometry
3. **Orient linkers**: Use Rodrigues rotation to align linker SBUs along each edge direction
4. **Build crystal**: Combine all atoms with fractional coordinates into a `CrystalStructure`

### Topologies

| Topology | Coordination | Nodes | Edges | Description |
|----------|-------------|-------|-------|-------------|
| **pcu** | 6 | 1 | 3 | Primitive cubic — octahedral nodes |
| **dia** | 4 | 2 | 4 | Diamond — tetrahedral nodes |
| **sql** | 4 | 1 | 2 | Square lattice — square planar nodes |

### SBU Types

**Metal Node**: A central metal atom with connection points arranged in the specified coordination geometry.

```rust
let node = Sbu::metal_node("Zn", 2.0, CoordinationGeometry::Tetrahedral);
```

**Linear Linker**: A bridging unit connecting two metal nodes.

```rust
let linker = Sbu::linear_linker(&["C", "C", "C"], 1.4, "BDC");
```

### Rodrigues Rotation

Linkers are oriented from the default $\hat{x}$ direction to the edge direction using the **Rodrigues rotation formula**:

$$
R(\vec{v}) = \vec{v} \cos\theta + (\hat{k} \times \vec{v}) \sin\theta + \hat{k}(\hat{k} \cdot \vec{v})(1 - \cos\theta)
$$

where $\hat{k}$ is the rotation axis (normalized cross product of $\hat{x}$ and the target direction) and $\theta$ is the angle between them.

## API

### Rust

```rust
use sci_form::{create_unit_cell, assemble_framework};

let cell = create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
let crystal = assemble_framework("Zn", "BDC", "pcu", &cell);
```

### CLI

```bash
# Unit cell info
sci-form cell --a 10 --b 10 --c 10 --alpha 90 --beta 90 --gamma 90

# Framework assembly
sci-form assemble --topology pcu --a 10.0 --metal Zn --geometry octahedral
```

### Python

```python
import sci_form
cell = sci_form.unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
print(cell.volume())       # 1000.0
print(cell.parameters())   # (a, b, c, alpha, beta, gamma)

crystal = sci_form.assemble_framework("pcu", 10.0, "Zn", "octahedral")
print(crystal.num_atoms())
```

## Validation

- **Cubic cell**: Volume $= a^3$, all angles 90°
- **Coordinate roundtrip**: frac → cart → frac recovers original
- **PBC wrapping**: All fractional coords in $[0, 1)$ after wrapping
- **Minimum image**: Distance is always ≤ half the cell diagonal
- **Supercell scaling**: $(2 \times 2 \times 2)$ supercell has 8× the atoms
