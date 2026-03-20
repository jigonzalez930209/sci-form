# E1: Conformal Geometric Algebra (CGA)

**Module:** `sci_form::experimental::cga`
**Feature flag:** `experimental-cga`

---

## Overview

Replaces the internal geometry engine (separate quaternions + vectors) with Conformal Geometric Algebra $G(4,1)$, unifying rotations and translations into a single multiplicative operator called a **Motor** ($M X \tilde{M}$). This eliminates gimbal lock and enables massively parallelizable geometric operations for conformer generation and MOF assembly.

<SvgDiagram src="/svg/experimental-cga.svg" alt="CGA Pipeline" />

---

## Theory

### G(4,1) Multivector

The conformal model embeds 3D Euclidean space into a 5D space with signature $(4,1)$:

$$e_1^2 = e_2^2 = e_3^2 = e_+^2 = 1, \quad e_-^2 = -1$$

A general multivector has $2^5 = 32$ components spanning all grade combinations. The **geometric product** encodes both the inner and outer products:

$$ab = a \cdot b + a \wedge b$$

### Motors

A **Motor** is the composition of a rotor (rotation) and translator:

- **Rotor:** $R = \cos(\theta/2) - \sin(\theta/2)\hat{L}$
- **Translator:** $T = 1 - \frac{1}{2}t \cdot e_\infty$
- **Motor:** $M = TR$

Objects transform via the sandwich product $X' = M X \tilde{M}$.

### Point Embedding

3D coordinates lift to conformal points:

$$P = p + \frac{1}{2}|p|^2 e_\infty + e_o$$

Coordinates are recovered from the transformed multivector by extracting the Euclidean components and normalizing.

---

## API

```rust
use sci_form::experimental::cga::*;

// Multivector operations
let a = Multivector::basis(1);           // e₁
let b = Multivector::basis(2);           // e₂
let ab = a.geometric_product(&b);        // e₁e₂
let rev = ab.reverse();                  // ẽ₁₂ = e₂e₁

// Motor from axis-angle
let motor = Motor::from_axis_angle([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);

// Transform a point
let point = cga_point([1.0, 0.0, 0.0]);
let rotated = motor.apply(&point);
let coords = extract_point(&rotated);   // ≈ [0, 1, 0]

// Dihedral rotation
let dihedral_motor = Motor::dihedral([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.047);
```

---

## Applications

- **Conformer generation:** Dihedral torsion scans using CGA motors instead of rotation matrices
- **MOF assembly:** Unit cell lattice vectors as translation motors, SBU placement via motor composition
- **Supercell construction:** Motor exponentiation $M_{super} = M_{cell}^{n_a \cdot n_b \cdot n_c}$

---

## Tests

```bash
cargo test --features experimental-cga --test regression -- test_cga
```

Covers: multivector algebra, motor construction, point embedding round-trips, dihedral rotations, and materials assembly comparison.
