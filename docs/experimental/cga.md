# E1: Conformal Geometric Algebra (CGA)

**Module:** `sci_form::cga` (core since 2026-10-14; `alpha::cga` re-export)
**Feature flag:** none (core); `alpha-cga` kept as no-op alias

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
use sci_form::cga::*;

// Motor from axis-angle
let motor = Motor::rotor([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);

// Transform a point
let rotated = motor.transform_point([1.0, 0.0, 0.0]);   // ≈ [0, 1, 0]

// Dihedral rotation about the a→b bond axis
let dm = dihedral_motor([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.047);
```

---

## Applications

- **Conformer generation:** Dihedral torsion scans using CGA motors instead of rotation matrices
- **MOF assembly:** Unit cell lattice vectors as translation motors, SBU placement via motor composition
- **Supercell construction:** Motor exponentiation $M_{super} = M_{cell}^{n_a \cdot n_b \cdot n_c}$

---

## Validation (F81, 2026-10-14)

External gate `tests/alpha_cga_validation.rs` (no feature needed) vs
`scripts/validation/oracle_cga.py` (Rodrigues closed form + `clifford`
Cl(4,1) sandwich, seed 42): motors 3.4e-15 on 16 rigid motions (validates
the bitmask blade encoding, incl. the e31 sign), dihedral rotation 6.6e-15,
subtree isolation bit-exact, dihedral angle exact up to the documented
mirror convention (ours ≡ −numpy mod 2π). Gaps: `refine_torsion_cga` vs
matrix method and materials assembly vs real MOF structures unmet.

## Tests

```bash
cargo test --features "beta-cpm alpha-cga" --test experimental_gates -- test_cga
```

20 integration tests (`tests/experimental/test_cga.rs`) covering: multivector algebra, motor construction, point embedding round-trips, dihedral rotations, and materials assembly comparison.

These tests live in the `experimental_gates` target — see [Overview](./overview.md#run-tests) for the full command matrix.
