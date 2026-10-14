//! Conformal Geometric Algebra (CGA) — core since 2026-10-14 (F81).
//!
//! Implements G(4,1) algebra with Motors for unified rotation-translation
//! operations on molecular coordinates and MOF assembly. Motors validated
//! vs Rodrigues + `clifford` to 3.4e-15; dihedral angle uses the mirror
//! convention (≡ −numpy mod 2π). Gaps: `refine_torsion_cga` vs matrix
//! method and materials assembly vs real MOF structures unmet.

mod conformer;
mod materials;
mod motor;
mod multivector;

pub use conformer::{
    apply_motor_to_subtree, dihedral_angle, dihedral_motor, embed_point, extract_point,
    refine_torsion_cga,
};
pub use materials::{assemble_framework_cga, place_sbu_cga, CgaFrame};
pub use motor::{conformal_point, extract_euclidean, Motor};
pub use multivector::Multivector;
