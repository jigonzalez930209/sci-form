//! Conformal Geometric Algebra (CGA) — Track E1
//!
//! Implements G(4,1) algebra with Motors for unified rotation-translation
//! operations on molecular coordinates and MOF assembly.

mod multivector;
mod motor;
mod conformer;
mod materials;

pub use multivector::Multivector;
pub use motor::Motor;
pub use conformer::{
    embed_point, extract_point, dihedral_motor, apply_motor_to_subtree,
    refine_torsion_cga,
};
pub use materials::{
    CgaFrame, place_sbu_cga, assemble_framework_cga,
};
