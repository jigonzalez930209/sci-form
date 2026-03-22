//! Conformal Geometric Algebra (CGA) — Track E1
//!
//! Implements G(4,1) algebra with Motors for unified rotation-translation
//! operations on molecular coordinates and MOF assembly.

mod conformer;
mod materials;
mod motor;
mod multivector;

pub use conformer::{
    apply_motor_to_subtree, dihedral_motor, embed_point, extract_point, refine_torsion_cga,
};
pub use materials::{assemble_framework_cga, place_sbu_cga, CgaFrame};
pub use motor::Motor;
pub use multivector::Multivector;
