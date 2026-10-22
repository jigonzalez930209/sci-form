//! Molecular force fields: UFF, MMFF94, ETKDG refinement, and minimization.
//!
//! # References
//!
//! - **UFF** — the Universal Force Field of Rappé, A. K. et al., "UFF, a full
//!   periodic table force field for molecular mechanics and molecular dynamics
//!   simulations", *J. Am. Chem. Soc.* **1992**, 114, 10024–10035. The
//!   parameter tables and the functional forms follow that paper and the
//!   canonical RDKit implementation (`ForceFieldHelpers/UFF`).
//! - **MMFF94** — the Merck Molecular Force Field of Halgren, "MMFF94: an
//!   extension of MMFF89", *J. Comput. Chem.* **1996**, 17, 490–519, with
//!   the published parameter data set.
//! - **ETKDG** — the distance-geometry embedding of Blainey, Todeschini &
//!   van Gunsteren, *Chem. Rev.* **2000**, 100, 2557, and the ETKDG refinement
//!   described by Weininger, *J. Chem. Inf. Comput. Sci.* **2008**, 48, 136.
//!
//! # Units
//!
//! Energies are in **kcal/mol**, gradients are `dE/dR` in **kcal/mol/Å**
//! (so the mechanical force is the negative of the returned gradient), bond
//! lengths and radii are in **Å**, and angles are in **radians** unless a
//! field name ends in `_deg` (degrees) or `_rad`.
//!
//! # Submodules
//!
//! - [`atom_typer`] — Element + topology → UFF atom type assignment.
//! - [`params`] — UFF parameter tables (bond radii, angles, torsions) and the
//!   Rappé formulas that turn them into force constants.
//! - [`uff`] / [`mmff94`] — Energy term implementations per force field.
//! - [`mmff94_params`] — Generated MMFF94 parameter tables (Halgren 1996).
//! - [`builder`] — Assembly of a UFF force field from a molecular graph.
//! - [`energy`] — Simplified UFF-style total-energy evaluation and the
//!   generic bond/angle/torsion/out-of-plane building blocks.
//! - [`gradients`] — Analytical gradient computation.
//! - [`minimizer`] — L-BFGS energy minimizer.
//! - [`traits`] — The `ForceFieldContribution` trait and the
//!   `MolecularForceField` container shared by every force field here.
//! - [`bounds_ff`] — Bounds-matrix force field for distance geometry.
//! - [`etkdg_3d`] — ETKDG 3D refinement with torsion terms.
//! - [`etkdg_lite`] — Lightweight ETKDG helpers (M6 empirical torsion).
//! - [`dg_terms`] — Distance-geometry specific energy contributions.
//! - [`torsion_scan`] — Systematic torsion angle scanning.
//! - [`reaxff`] — Reactive force field (bond orders, EEM, shielded nonbonded;
//!   validated vs LAMMPS `reaxff`, F78).
//!
//! # Caveats
//!
//! The **MMFF94** implementation covers bond stretching, angle bending,
//! torsion and van der Waals only. MMFF94's electrostatic, stretch–bend and
//! out-of-plane bending terms are *not* implemented, so the energies reported
//! here are the MMFF94 *valence + dispersion* part of the total, not the
//! complete MMFF94 energy.
//!
//! The **UFF** implementation uses a geometric-mean van der Waals combination
//! with no 1-4 scaling, and builds torsions only for non-terminal,
//! non-triple central bonds.

pub mod atom_typer;
pub mod bounds_ff;
pub mod builder;
pub mod dg_terms;
pub mod energy;
pub mod etkdg_3d;
pub mod etkdg_lite;
pub mod gradients;
pub mod minimizer;
pub mod mmff94;
pub mod mmff94_params;
pub mod params;
pub mod reaxff;
pub mod torsion_scan;
pub mod traits;
pub mod uff;

pub use bounds_ff::*;
pub use dg_terms::*;
pub use energy::*;
pub use etkdg_3d::*;
pub use gradients::*;
pub use minimizer::*;
pub use torsion_scan::*;
pub use traits::*;
