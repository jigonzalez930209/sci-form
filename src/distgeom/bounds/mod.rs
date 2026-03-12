//! Bounds matrix construction and triangle smoothing for distance geometry.

mod geometry;
mod matrix;
mod smoothing;

pub use geometry::get_bond_length;
pub use matrix::{calculate_bounds_matrix, calculate_bounds_matrix_opts};
pub use smoothing::{smooth_bounds_matrix, triangle_smooth, triangle_smooth_tol};
