//! Batch transport layer for conformer results.
//!
//! - [`arrow`] — Apache Arrow columnar encoding for zero-copy interop.
//! - [`chunked`] — Splitting large SMILES lists into balanced chunks.
//! - [`worker`] — Task subdivision and worker count estimation for Web Workers.

pub mod arrow;
pub mod chunked;
pub mod worker;

pub use arrow::*;
pub use chunked::*;
pub use worker::*;
