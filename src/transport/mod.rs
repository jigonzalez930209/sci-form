//! Batch transport layer for conformer results.
//!
//! - [`arrow`] — Arrow-*like* columnar buffers for zero-copy interop. This
//!   is **not** the Apache Arrow IPC format: the crate emits no IPC stream
//!   and has no `arrow` dependency, only a columnar byte layout that a
//!   consumer can turn into a real Arrow batch cheaply.
//! - [`chunked`] — Splitting large SMILES lists into balanced chunks.
//! - [`worker`] — Task subdivision and worker count estimation for Web Workers.

pub mod arrow;
pub mod chunked;
pub mod worker;

pub use arrow::*;
pub use chunked::*;
pub use worker::*;
