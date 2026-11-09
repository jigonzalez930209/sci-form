//! Backwards-compatible re-export of the core [`crate::render_bridge`] builders.
//!
//! Render-bridge domain builders were promoted to core on 2026-11-08 (F90).
//! Prefer `sci_form::render_bridge`; this path remains for existing
//! `render_bridge` imports.

pub use crate::render_bridge::*;
