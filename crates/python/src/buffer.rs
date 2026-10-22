//! A minimal, safe wrapper around CPython's buffer protocol.
//!
//! PyO3 ships a `pyo3::buffer::PyBuffer<T>` type, but it is only compiled
//! when the crate is *not* built against the limited API.  The Python
//! bindings are built as `abi3-py38` wheels, so `pyo3::buffer` and the
//! `pyo3::ffi` buffer declarations are both compiled out.
//!
//! This module therefore declares the two CPython entry points it needs
//! itself.  `PyObject_GetBuffer` and `PyBuffer_Release` have been exported
//! by every CPython build since 2.6, so the symbols resolve at load time
//! exactly like the rest of the CPython API that PyO3 already calls.
//!
//! The wrapper is deliberately small:
//!
//! * it requests a **writable** buffer and refuses read-only ones,
//! * it refuses any element type other than `float64`,
//! * it refuses non-contiguous or multi-dimensional buffers,
//! * it releases the buffer in `Drop`, so a Python exception can never leak
//!   the view.
//!
//! # Core change needed
//!
//! Building the crate with `abi3-py311` (or without `abi3` at all) makes
//! `pyo3::buffer::PyBuffer<f64>` available and lets this module be replaced
//! by the upstream type.  That change belongs to `crates/python/Cargo.toml`,
//! which is owned by the packaging work stream.

use std::ffi::{c_char, c_int, c_void};

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// CPython's `Py_buffer` view structure.
///
/// Layout-compatible with `CPython/Include/object.h`; the fields are the
/// ones the buffer protocol documents.
#[repr(C)]
struct CPyBuffer {
    buf: *mut c_void,
    obj: *mut c_void,
    len: isize,
    itemsize: isize,
    readonly: c_int,
    ndim: c_int,
    format: *mut c_char,
    shape: *mut isize,
    strides: *mut isize,
    suboffsets: *mut isize,
    internal: *mut c_void,
}

/// `PyBUF_WRITABLE`: the exporter must provide a writable buffer.
const PYBUF_WRITABLE: c_int = 0x0001;
/// `PyBUF_FORMAT`: the view must carry the element format string.
const PYBUF_FORMAT: c_int = 0x0004;

extern "C" {
    /// Request a buffer view of `obj`.
    ///
    /// Returns 0 on success and -1 with a Python exception set on failure.
    fn PyObject_GetBuffer(obj: *mut pyo3::ffi::PyObject, view: *mut CPyBuffer, flags: c_int)
        -> c_int;

    /// Release a buffer view previously obtained with `PyObject_GetBuffer`.
    fn PyBuffer_Release(view: *mut CPyBuffer);
}

/// A borrowed, writable, contiguous `float64` buffer.
pub(crate) struct F64Buffer {
    view: Option<CPyBuffer>,
    len: usize,
}

impl F64Buffer {
    /// Borrow `obj` as a writable, contiguous `float64` buffer.
    pub(crate) fn writable_f64(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        // SAFETY: `PyObject_GetBuffer` fills every field of the view on
        // success.  `MaybeUninit` avoids initialising a struct that holds
        // raw pointers.
        let mut view = unsafe {
            let mut uninit = std::mem::MaybeUninit::<CPyBuffer>::uninit();
            let ptr = uninit.as_mut_ptr();
            let flags = PYBUF_WRITABLE | PYBUF_FORMAT;
            if PyObject_GetBuffer(obj.as_ptr(), ptr, flags) != 0 {
                return Err(PyErr::fetch(obj.py()));
            }
            uninit.assume_init()
        };
        let itemsize = view.itemsize as usize;
        let len = (view.len as usize).checked_div(itemsize).unwrap_or(0);
        let format = if view.format.is_null() {
            String::new()
        } else {
            // SAFETY: `view.format` is a NUL-terminated string owned by the
            // exporter, valid until the view is released.
            unsafe { std::ffi::CStr::from_ptr(view.format) }
                .to_string_lossy()
                .trim_end_matches('\0')
                .to_string()
        };
        // Without `PyBUF_STRIDES` the exporter must hand out a C-contiguous
        // buffer, and without `PyBUF_ND` the view is one-dimensional.
        let layout_ok = view.ndim <= 1 && itemsize == 8 && format == "d";
        if !layout_ok {
            // SAFETY: the view was filled by `PyObject_GetBuffer` and has not
            // been released yet.
            unsafe {
                PyBuffer_Release(&mut view);
            }
            return Err(PyValueError::new_err(format!(
                "Invalid buffer: expected a writable, contiguous float64 (64-bit) array, \
                 got format '{format}' with {len} element(s)."
            )));
        }
        Ok(Self {
            view: Some(view),
            len,
        })
    }

    /// Number of `f64` elements in the buffer.
    pub(crate) fn len(&self) -> usize {
        self.len
    }

    /// The buffer contents as a mutable slice.
    pub(crate) fn as_mut_slice(&mut self) -> &mut [f64] {
        let ptr = self
            .view
            .as_ref()
            .expect("the view is released only in Drop")
            .buf as *mut f64;
        // SAFETY: the buffer is writable, contiguous and `len` elements long.
        unsafe { std::slice::from_raw_parts_mut(ptr, self.len) }
    }
}

impl Drop for F64Buffer {
    fn drop(&mut self) {
        if let Some(mut view) = self.view.take() {
            // SAFETY: `view` was filled by `PyObject_GetBuffer` and has not
            // been released yet.
            unsafe {
                PyBuffer_Release(&mut view);
            }
        }
    }
}
