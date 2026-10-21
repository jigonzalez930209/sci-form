//! Aggregated test target for the `tests/experimental/` suite.
//!
//! These files live in a sub-directory (so Cargo does not auto-discover them)
//! and cover the promoted-to-core experimental modules (EEQ, ALPB, D4, CPM)
//! plus the alpha/beta-gated ones (CGA, GSM, SDR, KPM, MBH, RandNLA,
//! Riemannian).
//!
//! The feature gates below use the CURRENT feature names. The historical
//! `experimental-*` flags are no-ops that do not enable the underlying
//! modules, so gating on them compiled these files to empty and silently
//! skipped every test inside them.

#[path = "experimental/test_alpb.rs"]
mod test_alpb;

#[path = "experimental/test_benchmarks.rs"]
mod test_benchmarks;

#[path = "experimental/test_cpm.rs"]
mod test_cpm;

#[path = "experimental/test_d4.rs"]
mod test_d4;

#[path = "experimental/test_eeq.rs"]
mod test_eeq;

#[cfg(feature = "alpha-cga")]
#[path = "experimental/test_cga.rs"]
mod test_cga;

#[cfg(feature = "alpha-gsm")]
#[path = "experimental/test_gsm.rs"]
mod test_gsm;

#[cfg(feature = "beta-kpm")]
#[path = "experimental/test_kpm.rs"]
mod test_kpm;

#[cfg(feature = "beta-mbh")]
#[path = "experimental/test_mbh.rs"]
mod test_mbh;

#[cfg(feature = "beta-randnla")]
#[path = "experimental/test_randnla.rs"]
mod test_randnla;

#[cfg(feature = "beta-riemannian")]
#[path = "experimental/test_riemannian.rs"]
mod test_riemannian;

#[cfg(feature = "alpha-sdr")]
#[path = "experimental/test_sdr.rs"]
mod test_sdr;
