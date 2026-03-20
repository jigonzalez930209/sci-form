//! Experimental algorithms — developed in parallel without interfering with
//! production code. Each sub-module is gated behind its own feature flag.

#[cfg(feature = "experimental-cga")]
pub mod cga;

#[cfg(feature = "experimental-randnla")]
pub mod rand_nla;

#[cfg(feature = "experimental-riemannian")]
pub mod riemannian;

#[cfg(feature = "experimental-kpm")]
pub mod kpm;

#[cfg(feature = "experimental-eeq")]
pub mod eeq;

#[cfg(feature = "experimental-alpb")]
pub mod alpb;

#[cfg(feature = "experimental-d4")]
pub mod d4;

#[cfg(feature = "experimental-sdr")]
pub mod sdr;

#[cfg(feature = "experimental-mbh")]
pub mod mbh;

#[cfg(feature = "experimental-cpm")]
pub mod cpm;

#[cfg(feature = "experimental-gsm")]
pub mod gsm;
