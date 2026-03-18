//! Smooth cutoff functions for distance-dependent interactions.
//!
//! Ensures that energy contributions vanish smoothly at the cutoff radius,
//! preventing discontinuities in the potential energy surface.

use std::f64::consts::PI;

/// Cosine cutoff function.
///
/// Returns 1.0 at r=0, smoothly decays to 0.0 at r=rc.
/// $$f_c(r) = \frac{1}{2}\cos\left(\frac{\pi r}{R_c}\right) + \frac{1}{2}$$
#[inline]
pub fn cosine_cutoff(r: f64, rc: f64) -> f64 {
    if r >= rc {
        0.0
    } else {
        0.5 * (PI * r / rc).cos() + 0.5
    }
}

/// Derivative of the cosine cutoff function with respect to r.
///
/// $$\frac{df_c}{dr} = -\frac{\pi}{2 R_c}\sin\left(\frac{\pi r}{R_c}\right)$$
#[inline]
pub fn cosine_cutoff_deriv(r: f64, rc: f64) -> f64 {
    if r >= rc {
        0.0
    } else {
        -0.5 * PI / rc * (PI * r / rc).sin()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cutoff_boundary_values() {
        let rc = 5.2;
        // At r=0, fc = 1.0
        assert!((cosine_cutoff(0.0, rc) - 1.0).abs() < 1e-12);
        // At r=rc, fc = 0.0
        assert!(cosine_cutoff(rc, rc).abs() < 1e-12);
        // Beyond rc, fc = 0.0
        assert_eq!(cosine_cutoff(rc + 1.0, rc), 0.0);
    }

    #[test]
    fn test_cutoff_monotonic_decrease() {
        let rc = 5.2;
        let mut prev = cosine_cutoff(0.0, rc);
        for i in 1..=100 {
            let r = rc * i as f64 / 100.0;
            let val = cosine_cutoff(r, rc);
            assert!(val <= prev + 1e-12, "Cutoff not monotonically decreasing");
            prev = val;
        }
    }

    #[test]
    fn test_cutoff_deriv_at_zero() {
        let rc = 5.2;
        assert!(cosine_cutoff_deriv(0.0, rc).abs() < 1e-12);
    }

    #[test]
    fn test_cutoff_deriv_numerical() {
        let rc = 5.2;
        let r = 2.5;
        let h = 1e-7;
        let numerical = (cosine_cutoff(r + h, rc) - cosine_cutoff(r - h, rc)) / (2.0 * h);
        let analytical = cosine_cutoff_deriv(r, rc);
        assert!(
            (numerical - analytical).abs() < 1e-6,
            "Numerical={numerical}, Analytical={analytical}"
        );
    }
}
