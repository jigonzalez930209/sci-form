//! Hermite tapering polynomial for smooth distance cutoff.
//!
//! Provides C³-continuous tapering from r=0 (tap=1) to r=cutoff (tap=0).
//! Used in ReaxFF for smooth non-bonded and bond-order cutoffs.

/// Seventh-order tapering polynomial (C³ continuity at cutoff).
///
/// tap(r) = 1 for r = 0
/// tap(r) = 0 for r ≥ r_cut
/// tap, tap', tap'', tap''' are continuous at r = r_cut
#[inline]
pub fn taper_function(r: f64, r_cut: f64) -> f64 {
    if r >= r_cut {
        return 0.0;
    }
    if r <= 0.0 {
        return 1.0;
    }

    let rr = r / r_cut;
    let rr3 = rr * rr * rr;
    let rr4 = rr3 * rr;
    let rr5 = rr4 * rr;
    let rr6 = rr5 * rr;
    let rr7 = rr6 * rr;

    // Coefficients for 7th-order polynomial ensuring C³ continuity
    20.0 * rr7 - 70.0 * rr6 + 84.0 * rr5 - 35.0 * rr4 + 1.0
}

/// Derivative of the taper function d(tap)/dr.
#[inline]
pub fn taper_derivative(r: f64, r_cut: f64) -> f64 {
    if r >= r_cut || r <= 0.0 {
        return 0.0;
    }

    let rr = r / r_cut;
    let rr2 = rr * rr;
    let rr3 = rr2 * rr;
    let rr4 = rr3 * rr;
    let rr5 = rr4 * rr;
    let rr6 = rr5 * rr;

    (140.0 * rr6 - 420.0 * rr5 + 420.0 * rr4 - 140.0 * rr3) / r_cut
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn taper_boundary_conditions() {
        let r_cut = 10.0;
        assert!((taper_function(0.0, r_cut) - 1.0).abs() < 1e-12);
        assert!(taper_function(r_cut, r_cut).abs() < 1e-12);
        assert!(taper_function(r_cut + 1.0, r_cut).abs() < 1e-12);
    }

    #[test]
    fn taper_is_monotone_decreasing() {
        let r_cut = 10.0;
        let mut prev = 1.0;
        for i in 1..100 {
            let r = i as f64 * r_cut / 100.0;
            let t = taper_function(r, r_cut);
            assert!(t <= prev + 1e-12);
            prev = t;
        }
    }
}
