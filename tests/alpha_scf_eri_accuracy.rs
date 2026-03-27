//! Integration tests for Obara-Saika electron repulsion integrals.
//!
//! Validates:
//! - Boys function F_0(0) = 1 (known limit)
//! - Boys function F_0(x→∞) ≈ √(π/x)/2 (asymptotic)
//! - (ss|ss) self-integral on same center has known analytical value
//! - Schwarz bound is an upper bound for |ERI|
//! - Screened ERI batch produces correct number of integrals
//! - ERI symmetry: (ab|cd) = (cd|ab)
//!
//! Reference: Helgaker, Jørgensen & Olsen, "Molecular Electronic-Structure Theory",
//! Chapter 9 (Coulomb integrals); McMurchie & Davidson (1978).

#![cfg(feature = "alpha-obara-saika")]

use sci_form::scf::obara_saika::{
    boys_function, compute_eris_screened, eri_ssss, schwarz_bound, ShellPairData,
};

// ─── Boys function ───────────────────────────────────────────────────────────

/// F_0(0) = 1 exactly.
#[test]
fn boys_f0_at_zero_is_one() {
    let f0 = boys_function(0, 0.0);
    assert!(
        (f0 - 1.0).abs() < 1e-10,
        "F_0(0) should be 1.0, got {:.12}",
        f0
    );
}

/// F_0(x) > 0 for all x > 0.
#[test]
fn boys_f0_positive_for_positive_x() {
    for &x in &[0.001, 0.1, 1.0, 5.0, 10.0, 50.0, 100.0] {
        let f = boys_function(0, x);
        assert!(f > 0.0, "F_0({}) should be positive, got {:.8}", x, f);
    }
}

/// F_0(x→∞) ≈ √(π/(4x)) (asymptotic limit).
#[test]
fn boys_f0_asymptotic_large_x() {
    let x = 100.0;
    let f = boys_function(0, x);
    let expected = (std::f64::consts::PI / x).sqrt() / 2.0;
    let rel_err = (f - expected).abs() / expected;
    assert!(
        rel_err < 0.01,
        "F_0(100) = {:.8}, expected ≈ {:.8}, rel error = {:.4}",
        f,
        expected,
        rel_err
    );
}

/// F_n(0) = 1/(2n+1).
#[test]
fn boys_fn_at_zero() {
    for n in 0..6 {
        let f = boys_function(n, 0.0);
        let expected = 1.0 / (2 * n + 1) as f64;
        assert!(
            (f - expected).abs() < 1e-10,
            "F_{}(0) = {:.8}, expected {:.8}",
            n,
            f,
            expected
        );
    }
}

// ─── (ss|ss) ERI ─────────────────────────────────────────────────────────────

/// Self-integral on same center: (sA sA | sA sA) with exponent α.
/// Known formula: (2α/π)^{3/2} × (π^{5/2})/(2α^2 √(4α)) × F_0(0)
///              = 2π^{5/2} / (4α² √(2α)) × 1
///
/// For α=1.0 on origin: prefactor = 2π^{5/2}/(4*sqrt(2)) ≈ 2*17.4932/(4*1.41421) ≈ 6.1685
/// But let's just check it's positive and finite.
#[test]
fn eri_ssss_self_integral_positive() {
    let sp = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);
    let eri = eri_ssss(&sp, &sp);
    assert!(
        eri > 0.0 && eri.is_finite(),
        "(ss|ss) self-integral should be positive and finite, got {:.8}",
        eri
    );
}

/// Known analytical (1s1s|1s1s) with α=β=γ=δ=1.0 on same center:
/// (ss|ss) = 2π^{5/2} / (p q √(p+q)) × K_ab × K_cd × F_0(T)
/// With α=β → p=2, K_ab = 1, same for cd → q=2, and T=0 (same center).
/// = 2π^{5/2} / (2×2×√4) = 2π^{5/2}/8 ≈ 2×17.4932/8 ≈ 4.3733
#[test]
fn eri_ssss_same_center_known_value() {
    let sp = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);
    let eri = eri_ssss(&sp, &sp);
    let pi = std::f64::consts::PI;
    let expected = 2.0 * pi.powi(2) * pi.sqrt() / 8.0;
    let rel_err = (eri - expected).abs() / expected;
    assert!(
        rel_err < 1e-6,
        "(ss|ss) same center = {:.8}, expected {:.8}, rel_err = {:.2e}",
        eri,
        expected,
        rel_err
    );
}

/// ERI symmetry: (ab|cd) = (cd|ab).
#[test]
fn eri_ssss_exchange_symmetry() {
    let sp_ab = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.5, [1.0, 0.0, 0.0]);
    let sp_cd = ShellPairData::new(0.8, [0.0, 1.0, 0.0], 1.2, [0.5, 0.5, 0.0]);

    let abcd = eri_ssss(&sp_ab, &sp_cd);
    let cdab = eri_ssss(&sp_cd, &sp_ab);

    assert!(
        (abcd - cdab).abs() < 1e-12,
        "(ab|cd) = {:.10}, (cd|ab) = {:.10} — should be equal",
        abcd,
        cdab
    );
}

/// ERI decreases with distance between shell pair centers.
#[test]
fn eri_decreases_with_distance() {
    let sp_ab = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);

    let sp_close = ShellPairData::new(1.0, [1.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]);
    let sp_far = ShellPairData::new(1.0, [5.0, 0.0, 0.0], 1.0, [5.0, 0.0, 0.0]);

    let eri_close = eri_ssss(&sp_ab, &sp_close);
    let eri_far = eri_ssss(&sp_ab, &sp_far);

    assert!(
        eri_close > eri_far,
        "ERI at 1Å ({:.6}) should be > ERI at 5Å ({:.6})",
        eri_close,
        eri_far
    );
}

// ─── Schwarz bounds ──────────────────────────────────────────────────────────

/// Schwarz bound should be a valid upper bound for |(ij|kl)|.
#[test]
fn schwarz_bound_is_upper_bound_for_eri() {
    let alpha = 1.0;
    let beta = 1.5;
    let gamma = 0.8;
    let delta = 1.2;
    let a = [0.0, 0.0, 0.0];
    let b = [1.0, 0.0, 0.0];
    let c = [0.0, 1.0, 0.0];
    let d = [0.5, 0.5, 0.0];

    let sp_ab = ShellPairData::new(alpha, a, beta, b);
    let sp_cd = ShellPairData::new(gamma, c, delta, d);

    let eri = eri_ssss(&sp_ab, &sp_cd).abs();
    let q_ab = schwarz_bound(alpha, a, beta, b);
    let q_cd = schwarz_bound(gamma, c, delta, d);

    let bound = q_ab * q_cd;
    assert!(
        eri <= bound + 1e-10,
        "|ERI| = {:.8} should be ≤ Q_ab×Q_cd = {:.8}",
        eri,
        bound
    );
}

/// Schwarz bound should be positive.
#[test]
fn schwarz_bound_positive() {
    let q = schwarz_bound(1.0, [0.0; 3], 1.0, [1.0, 0.0, 0.0]);
    assert!(q > 0.0, "Schwarz bound should be positive, got {:.8}", q);
}

// ─── Screened ERI batch ──────────────────────────────────────────────────────

/// Two shell pairs: should produce ≥1 significant integral.
#[test]
fn screened_eris_h2_produces_integrals() {
    let sp1 = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);
    let sp2 = ShellPairData::new(1.0, [0.74, 0.0, 0.0], 1.0, [0.74, 0.0, 0.0]);

    let q1 = schwarz_bound(1.0, [0.0; 3], 1.0, [0.0; 3]);
    let q2 = schwarz_bound(1.0, [0.74, 0.0, 0.0], 1.0, [0.74, 0.0, 0.0]);

    let eris = compute_eris_screened(&[sp1, sp2], &[q1, q2], 1e-10);
    assert!(
        eris.len() >= 2,
        "H₂ should have ≥2 significant ERIs, got {}",
        eris.len()
    );
}

/// Tight threshold should filter more integrals than loose threshold.
#[test]
fn screening_threshold_filters_integrals() {
    let sp1 = ShellPairData::new(1.0, [0.0, 0.0, 0.0], 1.0, [0.0, 0.0, 0.0]);
    let sp2 = ShellPairData::new(1.0, [3.0, 0.0, 0.0], 1.0, [3.0, 0.0, 0.0]);

    let q1 = schwarz_bound(1.0, [0.0; 3], 1.0, [0.0; 3]);
    let q2 = schwarz_bound(1.0, [3.0, 0.0, 0.0], 1.0, [3.0, 0.0, 0.0]);

    let eris_loose = compute_eris_screened(&[sp1.clone(), sp2.clone()], &[q1, q2], 1e-15);
    let eris_tight = compute_eris_screened(&[sp1, sp2], &[q1, q2], 1.0);

    assert!(
        eris_loose.len() >= eris_tight.len(),
        "Loose threshold should yield ≥ integrals than tight: {} vs {}",
        eris_loose.len(),
        eris_tight.len()
    );
}
