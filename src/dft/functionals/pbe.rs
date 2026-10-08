//! PBE exchange-correlation functional (GGA).
//!
//! Perdew-Burke-Ernzerhof generalized gradient approximation.
//! Reference: Phys. Rev. Lett. 77, 3865 (1996)

use std::f64::consts::PI;

/// PBE constants.
const KAPPA: f64 = 0.804;
const MU: f64 = 0.2195_1;
const BETA: f64 = 0.066_724_550_603_149_22;

/// Evaluate PBE exchange-correlation energy density and potential.
///
/// Input: electron density ρ and |∇ρ|² (gradient squared).
/// Returns: (exc, vxc_rho, vxc_sigma) where:
/// - exc: XC energy density (Hartree)
/// - vxc_rho: ∂(ρ·εxc)/∂ρ
/// - vxc_sigma: ∂(ρ·εxc)/∂(|∇ρ|²)
pub fn pbe(rho: f64, sigma: f64) -> (f64, f64, f64) {
    if rho < 1e-20 {
        return (0.0, 0.0, 0.0);
    }

    let (ex, vx_rho, vx_sigma) = pbe_exchange(rho, sigma);
    let (ec, vc_rho, vc_sigma) = pbe_correlation(rho, sigma);

    (ex + ec, vx_rho + vc_rho, vx_sigma + vc_sigma)
}

/// PBE exchange enhancement factor.
///
/// F_x(s) = 1 + κ - κ / (1 + μs²/κ)
/// where s = |∇ρ| / (2 k_F ρ), k_F = (3π²ρ)^{1/3}
fn pbe_exchange(rho: f64, sigma: f64) -> (f64, f64, f64) {
    let rho_third = rho.powf(1.0 / 3.0);
    let _rho_four_third = rho * rho_third;

    // LDA exchange energy density
    let cx = -0.75 * (3.0 / PI).powf(1.0 / 3.0);
    let ex_lda = cx * rho_third;

    // Fermi wavevector and reduced gradient
    let kf = (3.0 * PI * PI * rho).powf(1.0 / 3.0);
    let grad_rho = sigma.sqrt();
    let s = grad_rho / (2.0 * kf * rho + 1e-30);
    let s2 = s * s;

    // Enhancement factor
    let denom = 1.0 + MU * s2 / KAPPA;
    let fx = 1.0 + KAPPA - KAPPA / denom;

    let ex = ex_lda * fx;

    // ∂(ρεx)/∂ρ at fixed σ. s depends on ρ, so the chain rule is numerical.
    let delta = 1e-7 * rho.max(1e-12);
    let ex_plus = pbe_exchange_only(rho + delta, sigma);
    let ex_minus = pbe_exchange_only((rho - delta).max(1e-20), sigma);
    let vx_rho = ex + rho * (ex_plus - ex_minus) / (2.0 * delta);

    let dfx_ds2 = MU / (denom * denom);
    let ds2_dsigma = 1.0 / (4.0 * kf * kf * rho * rho + 1e-30);
    let vx_sigma = ex_lda * rho * dfx_ds2 * ds2_dsigma;

    (ex, vx_rho, vx_sigma)
}

fn pbe_exchange_only(rho: f64, sigma: f64) -> f64 {
    if rho < 1e-20 {
        return 0.0;
    }
    let cx = -0.75 * (3.0 / PI).powf(1.0 / 3.0);
    let ex_lda = cx * rho.powf(1.0 / 3.0);
    let kf = (3.0 * PI * PI * rho).powf(1.0 / 3.0);
    let s = sigma.sqrt() / (2.0 * kf * rho + 1e-30);
    let denom = 1.0 + MU * s * s / KAPPA;
    let fx = 1.0 + KAPPA - KAPPA / denom;
    ex_lda * fx
}

/// PBE correlation.
fn pbe_correlation(rho: f64, sigma: f64) -> (f64, f64, f64) {
    if rho < 1e-20 {
        return (0.0, 0.0, 0.0);
    }

    // Use PW correlation (simplified)
    let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);

    // Simplified PW92 correlation energy
    let a = 0.031091;
    let alpha1 = 0.21370;
    let beta1 = 7.5957;
    let beta2 = 3.5876;
    let beta3 = 1.6382;
    let beta4 = 0.49294;

    let rs_sqrt = rs.sqrt();
    let rs_32 = rs * rs_sqrt;

    let denom = 2.0 * a * (beta1 * rs_sqrt + beta2 * rs + beta3 * rs_32 + beta4 * rs * rs);
    let ec_lda = -2.0 * a * (1.0 + alpha1 * rs) * (1.0 + 1.0 / denom).ln();

    let kf = (3.0 * PI * PI * rho).powf(1.0 / 3.0);
    let ec = ec_lda + pbe_h(ec_lda, rho, sigma, kf);

    // Numerical derivatives for robustness
    let delta = 1e-7 * rho.max(1e-12);
    let (ec_plus, _, _) = pbe_correlation_energy(rho + delta, sigma);
    let (ec_minus, _, _) = pbe_correlation_energy(rho - delta, sigma);
    let vc_rho = ec + rho * (ec_plus - ec_minus) / (2.0 * delta);

    let delta_s = 1e-7 * sigma.max(1e-12);
    let (ec_sp, _, _) = pbe_correlation_energy(rho, sigma + delta_s);
    let (ec_sm, _, _) = pbe_correlation_energy(rho, sigma - delta_s);
    let vc_sigma = rho * (ec_sp - ec_sm) / (2.0 * delta_s);

    (ec, vc_rho, vc_sigma)
}

/// Helper for numerical derivatives.
fn pbe_correlation_energy(rho: f64, sigma: f64) -> (f64, f64, f64) {
    if rho < 1e-20 {
        return (0.0, 0.0, 0.0);
    }

    let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
    let a = 0.031091;
    let alpha1 = 0.21370;
    let beta1 = 7.5957;
    let beta2 = 3.5876;
    let beta3 = 1.6382;
    let beta4 = 0.49294;

    let rs_sqrt = rs.sqrt();
    let rs_32 = rs * rs_sqrt;
    let denom = 2.0 * a * (beta1 * rs_sqrt + beta2 * rs + beta3 * rs_32 + beta4 * rs * rs);
    let ec_lda = -2.0 * a * (1.0 + alpha1 * rs) * (1.0 + 1.0 / denom).ln();

    let kf = (3.0 * PI * PI * rho).powf(1.0 / 3.0);
    (ec_lda + pbe_h(ec_lda, rho, sigma, kf), 0.0, 0.0)
}

/// PBE correlation gradient correction H.
///
/// H = γ ln(1 + (β/γ) t² (1+At²)/(1+At²+A²t⁴)),
/// γ = (1−ln 2)/π², A = (β/γ) / (exp(−εc/γ) − 1).
fn pbe_h(ec_lda: f64, rho: f64, sigma: f64, kf: f64) -> f64 {
    let gamma = (1.0 - 2.0_f64.ln()) / (PI * PI);
    let ks = (4.0 * kf / PI).sqrt();
    let t2 = sigma / (4.0 * ks * ks * rho * rho + 1e-30);
    let denom = ((-ec_lda / gamma).exp() - 1.0).max(1e-30);
    let a_pbe = (BETA / gamma) / denom;
    let at2 = a_pbe * t2;
    let frac = t2 * (1.0 + at2) / (1.0 + at2 + at2 * at2);
    gamma * (1.0 + (BETA / gamma) * frac).ln()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pbe_zero_density() {
        let (exc, vr, vs) = pbe(0.0, 0.0);
        assert_eq!(exc, 0.0);
        assert_eq!(vr, 0.0);
        assert_eq!(vs, 0.0);
    }

    #[test]
    fn pbe_exchange_is_negative() {
        let (exc, _, _) = pbe(0.1, 0.0);
        assert!(exc < 0.0);
    }
}
