//! SVWN functional: Slater exchange + VWN-5 correlation (LDA).
//!
//! S: Dirac-Slater exchange (exact for uniform electron gas)
//! VWN: Vosko-Wilk-Nusair parametrization V (1980)

use std::f64::consts::PI;

/// Evaluate SVWN exchange-correlation energy density and potential.
///
/// Input: electron density ρ at a grid point.
/// Returns: (exc, vxc) where exc is XC energy density (Hartree) and vxc is XC potential.
pub fn svwn(rho: f64) -> (f64, f64) {
    if rho < 1e-20 {
        return (0.0, 0.0);
    }

    let (ex, vx) = slater_exchange(rho);
    let (ec, vc) = vwn5_correlation(rho);

    (ex + ec, vx + vc)
}

/// Slater (Dirac) exchange for the uniform electron gas.
///
/// εx = -3/4 (3/π)^{1/3} ρ^{1/3}
/// Vx = 4/3 εx
fn slater_exchange(rho: f64) -> (f64, f64) {
    let cx = -0.75 * (3.0 / PI).powf(1.0 / 3.0);
    let rho_third = rho.powf(1.0 / 3.0);

    let ex = cx * rho_third;
    let vx = 4.0 / 3.0 * ex;

    (ex, vx)
}

/// VWN-5 correlation functional (Vosko-Wilk-Nusair, parametrization V).
///
/// Reference: S.H. Vosko, L. Wilk, M. Nusair, Can. J. Phys. 58, 1200 (1980)
fn vwn5_correlation(rho: f64) -> (f64, f64) {
    // Wigner-Seitz radius: r_s = (3 / (4πρ))^{1/3}
    let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);

    // VWN-5 parameters for paramagnetic case
    let a = 0.0621814;
    let x0 = -0.10498;
    let b = 3.72744;
    let c = 12.9352;

    let x = rs.sqrt();
    let x_x0 = x - x0;
    let xx = x * x + b * x + c;
    let xx0 = x0 * x0 + b * x0 + c;
    let q = (4.0 * c - b * b).sqrt();

    let _ec = a
        * (x.powi(2).ln() / xx.ln() + 2.0 * b / q * (2.0 * x + b).atan2(q).atan()
            - b * x0 / xx0 * ((x_x0).powi(2) / xx).ln()
            + 2.0 * (b + 2.0 * x0) / q * (2.0 * x + b).atan2(q).atan());

    // Simplified VWN-5: use the standard parameterization
    // ε_c(r_s) = A/2 { ln(x²/X(x)) + 2b/Q arctan(Q/(2x+b))
    //            - bx0/X(x0) [ ln((x-x0)²/X(x)) + 2(b+2x0)/Q arctan(Q/(2x+b)) ] }
    let atan_term = (q / (2.0 * x + b)).atan();
    let ec_real = 0.5
        * a
        * (2.0 * (x * x / xx).ln() + 2.0 * b / q * atan_term
            - b * x0 / xx0
                * (2.0 * ((x - x0).powi(2) / xx).ln() + 2.0 * (b + 2.0 * x0) / q * atan_term));

    // Numerical derivative for V_c = ε_c - (r_s / 3) dε_c/dr_s
    let delta = 1e-6 * rs.max(1e-10);
    let rho_plus = 3.0 / (4.0 * PI * (rs + delta).powi(3));
    let rho_minus = 3.0 / (4.0 * PI * (rs - delta).powi(3));
    let (ec_plus, _) = vwn5_correlation_energy(rho_plus);
    let (ec_minus, _) = vwn5_correlation_energy(rho_minus);

    // Use the real ec value
    let vc = ec_real + rho * (ec_plus - ec_minus) / (rho_plus - rho_minus);

    (ec_real, vc)
}

/// Helper: VWN-5 energy density (without potential computation).
fn vwn5_correlation_energy(rho: f64) -> (f64, f64) {
    if rho < 1e-20 {
        return (0.0, 0.0);
    }
    let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);

    let a = 0.0621814;
    let x0 = -0.10498;
    let b = 3.72744;
    let c = 12.9352;

    let x = rs.sqrt();
    let xx = x * x + b * x + c;
    let xx0 = x0 * x0 + b * x0 + c;
    let q = (4.0 * c - b * b).sqrt();
    let atan_term = (q / (2.0 * x + b)).atan();

    let ec = 0.5
        * a
        * (2.0 * (x * x / xx).ln() + 2.0 * b / q * atan_term
            - b * x0 / xx0
                * (2.0 * ((x - x0).powi(2) / xx).ln() + 2.0 * (b + 2.0 * x0) / q * atan_term));

    (ec, 0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn svwn_zero_density() {
        let (exc, vxc) = svwn(0.0);
        assert_eq!(exc, 0.0);
        assert_eq!(vxc, 0.0);
    }

    #[test]
    fn slater_exchange_is_negative() {
        let (ex, vx) = slater_exchange(0.1);
        assert!(ex < 0.0);
        assert!(vx < 0.0);
    }
}
