//! Analytical potentials for reaction path benchmarking and testing.

/// 1D double well potential: V(x) = (x^2 - 1)^2
/// Symmetric, minima at +/- 1, barrier at x=0
pub fn double_well_1d(x: f64) -> f64 {
    (x * x - 1.0).powi(2)
}

/// 1D asymmetric potential: V(x) = (x^2 - 1)^2 - 0.3x
/// Exothermic with product at x=+1 being more stable
pub fn asymmetric_1d(x: f64) -> f64 {
    (x * x - 1.0).powi(2) - 0.3 * x
}

/// 1D SN2 model potential: V(x) = -2*exp(-2(x+1.5)^2) - 3*exp(-2(x-1.5)^2) + 1.5*exp(-8x^2)
pub fn sn2_model_1d(x: f64) -> f64 {
    -2.0 * (-(2.0 * (x + 1.5).powi(2))).exp()
        - 3.0 * (-(2.0 * (x - 1.5).powi(2))).exp()
        + 1.5 * (-(8.0 * x.powi(2))).exp()
}

/// 2D Müller-Brown potential
pub fn muller_brown_2d(x: f64, y: f64) -> f64 {
    let a = [-200.0_f64, -100., -170., 15.];
    let b = [-1.0, -1., -6.5, 0.7];
    let c = [0.0, 0., 11.0, 0.6];
    let d = [-10.0, -10., -6.5, 0.7];
    let x0 = [1.0, 0., -0.5, -1.0];
    let y0 = [0.0, 0.5, 1.5, 1.0];
    (0..4)
        .map(|k| {
            let dx = x - x0[k];
            let dy = y - y0[k];
            a[k] * (b[k] * dx * dx + c[k] * dx * dy + d[k] * dy * dy).exp()
        })
        .sum::<f64>()
        * 0.001
}

/// 1D Marcus/EVB adiabatic potential for proton transfer
/// s in [-1, 1], symmetric TS exactly at s=0
pub fn proton_transfer_evb(s: f64, k: f64, coupling: f64) -> f64 {
    let eps_d = k * (s + 1.0).powi(2);
    let eps_a = k * (s - 1.0).powi(2);
    let sum2 = (eps_d + eps_a) / 2.0;
    let diff2 = (eps_d - eps_a) / 2.0;
    sum2 - (diff2 * diff2 + coupling * coupling).sqrt()
}

/// 2D Morse-based potential for proton transfer A-H...B
pub fn proton_transfer_2d_morse(r_ah: f64, r_hb: f64) -> f64 {
    let de = 2.0;
    let alpha = 2.0;
    let r0 = 1.0;
    let morse_ah = de * (1.0 - (-alpha * (r_ah - r0)).exp()).powi(2);
    let morse_hb = de * (1.0 - (-alpha * (r_hb - r0)).exp()).powi(2);
    let constraint = 50.0 * (r_ah + r_hb - 2.5).powi(2);
    morse_ah + morse_hb + constraint
}
