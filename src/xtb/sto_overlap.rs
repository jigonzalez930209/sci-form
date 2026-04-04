//! STO-nG overlap integrals for tight-binding solvers.
//!
//! Expands Slater-type orbitals (STOs) in a small number of Gaussian
//! primitives (STO-nG), then computes overlap analytically using the
//! Obara-Saika recursion.  This gives proper angular-momentum–dependent
//! overlap integrals instead of the spherical 1s approximation.
//!
//! STO-3G coefficients from:
//!   Hehre, Stewart, Pople, JCP 51, 2657 (1969)
//!   Stewart, JCP 52, 431 (1970)

use std::f64::consts::PI;

// ── STO-3G expansion coefficients ──
// From: Robert F. Stewart, J. Chem. Phys. 52, 431-438 (1970), Table I-V.
// Independent exponents per (n, l) shell.
// Each entry: (exponent_ratio, contraction_coefficient)
// Actual Gaussian exponent = ratio × ζ²

/// 1s STO expanded in 3 Gaussians.
static STO3G_1S: [(f64, f64); 3] = [
    (2.227660584e+0, 1.543289673e-1),
    (4.057711562e-1, 5.353281423e-1),
    (1.098175104e-1, 4.446345422e-1),
];

/// 2s STO expanded in 3 Gaussians.
static STO3G_2S: [(f64, f64); 3] = [
    (2.581578398e+0, -5.994474934e-2),
    (1.567622104e-1, 5.960385398e-1),
    (6.018332272e-2, 4.581786291e-1),
];

/// 2p STO expanded in 3 Gaussians.
static STO3G_2P: [(f64, f64); 3] = [
    (9.192379002e-1, 1.623948553e-1),
    (2.359194503e-1, 5.661708862e-1),
    (8.009805746e-2, 4.223071752e-1),
];

/// 3s STO expanded in 3 Gaussians.
static STO3G_3S: [(f64, f64); 3] = [
    (5.641487709e-1, -1.782577972e-1),
    (6.924421391e-2, 8.612761663e-1),
    (3.269529097e-2, 2.261841969e-1),
];

/// 3p STO expanded in 3 Gaussians.
static STO3G_3P: [(f64, f64); 3] = [
    (2.692880368e+0, -1.061945788e-2),
    (1.489359592e-1, 5.218564264e-1),
    (5.739585040e-2, 5.450015143e-1),
];

/// 3d STO expanded in 3 Gaussians.
static STO3G_3D: [(f64, f64); 3] = [
    (5.229112225e-1, 1.686596060e-1),
    (1.639595876e-1, 5.847984817e-1),
    (6.386630021e-2, 4.056779523e-1),
];

/// 4s STO expanded in 3 Gaussians.
static STO3G_4S: [(f64, f64); 3] = [
    (2.267938753e-1, -3.349048323e-1),
    (4.448178019e-2, 1.056744667e+0),
    (2.195294664e-2, 1.256661680e-1),
];

/// 4p STO expanded in 3 Gaussians.
static STO3G_4P: [(f64, f64); 3] = [
    (4.859692220e-1, -6.147823411e-2),
    (7.430216918e-2, 6.604172234e-1),
    (3.653340923e-2, 3.932639495e-1),
];

/// 4d STO expanded in 3 Gaussians.
static STO3G_4D: [(f64, f64); 3] = [
    (1.777717219e-1, 2.308552718e-1),
    (8.040647350e-2, 6.042409177e-1),
    (3.949855551e-2, 2.595768926e-1),
];

/// 5s STO expanded in 3 Gaussians.
static STO3G_5S: [(f64, f64); 3] = [
    (1.080198458e-1, -6.617401158e-1),
    (4.408119382e-2, 7.467595004e-1),
    (2.610811810e-2, 7.146490945e-1),
];

/// 5p STO expanded in 3 Gaussians.
static STO3G_5P: [(f64, f64); 3] = [
    (2.127482317e-1, -1.389529695e-1),
    (4.729648620e-2, 8.076691064e-1),
    (2.604865324e-2, 2.726029342e-1),
];

/// 5d STO expanded in 3 Gaussians.
static STO3G_5D: [(f64, f64); 3] = [
    (4.913352950e-1, -2.010175008e-2),
    (7.329090601e-2, 5.899370608e-1),
    (3.594209290e-2, 4.658445960e-1),
];

// ── STO-4G expansion coefficients ──
// From: Robert F. Stewart, J. Chem. Phys. 52, 431-438 (1970), Table I-V.
// Independent exponents per (n, l) shell, same as tblite slater.f90 pAlpha4/pCoeff4.
// (exponent_ratio, contraction_coefficient)

/// 1s STO-4G.
static STO4G_1S: [(f64, f64); 4] = [
    (5.216844534e+0, 5.675242080e-2),
    (9.546182760e-1, 2.601413550e-1),
    (2.652034102e-1, 5.328461143e-1),
    (8.801862774e-2, 2.916254405e-1),
];

/// 2s STO-4G.
static STO4G_2S: [(f64, f64); 4] = [
    (1.161525551e+1, -1.198411747e-2),
    (2.000243111e+0, -5.472052539e-2),
    (1.607280687e-1, 5.805587176e-1),
    (6.125744532e-2, 4.770079976e-1),
];

/// 2p STO-4G.
static STO4G_2P: [(f64, f64); 4] = [
    (1.798260992e+0, 5.713170255e-2),
    (4.662622228e-1, 2.857455515e-1),
    (1.643718620e-1, 5.517873105e-1),
    (6.543927065e-2, 2.632314924e-1),
];

/// 3s STO-4G.
static STO4G_3S: [(f64, f64); 4] = [
    (1.513265591e+0, -3.295496352e-2),
    (4.262497508e-1, -1.724516959e-1),
    (7.643320863e-2, 7.518511194e-1),
    (3.760545063e-2, 3.589627317e-1),
];

/// 3p STO-4G.
static STO4G_3P: [(f64, f64); 4] = [
    (1.853180239e+0, -1.434249391e-2),
    (1.915075719e-1, 2.755177589e-1),
    (8.655487938e-2, 5.846750879e-1),
    (4.184253862e-2, 2.144986514e-1),
];

/// 3d STO-4G.
static STO4G_3D: [(f64, f64); 4] = [
    (9.185846715e-1, 5.799057705e-2),
    (2.920461109e-1, 3.045581349e-1),
    (1.187568890e-1, 5.601358038e-1),
    (5.286755896e-2, 2.432423313e-1),
];

/// 4s STO-4G.
static STO4G_4S: [(f64, f64); 4] = [
    (3.242212833e-1, -1.120682822e-1),
    (1.663217177e-1, -2.845426863e-1),
    (5.081097451e-2, 8.909873788e-1),
    (2.829066600e-2, 3.517811205e-1),
];

/// 4p STO-4G.
static STO4G_4P: [(f64, f64); 4] = [
    (1.492607880e+0, -6.035216774e-3),
    (4.327619272e-1, -6.013310874e-2),
    (7.553156064e-2, 6.451518200e-1),
    (3.706272183e-2, 4.117923820e-1),
];

/// 4d STO-4G.
static STO4G_4D: [(f64, f64); 4] = [
    (1.995825422e+0, -2.816702620e-3),
    (1.823461280e-1, 2.177095871e-1),
    (8.197240896e-2, 6.058047348e-1),
    (4.000634951e-2, 2.717811257e-1),
];

/// 5s STO-4G.
static STO4G_5S: [(f64, f64); 4] = [
    (8.602284252e-1, 1.103657561e-2),
    (1.189050200e-1, -5.606519023e-1),
    (3.446076176e-2, 1.179429987e+0),
    (1.974798796e-2, 1.734974376e-1),
];

/// 5p STO-4G.
static STO4G_5P: [(f64, f64); 4] = [
    (3.962838833e-1, -1.801459207e-2),
    (1.838858552e-1, -1.360777372e-1),
    (4.943555157e-2, 7.533973719e-1),
    (2.750222273e-2, 3.409304859e-1),
];

/// 5d STO-4G.
static STO4G_5D: [(f64, f64); 4] = [
    (4.230617826e-1, -2.421626009e-2),
    (8.293863702e-2, 3.937644956e-1),
    (4.590326388e-2, 5.489520286e-1),
    (2.628744797e-2, 1.190436963e-1),
];

// ── STO-6G expansion coefficients ──
// From: Robert F. Stewart, J. Chem. Phys. 52, 431-438 (1970), Table I-V.
// tblite slater.f90 pAlpha6s/pCoeff6s and pAlpha6p/pCoeff6p for n=6.

/// 6s STO-6G (special n=6 s).
static STO6G_6S: [(f64, f64); 6] = [
    (5.800292686e-1, 4.554359511e-3),
    (2.718262251e-1, 5.286443143e-2),
    (7.938523262e-2, -7.561016358e-1),
    (4.975088254e-2, -2.269803820e-1),
    (2.983643556e-2, 1.332494651e+0),
    (1.886067216e-2, 3.622518293e-1),
];

/// 6p STO-6G (special n=6 p).
static STO6G_6P: [(f64, f64); 6] = [
    (6.696537714e-1, 2.782723680e-3),
    (1.395089793e-1, -1.282887780e-1),
    (8.163894960e-2, -2.266255943e-1),
    (4.586329272e-2, 4.682259383e-1),
    (2.961305556e-2, 6.752048848e-1),
    (1.882221321e-2, 1.091534212e-1),
];

/// Look up STO-3G expansion for given principal quantum number and angular momentum.
fn sto3g_coefficients(n: u8, l: u8) -> &'static [(f64, f64); 3] {
    match (n, l) {
        (1, 0) => &STO3G_1S,
        (2, 0) => &STO3G_2S,
        (2, 1) => &STO3G_2P,
        (3, 0) => &STO3G_3S,
        (3, 1) => &STO3G_3P,
        (3, 2) => &STO3G_3D,
        (4, 0) => &STO3G_4S,
        (4, 1) => &STO3G_4P,
        (4, 2) => &STO3G_4D,
        (5, 0) => &STO3G_5S,
        (5, 1) => &STO3G_5P,
        (5, 2) => &STO3G_5D,
        _ => &STO3G_1S, // fallback
    }
}

/// Look up STO-4G expansion for given principal quantum number and angular momentum.
fn sto4g_coefficients(n: u8, l: u8) -> &'static [(f64, f64); 4] {
    match (n, l) {
        (1, 0) => &STO4G_1S,
        (2, 0) => &STO4G_2S,
        (2, 1) => &STO4G_2P,
        (3, 0) => &STO4G_3S,
        (3, 1) => &STO4G_3P,
        (3, 2) => &STO4G_3D,
        (4, 0) => &STO4G_4S,
        (4, 1) => &STO4G_4P,
        (4, 2) => &STO4G_4D,
        (5, 0) => &STO4G_5S,
        (5, 1) => &STO4G_5P,
        (5, 2) => &STO4G_5D,
        _ => &STO4G_1S, // fallback
    }
}

/// Map (l, m_index) to a list of Cartesian angular momentum components with coefficients.
///
/// For s and p orbitals, returns a single Cartesian component.
/// For d orbitals, returns the proper spherical→Cartesian decomposition:
///   m=0 (dxy): xy
///   m=1 (dxz): xz
///   m=2 (dyz): yz
///   m=3 (dx²-y²): sqrt(3)/2 * xx - sqrt(3)/2 * yy
///   m=4 (dz²): -0.5 * xx - 0.5 * yy + 1.0 * zz
///
/// The coefficients are for integrals computed with per-Cartesian normalization
/// (each Cartesian Gaussian individually normalized). The tblite `dtrafo` matrix
/// uses a common normalization for all d Cartesians (equal to the xx normalization),
/// so xy/xz/yz factors of sqrt(3) cancel with the normalization ratio sqrt(3).
fn spherical_to_cartesian_components(l: u8, m: u8) -> Vec<((u8, u8, u8), f64)> {
    match l {
        0 => vec![((0, 0, 0), 1.0)],
        1 => match m {
            0 => vec![((1, 0, 0), 1.0)], // px
            1 => vec![((0, 1, 0), 1.0)], // py
            2 => vec![((0, 0, 1), 1.0)], // pz
            _ => vec![((0, 0, 0), 1.0)],
        },
        2 => {
            let s3_4: f64 = (3.0_f64 / 4.0).sqrt(); // sqrt(3)/2
            match m {
                0 => vec![((1, 1, 0), 1.0)],                                       // dxy
                1 => vec![((1, 0, 1), 1.0)],                                       // dxz
                2 => vec![((0, 1, 1), 1.0)],                                       // dyz
                3 => vec![((2, 0, 0), s3_4), ((0, 2, 0), -s3_4)],                  // dx²-y²
                4 => vec![((2, 0, 0), -0.5), ((0, 2, 0), -0.5), ((0, 0, 2), 1.0)], // dz²
                _ => vec![((0, 0, 0), 1.0)],
            }
        }
        _ => vec![((0, 0, 0), 1.0)],
    }
}

/// Legacy single-Cartesian mapping for s and p orbitals (not used for d).
#[allow(dead_code)]
fn angular_to_cartesian(l: u8, m: u8) -> (u8, u8, u8) {
    match l {
        0 => (0, 0, 0),
        1 => match m {
            0 => (1, 0, 0), // px
            1 => (0, 1, 0), // py
            2 => (0, 0, 1), // pz
            _ => (0, 0, 0),
        },
        // For d-orbitals, use spherical_to_cartesian_components instead.
        // These are kept only as a fallback; d-orbitals should never reach here.
        2 => match m {
            0 => (1, 1, 0), // dxy
            1 => (1, 0, 1), // dxz
            2 => (0, 1, 1), // dyz
            3 => (2, 0, 0), // dx²-y² — WRONG for standalone use
            4 => (0, 0, 2), // dz² — WRONG for standalone use
            _ => (0, 0, 0),
        },
        _ => (0, 0, 0),
    }
}

/// Double factorial: (2n-1)!! = 1×3×5×...×(2n-1). Returns 1 for n≤0.
fn double_factorial(n: i32) -> f64 {
    if n <= 0 {
        return 1.0;
    }
    let mut result = 1.0;
    let mut k = n;
    while k > 0 {
        result *= k as f64;
        k -= 2;
    }
    result
}

/// Normalization constant for a Cartesian Gaussian primitive.
///
/// g(r) = N × x^lx × y^ly × z^lz × exp(-α|r-A|²)
fn gaussian_norm(lx: u8, ly: u8, lz: u8, alpha: f64) -> f64 {
    let l = (lx + ly + lz) as i32;
    let num = (4.0 * alpha).powi(l) * (2.0 * alpha / PI).powf(1.5);
    let den = double_factorial(2 * lx as i32 - 1)
        * double_factorial(2 * ly as i32 - 1)
        * double_factorial(2 * lz as i32 - 1);
    (num / den).sqrt()
}

/// Obara-Saika 1D overlap recursion.
///
/// Computes S(la, lb) for one Cartesian component.
/// `pa` = P_x - A_x, `pb` = P_x - B_x, where P is the Gaussian product center.
fn overlap_1d(la: u8, lb: u8, pa: f64, pb: f64, gamma: f64) -> f64 {
    let s00 = (PI / gamma).sqrt();
    if la == 0 && lb == 0 {
        return s00;
    }

    let la = la as usize;
    let lb = lb as usize;

    // Build overlap table s[i][j] using Obara-Saika recursion:
    // s[i+1][j] = pa × s[i][j] + 1/(2γ) × (i × s[i-1][j] + j × s[i][j-1])
    let inv2g = 0.5 / gamma;
    let mut s = vec![vec![0.0f64; lb + 1]; la + 1];
    s[0][0] = s00;

    // Fill first column: s[i][0]
    for i in 0..la {
        s[i + 1][0] = pa * s[i][0]
            + if i > 0 {
                inv2g * (i as f64) * s[i - 1][0]
            } else {
                0.0
            };
    }

    // Fill remaining columns
    for j in 0..lb {
        // s[0][j+1]
        s[0][j + 1] = pb * s[0][j]
            + if j > 0 {
                inv2g * (j as f64) * s[0][j - 1]
            } else {
                0.0
            };
        // s[i][j+1] for i > 0
        for i in 1..=la {
            s[i][j + 1] = pa * s[i - 1][j + 1]
                + inv2g
                    * ((i - 1) as f64 * if i >= 2 { s[i - 2][j + 1] } else { 0.0 }
                        + (j + 1) as f64 * s[i - 1][j]);
        }
    }

    s[la][lb]
}

/// Overlap between two normalized Cartesian Gaussian primitives.
fn primitive_overlap(
    lx_a: u8,
    ly_a: u8,
    lz_a: u8,
    alpha_a: f64,
    a: &[f64; 3],
    lx_b: u8,
    ly_b: u8,
    lz_b: u8,
    alpha_b: f64,
    b: &[f64; 3],
) -> f64 {
    let gamma = alpha_a + alpha_b;
    let mu = alpha_a * alpha_b / gamma;
    let r2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
    let k_ab = (-mu * r2).exp();

    let p = [
        (alpha_a * a[0] + alpha_b * b[0]) / gamma,
        (alpha_a * a[1] + alpha_b * b[1]) / gamma,
        (alpha_a * a[2] + alpha_b * b[2]) / gamma,
    ];

    let sx = overlap_1d(lx_a, lx_b, p[0] - a[0], p[0] - b[0], gamma);
    let sy = overlap_1d(ly_a, ly_b, p[1] - a[1], p[1] - b[1], gamma);
    let sz = overlap_1d(lz_a, lz_b, p[2] - a[2], p[2] - b[2], gamma);

    let na = gaussian_norm(lx_a, ly_a, lz_a, alpha_a);
    let nb = gaussian_norm(lx_b, ly_b, lz_b, alpha_b);

    na * nb * k_ab * sx * sy * sz
}

/// Compute STO overlap between two basis functions using STO-nG expansion.
///
/// # Arguments
/// * `n_a, n_b` — principal quantum number (1, 2, 3, ...)
/// * `l_a, l_b` — angular momentum (0=s, 1=p, 2=d)
/// * `m_a, m_b` — m-component index (0-based)
/// * `zeta_a, zeta_b` — Slater exponents (bohr⁻¹)
/// * `pos_a, pos_b` — atomic positions (bohr)
pub fn sto_overlap_ng(
    n_a: u8,
    l_a: u8,
    m_a: u8,
    zeta_a: f64,
    pos_a: &[f64; 3],
    n_b: u8,
    l_b: u8,
    m_b: u8,
    zeta_b: f64,
    pos_b: &[f64; 3],
) -> f64 {
    // Default: STO-3G
    sto_overlap_with_ng(
        n_a, l_a, m_a, zeta_a, pos_a, 3, n_b, l_b, m_b, zeta_b, pos_b, 3,
    )
}

/// Compute STO overlap between two basis functions with explicit number of Gaussians.
///
/// `ng_a` and `ng_b` specify 3, 4, or 6 Gaussians per center.
/// For d-orbitals, uses proper spherical harmonic → Cartesian decomposition.
pub fn sto_overlap_with_ng(
    n_a: u8,
    l_a: u8,
    m_a: u8,
    zeta_a: f64,
    pos_a: &[f64; 3],
    ng_a: u8,
    n_b: u8,
    l_b: u8,
    m_b: u8,
    zeta_b: f64,
    pos_b: &[f64; 3],
    ng_b: u8,
) -> f64 {
    let comps_a = spherical_to_cartesian_components(l_a, m_a);
    let comps_b = spherical_to_cartesian_components(l_b, m_b);
    let zeta_a_sq = zeta_a * zeta_a;
    let zeta_b_sq = zeta_b * zeta_b;

    // Collect primitives for A
    let prims_a: Vec<(f64, f64)> = match ng_a {
        4 => sto4g_coefficients(n_a, l_a).to_vec(),
        6 => sto6g_coefficients(n_a, l_a).to_vec(),
        _ => sto3g_coefficients(n_a, l_a).to_vec(),
    };
    // Collect primitives for B
    let prims_b: Vec<(f64, f64)> = match ng_b {
        4 => sto4g_coefficients(n_b, l_b).to_vec(),
        6 => sto6g_coefficients(n_b, l_b).to_vec(),
        _ => sto3g_coefficients(n_b, l_b).to_vec(),
    };

    let mut s = 0.0;
    for &((lx_a, ly_a, lz_a), coeff_a) in comps_a.iter() {
        for &((lx_b, ly_b, lz_b), coeff_b) in comps_b.iter() {
            let mut s_cart = 0.0;
            for &(ar, ca) in prims_a.iter() {
                let alpha_a = ar * zeta_a_sq;
                for &(br, cb) in prims_b.iter() {
                    let alpha_b = br * zeta_b_sq;
                    s_cart += ca
                        * cb
                        * primitive_overlap(
                            lx_a, ly_a, lz_a, alpha_a, pos_a, lx_b, ly_b, lz_b, alpha_b, pos_b,
                        );
                }
            }
            s += coeff_a * coeff_b * s_cart;
        }
    }
    s
}

/// Look up STO-6G expansion for given principal quantum number and angular momentum.
fn sto6g_coefficients(n: u8, l: u8) -> Vec<(f64, f64)> {
    match (n, l) {
        (6, 0) => STO6G_6S.to_vec(),
        (6, 1) => STO6G_6P.to_vec(),
        // For n < 6 or d-shells, use STO-4G as a fallback (these won't normally be called with ng=6)
        _ => sto4g_coefficients(n, l).to_vec(),
    }
}

// ── Multipole integrals (dipole + quadrupole) ─────────────────────────────────

/// Closed-form 1D overlap moment integral: ∫ x^l exp(-α x²) dx
/// Returns (0.5/alpha)^(l/2) * (l-1)!! for even l, 0 for odd l.
fn s1d_moment(l: usize, eab: f64) -> f64 {
    if !l.is_multiple_of(2) {
        return 0.0;
    }
    let half_inv = 0.5 / eab;
    let mut val = (PI / eab).sqrt(); // base case l=0
    let half_l = l / 2;
    for i in 0..half_l {
        val *= half_inv * (2 * i + 1) as f64;
    }
    val
}

/// Shift polynomial coefficients from Gaussian product center P to center A.
/// coeffs[0..=l] are polynomial coefficients; ae = P_k - A_k.
/// After shift, coeffs represent the expansion about P instead of A.
fn horizontal_shift(ae: f64, l: usize, coeffs: &mut [f64]) {
    match l {
        0 => {}
        1 => {
            coeffs[0] += ae * coeffs[1];
        }
        2 => {
            coeffs[0] += ae * ae * coeffs[2];
            coeffs[1] += 2.0 * ae * coeffs[2];
        }
        3 => {
            coeffs[0] += ae * ae * ae * coeffs[3];
            coeffs[1] += 3.0 * ae * ae * coeffs[3];
            coeffs[2] += 3.0 * ae * coeffs[3];
        }
        4 => {
            coeffs[0] += ae * ae * ae * ae * coeffs[4];
            coeffs[1] += 4.0 * ae * ae * ae * coeffs[4];
            coeffs[2] += 6.0 * ae * ae * coeffs[4];
            coeffs[3] += 4.0 * ae * coeffs[4];
        }
        _ => {} // l > 4 not needed for s/p/d
    }
}

/// Form product of two polynomial coefficient arrays.
/// a[0..=la], b[0..=lb] → d[0..=la+lb]
fn form_product(a: &[f64], b: &[f64], la: usize, lb: usize, d: &mut [f64]) {
    for i in 0..=la + lb {
        d[i] = 0.0;
    }
    for i in 0..=la {
        for j in 0..=lb {
            d[i + j] += a[i] * b[j];
        }
    }
}

/// Result of multipole integral computation for one Cartesian primitive pair.
#[derive(Clone, Debug)]
pub struct MultipoleResult {
    pub overlap: f64,
    pub dipole: [f64; 3],
    pub quadrupole: [f64; 6],
}

/// Compute overlap, dipole, and quadrupole integrals between two Cartesian Gaussian primitives.
/// Dipole is about center A (the first atom). Quadrupole is traceless about center A.
fn primitive_multipole(
    lx_a: u8,
    ly_a: u8,
    lz_a: u8,
    alpha_a: f64,
    a: &[f64; 3],
    lx_b: u8,
    ly_b: u8,
    lz_b: u8,
    alpha_b: f64,
    b: &[f64; 3],
) -> MultipoleResult {
    let gamma = alpha_a + alpha_b;
    let mu = alpha_a * alpha_b / gamma;
    let r2 = (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2);
    let pre = (-mu * r2).exp() * (PI / gamma).sqrt().powi(3);

    let p = [
        (alpha_a * a[0] + alpha_b * b[0]) / gamma,
        (alpha_a * a[1] + alpha_b * b[1]) / gamma,
        (alpha_a * a[2] + alpha_b * b[2]) / gamma,
    ];
    let rpi = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let rpj = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];

    let la = [lx_a as usize, ly_a as usize, lz_a as usize];
    let lb = [lx_b as usize, ly_b as usize, lz_b as usize];

    // Maximum angular momentum sum + 2 for quadrupole
    let max_l = la[0] + lb[0] + la[1] + lb[1] + la[2] + lb[2] + 4;
    let mut s1d_arr: Vec<f64> = (0..=max_l).map(|l| s1d_moment(l, gamma)).collect();
    // Normalize: s1d already includes sqrt(pi/gamma) from s1d_moment for l=0.
    // But s1d_moment gives full 1D integral; we need to factor out sqrt(pi/gamma)
    // since pre already includes (pi/gamma)^{3/2}.
    // Actually, tblite's `pre` = exp(-est) * sqrtpi3 * sqrt(oab)^3
    // and their s1d(l) = overlap_1d(l, eab) which is (0.5/eab)^(l/2)*(l-1)!!
    // So s1d(0) = 1.0 in tblite and the pi/gamma factor is in `pre`.
    // Let me normalize: divide each s1d by sqrt(pi/gamma) to get the tblite convention.
    let sqrt_pi_gamma = (PI / gamma).sqrt();
    for s in s1d_arr.iter_mut() {
        *s /= sqrt_pi_gamma;
    }

    let na = gaussian_norm(lx_a, ly_a, lz_a, alpha_a);
    let nb = gaussian_norm(lx_b, ly_b, lz_b, alpha_b);
    let cc = na * nb * pre;

    // For each Cartesian direction, compute v1d[k][0..3] = (overlap, dipole, quadrupole)
    let mut v1d = [[0.0f64; 3]; 3];

    for k in 0..3 {
        // Build polynomial expansions for center i and j
        let mut vi = vec![0.0f64; la[k] + lb[k] + 5];
        let mut vj = vec![0.0f64; la[k] + lb[k] + 5];
        vi[la[k]] = 1.0;
        vj[lb[k]] = 1.0;

        horizontal_shift(rpi[k], la[k], &mut vi);
        horizontal_shift(rpj[k], lb[k], &mut vj);

        let mut vv = vec![0.0f64; la[k] + lb[k] + 5];
        form_product(&vi, &vj, la[k], lb[k], &mut vv);

        for l in 0..=la[k] + lb[k] {
            v1d[k][0] += s1d_arr[l] * vv[l];
            v1d[k][1] += (s1d_arr[l + 1] + rpi[k] * s1d_arr[l]) * vv[l];
            v1d[k][2] +=
                (s1d_arr[l + 2] + 2.0 * rpi[k] * s1d_arr[l + 1] + rpi[k] * rpi[k] * s1d_arr[l])
                    * vv[l];
        }
    }

    let s3d = v1d[0][0] * v1d[1][0] * v1d[2][0];

    let d3d = [
        v1d[0][1] * v1d[1][0] * v1d[2][0], // x
        v1d[0][0] * v1d[1][1] * v1d[2][0], // y
        v1d[0][0] * v1d[1][0] * v1d[2][1], // z
    ];

    let q3d = [
        v1d[0][2] * v1d[1][0] * v1d[2][0], // xx
        v1d[0][1] * v1d[1][1] * v1d[2][0], // xy
        v1d[0][0] * v1d[1][2] * v1d[2][0], // yy
        v1d[0][1] * v1d[1][0] * v1d[2][1], // xz
        v1d[0][0] * v1d[1][1] * v1d[2][1], // yz
        v1d[0][0] * v1d[1][0] * v1d[2][2], // zz
    ];

    MultipoleResult {
        overlap: cc * s3d,
        dipole: [cc * d3d[0], cc * d3d[1], cc * d3d[2]],
        quadrupole: [
            cc * q3d[0],
            cc * q3d[1],
            cc * q3d[2],
            cc * q3d[3],
            cc * q3d[4],
            cc * q3d[5],
        ],
    }
}

/// Compute STO-NG multipole integrals (overlap, dipole[3], traceless quadrupole[6])
/// between two basis functions.
///
/// Dipole is about center A (pos_a). Quadrupole is traceless.
/// For d-orbitals, uses proper spherical harmonic → Cartesian decomposition.
/// Returns (overlap, dipole[3], quadrupole[6]).
pub fn sto_multipole_with_ng(
    n_a: u8,
    l_a: u8,
    m_a: u8,
    zeta_a: f64,
    pos_a: &[f64; 3],
    ng_a: u8,
    n_b: u8,
    l_b: u8,
    m_b: u8,
    zeta_b: f64,
    pos_b: &[f64; 3],
    ng_b: u8,
) -> (f64, [f64; 3], [f64; 6]) {
    let comps_a = spherical_to_cartesian_components(l_a, m_a);
    let comps_b = spherical_to_cartesian_components(l_b, m_b);
    let zeta_a_sq = zeta_a * zeta_a;
    let zeta_b_sq = zeta_b * zeta_b;

    let prims_a: Vec<(f64, f64)> = match ng_a {
        4 => sto4g_coefficients(n_a, l_a).to_vec(),
        6 => sto6g_coefficients(n_a, l_a).to_vec(),
        _ => sto3g_coefficients(n_a, l_a).to_vec(),
    };
    let prims_b: Vec<(f64, f64)> = match ng_b {
        4 => sto4g_coefficients(n_b, l_b).to_vec(),
        6 => sto6g_coefficients(n_b, l_b).to_vec(),
        _ => sto3g_coefficients(n_b, l_b).to_vec(),
    };

    let mut s_total = 0.0;
    let mut d_total = [0.0f64; 3];
    let mut q_total = [0.0f64; 6];

    for &((lx_a, ly_a, lz_a), coeff_sph_a) in comps_a.iter() {
        for &((lx_b, ly_b, lz_b), coeff_sph_b) in comps_b.iter() {
            let cc_sph = coeff_sph_a * coeff_sph_b;
            let mut s_cart = 0.0;
            let mut d_cart = [0.0f64; 3];
            let mut q_cart = [0.0f64; 6];

            for &(ar, ca) in prims_a.iter() {
                let alpha_a = ar * zeta_a_sq;
                for &(br, cb) in prims_b.iter() {
                    let alpha_b = br * zeta_b_sq;
                    let mp = primitive_multipole(
                        lx_a, ly_a, lz_a, alpha_a, pos_a, lx_b, ly_b, lz_b, alpha_b, pos_b,
                    );
                    let c = ca * cb;
                    s_cart += c * mp.overlap;
                    for i in 0..3 {
                        d_cart[i] += c * mp.dipole[i];
                    }
                    for i in 0..6 {
                        q_cart[i] += c * mp.quadrupole[i];
                    }
                }
            }

            s_total += cc_sph * s_cart;
            for i in 0..3 {
                d_total[i] += cc_sph * d_cart[i];
            }
            for i in 0..6 {
                q_total[i] += cc_sph * q_cart[i];
            }
        }
    }

    // Make quadrupole traceless: q' = 1.5*q - 0.5*tr*I
    let tr = 0.5 * (q_total[0] + q_total[2] + q_total[5]); // 0.5*(xx + yy + zz)
    q_total[0] = 1.5 * q_total[0] - tr; // xx
    q_total[1] *= 1.5; // xy
    q_total[2] = 1.5 * q_total[2] - tr; // yy
    q_total[3] *= 1.5; // xz
    q_total[4] *= 1.5; // yz
    q_total[5] = 1.5 * q_total[5] - tr; // zz

    (s_total, d_total, q_total)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_1s_1s_overlap_same_center() {
        // Same atom → overlap = 1.0
        let pos = [0.0, 0.0, 0.0];
        let s = sto_overlap_ng(1, 0, 0, 1.23, &pos, 1, 0, 0, 1.23, &pos);
        assert!(
            (s - 1.0).abs() < 0.02,
            "1s-1s self-overlap should be ~1.0, got {s}"
        );
    }

    #[test]
    fn test_1s_1s_overlap_h2() {
        // H-H at 1.398 bohr with ζ=1.23
        // Exact: exp(-1.720)(1+1.720+0.987) = 0.663
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 1.398];
        let s = sto_overlap_ng(1, 0, 0, 1.23, &a, 1, 0, 0, 1.23, &b);
        assert!(
            (s - 0.663).abs() < 0.02,
            "H-H 1s-1s overlap should be ~0.663, got {s}"
        );
    }

    #[test]
    fn test_2s_2s_overlap_cc() {
        // C-C at 2.91 bohr (1.54 Å) with ζ_s=2.096
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 2.91];
        let s = sto_overlap_ng(2, 0, 0, 2.096, &a, 2, 0, 0, 2.096, &b);
        // 2s-2s should be significantly larger than the 1s formula gives (0.044)
        assert!(s > 0.10, "2s-2s overlap should be > 0.10, got {s}");
        assert!(s < 0.50, "2s-2s overlap should be < 0.50, got {s}");
    }

    #[test]
    fn test_2p_2p_sigma_overlap() {
        // Two pz orbitals (both m=2→pz) along the z-axis (σ bond)
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 2.91];
        let s = sto_overlap_ng(2, 1, 2, 1.80, &a, 2, 1, 2, 1.80, &b);
        // pz-pz σ overlap should be significant
        assert!(
            s.abs() > 0.01,
            "2pσ-2pσ overlap should be significant, got {s}"
        );
    }

    #[test]
    fn test_2p_2p_pi_overlap() {
        // Two px orbitals (both m=0→px) along the z-axis (π bond)
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 2.91];
        let s = sto_overlap_ng(2, 1, 0, 1.80, &a, 2, 1, 0, 1.80, &b);
        // px-px π overlap should be smaller than σ
        assert!(
            s.abs() > 0.001,
            "2pπ-2pπ overlap should be nonzero, got {s}"
        );
    }

    #[test]
    fn test_sp_orthogonality() {
        // 2s and 2px on the same atom should be orthogonal
        let pos = [0.0, 0.0, 0.0];
        let s = sto_overlap_ng(2, 0, 0, 2.096, &pos, 2, 1, 0, 1.80, &pos);
        assert!(s.abs() < 0.01, "s-p on same center should be ~0, got {s}");
    }

    #[test]
    fn test_multipole_self_overlap() {
        // 1s on same center: overlap = 1 (approx with STO-3G), dipole = 0, quadrupole = 0
        let pos = [0.0, 0.0, 0.0];
        let (s, d, q) = sto_multipole_with_ng(1, 0, 0, 1.23, &pos, 3, 1, 0, 0, 1.23, &pos, 3);
        assert!((s - 1.0).abs() < 0.02, "1s self-overlap ~1.0, got {s}");
        assert!(d[0].abs() < 1e-10, "dipole x should be 0 on same center");
        assert!(d[1].abs() < 1e-10, "dipole y should be 0 on same center");
        assert!(d[2].abs() < 1e-10, "dipole z should be 0 on same center");
        // Traceless quadrupole on same center for s orbital should be zero
        let q_trace = q[0] + q[2] + q[5];
        assert!(
            q_trace.abs() < 1e-10,
            "traceless quadrupole trace should be 0, got {q_trace}"
        );
    }

    #[test]
    fn test_multipole_h2_dipole() {
        // H2 along z: 1s-1s at 1.398 bohr separation.
        // Dipole integral about atom A should have nonzero z component, zero x/y.
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 1.398];
        let (s, d, _q) = sto_multipole_with_ng(1, 0, 0, 1.23, &a, 3, 1, 0, 0, 1.23, &b, 3);
        assert!(s > 0.5, "overlap should be substantial, got {s}");
        assert!(d[0].abs() < 1e-10, "dipole x should be 0 by symmetry");
        assert!(d[1].abs() < 1e-10, "dipole y should be 0 by symmetry");
        assert!(
            d[2].abs() > 0.01,
            "dipole z should be nonzero, got {}",
            d[2]
        );
        // Dipole z should be positive (center of charge between A and B, relative to A)
        assert!(
            d[2] > 0.0,
            "dipole z about A should be positive for B at +z"
        );
    }

    #[test]
    fn test_multipole_overlap_matches_sto_overlap() {
        // The overlap from sto_multipole_with_ng should match sto_overlap_with_ng
        let a = [0.0, 0.0, 0.0];
        let b = [1.5, 0.8, 2.1];
        for ng in [3u8, 4, 6] {
            let s_ref = sto_overlap_with_ng(2, 0, 0, 2.096, &a, ng, 2, 0, 0, 1.80, &b, ng);
            let (s_mp, _, _) = sto_multipole_with_ng(2, 0, 0, 2.096, &a, ng, 2, 0, 0, 1.80, &b, ng);
            assert!(
                (s_ref - s_mp).abs() < 1e-8,
                "ng={ng}: overlap mismatch: ref={s_ref}, multipole={s_mp}"
            );
        }
    }

    #[test]
    fn test_multipole_traceless_quadrupole() {
        // Quadrupole should be traceless for any pair
        let a = [0.0, 0.0, 0.0];
        let b = [1.2, 0.5, 0.8];
        let (_s, _d, q) = sto_multipole_with_ng(2, 1, 0, 1.80, &a, 4, 2, 0, 0, 2.0, &b, 4);
        let trace = q[0] + q[2] + q[5]; // xx + yy + zz
        assert!(
            trace.abs() < 1e-10,
            "quadrupole trace should be 0, got {trace}"
        );
    }

    #[test]
    fn test_multipole_vs_overlap_consistency() {
        // Compare primitive_multipole overlap with primitive_overlap for s-s
        let a = [0.0, 0.0, 0.0];
        let b = [0.0, 0.0, 2.672];
        let alpha_a = 2.227660 * 2.363_f64.powi(2);
        let alpha_b = 3.242212833e-1 * 2.026128_f64.powi(2);

        let s_direct = primitive_overlap(0, 0, 0, alpha_a, &a, 0, 0, 0, alpha_b, &b);
        let mp = primitive_multipole(0, 0, 0, alpha_a, &a, 0, 0, 0, alpha_b, &b);

        println!("primitive_overlap:  {:.15e}", s_direct);
        println!("multipole.overlap:  {:.15e}", mp.overlap);
        println!("ratio: {:.15}", mp.overlap / s_direct);

        // Now compare full STO-nG overlap vs multipole
        let (s_mp, _, _) = sto_multipole_with_ng(1, 0, 0, 2.363, &a, 3, 4, 0, 0, 2.026128, &b, 4);
        let s_ov = sto_overlap_with_ng(1, 0, 0, 2.363, &a, 3, 4, 0, 0, 2.026128, &b, 4);
        println!("\nFull STO-nG:");
        println!("sto_overlap_with_ng:    {:.15e}", s_ov);
        println!("sto_multipole_with_ng:  {:.15e}", s_mp);
        println!("ratio: {:.15}", s_mp / s_ov);

        assert!(
            (s_direct - mp.overlap).abs() < 1e-12,
            "primitive overlap vs multipole differ: {} vs {}",
            s_direct,
            mp.overlap
        );
    }

    #[test]
    fn test_hbr_overlap_debug() {
        // HBr along z-axis: H at origin, Br at z=2.672 bohr (~1.414 Å)
        // GFN2 params for H: pqn=1, l=0, zeta=2.363, STO-4G
        // GFN2 params for Br:
        //   4s: pqn=4, l=0, zeta=2.026128, STO-4G
        //   4p: pqn=4, l=1, zeta=1.949257, STO-4G
        //   4d: pqn=4, l=2, zeta=1.040181, STO-3G
        let h_pos = [0.0, 0.0, 0.0];
        let br_pos = [0.0, 0.0, 2.672];

        println!("\n=== HBr overlap matrix (GFN2 params) ===");

        // H 1s - Br 4s
        let s_hs = sto_overlap_with_ng(1, 0, 0, 2.363, &h_pos, 4, 4, 0, 0, 2.026128, &br_pos, 4);
        println!("S(H_1s, Br_4s) = {:.8e}", s_hs);

        // H 1s - Br 4p (m=0 px, m=1 py, m=2 pz)
        for m in 0..3u8 {
            let s = sto_overlap_with_ng(1, 0, 0, 2.363, &h_pos, 4, 4, 1, m, 1.949257, &br_pos, 4);
            println!("S(H_1s, Br_4p_m{}) = {:.8e}", m, s);
        }

        // H 1s - Br 4d (m=0..4: dxy, dxz, dyz, dx2-y2, dz2)
        for m in 0..5u8 {
            let name = match m {
                0 => "dxy",
                1 => "dxz",
                2 => "dyz",
                3 => "dx2-y2",
                _ => "dz2",
            };
            let s = sto_overlap_with_ng(1, 0, 0, 2.363, &h_pos, 4, 4, 2, m, 1.040181, &br_pos, 3);
            println!("S(H_1s, Br_{}) = {:.8e}", name, s);
        }

        // Now compute the raw Cartesian d overlaps to understand the decomposition
        println!("\n=== Raw Cartesian d primitive overlaps ===");
        let zeta_h = 2.363;
        let zeta_d = 1.040181;

        let cart_labels = ["xx", "yy", "zz", "xy", "xz", "yz"];
        let cart_lxyz: [(u8, u8, u8); 6] = [
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1),
        ];

        for (idx, &(lx, ly, lz)) in cart_lxyz.iter().enumerate() {
            let mut s = 0.0;
            // H: STO-4G, 1s
            let prims_h = [
                (5.216844534e+0, 5.675242080e-2),
                (9.546182760e-1, 2.601413550e-1),
                (2.652034102e-1, 5.324461400e-1),
                (8.801862774e-2, 2.916254405e-1),
            ];
            // Br d: STO-3G, 4d
            let prims_d = [
                (0.338099, 0.023893),
                (0.085930, 0.505389),
                (0.038380, 0.562553),
            ];

            for &(ar, ca) in &prims_h {
                let alpha_a = ar * zeta_h * zeta_h;
                for &(br, cb) in &prims_d {
                    let alpha_b = br * zeta_d * zeta_d;
                    s += ca
                        * cb
                        * primitive_overlap(0, 0, 0, alpha_a, &h_pos, lx, ly, lz, alpha_b, &br_pos);
                }
            }
            println!("  <1s_H | d_{}_Br> = {:.8e}", cart_labels[idx], s);
        }

        // Verify: dz2 = -0.5*<s|xx> -0.5*<s|yy> + 1.0*<s|zz>
        // vs old: dz2 = <s|zz>
        println!("\n=== Decomposition check ===");
        let comps = spherical_to_cartesian_components(2, 4);
        println!("dz2 components: {:?}", comps);
        let comps2 = spherical_to_cartesian_components(2, 3);
        println!("dx2-y2 components: {:?}", comps2);

        // Br d-d self-overlaps (should be ~1 for same m)
        println!("\n=== Br 4d self-overlaps ===");
        for m in 0..5u8 {
            let s =
                sto_overlap_with_ng(4, 2, m, 1.040181, &br_pos, 3, 4, 2, m, 1.040181, &br_pos, 3);
            println!("S(Br_4d_m{}, Br_4d_m{}) = {:.8}", m, m, s);
        }
    }

    #[test]
    fn test_water_overlap_vs_tblite() {
        // Water geometry (Å → bohr)
        let aatoau = 1.8897259886;
        let o_pos = [0.0 * aatoau, 0.0 * aatoau, 0.117300 * aatoau];
        let h1_pos = [0.0 * aatoau, 0.757200 * aatoau, -0.469200 * aatoau];

        let zeta_o_s = 2.439742;
        let zeta_h = 1.23;

        // STO-4G for O 2s (n=2, l=0)
        let sto4g_2s = sto4g_coefficients(2, 0);
        // STO-3G for H 1s (n=1, l=0)
        let sto3g_1s = sto3g_coefficients(1, 0);

        println!("\nO pos bohr: {:?}", o_pos);
        println!("H1 pos bohr: {:?}", h1_pos);
        println!("zeta_O_s = {}, zeta_H = {}", zeta_o_s, zeta_h);
        println!("\nSTO-4G(2s) coefficients:");
        for (a, c) in sto4g_2s {
            println!("  alpha_raw={:.12e}, coeff={:.12e}", a, c);
        }
        println!("\nSTO-3G(1s) coefficients:");
        for (a, c) in sto3g_1s {
            println!("  alpha_raw={:.12e}, coeff={:.12e}", a, c);
        }

        let zeta_o_sq = zeta_o_s * zeta_o_s;
        let zeta_h_sq = zeta_h * zeta_h;

        let mut total = 0.0;
        for (i, &(ar_o, c_o)) in sto4g_2s.iter().enumerate() {
            let alpha_o = ar_o * zeta_o_sq;
            for (j, &(ar_h, c_h)) in sto3g_1s.iter().enumerate() {
                let alpha_h = ar_h * zeta_h_sq;
                let prim = primitive_overlap(0, 0, 0, alpha_o, &o_pos, 0, 0, 0, alpha_h, &h1_pos);
                total += c_o * c_h * prim;
                println!(
                    "prim[{},{}] alpha_o={:.8e} alpha_h={:.8e} overlap={:.15e}",
                    i, j, alpha_o, alpha_h, prim
                );
            }
        }
        println!("\nTotal primitive_overlap sum: {:.15e}", total);

        // Also compute via sto_overlap_with_ng
        let s_func = sto_overlap_with_ng(2, 0, 0, zeta_o_s, &o_pos, 4, 1, 0, 0, zeta_h, &h1_pos, 3);
        println!("sto_overlap_with_ng:       {:.15e}", s_func);

        // Also compute via sto_multipole_with_ng
        let (s_mp, _, _) =
            sto_multipole_with_ng(2, 0, 0, zeta_o_s, &o_pos, 4, 1, 0, 0, zeta_h, &h1_pos, 3);
        println!("sto_multipole_with_ng:     {:.15e}", s_mp);

        // Expected from tblite (Python): 4.446492247020593e-01
        let expected = 4.446492247020593e-01;
        println!("\nExpected (tblite/Python):   {:.15e}", expected);
        println!("Diff primitive_overlap:     {:.3e}", total - expected);
        println!("Diff sto_overlap_with_ng:   {:.3e}", s_func - expected);
        println!("Diff sto_multipole_with_ng: {:.3e}", s_mp - expected);
    }
}
