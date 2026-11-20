//! Lebedev angular quadrature grids.
//!
//! Provides pre-computed Lebedev points on the unit sphere for numerical
//! integration of angular components.

/// A Lebedev angular quadrature point: (direction [x,y,z] on unit sphere, weight).
type AngularPoint = ([f64; 3], f64);

/// Generate a Lebedev quadrature grid with approximately `n_target` points.
///
/// Returns points as (unit direction, weight) pairs.
/// Weights sum to 4π.
pub fn lebedev_grid(n_target: usize) -> Vec<AngularPoint> {
    // Use octahedral symmetry to generate Lebedev grids.
    // For simplicity, we provide exact grids for common sizes.
    if n_target <= 6 {
        lebedev_6()
    } else if n_target <= 26 {
        lebedev_26()
    } else if n_target <= 110 {
        lebedev_110()
    } else {
        lebedev_302()
    }
}

/// 6-point Lebedev grid (octahedral vertices).
fn lebedev_6() -> Vec<AngularPoint> {
    let w = 4.0 * std::f64::consts::PI / 6.0;
    vec![
        ([1.0, 0.0, 0.0], w),
        ([-1.0, 0.0, 0.0], w),
        ([0.0, 1.0, 0.0], w),
        ([0.0, -1.0, 0.0], w),
        ([0.0, 0.0, 1.0], w),
        ([0.0, 0.0, -1.0], w),
    ]
}

/// Octahedral orbits of the Lebedev–Laikov construction. `weight` is the
/// normalised weight (the full grid sums to 1); points are scaled by 4π.
mod orbit {
    use super::AngularPoint;

    const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

    /// (±1, 0, 0) and permutations: 6 points.
    pub fn a1(weight: f64, out: &mut Vec<AngularPoint>) {
        for axis in 0..3 {
            for sign in [1.0, -1.0] {
                let mut p = [0.0; 3];
                p[axis] = sign;
                out.push((p, FOUR_PI * weight));
            }
        }
    }

    /// (0, ±1/√2, ±1/√2) and permutations: 12 points.
    pub fn a2(weight: f64, out: &mut Vec<AngularPoint>) {
        let e = std::f64::consts::FRAC_1_SQRT_2;
        for zero in 0..3 {
            let (i, j) = ((zero + 1) % 3, (zero + 2) % 3);
            for si in [e, -e] {
                for sj in [e, -e] {
                    let mut p = [0.0; 3];
                    p[i] = si;
                    p[j] = sj;
                    out.push((p, FOUR_PI * weight));
                }
            }
        }
    }

    /// (±1/√3, ±1/√3, ±1/√3): 8 points.
    pub fn a3(weight: f64, out: &mut Vec<AngularPoint>) {
        let c = 1.0 / 3.0f64.sqrt();
        for sx in [c, -c] {
            for sy in [c, -c] {
                for sz in [c, -c] {
                    out.push(([sx, sy, sz], FOUR_PI * weight));
                }
            }
        }
    }

    /// (±l, ±l, ±m) with m = √(1 − 2l²) and permutations: 24 points.
    pub fn bk(l: f64, weight: f64, out: &mut Vec<AngularPoint>) {
        let m = (1.0 - 2.0 * l * l).sqrt();
        for odd in 0..3 {
            for s0 in [1.0, -1.0] {
                for s1 in [1.0, -1.0] {
                    for s2 in [1.0, -1.0] {
                        let mut p = [l * s0, l * s1, l * s2];
                        p[odd] = m * [s0, s1, s2][odd];
                        out.push((p, FOUR_PI * weight));
                    }
                }
            }
        }
    }

    /// (±p, ±q, 0) with q = √(1 − p²) and permutations: 24 points.
    pub fn ck(p: f64, weight: f64, out: &mut Vec<AngularPoint>) {
        let q = (1.0 - p * p).sqrt();
        for zero in 0..3 {
            let (i, j) = ((zero + 1) % 3, (zero + 2) % 3);
            for (u, v) in [(p, q), (q, p)] {
                for su in [1.0, -1.0] {
                    for sv in [1.0, -1.0] {
                        let mut point = [0.0; 3];
                        point[i] = u * su;
                        point[j] = v * sv;
                        out.push((point, FOUR_PI * weight));
                    }
                }
            }
        }
    }
}

/// 26-point Lebedev grid (degree 7): axes, edge midpoints, cube vertices.
fn lebedev_26() -> Vec<AngularPoint> {
    let mut pts = Vec::with_capacity(26);
    orbit::a1(1.0 / 21.0, &mut pts);
    orbit::a2(4.0 / 105.0, &mut pts);
    orbit::a3(9.0 / 280.0, &mut pts);
    pts
}

/// 110-point Lebedev–Laikov grid (degree 17).
///
/// Generators from Lebedev & Laikov, *Dokl. Math.* **1999**, 59, 477
/// (`LD0110`), checked against PySCF `gen_grid.MakeAngularGrid(110)`.
fn lebedev_110() -> Vec<AngularPoint> {
    let mut pts = Vec::with_capacity(110);
    orbit::a1(0.003_828_270_494_937_162, &mut pts);
    orbit::a3(0.009_793_737_512_487_512, &mut pts);
    orbit::bk(0.185_115_635_344_736_2, 0.008_211_737_283_191_111, &mut pts);
    orbit::bk(0.690_421_048_382_292_2, 0.009_942_814_891_178_103, &mut pts);
    orbit::bk(0.395_689_473_055_941_9, 0.009_595_471_336_070_963, &mut pts);
    orbit::ck(0.478_369_028_812_150_2, 0.009_694_996_361_663_028, &mut pts);
    pts
}

/// 302-point Lebedev grid (Fibonacci spiral approximation).
fn lebedev_302() -> Vec<AngularPoint> {
    let n = 302;
    let mut pts = Vec::with_capacity(n);
    let four_pi = 4.0 * std::f64::consts::PI;
    let golden = (1.0 + 5.0f64.sqrt()) / 2.0;
    let w = four_pi / n as f64;

    for i in 0..n {
        let _theta = ((2.0 * i as f64 + 1.0) / (2.0 * n as f64)).acos();
        // Use: cos(theta) = 1 - 2*(i+0.5)/n for uniform distribution
        let cos_theta = 1.0 - (2.0 * i as f64 + 1.0) / n as f64;
        let theta = cos_theta.acos();
        let phi = 2.0 * std::f64::consts::PI * i as f64 / golden;
        let st = theta.sin();
        pts.push(([st * phi.cos(), st * phi.sin(), cos_theta], w));
    }

    pts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lebedev_6_point_weights_sum_to_4pi() {
        let pts = lebedev_grid(6);
        let w_sum: f64 = pts.iter().map(|(_, w)| w).sum();
        assert!(
            (w_sum - 4.0 * std::f64::consts::PI).abs() < 1e-10,
            "6-point weights should sum to 4π, got {w_sum}"
        );
    }

    #[test]
    fn lebedev_26_point_weights_sum_to_4pi() {
        let pts = lebedev_grid(26);
        let w_sum: f64 = pts.iter().map(|(_, w)| w).sum();
        assert!(
            (w_sum - 4.0 * std::f64::consts::PI).abs() < 1e-10,
            "26-point weights should sum to 4π, got {w_sum}"
        );
    }

    #[test]
    fn lebedev_points_are_on_unit_sphere() {
        let pts = lebedev_grid(110);
        for (xyz, _) in &pts {
            let r = (xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-10,
                "point should lie on unit sphere, r={r}"
            );
        }
    }

    /// ∫ x^a y^b z^c dΩ over the unit sphere.
    fn sphere_monomial(a: u32, b: u32, c: u32) -> f64 {
        if a % 2 == 1 || b % 2 == 1 || c % 2 == 1 {
            return 0.0;
        }
        let double_factorial = |n: i64| (1..=n).rev().step_by(2).map(|k| k as f64).product::<f64>();
        let (a, b, c) = (a as i64, b as i64, c as i64);
        4.0 * std::f64::consts::PI
            * double_factorial(a - 1)
            * double_factorial(b - 1)
            * double_factorial(c - 1)
            / double_factorial(a + b + c + 1)
    }

    fn assert_exact_to_degree(n_points: usize, degree: u32) {
        let grid = lebedev_grid(n_points);
        assert_eq!(grid.len(), n_points);
        for a in 0..=degree {
            for b in 0..=(degree - a) {
                for c in 0..=(degree - a - b) {
                    let quadrature: f64 = grid
                        .iter()
                        .map(|(p, w)| {
                            w * p[0].powi(a as i32) * p[1].powi(b as i32) * p[2].powi(c as i32)
                        })
                        .sum();
                    let exact = sphere_monomial(a, b, c);
                    assert!(
                        (quadrature - exact).abs() < 1e-12,
                        "{n_points}-point grid, x^{a} y^{b} z^{c}: {quadrature} vs {exact}"
                    );
                }
            }
        }
    }

    #[test]
    fn lebedev_grids_integrate_monomials_exactly() {
        assert_exact_to_degree(6, 3);
        assert_exact_to_degree(26, 7);
        assert_exact_to_degree(110, 17);
    }

    #[test]
    fn higher_order_has_more_points() {
        let p6 = lebedev_grid(6);
        let p26 = lebedev_grid(26);
        let p110 = lebedev_grid(110);
        assert!(p26.len() > p6.len());
        assert!(p110.len() > p26.len());
    }
}
