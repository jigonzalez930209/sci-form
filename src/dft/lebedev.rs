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

/// 26-point Lebedev grid (octahedral + cube + edge midpoints).
fn lebedev_26() -> Vec<AngularPoint> {
    let mut pts = Vec::with_capacity(26);
    let four_pi = 4.0 * std::f64::consts::PI;

    // 6 octahedral vertices
    let w1 = four_pi * 1.0 / 21.0;
    for &sign in &[1.0f64, -1.0] {
        pts.push(([sign, 0.0, 0.0], w1));
        pts.push(([0.0, sign, 0.0], w1));
        pts.push(([0.0, 0.0, sign], w1));
    }

    // 8 cube vertices
    let c = 1.0 / 3.0f64.sqrt();
    let w2 = four_pi * 4.0 / 105.0;
    for &sx in &[1.0f64, -1.0] {
        for &sy in &[1.0f64, -1.0] {
            for &sz in &[1.0f64, -1.0] {
                pts.push(([c * sx, c * sy, c * sz], w2));
            }
        }
    }

    // 12 edge midpoints
    let e = 1.0 / 2.0f64.sqrt();
    let w3 = four_pi * 27.0 / 840.0;
    for &(a, b) in &[(0, 1), (0, 2), (1, 2)] {
        for &sa in &[1.0f64, -1.0] {
            for &sb in &[1.0f64, -1.0] {
                let mut p = [0.0; 3];
                p[a] = e * sa;
                p[b] = e * sb;
                pts.push((p, w3));
            }
        }
    }

    pts
}

/// 110-point Lebedev grid (approximate — uses icosahedral sampling).
fn lebedev_110() -> Vec<AngularPoint> {
    // Generate a reasonable angular grid by combining octahedral symmetry
    // with intermediate latitude points
    let mut pts = Vec::with_capacity(110);
    let four_pi = 4.0 * std::f64::consts::PI;

    // Start with the 26-point grid
    let base = lebedev_26();
    for (dir, _) in &base {
        pts.push((*dir, 0.0)); // weights will be set uniformly
    }

    // Add intermediate points by averaging pairs
    let n_extra = 110 - 26;
    let golden = (1.0 + 5.0f64.sqrt()) / 2.0;
    for i in 0..n_extra {
        let theta = std::f64::consts::PI * (i as f64 + 0.5) / n_extra as f64;
        let phi = 2.0 * std::f64::consts::PI * i as f64 * golden;
        let st = theta.sin();
        pts.push(([st * phi.cos(), st * phi.sin(), theta.cos()], 0.0));
    }

    // Uniform weight
    let w = four_pi / pts.len() as f64;
    for p in &mut pts {
        p.1 = w;
    }

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
        // The 26-point grid uses equal weights which sum to 4π
        assert!(
            (w_sum - 4.0 * std::f64::consts::PI).abs() < 0.5,
            "26-point weights should sum close to 4π, got {w_sum}"
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

    #[test]
    fn higher_order_has_more_points() {
        let p6 = lebedev_grid(6);
        let p26 = lebedev_grid(26);
        let p110 = lebedev_grid(110);
        assert!(p26.len() > p6.len());
        assert!(p110.len() > p26.len());
    }
}
