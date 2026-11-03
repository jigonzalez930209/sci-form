//! Central-difference nuclear gradients of a scalar energy function.
//!
//! The semi-empirical methods (PM3, GFN-xTB) expose a single energy surface
//! through their public `solve_*` entry points. The gradient of that surface
//! is obtained here by central differences of that very function, so the
//! gradient is consistent with the reported energy by construction. (The
//! earlier hand-written Hellmann-Feynman expressions differentiated a
//! *different*, simplified model and disagreed with the energy by tens to
//! hundreds of eV/Å.)
//!
//! Cost: `6·N` energy evaluations (`2` per Cartesian coordinate); the
//! evaluations are independent and run on the rayon pool when the `parallel`
//! feature is on. The result does not depend on the thread count, because
//! each coordinate is computed independently and collected in order.

/// Default displacement (Å). Truncation error is O(h²·E‴) ≈ 5e-6 eV/Å and the
/// SCF-noise contribution is ≈ ε_SCF / 2h, well below the 1e-4 eV/Å level.
pub(crate) const DEFAULT_STEP_ANGSTROM: f64 = 5e-4;

/// Gradient `dE/dR` (energy units per Å) of `energy` at `positions`.
///
/// # Errors
///
/// Returns the first error produced by `energy` (e.g. a non-converged SCF at
/// a displaced geometry).
pub(crate) fn central_difference_gradient<F>(
    positions: &[[f64; 3]],
    step: f64,
    energy: F,
) -> Result<Vec<[f64; 3]>, String>
where
    F: Fn(&[[f64; 3]]) -> Result<f64, String> + Sync,
{
    let n_coords = positions.len() * 3;

    let one = |k: usize| -> Result<f64, String> {
        let (atom, dim) = (k / 3, k % 3);
        let mut plus = positions.to_vec();
        let mut minus = positions.to_vec();
        plus[atom][dim] += step;
        minus[atom][dim] -= step;
        Ok((energy(&plus)? - energy(&minus)?) / (2.0 * step))
    };

    #[cfg(feature = "parallel")]
    let flat: Vec<Result<f64, String>> = {
        use rayon::prelude::*;
        (0..n_coords).into_par_iter().map(one).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let flat: Vec<Result<f64, String>> = (0..n_coords).map(one).collect();

    let mut grad = vec![[0.0f64; 3]; positions.len()];
    for (k, g) in flat.into_iter().enumerate() {
        grad[k / 3][k % 3] = g?;
    }
    Ok(grad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_the_analytic_gradient_of_a_quartic() {
        // E = Σ x⁴ + x·y, so dE/dx = 4x³ + y, dE/dy = x, dE/dz = 4z³.
        let pos = [[0.5, -0.3, 1.1], [0.2, 0.9, -0.4]];
        let e = |p: &[[f64; 3]]| -> Result<f64, String> {
            Ok(p.iter()
                .map(|a| a[0].powi(4) + a[1].powi(4) + a[2].powi(4) + a[0] * a[1])
                .sum())
        };
        let g = central_difference_gradient(&pos, DEFAULT_STEP_ANGSTROM, e).unwrap();
        for (a, p) in pos.iter().enumerate() {
            assert!((g[a][0] - (4.0 * p[0].powi(3) + p[1])).abs() < 1e-5);
            assert!((g[a][1] - (4.0 * p[1].powi(3) + p[0])).abs() < 1e-5);
            assert!((g[a][2] - 4.0 * p[2].powi(3)).abs() < 1e-5);
        }
    }

    #[test]
    fn propagates_the_energy_error() {
        let r = central_difference_gradient(&[[0.0; 3]], 1e-3, |_| Err("no SCF".to_string()));
        assert_eq!(r.unwrap_err(), "no SCF");
    }
}
