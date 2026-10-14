//! G(4,1) Multivector — the core algebraic type for Conformal Geometric Algebra.
//!
//! The conformal model embeds 3D Euclidean space into G(4,1) with basis vectors
//! e1, e2, e3 (Euclidean, square to +1) plus e+ (squares to +1) and e- (squares
//! to -1).  The null basis is e_o = (e- - e+)/2 (origin) and e_inf = e- + e+
//! (point at infinity).
//!
//! A general multivector has 2^5 = 32 components.

use std::ops::{Add, BitAnd, BitOr, Mul, Neg, Sub};

/// Number of basis blades in G(4,1).
const N_BLADES: usize = 32;

/// Basis blade labels for reference (bit-encoded: bit0=e1, bit1=e2, bit2=e3,
/// bit3=e+, bit4=e-).
const BASIS_LABELS: [&str; 32] = [
    "1", "e1", "e2", "e12", "e3", "e13", "e23", "e123", "e+", "e1+", "e2+", "e12+", "e3+", "e13+",
    "e23+", "e123+", "e-", "e1-", "e2-", "e12-", "e3-", "e13-", "e23-", "e123-", "e+-", "e1+-",
    "e2+-", "e12+-", "e3+-", "e13+-", "e23+-", "e123+-",
];

/// Metric signature: e1²=+1, e2²=+1, e3²=+1, e+²=+1, e-²=-1
const METRIC: [f64; 5] = [1.0, 1.0, 1.0, 1.0, -1.0];

/// A multivector in the G(4,1) Conformal Geometric Algebra.
///
/// Components are stored in a 32-element array indexed by blade bitmask.
#[derive(Debug, Clone, Copy)]
pub struct Multivector {
    pub data: [f64; N_BLADES],
}

impl Multivector {
    /// The zero multivector.
    pub const ZERO: Self = Self {
        data: [0.0; N_BLADES],
    };

    /// The scalar 1.
    pub fn one() -> Self {
        let mut m = Self::ZERO;
        m.data[0] = 1.0;
        m
    }

    /// Basis vector by index (0..4 → e1..e-).
    pub fn basis(i: usize) -> Self {
        assert!(i < 5);
        let mut m = Self::ZERO;
        m.data[1 << i] = 1.0;
        m
    }

    /// Construct from a single blade.
    pub fn blade(index: usize, value: f64) -> Self {
        let mut m = Self::ZERO;
        m.data[index] = value;
        m
    }

    /// The null origin point: e_o = (e- - e+) / 2
    pub fn origin() -> Self {
        let mut m = Self::ZERO;
        m.data[1 << 4] = 0.5; // e-
        m.data[1 << 3] = -0.5; // -e+
        m
    }

    /// The point at infinity: e_inf = e- + e+
    pub fn infinity() -> Self {
        let mut m = Self::ZERO;
        m.data[1 << 4] = 1.0; // e-
        m.data[1 << 3] = 1.0; // e+
        m
    }

    /// Scalar part.
    pub fn scalar(&self) -> f64 {
        self.data[0]
    }

    /// Get component of a specific blade.
    pub fn component(&self, blade: usize) -> f64 {
        self.data[blade]
    }

    /// Grade of a blade (popcount of its bitmask).
    fn blade_grade(blade: usize) -> u32 {
        (blade as u32).count_ones()
    }

    /// Extract grade-k part.
    pub fn grade(&self, k: u32) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            if Self::blade_grade(i) == k {
                m.data[i] = self.data[i];
            }
        }
        m
    }

    /// Reverse: reverse the order of factors in each blade.
    /// For a grade-k blade: rev = (-1)^{k(k-1)/2} * blade
    pub fn reverse(&self) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            let k = Self::blade_grade(i);
            // k*(k-1)/2 mod 2: for k=0 → 0, k=1 → 0, k=2 → 1, k=3 → 1, k=4 → 0, k=5 → 0
            let sign = if (k / 2) % 2 == 0 { 1.0 } else { -1.0 };
            m.data[i] = sign * self.data[i];
        }
        m
    }

    /// Grade involution: (-1)^k for each grade-k blade.
    pub fn involute(&self) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            let k = Self::blade_grade(i);
            let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
            m.data[i] = sign * self.data[i];
        }
        m
    }

    /// Conjugate = reverse ∘ grade involution.
    pub fn conjugate(&self) -> Self {
        self.reverse().involute()
    }

    /// Squared norm: M * ~M (scalar part).
    pub fn norm_squared(&self) -> f64 {
        (*self * self.reverse()).scalar()
    }

    /// Norm: sqrt(|M ~M|).
    pub fn norm(&self) -> f64 {
        self.norm_squared().abs().sqrt()
    }

    /// Normalize to unit multivector.
    pub fn normalized(&self) -> Self {
        let n = self.norm();
        if n < 1e-15 {
            return *self;
        }
        *self * (1.0 / n)
    }

    /// Geometric product sign and resulting blade for two basis blades.
    /// Uses the canonical reordering with metric signature.
    fn geo_sign_blade(a: usize, b: usize) -> (f64, usize) {
        let mut sign = 1.0;
        let mut result = a;

        for j in 0..5 {
            if b & (1 << j) == 0 {
                continue;
            }
            // Count how many bits in `result` above position j need to be swapped past
            let mut swaps = 0;
            for k in (j + 1)..5 {
                if result & (1 << k) != 0 {
                    swaps += 1;
                }
            }
            if swaps % 2 != 0 {
                sign = -sign;
            }

            if result & (1 << j) != 0 {
                // e_j * e_j = metric[j]
                sign *= METRIC[j];
                result ^= 1 << j;
            } else {
                result ^= 1 << j;
            }
        }

        (sign, result)
    }

    /// Geometric product.
    pub fn geo(&self, other: &Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..N_BLADES {
            if self.data[i] == 0.0 {
                continue;
            }
            for j in 0..N_BLADES {
                if other.data[j] == 0.0 {
                    continue;
                }
                let (sign, blade) = Self::geo_sign_blade(i, j);
                result.data[blade] += sign * self.data[i] * other.data[j];
            }
        }
        result
    }

    /// Outer (wedge) product: grade(a)+grade(b) part of geo product.
    pub fn outer(&self, other: &Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..N_BLADES {
            if self.data[i] == 0.0 {
                continue;
            }
            let gi = Self::blade_grade(i);
            for j in 0..N_BLADES {
                if other.data[j] == 0.0 {
                    continue;
                }
                let gj = Self::blade_grade(j);
                let (sign, blade) = Self::geo_sign_blade(i, j);
                if Self::blade_grade(blade) == gi + gj {
                    result.data[blade] += sign * self.data[i] * other.data[j];
                }
            }
        }
        result
    }

    /// Inner (left contraction) product.
    pub fn inner(&self, other: &Self) -> Self {
        let mut result = Self::ZERO;
        for i in 0..N_BLADES {
            if self.data[i] == 0.0 {
                continue;
            }
            let gi = Self::blade_grade(i);
            for j in 0..N_BLADES {
                if other.data[j] == 0.0 {
                    continue;
                }
                let gj = Self::blade_grade(j);
                if gj < gi {
                    continue;
                }
                let (sign, blade) = Self::geo_sign_blade(i, j);
                if Self::blade_grade(blade) == gj - gi {
                    result.data[blade] += sign * self.data[i] * other.data[j];
                }
            }
        }
        result
    }

    /// Sandwich product: M * X * ~M
    pub fn sandwich(&self, x: &Self) -> Self {
        self.geo(x).geo(&self.reverse())
    }

    /// Check if effectively zero.
    pub fn is_zero(&self, tol: f64) -> bool {
        self.data.iter().all(|&v| v.abs() < tol)
    }

    /// Label of the i-th basis blade (for debug).
    pub fn blade_label(i: usize) -> &'static str {
        BASIS_LABELS[i]
    }
}

// ── Operator overloads ───────────────────────────────────────────────────────

impl Add for Multivector {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            m.data[i] = self.data[i] + rhs.data[i];
        }
        m
    }
}

impl Sub for Multivector {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            m.data[i] = self.data[i] - rhs.data[i];
        }
        m
    }
}

impl Neg for Multivector {
    type Output = Self;
    fn neg(self) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            m.data[i] = -self.data[i];
        }
        m
    }
}

impl Mul<Multivector> for Multivector {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        self.geo(&rhs)
    }
}

impl Mul<f64> for Multivector {
    type Output = Self;
    fn mul(self, rhs: f64) -> Self {
        let mut m = Self::ZERO;
        for i in 0..N_BLADES {
            m.data[i] = self.data[i] * rhs;
        }
        m
    }
}

impl Mul<Multivector> for f64 {
    type Output = Multivector;
    fn mul(self, rhs: Multivector) -> Multivector {
        rhs * self
    }
}

/// Outer product operator: a ^ b
impl BitAnd for Multivector {
    type Output = Self;
    fn bitand(self, rhs: Self) -> Self {
        self.outer(&rhs)
    }
}

/// Inner product operator: a | b
impl BitOr for Multivector {
    type Output = Self;
    fn bitor(self, rhs: Self) -> Self {
        self.inner(&rhs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basis_anticommute() {
        // e_i e_j = -e_j e_i for i ≠ j
        for i in 0..5 {
            for j in (i + 1)..5 {
                let ei = Multivector::basis(i);
                let ej = Multivector::basis(j);
                let eij = ei * ej;
                let eji = ej * ei;
                let sum = eij + eji;
                assert!(
                    sum.is_zero(1e-14),
                    "e{} e{} + e{} e{} is not zero",
                    i,
                    j,
                    j,
                    i
                );
            }
        }
    }

    #[test]
    fn test_basis_metric() {
        // e_i^2 = metric[i]
        for i in 0..5 {
            let ei = Multivector::basis(i);
            let sq = ei * ei;
            assert!(
                (sq.scalar() - METRIC[i]).abs() < 1e-14,
                "e{}^2 = {}, expected {}",
                i,
                sq.scalar(),
                METRIC[i]
            );
        }
    }

    #[test]
    fn test_null_basis_properties() {
        // e+ e- + e- e+ = 0 (anticommute)
        let ep = Multivector::basis(3); // e+
        let em = Multivector::basis(4); // e-
        let sum = ep * em + em * ep;
        assert!(sum.is_zero(1e-14), "e+ e- + e- e+ should be 0");

        // e_o · e_inf = -1
        let eo = Multivector::origin();
        let einf = Multivector::infinity();
        let dot = eo | einf;
        assert!(
            (dot.scalar() - (-1.0)).abs() < 1e-14,
            "e_o · e_inf = {}, expected -1",
            dot.scalar()
        );
    }

    #[test]
    fn test_reverse_identity_for_scalars() {
        let s = Multivector::blade(0, std::f64::consts::PI);
        let rev = s.reverse();
        assert!((rev.scalar() - std::f64::consts::PI).abs() < 1e-14);
    }

    #[test]
    fn test_unit_motor_reverse_is_identity() {
        // For a unit rotor R, R ~R = 1
        let angle = std::f64::consts::FRAC_PI_4; // 45°
        let half = angle / 2.0;
        // Rotor in e12 plane: R = cos(θ/2) - sin(θ/2) e12
        let mut r = Multivector::ZERO;
        r.data[0] = half.cos(); // scalar
        r.data[3] = -half.sin(); // e12 blade (index = 0b0011 = 3)
        let product = r * r.reverse();
        assert!(
            (product.scalar() - 1.0).abs() < 1e-12,
            "R ~R scalar = {}, expected 1",
            product.scalar()
        );
        // All non-scalar components should be ~0
        for i in 1..32 {
            assert!(
                product.data[i].abs() < 1e-12,
                "R ~R blade {} = {}, expected 0",
                i,
                product.data[i]
            );
        }
    }
}
