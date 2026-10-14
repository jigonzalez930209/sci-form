//! CGA Motor — unified rotation + translation operator.
//!
//! A Motor M = T R combines a rotor R (rotation) and a translator T (translation)
//! into a single element of the even sub-algebra.  Applying M to a conformal point
//! P via the sandwich M P ~M performs both rotation and translation simultaneously.

use super::multivector::Multivector;

/// A Motor in CGA: M = T R (translator composed with rotor).
///
/// Internally stored as a general multivector restricted to even-grade blades.
#[derive(Debug, Clone, Copy)]
pub struct Motor {
    pub mv: Multivector,
}

impl Motor {
    /// Identity motor (no rotation, no translation).
    pub fn identity() -> Self {
        Self {
            mv: Multivector::one(),
        }
    }

    /// Build a rotor from an axis and angle (radians).
    ///
    /// Axis must be a unit 3D vector `[ax, ay, az]`.
    /// R = cos(θ/2) - sin(θ/2) (ax e23 + ay e31 + az e12)
    ///
    /// In the bivector encoding:
    /// - e23 = e2∧e3 → bitmask 0b00110 = 6
    /// - e31 = e3∧e1 → bitmask 0b00101 = 5  (with sign: e31 = -e13)
    /// - e12 = e1∧e2 → bitmask 0b00011 = 3
    pub fn rotor(axis: [f64; 3], angle: f64) -> Self {
        let half = angle / 2.0;
        let c = half.cos();
        let s = half.sin();

        let mut mv = Multivector::ZERO;
        mv.data[0] = c; // scalar

        // The bivector basis for 3D rotations:
        // e23 blade index = 0b00110 = 6
        // e13 blade index = 0b00101 = 5, but e31 = -e13, so ax component → +s*ax on e23,
        //   ay on e31 = -e13, az on e12
        mv.data[6] = -s * axis[0]; // e23
        mv.data[5] = s * axis[1]; // e13 (note: e31 = -e13, so negated)
        mv.data[3] = -s * axis[2]; // e12

        Self { mv }
    }

    /// Build a translator from a 3D translation vector.
    ///
    /// T = 1 - (1/2) t · e_inf
    /// where t · e_inf = tx e1∧e_inf + ty e2∧e_inf + tz e3∧e_inf
    ///
    /// e_inf = e+ + e-  (basis indices 3 and 4 → bitmask 0b01000 and 0b10000)
    /// e1∧e_inf has components on e1∧e+ and e1∧e- blades.
    pub fn translator(t: [f64; 3]) -> Self {
        let mut mv = Multivector::one();

        // e_inf = e+ + e-
        // e1∧e_inf = e1∧e+ + e1∧e-
        //   e1∧e+ → bitmask 0b01001 = 9
        //   e1∧e- → bitmask 0b10001 = 17
        // e2∧e_inf = e2∧e+ + e2∧e-
        //   e2∧e+ → bitmask 0b01010 = 10
        //   e2∧e- → bitmask 0b10010 = 18
        // e3∧e_inf = e3∧e+ + e3∧e-
        //   e3∧e+ → bitmask 0b01100 = 12
        //   e3∧e- → bitmask 0b10100 = 20

        let half = -0.5;
        mv.data[9] += half * t[0]; // e1 e+
        mv.data[17] += half * t[0]; // e1 e-
        mv.data[10] += half * t[1]; // e2 e+
        mv.data[18] += half * t[1]; // e2 e-
        mv.data[12] += half * t[2]; // e3 e+
        mv.data[20] += half * t[2]; // e3 e-

        Self { mv }
    }

    /// Compose a translator and rotor: M = T R (translate after rotate).
    pub fn from_rotation_translation(axis: [f64; 3], angle: f64, translation: [f64; 3]) -> Self {
        let r = Self::rotor(axis, angle);
        let t = Self::translator(translation);
        Self {
            mv: t.mv.geo(&r.mv),
        }
    }

    /// Compose two motors: self ∘ other = self * other.
    pub fn compose(&self, other: &Motor) -> Motor {
        Motor {
            mv: self.mv.geo(&other.mv),
        }
    }

    /// Apply this motor to a conformal point via sandwich: M P ~M.
    pub fn apply(&self, point: &Multivector) -> Multivector {
        self.mv.sandwich(point)
    }

    /// Apply to a 3D point given as `[x, y, z]`.
    /// Lifts to CGA, applies motor, extracts back.
    pub fn transform_point(&self, p: [f64; 3]) -> [f64; 3] {
        let cp = conformal_point(p);
        let result = self.apply(&cp);
        extract_euclidean(&result)
    }

    /// Reverse of the motor.
    pub fn reverse(&self) -> Motor {
        Motor {
            mv: self.mv.reverse(),
        }
    }

    /// Normalize the motor to unit magnitude.
    pub fn normalized(&self) -> Motor {
        Motor {
            mv: self.mv.normalized(),
        }
    }
}

/// Lift a 3D Euclidean point to a CGA null vector.
///
/// P = p + (1/2)|p|² e_inf + e_o
///   = x e1 + y e2 + z e3 + (1/2)(x²+y²+z²) e_inf + e_o
pub fn conformal_point(p: [f64; 3]) -> Multivector {
    let [x, y, z] = p;
    let sq = x * x + y * y + z * z;

    let mut mv = Multivector::ZERO;
    // Euclidean components
    mv.data[1] = x; // e1
    mv.data[2] = y; // e2
    mv.data[4] = z; // e3

    // e_inf = e+ + e- components (scaled by |p|²/2)
    mv.data[1 << 3] += 0.5 * sq; // e+
    mv.data[1 << 4] += 0.5 * sq; // e-

    // e_o = (e- - e+) / 2
    mv.data[1 << 3] -= 0.5; // -0.5 on e+
    mv.data[1 << 4] += 0.5; // +0.5 on e-

    mv
}

/// Extract 3D Euclidean coordinates from a CGA point.
///
/// A CGA point P has the form  α (x e1 + y e2 + z e3 + ... e_inf + e_o)
/// We read x, y, z from e1, e2, e3 components, normalized by the e_o inner product.
pub fn extract_euclidean(p: &Multivector) -> [f64; 3] {
    // The normalization factor: P · e_inf gives the weight.
    // For a properly normalized null vector:  -P · e_inf = 1
    // e_inf = e+ + e-
    // Inner product with grade-1 components
    let einf = Multivector::infinity();
    let dot = p.inner(&einf);
    let weight = -dot.scalar();

    if weight.abs() < 1e-15 {
        return [p.data[1], p.data[2], p.data[4]];
    }

    [
        p.data[1] / weight, // e1 component
        p.data[2] / weight, // e2 component
        p.data[4] / weight, // e3 component
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::FRAC_PI_2;

    #[test]
    fn test_identity_motor() {
        let m = Motor::identity();
        let p = [1.0, 2.0, 3.0];
        let result = m.transform_point(p);
        for i in 0..3 {
            assert!(
                (result[i] - p[i]).abs() < 1e-10,
                "Identity motor changed coordinate {}: {} → {}",
                i,
                p[i],
                result[i]
            );
        }
    }

    #[test]
    fn test_90_deg_rotation_z() {
        // Rotate (1, 0, 0) by 90° about Z → (0, 1, 0)
        let m = Motor::rotor([0.0, 0.0, 1.0], FRAC_PI_2);
        let result = m.transform_point([1.0, 0.0, 0.0]);
        assert!((result[0] - 0.0).abs() < 1e-10, "x = {}", result[0]);
        assert!((result[1] - 1.0).abs() < 1e-10, "y = {}", result[1]);
        assert!((result[2] - 0.0).abs() < 1e-10, "z = {}", result[2]);
    }

    #[test]
    fn test_translation() {
        let m = Motor::translator([1.0, 0.0, 0.0]);
        let result = m.transform_point([0.0, 0.0, 0.0]);
        assert!((result[0] - 1.0).abs() < 1e-10, "x = {}", result[0]);
        assert!((result[1] - 0.0).abs() < 1e-10, "y = {}", result[1]);
        assert!((result[2] - 0.0).abs() < 1e-10, "z = {}", result[2]);
    }

    #[test]
    fn test_rotation_then_translation() {
        // Rotate (1,0,0) 90° about Z → (0,1,0), then translate by (3,0,0) → (3,1,0)
        let m = Motor::from_rotation_translation([0.0, 0.0, 1.0], FRAC_PI_2, [3.0, 0.0, 0.0]);
        let result = m.transform_point([1.0, 0.0, 0.0]);
        assert!((result[0] - 3.0).abs() < 1e-10, "x = {}", result[0]);
        assert!((result[1] - 1.0).abs() < 1e-10, "y = {}", result[1]);
        assert!((result[2] - 0.0).abs() < 1e-10, "z = {}", result[2]);
    }

    #[test]
    fn test_conformal_round_trip() {
        let p = [1.234, -5.678, 9.012];
        let cp = conformal_point(p);
        let ep = extract_euclidean(&cp);
        for i in 0..3 {
            assert!(
                (ep[i] - p[i]).abs() < 1e-10,
                "Round-trip error at {}: {} vs {}",
                i,
                ep[i],
                p[i]
            );
        }
    }

    #[test]
    fn test_motor_composition() {
        // Two 90° rotations about Z should give 180°
        let r90 = Motor::rotor([0.0, 0.0, 1.0], FRAC_PI_2);
        let r180 = r90.compose(&r90);
        let result = r180.transform_point([1.0, 0.0, 0.0]);
        assert!((result[0] - (-1.0)).abs() < 1e-10, "x = {}", result[0]);
        assert!((result[1] - 0.0).abs() < 1e-10, "y = {}", result[1]);
    }
}
