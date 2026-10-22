//! Small-core effective core potentials for Z > 30.
//!
//! Core electrons are replaced by Z_eff = Z − n_core and a local Gaussian
//! potential on the nucleus. f electrons stay in the core. This is an alpha
//! ECP: it unblocks 4d/5d/p-block metals; it is not a fitted Stuttgart/LANL set.

use crate::scf::basis::{BasisFunction, BasisSet, ContractedShell};
use crate::scf::gaussian_integrals::overlap_primitive;

/// Small-core electrons removed from the explicit Hilbert space.
pub fn n_core(z: u8) -> u8 {
    match z {
        0..=30 => 0,
        31..=36 => 18,
        37..=38 => 28,
        39..=48 => 28,
        49..=54 => 36,
        55..=56 => 46,
        57..=71 => 54,
        72..=80 => 60,
        81..=86 => 68,
        87..=103 => 78,
        _ => 0,
    }
}

/// Effective nuclear charge `Z_eff = Z - n_core(Z)` in units of `e`.
///
/// The core electrons removed by [`n_core`] are replaced by a local
/// Gaussian potential (see [`local_gaussian`]), so the explicit valence
/// space sees a reduced nuclear charge. The result is clamped from below at
/// `1.0`, so a hypothetical element with more core electrons than valence
/// electrons still reports a positive charge.
pub fn z_eff(z: u8) -> f64 {
    (z as f64 - n_core(z) as f64).max(1.0)
}

/// Effective nuclear charges of a whole molecule, in units of `e`.
///
/// The returned vector has one entry per element in `elements` and follows
/// the same atom ordering as the input slice.
pub fn effective_charges(elements: &[u8]) -> Vec<f64> {
    elements.iter().map(|&z| z_eff(z)).collect()
}

/// Local Gaussian ECP: U(r) = A exp(−α r²).
///
/// Returns `(amplitude, exponent)` with the amplitude in Hartree and the
/// exponent in Bohr⁻² (the potential is evaluated with `r` in Bohr).
///
/// `None` is returned for elements that have no core electrons to remove
/// (`Z <= 30`), i.e. when no ECP is needed.
pub fn local_gaussian(z: u8) -> Option<(f64, f64)> {
    let n = n_core(z);
    if n == 0 {
        return None;
    }
    let nf = n as f64;
    Some((0.5 * nf, 0.15 * nf.powf(2.0 / 3.0)))
}

/// Add the local Gaussian ECP to the nuclear attraction matrix.
///
/// `v` is the nuclear attraction matrix in Hartree; the ECP contribution
/// (also in Hartree) is added to it in place. The basis centres are in
/// Bohr.
pub fn add_local_ecp(v: &mut nalgebra::DMatrix<f64>, basis: &BasisSet, elements: &[u8]) {
    for (atom, &z) in elements.iter().enumerate() {
        let Some((amp, alpha)) = local_gaussian(z) else {
            continue;
        };
        let center = basis
            .functions
            .iter()
            .find(|bf| bf.atom_index == atom)
            .map(|bf| bf.center)
            .unwrap_or([0.0; 3]);
        let n = basis.n_basis;
        for i in 0..n {
            for j in i..n {
                let vij = contracted_local(
                    &basis.functions[i],
                    &basis.functions[j],
                    alpha,
                    &center,
                    amp,
                );
                v[(i, j)] += vij;
                if i != j {
                    v[(j, i)] += vij;
                }
            }
        }
    }
}

fn contracted_local(
    a: &BasisFunction,
    b: &BasisFunction,
    gamma: f64,
    c: &[f64; 3],
    amp: f64,
) -> f64 {
    let mut s = 0.0;
    for pa in &a.primitives {
        let na = BasisFunction::normalization(pa.alpha, a.angular[0], a.angular[1], a.angular[2]);
        for pb in &b.primitives {
            let nb =
                BasisFunction::normalization(pb.alpha, b.angular[0], b.angular[1], b.angular[2]);
            let p = pa.alpha + pb.alpha;
            let px = (pa.alpha * a.center[0] + pb.alpha * b.center[0]) / p;
            let py = (pa.alpha * a.center[1] + pb.alpha * b.center[1]) / p;
            let pz = (pa.alpha * a.center[2] + pb.alpha * b.center[2]) / p;
            let cart = overlap_primitive(
                pa.alpha, &a.center, a.angular, pb.alpha, &b.center, b.angular,
            );
            let s_ss = overlap_primitive(p, &[px, py, pz], [0, 0, 0], gamma, c, [0, 0, 0]);
            let s_ab_ss = overlap_primitive(
                pa.alpha,
                &a.center,
                [0, 0, 0],
                pb.alpha,
                &b.center,
                [0, 0, 0],
            );
            let scale = if s_ab_ss.abs() > 1e-30 {
                s_ss / s_ab_ss
            } else {
                0.0
            };
            s += na * pa.coefficient * nb * pb.coefficient * amp * cart * scale;
        }
    }
    s
}

/// Principal valence quantum number `n` used to size the valence basis.
///
/// * `4` for Ga(31) – Xe(36)
/// * `5` for Rb(37) – Xe(54)
/// * `6` for Cs(55) – Rn(86)
/// * `7` for everything else (including the light elements and the
///   superheavy actinides, which fall back to the largest value).
pub fn valence_n(z: u8) -> u8 {
    match z {
        31..=36 => 4,
        37..=54 => 5,
        55..=86 => 6,
        _ => 7,
    }
}

/// Whether the element is given an explicit d-type valence shell.
///
/// True for the transition metals Ga(31) – Kr(36), Y(39) – Cd(48) and
/// Hf(72) – Hg(80) plus the lanthanides La(57) – Lu(71): the metals whose
/// valence space contains a filled or partially filled `(n-1)d` sub-shell
/// that must be represented explicitly. False for the main-group p-block
/// metals, the alkali and alkaline-earth metals, and the actinides.
pub fn needs_d_valence(z: u8) -> bool {
    matches!(z, 31..=36 | 39..=48 | 57..=80)
}

/// Re-export shell builder pieces used by the valence ECP basis.
pub fn ecp_shells(
    z: u8,
    atom_index: usize,
    center: [f64; 3],
    make_shell: impl Fn(usize, [f64; 3], u32, f64, &[(f64, f64)]) -> ContractedShell,
    s_val: &[(f64, f64)],
    p_val: &[(f64, f64)],
) -> Vec<ContractedShell> {
    let zeff = z_eff(z);
    let n = valence_n(z);
    let nstar = match n {
        4 => 3.7,
        5 => 4.0,
        6 => 4.2,
        _ => 4.4,
    };
    let zeta_s = (zeff / nstar * 0.35).clamp(0.25, 20.0);
    let zeta_p = (zeff / nstar * 0.30).clamp(0.20, 18.0);
    let zeta_d = (zeff / (nstar - 0.7) * 0.45).clamp(0.30, 22.0);
    let mut shells = vec![make_shell(atom_index, center, 0, zeta_s, s_val)];
    shells.push(make_shell(atom_index, center, 1, zeta_p, p_val));
    if needs_d_valence(z) {
        shells.push(make_shell(atom_index, center, 2, zeta_d, p_val));
    }
    shells
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zinc_has_no_core() {
        assert_eq!(n_core(30), 0);
        assert!((z_eff(30) - 30.0).abs() < 1e-12);
    }

    #[test]
    fn bromine_uses_ar_core() {
        assert_eq!(n_core(35), 18);
        assert!((z_eff(35) - 17.0).abs() < 1e-12);
        assert!(local_gaussian(35).is_some());
    }

    #[test]
    fn gold_uses_small_core() {
        assert_eq!(n_core(79), 60);
        assert!((z_eff(79) - 19.0).abs() < 1e-12);
    }
}
