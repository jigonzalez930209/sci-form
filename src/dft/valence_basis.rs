//! All-electron minimal basis for elements outside the Hehre STO-3G table.
//!
//! Exponents follow Slater's rules, ζ = (Z − σ) / n*. Each Slater orbital is
//! the same three-Gaussian fit used by STO-3G, with α_i = α_i(ζ=1) · ζ².
//! An empty valence p shell is kept so the Kohn–Sham problem has virtuals.
//! This is the all-electron path through Zn. Heavier metals need ECPs.

use crate::scf::basis::{get_sto3g_shells, sto3g_tabulated, BasisSet, ContractedShell, GaussianPrimitive};

const S_CORE: [(f64, f64); 3] = [
    (2.227660, 0.15432897),
    (0.405771, 0.53532814),
    (0.109818, 0.44463454),
];
const S_VAL: [(f64, f64); 3] = [
    (0.994203, -0.09996723),
    (0.231031, 0.39951283),
    (0.075139, 0.70011547),
];
const P_VAL: [(f64, f64); 3] = [
    (0.994203, 0.15591627),
    (0.231031, 0.60768372),
    (0.075139, 0.39195739),
];

pub fn ks_basis(elements: &[u8], positions_bohr: &[[f64; 3]]) -> Result<BasisSet, String> {
    if let Some(z) = elements.iter().copied().find(|z| *z > 30) {
        return Err(format!(
            "KS-DFT all-electron minimal basis covers Z <= 30 (through Zn). Z={z} needs an ECP, which is the next step of this cycle."
        ));
    }
    let mut shells = Vec::new();
    for (atom_idx, (&z, &center)) in elements.iter().zip(positions_bohr.iter()).enumerate() {
        if sto3g_tabulated(z) {
            shells.extend(get_sto3g_shells(z, atom_idx, center));
        } else {
            shells.extend(slater_shells(z, atom_idx, center));
        }
    }
    let basis = BasisSet::from_shells(shells);
    let n_electrons: usize = elements.iter().map(|z| *z as usize).sum();
    if basis.n_basis * 2 < n_electrons {
        return Err(format!(
            "minimal basis has {} functions for {} electrons",
            basis.n_basis, n_electrons
        ));
    }
    Ok(basis)
}

fn n_star(n: u8) -> f64 {
    match n {
        1 => 1.0,
        2 => 2.0,
        3 => 3.0,
        _ => 3.7,
    }
}

/// Slater groups: 1s, 2s+2p, 3s+3p, 3d, 4s.
fn groups(z: u8) -> [f64; 5] {
    let mut left = z as i32;
    let mut take = |n: i32| {
        let k = left.min(n).max(0);
        left -= k;
        k as f64
    };
    let g1 = take(2);
    let g2 = take(8);
    let g3 = take(8);
    if z >= 19 {
        let g4s = take(2);
        let gd = take(10);
        [g1, g2, g3, gd, g4s]
    } else {
        [g1, g2, g3, 0.0, 0.0]
    }
}

fn zeta(z: u8, shell: u8) -> f64 {
    let g = groups(z);
    let zf = z as f64;
    let sigma = match shell {
        0 => 0.30 * (g[0] - 1.0).max(0.0),
        1 => 0.85 * g[0] + 0.35 * (g[1] - 1.0).max(0.0),
        2 => 1.0 * g[0] + 0.85 * g[1] + 0.35 * (g[2] - 1.0).max(0.0),
        3 => 1.0 * (g[0] + g[1] + g[2]) + 0.35 * (g[3] - 1.0).max(0.0),
        _ => 1.0 * (g[0] + g[1]) + 0.85 * (g[2] + g[3]) + 0.35 * (g[4] - 1.0).max(0.0),
    };
    let n = match shell {
        0 => 1,
        1 => 2,
        2 | 3 => 3,
        _ => 4,
    };
    ((zf - sigma) / n_star(n)).clamp(0.4, 30.0)
}

fn shell(atom_index: usize, center: [f64; 3], l: u32, zeta: f64, fit: &[(f64, f64)]) -> ContractedShell {
    let z2 = zeta * zeta;
    ContractedShell {
        atom_index,
        center,
        l,
        primitives: fit
            .iter()
            .map(|&(alpha, coefficient)| GaussianPrimitive {
                alpha: alpha * z2,
                coefficient,
            })
            .collect(),
    }
}

fn slater_shells(z: u8, atom_index: usize, center: [f64; 3]) -> Vec<ContractedShell> {
    let mut shells = vec![shell(atom_index, center, 0, zeta(z, 0), &S_CORE)];
    if z >= 3 {
        shells.push(shell(atom_index, center, 0, zeta(z, 1), &S_VAL));
    }
    if z >= 4 {
        shells.push(shell(atom_index, center, 1, zeta(z, 1), &P_VAL));
    }
    if z >= 11 {
        shells.push(shell(atom_index, center, 0, zeta(z, 2), &S_VAL));
    }
    if z >= 12 {
        shells.push(shell(atom_index, center, 1, zeta(z, 2), &P_VAL));
    }
    if z >= 19 {
        shells.push(shell(atom_index, center, 0, zeta(z, 4), &S_VAL));
    }
    if z >= 21 {
        shells.push(shell(atom_index, center, 2, zeta(z, 3), &P_VAL));
    }
    if z >= 20 {
        shells.push(shell(atom_index, center, 1, zeta(z, 4), &P_VAL));
    }
    shells
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iron_basis_has_a_d_shell_and_virtuals() {
        let basis = ks_basis(&[26], &[[0.0, 0.0, 0.0]]).unwrap();
        assert!(basis.functions.iter().any(|f| f.l_total == 2));
        assert!(basis.n_basis * 2 > 26);
    }

    #[test]
    fn hydrogen_stays_on_the_hehre_table() {
        let basis = ks_basis(&[1, 1], &[[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]).unwrap();
        assert_eq!(basis.n_basis, 2);
    }

    #[test]
    fn past_zinc_asks_for_an_ecp() {
        let err = ks_basis(&[31], &[[0.0, 0.0, 0.0]]).unwrap_err();
        assert!(err.contains("ECP"));
    }
}
