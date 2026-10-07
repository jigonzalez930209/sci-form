//! Hydrogenic EHT parameters for metals absent from the literature table.
//!
//! The diagonal is the experimental first ionization potential. The Slater
//! exponent follows from that binding energy, ζ = sqrt(IP / 13.6057), which is
//! the hydrogenic single-zeta that reproduces the same ionization energy.
//! Lanthanides and actinides are treated with the f shell in the core: the
//! valence basis is ns, np and (n−1)d, because the overlap engine is s/p/d.

use std::collections::HashMap;
use std::sync::OnceLock;

use super::params::{EhtParams, OrbitalDef};

const HARTREE_EV: f64 = 13.605693;

pub fn is_derived(z: u8) -> bool {
    derived_params(z).is_some()
}

pub fn valence_electrons(z: u8) -> usize {
    match z {
        3 | 11 | 19 | 37 | 55 | 87 => 1,
        4 | 12 | 20 | 38 | 56 | 88 => 2,
        13 | 31 | 49 | 81 => 3,
        50 | 82 => 4,
        83 => 5,
        63 | 70 => 2, // Eu, Yb: 6s², empty d in the neutral atom
        57..=71 => 3,
        89 => 3,
        90 => 4,
        91 => 5,
        92 => 6,
        93 => 7,
        94 => 8,
        95 => 9,
        96 => 10,
        97 => 11,
        98 => 12,
        99 => 13,
        100 => 14,
        101 => 15,
        102 => 16,
        103 => 3,
        _ => 0,
    }
}

fn first_ip_ev(z: u8) -> f64 {
    match z {
        3 => 5.392,
        4 => 9.323,
        11 => 5.139,
        12 => 7.646,
        13 => 5.986,
        19 => 4.341,
        20 => 6.113,
        31 => 5.999,
        37 => 4.177,
        38 => 5.695,
        49 => 5.786,
        50 => 7.344,
        55 => 3.894,
        56 => 5.212,
        57 => 5.577,
        58 => 5.539,
        59 => 5.473,
        60 => 5.525,
        61 => 5.582,
        62 => 5.644,
        63 => 5.670,
        64 => 6.150,
        65 => 5.864,
        66 => 5.939,
        67 => 6.021,
        68 => 6.108,
        69 => 6.184,
        70 => 6.254,
        71 => 5.426,
        81 => 6.108,
        82 => 7.417,
        83 => 7.286,
        87 => 4.073,
        88 => 5.279,
        89 => 5.170,
        90 => 6.307,
        91 => 5.890,
        92 => 6.194,
        93 => 6.266,
        94 => 6.026,
        95 => 5.974,
        _ => 6.0,
    }
}

fn valence_n(z: u8) -> u8 {
    match z {
        3 | 4 => 2,
        11..=13 => 3,
        19 | 20 | 31 => 4,
        37 | 38 | 49 | 50 => 5,
        55..=83 => 6,
        _ => 7,
    }
}

fn shell_label(n: u8, l: u8) -> &'static str {
    match (n, l) {
        (2, 0) => "2s",
        (3, 0) => "3s",
        (3, 1) => "3p",
        (4, 0) => "4s",
        (4, 1) => "4p",
        (4, 2) => "3d",
        (5, 0) => "5s",
        (5, 1) => "5p",
        (5, 2) => "4d",
        (6, 0) => "6s",
        (6, 1) => "6p",
        (6, 2) => "5d",
        (7, 0) => "7s",
        (7, 1) => "7p",
        (7, 2) => "6d",
        _ => "val",
    }
}

fn symbol(z: u8) -> &'static str {
    match z {
        3 => "Li",
        4 => "Be",
        11 => "Na",
        12 => "Mg",
        13 => "Al",
        19 => "K",
        20 => "Ca",
        31 => "Ga",
        37 => "Rb",
        38 => "Sr",
        49 => "In",
        50 => "Sn",
        55 => "Cs",
        56 => "Ba",
        57 => "La",
        58 => "Ce",
        59 => "Pr",
        60 => "Nd",
        61 => "Pm",
        62 => "Sm",
        63 => "Eu",
        64 => "Gd",
        65 => "Tb",
        66 => "Dy",
        67 => "Ho",
        68 => "Er",
        69 => "Tm",
        70 => "Yb",
        71 => "Lu",
        81 => "Tl",
        82 => "Pb",
        83 => "Bi",
        87 => "Fr",
        88 => "Ra",
        89 => "Ac",
        90 => "Th",
        91 => "Pa",
        92 => "U",
        93 => "Np",
        94 => "Pu",
        95 => "Am",
        96 => "Cm",
        97 => "Bk",
        98 => "Cf",
        99 => "Es",
        100 => "Fm",
        101 => "Md",
        102 => "No",
        103 => "Lr",
        _ => "M",
    }
}

/// ζ in bohr⁻¹ from a hydrogenic ionization energy.
fn zeta_from_ip(ip_ev: f64) -> f64 {
    (ip_ev.max(1.0) / HARTREE_EV).sqrt().clamp(0.35, 6.0)
}

fn uses_d_valence(z: u8) -> bool {
    matches!(z, 57..=71 | 89..=103)
}

fn build_orbitals(z: u8) -> Vec<OrbitalDef> {
    let n = valence_n(z);
    let ip = first_ip_ev(z);
    let mut orbitals = vec![OrbitalDef {
        n,
        l: 0,
        label: shell_label(n, 0),
        vsip: -ip,
        zeta: zeta_from_ip(ip),
    }];
    let with_p = matches!(z, 13 | 31 | 49 | 50 | 81 | 82 | 83) || uses_d_valence(z);
    if with_p {
        let ip_p = (ip - 3.5).max(2.0);
        orbitals.push(OrbitalDef {
            n,
            l: 1,
            label: shell_label(n, 1),
            vsip: -ip_p,
            zeta: zeta_from_ip(ip_p),
        });
    }
    if uses_d_valence(z) {
        let ip_d = ip + 2.0;
        orbitals.push(OrbitalDef {
            n: n - 1,
            l: 2,
            label: shell_label(n, 2),
            vsip: -ip_d,
            zeta: zeta_from_ip(ip_d),
        });
    }
    orbitals
}

pub fn derived_params(z: u8) -> Option<&'static EhtParams> {
    if valence_electrons(z) == 0 {
        return None;
    }
    static CACHE: OnceLock<HashMap<u8, &'static EhtParams>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| {
        let mut map = HashMap::new();
        for z in 1u8..=103 {
            if valence_electrons(z) == 0 {
                continue;
            }
            let orbitals = Box::leak(build_orbitals(z).into_boxed_slice());
            let params = Box::leak(Box::new(EhtParams {
                z,
                symbol: symbol(z),
                orbitals,
            }));
            map.insert(z, &*params);
        }
        map
    });
    cache.get(&z).copied()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::eht::basis::sto3g_expansion;
    use crate::eht::params::{get_params, num_basis_functions};
    use crate::eht::solver::solve_eht;
    use crate::reactivity::collect_stda_excitations;

    #[test]
    fn lithium_zeta_matches_hydrogenic_ip() {
        let li = get_params(3).expect("Li");
        let expected = (5.392_f64 / HARTREE_EV).sqrt();
        assert!((li.orbitals[0].zeta - expected).abs() < 1e-6);
        assert!((li.orbitals[0].vsip + 5.392).abs() < 1e-9);
        assert_eq!(num_basis_functions(3), 1);
    }

    #[test]
    fn literature_iron_is_not_replaced() {
        let fe = get_params(26).unwrap();
        assert!(fe.orbitals.iter().any(|o| o.label == "3d"));
        assert!(fe.orbitals[0].vsip < -5.0);
    }

    #[test]
    fn sixth_shell_primitives_are_not_empty() {
        assert_eq!(sto3g_expansion(6, 0, 1.2).len(), 3);
        assert_eq!(sto3g_expansion(7, 2, 1.5).len(), 3);
        assert!(sto3g_expansion(6, 0, 1.2).iter().all(|g| g.alpha > 0.0));
    }

    #[test]
    fn lanthanum_and_uranium_solve_and_emit() {
        for metal in [57u8, 92] {
            let elements = [metal, 6, 1, 1, 1];
            let positions = [
                [0.0, 0.0, 0.0],
                [2.2, 0.0, 0.0],
                [2.6, 1.0, 0.0],
                [2.6, -0.5, 0.9],
                [2.6, -0.5, -0.9],
            ];
            let eht = solve_eht(&elements, &positions, None).expect("EHT");
            assert!(eht.n_electrons >= 7, "Z={metal}");
            assert!(eht.gap.is_finite());
            let excitations = collect_stda_excitations(&elements, &positions, 12.0).unwrap();
            assert!(
                excitations.iter().any(|e| e.oscillator_strength > 1e-6),
                "Z={metal} should have a bright excitation"
            );
        }
    }

    #[test]
    fn open_shell_counts_the_somo_as_occupied() {
        let elements = [7u8, 8];
        let positions = [[0.0, 0.0, 0.0], [1.15, 0.0, 0.0]];
        let eht = solve_eht(&elements, &positions, None).unwrap();
        assert_eq!(eht.n_electrons, 11);
        let excitations = collect_stda_excitations(&elements, &positions, 20.0).unwrap();
        let n_occ = eht.n_electrons.div_ceil(2);
        assert!(excitations.iter().any(|ex| ex.from_mo + 1 == n_occ));
    }
}
