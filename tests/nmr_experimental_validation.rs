//! Experimental ³J(H,H) anchors for the NMR coupling predictor.
//!
//! The predictor uses a Karplus family; F54 validates the *equation* against
//! the Stenutz HLA reference calculator. This file adds the missing check
//! against measured couplings for the cases where the molecular topology and
//! the observable are unambiguous.
//!
//! Benzene: ortho ³J(H,H) = 7.54 Hz (liquid/gas consensus; the sign is
//! positive and the magnitude is the conventional reference for aromatic
//! ortho coupling). The embedded D6h geometry gives one torsion per adjacent
//! H pair, so all six ³J values should be equal.

use sci_form::{embed, predict_nmr_couplings};

#[test]
fn benzene_ortho_3j_matches_experimental() {
    let conformer = embed("c1ccccc1", 42);
    assert!(
        conformer.error.is_none(),
        "embed failed: {:?}",
        conformer.error
    );
    let positions: Vec<[f64; 3]> = conformer
        .coords
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| [p[0], p[1], p[2]])
        .collect();

    let couplings = predict_nmr_couplings("c1ccccc1", &positions).expect("J-couplings");
    let ortho: Vec<f64> = couplings
        .iter()
        .filter(|c| c.n_bonds == 3 && c.coupling_type.contains("vicinal_3J_H-C-C-H"))
        .map(|c| c.j_hz)
        .collect();
    assert_eq!(ortho.len(), 6, "benzene has six adjacent H pairs");
    let mean = ortho.iter().sum::<f64>() / ortho.len() as f64;
    let spread = ortho
        .iter()
        .map(|j| (j - mean).abs())
        .fold(0.0f64, f64::max);
    eprintln!("Benzene ortho ³J: model {mean:.3} Hz (spread {spread:.2e}), experimental 7.54 Hz");
    assert!(spread < 1e-6, "D6h symmetry broken: spread {spread:.3e} Hz");
    assert!(
        (mean - 7.54).abs() <= 1.0,
        "benzene ³J {mean:.3} Hz vs experimental 7.54 Hz"
    );
}

/// Ethylene is the sp2 anchor: the planar H-C=C-H pairs are the cis
/// (φ≈0, 11.7 Hz) and trans (φ≈180, 19.0 Hz) couplings. The sp2-specific
/// Karplus parameters are fitted to those two measured values; the sp3
/// parameters used for single C-C bonds would give ~8-10 Hz for both.
#[test]
fn ethylene_cis_trans_3j_matches_experimental() {
    let conformer = embed("C=C", 42);
    assert!(
        conformer.error.is_none(),
        "embed failed: {:?}",
        conformer.error
    );
    let positions: Vec<[f64; 3]> = conformer
        .coords
        .as_chunks::<3>()
        .0
        .iter()
        .map(|p| [p[0], p[1], p[2]])
        .collect();
    let couplings = predict_nmr_couplings("C=C", &positions).expect("J-couplings");
    let mut vicinal: Vec<f64> = couplings
        .iter()
        .filter(|c| c.n_bonds == 3)
        .map(|c| c.j_hz)
        .collect();
    assert_eq!(vicinal.len(), 4, "ethylene has four vicinal H pairs");
    vicinal.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let cis = 0.5 * (vicinal[0] + vicinal[1]);
    let trans = 0.5 * (vicinal[2] + vicinal[3]);
    eprintln!("Ethylene ³J: cis {cis:.3} Hz (exp 11.7), trans {trans:.3} Hz (exp 19.0)");
    assert!((cis - 11.7).abs() <= 1.0, "cis ³J {cis:.3} Hz");
    assert!((trans - 19.0).abs() <= 1.0, "trans ³J {trans:.3} Hz");
}
