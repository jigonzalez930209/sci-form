//! Frontier-orbital reactivity descriptors derived from EHT molecular orbitals.

use serde::{Deserialize, Serialize};

/// One empirical pKa estimate for a specific atom-centered acid/base site.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmpiricalPkaSite {
    /// Atom index in the input graph.
    pub atom_index: usize,
    /// Site classification: `acidic` or `basic`.
    pub site_type: String,
    /// Textual environment label used by the heuristic rule.
    pub environment: String,
    /// Estimated pKa value from simple charge/environment heuristics.
    pub estimated_pka: f64,
    /// Rule confidence in [0,1].
    pub confidence: f64,
}

/// Empirical pKa summary for acidic and basic candidate sites.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EmpiricalPkaResult {
    /// Ranked acidic sites (lowest pKa first).
    pub acidic_sites: Vec<EmpiricalPkaSite>,
    /// Ranked basic sites (highest pKa first).
    pub basic_sites: Vec<EmpiricalPkaSite>,
    /// Human-readable caveats and guidance.
    pub notes: Vec<String>,
}

/// UFF energy enriched with aromaticity-aware correction metadata.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UffHeuristicEnergy {
    /// Raw UFF energy in kcal/mol.
    pub raw_energy_kcal_mol: f64,
    /// Aromatic stabilization correction in kcal/mol (negative lowers energy).
    pub aromatic_stabilization_kcal_mol: f64,
    /// Corrected heuristic energy in kcal/mol.
    pub corrected_energy_kcal_mol: f64,
    /// Number of aromatic bonds found in the parsed molecular graph.
    pub aromatic_bond_count: usize,
    /// Interpretation notes for downstream consumers.
    pub notes: Vec<String>,
}

/// Atom-resolved frontier-orbital descriptor summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FrontierDescriptors {
    /// Number of atoms in the system.
    pub num_atoms: usize,
    /// HOMO atom contributions, normalized to sum to ~1.
    pub homo_atom_contributions: Vec<f64>,
    /// LUMO atom contributions, normalized to sum to ~1.
    pub lumo_atom_contributions: Vec<f64>,
    /// Simple dual-descriptor proxy: LUMO contribution minus HOMO contribution.
    pub dual_descriptor: Vec<f64>,
    /// HOMO energy in eV.
    pub homo_energy: f64,
    /// LUMO energy in eV.
    pub lumo_energy: f64,
    /// HOMO-LUMO gap in eV.
    pub gap: f64,
}

/// One condensed Fukui descriptor row for an atom.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CondensedFukuiAtom {
    /// Atom index in the input order.
    pub atom_index: usize,
    /// Nucleophilic-attack susceptibility proxy (f+).
    pub f_plus: f64,
    /// Electrophilic-attack susceptibility proxy (f-).
    pub f_minus: f64,
    /// Radical-attack susceptibility proxy (f0).
    pub f_radical: f64,
    /// Dual descriptor (f+ - f-).
    pub dual_descriptor: f64,
}

/// Fukui workflow output, including atom-wise condensed descriptors.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FukuiDescriptors {
    /// Number of atoms in the system.
    pub num_atoms: usize,
    /// Nucleophilic-attack susceptibility proxy per atom.
    pub f_plus: Vec<f64>,
    /// Electrophilic-attack susceptibility proxy per atom.
    pub f_minus: Vec<f64>,
    /// Radical-attack susceptibility proxy per atom.
    pub f_radical: Vec<f64>,
    /// Dual descriptor (f+ - f-) per atom.
    pub dual_descriptor: Vec<f64>,
    /// Condensed atom-wise summary.
    pub condensed: Vec<CondensedFukuiAtom>,
    /// HOMO energy in eV.
    pub homo_energy: f64,
    /// LUMO energy in eV.
    pub lumo_energy: f64,
    /// HOMO-LUMO gap in eV.
    pub gap: f64,
    /// Domain and confidence caveats for interpretation.
    pub validity_notes: Vec<String>,
}

/// One ranked site entry for empirical local-reactivity helpers.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactivitySiteScore {
    /// Atom index in the input order.
    pub atom_index: usize,
    /// Composite empirical score.
    pub score: f64,
}

/// Empirical local-reactivity ranking summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactivityRanking {
    /// Ranked candidate sites for nucleophilic attack on the molecule.
    pub nucleophilic_attack_sites: Vec<ReactivitySiteScore>,
    /// Ranked candidate sites for electrophilic attack on the molecule.
    pub electrophilic_attack_sites: Vec<ReactivitySiteScore>,
    /// Ranked candidate sites for radical attack on the molecule.
    pub radical_attack_sites: Vec<ReactivitySiteScore>,
    /// Notes about how scores are constructed.
    pub notes: Vec<String>,
}

/// A single broad-band transition peak in an exploratory UV-Vis-like spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UvVisPeak {
    /// Transition energy in eV.
    pub energy_ev: f64,
    /// Wavelength in nm.
    pub wavelength_nm: f64,
    /// Relative intensity (arbitrary units).
    pub intensity: f64,
    /// Occupied MO index.
    pub from_mo: usize,
    /// Virtual MO index.
    pub to_mo: usize,
}

/// Exploratory UV-Vis-like spectrum built from low-cost EHT transitions.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UvVisSpectrum {
    /// Energy axis in eV.
    pub energies_ev: Vec<f64>,
    /// Intensity axis in arbitrary units.
    pub intensities: Vec<f64>,
    /// Dominant transitions used for interpretation.
    pub peaks: Vec<UvVisPeak>,
    /// Broadening sigma in eV.
    pub sigma: f64,
    /// Method caveats for use in UI/API consumers.
    pub notes: Vec<String>,
}

fn orbital_atom_contributions(
    basis: &[crate::eht::basis::AtomicOrbital],
    overlap: &nalgebra::DMatrix<f64>,
    coefficients: &[Vec<f64>],
    mo_index: usize,
    n_atoms: usize,
) -> Vec<f64> {
    let ao_to_atom = crate::population::population::ao_to_atom_map(basis);
    let mut contributions = vec![0.0; n_atoms];

    for mu in 0..basis.len() {
        let mut gross = 0.0;
        for nu in 0..basis.len() {
            gross += coefficients[mu][mo_index] * overlap[(mu, nu)] * coefficients[nu][mo_index];
        }
        contributions[ao_to_atom[mu]] += gross;
    }

    let total: f64 = contributions.iter().sum();
    if total.abs() > 1e-12 {
        for value in &mut contributions {
            *value /= total;
        }
    }

    contributions
}

/// Compute atom-resolved HOMO/LUMO descriptors from an EHT result.
pub fn compute_frontier_descriptors(
    elements: &[u8],
    positions: &[[f64; 3]],
    eht_result: &crate::eht::EhtResult,
) -> FrontierDescriptors {
    let basis = crate::eht::basis::build_basis(elements, positions);
    let overlap = crate::eht::build_overlap_matrix(&basis);

    let homo_atom_contributions = orbital_atom_contributions(
        &basis,
        &overlap,
        &eht_result.coefficients,
        eht_result.homo_index,
        elements.len(),
    );
    let lumo_atom_contributions = orbital_atom_contributions(
        &basis,
        &overlap,
        &eht_result.coefficients,
        eht_result.lumo_index,
        elements.len(),
    );
    let dual_descriptor = lumo_atom_contributions
        .iter()
        .zip(homo_atom_contributions.iter())
        .map(|(lumo, homo)| lumo - homo)
        .collect();

    FrontierDescriptors {
        num_atoms: elements.len(),
        homo_atom_contributions,
        lumo_atom_contributions,
        dual_descriptor,
        homo_energy: eht_result.homo_energy,
        lumo_energy: eht_result.lumo_energy,
        gap: eht_result.gap,
    }
}

fn validity_notes(elements: &[u8]) -> Vec<String> {
    let support = crate::eht::analyze_eht_support(elements);
    let mut notes = vec![
        "Fukui descriptors are computed from frontier-orbital atom contributions (EHT proxy), not from full finite-difference electron-addition/removal calculations.".to_string(),
        "Interpret values comparatively within related structures; absolute values are semi-quantitative.".to_string(),
    ];

    match support.level {
        crate::eht::SupportLevel::Supported => {
            notes.push(
                "Element set is in supported EHT domain; qualitative organic trend interpretation is usually reliable."
                    .to_string(),
            );
        }
        crate::eht::SupportLevel::Experimental => {
            notes.push(
                "Element set is in experimental EHT domain (typically transition metals); use rankings as exploratory guidance only."
                    .to_string(),
            );
        }
        crate::eht::SupportLevel::Unsupported => {
            notes.push(
                "Element set is outside supported EHT parameterization; descriptor reliability is low or undefined."
                    .to_string(),
            );
        }
    }

    notes.extend(support.warnings);
    notes
}

/// Compute Fukui-style local reactivity descriptors from EHT frontier orbitals.
pub fn compute_fukui_descriptors(
    elements: &[u8],
    positions: &[[f64; 3]],
    eht_result: &crate::eht::EhtResult,
) -> FukuiDescriptors {
    let frontier = compute_frontier_descriptors(elements, positions, eht_result);
    let f_plus = frontier.lumo_atom_contributions.clone();
    let f_minus = frontier.homo_atom_contributions.clone();
    let f_radical: Vec<f64> = f_plus
        .iter()
        .zip(f_minus.iter())
        .map(|(fp, fm)| 0.5 * (fp + fm))
        .collect();
    let dual_descriptor: Vec<f64> = f_plus
        .iter()
        .zip(f_minus.iter())
        .map(|(fp, fm)| fp - fm)
        .collect();

    let condensed: Vec<CondensedFukuiAtom> = (0..elements.len())
        .map(|idx| CondensedFukuiAtom {
            atom_index: idx,
            f_plus: f_plus[idx],
            f_minus: f_minus[idx],
            f_radical: f_radical[idx],
            dual_descriptor: dual_descriptor[idx],
        })
        .collect();

    FukuiDescriptors {
        num_atoms: elements.len(),
        f_plus,
        f_minus,
        f_radical,
        dual_descriptor,
        condensed,
        homo_energy: frontier.homo_energy,
        lumo_energy: frontier.lumo_energy,
        gap: frontier.gap,
        validity_notes: validity_notes(elements),
    }
}

fn sorted_scores(mut scores: Vec<ReactivitySiteScore>) -> Vec<ReactivitySiteScore> {
    scores.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    scores
}

/// Build empirical local-reactivity rankings from Fukui descriptors and Mulliken charges.
pub fn rank_reactivity_sites(
    fukui: &FukuiDescriptors,
    mulliken_charges: &[f64],
) -> ReactivityRanking {
    let n = fukui.num_atoms.min(mulliken_charges.len());

    let mut nucleophilic_attack_sites = Vec::with_capacity(n);
    let mut electrophilic_attack_sites = Vec::with_capacity(n);
    let mut radical_attack_sites = Vec::with_capacity(n);

    for i in 0..n {
        let q = mulliken_charges[i];
        let fp = fukui.f_plus[i];
        let fm = fukui.f_minus[i];
        let f0 = fukui.f_radical[i];

        // Positive charge reinforces nucleophilic attack susceptibility.
        let nuc_score = fp + 0.25 * q.max(0.0);
        // Negative charge reinforces electrophilic attack susceptibility.
        let elec_score = fm + 0.25 * (-q).max(0.0);
        let rad_score = f0 + 0.1 * q.abs();

        nucleophilic_attack_sites.push(ReactivitySiteScore {
            atom_index: i,
            score: nuc_score,
        });
        electrophilic_attack_sites.push(ReactivitySiteScore {
            atom_index: i,
            score: elec_score,
        });
        radical_attack_sites.push(ReactivitySiteScore {
            atom_index: i,
            score: rad_score,
        });
    }

    ReactivityRanking {
        nucleophilic_attack_sites: sorted_scores(nucleophilic_attack_sites),
        electrophilic_attack_sites: sorted_scores(electrophilic_attack_sites),
        radical_attack_sites: sorted_scores(radical_attack_sites),
        notes: vec![
            "Scores are empirical composites of condensed Fukui terms and Mulliken charge bias.".to_string(),
            "Use ranking order for exploratory prioritization; values are not calibrated kinetic barriers.".to_string(),
        ],
    }
}

fn gaussian(x: f64, mu: f64, sigma: f64) -> f64 {
    let s = sigma.max(1e-6);
    let norm = 1.0 / (s * (2.0 * std::f64::consts::PI).sqrt());
    let dx = x - mu;
    norm * (-0.5 * dx * dx / (s * s)).exp()
}

/// Build an exploratory UV-Vis-like spectrum from low-cost EHT occupied→virtual transitions.
pub fn compute_uv_vis_like_spectrum(
    eht_result: &crate::eht::EhtResult,
    sigma: f64,
    e_min: f64,
    e_max: f64,
    n_points: usize,
) -> UvVisSpectrum {
    let n_points = n_points.max(2);
    let span = (e_max - e_min).max(1e-6);
    let step = span / (n_points as f64 - 1.0);
    let energies_ev: Vec<f64> = (0..n_points).map(|i| e_min + step * i as f64).collect();
    let mut intensities = vec![0.0; n_points];

    let mut peaks = Vec::new();
    for occ in 0..=eht_result.homo_index {
        for virt in eht_result.lumo_index..eht_result.energies.len() {
            let delta_e = eht_result.energies[virt] - eht_result.energies[occ];
            if delta_e <= 1e-6 {
                continue;
            }

            let mut overlap_strength = 0.0;
            for ao in 0..eht_result.coefficients.len() {
                overlap_strength +=
                    (eht_result.coefficients[ao][occ] * eht_result.coefficients[ao][virt]).abs();
            }

            let intensity = overlap_strength * overlap_strength;
            if intensity <= 1e-12 {
                continue;
            }

            if peaks.len() < 24 {
                peaks.push(UvVisPeak {
                    energy_ev: delta_e,
                    wavelength_nm: 1239.841984 / delta_e,
                    intensity,
                    from_mo: occ,
                    to_mo: virt,
                });
            }

            for (idx, e) in energies_ev.iter().enumerate() {
                intensities[idx] += intensity * gaussian(*e, delta_e, sigma);
            }
        }
    }

    peaks.sort_by(|a, b| {
        b.intensity
            .partial_cmp(&a.intensity)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    UvVisSpectrum {
        energies_ev,
        intensities,
        peaks,
        sigma,
        notes: vec![
            "Exploratory UV-Vis-like spectrum from EHT MO energy differences and coefficient-overlap intensity proxy.".to_string(),
            "This is a qualitative visualization aid, not a calibrated excited-state method (no CI/TDDFT).".to_string(),
        ],
    }
}

fn has_bond_order(
    mol: &crate::graph::Molecule,
    a: usize,
    b: usize,
    order: crate::graph::BondOrder,
) -> bool {
    let ia = petgraph::graph::NodeIndex::new(a);
    let ib = petgraph::graph::NodeIndex::new(b);
    if let Some(edge_idx) = mol.graph.find_edge(ia, ib) {
        return mol.graph[edge_idx].order == order;
    }
    false
}

fn is_carboxylic_acid_oxygen(mol: &crate::graph::Molecule, atom_idx: usize) -> bool {
    let idx = petgraph::graph::NodeIndex::new(atom_idx);
    let atom = &mol.graph[idx];
    if atom.element != 8 {
        return false;
    }

    let neighbors: Vec<_> = mol.graph.neighbors(idx).collect();
    let has_h = neighbors.iter().any(|n| mol.graph[*n].element == 1);
    if !has_h {
        return false;
    }

    let carbon_neighbor = neighbors.iter().find(|n| mol.graph[**n].element == 6);
    if let Some(c_idx) = carbon_neighbor {
        let c_neighbors: Vec<_> = mol.graph.neighbors(*c_idx).collect();
        let carbonyl_o_count = c_neighbors
            .iter()
            .filter(|n| mol.graph[**n].element == 8)
            .filter(|n| {
                has_bond_order(
                    mol,
                    c_idx.index(),
                    n.index(),
                    crate::graph::BondOrder::Double,
                )
            })
            .count();
        return carbonyl_o_count >= 1;
    }

    false
}

fn is_phenol_oxygen(mol: &crate::graph::Molecule, atom_idx: usize) -> bool {
    let idx = petgraph::graph::NodeIndex::new(atom_idx);
    let atom = &mol.graph[idx];
    if atom.element != 8 {
        return false;
    }
    let neighbors: Vec<_> = mol.graph.neighbors(idx).collect();
    let has_h = neighbors.iter().any(|n| mol.graph[*n].element == 1);
    if !has_h {
        return false;
    }
    neighbors.iter().any(|n| {
        if mol.graph[*n].element != 6 {
            return false;
        }
        mol.graph
            .edges(*n)
            .any(|e| matches!(e.weight().order, crate::graph::BondOrder::Aromatic))
    })
}

fn is_thiol_sulfur(mol: &crate::graph::Molecule, atom_idx: usize) -> bool {
    let idx = petgraph::graph::NodeIndex::new(atom_idx);
    let atom = &mol.graph[idx];
    if atom.element != 16 {
        return false;
    }
    mol.graph.neighbors(idx).any(|n| mol.graph[n].element == 1)
}

fn is_neutral_amine_nitrogen(mol: &crate::graph::Molecule, atom_idx: usize) -> bool {
    let idx = petgraph::graph::NodeIndex::new(atom_idx);
    let atom = &mol.graph[idx];
    if atom.element != 7 || atom.formal_charge != 0 {
        return false;
    }
    if !matches!(atom.hybridization, crate::graph::Hybridization::SP3) {
        return false;
    }
    !mol.graph
        .edges(idx)
        .any(|e| matches!(e.weight().order, crate::graph::BondOrder::Aromatic))
}

fn is_aromatic_nitrogen(mol: &crate::graph::Molecule, atom_idx: usize) -> bool {
    let idx = petgraph::graph::NodeIndex::new(atom_idx);
    let atom = &mol.graph[idx];
    if atom.element != 7 {
        return false;
    }
    mol.graph
        .edges(idx)
        .any(|e| matches!(e.weight().order, crate::graph::BondOrder::Aromatic))
}

/// Estimate acidic/basic pKa sites from graph environment and Gasteiger charges.
pub fn estimate_empirical_pka(mol: &crate::graph::Molecule, charges: &[f64]) -> EmpiricalPkaResult {
    let n = mol.graph.node_count().min(charges.len());
    let mut acidic_sites = Vec::new();
    let mut basic_sites = Vec::new();

    for atom_idx in 0..n {
        let idx = petgraph::graph::NodeIndex::new(atom_idx);
        let atom = &mol.graph[idx];
        let q = charges[atom_idx];

        if is_carboxylic_acid_oxygen(mol, atom_idx) {
            acidic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "acidic".to_string(),
                environment: "carboxylic_acid_oxygen".to_string(),
                estimated_pka: (4.5 - 2.0 * q).clamp(-1.0, 14.0),
                confidence: 0.82,
            });
        } else if is_phenol_oxygen(mol, atom_idx) {
            acidic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "acidic".to_string(),
                environment: "phenol_oxygen".to_string(),
                estimated_pka: (10.0 - 1.5 * q).clamp(2.0, 16.0),
                confidence: 0.68,
            });
        } else if is_thiol_sulfur(mol, atom_idx) {
            acidic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "acidic".to_string(),
                environment: "thiol_sulfur".to_string(),
                estimated_pka: (10.5 - 1.2 * q).clamp(2.0, 16.0),
                confidence: 0.64,
            });
        } else if atom.element == 7 && atom.formal_charge > 0 {
            acidic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "acidic".to_string(),
                environment: "ammonium_like".to_string(),
                estimated_pka: (9.7 - 1.0 * q).clamp(4.0, 14.0),
                confidence: 0.62,
            });
        }

        if is_neutral_amine_nitrogen(mol, atom_idx) {
            let h_count = mol
                .graph
                .neighbors(idx)
                .filter(|n| mol.graph[*n].element == 1)
                .count();
            let base_pka = if h_count >= 2 {
                10.8
            } else if h_count == 1 {
                10.4
            } else {
                9.8
            };
            basic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "basic".to_string(),
                environment: "aliphatic_amine".to_string(),
                estimated_pka: (base_pka - 2.5 * q).clamp(2.0, 14.5),
                confidence: 0.75,
            });
        } else if is_aromatic_nitrogen(mol, atom_idx) && atom.formal_charge <= 0 {
            basic_sites.push(EmpiricalPkaSite {
                atom_index: atom_idx,
                site_type: "basic".to_string(),
                environment: "aromatic_nitrogen".to_string(),
                estimated_pka: (5.2 - 1.8 * q).clamp(-1.0, 10.0),
                confidence: 0.6,
            });
        }
    }

    acidic_sites.sort_by(|a, b| {
        a.estimated_pka
            .partial_cmp(&b.estimated_pka)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    basic_sites.sort_by(|a, b| {
        b.estimated_pka
            .partial_cmp(&a.estimated_pka)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    EmpiricalPkaResult {
        acidic_sites,
        basic_sites,
        notes: vec![
            "Empirical pKa estimates use simple graph environments plus Gasteiger-charge adjustments; values are coarse screening hints, not publication-grade predictions.".to_string(),
            "Best use is relative ranking within related congeneric series under similar conditions.".to_string(),
        ],
    }
}

/// Apply a lightweight aromatic stabilization correction on top of raw UFF energy.
pub fn apply_aromatic_uff_correction(
    mol: &crate::graph::Molecule,
    raw_energy_kcal_mol: f64,
) -> UffHeuristicEnergy {
    let aromatic_bond_count = mol
        .graph
        .edge_references()
        .filter(|e| matches!(e.weight().order, crate::graph::BondOrder::Aromatic))
        .count();
    let aromatic_stabilization_kcal_mol = -0.08 * aromatic_bond_count as f64;

    UffHeuristicEnergy {
        raw_energy_kcal_mol,
        aromatic_stabilization_kcal_mol,
        corrected_energy_kcal_mol: raw_energy_kcal_mol + aromatic_stabilization_kcal_mol,
        aromatic_bond_count,
        notes: vec![
            "Aromatic correction is an empirical post-UFF heuristic tied to aromatic bond count and should be used for ranking guidance only.".to_string(),
            "Raw UFF and corrected values are both reported so downstream workflows can choose their own policy.".to_string(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_frontier_contributions_are_symmetric() {
        let elements = [1u8, 1];
        let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
        let eht = crate::eht::solve_eht(&elements, &positions, None).unwrap();
        let descriptors = compute_frontier_descriptors(&elements, &positions, &eht);

        assert!(
            (descriptors.homo_atom_contributions[0] - descriptors.homo_atom_contributions[1]).abs()
                < 1e-6
        );
        assert!(
            (descriptors.lumo_atom_contributions[0] - descriptors.lumo_atom_contributions[1]).abs()
                < 1e-6
        );
    }

    #[test]
    fn test_frontier_contributions_are_normalized() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let eht = crate::eht::solve_eht(&elements, &positions, None).unwrap();
        let descriptors = compute_frontier_descriptors(&elements, &positions, &eht);

        let homo_sum: f64 = descriptors.homo_atom_contributions.iter().sum();
        let lumo_sum: f64 = descriptors.lumo_atom_contributions.iter().sum();
        assert!((homo_sum - 1.0).abs() < 1e-6);
        assert!((lumo_sum - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_fukui_descriptors_are_consistent_and_normalized() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let eht = crate::eht::solve_eht(&elements, &positions, None).unwrap();
        let fukui = compute_fukui_descriptors(&elements, &positions, &eht);

        let f_plus_sum: f64 = fukui.f_plus.iter().sum();
        let f_minus_sum: f64 = fukui.f_minus.iter().sum();
        let f0_sum: f64 = fukui.f_radical.iter().sum();
        assert!((f_plus_sum - 1.0).abs() < 1e-6);
        assert!((f_minus_sum - 1.0).abs() < 1e-6);
        assert!((f0_sum - 1.0).abs() < 1e-6);
        assert_eq!(fukui.condensed.len(), elements.len());
    }

    #[test]
    fn test_uv_vis_like_spectrum_has_requested_grid() {
        let elements = [6u8, 6, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [1.34, 0.0, 0.0],
            [-0.6, 0.92, 0.0],
            [-0.6, -0.92, 0.0],
            [1.94, 0.92, 0.0],
            [1.94, -0.92, 0.0],
        ];
        let eht = crate::eht::solve_eht(&elements, &positions, None).unwrap();
        let spec = compute_uv_vis_like_spectrum(&eht, 0.2, 0.5, 8.0, 300);
        assert_eq!(spec.energies_ev.len(), 300);
        assert_eq!(spec.intensities.len(), 300);
    }

    #[test]
    fn test_empirical_pka_detects_carboxylic_acid_site() {
        let mol = crate::graph::Molecule::from_smiles("CC(=O)O").unwrap();
        let charges = crate::compute_charges("CC(=O)O").unwrap().charges;
        let result = estimate_empirical_pka(&mol, &charges);

        assert!(!result.acidic_sites.is_empty());
        assert!(result
            .acidic_sites
            .iter()
            .any(|site| site.environment == "carboxylic_acid_oxygen"));
    }

    #[test]
    fn test_aromatic_uff_correction_is_negative_for_benzene() {
        let mol = crate::graph::Molecule::from_smiles("c1ccccc1").unwrap();
        let result = apply_aromatic_uff_correction(&mol, 10.0);
        assert!(result.aromatic_bond_count >= 6);
        assert!(result.aromatic_stabilization_kcal_mol < 0.0);
        assert!(result.corrected_energy_kcal_mol < result.raw_energy_kcal_mol);
    }
}
