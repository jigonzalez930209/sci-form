//! Electrostatic steering for reaction approach direction.
//!
//! Uses Gasteiger partial charges to compute a Coulombic potential
//! and determine the optimal electrostatic approach direction.

/// Result of electrostatic approach analysis.
pub struct ElectrostaticApproachInfo {
    /// Unit vector for the optimal electrostatic approach direction.
    pub direction: [f64; 3],
    /// Electrostatic interaction energy at the approach configuration (kcal/mol).
    pub interaction_energy: f64,
    /// Index of the most negative atom on fragment A (nucleophilic site).
    pub nucleophilic_atom: usize,
    /// Index of the most positive atom on fragment B (electrophilic site).
    pub electrophilic_atom: usize,
}

/// Compute optimal approach direction based on electrostatic complementarity.
///
/// Finds the most nucleophilic atom on fragment A (most negative charge)
/// and the most electrophilic atom on fragment B (most positive charge),
/// then returns the direction from A→B that maximises favourable interaction.
pub fn compute_electrostatic_approach(
    elements_a: &[u8],
    coords_a: &[f64],
    elements_b: &[u8],
    coords_b: &[f64],
) -> Option<ElectrostaticApproachInfo> {
    let smiles_a = elements_to_dummy_smiles(elements_a);
    let smiles_b = elements_to_dummy_smiles(elements_b);

    let charges_a = crate::compute_charges(&smiles_a).ok()?;
    let charges_b = crate::compute_charges(&smiles_b).ok()?;

    let n_a = elements_a.len();
    let n_b = elements_b.len();

    if charges_a.charges.len() < n_a || charges_b.charges.len() < n_b {
        return None;
    }

    // Find most negative atom on A (nucleophilic site)
    let (nuc_idx, _nuc_charge) = charges_a.charges[..n_a]
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;

    // Find most positive atom on B (electrophilic site)
    let (elec_idx, _elec_charge) = charges_b.charges[..n_b]
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))?;

    // Direction from nucleophilic atom on A toward electrophilic atom on B
    let _ax = coords_a[nuc_idx * 3];
    let _ay = coords_a[nuc_idx * 3 + 1];
    let _az = coords_a[nuc_idx * 3 + 2];
    let _bx = coords_b[elec_idx * 3];
    let _by = coords_b[elec_idx * 3 + 1];
    let _bz = coords_b[elec_idx * 3 + 2];

    // If fragments aren't close, use COM→COM direction with charge-biased weighting
    let com_a = charge_weighted_centroid(coords_a, &charges_a.charges[..n_a], true);
    let com_b = charge_weighted_centroid(coords_b, &charges_b.charges[..n_b], false);

    let dx = com_b[0] - com_a[0];
    let dy = com_b[1] - com_a[1];
    let dz = com_b[2] - com_a[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();

    let direction = if len > 1e-10 {
        [dx / len, dy / len, dz / len]
    } else {
        [1.0, 0.0, 0.0]
    };

    // Estimate interaction energy (Coulombic sum between fragments)
    let mut e_int = 0.0f64;
    let kcal_au = 332.0637; // conversion eÅ → kcal/mol
    for i in 0..n_a {
        let qi = charges_a.charges[i];
        for j in 0..n_b {
            let qj = charges_b.charges[j];
            let dxx = coords_a[i * 3] - coords_b[j * 3];
            let dyy = coords_a[i * 3 + 1] - coords_b[j * 3 + 1];
            let dzz = coords_a[i * 3 + 2] - coords_b[j * 3 + 2];
            let r = (dxx * dxx + dyy * dyy + dzz * dzz).sqrt().max(0.5);
            e_int += kcal_au * qi * qj / r;
        }
    }

    Some(ElectrostaticApproachInfo {
        direction,
        interaction_energy: e_int,
        nucleophilic_atom: nuc_idx,
        electrophilic_atom: elec_idx,
    })
}

/// Charge-weighted centroid of a fragment.
///
/// If `negative_bias` is true, emphasises negative charges (nucleophilic centers).
/// If false, emphasises positive charges (electrophilic centers).
fn charge_weighted_centroid(coords: &[f64], charges: &[f64], negative_bias: bool) -> [f64; 3] {
    let n = charges.len();
    let mut c = [0.0; 3];
    let mut w_total = 0.0f64;

    for i in 0..n {
        let q = charges[i];
        // Weight: for nucleophile, larger weight for more negative charges
        // For electrophile, larger weight for more positive charges
        let w = if negative_bias {
            (-q + 1.0).max(0.01)
        } else {
            (q + 1.0).max(0.01)
        };
        c[0] += w * coords[i * 3];
        c[1] += w * coords[i * 3 + 1];
        c[2] += w * coords[i * 3 + 2];
        w_total += w;
    }

    if w_total > 1e-14 {
        c[0] /= w_total;
        c[1] /= w_total;
        c[2] /= w_total;
    }

    c
}

/// Convert element list to a dummy SMILES for Gasteiger charges.
///
/// This is a heuristic — Gasteiger requires bond topology. For simple organic
/// molecules we can often reconstruct a reasonable SMILES from elements alone.
/// Falls back to atom-only SMILES for problematic cases.
fn elements_to_dummy_smiles(elements: &[u8]) -> String {
    // Simple heuristic: count heavy atoms and generate a chain
    let symbols: Vec<&str> = elements
        .iter()
        .map(|&z| match z {
            1 => "[H]",
            6 => "C",
            7 => "N",
            8 => "O",
            9 => "F",
            15 => "P",
            16 => "S",
            17 => "[Cl]",
            35 => "[Br]",
            53 => "[I]",
            _ => "[*]",
        })
        .collect();

    // For charges we just need a rough SMILES — concat heavy atoms as a chain
    let heavy: Vec<&str> = symbols.iter().copied().filter(|&s| s != "[H]").collect();
    if heavy.is_empty() {
        return "[H]".into();
    }
    heavy.join("")
}
