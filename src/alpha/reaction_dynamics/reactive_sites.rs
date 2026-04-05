//! Reactive site identification for 3D reaction dynamics.
//!
//! This module answers the two critical questions:
//! 1. **Which atoms are the reactive sites?**
//!    Uses multiple information sources ranked by reliability:
//!    a) SMIRKS atom mapping (explicit, highest confidence)
//!    b) Fukui/frontier MO descriptors (quantum, medium confidence)
//!    c) pKa/charge heuristics (empirical, lower confidence)
//!    d) Proximity heuristics (geometric, fallback)
//!
//! 2. **How should molecules approach?**
//!    Once reactive atoms are identified, computes the optimal approach
//!    vector from the 3D positions of reactive sites.
//!
//! This module bridges the gap between SMIRKS (graph-level, no 3D) and
//! the reaction dynamics pipeline (needs 3D approach directions).

use serde::{Deserialize, Serialize};

/// A reactive site on a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactiveSite {
    /// Atom index within the fragment.
    pub atom_index: usize,
    /// Confidence score [0, 1] — how certain we are this is the reactive atom.
    pub confidence: f64,
    /// Source of identification.
    pub source: SiteSource,
    /// Role in reaction: nucleophile, electrophile, radical, or bond-forming/breaking.
    pub role: SiteRole,
    /// 3D position [x, y, z] (Å).
    pub position: [f64; 3],
}

/// How the reactive site was identified.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SiteSource {
    /// Explicit SMIRKS atom mapping (highest confidence).
    SmirksMapping,
    /// Fukui/frontier MO analysis.
    FukuiFrontier,
    /// Charge and pKa heuristics.
    ChargeHeuristic,
    /// Geometric proximity fallback.
    ProximityFallback,
}

/// Role of the reactive site.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SiteRole {
    /// Nucleophilic site (electron donor, high f⁻).
    Nucleophile,
    /// Electrophilic site (electron acceptor, high f⁺).
    Electrophile,
    /// Radical site (high f_radical).
    Radical,
    /// Bond-forming site (from SMIRKS).
    BondForming,
    /// Bond-breaking site (from SMIRKS).
    BondBreaking,
    /// Unclassified reactive site.
    Unclassified,
}

/// Complete reactive site analysis for a set of fragments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactiveSiteAnalysis {
    /// Reactive sites per fragment, indexed by fragment index.
    pub fragment_sites: Vec<Vec<ReactiveSite>>,
    /// Optimal approach direction (unit vector), from fragment 0 to fragment 1.
    pub approach_direction: Option<[f64; 3]>,
    /// Distance between the two primary reactive sites (Å).
    pub reactive_distance: Option<f64>,
    /// SMIRKS-derived bond changes (if available).
    pub bond_changes: Vec<BondChangeInfo>,
    /// Diagnostic notes.
    pub notes: Vec<String>,
}

/// Bond change information with 3D context.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BondChangeInfo {
    /// Fragment index of first atom.
    pub frag_a: usize,
    /// Atom index within fragment A.
    pub atom_a: usize,
    /// Fragment index of second atom.
    pub frag_b: usize,
    /// Atom index within fragment B.
    pub atom_b: usize,
    /// Whether this bond is being formed (true) or broken (false).
    pub forming: bool,
}

/// Identify reactive sites using all available information.
///
/// Priority chain:
/// 1. If `smirks` is provided → use explicit atom mapping
/// 2. Always compute Fukui descriptors for electronic insight
/// 3. Fall back to proximity/charge heuristics
///
/// The result includes approach directions derived from the identified sites.
pub fn identify_reactive_sites(
    fragment_elements: &[&[u8]],
    fragment_coords: &[&[f64]],
    fragment_smiles: &[&str],
    smirks: Option<&str>,
) -> Result<ReactiveSiteAnalysis, String> {
    let n_frags = fragment_elements.len();
    if n_frags == 0 {
        return Err("No fragments provided".into());
    }

    let mut fragment_sites: Vec<Vec<ReactiveSite>> = vec![Vec::new(); n_frags];
    let mut bond_changes = Vec::new();
    let mut notes = Vec::new();

    // ── 1. SMIRKS-based identification (highest confidence) ──────────
    if let Some(smirks_str) = smirks {
        if let Ok(transform) = crate::smirks::parse_smirks(smirks_str) {
            notes.push(format!(
                "SMIRKS parsed: {} atom maps, {} bond changes",
                transform.atom_map.len(),
                transform.bond_changes.len()
            ));

            // Try to match SMIRKS against each fragment
            let smirks_sites = identify_sites_from_smirks(
                &transform,
                fragment_elements,
                fragment_coords,
                fragment_smiles,
            );

            if let Ok((sites, changes)) = smirks_sites {
                for (frag_idx, frag_sites) in sites.into_iter().enumerate() {
                    if frag_idx < n_frags {
                        fragment_sites[frag_idx].extend(frag_sites);
                    }
                }
                bond_changes = changes;
                notes.push("Reactive sites identified from SMIRKS atom mapping.".into());
            } else {
                notes.push("SMIRKS matching failed, falling back to Fukui analysis.".into());
            }
        }
    }

    // ── 2. Fukui/frontier MO descriptors (always, for electronic insight) ──
    for (frag_idx, (elems, coords)) in fragment_elements
        .iter()
        .zip(fragment_coords.iter())
        .enumerate()
    {
        let n_atoms = elems.len();
        let pos: Vec<[f64; 3]> = coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        // Only add Fukui sites if SMIRKS didn't already identify sites for this fragment
        let has_smirks_sites = !fragment_sites[frag_idx].is_empty();

        if let Ok(ranking) = crate::compute_reactivity_ranking(elems, &pos) {
            // Top nucleophilic site
            if let Some(nuc_site) = ranking.nucleophilic_attack_sites.first() {
                let a = nuc_site.atom_index;
                if a < n_atoms {
                    let site = ReactiveSite {
                        atom_index: a,
                        confidence: if has_smirks_sites { 0.5 } else { 0.8 },
                        source: SiteSource::FukuiFrontier,
                        role: SiteRole::Electrophile, // high f⁺ = attacked by nucleophile
                        position: pos[a],
                    };
                    if !has_smirks_sites {
                        fragment_sites[frag_idx].push(site);
                    }
                }
            }

            // Top electrophilic site
            if let Some(elec_site) = ranking.electrophilic_attack_sites.first() {
                let a = elec_site.atom_index;
                if a < n_atoms {
                    let site = ReactiveSite {
                        atom_index: a,
                        confidence: if has_smirks_sites { 0.5 } else { 0.8 },
                        source: SiteSource::FukuiFrontier,
                        role: SiteRole::Nucleophile, // high f⁻ = attacked by electrophile
                        position: pos[a],
                    };
                    if !has_smirks_sites {
                        fragment_sites[frag_idx].push(site);
                    }
                }
            }

            notes.push(format!(
                "Fragment {}: Fukui reactivity ranking computed.",
                frag_idx
            ));
        }
    }

    // ── 3. Fallback: proximity/charge heuristics ──────────────────────
    for (frag_idx, (elems, coords)) in fragment_elements
        .iter()
        .zip(fragment_coords.iter())
        .enumerate()
    {
        if !fragment_sites[frag_idx].is_empty() {
            continue; // Already have sites
        }

        let n_atoms = elems.len();

        // Use Gasteiger charges to find the most polar atom
        if let Ok(charge_result) = crate::compute_charges(fragment_smiles[frag_idx]) {
            let charges = &charge_result.charges;
            // Most negative atom → nucleophile
            if let Some((idx, _)) = charges
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            {
                if idx < n_atoms {
                    let pos: [f64; 3] = [coords[idx * 3], coords[idx * 3 + 1], coords[idx * 3 + 2]];
                    fragment_sites[frag_idx].push(ReactiveSite {
                        atom_index: idx,
                        confidence: 0.4,
                        source: SiteSource::ChargeHeuristic,
                        role: SiteRole::Nucleophile,
                        position: pos,
                    });
                }
            }
            // Most positive atom → electrophile
            if let Some((idx, _)) = charges
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            {
                if idx < n_atoms {
                    let pos: [f64; 3] = [coords[idx * 3], coords[idx * 3 + 1], coords[idx * 3 + 2]];
                    fragment_sites[frag_idx].push(ReactiveSite {
                        atom_index: idx,
                        confidence: 0.4,
                        source: SiteSource::ChargeHeuristic,
                        role: SiteRole::Electrophile,
                        position: pos,
                    });
                }
            }
        }

        // Ultimate fallback: use atom with highest atomic number (heaviest → most reactive)
        if fragment_sites[frag_idx].is_empty() && n_atoms > 0 {
            let (best_idx, _) = elems.iter().enumerate().max_by_key(|(_, &z)| z).unwrap();
            let pos: [f64; 3] = [
                coords[best_idx * 3],
                coords[best_idx * 3 + 1],
                coords[best_idx * 3 + 2],
            ];
            fragment_sites[frag_idx].push(ReactiveSite {
                atom_index: best_idx,
                confidence: 0.2,
                source: SiteSource::ProximityFallback,
                role: SiteRole::Unclassified,
                position: pos,
            });
            notes.push(format!(
                "Fragment {}: fallback to heaviest atom (Z={}).",
                frag_idx, elems[best_idx]
            ));
        }
    }

    // ── 4. Compute approach direction from reactive sites ─────────────
    let (approach_direction, reactive_distance) = if n_frags >= 2 {
        compute_approach_from_sites(&fragment_sites)
    } else {
        (None, None)
    };

    Ok(ReactiveSiteAnalysis {
        fragment_sites,
        approach_direction,
        reactive_distance,
        bond_changes,
        notes,
    })
}

/// Use SMIRKS atom mapping to identify reactive sites with 3D positions.
fn identify_sites_from_smirks(
    transform: &crate::smirks::SmirksTransform,
    fragment_elements: &[&[u8]],
    fragment_coords: &[&[f64]],
    fragment_smiles: &[&str],
) -> Result<(Vec<Vec<ReactiveSite>>, Vec<BondChangeInfo>), String> {
    let n_frags = fragment_elements.len();
    let mut fragment_sites: Vec<Vec<ReactiveSite>> = vec![Vec::new(); n_frags];
    let mut bond_changes = Vec::new();

    // Try to match SMIRKS against fragments
    let smirks_result = if n_frags == 1 {
        crate::smirks::apply_smirks(&transform.smirks, fragment_smiles[0])
    } else {
        crate::smirks::apply_smirks_multi(&transform.smirks, fragment_smiles)
    };

    let result = smirks_result?;
    if !result.success {
        return Err("SMIRKS pattern did not match".into());
    }

    // atom_mapping: map_num → molecule_atom_index
    // We need to figure out which fragment each mapped atom belongs to
    let mut frag_offsets = Vec::with_capacity(n_frags);
    let mut offset = 0usize;
    for elems in fragment_elements {
        frag_offsets.push(offset);
        offset += elems.len();
    }

    for (&map_num, &mol_atom_idx) in &result.atom_mapping {
        // Determine which fragment this atom belongs to
        let (frag_idx, local_idx) = find_fragment_for_atom(mol_atom_idx, fragment_elements);

        if frag_idx < n_frags && local_idx < fragment_elements[frag_idx].len() {
            let coords = fragment_coords[frag_idx];
            let pos = [
                coords[local_idx * 3],
                coords[local_idx * 3 + 1],
                coords[local_idx * 3 + 2],
            ];

            // Determine role from bond changes
            let role = determine_role_from_smirks(map_num, &transform.bond_changes);

            fragment_sites[frag_idx].push(ReactiveSite {
                atom_index: local_idx,
                confidence: 0.95,
                source: SiteSource::SmirksMapping,
                role,
                position: pos,
            });
        }
    }

    // Build 3D bond change info
    for bc in &transform.bond_changes {
        if let (Some(&atom_a_mol), Some(&atom_b_mol)) = (
            result.atom_mapping.get(&bc.atom1_map),
            result.atom_mapping.get(&bc.atom2_map),
        ) {
            let (frag_a, local_a) = find_fragment_for_atom(atom_a_mol, fragment_elements);
            let (frag_b, local_b) = find_fragment_for_atom(atom_b_mol, fragment_elements);
            let forming = bc.old_order.is_none() && bc.new_order.is_some();
            bond_changes.push(BondChangeInfo {
                frag_a,
                atom_a: local_a,
                frag_b,
                atom_b: local_b,
                forming,
            });
        }
    }

    Ok((fragment_sites, bond_changes))
}

/// Find which fragment an atom index belongs to (multi-fragment system).
fn find_fragment_for_atom(atom_idx: usize, fragment_elements: &[&[u8]]) -> (usize, usize) {
    let mut offset = 0;
    for (frag_idx, elems) in fragment_elements.iter().enumerate() {
        if atom_idx < offset + elems.len() {
            return (frag_idx, atom_idx - offset);
        }
        offset += elems.len();
    }
    // Fallback: return last fragment
    let last = fragment_elements.len().saturating_sub(1);
    (last, atom_idx.saturating_sub(offset))
}

/// Determine the role of a mapped atom from SMIRKS bond changes.
fn determine_role_from_smirks(
    map_num: usize,
    bond_changes: &[crate::smirks::BondChange],
) -> SiteRole {
    for bc in bond_changes {
        if bc.atom1_map == map_num || bc.atom2_map == map_num {
            if bc.old_order.is_none() && bc.new_order.is_some() {
                return SiteRole::BondForming;
            }
            if bc.old_order.is_some() && bc.new_order.is_none() {
                return SiteRole::BondBreaking;
            }
        }
    }
    SiteRole::Unclassified
}

/// Compute optimal approach direction from identified reactive sites on two fragments.
///
/// Strategy:
/// - Find the highest-confidence site on each fragment
/// - Approach direction = vector from fragment 0's site to fragment 1's site
/// - If multiple sites, use the pair with the most complementary roles
///   (nucleophile → electrophile preferred)
fn compute_approach_from_sites(
    fragment_sites: &[Vec<ReactiveSite>],
) -> (Option<[f64; 3]>, Option<f64>) {
    if fragment_sites.len() < 2 {
        return (None, None);
    }

    let sites_a = &fragment_sites[0];
    let sites_b = &fragment_sites[1];

    if sites_a.is_empty() || sites_b.is_empty() {
        return (None, None);
    }

    // Score each (site_a, site_b) pair by:
    // 1. Complementary roles (nuc→elec or elec→nuc)
    // 2. Sum of confidence scores
    let mut best_score = -1.0f64;
    let mut best_a = 0;
    let mut best_b = 0;

    for (ia, sa) in sites_a.iter().enumerate() {
        for (ib, sb) in sites_b.iter().enumerate() {
            let mut score = sa.confidence + sb.confidence;

            // Bonus for complementary electronic roles
            score += role_complementarity(&sa.role, &sb.role);

            if score > best_score {
                best_score = score;
                best_a = ia;
                best_b = ib;
            }
        }
    }

    let pa = sites_a[best_a].position;
    let pb = sites_b[best_b].position;

    let dx = pb[0] - pa[0];
    let dy = pb[1] - pa[1];
    let dz = pb[2] - pa[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();

    if dist < 1e-6 {
        return (None, None);
    }

    let direction = [dx / dist, dy / dist, dz / dist];
    (Some(direction), Some(dist))
}

/// Score how complementary two reactive roles are.
fn role_complementarity(a: &SiteRole, b: &SiteRole) -> f64 {
    match (a, b) {
        (SiteRole::Nucleophile, SiteRole::Electrophile)
        | (SiteRole::Electrophile, SiteRole::Nucleophile) => 1.0,
        (SiteRole::BondForming, SiteRole::BondForming) => 0.8,
        (SiteRole::BondBreaking, SiteRole::BondBreaking) => 0.6,
        (SiteRole::Radical, SiteRole::Radical) => 0.5,
        _ => 0.0,
    }
}

/// Identify reactive atom pairs between two fragments for constrained optimization.
///
/// Returns (atom_in_frag_a, atom_in_frag_b) pairs that should form bonds,
/// derived from reactive site analysis.
pub fn reactive_atom_pairs(analysis: &ReactiveSiteAnalysis) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();

    // From bond change info (SMIRKS-derived, highest confidence)
    for bc in &analysis.bond_changes {
        if bc.forming && bc.frag_a == 0 && bc.frag_b == 1 {
            pairs.push((bc.atom_a, bc.atom_b));
        } else if bc.forming && bc.frag_a == 1 && bc.frag_b == 0 {
            pairs.push((bc.atom_b, bc.atom_a));
        }
    }

    // If no SMIRKS bond changes, use highest-confidence sites
    if pairs.is_empty() && analysis.fragment_sites.len() >= 2 {
        let sites_a = &analysis.fragment_sites[0];
        let sites_b = &analysis.fragment_sites[1];

        if let (Some(sa), Some(sb)) = (
            sites_a.iter().max_by(|a, b| {
                a.confidence
                    .partial_cmp(&b.confidence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
            sites_b.iter().max_by(|a, b| {
                a.confidence
                    .partial_cmp(&b.confidence)
                    .unwrap_or(std::cmp::Ordering::Equal)
            }),
        ) {
            pairs.push((sa.atom_index, sb.atom_index));
        }
    }

    pairs
}
