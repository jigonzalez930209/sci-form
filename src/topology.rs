//! Coordination-topology analysis of transition-metal centres.
//!
//! Given the atomic numbers and Cartesian coordinates of a molecule, this
//! module finds the transition-metal centres, determines which atoms are
//! coordinated to them from a covalent-radius criterion, and classifies the
//! resulting coordination geometry by comparing the ligand directions with
//! the idealised directions of the common polyhedra.
//!
//! The geometry classification is purely geometric: the metal–ligand
//! vectors are matched against the ideal directions of the linear,
//! trigonal, tetrahedral, square-planar, trigonal-bipyramidal and octahedral
//! polyhedra by minimising `1 − v̂ · î` over all assignments, and the winning
//! geometry is reported together with a dimensionless quality score.

use serde::{Deserialize, Serialize};

/// Best-fit coordination geometry of a metal centre.
///
/// The variant is selected by comparing the observed metal–ligand directions
/// with the ideal directions of the candidate polyhedra; it is a guess from
/// the geometry alone and carries no energetic information.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CoordinationGeometryGuess {
    /// Two-coordinate, 180° (a linear MCl₂ fragment).
    Linear,
    /// Three-coordinate, 120° planar.
    Trigonal,
    /// Four-coordinate, 109.5° tetrahedral (e.g. ZnCl₄²⁻).
    Tetrahedral,
    /// Four-coordinate, square planar with 90° angles (e.g. PtCl₄²⁻).
    SquarePlanar,
    /// Five-coordinate, trigonal bipyramidal (e.g. Fe(CO)₅).
    TrigonalBipyramidal,
    /// Six-coordinate, octahedral (e.g. W(CO)₆).
    Octahedral,
    /// Geometry could not be assigned (unsupported coordination number).
    Unknown,
}

/// A transition-metal centre found by [`analyze_topology`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetalCoordinationCenter {
    /// Index of the metal atom in the input arrays.
    pub atom_index: usize,
    /// Atomic number (Z) of the metal.
    pub element: u8,
    /// Indices of the atoms classified as ligands of this metal.
    pub ligand_indices: Vec<usize>,
    /// Number of coordinated ligands (the coordination number).
    pub coordination_number: usize,
    /// Best-fit coordination geometry of the ligand arrangement.
    pub geometry: CoordinationGeometryGuess,
    /// Dimensionless quality of the geometry fit, in the range 0.0 (poor) to
    /// 1.0 (perfect match with the ideal polyhedron).
    pub geometry_score: f64,
}

/// Result of a coordination-topology analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologyAnalysisResult {
    /// The transition-metal centres that were located.
    pub metal_centers: Vec<MetalCoordinationCenter>,
    /// Human-readable caveats, e.g. when no metal centre was found.
    pub warnings: Vec<String>,
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn normalize(v: [f64; 3]) -> [f64; 3] {
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag <= 1e-12 {
        [0.0, 0.0, 1.0]
    } else {
        [v[0] / mag, v[1] / mag, v[2] / mag]
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn is_likely_coordinated(metal_z: u8, ligand_z: u8, dist: f64) -> bool {
    if ligand_z == 1 {
        return false;
    }
    let cutoff = 1.3
        * (crate::graph::get_covalent_radius(metal_z)
            + crate::graph::get_covalent_radius(ligand_z))
        + 0.25;
    dist <= cutoff
}

fn ideal_directions(geometry: CoordinationGeometryGuess) -> &'static [[f64; 3]] {
    const LINEAR: [[f64; 3]; 2] = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    const TRIGONAL: [[f64; 3]; 3] = [
        [1.0, 0.0, 0.0],
        [-0.5, 0.8660254037844386, 0.0],
        [-0.5, -0.8660254037844386, 0.0],
    ];
    const TETRAHEDRAL: [[f64; 3]; 4] = [
        [0.5773502691896258, 0.5773502691896258, 0.5773502691896258],
        [0.5773502691896258, -0.5773502691896258, -0.5773502691896258],
        [-0.5773502691896258, 0.5773502691896258, -0.5773502691896258],
        [-0.5773502691896258, -0.5773502691896258, 0.5773502691896258],
    ];
    const SQUARE_PLANAR: [[f64; 3]; 4] = [
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, -1.0, 0.0],
    ];
    const TRIGONAL_BIPYRAMIDAL: [[f64; 3]; 5] = [
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
        [1.0, 0.0, 0.0],
        [-0.5, 0.8660254037844386, 0.0],
        [-0.5, -0.8660254037844386, 0.0],
    ];
    const OCTAHEDRAL: [[f64; 3]; 6] = [
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];

    match geometry {
        CoordinationGeometryGuess::Linear => &LINEAR,
        CoordinationGeometryGuess::Trigonal => &TRIGONAL,
        CoordinationGeometryGuess::Tetrahedral => &TETRAHEDRAL,
        CoordinationGeometryGuess::SquarePlanar => &SQUARE_PLANAR,
        CoordinationGeometryGuess::TrigonalBipyramidal => &TRIGONAL_BIPYRAMIDAL,
        CoordinationGeometryGuess::Octahedral => &OCTAHEDRAL,
        CoordinationGeometryGuess::Unknown => &[],
    }
}

fn assignment_cost(vectors: &[[f64; 3]], ideals: &[[f64; 3]]) -> f64 {
    fn recurse(
        vectors: &[[f64; 3]],
        ideals: &[[f64; 3]],
        used: &mut [bool],
        idx: usize,
        current: f64,
        best: &mut f64,
    ) {
        if idx == vectors.len() {
            *best = best.min(current);
            return;
        }
        if current >= *best {
            return;
        }

        for ideal_idx in 0..ideals.len() {
            if used[ideal_idx] {
                continue;
            }
            used[ideal_idx] = true;
            let cost = 1.0 - dot(vectors[idx], ideals[ideal_idx]);
            recurse(vectors, ideals, used, idx + 1, current + cost, best);
            used[ideal_idx] = false;
        }
    }

    let mut used = vec![false; ideals.len()];
    let mut best = f64::INFINITY;
    recurse(vectors, ideals, &mut used, 0, 0.0, &mut best);
    best / (vectors.len() as f64)
}

fn classify_geometry(vectors: &[[f64; 3]]) -> (CoordinationGeometryGuess, f64) {
    let candidates: &[CoordinationGeometryGuess] = match vectors.len() {
        2 => &[CoordinationGeometryGuess::Linear],
        3 => &[CoordinationGeometryGuess::Trigonal],
        4 => &[
            CoordinationGeometryGuess::Tetrahedral,
            CoordinationGeometryGuess::SquarePlanar,
        ],
        5 => &[CoordinationGeometryGuess::TrigonalBipyramidal],
        6 => &[CoordinationGeometryGuess::Octahedral],
        _ => return (CoordinationGeometryGuess::Unknown, 0.0),
    };

    let mut best_geometry = CoordinationGeometryGuess::Unknown;
    let mut best_cost = f64::INFINITY;

    for &geometry in candidates {
        let ideals = ideal_directions(geometry);
        let cost = assignment_cost(vectors, ideals);
        if cost < best_cost {
            best_cost = cost;
            best_geometry = geometry;
        }
    }

    let score = (1.0 - (best_cost / 2.0).sqrt()).clamp(0.0, 1.0);
    (best_geometry, score)
}

/// Locate the transition-metal centres of a molecule and classify their
/// coordination geometry.
///
/// `elements` holds the atomic numbers (Z) of the atoms and `positions` their
/// Cartesian coordinates in **Ångström**; the two slices are walked in
/// parallel and the shorter one bounds the analysis, so a missing trailing
/// coordinate is treated as "no atom" rather than causing a panic.
///
/// An atom counts as a ligand of a metal when its distance to the metal is
/// at most `1.3 · (r_covalent(metal) + r_covalent(ligand)) + 0.25 Å` and it
/// is not a hydrogen. Only elements with Z in the d-block ranges 21–30,
/// 39–48 and 72–80 are treated as metal centres.
///
/// The returned geometry is the candidate polyhedron whose ideal directions
/// best match the observed metal–ligand vectors; the accompanying
/// `geometry_score` (0.0–1.0) reports how good that match is. The result is a
/// geometric guess: it does not consider oxidation state, d-electron count or
/// ligand-field stabilisation.
///
/// # Example
///
/// ```
/// use sci_form::topology::analyze_topology;
///
/// // A square-planar Pt(II) centre with four chloride ligands.
/// let elements = [78u8, 17, 17, 17, 17];
/// let positions = [
///     [0.0, 0.0, 0.0],
///     [2.32, 0.0, 0.0],
///     [-2.32, 0.0, 0.0],
///     [0.0, 2.05, 0.0],
///     [0.0, -2.05, 0.0],
/// ];
/// let result = analyze_topology(&elements, &positions);
/// assert_eq!(result.metal_centers.len(), 1);
/// assert_eq!(result.metal_centers[0].coordination_number, 4);
/// ```
pub fn analyze_topology(elements: &[u8], positions: &[[f64; 3]]) -> TopologyAnalysisResult {
    let mut metal_centers = Vec::new();

    // Be tolerant of mismatched input lengths: a missing coordinate is
    // treated as "not a ligand" rather than an out-of-bounds panic.
    let n_atoms = elements.len().min(positions.len());

    for (metal_idx, &metal_z) in elements.iter().take(n_atoms).enumerate() {
        if !crate::eht::is_transition_metal(metal_z) {
            continue;
        }

        let mut ligand_indices = Vec::new();
        let mut vectors = Vec::new();
        for ligand_idx in 0..n_atoms {
            if ligand_idx == metal_idx {
                continue;
            }
            let dist = distance(positions[metal_idx], positions[ligand_idx]);
            if is_likely_coordinated(metal_z, elements[ligand_idx], dist) {
                ligand_indices.push(ligand_idx);
                vectors.push(normalize([
                    positions[ligand_idx][0] - positions[metal_idx][0],
                    positions[ligand_idx][1] - positions[metal_idx][1],
                    positions[ligand_idx][2] - positions[metal_idx][2],
                ]));
            }
        }

        let (geometry, geometry_score) = classify_geometry(&vectors);
        metal_centers.push(MetalCoordinationCenter {
            atom_index: metal_idx,
            element: metal_z,
            coordination_number: ligand_indices.len(),
            ligand_indices,
            geometry,
            geometry_score,
        });
    }

    let warnings = if metal_centers.is_empty() {
        vec!["No transition-metal centers detected for topology analysis.".to_string()]
    } else {
        Vec::new()
    };

    TopologyAnalysisResult {
        metal_centers,
        warnings,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_square_planar_pt() {
        let elements = [78u8, 17, 17, 7, 7, 1, 1, 1, 1, 1, 1];
        let positions = [
            [0.0, 0.0, 0.0],
            [2.32, 0.0, 0.0],
            [-2.32, 0.0, 0.0],
            [0.0, 2.05, 0.0],
            [0.0, -2.05, 0.0],
            [0.90, 2.65, 0.0],
            [-0.90, 2.65, 0.0],
            [0.00, 2.05, 0.95],
            [0.90, -2.65, 0.0],
            [-0.90, -2.65, 0.0],
            [0.00, -2.05, -0.95],
        ];

        let result = analyze_topology(&elements, &positions);
        assert_eq!(result.metal_centers.len(), 1);
        assert_eq!(result.metal_centers[0].coordination_number, 4);
        assert_eq!(
            result.metal_centers[0].geometry,
            CoordinationGeometryGuess::SquarePlanar
        );
    }

    #[test]
    fn test_detect_octahedral_fe() {
        let elements = [26u8, 17, 17, 17, 17, 17, 17];
        let positions = [
            [0.0, 0.0, 0.0],
            [2.30, 0.0, 0.0],
            [-2.30, 0.0, 0.0],
            [0.0, 2.30, 0.0],
            [0.0, -2.30, 0.0],
            [0.0, 0.0, 2.30],
            [0.0, 0.0, -2.30],
        ];

        let result = analyze_topology(&elements, &positions);
        assert_eq!(
            result.metal_centers[0].geometry,
            CoordinationGeometryGuess::Octahedral
        );
    }

    #[test]
    fn test_detect_trigonal_bipyramidal_center() {
        let elements = [25u8, 17, 17, 17, 17, 17];
        let positions = [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 2.20],
            [0.0, 0.0, -2.20],
            [2.10, 0.0, 0.0],
            [-1.05, 1.8186533479473213, 0.0],
            [-1.05, -1.8186533479473213, 0.0],
        ];

        let result = analyze_topology(&elements, &positions);
        assert_eq!(
            result.metal_centers[0].geometry,
            CoordinationGeometryGuess::TrigonalBipyramidal
        );
    }
}
