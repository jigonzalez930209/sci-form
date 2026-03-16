use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CoordinationGeometryGuess {
    Linear,
    Trigonal,
    Tetrahedral,
    SquarePlanar,
    TrigonalBipyramidal,
    Octahedral,
    Unknown,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetalCoordinationCenter {
    pub atom_index: usize,
    pub element: u8,
    pub ligand_indices: Vec<usize>,
    pub coordination_number: usize,
    pub geometry: CoordinationGeometryGuess,
    pub geometry_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TopologyAnalysisResult {
    pub metal_centers: Vec<MetalCoordinationCenter>,
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

pub fn analyze_topology(elements: &[u8], positions: &[[f64; 3]]) -> TopologyAnalysisResult {
    let mut metal_centers = Vec::new();

    for (metal_idx, &metal_z) in elements.iter().enumerate() {
        if !crate::eht::is_transition_metal(metal_z) {
            continue;
        }

        let mut ligand_indices = Vec::new();
        let mut vectors = Vec::new();
        for ligand_idx in 0..elements.len() {
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
