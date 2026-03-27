//! Verlet neighbor list for O(N) non-bonded force evaluation with PBC.

use super::pbc::SimulationBox;

/// Verlet neighbor list entry: pairs of atom indices within cutoff + skin.
#[derive(Debug, Clone)]
pub struct NeighborList {
    /// Neighbor pairs (i, j) with i < j.
    pub pairs: Vec<(usize, usize)>,
    /// Cutoff distance used (Å).
    pub cutoff: f64,
    /// Skin distance (Å).
    pub skin: f64,
    /// Positions when list was last built (flat).
    positions_at_build: Vec<f64>,
    /// Maximum displacement since last build.
    max_displacement: f64,
}

impl NeighborList {
    /// Build a neighbor list for all pairs within cutoff + skin.
    pub fn build(positions_flat: &[f64], sim_box: &SimulationBox, cutoff: f64, skin: f64) -> Self {
        let n = positions_flat.len() / 3;
        let r_list = cutoff + skin;
        let r_list_sq = r_list * r_list;
        let mut pairs = Vec::new();

        for i in 0..n {
            let ri = [
                positions_flat[3 * i],
                positions_flat[3 * i + 1],
                positions_flat[3 * i + 2],
            ];
            for j in (i + 1)..n {
                let rj = [
                    positions_flat[3 * j],
                    positions_flat[3 * j + 1],
                    positions_flat[3 * j + 2],
                ];
                let d2 = sim_box.distance_sq(&ri, &rj);
                if d2 < r_list_sq {
                    pairs.push((i, j));
                }
            }
        }

        Self {
            pairs,
            cutoff,
            skin,
            positions_at_build: positions_flat.to_vec(),
            max_displacement: 0.0,
        }
    }

    /// Check whether the neighbor list needs rebuilding based on atomic displacements.
    ///
    /// Rebuild is needed when twice the maximum displacement exceeds the skin distance.
    pub fn needs_rebuild(&mut self, positions_flat: &[f64]) -> bool {
        let n = positions_flat.len() / 3;
        let mut max_disp_sq = 0.0f64;
        for i in 0..n {
            let dx = positions_flat[3 * i] - self.positions_at_build[3 * i];
            let dy = positions_flat[3 * i + 1] - self.positions_at_build[3 * i + 1];
            let dz = positions_flat[3 * i + 2] - self.positions_at_build[3 * i + 2];
            max_disp_sq = max_disp_sq.max(dx * dx + dy * dy + dz * dz);
        }
        self.max_displacement = max_disp_sq.sqrt();
        // Rebuild needed when any atom has moved more than skin/2
        2.0 * self.max_displacement > self.skin
    }

    /// Rebuild the list at new positions.
    pub fn rebuild(&mut self, positions_flat: &[f64], sim_box: &SimulationBox) {
        *self = Self::build(positions_flat, sim_box, self.cutoff, self.skin);
    }

    /// Number of neighbor pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Whether the list is empty.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dynamics_live::pbc::SimulationBox;

    #[test]
    fn build_includes_close_pairs() {
        // Two atoms within cutoff
        let positions = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let sim_box = SimulationBox::none();
        let nl = NeighborList::build(&positions, &sim_box, 3.0, 0.5);
        assert_eq!(nl.pairs.len(), 1);
        assert_eq!(nl.pairs[0], (0, 1));
    }

    #[test]
    fn build_excludes_distant_pairs() {
        // Two atoms beyond cutoff+skin
        let positions = [0.0, 0.0, 0.0, 10.0, 0.0, 0.0];
        let sim_box = SimulationBox::none();
        let nl = NeighborList::build(&positions, &sim_box, 3.0, 0.5);
        assert!(nl.is_empty());
    }

    #[test]
    fn needs_rebuild_after_displacement() {
        let positions = [0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let sim_box = SimulationBox::none();
        let mut nl = NeighborList::build(&positions, &sim_box, 3.0, 0.5);
        // Small displacement — no rebuild
        let pos2 = [0.1, 0.0, 0.0, 1.5, 0.0, 0.0];
        assert!(!nl.needs_rebuild(&pos2));
        // Large displacement — triggers rebuild
        let pos3 = [0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        assert!(nl.needs_rebuild(&pos3));
    }

    #[test]
    fn len_and_is_empty() {
        let positions = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 10.0, 0.0, 0.0];
        let sim_box = SimulationBox::none();
        let nl = NeighborList::build(&positions, &sim_box, 3.0, 0.5);
        assert_eq!(nl.len(), 1); // only (0,1)
        assert!(!nl.is_empty());
    }
}
