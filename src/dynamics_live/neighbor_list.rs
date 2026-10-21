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
    ///
    /// The pairs are found with a cell list whose cell edge is `cutoff + skin`,
    /// so the search is O(N) in the number of atoms for a homogeneous system.
    /// The 27 neighbouring cells are scanned for every atom and the
    /// minimum-image distance decides the membership of a pair.
    ///
    /// The cell-list shortcut is only valid when the shortest box vector is at
    /// least `2 * (cutoff + skin)` long; otherwise a pair can be closer than
    /// `cutoff + skin` across a cell boundary that is not in the 27-cell
    /// neighbourhood, and the exact O(N^2) scan is used instead.
    pub fn build(positions_flat: &[f64], sim_box: &SimulationBox, cutoff: f64, skin: f64) -> Self {
        let n = positions_flat.len() / 3;
        let r_list = cutoff + skin;
        let r_list_sq = r_list * r_list;
        let mut pairs = Vec::new();

        let use_cells = sim_box.periodic && sim_box.shortest_vector() >= 2.0 * r_list;

        if use_cells {
            // Bin the atoms into cubic cells of edge r_list.
            let mut cell_min = [f64::INFINITY; 3];
            let mut cell_max = [f64::NEG_INFINITY; 3];
            for i in 0..n {
                for d in 0..3 {
                    let v = positions_flat[3 * i + d];
                    if v < cell_min[d] {
                        cell_min[d] = v;
                    }
                    if v > cell_max[d] {
                        cell_max[d] = v;
                    }
                }
            }
            let n_cells = [
                (((cell_max[0] - cell_min[0]) / r_list).ceil() as usize).max(1),
                (((cell_max[1] - cell_min[1]) / r_list).ceil() as usize).max(1),
                (((cell_max[2] - cell_min[2]) / r_list).ceil() as usize).max(1),
            ];
            let total_cells = n_cells[0] * n_cells[1] * n_cells[2];
            let mut bins: Vec<Vec<usize>> = vec![Vec::new(); total_cells];

            let cell_index = |pos: &[f64; 3]| -> usize {
                let cx = (((pos[0] - cell_min[0]) / r_list).floor().max(0.0) as usize)
                    .min(n_cells[0] - 1);
                let cy = (((pos[1] - cell_min[1]) / r_list).floor().max(0.0) as usize)
                    .min(n_cells[1] - 1);
                let cz = (((pos[2] - cell_min[2]) / r_list).floor().max(0.0) as usize)
                    .min(n_cells[2] - 1);
                cx * n_cells[1] * n_cells[2] + cy * n_cells[2] + cz
            };

            let mut atom_cell = vec![0usize; n];
            for i in 0..n {
                let pos = [
                    positions_flat[3 * i],
                    positions_flat[3 * i + 1],
                    positions_flat[3 * i + 2],
                ];
                let ci = cell_index(&pos);
                atom_cell[i] = ci;
                bins[ci].push(i);
            }

            for i in 0..n {
                let ri = [
                    positions_flat[3 * i],
                    positions_flat[3 * i + 1],
                    positions_flat[3 * i + 2],
                ];
                let ci = atom_cell[i];
                let cx = ci / (n_cells[1] * n_cells[2]);
                let cy = (ci / n_cells[2]) % n_cells[1];
                let cz = ci % n_cells[2];
                for dx in -1isize..=1 {
                    for dy in -1isize..=1 {
                        for dz in -1isize..=1 {
                            let nx = cx as isize + dx;
                            let ny = cy as isize + dy;
                            let nz = cz as isize + dz;
                            if nx < 0
                                || ny < 0
                                || nz < 0
                                || nx >= n_cells[0] as isize
                                || ny >= n_cells[1] as isize
                                || nz >= n_cells[2] as isize
                            {
                                continue;
                            }
                            let ni = (nx as usize) * n_cells[1] * n_cells[2]
                                + (ny as usize) * n_cells[2]
                                + nz as usize;
                            for &j in &bins[ni] {
                                if j <= i {
                                    continue;
                                }
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
                    }
                }
            }
            pairs.sort_unstable();
        } else {
            pairs = Self::brute_force_pairs(positions_flat, sim_box, r_list_sq);
        }

        Self {
            pairs,
            cutoff,
            skin,
            positions_at_build: positions_flat.to_vec(),
            max_displacement: 0.0,
        }
    }

    /// Exact O(N^2) scan of all pairs within `r_list_sq`.
    ///
    /// This is the reference implementation the cell list is checked against.
    pub fn brute_force_pairs(
        positions_flat: &[f64],
        sim_box: &SimulationBox,
        r_list_sq: f64,
    ) -> Vec<(usize, usize)> {
        let n = positions_flat.len() / 3;
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
        pairs
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

    /// The cell-list path and the exact O(N^2) scan must produce the same
    /// pair set for the same simulation box.
    #[test]
    fn cell_list_matches_brute_force() {
        // 12 atoms on a simple cubic lattice.
        let mut positions = Vec::new();
        for i in 0..12 {
            let x = (i % 4) as f64 * 2.0;
            let y = ((i / 4) % 3) as f64 * 2.0;
            let z = (i / 12) as f64 * 2.0;
            positions.extend_from_slice(&[x, y, z]);
        }

        // A box whose shortest vector is at least 2*(cutoff+skin) long takes
        // the cell-list path.
        let sim = SimulationBox::orthorhombic(20.0, 20.0, 20.0);
        assert!(
            sim.shortest_vector() >= 2.0 * 3.5,
            "the test box must exercise the cell-list path"
        );
        let cells = NeighborList::build(&positions, &sim, 3.0, 0.5);

        // The reference is the exact scan on the very same box.
        let r_list_sq = 3.5 * 3.5;
        let brute = NeighborList::brute_force_pairs(&positions, &sim, r_list_sq);

        assert!(
            !brute.is_empty(),
            "the reference configuration should have neighbours"
        );
        assert_eq!(
            brute, cells.pairs,
            "the cell list and the exact scan must agree"
        );
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
