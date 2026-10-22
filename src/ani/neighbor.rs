//! Cell-list spatial partitioning for efficient neighbor search.
//!
//! Bins atoms into a 3D grid of cells with side length equal to the cutoff
//! radius. Only atoms in adjacent cells need distance evaluation, giving
//! O(N) scaling for uniformly distributed systems.

/// A neighbor pair: (atom_i, atom_j, squared_distance).
#[derive(Debug, Clone, Copy)]
pub struct NeighborPair {
    pub i: usize,
    pub j: usize,
    pub dist_sq: f64,
}

/// Spatial cell-list for fast cutoff-radius neighbor search.
pub struct CellList {
    /// Cutoff radius.
    pub cutoff: f64,
    /// Cutoff radius squared (for comparison without sqrt).
    cutoff_sq: f64,
    /// Cell side length (= cutoff).
    _cell_size: f64,
    /// Number of cells along each axis.
    n_cells: [usize; 3],
    /// Origin (minimum corner) of the bounding box.
    _origin: [f64; 3],
    /// Atom indices in each cell, indexed by flat cell index.
    cells: Vec<Vec<usize>>,
}

impl CellList {
    /// Build a cell list from atom positions.
    ///
    /// `positions` is a slice of `[x, y, z]` arrays.
    ///
    /// This constructor never panics: the number of cells per axis and their
    /// product are computed with checked arithmetic and clamped, so a
    /// pathological geometry (or a non-finite coordinate) cannot overflow the
    /// allocation.  Use [`CellList::try_new`] to get a proper error message
    /// for a non-finite coordinate.
    pub fn new(positions: &[[f64; 3]], cutoff: f64) -> Self {
        match Self::try_new(positions, cutoff) {
            Ok(list) => list,
            Err(_) => {
                // Fall back to a single cell: correctness (every pair is
                // checked) is preserved, only the speed-up is lost.
                let mut cells = vec![Vec::new(); 1];
                cells[0].extend(0..positions.len());
                CellList {
                    cutoff,
                    cutoff_sq: cutoff * cutoff,
                    _cell_size: cutoff,
                    n_cells: [1, 1, 1],
                    _origin: [0.0, 0.0, 0.0],
                    cells,
                }
            }
        }
    }

    /// Build a cell list from atom positions, reporting invalid input.
    ///
    /// Returns `Err` when a coordinate is not finite (a `NaN`/`Inf` used to
    /// blow up the cell count and abort the process with a capacity
    /// overflow) or when the required number of cells does not fit in a
    /// `usize`.
    pub fn try_new(positions: &[[f64; 3]], cutoff: f64) -> Result<Self, String> {
        if !cutoff.is_finite() || cutoff <= 0.0 {
            return Err(format!(
                "cutoff must be a finite positive number (got {})",
                cutoff
            ));
        }
        let padding = 0.01;
        let mut min = [f64::INFINITY; 3];
        let mut max = [f64::NEG_INFINITY; 3];

        for (atom, pos) in positions.iter().enumerate() {
            for d in 0..3 {
                if !pos[d].is_finite() {
                    let axis_name = ["x", "y", "z"][d];
                    return Err(format!(
                        "atom {} has a non-finite {}-coordinate ({}) (NaN and Inf are not valid Cartesian coordinates)",
                        atom, axis_name, pos[d]
                    ));
                }
                if pos[d] < min[d] {
                    min[d] = pos[d];
                }
                if pos[d] > max[d] {
                    max[d] = pos[d];
                }
            }
        }

        let origin = [min[0] - padding, min[1] - padding, min[2] - padding];
        let cell_size = cutoff;
        let mut n_cells = [1usize; 3];
        for d in 0..3 {
            let extent = (max[d] - origin[d] + padding) / cell_size;
            if !extent.is_finite() {
                return Err(format!(
                    "cell count along axis {} is not finite (the geometry does not fit in a cell list with cutoff {})",
                    d, cell_size
                ));
            }
            n_cells[d] = (extent.ceil().max(1.0) as usize).max(1);
        }

        let total = n_cells[0]
            .checked_mul(n_cells[1])
            .and_then(|v| v.checked_mul(n_cells[2]))
            .ok_or_else(|| {
                format!(
                    "cell count {}x{}x{} overflows usize (cutoff {} is too small for this geometry)",
                    n_cells[0], n_cells[1], n_cells[2], cell_size
                )
            })?;
        let mut cells = vec![Vec::new(); total];

        for (idx, pos) in positions.iter().enumerate() {
            let ci = Self::cell_index_static(pos, &origin, cell_size, &n_cells);
            cells[ci].push(idx);
        }

        Ok(CellList {
            cutoff,
            cutoff_sq: cutoff * cutoff,
            _cell_size: cell_size,
            n_cells,
            _origin: origin,
            cells,
        })
    }

    fn cell_index_static(
        pos: &[f64; 3],
        origin: &[f64; 3],
        cell_size: f64,
        n_cells: &[usize; 3],
    ) -> usize {
        let cx = ((pos[0] - origin[0]) / cell_size) as usize;
        let cy = ((pos[1] - origin[1]) / cell_size) as usize;
        let cz = ((pos[2] - origin[2]) / cell_size) as usize;
        let cx = cx.min(n_cells[0] - 1);
        let cy = cy.min(n_cells[1] - 1);
        let cz = cz.min(n_cells[2] - 1);
        cx * n_cells[1] * n_cells[2] + cy * n_cells[2] + cz
    }

    /// Find all neighbor pairs within the cutoff radius.
    ///
    /// Returns pairs with i < j to avoid duplicates.
    pub fn find_neighbors(&self, positions: &[[f64; 3]]) -> Vec<NeighborPair> {
        let mut pairs = Vec::new();
        let nc = &self.n_cells;

        for cx in 0..nc[0] {
            for cy in 0..nc[1] {
                for cz in 0..nc[2] {
                    let ci = cx * nc[1] * nc[2] + cy * nc[2] + cz;
                    self.pairs_within_cell(&self.cells[ci], positions, &mut pairs);

                    // Check 13 forward-neighbors (half of 26 neighbors)
                    for &(dx, dy, dz) in &HALF_NEIGHBOR_OFFSETS {
                        let nx = cx as isize + dx;
                        let ny = cy as isize + dy;
                        let nz = cz as isize + dz;
                        if nx < 0 || ny < 0 || nz < 0 {
                            continue;
                        }
                        let (nx, ny, nz) = (nx as usize, ny as usize, nz as usize);
                        if nx >= nc[0] || ny >= nc[1] || nz >= nc[2] {
                            continue;
                        }
                        let ni = nx * nc[1] * nc[2] + ny * nc[2] + nz;
                        self.pairs_between_cells(
                            &self.cells[ci],
                            &self.cells[ni],
                            positions,
                            &mut pairs,
                        );
                    }
                }
            }
        }
        pairs
    }

    fn pairs_within_cell(
        &self,
        cell: &[usize],
        positions: &[[f64; 3]],
        pairs: &mut Vec<NeighborPair>,
    ) {
        for a in 0..cell.len() {
            for b in (a + 1)..cell.len() {
                let i = cell[a];
                let j = cell[b];
                let dsq = dist_sq(&positions[i], &positions[j]);
                if dsq < self.cutoff_sq {
                    pairs.push(NeighborPair { i, j, dist_sq: dsq });
                }
            }
        }
    }

    fn pairs_between_cells(
        &self,
        cell_a: &[usize],
        cell_b: &[usize],
        positions: &[[f64; 3]],
        pairs: &mut Vec<NeighborPair>,
    ) {
        for &i in cell_a {
            for &j in cell_b {
                let dsq = dist_sq(&positions[i], &positions[j]);
                if dsq < self.cutoff_sq {
                    let (lo, hi) = if i < j { (i, j) } else { (j, i) };
                    pairs.push(NeighborPair {
                        i: lo,
                        j: hi,
                        dist_sq: dsq,
                    });
                }
            }
        }
    }
}

#[inline]
fn dist_sq(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}

/// 13 forward-neighbor cell offsets (avoids double-counting).
const HALF_NEIGHBOR_OFFSETS: [(isize, isize, isize); 13] = [
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
    (1, 1, 0),
    (1, -1, 0),
    (1, 0, 1),
    (1, 0, -1),
    (0, 1, 1),
    (0, 1, -1),
    (1, 1, 1),
    (1, 1, -1),
    (1, -1, 1),
    (1, -1, -1),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_water_neighbors() {
        // Water geometry: O at origin, two H at ~0.96 Å
        let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
        let cl = CellList::new(&positions, 5.2);
        let pairs = cl.find_neighbors(&positions);
        // All 3 pairs should be within 5.2 Å
        assert_eq!(pairs.len(), 3, "Water should have 3 pairs");
    }

    #[test]
    fn test_distant_atoms_excluded() {
        let positions = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];
        let cl = CellList::new(&positions, 5.2);
        let pairs = cl.find_neighbors(&positions);
        assert_eq!(pairs.len(), 0, "Atoms 10 Å apart should not be neighbors");
    }
}
