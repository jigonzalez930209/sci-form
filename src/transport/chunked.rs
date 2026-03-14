//! Chunked streaming for large datasets to avoid JSON bottlenecks.
//!
//! Instead of serializing entire grids or coordinate arrays as one JSON blob,
//! this module provides chunk-based iteration that can be streamed over
//! WebSocket, Server-Sent Events, or consumed incrementally.

use serde::{Deserialize, Serialize};

/// A chunk of data from a streamed computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataChunk {
    /// Chunk index (0-based).
    pub index: usize,
    /// Total number of chunks (if known).
    pub total: Option<usize>,
    /// Data type identifier.
    pub kind: ChunkKind,
    /// Float values in this chunk.
    pub values: Vec<f64>,
    /// Number of logical items in this chunk.
    pub count: usize,
}

/// Type of data in a chunk.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum ChunkKind {
    /// 3D coordinates (flat xyz triples).
    Coordinates,
    /// ESP grid values.
    EspValues,
    /// DOS curve data.
    DosValues,
    /// Generic numeric data.
    Generic,
}

/// Iterator that produces chunks from a large float buffer.
pub struct ChunkedIterator {
    data: Vec<f64>,
    chunk_size: usize,
    position: usize,
    kind: ChunkKind,
    total_chunks: usize,
    current_index: usize,
}

impl ChunkedIterator {
    /// Create a new chunked iterator over data.
    ///
    /// `chunk_size`: number of f64 values per chunk.
    pub fn new(data: Vec<f64>, chunk_size: usize, kind: ChunkKind) -> Self {
        let total = (data.len() + chunk_size - 1) / chunk_size.max(1);
        Self {
            data,
            chunk_size: chunk_size.max(1),
            position: 0,
            kind,
            total_chunks: total,
            current_index: 0,
        }
    }

    /// Total number of chunks.
    pub fn total_chunks(&self) -> usize {
        self.total_chunks
    }

    /// Whether all data has been consumed.
    pub fn is_done(&self) -> bool {
        self.position >= self.data.len()
    }
}

impl Iterator for ChunkedIterator {
    type Item = DataChunk;

    fn next(&mut self) -> Option<DataChunk> {
        if self.position >= self.data.len() {
            return None;
        }

        let end = (self.position + self.chunk_size).min(self.data.len());
        let values = self.data[self.position..end].to_vec();
        let count = values.len();
        let chunk = DataChunk {
            index: self.current_index,
            total: Some(self.total_chunks),
            kind: self.kind,
            values,
            count,
        };

        self.position = end;
        self.current_index += 1;
        Some(chunk)
    }
}

/// Split an ESP grid into chunks for streaming.
pub fn chunk_esp_grid(grid: &crate::esp::EspGrid, chunk_size: usize) -> ChunkedIterator {
    ChunkedIterator::new(grid.values.clone(), chunk_size, ChunkKind::EspValues)
}

/// Split coordinate data into chunks (aligned to xyz triples).
pub fn chunk_coordinates(coords: &[f64], atoms_per_chunk: usize) -> ChunkedIterator {
    ChunkedIterator::new(coords.to_vec(), atoms_per_chunk * 3, ChunkKind::Coordinates)
}

/// Split DOS data into chunks.
pub fn chunk_dos(dos: &crate::dos::DosResult, points_per_chunk: usize) -> Vec<DataChunk> {
    let mut chunks = Vec::new();
    let n = dos.energies.len();
    let total = (n + points_per_chunk - 1) / points_per_chunk.max(1);

    for i in 0..total {
        let start = i * points_per_chunk;
        let end = (start + points_per_chunk).min(n);

        // Interleave energy and DOS values: [e0, d0, e1, d1, ...]
        let mut values = Vec::with_capacity((end - start) * 2);
        for j in start..end {
            values.push(dos.energies[j]);
            values.push(dos.total_dos[j]);
        }

        chunks.push(DataChunk {
            index: i,
            total: Some(total),
            kind: ChunkKind::DosValues,
            values,
            count: end - start,
        });
    }
    chunks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunked_iterator_basic() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut iter = ChunkedIterator::new(data, 2, ChunkKind::Generic);

        let c1 = iter.next().unwrap();
        assert_eq!(c1.index, 0);
        assert_eq!(c1.values, vec![1.0, 2.0]);
        assert_eq!(c1.total, Some(3));

        let c2 = iter.next().unwrap();
        assert_eq!(c2.index, 1);
        assert_eq!(c2.values, vec![3.0, 4.0]);

        let c3 = iter.next().unwrap();
        assert_eq!(c3.index, 2);
        assert_eq!(c3.values, vec![5.0]);

        assert!(iter.next().is_none());
        assert!(iter.is_done());
    }

    #[test]
    fn test_chunked_iterator_exact_division() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let iter = ChunkedIterator::new(data, 3, ChunkKind::Generic);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].count, 3);
        assert_eq!(chunks[1].count, 3);
    }

    #[test]
    fn test_chunk_coordinates_alignment() {
        // 4 atoms = 12 coords; chunk by 2 atoms = 6 values
        let coords = vec![0.0; 12];
        let iter = chunk_coordinates(&coords, 2);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 2);
        assert_eq!(chunks[0].count, 6);
        assert_eq!(chunks[0].kind, ChunkKind::Coordinates);
    }

    #[test]
    fn test_chunk_esp() {
        let grid = crate::esp::EspGrid {
            origin: [0.0; 3],
            spacing: 0.5,
            dims: [3, 3, 3],
            values: vec![0.1; 27],
        };
        let iter = chunk_esp_grid(&grid, 10);
        let chunks: Vec<_> = iter.collect();
        assert_eq!(chunks.len(), 3); // 27 / 10 = 3 chunks
        assert_eq!(chunks[0].count, 10);
        assert_eq!(chunks[1].count, 10);
        assert_eq!(chunks[2].count, 7);
    }

    #[test]
    fn test_chunk_dos() {
        let dos = crate::dos::DosResult {
            energies: vec![1.0, 2.0, 3.0, 4.0, 5.0],
            total_dos: vec![0.1, 0.2, 0.3, 0.4, 0.5],
            pdos: vec![],
            sigma: 0.3,
        };
        let chunks = chunk_dos(&dos, 2);
        assert_eq!(chunks.len(), 3); // 5 / 2 = 3
                                     // First chunk: [e0, d0, e1, d1]
        assert_eq!(chunks[0].values, vec![1.0, 0.1, 2.0, 0.2]);
        assert_eq!(chunks[0].count, 2);
    }

    #[test]
    fn test_total_chunks() {
        let iter = ChunkedIterator::new(vec![0.0; 100], 30, ChunkKind::Generic);
        assert_eq!(iter.total_chunks(), 4);
    }
}
