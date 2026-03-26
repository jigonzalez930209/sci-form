//! Web Worker dispatch strategy for WASM-safe parallelism.
//!
//! In browsers, heavy computations block the main thread. Web Workers provide
//! separate threads, but SharedArrayBuffer requires cross-origin isolation.
//!
//! This module provides:
//! - Task descriptors for dispatching work to workers
//! - Result aggregation for collecting worker outputs
//! - Batch splitting for dividing work across N workers
//!
//! The actual Web Worker creation happens in JavaScript; this module provides
//! the data structures and splitting logic.

use serde::{Deserialize, Serialize};

/// A task to be dispatched to a Web Worker.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkerTask {
    /// Unique task identifier.
    pub id: usize,
    /// Task type.
    pub kind: TaskKind,
    /// SMILES strings to process (for batch embedding tasks).
    pub smiles: Vec<String>,
    /// Numeric parameters (e.g., seed, spacing, etc.).
    pub params: Vec<f64>,
}

/// Types of tasks that can be dispatched to workers.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum TaskKind {
    /// Batch conformer generation.
    EmbedBatch,
    /// ESP grid computation.
    ComputeEsp,
    /// DOS computation.
    ComputeDos,
    /// Population analysis batch.
    ComputePopulation,
    /// UFF energy evaluation batch.
    ComputeUff,
}

/// Result from a completed worker task.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorkerResult {
    /// Task ID this result corresponds to.
    pub task_id: usize,
    /// Whether the task completed successfully.
    pub success: bool,
    /// Error message if failed.
    pub error: Option<String>,
    /// JSON-encoded result data.
    pub data: String,
}

/// Split a batch of SMILES into N worker tasks.
///
/// `smiles`: all SMILES to process
/// `n_workers`: number of workers to distribute across
/// `seed`: RNG seed
///
/// Each chunk is capped at 10 000 SMILES to prevent unbounded memory
/// allocation when few workers are requested for very large inputs.
///
/// Returns a vector of tasks, one per worker.
pub fn split_batch(smiles: &[String], n_workers: usize, seed: u64) -> Vec<WorkerTask> {
    let n = smiles.len();
    const MAX_CHUNK: usize = 10_000;
    // Ensure enough workers so no single chunk exceeds MAX_CHUNK
    let min_workers_for_cap = n.div_ceil(MAX_CHUNK);
    let workers = n_workers.max(1).max(min_workers_for_cap).min(n);
    let chunk_size = n.div_ceil(workers);

    smiles
        .chunks(chunk_size)
        .enumerate()
        .map(|(i, chunk)| WorkerTask {
            id: i,
            kind: TaskKind::EmbedBatch,
            smiles: chunk.to_vec(),
            params: vec![seed as f64],
        })
        .collect()
}

/// Merge worker results back into order (by task_id).
pub fn merge_results(mut results: Vec<WorkerResult>) -> Vec<WorkerResult> {
    results.sort_by_key(|r| r.task_id);
    results
}

/// Estimate the optimal number of workers based on data size.
///
/// Heuristic: 1 worker per 100 molecules, max 8 for browser environments.
pub fn estimate_workers(n_items: usize, max_workers: usize) -> usize {
    let ideal = (n_items / 100).max(1);
    ideal.min(max_workers).min(8)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_batch_even() {
        let smiles: Vec<String> = (0..10).map(|i| format!("mol{}", i)).collect();
        let tasks = split_batch(&smiles, 2, 42);
        assert_eq!(tasks.len(), 2);
        assert_eq!(tasks[0].smiles.len(), 5);
        assert_eq!(tasks[1].smiles.len(), 5);
        assert_eq!(tasks[0].id, 0);
        assert_eq!(tasks[1].id, 1);
    }

    #[test]
    fn test_split_batch_uneven() {
        let smiles: Vec<String> = (0..7).map(|i| format!("mol{}", i)).collect();
        let tasks = split_batch(&smiles, 3, 42);
        assert_eq!(tasks.len(), 3);
        assert_eq!(tasks[0].smiles.len(), 3);
        assert_eq!(tasks[1].smiles.len(), 3);
        assert_eq!(tasks[2].smiles.len(), 1);
    }

    #[test]
    fn test_split_batch_more_workers_than_items() {
        let smiles: Vec<String> = vec!["C".to_string(), "CC".to_string()];
        let tasks = split_batch(&smiles, 10, 42);
        assert_eq!(tasks.len(), 2); // min(10, 2) = 2 workers
    }

    #[test]
    fn test_merge_results_ordered() {
        let results = vec![
            WorkerResult {
                task_id: 2,
                success: true,
                error: None,
                data: "r2".to_string(),
            },
            WorkerResult {
                task_id: 0,
                success: true,
                error: None,
                data: "r0".to_string(),
            },
            WorkerResult {
                task_id: 1,
                success: true,
                error: None,
                data: "r1".to_string(),
            },
        ];
        let merged = merge_results(results);
        assert_eq!(merged[0].task_id, 0);
        assert_eq!(merged[1].task_id, 1);
        assert_eq!(merged[2].task_id, 2);
    }

    #[test]
    fn test_estimate_workers() {
        assert_eq!(estimate_workers(50, 4), 1);
        assert_eq!(estimate_workers(500, 8), 5);
        assert_eq!(estimate_workers(10000, 8), 8);
        assert_eq!(estimate_workers(10000, 4), 4);
    }

    #[test]
    fn test_worker_task_serialization() {
        let task = WorkerTask {
            id: 0,
            kind: TaskKind::EmbedBatch,
            smiles: vec!["C".to_string(), "CC".to_string()],
            params: vec![42.0],
        };
        let json = serde_json::to_string(&task).unwrap();
        let back: WorkerTask = serde_json::from_str(&json).unwrap();
        assert_eq!(back.id, 0);
        assert_eq!(back.smiles.len(), 2);
    }
}
