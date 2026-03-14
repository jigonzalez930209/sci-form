# Web Data Transport and Parallelism

sci-form includes a transport layer for efficient data transfer between WASM modules and the browser, avoiding JSON serialization overhead for large datasets.

## Arrow Columnar Layout

<img src="/svg/transport-arrow.svg" alt="Arrow columnar memory layout" class="svg-diagram" />

### Why Columnar?

| Approach | Serialization | Memory Layout | GPU Transfer |
|----------|---------------|---------------|-------------|
| JSON rows | Parse each field | Scattered | Requires copy+repack |
| Arrow columns | Zero-copy possible | Contiguous typed arrays | Direct `Float32Array` handoff |

### RecordBatch Structure

A `RecordBatch` organizes data into typed columns:

```rust
pub struct RecordBatch {
    columns: Vec<Column>,
    schema: Vec<ColumnSchema>,
    num_rows: usize,
}

pub enum Column {
    Float64(Vec<f64>),
    Int32(Vec<i32>),
    Uint8(Vec<u8>),
}
```

### Pack Functions

**Conformer data**:
```rust
let batch = pack_conformers(&conformer_results);
// Columns: x (f64), y (f64), z (f64), energy (f64), element (u8)
```

**ESP grid**:
```rust
let batch = pack_esp_grid(&esp_grid);
// Columns: esp_values (f64), grid_x (f64), grid_y (f64), grid_z (f64)
```

**DOS data**:
```rust
let batch = pack_dos(&dos_result);
// Columns: energies (f64), total_dos (f64), pdos per atom (f64 each)
```

## Chunked Streaming

<img src="/svg/transport-chunked.svg" alt="Chunked streaming and Web Workers" class="svg-diagram" />

### ChunkedIterator

For large datasets that shouldn't be sent as a single allocation, `ChunkedIterator` lazily yields fixed-size slices:

```rust
pub struct ChunkedIterator {
    data: Vec<f64>,
    chunk_size: usize,
    position: usize,
    total_chunks: usize,
    kind: ChunkKind,
}

impl Iterator for ChunkedIterator {
    type Item = DataChunk;
    // yields DataChunk { index, total, kind, values, count }
}
```

### Chunk Functions

**ESP grid chunking**:
$$
N_{\text{chunks}} = \lceil N_{\text{values}} / \text{chunk\_size} \rceil
$$

```rust
let chunks = chunk_esp_grid(&esp_grid, 10000);
for chunk in chunks {
    // Send chunk.values to GPU or Web Worker
}
```

**Coordinate chunking** (3 values per atom):
```rust
let chunks = chunk_coordinates(&positions, 1000);
// Each chunk has up to 1000 atoms (3000 f64 values)
```

## Web Worker Split/Merge

### Task Distribution

`split_batch` distributes SMILES across workers using round-robin:

```rust
let tasks = split_batch(&smiles_list, n_workers, seed);
// Each WorkerTask has: id, kind, smiles subset, params
```

### Worker Estimation

The optimal number of workers is:

$$
n_{\text{workers}} = \min\left(\text{max\_workers}, \; \left\lceil \frac{n_{\text{items}}}{64} \right\rceil\right)
$$

At least 1 worker, capped at `max_workers` (typically `navigator.hardwareConcurrency`).

### Result Merging

After workers complete, merge results ordered by task ID:

```rust
let merged = merge_results(worker_results);
// Sorted by task_id, preserves original ordering
```

## API

### WASM

```typescript
import { pack_batch_arrow, split_worker_tasks, estimate_workers } from 'sci-form';

// Pack conformer results into columnar format
const batch = pack_batch_arrow(smilesArray);

// Distribute work across Web Workers
const nWorkers = estimate_workers(smiles.length, navigator.hardwareConcurrency);
const tasks = split_worker_tasks(smiles, nWorkers, 42);

// Each worker processes its task subset
for (const task of tasks) {
    worker.postMessage(task);
}
```

### Python

```python
import sci_form

# Pack conformers
batch = sci_form.pack_conformers(["CCO", "c1ccccc1"])
print(batch.num_columns())  # 5
print(batch.num_rows())     # total atoms across all molecules

# Estimate parallelism
n = sci_form.estimate_workers(10000, 8)  # → 8

# Split work
tasks = sci_form.split_worker_tasks(smiles_list, n, 42)
```

## Performance Considerations

| Data Size | Recommended Approach |
|-----------|---------------------|
| < 1 MB | Direct JSON or single RecordBatch |
| 1–50 MB | Chunked streaming (avoid single allocation) |
| > 50 MB | Web Workers + chunked streaming |

The transport layer eliminates the JSON bottleneck for datasets with thousands of atoms or large volumetric grids, enabling real-time browser visualization.
