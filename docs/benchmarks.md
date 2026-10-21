# Benchmarks

This page reports how `sci-form` compares against RDKit, the reference implementation for 3D molecular conformer generation.

## How these numbers were produced

There are **three different tiers** of evidence in this page. They are not interchangeable, so each section below states which tier it belongs to.

### Tier 1 — committed fixture + a test that runs in CI

The molecule corpus is committed to git and a Rust test recomputes the numbers on every CI run. These are the numbers you can trust to still hold.

| What | Fixture | Test | CI job |
|------|---------|------|--------|
| 131-molecule embed / geometry-quality sweep | `tests/fixtures/diverse_molecules.smi` | `tests/regression/test_diverse_molecules.rs::test_diverse_molecules` | `validation` → `regression` |
| GDB-20 ensemble RMSD vs RDKit | `tests/fixtures/gdb20_ensemble_300.json` + `tests/fixtures/gdb20_reference_1k.json` | `tests/regression/test_ensemble_rmsd.rs::test_ensemble_rmsd` | `validation` → `regression` |
| Bounds / embedding / 3D-FF / BFGS / SMARTS parity | `tests/fixtures/validation/*.json` | `tests/rdkit_parity.rs` | `validation` → `rdkit-parity` |
| CIP / ML descriptors / canonical SMILES / ALPB / variable cell | `tests/fixtures/validation/*.json` | `tests/rdkit_chem_parity.rs` | `validation` → `rdkit-chem-parity` |

Run them yourself:

```bash
cargo test --release --test regression test_diverse_molecules -- --nocapture
cargo test --release --test regression test_ensemble_rmsd -- --nocapture
```

### Tier 2 — Python scripts that need RDKit and are **not** re-run by CI

These scripts live in `scripts/` and need `rdkit` (and, for some of them, `pyscf`, `torch`/`torchani` or `openbabel`) installed. They also need a release build of the CLI at `target/release/sci-form`. **No CI job installs RDKit or executes these scripts**, so their output is not refreshed automatically and can drift from the current code.

| Script | What it does | Input |
|--------|--------------|-------|
| `scripts/benchmark_diverse.py` | sci-form vs RDKit on the 131-molecule corpus, heavy-atom pairwise-distance RMSD, minimum over all seed pairs | `tests/fixtures/diverse_molecules.smi` |
| `scripts/benchmark_rdkit.py` | RDKit ETKDG reference conformers for a random sample of 100 SMILES | `scripts/10k_smiles.smi` |
| `scripts/benchmark_100.py` | Three-way comparison (sci-form / RDKit / PySCF, optionally TorchANI) | `scripts/10k_smiles.smi` |
| `scripts/three_way_comparison.py` | sci-form vs RDKit vs OpenBabel on GDB-20 | GDB-20 reference fixture |
| `scripts/generate_chembl_reference.py` | RDKit reference conformers for the ChemBL 10K set | `data/chembl_10k_largest.smi` |

### Tier 3 — numbers in this page that cannot be traced to a committed artifact

Some of the figures below (notably the RDKit RMSD averages and the per-category breakdown) were recorded from a local run of a Tier 2 script. **The output of that run is not committed**, so those numbers cannot be re-derived from the repository alone, and no CI job re-runs them. They are reported here as historical measurements, with the script that produces them named so they can be regenerated on demand.

---

## Methodology

All RDKit comparisons use **heavy-atom pairwise-distance RMSD** — the root-mean-square deviation of all pairwise distances between non-hydrogen atoms. This metric is alignment-free (no superposition needed) and focuses on the chemically important scaffold geometry.

$$\text{RMSD}_{\text{pairwise}} = \sqrt{\frac{1}{|\mathcal{P}|}\sum_{(i,j) \in \mathcal{P}} (d_{ij}^{\text{sci-form}} - d_{ij}^{\text{RDKit}})^2)}$$

where $\mathcal{P}$ is the set of all heavy-atom pairs.

---

## Diverse Molecule Benchmark

**Tier 1.** 131 curated molecules spanning 27 chemical categories, from simple alkanes to macrocycles and metal-containing compounds. The corpus is `tests/fixtures/diverse_molecules.smi`; the numbers below are what `test_diverse_molecules` prints on the current code (single seed, seed 42).

### Overall Results

| Metric | Value |
|--------|-------|
| Total molecules | 131 |
| Categories | 27 |
| Parse success | **100%** (131/131) |
| Embed success | **99.2%** (130/131) |
| Geometry quality | **99.2%** (130/131) |

Throughput is machine-dependent and is printed by the test itself (138.9 mol/s on the machine that produced the table above; the test previously reported ≈60 mol/s on an 8-core consumer CPU).

### Per-Category Results

<SvgDiagram src="/svg/benchmarks-rmsd.svg" alt="benchmarks-rmsd" />

### Embed Failures

Only **one** molecule currently fails to embed (out of 131):

| Molecule | Category | Reason |
|----------|----------|--------|
| Cubane (`C12C3C4C1C5C3C2C45`) | strained | Extreme 90° angles in 4-rings — "Failed to generate valid conformer after 160 iterations" |

Every other category embeds at 100%.

> **Historical note:** earlier revisions of this page reported 97.7% embed success with three failures (pyrene, cubane, fluoranthene). Two of those three now embed successfully; the numbers above are the current output of the Tier 1 test.

---

## RDKit Comparison

**Tier 3.** Heavy-atom pairwise-distance RMSD between sci-form and RDKit conformers, multi-seed ensemble comparison (5 seeds per molecule: 42, 0, 1, 2, 3), minimum RMSD over all 5 × 5 seed pairs.

Regenerate with:

```bash
pip install rdkit numpy
cargo build --release -p sci-form-cli
python scripts/benchmark_diverse.py tests/fixtures/diverse_molecules.smi 5
```

The recorded run reported:

### Overall Results

| Metric | Value |
|--------|-------|
| Average RMSD | **0.064 Å** |
| Median RMSD | **0.011 Å** |
| < 0.1 Å | **82.8%** |
| < 0.3 Å | **94.4%** |
| < 0.5 Å | **98.4%** |
| < 1.0 Å | **100%** |

### RMSD Distribution

<SvgDiagram src="/svg/benchmarks-timing.svg" alt="benchmarks-timing" />

### Hardest Categories

| Category | Avg RMSD | Description |
|----------|---------|-------------|
| silicon | 0.543 Å | Si atom typing differences |
| selenium | 0.507 Å | Se parameter approximations |
| strained | 0.182 Å | Cubane, cyclopropane |
| polycyclic | 0.112 Å | Fused aromatic systems |

These figures are **not** produced by any test in the repository and are not regenerated by CI. Treat them as a snapshot of one local run.

---

## GDB-20 Ensemble Comparison

**Tier 1 for the corpus, Tier 3 for the statistics.** Large-scale validation against RDKit on the GDB-20 database (molecules with up to 20 heavy atoms).

The committed corpus is `tests/fixtures/gdb20_ensemble_300.json`: **300 molecules**, each with **21 RDKit conformers** (seeds 42 and 0–19). `tests/fixtures/gdb20_reference_1k.json` holds the single-seed (seed 42) RDKit reference for 1000 GDB-20 molecules and is used for the torsion data.

The test is:

```bash
cargo test --release --test regression test_ensemble_rmsd -- --nocapture
```

It prints, for the ensemble comparison, the minimum RMSD across all RDKit seeds, both all-atom and heavy-atom, together with the fraction above 0.1 / 0.3 / 0.5 Å.

By default the test uses **one** sci-form seed (42). Increase it with environment variables:

| Variable | Default | Meaning |
|----------|---------|---------|
| `OUR_SEEDS` | `1` | sci-form seeds to try (`21` → seeds 42, 0–19, matching the RDKit ensemble) |
| `GDB20_LIMIT` | `300` | Number of molecules to process (`0` → all) |
| `OWN_TORSIONS` | `0` | `1` uses sci-form's own SMARTS matcher instead of the RDKit oracle torsions |

> **Historical note:** earlier revisions of this page reported a 500-molecule run with a 5-seed sci-form ensemble and a min-RMSD distribution (83.80% below 0.05 Å) plus an "88.4% ensemble rescue rate". Those figures correspond to a local run with `OUR_SEEDS=5` over a 500-molecule sample. The committed fixture holds 300 molecules, and the numbers are **not** pinned by an assertion in the test — the test prints them rather than asserting on them. Re-run the test with the environment variables above to reproduce the current values.

---

## ChemBL 10K Benchmark

**Tier 2.** Stress test on 10,000 molecules from the ChemBL database with practical pharmaceutical relevance (up to 100 atoms). The corpus is committed as `data/chembl_10k_practical.smi` (and `data/chembl_10k_practical_asc.smi`). The RDKit reference conformers are produced by `scripts/generate_chembl_reference.py`, which writes `tests/fixtures/chembl_10k_reference.json` — that file is **not committed**, so the sweep has to be regenerated locally.

A recorded run reported:

| Metric | Value |
|--------|-------|
| Parse success | **100%** |
| Embed success | **97.54%** |
| Geometry quality | **97.18%** |
| Throughput | **2.1 mol/s** |

Lower throughput is expected for larger molecules due to $O(N^3)$ scaling of Floyd-Warshall and eigendecomposition.

No CI job re-runs this sweep, and no committed artifact stores its result.

---

## Performance Scaling

<SvgDiagram src="/svg/benchmarks-success.svg" alt="benchmarks-success" />

The dominant cost is the Floyd-Warshall triangle smoothing ($O(N^3)$) and the BFGS optimization (each iteration is $O(N^2)$ for the inverse Hessian update).

---

## Property Calculation Performance

### Conformer Generation

| Dataset | Molecules | Success | Throughput |
|---------|-----------|---------|------------|
| Diverse (131 molecules) | 131 | 99.2% | machine-dependent; printed by the test |
| ChemBL 10K | 10,000 | 97.5% | 2.1 mol/s (recorded run) |
| GDB-20 (300 sample) | 300 | 100% | machine-dependent |

Single-seed mode, no ensemble.

### Electronic Structure (EHT)

| Molecule | Atoms | Basis Functions | Time |
|----------|-------|----------------|------|
| H₂O | 3 | 5 | < 1 ms |
| Benzene | 12 | 18 | < 2 ms |
| Naphthalene | 18 | 28 | < 5 ms |
| Drug-like (~30 heavy) | ~40 | ~60 | ~10 ms |

EHT cost scales as $O(N_\text{AO}^3)$ for diagonalization. The STO-3G minimal basis keeps $N_\text{AO}$ small even for medium molecules.

### ESP Grid

| Grid Resolution | Spacing | Typical Size | Time |
|----------------|---------|-------------|------|
| Coarse | 1.0 Å | 10³ grid | < 5 ms |
| Standard | 0.5 Å | 20³ grid | ~20 ms |
| Fine | 0.2 Å | 50³ grid | ~300 ms |

ESP evaluation is $O(N_\text{atoms} \times N_\text{grid})$. The parallel evaluator (`compute_esp_grid_parallel`) gives near-linear speedup on multi-core systems.

### Complete Property Pipeline (single molecule)

Full pipeline (embed + charges + EHT + ESP + DOS + SASA + dipole) on a drug-like molecule (~30 heavy atoms):

| Step | Time |
|------|------|
| Conformer generation | ~10 ms |
| Gasteiger charges | < 1 ms |
| EHT calculation | ~5 ms |
| ESP grid (0.5 Å) | ~20 ms |
| DOS computation | < 1 ms |
| SASA | < 1 ms |
| Dipole | < 1 ms |
| **Total** | **~38 ms** |

---

## Experimental Quantum Engine Benchmarks

The isolated `experimental_2` stack adds a Roothaan-Hall RHF/STO-3G implementation with explicit two-electron integrals and rayon-parallel ERI construction.

### Regression Coverage

| Suite | Scope | Result |
|-------|-------|--------|
| `test_experimental_comparison` | SCF correctness, legacy-vs-experimental comparison, NIST references, spectroscopy smoke checks | **54 passed, 0 failed** |
| `test_extended_molecules` | Complex molecules, NMR/IR/UV-Vis reference checks, ML properties, drug-like molecules | **14 passed, 0 failed, 7 ignored** |

Both suites live in `tests/regression/` and run in the CI `regression` gate.

### Parallel ERI Speedups

Observed on an Intel i5-10500H (12 logical cores):

| System size | Typical behavior |
|------------|------------------|
| Small systems (< 15 basis functions) | ~0.8–1.1×, parallel overhead dominates |
| Benzene (~36 basis functions) | ~3–5× speedup |
| Pyridine (~42 basis functions) | ~4–6× speedup |

The heavy speedup tables are intentionally marked `#[ignore]` in debug-mode regression runs because explicit ERIs scale as $O(N^4)$ and distort CI time on larger aromatic systems. Run them in release mode:

```bash
cargo test --release --test regression test_extended_molecules -- --include-ignored
```

### Interpretation

- The current acceleration path is CPU-only via rayon
- GPU execution is planned, but the current `phase1_gpu_infrastructure` backend is still a CPU fallback with WGSL-ready interfaces
- The experimental engine is best viewed today as a validated research track rather than the default production quantum backend
