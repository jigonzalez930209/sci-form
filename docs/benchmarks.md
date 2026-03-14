# Benchmarks

Comprehensive comparison of sci-form against RDKit, the gold standard for 3D molecular conformer generation.

## Methodology

All comparisons use **heavy-atom pairwise-distance RMSD** — the root-mean-square deviation of all pairwise distances between non-hydrogen atoms. This metric is alignment-free (no superposition needed) and focuses on the chemically important scaffold geometry.

$$\text{RMSD}_{\text{pairwise}} = \sqrt{\frac{1}{|\mathcal{P}|}\sum_{(i,j) \in \mathcal{P}} (d_{ij}^{\text{sci-form}} - d_{ij}^{\text{RDKit}})^2}$$

where $\mathcal{P}$ is the set of all heavy-atom pairs.

## Diverse Molecule Benchmark

131 curated molecules spanning 27 chemical categories, from simple alkanes to macrocycles and metal-containing compounds.

### Overall Results

| Metric | Value |
|--------|-------|
| Total molecules | 131 |
| Parse success | **100%** |
| Embed success | **97.7%** (128/131) |
| Geometry quality | **97.7%** |
| Throughput | **60 mol/s** |

### Per-Category Results

<SvgDiagram src="/svg/benchmarks-rmsd.svg" alt="benchmarks-rmsd" />

### Embed Failures

Only 3 molecules fail to embed (out of 131):

| Molecule | Category | Reason |
|----------|----------|--------|
| Pyrene | polycyclic | 4-ring fused polyaromatic |
| Cubane | strained | Extreme 90° angles in 4-rings |
| Fluoranthene | polycyclic | Fused 5-6-6-5 ring system |

These are well-known hard cases for distance geometry due to their extreme geometric constraints.

## RDKit Comparison

Heavy-atom pairwise-distance RMSD between sci-form and RDKit conformers. Multi-seed ensemble comparison (5 seeds per molecule, minimum RMSD reported).

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

## GDB-20 Ensemble Comparison

Large-scale validation on 500 molecules from the GDB-20 database (molecules with up to 20 heavy atoms), using an ensemble of 5 sci-form seeds compared against 21 RDKit seeds. The minimum RMSD across all seed combinations is reported.

### Results

| Metric | All-atom | Heavy-atom |
|--------|----------|------------|
| Embed success | 100% | 100% |
| Average min-RMSD | 0.063 Å | 0.024 Å |
| > 0.1 Å | 13.00% | 9.20% |
| > 0.3 Å | 6.80% | 1.00% |
| > 0.5 Å | 1.60% | **0.00%** |
| > 0.7 Å | 0.00% | 0.00% |

### min-RMSD Distribution (all atoms)

| Range | Count | Share |
|-------|-------|-------|
| 0.00–0.05 Å | 419 | **83.80%** |
| 0.05–0.10 Å | 16 | 3.20% |
| 0.10–0.20 Å | 16 | 3.20% |
| 0.20–0.30 Å | 15 | 3.00% |
| 0.30–0.50 Å | 26 | 5.20% |
| 0.50–0.70 Å | 8 | 1.60% |
| > 0.70 Å | 0 | 0.00% |

### Ensemble Rescue Rate

Of molecules with single-seed RMSD > 0.5 Å, the multi-seed ensemble rescued **88.4%** (61/69) to below 0.5 Å. Only 8 molecules remain above the threshold after ensemble selection.

## ChemBL 10K Benchmark

Stress test on 10,000 molecules from the ChemBL database with practical pharmaceutical relevance (up to 100 atoms).

| Metric | Value |
|--------|-------|
| Parse success | **100%** |
| Embed success | **97.54%** |
| Geometry quality | **97.18%** |
| Throughput | **2.1 mol/s** |

Lower throughput is expected for larger molecules due to $O(N^3)$ scaling of Floyd-Warshall and eigendecomposition.

## Performance Scaling

<SvgDiagram src="/svg/benchmarks-success.svg" alt="benchmarks-success" />

The dominant cost is the Floyd-Warshall triangle smoothing ($O(N^3)$) and the BFGS optimization (each iteration is $O(N^2)$ for the inverse Hessian update).

## Property Calculation Performance

### Conformer Generation

| Dataset | Molecules | Success | Throughput |
|---------|-----------|---------|------------|
| Diverse (131 molecules) | 131 | 97.7% | 60 mol/s |
| ChemBL 10K | 10,000 | 97.5% | 2.1 mol/s |
| GDB-20 (500 sample) | 500 | 100% | ~50 mol/s |

Single-seed mode, no ensemble. Throughput measured on consumer hardware (8-core).

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
