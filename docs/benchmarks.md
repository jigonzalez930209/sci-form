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

<img src="/svg/benchmarks-rmsd.svg" alt="benchmarks-rmsd" class="svg-diagram" />

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

<img src="/svg/benchmarks-timing.svg" alt="benchmarks-timing" class="svg-diagram" />

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

<img src="/svg/benchmarks-success.svg" alt="benchmarks-success" class="svg-diagram" />

The dominant cost is the Floyd-Warshall triangle smoothing ($O(N^3)$) and the BFGS optimization (each iteration is $O(N^2)$ for the inverse Hessian update).
