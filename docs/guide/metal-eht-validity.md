# Metal EHT Validity and Confidence

This page defines the current validity domain for transition-metal calculations in sci-form.

## Scope

Transition-metal EHT support is available for:

- First-row metals: Sc-Zn
- Second-row metals: Y-Cd
- Third-row subset: Hf-Hg

These systems use an extended valence basis with real d orbitals and are reported with experimental confidence metadata.

## Confidence Levels

- High: calibrated and benchmarked workflows (current default for stable organic main-group EHT workflows)
- Experimental: operational workflow with regression coverage but without full literature calibration (current level for transition-metal EHT)
- Unsupported: no reliable EHT parameterization for requested elements or workflow

For transition-metal systems, treat absolute orbital energies and fine MO ordering as provisional until full calibration is complete.

## Quantified Regression Baseline

The project includes fixture-driven regression tests for representative metal systems:

- Ferrocene-like
- Cisplatin-like square-planar Pt complex
- PdCl4-like square-planar Pd complex
- FeCl6 octahedral complex

Validation checks currently enforce:

- Orbital count consistency
- Valence-electron count consistency
- Support-level metadata consistency
- HOMO, LUMO, and gap stability against fixture baselines
- Sub-0.1% variation threshold for HOMO/LUMO (and for gap when reference gap is not near zero)

Reference fixture and test implementation:

- Fixture: `tests/fixtures/eht_metal_reference.json`
- Test: `tests/test_eht_metal_references.rs`

Run locally:

```bash
cargo test --release --test test_eht_metal_references
```

## Experimental Geometry Benchmark (<1% Error)

To validate representative transition-metal systems against experimental structural data, sci-form includes a geometry benchmark with light and heavy complexes:

- Light: ferrocene (Fe)
- Heavy: cisplatin (Pt), hexachloroplatinate motif (Pt)

The benchmark fixture stores target metal-ligand distances and a strict maximum relative error of 1%:

- Fixture: `tests/fixtures/metal_experimental_geometry.json`
- Test: `tests/test_metal_experimental_geometry.rs`

Run locally:

```bash
cargo test --release --test test_metal_experimental_geometry
```

The test computes average metal-ligand distances from platform geometry outputs and enforces:

$$
	ext{error\_\%} = \frac{|d_{pred} - d_{exp}|}{d_{exp}} \times 100 \le 1.0
$$

## Current Practical Guidance

- Use EHT metal results for qualitative analysis first:
  - frontier-orbital shape inspection
  - relative trends across related structures
  - exploratory DOS and population analysis
- Prefer UFF-backed geometry and energy workflows when the target task is force-field centric.
- For publication-grade quantitative metal energetics, cross-check against calibrated external methods.

## Reactivity Descriptor Validity

New reactivity workflows are available through:

- `compute_fukui_descriptors(elements, positions)`
- `compute_reactivity_ranking(elements, positions)`

These workflows use EHT-derived frontier atom contributions as low-cost proxies:

- $f^+$ uses LUMO atom contributions
- $f^-$ uses HOMO atom contributions
- $f^0 = \frac{f^+ + f^-}{2}$
- dual descriptor uses $f^+ - f^-$

Interpretation scope:

- Organic main-group systems: useful qualitative site-priority trends in related series.
- Transition-metal systems: exploratory only; use as ranking hints and cross-check with calibrated external methods.
- Unsupported EHT elements: descriptor reliability is not guaranteed.

Empirical ranking output mixes condensed Fukui terms with Mulliken-charge bias. It is intended for triage and visualization, not as a substitute for kinetics or high-level reactivity models.

## Exploratory UV-Vis-Like Output

An exploratory spectral helper is available via:

- `compute_uv_vis_spectrum(elements, positions, sigma, e_min, e_max, n_points)`

This output is generated from occupied→virtual EHT MO energy differences with a coefficient-overlap intensity proxy and Gaussian broadening.

Limits:

- Not a calibrated excited-state method.
- No CI, TDDFT, or solvent/environment effects.
- Best used as an interactive qualitative view for trend comparison.

## Empirical pKa Heuristic

An empirical pKa helper is available via:

- `compute_empirical_pka(smiles)`

This workflow uses graph-environment rules (for example carboxylic-acid oxygen, phenol oxygen, aliphatic amine, aromatic nitrogen) combined with Gasteiger-charge adjustments to estimate coarse acidic/basic site trends.

Limits:

- Intended for ranking and triage, not calibrated thermodynamic pKa prediction.
- Works best on common organic functional groups.
- Treat transition-metal and unusual ionic systems as low confidence.

## Aromatic UFF Heuristic

An aromaticity-informed UFF helper is available via:

- `compute_uff_energy_with_aromatic_heuristics(smiles, coords)`

It reports raw UFF energy plus a lightweight aromatic stabilization correction based on aromatic bond count. This is a downstream ranking heuristic for comparative workflows; both raw and corrected energies are returned so callers can choose policy.

## Programmatic Method Planning

sci-form now exposes a structured method-planning API that reports:

- recommended method per property domain
- fallback path when applicable
- confidence level and numeric confidence score
- structured limitations and warnings

Available entry points:

- Rust: `get_system_method_plan(elements)`
- Python: `system_method_plan(elements)`
- WASM: `system_method_plan(elements_json)`

This is intended to help UI and workflow layers choose between embedding, UFF-backed energy workflows, and EHT-backed orbital workflows without hard-coding metal-specific rules outside the library.

## Multi-Method Comparison Workflow

sci-form now exposes a comparison workflow that runs available methods on the same input geometry and returns per-method status, confidence metadata, warnings, limitations, and compact outputs.

Currently compared methods:

- UFF: force-field energy (kcal/mol)
- EHT: HOMO, LUMO, and HOMO-LUMO gap (eV)

Available entry points:

- Rust: `compare_methods(smiles, elements, positions, allow_experimental_eht)`
- Python: `compare_methods(smiles, elements, coords, allow_experimental_eht=False)`
- WASM: `compare_methods(smiles, elements_json, coords_json, allow_experimental_eht)`

## Structured Topology Analysis

sci-form now also exposes machine-readable topology analysis for transition-metal centers.

Current output includes:

- detected metal centers
- ligand atom indices assigned to each center
- coordination number
- inferred geometry and fit score

Supported motif detection currently includes:

- linear
- trigonal
- tetrahedral
- square-planar
- trigonal-bipyramidal
- octahedral

Available entry points:

- Rust: `compute_topology(elements, positions)`
- Python: `topology_analysis(elements, coords)`
- WASM: `compute_topology(elements_json, coords_json)`

## Known Limits

- Transition-metal parameters remain provisional and are not yet literature-calibrated across a broad benchmark set.
- The current regression threshold protects consistency over time, not absolute agreement with experiment.
- Target-level calibration for the benchmark set is still an open roadmap item.

## Roadmap Link

This page satisfies the roadmap requirement to document confidence and validity limits for metal chemistry while calibration is still in progress.