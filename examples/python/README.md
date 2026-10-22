# Python Example / Test Suite

Comprehensive examples for the `sci_form` Python bindings, ordered from
basic operations to advanced and experimental methods.

## Prerequisites

Install the all-features wheel (Python 3.14 uses the ABI3 forward-compatible path):

```bash
cd crates/python
PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1 maturin build --release --all-features
pip install --user --break-system-packages --force-reinstall \\
  target/wheels/sciforma-1.0.0-cp314-cp314-linux_x86_64.whl
```

## Run

```bash
cd examples/python
python3 run_all.py
```

Expected result: all 11 modules run without raising.

## Modules

- 01_embed.py — conformer generation, batch, SMILES parsing
- 02_charges_sasa_dipole.py — charges, SASA, dipole, population, DOS, RMSD
- 03_eht_forcefields.py — EHT, orbital meshes, bond orders, frontier, UFF/MMFF94, ANI, ESP
- 04_quantum_methods.py — PM3, xTB, GFN1, GFN2, HF-3c
- 05_reactivity_ml_stereo.py — Fukui, reactivity, pKa, ML descriptors/predict, stereo, rings, fingerprints
- 06_spectroscopy.py — sTDA UV-Vis, vibrational analysis, IR
- 07_nmr.py — NMR shifts, couplings, spectrum, GIAO, HOSE codes
- 08_solvation_materials.py — solvation, unit cells, frameworks
- 09_clustering_transport.py — clustering, fingerprints, transport, descriptors
- 10_alpha.py — DFT, ReaxFF, MLFF AEVs, ERI, CGA, GSM, SDR, reaction dynamics
- 11_beta.py — KPM, MBH, RandNLA, CPM, PSD geometry

Each script prints labelled inputs/outputs for the corresponding API.
