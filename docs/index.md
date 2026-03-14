---
layout: home

hero:
  name: sci-form
  text: Computational Chemistry in Rust
  tagline: 3D conformer generation, EHT, ESP, DOS, MMFF94, alignment, materials — Rust, Python, TypeScript, CLI
  image:
    src: /logo.svg
    alt: sci-form
  actions:
    - theme: brand
      text: Get Started
      link: /guide/getting-started
    - theme: alt
      text: Algorithm Deep Dive
      link: /algorithm/overview
    - theme: alt
      text: View on GitHub
      link: https://github.com/jigonzalez930209/sci-form

features:
  - icon: ⚡
    title: High Performance
    details: 60+ conformers/second in native Rust. Parallel batch processing via rayon. Zero runtime dependencies.
  - icon: 🎯
    title: RDKit-Quality Accuracy
    details: 0.00% heavy-atom RMSD above 0.5 Å vs RDKit on GDB-20. 846 CSD torsion patterns for realistic geometry.
  - icon: 🔬
    title: Extended Hückel Theory
    details: EHT Hamiltonian (Wolfsberg-Helmholtz), Löwdin orthogonalization, HOMO/LUMO gaps, Mulliken & Löwdin population analysis, dipole moments.
  - icon: 🌊
    title: Electrostatic Potential
    details: Coulomb ESP grid from Mulliken charges, red/white/blue color mapping, parallel rayon evaluation, Gaussian Cube file export.
  - icon: 📊
    title: Density of States
    details: Total DOS and per-atom PDOS with Gaussian smearing from EHT orbital energies. MSE metric, JSON export.
  - icon: 🧲
    title: Force Fields
    details: UFF (50+ element types including transition metals) and MMFF94 (Halgren 14-7 vdW, quartic stretch, cubic bend, 3-term torsion).
  - icon: 📐
    title: Molecular Alignment
    details: Kabsch SVD alignment and quaternion-based Coutsias 2004 method. Optimal rotation, RMSD computation after superposition.
  - icon: 🏗️
    title: Materials Assembly
    details: Periodic unit cells from lattice parameters, secondary building unit (SBU) topology, MOF-type framework crystal structure generation.
  - icon: 🌐
    title: Multi-Platform
    details: Native Rust library, Python (PyO3), TypeScript/JS (WASM) with typed-array APIs, and cross-platform CLI for Linux, macOS, Windows.
  - icon: 📖
    title: Documented Algorithms
    details: Complete theoretical foundations with mathematical derivations, SVG diagrams, and step-by-step pipeline explanations for every module.
---
