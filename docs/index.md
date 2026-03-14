---
layout: home

hero:
  name: sci-form
  text: 3D Molecular Conformer Generation
  tagline: High-performance ETKDG distance geometry from SMILES — Rust, Python, TypeScript, CLI
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
    details: 60+ molecules/second in native Rust. Parallel batch processing via rayon. Zero runtime dependencies.
  - icon: 🎯
    title: RDKit-Quality Accuracy
    details: 0.00% heavy-atom RMSD above 0.5 Å vs RDKit on GDB-20. 837 CSD torsion patterns for realistic geometry.
  - icon: 🌐
    title: Multi-Platform
    details: Native Rust library, Python (PyO3), TypeScript/JS (WASM), and cross-platform CLI for Linux, macOS, and Windows.
  - icon: 🧬
    title: Full Chemical Coverage
    details: Handles all functional groups, stereo centers, macrocycles, fused rings, metals, halogens, and heavy atoms.
  - icon: 📐
    title: Rigorous Validation
    details: Tetrahedral center checks, planarity enforcement, double-bond geometry verification, and chiral volume constraints.
  - icon: 📖
    title: Documented Algorithms
    details: Complete theoretical foundation with mathematical derivations, SVG diagrams, and step-by-step pipeline explanations.
---
