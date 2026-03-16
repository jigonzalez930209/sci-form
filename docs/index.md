---
layout: home

hero:
  name: sci-form
  text: Computational Chemistry in Rust
  tagline: Conformer generation, EHT, PM3, GFN-xTB, ML properties, force fields — Rust, Python, TypeScript, CLI (v0.4.3)
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
    details: "60+ conformers/second, parallel batch via rayon, zero C++ dependencies. ~2 MB binary."
  - icon: 🎯
    title: RDKit-Quality Accuracy
    details: "0.064 Å avg RMSD vs RDKit. 846 CSD torsion patterns. 100% stereo validation."
  - icon: 🔬
    title: Three QM Methods
    details: "**NEW v0.4** — EHT (all elements), PM3 (NDDO SCF, thermochemistry), GFN-xTB (25 elements, ultra-fast)."
  - icon: ⚛️
    title: Extended Hückel Theory
    details: "Wolfsberg-Helmholtz Hamiltonian, Löwdin orthogonalization, HOMO/LUMO gaps, population analysis, dipole moments."
  - icon: 🌊
    title: Electrostatic Potential & DOS
    details: "Coulomb ESP grids (red/white/blue mapping), total/per-atom DOS with Gaussian smearing, volumetric orbital grids, Marching Cubes isosurfaces."
  - icon: 🧬
    title: ML Descriptors & Properties
    details: "**NEW v0.4** — 17 descriptors (no 3D), LogP, solubility, Lipinski Ro5, druglikeness. SMILES → results in ~1 μs."
  - icon: 🧲
    title: Force Fields
    details: "UFF (50+ elements + TM) and MMFF94 (Halgren 14-7 vdW, quartic/cubic bends, Fourier torsions)."
  - icon: 📐
    title: Molecular Alignment
    details: "Kabsch SVD and quaternion-based Coutsias 2004. RMSD after optimal superposition."
  - icon: 🏗️
    title: Materials Assembly
    details: "Periodic unit cells, SBU topology, MOF-type framework crystal structure generation."
  - icon: 🌐
    title: Multi-Platform Bindings
    details: "Native Rust, PyO3 (Python), WASM (TypeScript/JS, browser/Node.js/Deno/Bun), and CLI."
  - icon: 📖
    title: Documented Algorithms
    details: "Complete theoretical foundations with mathematical derivations, SVG diagrams, and step-by-step pipeline explanations for every module."
---
