---
layout: home

hero:
  name: sci-form
  text: Computational Chemistry Toolkit in Rust
  tagline: ETKDG conformers, PM3 and GFN-xTB, HF-3c, ANI, spectroscopy, descriptors, and materials workflows across Rust, Python, TypeScript, and CLI
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
    title: Multi-method Quantum Stack
    details: "EHT, PM3/PM3(tm), GFN0/GFN1/GFN2-xTB, HF-3c, and CISD cover fast screening through corrected ab-initio workflows."
  - icon: 🌈
    title: Spectroscopy (Track D)
    details: "UV-Vis sTDA, vibrational analysis + broadened IR, and NMR shifts/couplings are available from the same molecular pipeline."
  - icon: ⚛️
    title: Extended Hückel Theory
    details: "Wolfsberg-Helmholtz Hamiltonian, Löwdin orthogonalization, HOMO/LUMO gaps, population analysis, dipole moments."
  - icon: 🌊
    title: Electrostatic Potential & DOS
    details: "Coulomb ESP grids (red/white/blue mapping), total/per-atom DOS with Gaussian smearing, volumetric orbital grids, Marching Cubes isosurfaces."
  - icon: 🧬
    title: ML Potentials & Descriptors
    details: "ANI-2x and ANI-TM potentials, topological descriptors, WHIM/RDF/GETAWAY, and built-in ML property models support fast screening."
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
    details: "Algorithm notes, API references, benchmark writeups, and WebGPU validation docs track the current repository surface."
---
