---
layout: home

hero:
  name: sci-form
  text: Computational Chemistry in Rust
  tagline: Conformer generation, EHT, PM3, GFN0/1/2-xTB, HF-3c, ANI-2x/TM, UV-Vis, IR, NMR, ML, stereochemistry, crystallography — Rust, Python, TypeScript, CLI (v0.10.6)
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
    details: "60+ conformers/second, parallel batch via rayon, zero C++ dependencies. ETKDGv2 with 846 CSD torsion patterns. ~2 MB binary."
  - icon: 🎯
    title: RDKit-Quality Accuracy
    details: "0.064 Å avg RMSD vs RDKit. 0.00% heavy-atom RMSD > 0.5 Å on GDB-20 (2000 molecules). 97.7% embed success on diverse molecules."
  - icon: 🔬
    title: Seven QM Methods
    details: "EHT (all elements + band structure + gradients), PM3 + PM3(tm) (transition metals), GFN0/GFN1/GFN2-xTB (25 elements, D3/D4), HF-3c (D3+gCP+SRB), CISD excited states."
  - icon: 🧠
    title: Neural Network Potentials
    details: "ANI-2x for H/C/N/O/F/S/Cl energies and forces. ANI-TM for 24 elements: H,C,N,O,F,Si,P,S,Cl,Ti,Cr,Mn,Fe,Co,Ni,Cu,Zn,Br,Ru,Pd,Ag,I,Pt,Au."
  - icon: 🌈
    title: Spectroscopy
    details: "UV-Vis sTDA with oscillator strengths. IR via numerical Hessian (3N modes, ZPVE, peak assignment, functional groups). ¹H/¹³C/¹⁹F/³¹P NMR shifts (HOSE codes), J-coupling (Karplus), ensemble averaging."
  - icon: 🤖
    title: Machine Learning
    details: "17+ 2D descriptors. 3D descriptors: WHIM (5 weightings), RDF, GETAWAY. LogP, solubility, Lipinski Ro5. Random Forest with variance, Gradient Boosting, cross-validation."
  - icon: 🧬
    title: Stereochemistry & Cheminformatics
    details: "CIP R/S, E/Z, atropisomeric M/P, helical chirality. SSSR ring perception. ECFP fingerprints + Tanimoto. Butina RMSD clustering. SMIRKS reaction transforms. Canonical SMILES."
  - icon: 💧
    title: Solvation & Reactivity
    details: "Non-polar SASA solvation. Generalized Born (HCT). ALPB. D4 dispersion + gradients. EEQ geometry-dependent charges. Fukui, frontier orbital descriptors, empirical pKa."
  - icon: 🏗️
    title: Materials & Crystallography
    details: "All 230 ITC space groups with symmetry operations. MOF-type framework assembly (pcu/dia/sql). BFGS+SD geometry optimization with PBC. Hapticity and metallocene detection."
  - icon: 🧲
    title: Force Fields & Alignment
    details: "UFF (50+ elements, TMs) and MMFF94 (Halgren 14-7 vdW). Kabsch SVD and quaternion alignment (Coutsias 2004). NPA/NBO population analysis. Bond orders (Wiberg, Mayer)."
  - icon: 🌐
    title: Multi-Platform Bindings
    details: "Native Rust crate (sci-form). PyO3 Python package (sciforma on PyPI). WASM TypeScript/JS module (sci-form-wasm on npm, browser/Node/Deno/Bun). Cross-platform CLI (sci-form-cli)."
  - icon: 🔭
    title: Experimental Engines
    details: "GPU acceleration (experimental-gpu): MMFF94 non-bonded + EEQ WGSL kernels. Alpha: CGA, GSM (reaction paths), SDR. Beta: KPM, MBH, CPM, RandNLA, Riemannian optimization."
---
