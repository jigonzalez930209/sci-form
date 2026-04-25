# sci-form — Estado del Proyecto y Tareas Pendientes

> Última revisión: 2026-04-07
> Consolidación de: ROADMAP.md, ROADMAP_EXPERIMENTAL.md, ROADMAP_ALPHA_*.md, ROADMAP_REACTION_DYNAMICS_3D.md, todo.md, SESSION_5_*.md, SMIRKS_*.md, sci-form-exportation.md, dynamics-roadmap-implementation.md, uv-vis-spectroscopy-findings.md

---

## Implementado y Funcionando

### Core — Conformer Engine (Track A)
- SMILES parser completo (60+ bracket elements, implicit H, stereo, isotopes)
- Distance geometry: bounds matrix 1-2/1-3/1-4, shortest-path smoothing, metric matrix → 3D
- ETKDG refinement: ETKDGv2/v3 torsion knowledge, macrocycles, amide/sulfonamide/phosphonate patterns
- UFF + MMFF94 force fields con gradientes analíticos y L-BFGS minimization
- Batch embedding paralelo (`rayon`), typed-array WASM, Arrow transport, Web Workers
- ChemBL 10K: 97.54% success, 97.18% perfect geometry, ~2.1 mol/s
- GDB-20 validation: 0.024 Å avg RMSD vs RDKit, 0.00% > 0.5 Å

### Electronic Structure
- EHT (Extended Hückel Theory): STO-3G basis, Wolfsberg-Helmholtz, volumetric orbital grids, Marching Cubes isosurfaces, gradients + geometry optimization
- EHT Band Structure: k-path periodic Bloch Hamiltonian
- PM3 (NDDO SCF): heat of formation, HOMO/LUMO, Mulliken charges, Gaussian core-core corrections (main-group + PM3(tm) transition metals Ti–Au)
- GFN0-xTB: SCC tight-binding, Broyden mixing, 25+ elements incl. transition metals
- GFN1-xTB: shell-resolved charges + D3 dispersion, Broyden SCC
- GFN2-xTB: multipole electrostatics + D4 + halogen bonding
- HF-3c: minimal-basis HF + D3-BJ + gCP + SRB corrections, CIS UV-Vis
- UHF / ROHF: level shift, damping, spin contamination analysis, α/β orbitals
- CISD: excited states from CIS+D, AO→MO 4-index integral transform
- ANI-2x: neural network potential (H,C,N,O,F,S,Cl), AEV + backprop forces
- ANI-TM: extended to 24 elements including transition metals

### Properties & Analysis
- Gasteiger-Marsili charges, EEQ charges (Z=1–86, improved Gaussian damping)
- Mulliken & Löwdin population (parallel, Z=1–86), bond orders (Wiberg/Mayer)
- NPA/NBO: Natural Population & Bond Orbital analysis
- Dipole moment (Debye), ESP grid (parallel), DOS + PDOS (multi-method)
- SASA (Shrake-Rupley), RMSD (Kabsch SVD + quaternion)
- Fukui descriptors, reactivity ranking, empirical pKa, frontier descriptors
- ML descriptors (17+), LogP/MR/solubility/druglikeness prediction
- WHIM, RDF, GETAWAY 3D molecular descriptors
- Random Forest, Gradient Boosting, cross-validation
- Stereochemistry: CIP R/S, E/Z, atropisomeric M/P, helical chirality
- SSSR (Horton's algorithm), ECFP fingerprints + Tanimoto
- Butina RMSD clustering, diversity filtering

### Spectroscopy
- UV-Vis: sTDA-xTB vertical excitations + spectral broadening
- IR: numerical Hessian, vibrational frequencies, dipole intensities, peak assignment, RRHO thermochemistry
- NMR: HOSE-code 1H/13C shifts, J-coupling (Karplus, 2J–5J incl. long-range), ensemble averaging

### Materials & Periodics
- Unit cell, MOF framework assembly (pcu/dia/sql topologies)
- 230 ITC space groups, equivalent position generation
- CIF import/export (parser + writer, uncertainty notation)
- Periodic molecular graphs, hapticity/metallocene detection
- Framework geometry optimization (BFGS + steepest descent, PBC)

### Reactions
- SMIRKS: parse + apply single/multi-component reactions
- NEB multi-method (UFF/MMFF94/PM3/xTB/GFN1/GFN2/HF-3c)
- MD: Velocity Verlet, Nosé-Hoover thermostat

### Solvation
- Non-polar SASA solvation, Generalized Born (HCT) electrostatic

### Bindings
- Rust crate, Python (PyO3/maturin), WASM (wasm-bindgen), CLI
- EDL module: 100% across all bindings (6/6 functions × 4 bindings)

### GPU (behind `experimental-gpu`)
- GPU MMFF94 non-bonded, GPU sTDA J-matrix, GPU Hessian, GPU EEQ, GPU CPM, GPU ALPB Born radii

### Alpha/Beta Infrastructure
- CI-NEB (climbing image), orbital-guided approach, constrained optimization, electrostatic steering, IRC, angular sampling (`src/alpha/reaction_dynamics/`)
- GSM (Growing String Method), MBH (Mobile Block Hessian), HTST kinetics
- KPM (Kernel Polynomial Method), RandNLA eigensolvers
- EEQ dynamic charges, ALPB solvation, D4 dispersion, CPM (Constant Potential)
- CGA (Conformal Geometric Algebra), Riemannian optimization, SDR embedding
- EDL (Electrical Double Layer) profiles + capacitance charts
- Render bridge infrastructure, transport/columnar data

---

## Tareas Pendientes

### 1. Core Engine — Rendimiento
- [ ] Reducir latencia single-conformer de ~106 ms/mol al target de 11 ms/mol
- [ ] CLI batch stdin: implementar streaming real (hoy materializa toda la entrada antes de procesar)

### 2. HF / Post-HF
- [ ] Completar CISD singles-doubles coupling
- [ ] Reducir coste de memoria GPU para ERI O(N⁴)
- [ ] Mejorar momentos de transición más allá de la aproximación de monopolo

### 3. PM3
- [ ] Completar parametrización residual: resonance integrals, nuclear repulsion, referencias atómicas para heat of formation
- [ ] Delimitar claramente PM3 vs PM3(tm) en la superficie pública

### 4. xTB
- [ ] Corrección completa D3-BJ C8 (actualmente simplificada)
- [ ] Refinamientos: gamma shell-resolved, SCC damping, criterio de convergencia de cargas
- [ ] Repulsión dependiente de coordinación

### 5. EHT
- [ ] Soporte real de orbitales 4d/5d (hoy usa tabla 3d STO-3G)
- [ ] Revisión de parámetros provisionales para metales de transición

### 6. NMR
- [ ] Ampliación fuerte de HOSE database
- [ ] Constantes de Karplus más específicas por entorno químico
- [ ] Aromaticidad real con SSSR/Hückel para shifts
- [ ] Salida más allá del supuesto first-order (strong coupling)

### 7. GPU
- [ ] Fallback CPU transparente cuando GPU no disponible
- [ ] Tiling/sparse batching para MMFF94
- [ ] Tuning dinámico de batch/workgroup según device
- [ ] Mensajes útiles para hardware con poca VRAM

### 8. Espectroscopía
- [ ] Prescreening real (Cauchy-Schwarz o similar) para reducir cuello de botella O(N⁴) en sTDA
- [ ] Clasificar excitaciones dark/bright de forma más explícita

### 9. SMARTS / SMIRKS
- [ ] Recursive SMARTS (`$()`)
- [ ] Endurecer validación de atom maps antes de aplicar transformaciones

### 10. Solvación
- [ ] Unificar SASA entre módulo principal y ALPB
- [ ] Separación clara entre GB y ALPB en la API

### 11. Materiales
- [ ] Normalización automática de coordenadas fraccionales
- [ ] Expansión más allá de topologías homogéneas

### 12. Force Fields
- [ ] MMFF94 cobertura completa para metales de transición
- [ ] Detector aromático independiente
- [ ] Breakdown por términos energéticos en la API pública

### 13. Dispersion
- [ ] Completar parametrización/tablas D4 C6 hasta cobertura amplia
- [ ] Descomposición energética visible cuando se combine con otros métodos

### 14. Otros Módulos Pendientes Menores
- [ ] ESP: hacer configurable la distancia mínima (0.01 Å) y endurecer validación spacing/padding
- [ ] Distance Geometry: detectar embeddings 2D/planos y re-perturbar; detección temprana de bounds contradictorios
- [ ] EEQ: implementar cEEQ y limpiar restos de experimentalidad en nombres de API
- [ ] Population: emitir warning cuando falte cobertura; exponer homo/lumo/gap en resultado público
- [ ] Fukui: evaluar Hirshfeld/CM5 para mayor robustez (actualmente basis-set-dependent via Mulliken)
- [ ] Surface: documentar tabla de fallback radii; valorar opción tipo Connolly surface
- [ ] Optimization: unificar optimizadores entre EHT/materials; variable-cell BFGS

### 15. Reaction Dynamics 3D (alpha)
- [ ] Deprecar `dynamics.rs` legacy como assembly primario
- [ ] Rotación Kabsch para >2 fragmentos
- [ ] Posicionamiento de moléculas extra (≥2) según átomos equivalentes en el producto
- [ ] Optimizador global de orientación multi-fragmento
- [ ] Comparación sistemática de ángulos de ataque, barreras y trayectorias vs referencias externas
- [ ] Frontend: eliminar V-path X-axis fallback (rotateToAlign ±X hardcodeados)

### 16. Alpha Export Gaps (bindings faltantes)
- [ ] Render Bridge: 13 funciones faltantes (arrhenius_chart, kpoint_path_chart, band_structure_chart, trajectory_chart, dos_chart, etc.)
- [ ] Periodic Linear: 10 funciones faltantes (band_structure_eht, fermi_surface_from_bands, kpath_standard, etc.)
- [ ] Kinetics CLI: 3 funciones faltantes

### 17. Validación a Gran Escala (Track G)
- [ ] Curar banco de ≥100 moléculas representativas diversas
- [ ] Pipeline automatizado de benchmarking vs PySCF/RDKit
- [ ] Threshold de <1% desviación en energías, gradientes y propiedades
- [ ] Documentar límites de escalado (tiempo + memoria)

### 18. Experimental Tracks Pendientes de Validación para Promoción
- [ ] Validación cuantitativa por submódulo (alpha + beta)
- [ ] Benchmarking frente a rutas estables
- [ ] Criterios claros de promoción experimental → producción
- [ ] Tests de integración más representativos para beta

---

## Tracks Experimentales (Resumen de Estado)

| Track | Módulo | Estado |
|-------|--------|--------|
| E1 CGA | `experimental::cga` | Implementado, pendiente validación |
| E2 RandNLA | `beta::rand_nla` | Implementado, pendiente benchmarks a escala |
| E3 Riemannian | `experimental::riemannian` | Implementado, baja prioridad (E8 resuelve lo mismo) |
| E4 KPM | `beta::kpm` | Implementado, pendiente validación vs exact para N>100 |
| E5 EEQ Dynamic | `experimental::eeq` / `charges_eeq` | Implementado, parcialmente integrado |
| E6 ALPB | `experimental::alpb` / `solvation_alpb` | Implementado, pendiente unificación con GB |
| E7 D4 Dispersion | `experimental::d4` | Implementado, pendiente cobertura de tablas |
| E8 SDR Embedding | `experimental::sdr` | Implementado, pendiente benchmark 10K |
| E9 MBH | `beta::mbh` | Implementado, pendiente validación completa de frecuencias |
| E10 CPM | `beta::cpm` | Implementado, pendiente benchmarks electroquímicos |
| E11 GSM | `alpha::gsm` | Implementado, pendiente validación vs NEB en 5+ reacciones |

---

## Módulos Cerrados (sin deuda relevante)

Alignment, ANI, Dipole, DOS, ML descriptors, Rings/SSSR, Transport, SCF (UHF/ROHF), CIF, AO→MO transform, EEQ damping Z=86, Population Z=86+parallel, NMR 5J coupling, SMIRKS multi-component, PM3 Gaussian corrections, xTB Broyden mixing.
