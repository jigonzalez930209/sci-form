# sci-form — Estado del Proyecto y Tareas Pendientes

## Ciclo actual — Kohn–Sham DFT en la librería

El DFT no se importa de otro paquete. El solver es `src/dft/` (`alpha-dft`): F = H + J + Vxc, SVWN y PBE, grilla de Becke.

- [x] Base mínima all-electron hasta Zn para los elementos que no están en la tabla Hehre. Antes caían en un solo 1s.
- [x] ERI cartesianos s/p (Obara–Saika) para el Coulomb KS; gate H₂/H₂O vs PySCF a 1 mHa.
- [x] ERI d-shell para el Coulomb de Sc–Zn (Obara–Saika cartesiano con memo; humo Zn átomo).
- [x] ECP de core pequeño para Z > 30 (Z_eff + potencial local Gaussiano; f en el core).
- [x] UKS para capas abiertas (Slater polarizado + VWN sobre ρ total).
- [x] Split-valence 3-21G (H, C, O) y gradiente KS frozen-P + Pulay.

El detalle está en `docs/algorithm/ks-dft.md`.

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
- UV-Vis: sTDA vertical excitations + spectral broadening
- Fluorescence: Kasha selection, spontaneous-emission lineshape, Einstein A lifetime (vertical; no S1 optimization)
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
- [ ] Reducir latencia single-conformer al target de 11 ms/mol — medición 2026-10-12 (release, llamada única): CCO 9.4 ms, benceno 5.7 ms, aspirina 15.2 ms, cafeína 25–31 ms (parse 0.05 + bounds 3–6 + embed ~24), colesterol 108 ms. El batch paralelo ya va a ~3.7 ms/mol (cafeína ×50). Perfilado bloqueado en esta máquina (`perf_event_paranoid`); siguiente paso: sesión dedicada con instrumentación por fases del pipeline ETKDG (bounds 15–20 %, embedding+refinamiento 80 %)
- [x] CLI batch stdin: streaming real en chunks de 512 (`embed_cmds::cmd_batch`), JSON/XYZ/SDF incrementales

### 2. HF / Post-HF
- [ ] Completar CISD singles-doubles coupling — dedicado; ruta de validación: `pyscf.ci.cisd` (disponible en el venv) con fixture nuevo
- [ ] Reducir coste de memoria GPU para ERI O(N⁴) — batching/tiling del tensor (ver plan §14.3: >200 funciones ≈ 3.2 GB en CPU)
- [ ] Mejorar momentos de transición más allá de la aproximación de monopolo: el fixture F41 (`pyscf_tda.json`) ya trae las fuerzas de oscilador de PySCF y el test las ignora (campo `oscillator` sin leer), así que el gate está listo; falta añadir integrales de dipolo AO al módulo `scf` (hoy no existen) y usarlas en `transition_dipole_from_ci``

### 3. PM3
- [x] MOPAC 23.2.5 descargado y usado como árbitro de F42 (fixture `tests/fixtures/validation/mopac_pm3.json`, 24 moléculas)
- [x] Núcleo PM3 fiel (port de MOPAC): tabla de parámetros de la fuente (42 elementos + eisol), overlaps analíticos `diat2` (casos 1–6), `reppd`/`rotmat`/`spcore`/`ccrep` (con la corrección N–H/O–H), `fock1`/`fock2` y HoF con eisol. Validado exacto (<0.02 kcal/mol) en hidruros X–H (C/N/O/F/Si/P/S/Cl), matriz H de agua y etileno (1e-4 eV), SCF de agua (eigenvalores y cargas exactos) y E_nuc por par (1e-4 eV)
- [ ] PM3 residual documentado: las integrales 2e "polares" (π) de la rotación `tx`/`rotatd` salen con el signo opuesto a MOPAC (error hasta ~130 kcal/mol en sistemas π: benceno, aromáticos). Siguiente paso: depurar contra el dump `HCORE` de MOPAC usando el bloque C–C alineado en z
- [x] Delimitar claramente PM3 vs PM3(tm) en la superficie pública (`Pm3Method` + `pm3_method`, flag `tm` en la tabla)

### 4. xTB
- [x] GFN1-xTB fiel a tblite 0.7.0 (H0 Slater con kpair/kshell/enscale, electrostática shell-resolved con tercer orden, repulsión efectiva, D3-BJ y corrección de halógeno); 8.1e-5 eV vs tblite en 18 moléculas (F44)
- [x] D3-BJ con interpolación de CN, C8 completo y ATM (validado vs s-dftd3, 2.8e-11 Ha; F46)
- [ ] UKS/unrestricted para capas abiertas — dedicado (spin-resolved SCC con densidad α/β y test vs tblite `uhf`); hoy el SCC de GFN1/GFN2 es restringido
- [ ] Decidir el destino de `solve_xtb` legacy: conservado solo para el Hessiano numérico y tests de DOS; la batería aromática ya usa GFN0 (F43). Reemplazo completo pendiente de migrar el Hessiano a GFN2
- [ ] Repulsión dependiente de coordinación más allá del esquema efectivo de tblite (sin referencia externa por ahora)

### 5. EHT
- [x] Soporte real de orbitales 4d/5d (contracciones STO-3G propias 4d/5d en `eht::basis`, tests de forma y de Mo)
- [ ] Revisión de parámetros provisionales para metales de transición — pendiente de un benchmark contra tablas EHT de literatura (hoy: tabla d-block + fallback hydrogenic con warnings, ζ vs IP testeado)

### 6. NMR
- [ ] Ampliación fuerte de HOSE database — bloqueado: requiere una base externa de shifts (NMRShiftDB/SDBS) con licencia y validación
- [ ] Constantes de Karplus más específicas por entorno — bloqueado: requiere implementar Haasnoot-de Leeuw-Altona y un set experimental de ³J para validar (hoy hay Karplus estándar por camino C-C/C-N/C-O/C-S, con tests cualitativos)
- [x] Aromaticidad real con SSSR/Hückel para shifts (percepción antes de HOSE; Kekulé == aromático en los tests)
- [ ] Salida más allá del supuesto first-order (strong coupling) — requiere formalismo de segundo orden y espectros de referencia para gate

### 7. GPU
- [x] Fallback CPU transparente cuando GPU no disponible (`compute_mmff94_nonbonded_with_report` + reportes de backend)
- [x] Tiling/sparse batching para MMFF94 (`mmff94_row_tiles`, evaluación tiled verificada contra la no tiled)
- [x] Tuning dinámico de batch/workgroup según device (`optimal_1d_dispatch`/`optimal_grid_dispatch` en MMFF94 y ESP)
- [x] Mensajes útiles para hardware con poca VRAM (`MemoryError::hint` con acciones sugeridas)

### 8. Espectroscopía
- [x] Prescreening real (Cauchy-Schwarz sobre normas de cargas de transición + ventanas de espacio activo; test de que el umbral anula el acoplamiento no diagonal)
- [x] Clasificar excitaciones dark/bright de forma más explícita (fluorescencia: umbral `dark_threshold` y bandera de Kasha)

### 9. SMARTS / SMIRKS
- [x] Recursive SMARTS (`$()`): parser y matcher implementados y validados contra RDKit (anillos, carbonilos, negaciones y patrones ramificados en `smarts::matcher`)
- [x] Endurecer validación de atom maps antes de aplicar transformaciones (mapas duplicados, `:0`, componentes vacíos y patrones SMARTS inválidos en `parse_smirks`; elementos de dos letras `Mg`/`As`/`se` en SMARTS)

### 10. Solvación
- [x] Unificar SASA entre módulo principal y ALPB (ambos usan `surface::sasa`, test de equivalencia)
- [x] Separación clara entre GB y ALPB en la API (`compute_gb_solvation` vs `compute_alpb_solvation`, documentados)

### 11. Materiales
- [x] Normalización automática de coordenadas fraccionales (el parser CIF envuelve a [0,1) con `UnitCell::wrap_frac`)
- [ ] Expansión más allá de topologías homogéneas (nets con aristas no equivalentes) — requiere más datos RCSR/pormake y ampliar el oráculo F67

### 12. Force Fields
- [x] MMFF94 cobertura para metales de transición: evaluado — el set oficial de MMFF94 no parametriza TMs; los metales se enrutan a UFF/GFN (documentado en el módulo)
- [x] Detector aromático independiente (Hückel en `rings::perceive_aromaticity`, aplicado antes del force field)
- [x] Breakdown por términos energéticos en la API pública (`compute_uff_energy_breakdown`, `compute_mmff94_energy_breakdown`)

### 13. Dispersion
- [ ] Completar parametrización/tablas D4 C6 hasta cobertura amplia
- [x] Descomposición energética visible cuando se combine con otros métodos (GFN0 2/3-cuerpos, GFN1 e6/e8, GFN2 2-cuerpos vs ATM, D4 standalone e2/e3)

### 14. Otros Módulos Pendientes Menores
- [x] ESP: distancia mínima configurable (`EspGridConfig::min_distance`) y validación estricta de spacing/padding (`compute_esp_grid_with_config`, `compute_esp_with_config`)
- [x] Distance Geometry: embeddings 2D/planos se re-perturban y re-minimizan en el pipeline; bounds contradictorios detectados por `triangle_smooth_tol` (escalera de reintentos + tests)
- [x] EEQ: limpiar restos de experimentalidad en nombres de API (comando CLI y bindings WASM sin feature gate)
- [ ] EEQ: cEEQ/EEQ-BC — bloqueado por falta de referencia (no existe en xtb 6.7.1 ni tblite; requeriría el paper del modelo y sus tablas propias)
- [x] Population: warning de cobertura (elementos sin modelo de valencia / sin parámetros EHT) y homo/lumo/gap expuestos en `PopulationResult` (Rust y binding Python)
- [x] Fukui: Hirshfeld/CM5 evaluado (requieren densidades promoleculares/cargas CM5 ausentes; NPA-condensado requiere EHT con carga, hoy no expuesto) — se mantiene Mulliken con la caveat documentada
- [x] Surface: tabla de radios documentada (`VDW_FALLBACK_RADIUS`, `vdw_radius_source`) y evaluación Connolly/SES documentada (diferida hasta que un consumidor necesite la malla)
- [x] Optimization: unificar optimizadores entre EHT/materials (EHT usa el BFGS compartido de `materials::geometry_opt`) y repulsión nuclear incluida en el gradiente EHT
- [ ] Optimization: variable-cell BFGS — pendiente de un tensor de stress en el modelo periódico (el callback energía/fuerzas no lo expone)

### 15. Reaction Dynamics 3D (alpha)
- [x] Deprecar `dynamics.rs` legacy como assembly primario — `dynamics.rs` conserva solo energía/NEB; el assembly vive en `alpha::reaction_dynamics::complex` (documentado en la cabecera del módulo)
- [x] Rotación Kabsch para >2 fragmentos (`assemble_fragments_at_product_positions`)
- [x] Posicionamiento de moléculas extra (≥2) según átomos equivalentes en el producto
- [x] Optimizador global de orientación multi-fragmento (`optimize_reactive_complex_global`: Kabsch + rotaciones aleatorias + relajación, se guarda la mejor; tests de cota)
- [ ] Comparación sistemática de ángulos de ataque, barreras y trayectorias vs referencias externas (requiere set de reacciones con barreras de literatura)
- [ ] Frontend: eliminar V-path X-axis fallback (código del frontend, fuera de este repo)

### 16. Alpha Export Gaps (bindings faltantes)
- [ ] Render Bridge: 13 funciones faltantes (arrhenius_chart, kpoint_path_chart, band_structure_chart, trajectory_chart, dos_chart, etc.) — dedicado; son wrappers de bindings sobre APIs Rust existentes
- [ ] Periodic Linear: 10 funciones faltantes (band_structure_eht, fermi_surface_from_bands, kpath_standard, etc.) — dedicado
- [ ] Kinetics CLI: 3 funciones faltantes — dedicado

### 17. Validación a Gran Escala (Track G)
- [x] Banco ≥100 moléculas (180 diversas + 1000 ChEMBL + 400 UFF; guard `track_g_bank_has_enough_molecules`)
- [x] Pipeline automatizado de benchmarking vs PySCF/RDKit (oráculos + `run_validation.sh` + `collect.py`)
- [x] Threshold <1 % documentado por observable (plan §14.2; los gated quedan 2–5 órdenes por debajo)
- [x] Límites de escalado documentados (plan §14.3: tiempos por método y memoria ERI O(N⁴))

### 18. Experimental Tracks Pendientes de Validación para Promoción
- [ ] Validación cuantitativa por submódulo (alpha + beta) — parcial: los promovidos tienen fichas (EEQ/D4/ALPB/sTDA); E1–E11 siguen con evidencia parcial
- [ ] Benchmarking frente a rutas estables (KPM vs diagonalización, GSM vs NEB, SDR vs Butina)
- [x] Criterios claros de promoción experimental → producción (`src/experimental_status.rs`, actualizado con la evidencia actual)
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

Alignment, ANI, Dipole, DOS, ML descriptors, Rings/SSSR, Transport, SCF (UHF/ROHF), CIF, AO→MO transform, EEQ damping Z=86, Population Z=86+parallel, NMR 5J coupling, SMIRKS multi-component, PM3 Gaussian corrections, xTB Broyden mixing, GFN1-xTB (F44), D3-BJ (F46), tablas STO-6G completas.
