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
- [x] Latencia single-conformer: perfilado desbloqueado (`perf_event_paranoid=1`) y ejecutado con `examples/latency_probe.rs` + `perf record` (release con símbolos). **Fases medidas**: `minimize_bfgs_rdkit` (bounds) 20 %, `minimize_etkdg_3d_bfgs` 21 %, energía+gradiente ETKDG 17 %, SMARTS `atom_matches` 8 %, SSSR 2.5 %, `min_path_excluding` 2.5 % — confirma el reparto hipotetizado (bounds ~20 %, embedding+refinamiento ~60 %). **Optimizaciones aplicadas**: visited de los BFS de `graph.rs` con `Vec<bool>` (era ~10 % del tiempo en hashing SipHash) y filtro de relevancia de anillos con merge de claves ordenadas. **Resultado (release, CLI proceso por llamada)**: CCO 8 ms, benceno 7 ms, aspirina 11 ms, cafeína 13 ms (media 9 ms) — objetivo de 11 ms/mol alcanzado; warm in-process cafeína 5.4 ms. Próxima palanca si se quiere más: los dos minimizadores BFGS (43 % self + 17 % su energía/gradiente)
- [x] CLI batch stdin: streaming real en chunks de 512 (`embed_cmds::cmd_batch`), JSON/XYZ/SDF incrementales

### 2. HF / Post-HF
- [x] CISD completo — `hf::cisd::compute_cisd` es ahora un CISD determinantal (orbitales de espín) sobre la referencia RHF: incluye los dobles de mismo orbital (i²→a² vía (iα iβ)→(aα aβ), el determinante dominante de H2/STO-3G), el acoplamiento completo con la referencia (reglas de Slater–Condon con las integrales monoelectrónicas derivadas de la condición de Fock canónica) y la energía de correlación real (`λ_min − E_ref`, antes el menor autovalor CIS). Gate `f41_cisd_matches_pyscf` en `cross_validation.rs` (fixture `pyscf_tda.json::cisd_correlation_hartree`, generado por `oracle_pyscf.py`): 8 moléculas STO-3G a ≤ 2.1e-8 Ha (tolerancia 1e-6). La matriz se valida además contra una referencia de segunda cuantización explícita en `tests/cisd_brute.rs` (0 discrepancias > 1e-9 en 361 determinantes sintéticos) y el camino perturbativo (CIS+MP2) contra `pyscf.mp.MP2`. Hallazgo clave: las creaciones de Slater–Condon van en orden descendente (`a†_a a†_b a_j a_i`); el signo invertido era invisible en H2O pero rompía CH2O/C2H4 (8.5e-3 Ha)
- [x] Memoria GPU para ERI O(N⁴): dispatch por tiles (`eri_quartet_tiles` + `compute_eris_gpu_tiled`) con salida *packed* (un f32 por cuarteto único) y desempaquetado 8-fold en el host; los buffers GPU quedan acotados por `GpuMemoryLimits` y bases que antes excedían el límite ya se pueden procesar. Tests: planificador (cobertura contigua, límite por tile) y round-trip packing↔tensor CPU (1e-12). Nota: el almacenamiento CPU del tensor completo (plan §14.3) sigue siendo n⁴ f64
- [x] Momentos de transición con integrales de dipolo AO: `scf::dipole_matrix` (recursión 1D), camino exacto en `hf::cis::compute_cis_with_ao_dipoles` y en `spectroscopy::compute_stda_with_dipoles`; gate F41 extendido: |ΔD| = 1.1e-9 a.u. vs `int1e_r` de PySCF y |Δf| < 1e-5 en 27/29 estados (H2O 0/1 documentado: intercambio de carácter entre configuraciones casi degeneradas)`

### 3. PM3
- [x] MOPAC 23.2.5 descargado y usado como árbitro de F42 (fixture `tests/fixtures/validation/mopac_pm3.json`, 24 moléculas)
- [x] Núcleo PM3 fiel (port de MOPAC): tabla de parámetros de la fuente (42 elementos + eisol), overlaps analíticos `diat2` (casos 1–6), `reppd`/`rotmat`/`spcore`/`ccrep` (con la corrección N–H/O–H), `fock1`/`fock2` y HoF con eisol. Validado exacto (<0.02 kcal/mol) en hidruros X–H (C/N/O/F/Si/P/S/Cl), matriz H de agua y etileno (1e-4 eV), SCF de agua (eigenvalores y cargas exactos) y E_nuc por par (1e-4 eV)
- [x] PM3 residual π — **resuelto** (2026-10-13): la rama heavy–heavy de `reppd` no reasignaba `qa`/`qb` a las separaciones cuadrupolares *simples* antes de `arg(54..72)` (MOPAC `reppd`: `qa = qq(ni); qb = qq(nj)`); los valores doblados corrompían exactamente las integrales polares (dxqxz/qxzdx/qxzqxz, `ri(10/15/20)`). Gate `f42_mopac_pm3` endurecido: las 24 moléculas del fixture a <0.02 kcal/mol (máx 8.4e-3, antes cota 130)
- [x] Delimitar claramente PM3 vs PM3(tm) en la superficie pública (`Pm3Method` + `pm3_method`, flag `tm` en la tabla)

### 4. xTB
- [x] GFN1-xTB fiel a tblite 0.7.0 (H0 Slater con kpair/kshell/enscale, electrostática shell-resolved con tercer orden, repulsión efectiva, D3-BJ y corrección de halógeno); 8.1e-5 eV vs tblite en 18 moléculas (F44)
- [x] D3-BJ con interpolación de CN, C8 completo y ATM (validado vs s-dftd3, 2.8e-11 Ha; F46)
- [x] UKS/unrestricted para capas abiertas — **implementado** (2026-10-13): `xtb::uks::{solve_gfn2_unrestricted, solve_gfn1_unrestricted}` con SCC spin-resolved (densidades α/β, magnetización por shell, interacción on-site de tblite `spin.f90`/`data/spin.f90` transcrita a `spin_params.rs` para Z=1–86). Hallazgo: tblite 0.7.0 no conecta el contenedor de espín (su `uhf` solo desdobla ocupaciones; ciclos idénticos al restringido en dobletes), así que el gate compara `E_uks − E_spin` vs `tblite uhf` (física no-espín exacta salvo relajación de 2º orden). Gates: `f45_uks_matches_tblite_uhf` (GFN2, máx 3.5e-2 eV), `f44_uks_matches_tblite_uhf` (GFN1, máx 3.4e-2 eV), `uks_uhf0_matches_restricted` (uhf=0 ≡ restringido a 1e-4 eV) + conservación exacta Σm=uhf y E_spin<0. Oráculo `oracle_uks_tblite.py` (¡posiciones tblite-Python en BOHR!) y fixture `uks_tblite.json` (OH/CH3/NO/O2/H2O)
- [x] Destino de `solve_xtb` legacy — **decidido y migrado** (2026-10-13): los consumidores de energía usan ya los solvers validados (`HessianMethod::Xtb` → `solve_gfn2`, coherente con el `compute_xtb` público = GFN2 y con los Hessianos de referencia F48 «xtb 6.7.1 GFN2-xTB»; dipolo vibracional → GFN2; `DosMethod::Xtb` → GFN2; mesh orbital → GFN2; backends GSM ya usaban `compute_xtb` = GFN2). El solver legacy (`solver.rs`) se conserva **solo** como proveedor del gradiente analítico de los backends de fuerza `"xtb"` (`dynamics.rs` NEB/MD y `compute_xtb_gradient` público, sin gradientes GFN0/GFN2 analíticos en el repo); su eliminación total queda pendiente de esos gradientes (dedicado)
- [x] Repulsión dependiente de coordinación — **cerrado como limitación documentada**: nuestro GFN1/GFN2 reproduce el esquema efectivo de tblite 0.7.0 (validado término a término, F44/F45); ir «más allá» exigiría una referencia externa que hoy no existe (tblite/xtb no exponen una repulsión CN-dependiente explícita). Ruta de desbloqueo: si una versión futura de tblite/xtb publica el esquema explícito, portarlo y extender el gate F44/F45

### 5. EHT
- [x] Soporte real de orbitales 4d/5d (contracciones STO-3G propias 4d/5d en `eht::basis`, tests de forma y de Mo)
- [x] Revisión de parámetros para metales de transición — **benchmark hecho** contra la tabla de referencia de YAeHMOP 3.0.3 (`tightbind/eht_parms.dat`, compilación de Alvarez 1993) con `scripts/validation/check_eht_metal_params.py`: 33 valores no coincidían con la fuente citada (ζ de 4d/5d en Zr/Nb/Mo/Tc/Ru/W/Re/Os/Ta/Ir/Cd, VSIP de Ni 3d, Zr, Tc, Cd, ζ de Mn/Zn) y se corrigieron a la literatura; gate `f31_metal_params_match_yaehmop` con Δ = 0 en 76 filas (tolerancia 0.02). Límites: Y/Ag/Hf conservan nuestros valores (filas incompletas en la referencia) y la base sigue siendo single-zeta (primer exponente de las d double-zeta). Siguen siendo «provisional» en el sentido de que el *rendimiento* químico no está validado contra datos experimentales, solo los parámetros contra la tabla

### 6. NMR
- [x] HOSE database — **cerrado como bloqueado con ruta concreta**: la implementación actual cubre la base de shifts integrada y está validada con tests cualitativos; una ampliación fuerte requiere una base externa (NMRShiftDB2 o SDBS) más una revisión de licencia y un protocolo de validación cruzada (shifts experimentales por entorno). Ruta de desbloqueo: descargar el dump de NMRShiftDB2, revisar licencia (CC-BY-SA para los datos), generar códigos HOSE de la base y validar contra shifts experimentales con tolerancia documentada
- [x] Constantes de Karplus específicas por entorno — **HLA implementado** (`nmr::coupling::{hla_coupling, huggins_dchi}`: `P1cos²φ + P2cosφ + P3 + Σ Δχ_i[P4 + P5cos²(ξ_iφ + P6|Δχ_i|)]` con P1=13.893, P2=−0.981, P4=0.637, P5=−2.243, P6=5.0 y Δχ de Huggins ajustados) validado contra la implementación de referencia (calculadora HLA de Stenutz, que implementa Haasnoot et al. 1980) con `scripts/validation/oracle_hla.py` (84 casos) y gate F54: base a 4e-3 Hz, sustituyentes C/S/N/Cl a 0.42 Hz. Convención ξ calibrada (S1/S3→+1, S2/S4→−1). Límite documentado: O/F/Si/Br/I/P requieren las ecuaciones específicas por patrón (C/D/E) de la referencia, y sigue pendiente la validación contra un set *experimental* de ³J (el gate usa una implementación de referencia, no datos de medida)
- [x] Aromaticidad real con SSSR/Hückel para shifts (percepción antes de HOSE; Kekulé == aromático en los tests)
- [x] Salida más allá del supuesto first-order: `nmr::second_order::SpinSystem` diagonaliza el Hamiltoniano de espín (base producto, hasta 6 espines) y da el espectro exacto de líneas; `strong_coupling_report` marca pares con |Δν|/|J| < 10 y `compute_nmr_spectrum` los anota. Gate: solución analítica AB (posiciones ν_c ± (C±|J|)/2 e intensidades 1 ± sin2θ) a 1e-8, límite débil (dobletes 1:1), A2 (singlete) y roof effect. La fórmula de referencia se corrigió contra la diagonalización (intensidades 1 ± J/C, no 1 ± Δν/C)

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
- [x] Topologías más allá de las homogéneas: `materials::topology_data` importa 6 nets de la base RCSR/pormake (acs, afw, sod, pts, xah, tfy) con tipos de nodo y — el punto del ítem — **aristas no equivalentes** (afw/xah/tfy con 2 tipos de arista, vía `# EDGE_CENTER`), generadas con `gen_topologies.py` desde los `.cgd`. Oráculo F67 ampliado (fixture con 15 nets) y gate nuevo `f67_extended_nets_match_pormake` (nodos, coordinaciones, coordenadas y tipos de arista vs pormake). Pendiente: que `assemble_framework` consuma nets con múltiples tipos de arista (hoy el ensamblado usa los 4 built-ins homogéneos)

### 12. Force Fields
- [x] MMFF94 cobertura para metales de transición: evaluado — el set oficial de MMFF94 no parametriza TMs; los metales se enrutan a UFF/GFN (documentado en el módulo)
- [x] Detector aromático independiente (Hückel en `rings::perceive_aromaticity`, aplicado antes del force field)
- [x] Breakdown por términos energéticos en la API pública (`compute_uff_energy_breakdown`, `compute_mmff94_energy_breakdown`)

### 13. Dispersion
- [x] D4 C6 hasta cobertura amplia: tablas regeneradas desde dftd4 v4.3.0 (`reference.inc` + `data/*.f90`) hasta Z=118 (`gen_d4_refs.py`, verificando 1–86 contra la tabla previa) y EEQ2019 extendida a Z=103 desde multicharge (`gen_eeq2019_params.py`). Gate F46 ampliado con 7 casos pesados (PbI2, Au2, HgCl2, LaH3, Rn2, UO2, LrH): exactos vs dftd4 (≤1e-9 Ha). Límite documentado: Z=104–111 no soportado por dftd4 (nosotros devolvemos 0) y Z=112–118 usan el modelo de carga nuevo de dftd4 (EEQ-BC, fuera de alcance)
- [x] Descomposición energética visible cuando se combine con otros métodos (GFN0 2/3-cuerpos, GFN1 e6/e8, GFN2 2-cuerpos vs ATM, D4 standalone e2/e3)

### 14. Otros Módulos Pendientes Menores
- [x] ESP: distancia mínima configurable (`EspGridConfig::min_distance`) y validación estricta de spacing/padding (`compute_esp_grid_with_config`, `compute_esp_with_config`)
- [x] Distance Geometry: embeddings 2D/planos se re-perturban y re-minimizan en el pipeline; bounds contradictorios detectados por `triangle_smooth_tol` (escalera de reintentos + tests)
- [x] EEQ: limpiar restos de experimentalidad en nombres de API (comando CLI y bindings WASM sin feature gate)
- [x] EEQ: cEEQ/EEQ-BC — port término a término en `src/charges_eeq/eeqbc.rs` (CN erf + CN dirigida ponderada por EN de Pauling, matriz de capacitancias de enlace con vdW por pares, A con screening `C`, `qloc` con overrides de actínidos) con tablas generadas por `scripts/validation/gen_eeqbc_params.py` (`eeqbc2025.f90` + mctc-lib `vdwrad`/`paulingen`). Oráculo **directo** (mejor que la referencia indirecta prevista): CLI de multicharge compilado con meson/ninja (`--model eeqbc`), fixture de 23 moléculas (`tests/fixtures/validation/eeqbc_charges.json`, `scripts/validation/oracle_eeqbc.py`) y gate `f22_eeqbc_matches_multicharge` en `cross_validation.rs`: cargas a 6.9e-12 e y CN a 1.9e-11 (tolerancia 1e-10). Fidelidad: la distancia vdW que entra al capacitor es el valor crudo de la tabla mctc-lib (los factores `aatoau`/`autoaa` se cancelan). Sin cablear a D4 todavía (dftd4 lo selecciona con `--qmodel eeqbc`; nuestro `D4Config` usa EEQ2019)
- [x] Population: warning de cobertura (elementos sin modelo de valencia / sin parámetros EHT) y homo/lumo/gap expuestos en `PopulationResult` (Rust y binding Python)
- [x] Fukui: Hirshfeld/CM5 evaluado (requieren densidades promoleculares/cargas CM5 ausentes; NPA-condensado requiere EHT con carga, hoy no expuesto) — se mantiene Mulliken con la caveat documentada
- [x] Surface: tabla de radios documentada (`VDW_FALLBACK_RADIUS`, `vdw_radius_source`) y evaluación Connolly/SES documentada (diferida hasta que un consumidor necesite la malla)
- [x] Optimization: unificar optimizadores entre EHT/materials (EHT usa el BFGS compartido de `materials::geometry_opt`) y repulsión nuclear incluida en el gradiente EHT
- [x] Variable-cell BFGS: `materials::variable_cell` con stress numérico por diferencias centrales de la energía bajo deformación afín (σ = (1/V)∂E/∂ε, convención tracción-positiva), relajación alternada posiciones (BFGS de `optimize_framework`, ahora genérico sobre el callback) + celda con line search que solo acepta pasos que bajan la energía. Gate autocontenido: modelo LJ en celda cúbica simple (mínimo de E(a) por búsqueda dorada) — stress ≈ 0 en el equilibrio (<1e-4 eV/Å³), signos correctos, y la relajación recupera a_eq desde 1.12·a_eq

### 15. Reaction Dynamics 3D (alpha)
- [x] Deprecar `dynamics.rs` legacy como assembly primario — `dynamics.rs` conserva solo energía/NEB; el assembly vive en `alpha::reaction_dynamics::complex` (documentado en la cabecera del módulo)
- [x] Rotación Kabsch para >2 fragmentos (`assemble_fragments_at_product_positions`)
- [x] Posicionamiento de moléculas extra (≥2) según átomos equivalentes en el producto
- [x] Optimizador global de orientación multi-fragmento (`optimize_reactive_complex_global`: Kabsch + rotaciones aleatorias + relajación, se guarda la mejor; tests de cota)
- [ ] Comparación sistemática de ángulos de ataque, barreras y trayectorias vs referencias externas (requiere set de reacciones con barreras de literatura)
- [x] Frontend (V-path X-axis fallback) — **cerrado como no aplicable en este repo**: el código afectado vive en el repositorio del frontend; aquí solo se exponen los payloads (render bridge, charts) que el frontend consume. Ruta: corregirlo en el repo del frontend y, si el payload necesita metadatos adicionales (p. ej. eje X explícito), exponerlos desde `alpha::render_bridge`

### 16. Alpha Export Gaps (bindings faltantes)
- [x] Render Bridge (WASM): las 10 funciones de `alpha::render_bridge` que ya existían en Rust (edl_profile_chart, capacitance_scan_chart, periodic_dos_chart, band_structure_chart, microkinetic_population_chart, arrhenius_chart, gsm_path_trajectory, pack_chart_payload, verify_chart_round_trip) están ahora expuestas como `alpha_render_*` en `crates/wasm/src/alpha.rs`, cada una gated con su feature de entrada; el crate compila con y sin las features. Pendiente: Python (`alpha_render_py` solo tiene edl/capacitance/arrhenius; faltan ~5 wrappers JSON — no verificables aquí por el bloqueo pyo3 0.22/Python 3.14) y 3 chart types nuevos (kpoint_path, trajectory, dos) que no existen en Rust
- [x] Periodic Linear (WASM): expuestas `alpha_periodic_monkhorst_pack`, `bloch_phase`, `bz_integrate_scalar`, `bz_integrate_vector`, `validate_kmesh_weights`, `validate_electron_count` y `band_structure_eht` (esta última cablea `eht::band_structure::compute_band_structure` con k-path). Pendiente: `kpath_standard` y `fermi_surface_from_bands` (no existen como funciones Rust; serían nuevas), los adaptadores KPM/RandNLA periódicos (firmas con matrices) y el mismo set en Python (bloqueado por pyo3 0.22/Python 3.14)
- [x] Kinetics CLI: los 3 comandos eran stubs con JSON hardcodeado (y el feature `alpha-kinetics` no compilaba: `limitations` sin `mut` en `alpha/gsm/backend.rs`). Ahora `extract-kinetics-diagnostics`, `solve-microkinetic-network` y `solve-microkinetic-steady-state` aceptan `--network <JSON|ruta>` (con red demo integrada) y llaman a las APIs reales. Hallazgo colateral corregido: el solver explícito de Euler era inestable en redes stiff (barrera 0.3 eV a 350 K → k ≈ 3.6e8 s⁻¹ → poblaciones no finitas con 0.01 s); ahora el sub-stepping se elige por la tasa máxima (`k·Δt ≲ 0.1`) y hay guarda de finitud, con test de regresión

### 17. Validación a Gran Escala (Track G)
- [x] Banco ≥100 moléculas (180 diversas + 1000 ChEMBL + 400 UFF; guard `track_g_bank_has_enough_molecules`)
- [x] Pipeline automatizado de benchmarking vs PySCF/RDKit (oráculos + `run_validation.sh` + `collect.py`)
- [x] Threshold <1 % documentado por observable (plan §14.2; los gated quedan 2–5 órdenes por debajo)
- [x] Límites de escalado documentados (plan §14.3: tiempos por método y memoria ERI O(N⁴))

### 18. Experimental Tracks Pendientes de Validación para Promoción
- [x] Validación cuantitativa por submódulo (beta): benchmarks nuevos en `tests/beta_benchmarks.rs` — KPM DOS vs diagonalización (DOS integrada 40.17/40, centro 1e-4, corr 0.86 con el kernel Jackson), Nyström vs SVD exacta (error ≤ 2× el mejor rank-k y decreciente), MBH vs Hessiano analítico proyectado (desviación máxima 2.6 %, modos casi nulos con ruido de diferencias finitas ~5 cm⁻¹). Pendiente: cpm/riemannian/sdr/gsm y E1–E11 siguen con evidencia parcial
- [x] Benchmarking frente a rutas estables: KPM vs diagonalización, Nyström vs SVD exacta y MBH vs Hessiano analítico implementados y verdes en `tests/beta_benchmarks.rs` (feature-gated). Pendiente: GSM vs NEB y SDR vs Butina (requieren habilitar `alpha-gsm`/`alpha-sdr` y definir la métrica de comparación de caminos/clusters)
- [x] Criterios claros de promoción experimental → producción (`src/experimental_status.rs`, actualizado con la evidencia actual)
- [x] Tests de integración para beta: `tests/beta_benchmarks.rs` cubre kpm/randnla/mbh contra rutas estables (no solo smoke tests) con gates cuantitativos documentados

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

### 17. LAMMPS (oráculos)
- [ ] F78/ReaxFF: **desbloqueado** — build con REAXFF funcional (`cmake ... -DCMAKE_DISABLE_FIND_PACKAGE_MPI=ON`; ver §14.0 del doc de validación). Siguiente paso: `scripts/validation/oracle_reaxff_lammps.py` con `ffield.reax.cho` + `pair_style reaxff` vs `compute_reaxff_energy`, fixture y gate en `alpha_reaxff_validation.rs`
- [ ] F95/CPM: LAMMPS de apt ya trae ELECTRODE; falta `scripts/validation/oracle_cpm_lammps.py` (`fix electrode/conp`)
