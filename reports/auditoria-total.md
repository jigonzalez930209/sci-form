# AUDITORÍA TOTAL — sci-form v0.15.2

**Fecha:** 2026-10-20 · **Alcance:** 100% del repositorio (357 archivos Rust en `src/`, 4 crates, 89 archivos de test, 63 docs, scripts, CI, packaging).
**Método:** lectura de código fuente por 7 subagentes de dominio + verificación mecánica (LOC, API pública, gates, cobertura) + ejecución real de la suite + verificación directa de afirmaciones documentales contra el código.

---

## 0. RESUMEN EJECUTIVO

### Números

| Métrica | Valor |
|---|---|
| Archivos Rust en `src/` (excl. `src/bin/`) | 330 |
| LOC en `src/` | 160 768 |
| LOC en `crates/` (python + wasm + cli) | 15 950 |
| LOC en `tests/` | 44 189 |
| Archivos de test | 89 |
| Tests unitarios (`cargo test --lib --all-features`) | **961 pass / 1 FAIL / 2 ignored** |
| Gates de validación cruzada (`cross_validation.rs`) | 38 tests, **F23/F64/F65/F66/F69 sin implementar** |
| Gates alpha/beta (F77–F95) | ~44 tests — **nunca ejecutados por CI** |
| Tests Python | 0 ejecutables (script sin pytest, no en CI) |
| Tests WASM/JS | 0 ejecutables (ruta de import rota, no en CI) |
| Funciones API que devuelven datos falsos (hard-coded) | ~26 en Python + ~26 en WASM + ~11 en CLI |
| Subcomandos CLI | 68 (7 handlers inalcanzables) |
| Claimed-vs-real en `docs/validation/cross-library-validation.es.md` | 5 fichas PASS sin test (F23/F64/F65/F66/F69) |
| Parches de CI que no ejecutan la science-gate suite | **44** |

### Veredicto por dominio

| Dominio | Estado | Gate externo real |
|---|---|---|
| SCF / HF-3c / CISD | **Sólido** — paridad PySCF 1e-9/1e-10 | F36, F37, F38, F39, F41 |
| PM3 | **Sólido** (ΔHf vs MOPAC 23.2.5 <0.02 kcal/mol) | F42 |
| GFN0/GFN1/GFN2 + UKS | **Sólido** (vs xtb 6.7.1 / tblite 0.7.0) | F43, F44, F45, F44-UKS, F45-UKS |
| D3 / D4 | **Sólido** (vs s-dftd3, dftd4 v4.3.0) | F46 |
| DFT (KS) | **Aceptable** — energía vs PySCF, pero `coulomb_energy=0.0` hard-coded | F77 |
| UFF / atom typing | **Aceptable** — energía agregada vs RDKit | F16, F17 |
| EHT | **Débil** — solo tabla de parámetros validada, no la energía | F31 (parcial) |
| MMFF94 | **FALSO** — parámetros placeholder, sin paridad posible | — |
| ETKDG / distance geometry | **Débil** — sin paridad vs RDKit `EmbedMolecule` | — |
| Gasteiger / SASA / dipole / Kabsch / SSSR / ECFP / Butina | **Aceptable** | F21, F28, F26, F13, F05, F59, F14 |
| ANI | **Roto** — CLI falla en runtime, AEV con alias, fuerzas ≠ gradiente | — |
| ESP / NPA / NMR shifts / descriptors ML | **Débil** — sin paridad externa | — |
| GB / ALPB | **Débil** — solo aserciones de signo | — |
| Materiales (XRD, CIF, space groups) | **Débil** — fixtures generados, ningún test los consume | — |
| MD / Langevin | **Roto** — inconsistencia de unidades kcal vs eV | F84 (fuerzas) |
| GPU | **Débil** — `best_available()` hace panic sin display; `mmff94_gpu` numéricamente inválido | — |
| Alpha/beta | **Aceptable** donde hay gate (F81–F95); KPM/RandNLA/MBH no cumplen criterios de promoción | F81–F95 |

### P0 — bloqueantes (rompen ciencia o builds)

| # | Hallazgo | Ubicación | Impacto |
|---|---|---|---|
| P0-1 | `GpuContext::best_available()` hace **panic** en `wgpu-hal` (`BadDisplay`) en vez de caer a CPU. Flaky: 961/1 fail en un run, 962/0 en otro. | `src/gpu/context.rs` | CI inestable; test `test_compute_esp_grid_water` falla |
| P0-2 | **MMFF94 usa parámetros placeholder** (`kb=5.0`, `k_a=0.5`, θ₀=109.5°, `r*=3.5/ε=0.05`). No hay tablas MMFF94 reales. | `src/forcefield/mmff94.rs` | `compute_mmff94_energy` NO devuelve energía MMFF94. Falsa confianza alta |
| P0-3 | **~26 funciones por lenguaje devuelven datos falsos**: `total_dos: vec![1.0;100]`, `homo_energy_ev:-5.0`, `report_generator → "Generated: 2026-03-28"`, `validate_experimental_result → always true` | `crates/python/src/experimental.rs`, `crates/wasm/src/experimental.rs`, `crates/cli/src/experimental_cmds.rs` | API pública que miente al usuario |
| P0-4 | `GpuContext` + `mmff94_gpu`: `atomicAdd` sobre `bitcast<u32>` de floats (suma entera de bit-patterns) | `src/gpu/mmff94_gpu.rs` | El path GPU no puede reproducir la energía CPU |
| P0-5 | `ANI` CLI usa `make_test_model(384)` cuando la longitud AEV real es 4032 → error en toda llamada | `crates/cli/src/calc_cmds.rs` | Subcomando `ani` 100% roto |
| P0-6 | `ani/aev.rs`: índice angular `lo*(2N-lo-1)/2 + (hi-lo)` **aliasa 6 de 28 bloques**; `ani/gradients.rs` backpropaga solo el bloque radial | `src/ani/aev.rs`, `src/ani/gradients.rs` | Las fuerzas ANI **no son** −dE/dR |
| P0-7 | Langevin BAOAB usa factor de aceleración en eV (`0.009648533`) y `k_B` en eV/K, pero las fuerzas están en kcal/(mol·Å) → error dimensional 10²–10⁵ | `src/dynamics_live/langevin.rs` | Termodinámica de Langevin cuantitativamente incorrecta |
| P0-8 | `docs/validation/cross-library-validation.es.md` §13.1 declara PASS con valores medidos para F23, F64, F65, F66, F69 — **no existe ningún `fn f23_*`/`f64_*`/`f65_*`/`f66_*`/`f69_*` en `tests/`** | docs + `tests/` | 5 validaciones ficticias |
| P0-9 | `pkg` y `pkg-node` son **symlinks absolutos** versionados a `/home/kuran/sci-form/...` | raíz del repo | Clone de cualquier otro usuario → symlinks colgantes; artefactos de build en git |
| P0-10 | `cargo clippy --all-targets` en CI solo corre para el paquete raíz → `crates/cli` y `crates/python` nunca se lintean | `.github/workflows/ci.yml` | Regresiones de lint en bindings |

### P1 — alto (validación ausente donde se afirma haberla)

| # | Hallazgo |
|---|---|
| P1-1 | **Sin paridad ETKDG vs RDKit.** F70 usa la geometría ETKDGv3 de RDKit como *entrada*, nunca como *referencia*. Sin fixture de bounds, RMSD ni tasa de rechazo vs `AllChem.EmbedMolecule`. |
| P1-2 | **Sin paridad de la matriz de bounds.** `tests/regression/test_bounds_compare.rs` está `#[ignore]`d y requiere un `/tmp/rdkit_bounds_*.json` manual. Es el control de mayor valor de todo el dominio geométrico. |
| P1-3 | **Sin paridad de force field ETKDG-3D** (`construct3DForceField`): ni conteo de constraints, ni longitudes min/max, ni energías por término. |
| P1-4 | **Sin paridad de BFGS** (`minimize_bfgs_rdkit` afirma port exacto de `BFGSOpt.cpp`: FUNCTOL/MOVETOL/EPS, line search cúbico, gradient-scaling) — nunca comparado. |
| P1-5 | **SMARTS sin paridad semántica**: solo "parsea sin error". Sin comparación de `GetSubstructMatches` (ni conteo ni mapeo) → una divergencia silenciosa corrompe las preferencias de torsión CSD. |
| P1-6 | **Tabla CSD de torsiones sin hash** contra `torsionPreferences_v2.in` de RDKit → deriva indetectable. |
| P1-7 | **Sin paridad de gradients analíticos** en ningún método cuántico (SCF, PM3, xTB, DFT, EHT). Solo se verifica la 3ª ley de Newton. |
| P1-8 | **HF-3c sin paridad end-to-end**: `gcp.rs` y `srb.rs` usan parameterizaciones aproximadas (proxy de `n_virtual`; damping SRB publicado no implementado) sin referencia. |
| P1-9 | **EHT sin paridad de energía/overlap/eigenvalores** — solo la tabla de parámetros (F31). |
| P1-10 | **ESP sin paridad** (APBS/pdb2pqr/Multiwfn/Psi4) y solo point-charge; el argumento `elements` se ignora. |
| P1-11 | **NPA/NBO sin paridad** (Psi4/GenNBO/Multiwfn) — 3 tests con matrices 2×2 de juguete. |
| P1-12 | **Población EHT sin paridad**: F36 valida `scf/mulliken.rs`, no `population/population.rs`. Wiberg/Mayer no comparados contra el estándar `ΣPμνSμν`. |
| P1-13 | **Estereoquímica sin referencia CIP**: sin `AssignStereochemistry` de RDKit; la convención de signo R/S parece **invertida** vs OpenSMILES; sin reglas CIP 4/5; ignora chiral tags del SMILES de entrada. |
| P1-14 | **GB/ALPB sin paridad** contra xTB/ALPB, AmberTools GB, OpenMM `CustomGBForce`; la integral de descreening no es la forma HCT publicada. |
| P1-15 | **XRD sin paridad**: fixture `pymatgen_xrd.json` existe, ningún test la carga. Sin extinciones sistemáticas ni multiplicidad. |
| P1-16 | **230 space groups sin paridad**: fixture `space_groups_230.json` (412 KB) existe, ningún test la carga. Solo se verifica existencia. |
| P1-17 | **CIF sin paridad** contra pymatgen/ASE (`pymatgen_cif.json` sin usar); sin corpus COD; elementos >Bi → Z=0. |
| P1-18 | **Dinámica de reacción estable sin ningún test**: `compute_reaction_dynamics` (expuesto como API pública) + `map_atoms_greedy`, `reorder_coords`, `slide_molecules`, `com_flat`, `centre_at_origin` → 0 tests. |
| P1-19 | **CI-NEB / IRC / complex-assembly sin paridad** vs ASE NEB (criterio ±1 kcal/mol y ≤0.05 Å documentado, no implementado). |
| P1-20 | **Variable-cell sin validación externa**: solo modelo LJ interno; `numerical_stress` tiene ramas diagonal/off-diagonal idénticas (falta el factor ½ documentado); solo strain isotrópico. |
| P1-21 | **Descriptors ML sin paridad**: F56/F57/F58 planificados, `scripts/compare_descriptors_rdkit.py` no existe. `balaban_j` es un placeholder `|ln µ|+1`. |
| P1-22 | **NMR sin paridad** contra nmrdb2/NIST SDBS: solo spot-check de 3 moléculas ±1 ppm. `nmr/hose.rs` es código muerto (`predict_shift_from_hose` nunca se llama). |
| P1-23 | **"GIAO" NMR no es GIAO**: sin base gauge-including, sin respuesta RPA/RFIA, anisotropía hard-coded 0.0. |
| P1-24 | **sTDA UV-Vis sin paridad** contra PySCF `tdscf.TDA`; `config.ax` declarado y nunca usado (no hay término de intercambio). |
| P1-25 | **IR intensities sin paridad** NIST/CCCBDB; `spectroscopy/ir_intensities.rs` tiene 0 tests y campos de config sin usar. |
| P1-26 | **EEQ2019 completamente sin test** (0 unit, 0 integración) y ya usado en producción por el D4. |
| P1-27 | **EEQ energy ≠ energía de los charges de EEQ** (σ distinto en `energy.rs` vs `charges.rs`); γ en 1/Å mezclado con η en eV. |
| P1-28 | **Descriptores ANI sin paridad** contra TorchANI/ASE; sin pesos ANI-2x oficiales y formato de pesos bespoke. |
| P1-29 | **Metal EHT sin paridad de orbitales** contra yaehmop; solo la tabla de parámetros. |
| P1-30 | **KPM/RandNLA/MBH no cumplen criterios de promoción** (>1% error en rank N/4, sin speedup N>500, sin ±2 cm⁻¹ experimental). |

### P2 — medio (calidad, consistencias, docs)

| # | Hallazgo |
|---|---|
| P2-1 | CI ejecuta solo `cargo test --lib` + `cargo test --test ci`. **Los 44 gates alpha/beta y los 38 de `cross_validation.rs` nunca corren en CI.** |
| P2-2 | `RUST_EXPERIMENTAL_FEATURES` omite `alpha-dynamics-live` → los 15 tests de `dynamics_live` nunca se ejecutan. |
| P2-3 | `tests/experimental/test_{cga,sdr,riemannian}.rs` están gated por features **legacy** (`experimental-*`) → compilan a módulos vacíos en CI. |
| P2-4 | `docs/api/python.md` documenta 3 funciones inexistentes (`charges_configured`, `ir_spectrum_broadened`, `ensemble_j_couplings`). |
| P2-5 | `docs/api/typescript.md` documenta 2 exports WASM inexistentes. |
| P2-6 | `docs/guide/python.md` usa nombres de función **incorrectos en todos los bloques de código** (`sci_form.compute_charges` etc. → `AttributeError`). |
| P2-7 | `docs/api/cli.md` documenta 22 de 68 subcomandos. `docs/guide/cli.md` cubre 12. |
| P2-8 | `docs/index.md` dice "(v0.5.0)" vs crate 0.15.2; `docs/guide/{rust,typescript,python}.md` imprimen "0.9.1". |
| P2-9 | `docs/algorithm/webgpu-validation.md` documenta un harness (`crates/wasm/examples/vite-webgpu-hybrid/`) **inexistente en git**. |
| P2-10 | 12 páginas docs son **orphans** (fuera del nav de VitePress), incluida `complete-supported-reference.md`. |
| P2-11 | `src/nmr/shifts.rs:1125` — `test_benzene_h_shifts` **le falta `#[test]`**, nunca se ejecuta. |
| P2-12 | `src/materials/mod.rs` declara un **segundo tipo `SpaceGroup`** (cell.rs legacy vs space_groups.rs) y no re-exporta `geometry_opt`/`space_groups`/`variable_cell`/`topology_data`. |
| P2-13 | `crates/wasm/js/alpha/index.d.ts` **deriva de las firmas Rust** (`apply_thermostat`, `apply_steering_force`, `set_temperature`, `num_atoms` ausente). |
| P2-14 | `crates/wasm/src/{alpha,beta}.rs` reportan versiones `0.14.4` / `0.11.2` vs crate 0.15.2. |
| P2-15 | `crates/cli/src/format.rs::element_symbol` devuelve `"X"` fuera de 13 elementos → SDF inválido para cualquier metal. |
| P2-16 | `crates/cli/src/embed_cmds.rs::cmd_info` anuncia "IR, NMR, UV-Vis" sin subcomandos correspondientes. |
| P2-17 | `crates/cli/src/calc_cmds.rs::cmd_esp` descarta los valores del grid (solo metadatos). |
| P2-18 | `crates/cli/src/calc_cmds.rs::cmd_hf3c` usa `n_cis_states: 0` vs default 5 en Python. |
| P2-19 | `crates/cli/src/calc_cmds.rs::cmd_ani` usa red de test (384) vs default de librería (4032) en Python → energías no comparables. |
| P2-20 | `crates/python/Cargo.toml` pinnea `pyo3 = "0.22"` (máx Python 3.13) mientras `pyproject.toml` declara `requires-python = ">=3.8"`; en Python 3.14 requiere `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1`. |
| P2-21 | `crates/wasm/Cargo.toml` usa `wgpu 29` vs `wgpu 0.19` del core → skew de versión. |
| P2-22 | `scripts/post-wasm-build.mjs` parchea `pkg-node/package.json` a `sci_form.js` pero wasm-pack emite `sci_form_wasm.js` → `npm run wasm:build:node` produce un pkg inservible. |
| P2-23 | `package.json` declara `"main": "index.js"` sin `index.js` y **no tiene script `test`**. |
| P2-24 | `tests/integration/test_wasm_integration.js` resuelve `../crates/wasm/pkg/...` desde `tests/` → ruta inexistente → `MODULE_NOT_FOUND`. |
| P2-25 | `src/experimental_status.rs` referencia paths `experimental/{kpm,mbh,...}` **borrados** y features legacy; omite `alpha/{reaction_dynamics,periodic_linear,render_bridge,kinetics}` y todo `src/gpu/`. |
| P2-26 | `src/gpu/shader_registry.rs` no registra 4 shaders existentes (ANI AEV, HCT Born, MMFF94 NB, sTDA matmul) → drift registry/kernel. |
| P2-27 | `src/distgeom/parallel.rs::compute_bounds_parallel` calcula bounds **completamente distintos** de `calculate_bounds_matrix_opts` (API divergente + GPU stub). |
| P2-28 | Dos implementaciones de SSSR (`distgeom/validation.rs::find_sssr` vs `rings/sssr.rs::compute_sssr`) y dos superficies de energía FF (`forcefield/energy.rs` vs `forcefield/uff.rs`) sin cross-check. |
| P2-29 | `forcefield/etkdg_lite.rs` reproduce deliberadamente un **bug de RDKit** (6º armónico usa V[4]/s[4]) sin fixture que lo demuestre. |
| P2-30 | `src/sc f/screening.rs::SchwarzScreener` y `src/scf/validation.rs` son **código muerto** (0 callers en src/ y tests/). |
| P2-31 | `src/beta/riemannian.rs` + `src/alpha/{cga,edl,sdr}.rs` son shims de back-compat sin callers. |
| P2-32 | `src/alpha/kinetics/gpu.rs`, `src/alpha/periodic_linear/gpu.rs`, `src/alpha/render_bridge/gpu.rs`, `src/edl/gpu.rs` son **stubs CPU** con tests "GPU/CPU parity" que comparan CPU contra CPU. |
| P2-33 | `src/materials/topology_data.rs`: `edges.len() != edge_types.len()` en 4 de 6 redes; `edges` nunca se consume. |
| P2-34 | `src/materials/variable_cell.rs::strained_lattice` contiene `let eps = …; let _ = eps;` muerto. |
| P2-35 | `src/spectroscopy/stda_uvvis.rs` `config.ax` sin usar; `src/spectroscopy/ir_intensities.rs` `broadening`/`freq_min`/`freq_max` sin usar. |
| P2-36 | `src/ml/symmetry_functions.rs`: el bucle radial hace `zip` de shifts×etas → solo 16 de 256 slots se llenan; el oráculo F79 reproduce la misma truncación. |
| P2-37 | `src/ml/advanced_models.rs`: `oob_score` siempre `None`; subsampling de features sembrado por camino de recursión, no por nodo. |
| P2-38 | `src/charges_eeq/energy.rs` × 23.06 eV→kcal/mol sin justificación (γ en 1/Å). |
| P2-39 | `src/riemannian/manifold.rs`: `tangent_projection` es solo simetrización (no es proyección tangente real); `psd_distance` es proxy Frobenius, no la métrica afín-invariantemente invariante anunciada. |
| P2-40 | `src/beta/riemannian.rs` `config.fallback_enabled` / `used_fallback` muertos. |
| P2-41 | `src/dft/ks_fock.rs`: `coulomb_energy` hard-coded a 0.0 ("approximate"). `build_vxc_uks` clona el resultado closed-shell para PBE. |
| P2-42 | `src/dft/lebedev.rs`: las grids de "110" y "302" puntos son **aproximaciones** (latitud/icosaedro, espiral de Fibonacci), no Lebedev — y `GridQuality::Fine/VeryFine` dependen de ellas. |
| P2-43 | `src/dft/functionals/{svwn,pbe}.rs`: ∂/∂ρ y ∂/∂σ por diferencias numéricas, sin check analítico-vs-numérico. |
| P2-44 | `src/dft/ecp.rs`: ECP "no fitted Stuttgart/LANL" (lo dice el docstring) → energías Z>30 sin validar. |
| P2-45 | `src/riemannian`, `src/beta/mbh/hessian.rs`: factor de conversión 108.59 hard-coded. |
| P2-46 | `src/xtb/gfn2.rs` tolerancia F45 (0.0217 eV) más laxa que F44 (0.01 eV). |
| P2-47 | `src/forcefield/builder.rs`: torsiones divididas por número de caminos, sin escalado 1-4 en LJ → paridad solo agregado. |
| P2-48 | `src/esp/esp.rs`: `.cube` escribe cargas, no la densidad ESP verdadera. |
| P2-49 | `src/transport/arrow.rs` no tiene dependencia `arrow` ni IPC → "Arrow-compatible" es overstated. |
| P2-50 | `src/dynamics.rs::atomic_mass_amu` devuelve 12.0 amu para ~100 elementos. |
| P2-51 | `src/dynamics_live/neighbor_list.rs` es O(N²) pese al doc "O(N)". |
| P2-52 | `src/dynamics_live/pbc.rs` + `materials/cell.rs::minimum_image_distance`: `round()` fraccional no es minimum-image real para celdas triclínicas. |
| P2-53 | `src/dynamics_live/barostat.rs`: el doc dice `N·k_B·T + virial`, el código usa `(2·KE + virial)/(3V)`; y no existe virial en el repo. |
| P2-54 | `src/pm3/solver.rs::test_pm3_h2` está `#[ignore]`d ("SCF rompe la simetría de inversión del H2"). |
| P2-55 | `src/smarts/torsion_data.rs` 3247 LOC de datos generados sin hash de verificación. |
| P2-56 | `src/hf/d3_data.rs` (24 458 LOC) y `src/xtb/d4_data.rs` (4676 LOC) sin test que compare contra el archivo upstream byte a byte. |
| P2-57 | `src/materials/space_group_table.rs` (6751 LOC) sin script generador versionado → deriva silenciosa. |
| P2-58 | `src/xtb/uks.rs`: el término de spin no tiene equivalente en tblite 0.7.0 → validado solo por signo y conservación. |
| P2-59 | `src/beta/kpm/density.rs`: KPM Mulliken usa simplificación "1 BF por núcleo". |
| P2-60 | `src/alpha/reaction_dynamics/per_frame.rs`: "Wiberg" es aproximación P² con base→átomo heurística. `orbital_approach.rs` usa tabla `basis_functions_for_element` en vez del basis EHT real. |

---

## 1. TABLA MAESTRA — `src/` (330 archivos)

Leyenda Testing: `UT n` = n tests unitarios; `INT` = tests de integración; `—` = sin cobertura.
Leyenda Externo: `✔` = existe gate de paridad; `✘` = NO existe gate de paridad.

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/lib.rs` | 3004 | Raíz del crate: 53 `pub mod` + 17 tipos públicos + 87 funciones públicas. Declara los tipos de contrato (`ConformerResult`, `ConformerConfig`, `ConformerSearchResult`, `MethodCapability`, `ScientificMethod`, `PropertyRequest`, `SystemMethodPlan`, `ElectronicWorkflowResult`, `MethodComparisonResult`, `AromaticityAnalysis`, `StereocenterAnalysis`, `GraphFeatureAnalysis`, `GiaoNmrConfig/Result`) y los entry points de alto nivel (`embed`, `embed_batch`, `embed_diverse`, `parse`, `compute_*` × 25, `compare_methods`, `get_system_capabilities`, `get_system_method_plan`). Gatea `alpha`/`beta`/`dft`/`mlff`/`dynamics_live` por feature | 171 items públicos (87 fn, 18 struct, 5 enum, 58 mod) | UT 23 ✔ · INT F01–F70 vía `lib::` — PASS | — es la fachada; su superficie es la que debe replicar bindings |

### 1.1 Geometría, conformeros, force fields

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/graph.rs` | 566 | Grafo molecular sobre `petgraph`; radios covalentes/VdW (derivados RDKit); BFS `min_path_excluding`; `get_corrected_ideal_angle` portado de RDKit (3/4-anillos 60°/90°, exocíclicos 116°/112°) | `Molecule`,`Atom`,`Bond`,`Hybridization`,`ChiralType`,`BondStereo`,`BondOrder`,`get_covalent_radius`,`get_vdw_radius`,`get_ideal_angle`,`get_corrected_ideal_angle`,`min_path_excluding` | UT 0 · INT F01, tet_centers, algorithm_validation — PASS | RDKit `ForceFieldHelpers` ✔ (F01 parcial) · **✘ sin paridad de la tabla de ángulos de anillo** |
| `src/smiles.rs` | 820 | Parser SMILES descendente recursivo: ramas, anillos, brackets/isótopos/cargas, `@`/`@@`, `/` `\`, H implícitos, hibridación | `SmilesParser::{new,parse}` | UT 6 ✔ · INT F01, algorithm_validation, smirks — PASS | RDKit `SmilesParse` ✔ (F01 parcial) · **✘ sin paridad de stereo/isótopos/H explícito** |
| `src/canonical_smiles.rs` | 384 | Ranking canónico Morgan (hash FNV-1a) + DFS iterativo para dígitos de anillo | `to_canonical_smiles` | UT 5 ✔ | RDKit `MolToSmiles` — **✘ NO HAY PARITY TEST** |
| `src/topology.rs` | 288 | Clasificación de geometría de coordinación de metales por cutoff de radio + brute-force vs poliedros ideales | `analyze_topology`,`CoordinationGeometryGuess`,`MetalCoordinationCenter` | UT 3 ✔ · INT ninguna | Open Babel / ChemAxon — **✘ sin paridad** |
| `src/distgeom/mod.rs` | 24 | Fachada de re-exports | — | — | — |
| `src/distgeom/bounds/mod.rs` | 10 | Fachada bounds | `get_bond_length`,`calculate_bounds_matrix`,`smooth_bounds_matrix`,`triangle_smooth*` | — | RDKit `DistGeom` |
| `src/distgeom/bounds/geometry.rs` | 442 | Longitudes 1-2 UFF (Rappé eq. 2-4) con radios conjugados `_R`; detección cis/trans 1-4/1-5 amide/ester portada de RDKit | `get_bond_length` (+ `pub(super)` resto) | UT 0 · INT bounds_compare (`#[ignore]`d), debug | RDKit `set12/14/15Bounds` — **✘ sin paridad en CI** |
| `src/distgeom/bounds/matrix.rs` | 459 | Matriz de bounds completa: FW topológico, 1-2/1-3/1-4/1-5 con semántica de unión RDKit, regla `squishAtoms` 5-anillos | `calculate_bounds_matrix`,`calculate_bounds_matrix_opts` | UT 0 · INT (ignorada) | RDKit `getDistanceBounds` — **✘ EL GAP #1 DEL DOMINIO** |
| `src/distgeom/bounds/smoothing.rs` | 101 | Suavizado por desigualdad triangular (upper/lower) con clamp de tolerancia | `triangle_smooth`,`triangle_smooth_tol`,`smooth_bounds_matrix` | UT 2 ✔ | RDKit `smoothBounds` — **✘ sin paridad** |
| `src/distgeom/chirality.rs` | 246 | Conjuntos quirales: orden canónico de vecinos, paridad por permutación, volumen firmado (±2..±50) | `calc_chiral_volume`,`identify_chiral_sets` | UT 3 ✔ | RDKit `findChiralSets` — **✘ sin paridad** |
| `src/distgeom/embedding/mod.rs` | 8 | Fachada | — | — | — |
| `src/distgeom/embedding/coordinates.rs` | 604 | Núcleo DG: muestreo en orden triangular RDKit, matriz métrica Cayley-Menger, `power_eigen_solver` RDKit-exacto | `pick_rdkit_distances`,`compute_metric_matrix`,`generate_3d_coordinates`,`compute_initial_coords_rdkit`,`power_eigen_solver` | UT 0 · INT embedding_trace, differential_testing | RDKit `computeInitialCoords` — **✘ sin paridad de métrica/eigen** |
| `src/distgeom/embedding/rng.rs` | 107 | `MinstdRand` (48271/2³¹−1) y `Mt19937` (624 palabras) compatibles con boost | `DistGeomRng`,`MinstdRand`,`Mt19937` | UT 0 | boost/RDKit RNG — **✘ sin vectores de referencia** |
| `src/distgeom/parallel.rs` | 294 | Bounds paralelos (radios crudos, **distintos** de los canónicos) + FW CPU/paralelo + stub GPU | `compute_bounds_parallel`,`smooth_bounds_gpu`,`floyd_warshall_*` | UT 2 ✔ | — **✘ API divergente sin paridad** |
| `src/distgeom/validation.rs` | 661 | Validación parity-RDKit: SSSR, centros tetraédricos, volumen normalizado, signo de volumen quiral, planaridad SP2, linealidad de dobles | `identify_tetrahedral_centers`,`check_tetrahedral_centers`,`check_chiral_centers`,`check_planarity`,`check_double_bond_geometry` | UT 2 (solo `perturb_if_planar`) | RDKit `checkTetrahedralCenters` — **✘ 6 predicados核心 sin test** |
| `src/conformer.rs` | 768 | Pipeline ETKDG end-to-end: bounds→suavizado, quiral/tetra, 3D/4D, BFGS stage1/2, FF ETKDG-3D, fallback | `generate_3d_conformer`,`generate_3d_conformer_best_of_k`,`compute_total_bounds_energy_f64` | UT 1 · INT new_pipeline, geometry_quality, gdb20, ensemble_rmsd, ci — PASS | RDKit `EmbedMolecule`/ETKDGv3 — **✘ EL GAP #1 DEL DOMINIO** |
| `src/etkdg.rs` | 111 | Dihedral por atan2 + penalización de torsión armónica | `calculate_dihedral_angle`,`calculate_torsion_penalty` | UT 3 ✔ | RDKit `computeSignedDihedralAngle` — **✘ sin paridad de signo** |
| `src/alignment/mod.rs` | 5 | Fachada | `align_coordinates`,`compute_rmsd` | — | — |
| `src/alignment/kabsch.rs` | 431 | Kabsch: centroides, covarianza 3×3, SVD, rotación corregida por reflexión, RMSD; variante cuaternión | `AlignmentResult`,`compute_rmsd`,`align_coordinates`,`align_quaternion` | UT 10 ✔ · INT F13 — PASS | RDKit/NumPy **✔ F13** · `align_quaternion` sin fixture |
| `src/forcefield/mod.rs` | 42 | Fachada del árbol FF | re-exports | — | — |
| `src/forcefield/builder.rs` | 370 | Constructor UFF espejo de `UFF/Builder.cpp`: re-percepción aromática, atom typing, ángulos con hacks 3/4-anillos, torsiones, 3 inversiones por SP2 | `build_uff_force_field` | UT 0 · INT F16, F70, F84 — PASS | RDKit UFF builder **✔ F16/F70/F84 (agregado)** · **✘ sin descomp. por término** |
| `src/forcefield/atom_typer.rs` | 248 | Atom typing UFF (elemento+dígito hibridación+`_R`) + modelo de conjugación estilo RDKit | `assign_uff_type`,`is_atom_aromatic`,`conjugated_bond_flags` | UT 0 · INT F17 — PASS | RDKit `AtomTyper.cpp` **✔ F17 (labels)** · `conjugated_bond_flags` sin comparar |
| `src/forcefield/params.rs` | 345 | Tablas UFF (Rappé 1992) + fórmulas de longitud/k, k angular, `calcTorsionParams` (caso grupo-6, eq. 17), coeficientes de inversión | `UffAtomParams`,`get_uff_params`,`get_uff_bond_length`,`get_uff_angle_force_constant` | UT 0 (indirecto F16/F70) | RDKit `Params.cpp` — **✘ sin paridad de valores individuales** |
| `src/forcefield/uff.rs` | 453 | Términos UFF: stretch armónico, bend Fourier orden 0-4, torsión `0.5V(1−cosφ₀cosnφ)`, inversión UFF, LJ 12-6 | `UffHarmonicBondStretch`,`UffAngleBend`,`UffTorsion`,`UffInversion`,`UffLennardJones` | UT 0 · INT F16/F70/F84 — PASS | RDKit `ForceField.cpp` — **✘ sin breakdown ni FD-check angular/torsión** |
| `src/forcefield/mmff94.rs` | 1228 | Atom typing MMFF94 (75 tipos) + términos, **pero `build` usa parámetros placeholder** | `Mmff94AtomType`,`Mmff94Builder::{build,total_energy}`,`Mmff94Variant` | UT 8 ✔ | RDKit MMFF94 — **✘ IMPOSIBLE: no hay tablas reales (P0-2)** |
| `src/forcefield/traits.rs` | 116 | Trait `ForceFieldContribution` + registro de términos + breakdown agregado + energía/gradiente paralelo | `ForceFieldContribution`,`MolecularForceField`,`aggregate_breakdown` | UT 0 | RDKit `GetEnergyTerms` — **✘ sin paridad de breakdown** |
| `src/forcefield/energy.rs` | 481 | Modelo f32 legacy de energía total (bond/angle/bounds/OOP/torsión/vdW) — **superficie paralela** | `calculate_total_energy`,`vdw_energy`,`bounds_energy`,`FFParams` | UT 0 | — **✘ superficie duplicada sin validar** |
| `src/forcefield/gradients.rs` | 554 | Gradientes analíticos f32 espejo de `energy.rs` | `compute_analytical_gradient`,`analytical_grad_*` | UT 0 · INT gradient_check, gradient_perterm, grad_torsion_detail — PASS | RDKit `calcGrad` — **✘ FD-check solo en regression/, no en CV** |
| `src/forcefield/minimizer.rs` | 210 | L-BFGS (m=10) con cap de gradiente 10.0 + Armijo; **RMSD Kabsch duplicado** | `minimize_energy_lbfgs`,`calculate_rmsd_kabsch` | UT 0 | RDKit `ForceField::minimize` — **✘ sin paridad** |
| `src/forcefield/dg_terms.rs` | 154 | Implementaciones `ForceFieldContribution` para `DistanceViolation`/`ChiralViolation` (no usadas por el pipeline vivo) | `DistanceViolation`,`ChiralViolation` | UT 0 | RDKit `DistViolationContribs` — **✘ sin paridad,疑似 código muerto** |
| `src/forcefield/etkdg_lite.rs` | 290 | Torsión M6 Fourier con Chebyshev; **reproduce el bug copy-paste de RDKit** (6º armónico usa V[4]/s[4]) | `calc_torsion_energy_m6`,`calc_torsion_grad_m6`,`infer_etkdg_parameters` | UT 0 | RDKit `TorsionAngleM6` — **✘ sin fixture que demuestre el bug reproducido** |
| `src/forcefield/torsion_scan.rs` | 678 | Detección de enlaces rotables, rotación rígida, optimización de torsiones (greedy/bounds/MC/annealing/systematic) | `find_rotatable_bonds`,`optimize_torsions_*`,`torsional_sampling_diverse` | UT 1 ✔ · INT algorithm_validation, test_cga | RDKit `EmbedMultipleConfs` — **✘ sin paridad de ensemble** |
| `src/forcefield/bounds_ff/mod.rs` | 27 | Fachada + `ChiralSet`,`EmbedTorsion` | `ChiralSet`,`EmbedTorsion` | — | — |
| `src/forcefield/bounds_ff/bfgs.rs` | 383 | `BFGSOpt::minimize` RDKit-exacto (Hessiana N×N inversa, line search cúbico, hack de escalado de gradiente, TOLX/FUNCTOL/EPS) | `minimize_bfgs_rdkit` | UT 0 · INT step_by_step, energy_comparison, debug | RDKit `BFGSOpt.cpp` — **✘ sin paridad de trayectoria** |
| `src/forcefield/bounds_ff/lbfgs.rs` | 589 | Variantes L-BFGS sobre bounds FF con basin threshold, 4D y quiral | `minimize_bounds_lbfgs*`,`minimize_embedding_lbfgs` | UT 0 · INT bounds_compare | RDKit `ForceField::minimize` — **✘ sin paridad** |
| `src/forcefield/bounds_ff/energy.rs` | 454 | Energía/gradiente de violación de bounds + volumen quiral (f32/f64) + recolección de torsiones 4D | `bounds_violation_energy*`,`chiral_violation_*` | UT 0 · INT gradient_check | RDKit `DistViolationContribs` — **✘ pre-factor f32 difiere del de RDKit, sin cuantificar** |
| `src/forcefield/etkdg_3d/mod.rs` | 89 | Fachada + tipos `DistConstraint`,`M6TorsionContrib`,`UFFInversionContrib`,`AngleConstraint`,`Etkdg3DFF` | los anteriores | — | — |
| `src/forcefield/etkdg_3d/builder.rs` | 741 | Constructor FF ETKDG-3D parity-RDKit: marking 1-4, torsiones planas 4-6 anillos, heurísticas cadena (ester/amide V1=100, Ar–C(=O) V2=8…), 3 permutaciones de inversión, constraints 1-2/1-3/largo | `build_etkdg_3d_ff`,`build_etkdg_3d_ff_with_torsions` | UT 0 · INT debug_ff_energy, regression ff | RDKit `construct3DForceField` — **✘ sin paridad término a término** |
| `src/forcefield/etkdg_3d/energy.rs` | 773 | Energía FF ETKDG-3D f32/f64: flat-bottom, angle constraints, inversión UFF, OOP, M6, + `planarity_check_energy` | `etkdg_3d_energy*`,`planarity_check_energy*` | UT 0 · INT debug_ff_energy, energy_comparison | RDKit `calcEnergy` — **✘ sin paridad del orden de acumulación** |
| `src/forcefield/etkdg_3d/gradient.rs` | 629 | Gradientes FF ETKDG-3D f32/f64 | `etkdg_3d_gradient`,`etkdg_3d_gradient_f64` | UT 0 · INT gradient_check, gradient_perterm | RDKit `calcGrad` — **✘ sin FD-check del path f64** |
| `src/forcefield/etkdg_3d/optimizer.rs` | 661 | Minimizadores L-BFGS/BFGS del FF ETKDG-3D + `correct_torsion_angles` | `minimize_etkdg_3d*`,`correct_torsion_angles` | UT 0 · INT geometry_quality, step_by_step | RDKit `ForceField::minimize` — **✘ sin paridad de geometría final** |
| `src/forcefield/reaxff/mod.rs` | 30 | Fachada ReaxFF | re-exports | — | — |
| `src/forcefield/reaxff/params.rs` | 260 | Parámetros ReaxFF + parser `ffield.reax.cho` (De, p_be1/2, p_bo1..6, γ, r_sigma) | `ReaxffParams`,`get_bond_de`,`get_equilibrium_angle` | UT 5 ✔ | LAMMPS `reaxff` — **✔ parcial F78** · `get_equilibrium_angle`/`get_angle_force_constant` son stubs hard-coded sin paridad |
| `src/forcefield/reaxff/bond_order.rs` | 305 | Órdenes de enlace continuos BO_σ+BO_π+BO_ππ con taper Hermite | `BondOrder`,`compute_bond_orders` | UT 3 ✔ · INT F78 — PASS | LAMMPS `reaxff` **✔ F78** (solo sistemas H₂) |
| `src/forcefield/reaxff/eem.rs` | 162 | EEM (electronegativity equalization) para ReaxFF | `EemAtomParams`,`solve_eem` | UT 4 ✔ | LAMMPS `reaxff` — **✘ sin paridad de cargas EEM** |
| `src/forcefield/reaxff/energy.rs` | 210 | Términos bonded ReaxFF vs orden continuo | `ReaxffEnergy`,`compute_bonded_energy` | UT 2 ✔ | LAMMPS `reaxff` — **✔ parcial F78**, sin paridad por componente |
| `src/forcefield/reaxff/gradients.rs` | 174 | Gradientes por diferencias centrales de la energía bonded | `compute_reaxff_gradient` | UT 2 ✔ | LAMMPS `reaxff` (analítico) — **✘ sin paridad vs gradiente analítico** |
| `src/forcefield/reaxff/nonbonded.rs` | 188 | vdW + Coulomb blindado con taper | `compute_nonbonded_energy*` | UT 3 ✔ | LAMMPS `reaxff` — **✔ parcial F78** |
| `src/forcefield/reaxff/taper.rs` | 72 | Polinomio de taper Hermite + derivada | `taper_function`,`taper_derivative` | UT 2 ✔ | LAMMPS `reaxff` — **✘ sin paridad de coeficientes** |
| `src/optimization/mod.rs` | 8 | Fachada | `RustBfgsEngine` | — | — |
| `src/optimization/bfgs.rs` | 146 | BFGS denso + Armijo + hack de escalado RDKit + clamp 0.3 | `RustBfgsEngine::execute_minimization` | UT 0 · **0 referencias en todo el repo** | RDKit `BFGSOpt` — **✘ sin paridad** |
| `src/smarts/mod.rs` | 21 | Fachada | `parse_smarts`,`substruct_match*`,`match_experimental_torsions` | — | — |
| `src/smarts/parser.rs` | 759 | Parser SMARTS: primitivas de átomo (`#n`,`Hn`,`Xn`,`Dn`,`R`,`rn`,`r{n-m}`,`^n`, quiralidad) y de enlace, recursivo, AND/OR/NOT, anillos, atom maps | `SmartsPattern`,`SmartsAtom`,`SmartsBond`,`parse_smarts` | UT 7 ✔ | RDKit `MolFromSmarts` — **✘ solo "parsea", sin paridad semántica** |
| `src/smarts/matcher.rs` | 668 | Matching por backtracking ordenado por conectividad, queries con anillos, SMARTS recursivo, batch paralelo | `substruct_match`,`substruct_match_batch`,`has_substruct_match` | UT 5 ✔ · INT new_features, smirks | RDKit `GetSubstructMatches` — **✘ sin paridad de conteos/mapeos** |
| `src/smarts/torsion_data.rs` | 3247 | Tablas CSD/RDKit `torsionPreferences_v2.in` (837 patrones, 12 coeficientes Fourier c/u) | `TORSION_PREFS_V2`,`TORSION_PREFS_SMALL_RINGS`,`TORSION_PREFS_MACROCYCLES` | UT 0 (datos) | RDKit `torsionPreferences_v2.in` — **✘ sin hash de contenido** |
| `src/smarts/torsion_matcher.rs` | 230 | Match de patrones CSD → `M6TorsionContrib` :1-:4, exclusión de anillos, caché thread-local | `match_experimental_torsions` | UT 1 ✔ · INT ensemble_rmsd, debug | RDKit `getExperimentalTorsions` — **✘ sin paridad del set matcheado** |
| `src/rings/mod.rs` | 14 | Fachada | `compute_sssr`,`perceive_aromaticity`,`compute_ecfp`,`compute_tanimoto` | — | — |
| `src/rings/sssr.rs` | 569 | SSSR estilo Horton (ciclo más corto por arista + ciclos fundamentales + filtro de independencia lineal) | `RingInfo`,`SssrResult`,`compute_sssr` | UT 8 ✔ · INT F05 — PASS | RDKit `GetSymmSSSR` **✔ F05** (multiset, no orden) · `distgeom/validation.rs::find_sssr` es una 2ª implementación sin cross-check |
| `src/rings/aromaticity.rs` | 220 | Percepción Hückel 4n+2 en input Kekulé; conteo π por elemento; reescribe enlaces a Aromatic y promueve a sp² | `perceive_aromaticity`,`perceived_aromaticity` | UT 6 ✔ | RDKit `setAromaticity` — **✘ doc afirma paridad, sin fixture** |
| `src/rings/ecfp.rs` | 297 | Morgan/ECFP circular: invariantes, radio 2/3, folding a bitvector, Tanimoto | `ECFPFingerprint`,`compute_ecfp`,`compute_tanimoto` | UT 7 ✔ · INT F59 — PASS | RDKit `GetMorganFingerprint` **✔ F59 (solo ranking Tanimoto, no igualdad de bits)** |
| `src/clustering.rs` | 273 | Butina: matriz RMSD todos-contra-todos (serial + rayon), vecinos, centroides voraces | `ClusterResult`,`compute_rmsd_matrix`,`butina_cluster`,`filter_diverse_conformers` | UT 5 ✔ · INT F14 — PASS | RDKit `Butina.ClusterData` **✔ F14** · `filter_diverse_conformers` sin referencia |
| `src/smirks.rs` | 594 | Parser/aplicador SMIRKS: biyección de atom maps, cambios de orden de enlace, adiciones/eliminaciones, multi-reactivo | `SmirksTransform`,`BondChange`,`SmirksResult`,`parse_smirks`,`apply_smirks`,`apply_smirks_multi` | UT 26 ✔ (el más testado) · INT test_smirks_reactions — PASS | RDKit `rdChemReactions` — **✘ sin paridad de productos** |
| `src/periodic.rs` | 400 | Grafo molecular periódico: detección de enlaces con PBC ±1 imagen, hapticity (η) para metalocenos | `PeriodicMolecule`,`HapticInteraction`,`HapticityAnalysis`,`build_periodic_molecule`,`detect_hapticity` | UT 2 ✔ · INT ninguna | Open Babel / ASE / pymatgen — **✘ sin paridad** |

### 1.2 Química cuántica y estructura electrónica

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/scf/mod.rs` | 44 | Índice de 22 submódulos SCF | — | INT F36/37/38/41 — PASS | — |
| `src/scf/constants.rs` | 31 | Constantes (Å↔Bohr, Hartree↔eV/kcal) + umbrales SCF | `ANGSTROM_TO_BOHR`,`HARTREE_TO_EV`,`SCF_MAX_ITER` | UT 0 | CODATA 2018 — **✘ sin test de valores** |
| `src/scf/types.rs` | 130 | `MolecularSystem` (conversión a Bohr, n_electrons) + `ScfResult` | `MolecularSystem`,`ScfResult`,`IrResult` | UT 0 · INT F36/37/38/41 — PASS | PySCF (implícito vía F36) |
| `src/scf/core_matrices.rs` | 166 | S, T, V, H⁰=T+V + `nuclear_repulsion_energy` | `CoreMatrices`,`nuclear_repulsion_energy*` | UT 5 ✔ · INT F36 — PASS | **PySCF ✔ F36 (≤1e-9 a.u.)** |
| `src/scf/basis.rs` | 940 | STO-3G HSP con renormalización estilo PySCF; primitivas, shells, tabla H..Cl | `BasisSet`,`BasisFunction`,`ContractedShell`,`sto3g_tabulated` | UT 5 ✔ · INT F36/F41 — PASS | **PySCF ✔ F36** · sin paridad para elementos fuera de la tabla |
| `src/scf/gaussian_integrals.rs` | 229 | Integrales primitivas Obara-Saika: teorema del producto, overlap 1D, `boys_function` (serie + asintótico) | `overlap_primitive`,`overlap_1d`,`boys_function` | UT 6 ✔ · INT F36 — PASS | **PySCF/libcint ✔ indirecto F36** · sin paridad a nivel primitiva |
| `src/scf/overlap_matrix.rs` | 110 | Matriz S contractual | `build_overlap_matrix` | UT 3 ✔ · INT F36 — PASS | **PySCF ✔ F36** |
| `src/scf/kinetic_matrix.rs` | 139 | T por recursión Schlegel/HBS | `build_kinetic_matrix` | UT 2 ✔ · INT F36 — PASS | **PySCF ✔ F36** |
| `src/scf/nuclear_matrix.rs` | 286 | V por Boys + expansión Hermite cartesiana; soporta cargas fraccionales | `build_nuclear_matrix`,`build_nuclear_matrix_charges` | UT 3 ✔ · INT F36 — PASS | **PySCF ✔ F36** |
| `src/scf/two_electron.rs` | 672 | ERIs (μν\|λσ) por recursión Obara-Saika completa, caché memoizada, simetría 8-fold, variante rayon | `TwoElectronIntegrals::{compute,compute_parallel,get}` | UT 8 ✔ · INT F37 — PASS | **PySCF/libcint ✔ F37 (≤1e-10)** · path rayon sin gate |
| `src/scf/obara_saika.rs` | 231 | Motor ERI s-only optimizado: `erf_accurate` (Taylor+Lentz), Boys, `eri_ssss`, `schwarz_bound` | `boys_function`,`eri_ssss`,`schwarz_bound`,`compute_eris_screened` | UT 2 ✔ · INT F80 (3/3), alpha_scf_eri_accuracy (12/12) — PASS | **scipy + PySCF/libcint ✔ F80** · solo shells s |
| `src/scf/screening.rs` | 135 | Schwarz pre-screening (Q_ij) | `SchwarzScreener` | UT 4 ✔ · **código muerto (0 callers)** | libcint — **✘ wrapper sin paridad; no está en el pipeline SCF** |
| `src/scf/fock_matrix.rs` | 167 | F = H⁰ + Σ P[(μν\|λσ)−½(μλ\|νσ)]; variante DFTB-SCC + `build_gamma_matrix` | `build_fock_matrix`,`build_fock_dftb`,`build_gamma_matrix` | UT 3 ✔ | PySCF — **✘ sin paridad directa de Fock**; tblite — **✘ sin paridad DFTB** |
| `src/scf/density_matrix.rs` | 241 | P closed-shell, spin density, FON con Fermi smearing, `density_rms_change`, `electron_count` | `build_density_matrix`,`build_density_matrix_fon` | UT 6 ✔ | PySCF — **✘ sin paridad directa de densidad** |
| `src/scf/diis.rs` | 204 | DIIS de Pulay: vector de error, B con Lagrange, LU, extrapolación | `DiisAccelerator` | UT 4 ✔ | PySCF DIIS — **✘ sin paridad directa** |
| `src/scf/energy.rs` | 106 | Descomposición: E_elec=½Tr[P(H⁰+F)], E_1e, E_2e, E_total | `electronic_energy`,`total_energy`,`one_electron_energy` | UT 2 ✔ | PySCF — **✘ sin paridad por término** |
| `src/scf/orthogonalization.rs` | 156 | Löwdin X=UΛ^{-1/2}Uᵀ con threshold | `lowdin_orthogonalization`,`transform_to_orthogonal` | UT 4 ✔ | PySCF — **✘ sin paridad directa de X** |
| `src/scf/mulliken.rs` | 150 | Mulliken (q=Z−Σ(PS)μμ) y Löwdin (S^½PS^½) | `mulliken_analysis`,`lowdin_analysis` | UT 3 ✔ · INT F23/F36 — PASS | **PySCF ✔ F23/F36 (≤0.005 e)** · sin paridad de Löwdin |
| `src/scf/dipole_matrix.rs` | 138 | Integrales de dipolo AO D^α_μν por recursión | `dipole_1d`,`dipole_matrices` | UT 2 ✔ · INT F41 — PASS | **PySCF `int1e_r` ✔ F41 (<1e-8 a.u.)** |
| `src/scf/gradients.rs` | 176 | Gradiente numérico por FD central, densidad energy-weighted W, gradiente analítico de repulsión nuclear | `numerical_gradient`,`energy_weighted_density` | UT 2 ✔ | PySCF/Psi4 — **✘ NO HAY PARITY TEST** · `numerical_gradient` sin caller |
| `src/scf/scf_loop.rs` | 459 | Driver RHF Roothaan-Hall con DIIS, level shift, damping, path paralelo | `ScfConfig`,`run_scf`,`run_scf_hardened` | UT 2 ✔ · INT F36 — PASS | **PySCF ✔ F36 (E≤1e-7 Ha)** · sin paridad de `run_scf_hardened` |
| `src/scf/uhf.rs` | 612 | UHF/ROHF: Fock α/β, densidades de spin, `<S²>`, proyección ROHF, escalamiento | `UhfConfig`,`run_uhf`,`run_rohf` | UT 0 · INT F38 — PASS | **PySCF ✔ F38 (E≤1e-6, ⟨S²⟩≤1e-6)** · ROHF solo H₂⁺ |
| `src/scf/validation.rs` | 137 | Comparadores CPU/GPU de matrices | `compare_matrices`,`validate_overlap_matrix` | UT 4 ✔ · **0 callers (código muerto)** | — **✘** |
| `src/scf/extended_basis.rs` | 386 | 3-21G y 6-31G (H,C,N,O), **sin renormalización de self-overlap** | `BasisSetType`,`build_basis_set` | UT 2 ✔ | PySCF — **✘ sin paridad para 3-21G/6-31G** |
| `src/hf/mod.rs` | 34 | Índice HF-3c | re-exports | — | — |
| `src/hf/api.rs` | 384 | API pública HF-3c: STO-3G → ERI → SCF DIIS → D3+gCP+SRB → CIS → Mulliken | `HfConfig`,`Hf3cResult`,`solve_hf3c`,`solve_hf3c_batch` | UT 3 ✔ | PySCF no tiene HF-3c — **✘ sin paridad end-to-end** |
| `src/hf/basis.rs` | 492 | Basis HF-3c con `Shell::normalize` + `to_scf_basis` | `Shell`,`BasisSet`,`build_sto3g_basis`,`ao_to_atom_map` | UT 2 ✔ · INT F41 — PASS | PySCF ✔ indirecto F41 |
| `src/hf/integrals.rs` | 79 | Repack de ERIs a layout 8-fold compacto | `compute_eris`,`eri_index`,`get_eri` | UT 2 ✔ | PySCF ✔ indirecto F37/F39 |
| `src/hf/fock.rs` | 118 | Fock HF-3c + energía electrónica | `build_fock`,`electronic_energy` | UT 2 ✔ | PySCF — **✘ sin paridad; implementación duplicada** |
| `src/hf/overlap_kin.rs` | 45 | Wrappers sobre `scf::{overlap,kinetic}_matrix` | `compute_overlap_matrix`,`compute_kinetic_matrix` | UT 2 ✔ | PySCF ✔ indirecto F36 |
| `src/hf/nuclear.rs` | 79 | Wrapper de V + **2ª implementación local de Boys** | `compute_nuclear_matrix`,`boys_function` | UT 2 ✔ | PySCF ✔ indirecto F36 · Boys local sin paridad propia |
| `src/hf/scf.rs` | 375 | **2º motor SCF independiente** (Löwdin, DIIS con SVD, level shift) | `ScfResult`,`ScfConfig`,`solve_scf` | UT 1 ✔ | PySCF — **✘ sin paridad: SCF duplicado sin test externo** |
| `src/hf/scf_trait.rs` | 180 | Trait `ScfSolver` + adaptadores HF/PM3/xTB; flags ADIIS **nunca implementados** | `ScfSolver`,`HfScfSolver`,`Pm3ScfSolver`,`XtbScfSolver` | UT 0 · **0 callers** | — **✘** |
| `src/hf/cis.rs` | 371 | CIS con Hamiltoniano singlet, MO-ERI on-the-fly, osciladores con dipolos AO exactos | `Excitation`,`CisResult`,`compute_cis*` | UT 1 ✔ · INT F41 — PASS | **PySCF `tdscf` ✔ F41 (ΔE≤0.01 eV, f≤5e-3)** |
| `src/hf/cisd.rs` | 671 | CISD por determinantes spin-orbital con reglas Slater-Condon; fallback perturbativo CIS+MP2 >5000 det | `CisdResult`,`compute_cisd`,`debug_space` | UT 4 ✔ · INT F41-CISD + cisd_brute — PASS | **PySCF `ci.cisd` ✔ F41** · sin paridad del fallback MP2 |
| `src/hf/mo_transform.rs` | 319 | Transformación AO→MO 4-index (Yoshimine, O(N⁵)) | `MoIntegrals`,`ao_to_mo_transform`,`ao_to_mo_transform_parallel` | UT 1 ✔ · INT F39 — PASS | **PySCF `ao2mo` ✔ F39 (≤1e-10)** · path paralelo sin gate |
| `src/hf/d3.rs` | 389 | D3-BJ: CN logístico, pesos Gaussian, C6 interpolado, damping BJ, término ATM | `D3Result`,`D3Params`,`compute_d3_energy`,`compute_d3_atm_energy` | UT 2 ✔ · INT F46 — PASS | **simple-dftd3 v3 ✔ F46 (≤1.6e-6 Ha)** |
| `src/hf/d3_data.rs` | 24 458 | Tablas D3 generadas (D3_C6, D3_R4R2, radios) | `D3_NREFS`,`D3_REFCN`,`D3_C6`,`D3_VDWRAD_BOHR` | UT 0 (indirecto F46) | simple-dftd3 — **✘ sin comparación byte a byte** |
| `src/hf/gcp.rs` | 99 | gCP: σ=0.1825, α=1.0265, β=1.2148; `n_virtual` es un **proxy** | `compute_gcp` | UT 2 ✔ | HF-3c gCP (Kruse-Grimme) — **✘ sin paridad; forma no publicada** |
| `src/hf/srb.rs` | 89 | SRB: s=0.03, γ=0.7, `exp(−γ(R−R_cov)²)` — **sin `f_damp` publicado** | `compute_srb` | UT 2 ✔ | HF-3c SRB — **✘ sin paridad** |
| `src/dft/mod.rs` | 30 | Índice DFT (gate `alpha-dft`) | re-exports | — | — |
| `src/dft/ks_fock.rs` | 657 | Driver KS-DFT: RKS (F=H⁰+J+V_XC, **sin exchange exacto**) y UKS con DIIS; ECP en core | `DftMethod`,`DftConfig`,`DftResult`,`solve_ks_dft` | UT 9 ✔ · INT F77 — PASS | **PySCF ✔ F77 (E≤1 mHa, gap≤0.05 eV)** · `coulomb_energy=0.0` hard-coded (P2-41) · sin paridad UKS/ECP |
| `src/dft/functionals/svwn.rs` | 147 | SVWN: Slater exchange + VWN-5, ∂V_c por diferencias numéricas | `svwn`,`svwn_spin` | UT 2 ✔ | PySCF `dft.LDA` — **✘ sin paridad a nivel funcional** |
| `src/dft/functionals/pbe.rs` | 177 | PBE: Slater + F_x(s) + PW92 simplificado + corrección de gradiente, derivadas numéricas | `pbe` | UT 2 ✔ | PySCF `dft.GGA` — **✘ sin paridad; correlación simplificada** |
| `src/dft/functionals/mod.rs` | 5 | Índice | — | — | — |
| `src/dft/grid.rs` | 191 | Grid: radios Bragg-Slater, mapeo radial Gauss-Chebyshev, Lebedev, pesos Becke | `MolecularGrid`,`GridQuality`,`GridPoint` | UT 4 ✔ | PySCF grids — **✘ sin estudio de convergencia de grid** |
| `src/dft/lebedev.rs` | 180 | Angular: 6 y 26 exactos; **"110" y "302" son aproximaciones** | `lebedev_grid` | UT 4 ✔ | Lebedev real — **✘ sin paridad; Fine/VeryFine dependen de aproximaciones** |
| `src/dft/becke.rs` | 124 | Partición de Becke 1988 con ajuste heteronuclear | `becke_weights` | UT 2 ✔ | Becke 1988 — **✘ sin paridad de pesos** |
| `src/dft/vxc_matrix.rs` | 265 | V_XC: eval basis+∇φ en grid, ρ y σ, LDA y GGA; `build_vxc_uks` | `build_vxc_matrix`,`build_vxc_uks` | UT 3 ✔ | PySCF — **✘ sin paridad de V_XC** · UKS/PBE clona closed-shell (P2-41) |
| `src/dft/valence_basis.rs` | 233 | Base mínima KS: HSP tabulada + extensión Slater-rule no tabulada; split-valence 3-21G | `DftBasisKind`,`ks_basis` | UT 4 ✔ | PySCF ✔ F77 (solo subset tabulado) · **✘ extensión Slater sin paridad** |
| `src/dft/ecp.rs` | 180 | ECP small-core Z>30 con potencial local **no fitted** | `n_core`,`z_eff`,`local_gaussian`,`ecp_shells` | UT 3 ✔ | Stuttgart/LANL — **✘ sin paridad; docstring admite "not fitted"** |
| `src/dft/gradient.rs` | 175 | Gradiente KS frozen-density con Pulay analítico | `ks_nuclear_gradient` | UT 1 ✔ | PySCF/Psi4 — **✘ sin paridad; solo 3ª ley de Newton** |
| `src/eht/mod.rs` | 36 | Índice EHT | re-exports | INT test_eht, F31, F47 | — |
| `src/eht/params.rs` | 1337 | Tablas Hoffmann: VSIP + exponentes Slater s/p/d; clasificador de soporte por elemento | `SupportLevel`,`EhtSupport`,`EhtParams`,`analyze_eht_support`,`get_params` | UT 8 ✔ · INT F31 — PASS | **yaehmop ✔ F31 (solo metales)** · **✘ sin paridad del grupo principal** |
| `src/eht/basis.rs` | 470 | STO-nG (HSP 1969 / Stewart 1970) con términos armónicos reales d | `AtomicOrbital`,`sto3g_expansion`,`build_basis` | UT 10 ✔ | xtb `slater.f90` — **✘ sin paridad de overlap/energía EHT** |
| `src/eht/overlap.rs` | 227 | S de EHT con recursión 1D estilo Obara-Saika | `build_overlap_matrix` | UT 5 ✔ | libcint — **✘ sin paridad** |
| `src/eht/hamiltonian.rs` | 119 | H: diagonal VSIP, off-diagonal Wolfsberg-Helmholtz (K=1.75) | `build_hamiltonian` | UT 4 ✔ | Hoffmann EHT — **✘ sin paridad de matriz** |
| `src/eht/solver.rs` | 359 | HC=SCE por Löwdin + conteo de electrones de valencia (incl. series 3d/4d/5d) | `EhtResult`,`solve_eht`,`solve_generalized_eigenproblem` | UT 12 ✔ | yaehmop — **✘ sin paridad de eigenvalores/gaps** |
| `src/eht/metal_fallback.rs` | 313 | Parámetros hidrogénicos derivados: ζ=√(IP/13.6057), cacheados en `OnceLock` | `derived_params`,`is_derived`,`valence_electrons` | UT 5 ✔ | NIST IP — **✘ sin paridad de la tabla IP** |
| `src/eht/gradients.rs` | 262 | Gradiente de superficie EHT por FD central + `optimize_geometry_eht` | `EhtGradient`,`EhtOptResult`,`compute_eht_gradient` | UT 2 ✔ | — **✘ sin paridad; doc admite superficie cualitativa** |
| `src/eht/band_structure.rs` | 403 | Band structure periódica: k-path, H(k)/S(k) con 6 imágenes y fase `cos` real | `KPoint`,`BandStructure`,`compute_band_structure` | UT 2 ✔ | QE/GPAW/CASTEP — **✘ sin paridad; sin parte imaginaria** |
| `src/eht/volume.rs` | 417 | Mapeo volumétrico de orbitales, densidad, variantes paralela y por slabs | `VolumetricGrid`,`evaluate_orbital_on_grid*` | UT 7 ✔ | Avogadro/Multiwfn — **✘ sin paridad numérica** |
| `src/eht/marching_cubes.rs` | 977 | Marching Cubes estándar + variante dual, simplificación, normales ponderadas | `IsosurfaceMesh`,`marching_cubes`,`simplify_mesh` | UT 8 ✔ | VTK/skimage — **✘ sin paridad de topología** |
| `src/pm3/mod.rs` | 35 | Índice PM3 | re-exports | — | — |
| `src/pm3/faithful.rs` | 1794 | Port fiel de MOPAC 23.2.5 término a término: overlaps diat/set, 22 integrales locales (`reppd`), rotmat, spcore, ccrep, calpar, Fock NDDO (`fock2`/`jab`/`kab`), ΔHf | `Pm3Derived`,`FaithfulResult`,`solve` | UT 4 ✔ · INT F42 — PASS | **MOPAC 23.2.5 ✔ F42 (ΔHf <0.02 kcal/mol, 24 mols)** · solo ΔHf, Z≤17, closed-shell |
| `src/pm3/solver.rs` | 846 | Solver PM3 NDDO genérico/legacy; **overlap STO s-p = 0.5×s-s** (aproximación) | `Pm3Result`,`solve_pm3` | UT 4 (1 `#[ignore]`d) · INT F42 | MOPAC ✔ F42 (path fiel) · **✘ path legacy Z>17 sin referencia** |
| `src/pm3/params.rs` | 371 | Accesores de parámetros PM3 + set PM3(tm) (Cundari 1998) | `Pm3Params`,`get_pm3_params`,`pm3_method` | UT 6 ✔ | MOPAC ✔ indirecto F42 · sin paridad por parámetro |
| `src/pm3/params_table.rs` | 1022 | Tabla PM3 main-group generada | `PM3_MAIN_GROUP` | UT 0 (indirecto F42) | MOPAC — **✘ sin comparación byte a byte** |
| `src/pm3/gradients.rs` | 321 | Gradientes PM3: W=2ΣεCC, HF+Pulay, derivadas analíticas de γ y core-core, micro-numéricas del overlap | `Pm3GradientResult`,`compute_pm3_gradient` | UT 3 ✔ | MOPAC — **✘ sin paridad de fuerzas** |
| `src/pm3/gpu.rs` | 145 | GPU del G diagonal de dos centros (WGSL f32) | `build_pm3_g_matrix_gpu` | UT 1 ✔ | — **✘ sin paridad GPU/CPU** |
| `src/xtb/mod.rs` | 36 | Índice xTB | re-exports | — | — |
| `src/xtb/gfn0.rs` | 1026 | Port GFN0-xTB de xtb 6.7.1: EEQ con erf, CN mctc, H⁰ con `dtrf2`, repulsión, SRB, D4(BJ) | `Gfn0Result`,`gfn0_*`,`solve_gfn0` | UT 1 ✔ · INT F43 — PASS | **xtb 6.7.1 ✔ F43 (≤0.01 eV)** · solo closed-shell; sin paridad por término |
| `src/xtb/gfn0_params.rs` | 2465 | Parámetros GFN0 de `param_gfn0-xtb.txt` | `Gfn0ElementParams`,`GFN0_GLOBALS`,`get_gfn0_params` | UT 0 (indirecto F43) | xtb — **✘ sin aserción por parámetro** |
| `src/xtb/gfn1.rs` | 886 | Port GFN1-xTB de tblite 0.7.0: kpair/enscale, electrostática 2º orden por shell, 3º orden, halogen bond, Broyden | `Gfn1Result`,`solve_gfn1` | UT 1 ✔ · INT F44 — PASS | **tblite 0.7.0 ✔ F44 (≤0.01 eV)** · sin paridad por término |
| `src/xtb/gfn1_params.rs` | 1720 | Parámetros GFN1 de `gfn1.f90` | `Gfn1ElementParams`,`GFN1_GLOBALS` | UT 0 (indirecto F44) | tblite — **✘ sin aserción por parámetro** |
| `src/xtb/gfn2.rs` | 1671 | GFN2-xTB: AES multipolo, 3º orden Hubbard, Fermi filling con Newton, D4 self-consistent, halogen bond | `Gfn2Result`,`solve_gfn2` | UT 6 ✔ · INT F45 — PASS | **tblite 0.7.0 ✔ F45 (≤0.0217 eV)** · tolerancia más laxa que F44 |
| `src/xtb/gfn2_params.rs` | 871 | Parámetros GFN2 oficiales | `Gfn2ElementParams`,`get_gfn2_params`,`effective_hubbard` | UT 4 ✔ | tblite ✔ indirecto F45 |
| `src/xtb/solver.rs` | 548 | **Solver "GFN0-inspired" legacy** (expone `solve_xtb`) | `XtbResult`,`solve_xtb` | UT 4 ✔ | — **✘ sin paridad: duplicado legacy exportado** |
| `src/xtb/params.rs` | 564 | Parámetros simplificados nivel GFN0 legacy | `XtbParams`,`get_xtb_params` | UT 3 ✔ | — **✘ sin paridad** |
| `src/xtb/d4.rs` | 837 | D4 self-consistente: pesos con ζ-scaling, C6 Casimir-Polder, damping BJ, potencial para Fock, ATM | `D4Model`,`D4Weights`,`D4Damping` | UT 2 ✔ · INT F46 — PASS | **dftd4 v4.3.0 ✔ F46 (≤1.6e-6 Ha)** · `get_potential` sin paridad |
| `src/xtb/d4_data.rs` | 4676 | Datos dftd4 v4.3.0 (Z=1..118) | `R4R2`,`CP_FREQ`,`REFSYS`,… | UT 0 (indirecto F46) | dftd4 — **✘ sin comparación byte a byte** |
| `src/xtb/broyden.rs` | 331 | Mixer de Broyden modificado para el loop SCC | `BroydenMixer` | UT 2 ✔ | tblite `scf/broyden` — **✘ sin paridad de trayectoria** |
| `src/xtb/spin_params.rs` | 115 | Constantes W de `data/spin.f90` (Hartree) | `SPIN_W`,`spin_w` | UT 0 (indirecto vía uks) | **tblite ✔ indirecto** |
| `src/xtb/sto_overlap.rs` | 1396 | Expansión STO-nG (Stewart 1970) + overlaps cartesianos + momentos multipolo (monopolo/dipolo/cuadrupolo traceless) | `sto_overlap_ng`,`sto_multipole_with_ng`,`MultipoleResult` | UT 13 ✔ | **tblite ✔ parcial (1 test de overlap)** · **✘ sin paridad de multipolos** |
| `src/xtb/uks.rs` | 1291 | GFN1/GFN2 UKS colineal con magnetización por shell y E_spin=½ΣmWm | `UksResult`,`solve_gfn2_unrestricted`,`solve_gfn1_unrestricted` | UT 1 ✔ · INT F44-UKS, F45-UKS — PASS | **tblite `uhf` ✔ F44/F45-UKS (≤0.05 eV)** · término de spin sin equivalente upstream |
| `src/xtb/gradients.rs` | 258 | Gradientes analíticos xTB (path legacy) | `XtbGradientResult`,`compute_xtb_gradient` | UT 3 ✔ | xtb/tblite — **✘ sin paridad de fuerzas** |
| `src/xtb/gpu.rs` | 215 | GPU de las 3 γ-matrices xTB (GFN0/1/2) | `build_xtb_gamma_gpu` | UT 2 ✔ | — **✘ sin paridad GPU/CPU** |
| `src/dispersion/mod.rs` | 26 | Índice D4 | re-exports | — | — |
| `src/dispersion/dispersion.rs` | 210 | D4 standalone vía `xtb::d4` con cargas EEQ2019; gradiente por FD | `D4Config`,`D4Result`,`compute_d4_energy`,`compute_d4_gradient` | UT 4 ✔ · INT F46 — PASS | **dftd4 v4.3.0 ✔ F46** · gradiente numérico sin paridad |
| `src/dispersion/params.rs` | 516 | Helpers D4 standalone (valores aproximados a mano) | `D4Params`,`get_d4_params`,`dynamic_c6` | UT 3 ✔ | dftd4 — **✘ el path de producción no usa estos helpers** |
| `src/dos/mod.rs` | 10 | Índice DOS | re-exports | — | — |
| `src/dos/dos.rs` | 519 | DOS Gauss + PDOS por pesos orbital-átomo (Mulliken-like) | `DosResult`,`compute_dos`,`compute_pdos` | UT 9 ✔ | xtb/ORCA/GPAW — **✘ sin paridad; ni integral vs conteo electrónico** |
| `src/dos/multi_method.rs` | 151 | DOS multi-método (EHT/PM3/GFN0/GFN1/GFN2/HF-3c) | `DosMethod`,`compute_dos_multimethod` | UT 1 (solo EHT) | — **✘ sin paridad del broadening; branches sin cobertura** |
| `src/mesh.rs` | 298 | Mesh orbital genérico; **usa base EHT para PM3/xTB/HF-3c** | `MeshMethod`,`OrbitalMeshResult`,`compute_orbital_mesh` | UT 4 ✔ | Avogadro/Multiwfn/VMD — **✘ sin paridad** · forma espacial usa base EHT |

### 1.3 Cargas, superficie, electrostática, población, reactividad, estereo

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/charges/mod.rs` | 13 | Fachada | re-exports | — | — |
| `src/charges/gasteiger.rs` | 747 | Gasteiger-Marsili PEOE: χ=a+bq+cq², Δq por enlace, damping 0.5, pre-paso de carboxilato | `GasteigerParams`,`GasteigerMode`,`gasteiger_marsili_charges_full`,`ChargeResult` | UT 8 ✔ · INT F21, test_charges, algorithm_validation — PASS | **RDKit ✔ F21 (max \|Δq\|=8.3e-17 e)** · OpenBabel sin paridad |
| `src/charges_eeq/mod.rs` | 33 | Fachada | re-exports 12 items | — | — |
| `src/charges_eeq/charges.rs` | 441 | EEQ legacy: CN fraccional tipo Fermi, Coulomb con erf, CN-shift de χ, solve (n+1)×(n+1) | `EeqParams`,`EeqConfig`,`compute_eeq_charges` | UT 4 ✔ · INT test_eeq, e3_justification — PASS | multicharge/dftd4 — **✘ sin paridad** · **γ en 1/Å mezclado con η en eV** |
| `src/charges_eeq/energy.rs` | 126 | E_EEQ + gradiente FD; **σ distinto al de `charges.rs`** | `EeqEnergyResult`,`compute_eeq_energy`,`compute_eeq_gradient` | UT 2 ✔ (solo finitud) | dftd4 — **✘ sin paridad; E ≠ energía de los charges** |
| `src/charges_eeq/eeq2019.rs` | 151 | Port fiel de multicharge `eeq2019` (a.u.) | `compute_eeq2019_charges` | UT 0 · INT ninguna | multicharge — **✘ CERO tests; usado en producción por D4** |
| `src/charges_eeq/eeq2019_data.rs` | 123 | Tablas EEQ2019 (χ,η,rad,kcnchi) Z=1..103 | `EEQ2019_CHI/ETA/RAD/KCNCHI` | UT 0 | multicharge — **✘ sin paridad de tablas** |
| `src/charges_eeq/eeqbc.rs` | 327 | Port EEQ-BC 2025 (JCP 162, 214109): erf exacto, qloc, matriz de capacitancia | `compute_eeqbc_charges`,`eeqbc_coordination_numbers` | UT 0 · INT F22 — PASS | **multicharge `eeqbc` ✔ F22 (max \|Δq\|=6.9e-12 e)** |
| `src/charges_eeq/eeqbc_data.rs` | 1669 | Parámetros EEQ-BC (5356 radios vdw + EN Pauling) | 10 tablas `EEQBC_*` | UT 0 (indirecto F22) | multicharge/mctc-lib — **✘ sin paridad de actinidos** |
| `src/dipole/mod.rs` | 5 | Fachada | re-exports | — | — |
| `src/dipole/dipole.rs` | 184 | μ=ΣqR (e·Å→Debye); split nuclear/electrónico desde EHT | `DipoleResult`,`compute_dipole`,`compute_dipole_from_eht` | UT 4 ✔ · INT F26, phase_c — PASS | **PySCF ✔ F26 (≤1.8e-7 D)** · `compute_dipole_from_eht` sin paridad · NIST `dipole_debye` sin usar |
| `src/surface/mod.rs` | 28 | Fachada | re-exports | — | — |
| `src/surface/sasa.rs` | 335 | Shrake-Rupley: esfera Fibonacci (960 pts), radio VdW+probe, área 4πR²·n_acc/n | `SasaResult`,`compute_sasa`,`compute_sasa_parallel` | UT 8 ✔ · INT F28, test_sasa — PASS | **FreeSASA/RDKit ✔ F28 (max rel 1.53 %)** · sin test parallel==serial |
| `src/esp/mod.rs` | 12 | Fachada | re-exports | — | — |
| `src/esp/esp.rs` | 835 | ESP point-charge en grid 3D con clamp nuclear, mapa RGB, escritura/lectura `.cube` | `EspGrid`,`EspGridConfig`,`compute_esp_grid*`,`export_cube`,`read_cube` | UT 12 ✔ · INT phase_c — PASS | APBS/pdb2pqr, Multiwfn, Psi4 — **✘ SIN PARIDAD** · `elements` ignorado · `.cube` escribe cargas |
| `src/population/mod.rs` | 15 | Fachada | re-exports | — | — |
| `src/population/population.rs` | 961 | Mulliken/Löwdin desde EHT + Wiberg-like (P² Löwdin) y Mayer-like (PS·PS) | `PopulationResult`,`BondOrderResult`,`compute_population*`,`compute_bond_orders*` | UT 19 ✔ · INT phase_c — PASS | PySCF — **✘ sin paridad (F36 cubre `scf/mulliken.rs`)** · Wiberg no es el estándar ΣPμνSμν |
| `src/population/npa.rs` | 403 | NPA/NBO simplificado: S^{-1/2}DS^{-1/2}, block-diagonalización, ocupaciones naturales, lone pairs | `NpaResult`,`NboResult`,`compute_npa`,`compute_nbo` | UT 3 (2×2 de juguete) · INT ninguna | **Psi4/GenNBO/NBO 6.0, Multiwfn — ✘ SIN PARIDAD** · Z nuclear en vez de valencia |
| `src/reactivity.rs` | 945 | Descriptores de frontera EHT, Fukui condensados, dual descriptor, ranking de sitios, pKa empírico, corrección aromática UFF, UV-Vis estilo sTDA | `FrontierDescriptors`,`FukuiDescriptors`,`ReactivityRanking`,`EmpiricalPkaResult`,`StdaUvVisSpectrum` | UT 6 ✔ · INT F47, roadmap_validation, reaction_battery — PASS | **F47 ✔ pero contra reglas, no contra QM** · **✘ sin paridad ORCA/Psi4 para Fukui** · pKa sin tabla experimental |
| `src/potentials.rs` | 58 | Potenciales analíticos: doble well 1D, SN2, Müller-Brown 2D, EVB, Morse 2D | `double_well_1d`,`muller_brown_2d`,`proton_transfer_2d_morse` | UT 0 · INT reaction_battery, alpha_reaction_dynamics_3d — PASS | F82 ✔ (parcial, cubre el algoritmo GSM) · `proton_transfer_2d_morse` sin test |
| `src/stereo.rs` | 1115 | CIP: árboles recursivos con phantom atoms, R/S por volumen firmado, E/Z, atropisomería, proquiral, helicoidal (aleno/ heliceno/metaloceno) | `StereoAnalysis`,`analyze_stereo`,`stereo_descriptors` | UT 6 ✔ · INT algorithm_validation, new_features, ci — PASS | **RDKit `AssignStereochemistry` — ✘ SIN PARIDAD** · **signo R/S aparentemente invertido vs OpenSMILES** · sin CIP 4/5 |
| `src/cpm/mod.rs` | 11 | Fachada CPM | re-exports | — | — |
| `src/cpm/grand_potential.rs` | 300 | CPM: Ω=Σχq+½Σηq²+ΣJqJq−μΣq, J=14.3996/(εr); solve lineal LU | `CpmConfig`,`CpmResult`,`compute_cpm_charges` | UT 3 ✔ · INT test_cpm (F95) — PASS | **LAMMPS `fix electrode/conp` ✔ F95 (ratio 1.76×, cos 0.44)** · gap de modelo documentado |
| `src/cpm/surface.rs` | 140 | Barrido μ → Q(μ), Ω(μ), C=dQ/dμ | `CpmSurface`,`compute_cpm_surface` | UT 3 ✔ · INT test_cpm — PASS | — **✘ sin paridad de la curva C(μ)** |
| `src/beta/cpm/mod.rs` | 11 | Shim de re-export del CPM core | re-exports | UT 0 | — shim sin paridad propia |

### 1.4 ML, ANI, descriptores

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/ml/mod.rs` | 31 | Índice ML; `inference`/`mlff`/`symmetry_functions` tras gate `alpha-mlff` | re-exports | — | — |
| `src/ml/descriptors.rs` | 465 | Descriptores 2D (MW, HBD/HBA, rotatable, Fsp3, Wiener BFS, anillos, **Balaban placeholder**) + 3D (giroCardano: asphericidad, NPR, esfericidad) | `MolecularDescriptors`,`Descriptors3D`,`compute_descriptors`,`compute_3d_descriptors` | UT 2 ✔ · INT new_features, extended_molecules — PASS | **RDKit Descriptors / mordred — ✘ SIN PARIDAD (F56 sin fixture)** · `balaban_j` = `\|ln µ\|+1` |
| `src/ml/models.rs` | 283 | Surrogados lineales: logP estilo Crippen, refractividad, logS Delaney, druglikeness | `MlPropertyResult`,`predict_properties` | UT 3 ✔ · INT extended_molecules — PASS | **RDKit Crippen / ESOL — ✘ sin paridad; coeficientes a mano** |
| `src/ml/inference.rs` | 259 | Red feed-forward: GELU/CELU/ReLU, forward, backprop manual reverse-mode | `Activation`,`DenseLayer`,`InferenceNet`,`MlffResult` | UT 3 ✔ · INT alpha_mlff_consistency (10), alpha_mlff_validation (2) — PASS | PyTorch — **✘ sin paridad** · F79 vs oráculo numpy (5.1e-14) |
| `src/ml/mlff.rs` | 227 | Pipeline MLFF: AEV → NN → energías atómicas; **fuerzas por FD central, no backprop** | `MlffConfig`,`compute_mlff`,`compute_mlff_energy_and_gradient` | UT 3 ✔ · INT F79 — PASS | **TorchANI / ASE — ✘ sin paridad; sin pesos entrenados ni benchmark MD17** |
| `src/ml/symmetry_functions.rs` | 260 | Behler-Parrinello generalizadas (species-agnostic) radiales+angulares, paralelo | `SymmetryFunctionParams`,`Aev`,`compute_aevs` | UT 5 ✔ · INT F79 — PASS | **TorchANI — ✘ sin paridad** · **bug: solo 16 de 256 slots radiales se llenan (P2-36)** |
| `src/ml/advanced_models.rs` | 714 | CART puro Rust: RF con bagging, GBM L2, CV k-fold, recalibración isotónica PAVA | `RandomForest`,`GradientBoosting`,`cross_validate`,`IsotonicCalibrator` | UT 5 ✔ (datos sintéticos) | **scikit-learn — ✘ SIN PARIDAD; sin benchmark MoleculeNet; `oob_score` siempre None** |
| `src/ml/ensemble.rs` | 493 | Consenso: 3 modelos de logP, TPSA Ertl, pKa heurístico, Veber, BBB | `EnsembleResult`,`VeberResult`,`compute_tpsa`,`predict_ensemble` | UT 7 ✔ | **RDKit `Descriptors.TPSA` — ✘ sin paridad** |
| `src/ml/getaway.rs` | 298 | GETAWAY: matriz influencia H=X(XᵀX)⁻¹Xᵀ, leverage, HATk/Rk, contenido de información | `GetawayDescriptors`,`compute_getaway` | UT 2 ✔ | **mordred/PaDEL/CDK — ✘ sin paridad (F57 sin fixture)** |
| `src/ml/pharmacophore.rs` | 340 | Percepción de farmacóforos + fingerprint hashed de 2 puntos con bins de 2 Å | `PharmFeature`,`PharmacophoreFingerprint`,`detect_features` | UT 3 ✔ · INT new_features — PASS | **RDKit `Pharm2D` — ✘ sin paridad (F58 sin fixture); layout de bits bespoke** |
| `src/ml/rdf_descriptors.rs` | 203 | RDF suavizado gaussiano (β=100 Å⁻²) sobre 128 radios a 12 Å, 5 ponderaciones | `RdfDescriptors`,`compute_rdf` | UT 2 ✔ | **mordred/PaDEL — ✘ sin paridad; normalización 4πr² no es la de Mordred** |
| `src/ml/whim.rs` | 405 | WHIM: covarianza ponderada, Jacobi, λ/ν/γ/η por eje, T/A/D/K, 5 esquemas | `WhimDescriptors`,`WhimWeighting`,`compute_whim` | UT 3 ✔ | **mordred/PaDEL — ✘ sin paridad; γ/η con momentos crudos** |
| `src/ani/mod.rs` | 21 | Índice ANI | re-exports | — | — |
| `src/ani/api.rs` | 197 | Driver ANI: validación de especies → cell list → AEV → NN → energía → fuerzas | `AniConfig`,`AniResult`,`compute_ani`,`compute_ani_test` | UT 2 ✔ · INT ninguna | **TorchANI/ASE — ✘ sin paridad** · formato de pesos bespoke |
| `src/ani/aev.rs` | 248 | AEV Behler-Parrinello con split por especie (radial + pares angulares) | `compute_aevs` | UT 2 ✔ (débiles) | **TorchANI — ✘ sin paridad** · **BUG: alias de 6 de 28 bloques angulares (P0-6)** |
| `src/ani/aev_params.rs` | 112 | Tablas ANI-2x de hiperparámetros (valores "representativos", no los publicados) | `ANI_ELEMENTS`,`AevParams`,`default_ani2x_params` | UT 2 ✔ | Devereux JCTC 2020 — **✘ sin paridad** |
| `src/ani/ani_tm.rs` | 287 | ANI-TM a 24 elementos, cutoffs 6.5/4.5, packing correcto de pares | `ANI_TM_ELEMENTS`,`compute_aevs_tm`,`is_ani_tm_supported` | UT 3 ✔ | torchani TM — **✘ sin paridad** · "transfer learning" solo documentado |
| `src/ani/cutoff.rs` | 78 | Cutoff coseno f_c(r) y su derivada analítica | `cosine_cutoff`,`cosine_cutoff_deriv` | UT 4 ✔ (incl. FD) | Behler-Parrinello — fórmula idéntica, sin fixture externo necesario |
| `src/ani/gradients.rs` | 204 | "Fuerzas analíticas": NN backward → dE/dAEV → cadena sobre funciones radiales | `compute_forces` | UT 1 (solo ΣF≈0) | **TorchANI/ASE — ✘ sin paridad** · **BUG: solo bloque radial (P0-6)** |
| `src/ani/neighbor.rs` | 218 | Cell-list O(N) de vecinos | `CellList`,`NeighborPair` | UT 2 ✔ | ASE `neighborlist`/scipy cKDTree — **✘ sin paridad** |
| `src/ani/nn.rs` | 226 | Red nalgebra con GELU/CELU + backprop | `FeedForwardNet`,`DenseLayer` | UT 2 ✔ (backward vs FD) | PyTorch — **✘ sin paridad; GELU con erf aproximado (1.5e-7)** |
| `src/ani/weights.rs` | 174 | Serializador binario de redes por elemento + generador de red de test | `load_weights`,`save_weights`,`make_test_model` | UT 1 ✔ | **ANI-2x oficial `.netFormat`/HDF5 — ✘ formato bespoke, pesos oficiales no cargables** |

### 1.5 Espectroscopia, IR, NMR

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/spectroscopy/mod.rs` | 36 | Fachada | re-exports | INT test_spectroscopy (34) — PASS | — |
| `src/spectroscopy/types.rs` | 102 | Tipos compartidos: `ScfInput`, `TransitionInfo`, `ShieldingTensor` | `ScfInput`,`TransitionInfo`,`SpectroscopyResult` | UT 0 · INT ninguna | — **UT NONE** |
| `src/spectroscopy/hessian.rs` | 386 | Hessiana semi-numérica 6N, simetrización, mass-weighting, eigendecomposición, cm⁻¹ | `HessianResult`,`compute_hessian*` | UT 2 ✔ | Psi4/PySCF/NIST — **✘ sin paridad** |
| `src/spectroscopy/ir_intensities.rs` | 197 | I_k=\|dμ/dQ\|² con 42.256 a.u.→km/mol; detección de transl/rot | `IrConfig`,`IrPeak`,`compute_ir_intensities` | UT 0 (solo smoke bench) | **NIST/CCCBDB — ✘ sin paridad; `broadening`/`freq_min`/`freq_max` sin usar** |
| `src/spectroscopy/stda_uvvis.rs` | 622 | sTDA: espacio activo, A=Δε+2×Coulomb monopolar (η=0.3), prescreening C-S, osciladores | `StdaConfig`,`compute_stda`,`compute_stda_with_dipoles` | UT 7 ✔ | **PySCF `tdscf.TDA` / Grimme sTDA — ✘ sin paridad** · `config.ax` sin usar (sin exchange) |
| `src/spectroscopy/transition_dipoles.rs` | 190 | Dipolo AO→MO con overlap ponderado; f=(2/3)ΔE\|μ\|² | `DipoleIntegrals`,`oscillator_strength`,`transition_dipole_from_ci_vector` | UT 3 ✔ | PySCF `ao2mo` — **✘ sin paridad de este módulo** |
| `src/spectroscopy/fluorescence.rs` | 308 | Fluorescencia desde sTDA-EHT: regla de Kasha, I(E)∝n²E³\|μ\|²g(E), A de Einstein, vida radiativa | `FluorescenceConfig`,`FluorescenceSpectrum`,`einstein_a_coefficient` | UT 6 ✔ | TD-DFT + vidas experimentales — **✘ sin paridad; sin opt. de geometría excitada** |
| `src/spectroscopy/giao_nmr.rs` | 307 | "GIAO": constante diamagnética de átomo libre + corrección empírica + SOS de solapamiento (×100) | `compute_nmr_shieldings*`,`shieldings_to_shifts` | UT 6 ✔ (estructurales) | **Psi4/PySCF NMR — ✘ sin paridad; no es GIAO (P1-23)** |
| `src/ir/mod.rs` | 20 | Fachada IR | re-exports | — | — |
| `src/ir/hessian.rs` | 588 | Hessiana numérica por FD de energías (EHT/PM3/GFN2) + semi-analítica + UFF analítica | `HessianMethod`,`compute_numerical_hessian*` | UT 2 ✔ | **xtb `vibspectrum` ✔ parcial (F48 cubre diagonalización, no estos stencils)** · sin paridad vs Psi4/PySCF |
| `src/ir/vibrations.rs` | 1225 | Workflow vibracional: base de modos rígidos + proyector, eigensolve, intensidades IR, RRHO a 298.15 K, ensanchado | `VibrationalAnalysis`,`Thermochemistry`,`compute_vibrational_analysis*`,`compute_ir_spectrum*` | UT 6 ✔ · INT F48, test_spectroscopy, algorithm_validation — PASS | **xtb ✔ F48 / NIST ✔ F49** · **✘ sin paridad de intensidades km/mol** |
| `src/ir/peak_assignment.rs` | 870 | Asignación de picos IR por tabla de 60 grupos con gating por elementos y amplitud dominante | `GroupFrequency`,`PeakAssignment`,`assign_peaks` | UT 3 ✔ · INT F49 — PASS | **NIST/CCCBDB ✔ F49 (recall 1.00, 10 mols)** |
| `src/nmr/mod.rs` | 29 | Fachada NMR | re-exports | — | — |
| `src/nmr/shifts.rs` | 1195 | Shifts empíricos ¹H/¹³C/¹⁹F/³¹P/¹⁵N/¹¹B/²⁹Si/⁷⁷Se/¹⁷O/³³S por reglas de entorno; offsets de isótopo | `NmrShiftResult`,`predict_chemical_shifts*` | UT 5 (1 sin `#[test]`, P2-11) · INT algorithm_validation — PASS | **nmrdb2 / NIST SDBS — ✘ sin paridad (spot-check 3 mols)** |
| `src/nmr/nucleus.rs` | 1467 | Registro de 85 isótopos + catálogo de 33 metales (frecuencia, spin, FWHM, rangos) | `NmrNucleus` (16 métodos) | UT 0 · INT indirecta | **CODATA/CRC — ✘ sin paridad de frecuencias/FWHM** |
| `src/nmr/coupling.rs` | 633 | Karplus con coeficientes por vía (Altona-Sundaralingam, Bystrov), reglas topológicas ²J/⁴J/⁵J, HLA generalized + Δχ de Huggins, promedio Boltzmann | `JCoupling`,`KarplusParams`,`predict_j_couplings`,`ensemble_averaged_j_couplings` | UT 6 ✔ · INT F54 — PASS | **Haasnoot/de Leeuw/Altona ✔ F54 (base 4e-3 Hz, 49 casos)** · `hla_coupling` no está en el path de producción |
| `src/nmr/spectrum.rs` | 700 | Espectro 1D sintético: agrupación, multiplets Pascal, Lorentzian, integraciones, flag de 2º orden | `NmrSpectrum`,`compute_nmr_spectrum*` | UT 7 ✔ | nmrdb2/experimental — **✘ sin paridad; solo 1º orden** |
| `src/nmr/hose.rs` | 390 | Códigos HOSE por esferas concéntricas + base de 51 patrones H→shift | `HoseCode`,`HoseShiftLookup`,`generate_hose_codes` | UT 3 ✔ | nmrdb2 HOSE — **✘ sin paridad; `predict_shift_from_hose` es código muerto** |
| `src/nmr/second_order.rs` | 385 | NMR exacto 2º orden: H de spin-½ en base producto, eigen-decomp, transiciones 1Q, efecto techo AB | `SpinSystem`,`ExactLine`,`ab_reference` | UT 6 ✔ | nmrglue/DAST — **✘ sin paridad externa (solo referencia analítica interna)** |

### 1.6 Materiales, periódico, CGA, SDR, Riemanniano

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/materials/mod.rs` | 28 | Índice; declara un **2º tipo `SpaceGroup`**; no re-exporta geometry_opt/space_groups/variable_cell | `pub use assembly,cell,cif,sbu` | INT F67/F68/F70, ci — PASS | — |
| `src/materials/cell.rs` | 923 | `UnitCell` con base en filas, frac↔cart, wrap, min-image, supercel; XRD en polvo (Bragg + Waasmaier-Kirfel + LP); porosidad (voxel + BFS percolante + PLD) | `UnitCell`,`PowderXrdResult`,`simulate_powder_xrd`,`PorosityResult`,`compute_porosity` | UT 7 ✔ · INT F68 (porosidad) — PASS · XRD solo smoke | **pymatgen `XRDCalculator` — ✘ fixture `pymatgen_xrd.json` sin test** · GSAS-II, pyzeo ✔ F68 |
| `src/materials/assembly.rs` | 644 | 10 redes RCSR built-in (pcu/dia/sql/bcu/fcu/nbo/pts/kgm/hxl/pillared_sql) + ensamblaje con rotación de linker (Rodrigues) | `Topology`,`CrystalStructure`,`assemble_framework` | UT 16 ✔ · INT F67 — PASS | **pormake ✔ F67 (MOF-5 + coordinación RCSR)** · ignora connection points de SBU |
| `src/materials/sbu.rs` | 247 | Modelo SBU: átomos + puntos de conexión + geometría de coordinación; `ideal_directions` | `Sbu`,`ConnectionPoint`,`CoordinationGeometry` | UT 6 ✔ | **pormake SBU library — ✘ sin paridad** |
| `src/materials/cif.rs` | 697 | Lector CIF 1.1 (celda, space group, loop `_atom_site_*`), quita incertidumbres `5.640(1)`, escribe CIF mínimo | `CifStructure`,`CifAtomSite`,`parse_cif`,`write_cif` | UT 6 ✔ · INT F68 — PASS | **pymatgen/ASE — ✘ fixture `pymatgen_cif.json` sin test** · elementos >Bi → Z=0 |
| `src/materials/space_groups.rs` | 265 | 230 grupos ITC: `SymmetryOp`, lookup por número/símbolo, posiciones equivalentes, multiplicidad | `SpaceGroup`,`SymmetryOp`,`space_group_by_number`,`space_group_by_symbol` | UT 7 ✔ (solo existencia) | **spglib/pymatgen — ✘ fixture `space_groups_230.json` (412 KB) sin test** |
| `src/materials/space_group_table.rs` | 6751 | Tabla generada de 230 grupos con matrices i8 + traslaciones en duodécimas | `group_table(n)` (crate) | UT 0 (datos) | spglib/pymatgen — **✘ sin generador versionado → deriva silenciosa** |
| `src/materials/topology_data.rs` | 122 | 6 redes RCSR (acs/afw/sod/pts/xah/tfy) | `TOPOLOGY_NETS`,`NetData` | UT 0 · INT F67 — PASS | **pormake ✔ F67** · `edges` nunca consumido, `edges.len()!=edge_types.len()` en 4/6 (P2-33) |
| `src/materials/geometry_opt.rs` | 627 | Relajación de celda fija: BFGS con Sherman-Morrison-Woodbury + Armijo, o steepest descent adaptativo; constraints PBC | `OptMethod`,`FrameworkOptResult`,`optimize_framework` | UT 4 ✔ · INT F70 — PASS | **RDKit UFF ✔ F70 (RMSD≤0.01 Å)** · solo molecular, path PBC sin test |
| `src/materials/variable_cell.rs` | 505 | Relajación variable: `numerical_stress`, deformación afín, paso de strain isotropico con backtracking | `VariableCellResult`,`numerical_stress`,`optimize_variable_cell` | UT 3 ✔ (modelo LJ interno) | **VASP ISIF / CP2K / QE vc-relax / LAMMPS fix box/relax — ✘ SIN PARIDAD** · ramas de stress idénticas (P2) |
| `src/materials/scattering_factors.rs` | 392 | Cromer-Mann / Waasmaier-Kirfel (Z=1..98) | `atomic_scattering_factor` (crate) | UT 0 · INT smoke NaCl | **pymatgen XRDCalculator — ✘ fixture existe sin test; summary.json dice value:null** |
| `src/kmesh.rs` | 795 | Mesh recíproco: Monkhorst-Pack (MP + Γ), fases de Bloch, H(k) con `cos` real, integración BZ, validación de electrones | `KMesh`,`KPoint`,`monkhorst_pack_mesh`,`bloch_phase`,`build_bloch_hamiltonian` | UT 15 ✔ · INT F88 — PASS | **ASE `monkhorst_pack` ✔ F88** · sin paridad de eigenvalores vs QE/VASP |
| `src/cga/mod.rs` | 21 | Fachada CGA | re-exports | INT F81, test_cga — PASS | — |
| `src/cga/multivector.rs` | 436 | Álgebra G(4,1): 32 componentes, métrica [1,1,1,1,−1], producto geométrico, reverse, involute, grade | `Multivector` (23 métodos) | UT 5 ✔ · INT F81 — PASS | **clifford Cl(4,1) ✔ F81 (3.4e-15)** |
| `src/cga/motor.rs` | 254 | Motor = translator∘rotor; `conformal_point`/`extract_euclidean` | `Motor`,`conformal_point`,`extract_euclidean` | UT 6 ✔ · INT F81 — PASS | **clifford + scipy Rotation ✔ F81** |
| `src/cga/conformer.rs` | 323 | Conformeros CGA: lift/extract, `dihedral_motor`, rotación de subárbol, `refine_torsion_cga` | `dihedral_motor`,`apply_motor_to_subtree`,`refine_torsion_cga` | UT 7 ✔ · INT F81 — PASS | **RDKit `GetDihedralRad` ✔ F81 (convención espejo documentada)** |
| `src/cga/materials.rs` | 300 | Frame cristalográfico CGA, frac↔cart, supercel, ensamblaje **solo PCU** | `CgaFrame`,`CgaCrystalStructure`,`assemble_framework_cga` | UT 6 ✔ | **pormake/CoRE MOF — ✘ sin paridad real (doc admite gap)** |
| `src/sdr/mod.rs` | 15 | Fachada SDR | re-exports | INT F83 — PASS | — |
| `src/sdr/projections.rs` | 301 | Operadores de relajación convexa: `project_psd`, `project_distances`, `svt_step`, proyecciones alternadas | `SdrConfig`,`project_psd`,`project_distances`,`svt_step` | UT 5 ✔ · INT F83 — PASS | **numpy eigh-clip + SVT ✔ F83** · `project_distances` heurístico |
| `src/sdr/embedding.rs` | 213 | Pipeline SDR: `warm_start_gram` (double-centering), `extract_coordinates` (top-3 eigenpairs), `sdr_embed` | `SdrResult`,`sdr_embed`,`warm_start_gram` | UT 4 ✔ · INT F83 — PASS | **CVXPY/SCS ✔ F83** · `warm_start_gram` O(n³) con medias recomputadas |
| `src/riemannian/mod.rs` | 13 | Fachada | re-exports | INT F94 — PASS | — |
| `src/riemannian/manifold.rs` | 275 | Cono PSD: `PsdManifold`, `psd_projection`, `psd_retraction`, `tangent_projection` (solo simetría), `psd_distance` (proxy Frobenius) | `PsdManifold`,`psd_projection`,`psd_retraction` | UT 6 ✔ · INT F94 — PASS | **numpy ✔ F94** · **geomstats/pymanopt — ✘ sin paridad; `psd_distance` no es la métrica anunciada** |
| `src/riemannian/lbfgs.rs` | 534 | L-BFGS Riemanniano: two-loop, γ=⟨s,y⟩/⟨y,y⟩, proyección tangencial, Armijo a lo largo de la retracción | `RiemannianLbfgs`,`RiemannianConfig` | UT 5 ✔ (FD gradient) · INT F94 — PASS | **scipy L-BFGS-B ✔ F94** · `fallback_enabled`/`used_fallback` muertos |

### 1.7 Dinámica, solvatación, EDL, IM(D)

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/dynamics.rs` | 1912 | MD gas-fase (Velocity Verlet con Berendsen/Nosé-Hoover), NEB simplificado con clamp 0.05 Å (**sin climbing image**), `compute_reaction_dynamics` completo, `atomic_mass_amu` (**12.0 amu por defecto**) | `MdBackend`,`NebBackend`,`MdTrajectory`,`NebPathResult`,`compute_reaction_dynamics`,`atomic_mass_amu` | UT 0 · INT algorithm_validation, reaction_dft_mesh, literature_barriers, reaction_battery — PASS | **LAMMPS/OpenMM/GROMACS/NAMD — ✘ SIN PARIDAD** · `compute_reaction_dynamics` sin ningún test (P1-18) |
| `src/dynamics_live/mod.rs` | 31 | Índice (gate `alpha-dynamics-live`) | re-exports | INT F84 — PASS · **CI no incluye el gate (P2-2)** | — |
| `src/dynamics_live/state.rs` | 419 | Estado MD zero-copy: buffers planos, init de velocidades MB, Verlet, termostatos | `LiveMolecularSystem`,`DynamicAtom` | UT 6 ✔ · INT F84 — PASS | **RDKit UFF ✔ F84 (fuerzas, drift NVE≤1e-3, reversibilidad ≤1e-6)** · sin paridad de trayectoria |
| `src/dynamics_live/pbc.rs` | 199 | PBC: `SimulationBox` con wrap/min-image/volume | `SimulationBox` | UT 2 ✔ · INT alpha_dynamics_conservation — PASS | **ASE `get_distances` — ✘ sin paridad; `pymatgen_pbc.json` sin test; min-image inválido para triclínico** |
| `src/dynamics_live/neighbor_list.rs` | 135 | Verlet neighbor list (**O(N²) real**) | `NeighborList` | UT 4 ✔ | LAMMPS `neigh`/ASE — **✘ sin test de completitud** |
| `src/dynamics_live/barostat.rs` | 109 | Berendsen NPT isotropico; **doc y código usan fórmulas distintas**; sin virial en el repo | `BerendsenBarostat`,`compute_pressure_bar` | UT 2 (triviales) | **LAMMPS `fix press/berendsen` — ✘ sin paridad** |
| `src/dynamics_live/langevin.rs` | 180 | Langevin BAOAB con PRNG xorshift64; **unidades eV vs kcal/mol incompatibles (P0-7)** | `LangevinConfig`,`langevin_baoa_step` | UT 1 (solo "velocidades no nulas") | **OpenMM/LAMMPS/GROMACS — ✘ sin paridad; cuantitativamente incorrecto** |
| `src/imd.rs` | 119 | Steering interactivo: `SteeringForces` acumula energía armónica y fuerzas | `SteeringForce`,`SteeringForces` | UT 3 ✔ · INT F85 — PASS | **F85 ✔ vs oráculo analítico numpy (1e-12/1e-9)** · sin sesión VMD/NAMD real |
| `src/solvation.rs` | 312 | GB flavour: ASP de Wesson-Eisenberg + SASA, radios Born Still-HCT, ecuación de Still con f_GB | `NonPolarSolvation`,`GbSolvation`,`compute_nonpolar_solvation`,`compute_born_radii`,`compute_gb_solvation` | UT 4 ✔ (solo signos) · INT ci, experimental/test_benchmarks — PASS | **AmberTools GB / OpenMM / xTB ALPB — ✘ SIN PARIDAD** · integral de descreening ≠ forma HCT publicada · ASP solo 11 elementos |
| `src/solvation_alpb/mod.rs` | 31 | Fachada ALPB | re-exports | INT test_alpb (8) — PASS | — |
| `src/solvation_alpb/born.rs` | 147 | Radios Born HCT + OBC-II (tanh 1/0.8/4.85) + `gb_kernel` | `AlpbBornRadii`,`compute_born_radii`,`gb_kernel` | UT 3 ✔ · INT test_alpb, gpu_candidates — PASS | **xTB/ALPB, AmberTools — ✘ sin paridad** · `ρ = r + 0.1·probe` ad-hoc |
| `src/solvation_alpb/solvation.rs` | 215 | Energía ALPB: Still-HCT + factor Klamt `(ε−1)/(ε+0.571412x)` + no polar γ·SASA | `AlpbConfig`,`AlpbResult`,`compute_alpb_solvation` | UT 5 ✔ · INT test_alpb, test_benchmarks — PASS | **xTB ALPB (`xtb --alpb`) — ✘ SIN PARIDAD pese a declararse "xTB-compatible"** · no gradientes |
| `src/edl/mod.rs` | 1554 | EDL: Gouy-Chapman, GCS, Helmholtz, ALPB-Stern, adaptadores EEQ/CPM, diagnósticos de convergencia | `EdlConfig`,`EdlProfileResult`,`compute_edl_profile`,`scan_edl_capacitance`,`couple_cpm_eeq_edl` | UT 20 ✔ · INT F86 — PASS | **F86 ✔ vs oráculo analítico numpy (2e-14)** · sin paridad vs solver PN (FEM) |
| `src/edl/gpu.rs` | 108 | Wrappers GPU de barridos EDL — **stub CPU sin kernel WGSL** | `scan_edl_capacitance_gpu`,`compute_edl_profiles_gpu` | UT 2 (CPU vs CPU) | — **✘** (P2-32) |

### 1.8 Alpha, beta, GPU, transport, misceláneo

| File | LOC | Qué hace | API pública | Testing | Comparación externa |
|---|---|---|---|---|---|
| `src/alpha/mod.rs` | 53 | Namespace alpha con 8 subsistemas tras gate | `pub mod cga,edl,gsm,kinetics,periodic_linear,render_bridge,sdr,reaction_dynamics` | INT varios — PASS | — |
| `src/alpha/cga.rs` | 8 | Shim → `crate::cga` | re-exports | INT F81 (vía core) | — |
| `src/alpha/edl.rs` | 8 | Shim → `crate::edl` | re-exports | INT F86 (vía core) | — |
| `src/alpha/sdr.rs` | 8 | Shim → `crate::sdr` | re-exports | INT F83 (vía core) | — |
| `src/alpha/gsm/mod.rs` | 18 | Fachada GSM | re-exports | INT F82 — PASS | — |
| `src/alpha/gsm/backend.rs` | 724 | Router de backends de energía (UFF/MMFF94/EHT/PM3/xTB/GFN1/GFN2/HF-3c/ReaxFF) + planificación por elemento | `GsmEnergyBackend`,`plan_gsm_backends`,`compare_gsm_backends` | UT 5 ✔ · INT F82 (no la matriz) | **ASE NEB/FIRE + barreras publicadas — ✘ sin paridad de barreras PM3/xTB/GFN2/HF-3c** |
| `src/alpha/gsm/saddle.rs` | 248 | Búsqueda de saddle: `refine_saddle` minimiza el gradiente perpendicular | `GsmResult`,`find_transition_state`,`refine_saddle` | UT 3 ✔ · INT F82 — PASS | **ASE `ase.mep.NEB` ✔ F82 (barrera 1.6e-3)** · sin check de modo imaginario |
| `src/alpha/gsm/string.rs` | 350 | Crecimiento de string bidireccional + reparametrización por longitud de arco | `GsmConfig`,`GsmPath`,`grow_string`,`interpolate_node` | UT 4 ✔ · INT F82 — PASS | pyGSM/ASE — **✘ parcial** |
| `src/alpha/kinetics/mod.rs` | 1076 | HTST/Eyring, K_eq, integrador microcinético Euler explícito con sub-stepping adaptativo, steady-state, cadena GSM→MBH→HTST | `evaluate_htst_rate`,`solve_microkinetic_network`,`analyze_gsm_mbh_htst_step*` | UT 14 ✔ · INT F87 — PASS | **Eyring + scipy `solve_ivp` ✔ F87 (1e-12/1e-4)** · Cantera sin paridad |
| `src/alpha/kinetics/gpu.rs` | 59 | Wrapper GPU de barrido HTST — **stub CPU** | `evaluate_htst_sweep_gpu` | UT 1 (CPU vs CPU) | — **✘** |
| `src/alpha/periodic_linear/mod.rs` | 195 | Adaptadores periódicos sobre KPM/RandNLA | `compute_periodic_kpm_dos_from_operators`,`solve_periodic_randnla` | UT 2 ✔ | — **✘ sin paridad de los adaptadores** |
| `src/alpha/periodic_linear/gpu.rs` | 59 | Wrapper GPU de H(k) — **stub CPU** | `assemble_periodic_operators_gpu` | UT 1 | — **✘** |
| `src/alpha/reaction_dynamics/mod.rs` | 1105 | Pipeline 3D: embed de fragmentos, sitios reactivos, complejo guiado por producto, mapeo greedy, CI-NEB, frames restringidos | `compute_reaction_dynamics_3d`,`compute_approach_dynamics` | UT 0 · INT 35/35 — PASS (auto-consistencia) | **ASE NEB — ✘ sin paridad (criterio ±1 kcal/mol no implementado)** |
| `src/alpha/reaction_dynamics/ci_neb.rs` | 468 | Climbing-Image NEB con tangente mejorada, FIRE, IDPP | `CiNebConfig`,`CiNebResult`,`compute_ci_neb_path` | INT 3 — PASS | **ASE `ase.mep.NEB` — ✘ sin paridad** |
| `src/alpha/reaction_dynamics/complex.rs` | 654 | Complejo 3D sin forzar eje X, guiado por producto con Kabsch, búsqueda global con `MinstdRand` | `build_3d_reaction_complex_guided`,`optimize_reactive_complex_global` | UT 2 ✔ · INT 3 — PASS | RDKit `AlignMol` — **✘ sin paridad** |
| `src/alpha/reaction_dynamics/constrained_opt.rs` | 591 | L-BFGS + SHAKE + proyección de gradiente en constraints | `optimize_with_constraints`,`optimize_frame_sequence` | INT 3 — PASS | ASE `ase.constraints` / PySCF `pybfgs` — **✘ sin paridad** |
| `src/alpha/reaction_dynamics/electrostatics.rs` | 164 | Dirección electrostática con cargas Gasteiger por fragmento | `compute_electrostatic_approach` | INT 1 — PASS | RDKit Gasteiger — **✘ sin paridad** |
| `src/alpha/reaction_dynamics/idpp.rs` | 151 | IDPP (Smidstrup 2014) con pesos 1/d⁴ | `generate_idpp_path` | INT 2 — PASS | ASE `ase.mep.idpp` — **✘ sin paridad** |
| `src/alpha/reaction_dynamics/irc.rs` | 343 | IRC: Hessiana FD, modo TS por power iteration, descenso steepest mass-weighted | `IrcResult`,`compute_irc` | INT 1 — PASS | pysisyphus/Gaussian IRC — **✘ sin paridad; no es un IRC real** |
| `src/alpha/reaction_dynamics/orbital_approach.rs` | 215 | Dirección por FMO: EHT por fragmento, HOMO/LUMO, centroide orbital | `compute_orbital_approach_direction` | INT 1 — PASS | FMO/DFT — **✘ sin paridad; base→átomo heurística** |
| `src/alpha/reaction_dynamics/per_frame.rs` | 154 | Propiedades por frame: EHT, Mulliken, dipolo, órdenes de enlace "Wiberg" (P²) | `compute_frame_properties` | INT 2 — PASS | PySCF — **✘ sin paridad** |
| `src/alpha/reaction_dynamics/reactive_sites.rs` | 513 | Sitios reactivos: SMIRKS (0.95) > Fukui (0.8) > Gasteiger (0.4) > pesado (0.2) | `identify_reactive_sites`,`reactive_atom_pairs` | INT 5 — PASS | RDKit — **✘ parcial (solo contra el propio matcher SMIRKS)** |
| `src/alpha/reaction_dynamics/sampling.rs` | 264 | Muestreo multi-angular Fibonacci con "probe" NEB de 5 imágenes | `sample_approach_orientations` | INT 1 — PASS | ASE multi-start / CREST — **✘ sin paridad; doc dice SO(3) pero muestrea S²** |
| `src/alpha/render_bridge/mod.rs` | 593 | Constructores de charts de dominio (EDL, DOS periódico, microcinética, capacitancia, Arrhenius, bandas) + packers | `edl_profile_chart`,`arrhenius_chart`,`band_structure_chart`,`build_full_render_payload` | UT 10 ✔ | JSON Schema/Vega-Lite — **✔ parcial (F90 cubre el core, no los builders)** |
| `src/alpha/render_bridge/gpu.rs` | 55 | Wrapper GPU de packing — **stub CPU** | `pack_charts_gpu` | UT 1 | — **✘** |
| `src/beta/mod.rs` | 49 | Namespace beta | `pub mod kpm,mbh,cpm,rand_nla,riemannian` | INT beta_* — PASS | — |
| `src/beta/kpm/mod.rs` | 14 | Fachada KPM | re-exports | INT F91 — PASS | — |
| `src/beta/kpm/chebyshev.rs` | 393 | Expansión de Chebyshev: bounds Gershgorin, rescale, kernel Jackson, momentos exactos + estimates estocásticos, DOS | `ChebyshevExpansion`,`jackson_kernel`,`rescale_matrix` | UT 5 ✔ · INT F91 — PASS | **numpy eigh ✔ F91 (3.5e-15)** |
| `src/beta/kpm/density.rs` | 386 | DOS KPM: coeficientes Fermi-Dirac por Gauss-Chebyshev, Löwdin, bisección de μ, Mulliken | `KpmConfig`,`compute_kpm_dos`,`compute_kpm_mulliken` | UT 4 ✔ · INT F91 — PASS | **numpy eigh ✔ F91** · criteria de promoción no cumplidos |
| `src/beta/mbh/mod.rs` | 11 | Fachada MBH | re-exports | INT F92 — PASS | — |
| `src/beta/mbh/blocks.rs` | 306 | Detección de bloques rígidos desde SSSR + matriz de proyección L | `detect_rigid_blocks`,`build_projection_matrix` | UT 4 ✔ · INT F92 — PASS | MBH literatura — **✘ parcial (solo estructura de bloques)** |
| `src/beta/mbh/hessian.rs` | 324 | Hessiana reducida 4-puntos, M_red=LᵀML, eigen en forma estándar, cm⁻¹ | `MbhConfig`,`compute_mbh_frequencies` | UT 3 ✔ · INT F92 — PASS | **Hessian cartesiana completa L0 ✔ parcial** · ±2 cm⁻¹ experimental no cumplido |
| `src/beta/rand_nla/mod.rs` | 10 | Fachada RandNLA | re-exports | INT F93 — PASS | — |
| `src/beta/rand_nla/nystrom.rs` | 283 | Gaussian sketch Ω, aproximación Nyström, S^{-1/2} proxy | `GaussianSketch`,`NystromApprox` | UT 6 ✔ · INT F93 — PASS | numpy/scipy ✔ F93 · sin `scikit-learn` randomized-eigensolver |
| `src/beta/rand_nla/solver.rs` | 320 | EHT randomizado: Nyström → H' → eig → back-transform, residual a posteriori, fallback exacto | `RandNlaConfig`,`solve_eht_randnla` | UT 4 ✔ · INT F93 — PASS | **numpy/scipy ✔ F93** · criterio "<1% en rank N/4" no cumplido (~0.13) |
| `src/beta/riemannian.rs` | 11 | Shim → `crate::riemannian` | re-exports | INT F94 (vía core) | — |
| `src/gpu/mod.rs` | 56 | Árbol GPU de 26 submódulos | `pub mod …` | INT gpu_candidates (40) — PASS | — |
| `src/gpu/context.rs` | 634 | Contexto wgpu: instance/adapter/device/queue, `run_compute` con staging readback, `best_available` | `GpuContext`,`ComputeCapabilities`,`run_compute` | UT 7 ✔ · **FLAKY: panic EGL en headless (P0-1)** | WebGPU/wgpu 0.19 — **✘ sin paridad; contrato de fallback violado** |
| `src/gpu/aligned_types.rs` | 172 | PODs `repr(C)` alineados para VRAM | `GpuAtom`,`GpuBasisFunction`,`GpuDenseMatrix` | UT 5 ✔ | WGSL layout — **✘ sin paridad** |
| `src/gpu/backend_report.rs` | 101 | Telemetría de backend GPU | `GpuActivationReport`,`OrbitalGridReport` | UT 3 ✔ | — |
| `src/gpu/buffer_manager.rs` | 309 | Ciclo de vida de buffers + plan SCF de 11 buffers + límites 128 MB | `BufferManager`,`ScfBufferPlan` | UT 5 ✔ | WebGPU limits — **✘ sin paridad con dispositivo real** |
| `src/gpu/memory_budget.rs` | 463 | Límites W3C/WebGPU + estimadores de memoria + dispatch óptimo | `GpuMemoryBudget`,`MemoryError` | UT 11 ✔ | Spec WebGPU — codificado, sin paridad con hardware |
| `src/gpu/shader_registry.rs` | 350 | Catálogo de 13 shaders WGSL con tiers | `ShaderDescriptor`,`ALL_SHADERS` | UT 8 ✔ | — **4 shaders existentes sin registrar (P2-26)** |
| `src/gpu/pipeline_coordinator.rs` | 421 | Scheduler multi-kernel con decisiones GPU/CPU | `PipelineCoordinator`,`plan_pipeline` | UT 11 ✔ | — **✘ no cableado a ningún path de producción** |
| `src/gpu/orbital_grid.rs` | 574 | MO y densidad en grid 3D (CPU + packing GPU) | `GridParams`,`evaluate_orbital_with_report` | UT 7 ✔ | Cubegen/Gaussian cube — **✘ sin paridad** |
| `src/gpu/density_grid_gpu.rs` | 233 | Densidad electrónica en grid por WGSL | `evaluate_density_with_report` | UT 0 · INT GPU vs CPU — PASS | PySCF cubegen — **✘ sin paridad** |
| `src/gpu/esp_grid_gpu.rs` | 216 | ESP point-charge en grid por WGSL | `compute_esp_grid_with_report` | UT 0 · INT GPU vs CPU — PASS | Psi4/CP2K — **✘ sin paridad** |
| `src/gpu/marching_cubes.rs` | 697 | Marching cubes CPU + GPU con contador atómico | `McOutput`,`marching_cubes_cpu`,`marching_cubes_gpu` | UT 4 ✔ | VTK/skimage — **✘ sin paridad; mesh GPU nunca comparado con CPU** |
| `src/gpu/isosurface.rs` | 228 | Pipeline de isosuperficie orbital (HOMO/LUMO) | `generate_orbital_isosurface` | UT 4 ✔ · **FLAKY (P0-1)** | VMD/PyMOL — **✘ sin paridad** |
| `src/gpu/two_electron_gpu.rs` | 515 | ERIs O(N⁴) en GPU con simetría 8-fold y tiling | `compute_eris_gpu`,`eri_quartet_tiles` | UT 6 ✔ (round-trip CPU 1e-12) | libcint — **✘ referencia solo interna** |
| `src/gpu/one_electron_gpu.rs` | 415 | S/T/V en GPU (solo s-type, fallback CPU si hay p/d) | `compute_one_electron_gpu` | UT 3 ✔ | libcint — **✘ sin paridad** |
| `src/gpu/fock_build_gpu.rs` | 202 | Fock en GPU | `build_fock_gpu` | UT 3 ✔ | PySCF — **✘ sin paridad** |
| `src/gpu/gamma_matrix_gpu.rs` | 136 | γ-matrix SCC-DFTB en GPU | `build_gamma_gpu` | UT 1 ✔ | tblite — **✘ sin paridad; módulo no cableado** |
| `src/gpu/eeq_gpu.rs` | 213 | EEQ: γ con erf en GPU + solve lineal CPU | `compute_eeq_charges_gpu` | UT 1 ✔ | multicharge — **✘ sin paridad GPU vs CPU** |
| `src/gpu/cpm_gpu.rs` | 170 | CPM: J en GPU + update SCF damping 0.5 en CPU | `compute_cpm_charges_gpu` | UT 0 | LAMMPS `electrode/conp` — **✘ sin paridad** |
| `src/gpu/d4_dispersion_gpu.rs` | 178 | D4 dos cuerpos en GPU (e3 en CPU) | `compute_d4_energy_gpu` | UT 0 | dftd4 — **✘ sin paridad** |
| `src/gpu/alpb_born_gpu.rs` | 172 | Radios Born ALPB en GPU | `compute_born_radii_gpu` | UT 1 ✔ | Amber ALPB — **✘ sin paridad con CPU** |
| `src/gpu/hct_born_gpu.rs` | 192 | Radios Born HCT en GPU | `compute_hct_born_radii_gpu` | UT 1 ✔ | HCT — **✘ sin paridad con CPU** |
| `src/gpu/ani_aev_gpu.rs` | 189 | AEV radial en GPU | `compute_radial_aev_gpu` | UT 1 ✔ | nANI/ANI-2x — **✘ sin paridad con `ani/aev.rs`** |
| `src/gpu/mmff94_gpu.rs` | 555 | MMFF94 no-bonded GPU (buffered 14-7 + DDD) | `compute_mmff94_nonbonded_gpu` | UT 4 (tiled CPU vs untiled CPU) | OpenMM/Tinker — **✘ P0-4: `atomicAdd` sobre floats bitcast es inválido** |
| `src/gpu/stda_gpu.rs` | 235 | Matriz de Coulomb sTDA en GPU (2 matmuls) | `compute_stda_j_matrix_gpu` | UT 0 | Turbomole sTDA — **✘ sin paridad; módulo no cableado** |
| `src/gpu/hessian_gpu.rs` | 149 | Batch de desplazamientos para Hessiana numérica (**sin kernel WGSL pese al nombre**) | `generate_hessian_displacements`,`assemble_hessian_from_energies` | UT 2 ✔ | — **✘ off-diagonal aproximado sin test vs 4-puntos** |
| `src/render_bridge.rs` | 134 | Contratos de payload de charts (serde JSON) + packing a `RecordBatch` + round-trip | `ChartPayload`,`RenderPayload`,`pack_chart_payload`,`verify_chart_round_trip` | UT 0 · INT F90 — PASS | **JSON Schema 2020-12 ✔ F90** · sin IPC Arrow real |
| `src/transport/mod.rs` | 14 | Fachada transport | `pub mod arrow,chunked,worker` | — | — |
| `src/transport/arrow.rs` | 307 | Layout columnar en memoria (**no es Arrow IPC**), packers conformers/ESP/DOS | `RecordBatch`,`ColumnSchema`,`pack_conformers` | UT 5 ✔ · INT ninguna | **Apache Arrow — ✘ SIN PARIDAD (sin dependencia `arrow`)** |
| `src/transport/chunked.rs` | 221 | Streaming por chunks de arrays grandes | `ChunkedIterator`,`DataChunk`,`chunk_*` | UT 6 ✔ | NDJSON/Arrow IPC — **✘ sin test de formato de cable** |
| `src/transport/worker.rs` | 183 | Despacho a Web Workers: split/merge/estimate | `WorkerTask`,`split_batch`,`merge_results`,`estimate_workers` | UT 6 ✔ | Web Worker API — **✘ sin consumidor JS** |
| `src/experimental_status.rs` | 194 | Módulo **solo-documentación**: criterios de promoción por módulo | (ninguna) | UT 0 | — **Stale (P2-25)** |
| `src/fixture_io.rs` | 56 | Resolución de fixtures con soporte `.gz` | `resolve_fixture_path`,`read_text_fixture` | UT 0 · INT varios — PASS | — **✘ sin test de la rama .gz** |
| `src/bin/*.rs` (27 binarios) | 3060 | Herramientas dev/debug/bench: `benchmark`, `benchmark_10k`, `full_benchmark`, `fast_bench`, `eht_metal_reference`, `check_grad`, `check_c_ints`, `print_c_eri`, `debug_*` (11), `dump_*` (3), `diagnose_failures`, `find_attempt`, `rmsd_per_mol`, `test_eigen`, `test_failures`, `trace_pipeline`, `verify_mt` (KAT de MT19937) | `fn main()` (sin API reutilizable) | 0 tests; no corren en `cargo test` | Fixtures de referencia RDKit existen como **datos** pero **sin test automatizado** |

---

## 2. TABLA — BINDINGS (Python / WASM / CLI)

### 2.1 Python (`crates/python/src/`, PyO3 0.22)

| File | LOC | Qué hace | API expuesta | Testing | Comparación externa |
|---|---|---|---|---|---|
| `crates/python/src/lib.rs` | 69 | Módulo raíz; registra 26 submódulos | `sci_form` | compila solo con `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1` en Python 3.14 | RDKit `rdkit.Chem` — sci-form no tiene modelo `Mol`; todo son funciones libres |
| `crates/python/src/embed.rs` | 116 | `embed`/`embed_batch`/`parse`/`version` → `ConformerResult` | `ConformerResult`,`embed`,`embed_batch`,`parse` | script manual, **no en CI** | RDKit: `Mol` + numpy (N,3). sci-form: lista plana + `get_positions()` |
| `crates/python/src/system.rs` | 186 | capabilities / method plan | `SystemCapabilitiesPy`,`system_capabilities`,`system_method_plan` | ninguna | sin análogo en el ecosistema |
| `crates/python/src/compare.rs` | 129 | comparación multi-método | `MethodComparisonResultPy`,`compare_methods` | ninguna | sin análogo |
| `crates/python/src/eht.rs` | 222 | EHT + fallback + support + marching cubes | `EhtResultPy`,`eht_calculate`,`eht_or_uff_fallback` | ninguna | **YAeHMOP/PySCF devuelven numpy; sci-form no expone `.mo_energy` ni eigenvectores en Python** |
| `crates/python/src/properties.rs` | 202 | charges, sasa, dipole, dos, rmsd | 5 pyclasses | ninguna | RDKit devuelve numpy; sci-form pyclass sin numpy |
| `crates/python/src/population.rs` | 156 | población, órdenes de enlace, frontier | 3 pyclasses | ninguna | PySCF `mf.mulliken_pop()` → numpy |
| `crates/python/src/reactivity.rs` | 135 | Fukui + ranking | 3 pyclasses | ninguna | Multiwfn: tabla por átomo; sci-form pierde índices |
| `crates/python/src/analysis.rs` | 171 | uvvis, graph features, topología | 5 pyclasses | ninguna | PySCF TD-DFT → numpy |
| `crates/python/src/forcefield.rs` | 103 | UFF, UFF+heurísticas, pKa, MMFF94 | 4 pyclasses | ninguna | RDKit devuelve float; sci-form float **sin descomp. por término** |
| `crates/python/src/materials.rs` | 122 | unit cell, ensamblaje de framework | `UnitCellPy`,`CrystalStructurePy` | ninguna | pymatgen/ASE: sin interop numpy |
| `crates/python/src/transport.rs` | 104 | pack conformers, split workers | `RecordBatchPy` | ninguna | **PyArrow: no hay IPC real; "smiles" es un string JSON, no lista** |
| `crates/python/src/electronic.rs` | 252 | PM3, xTB, GFN1, GFN2 | 4 pyclasses | ninguna | **xtb-python devuelve `.get_energy()` en eV; sci-form no etiqueta unidades** |
| `crates/python/src/ml.rs` | 119 | descriptores + predicción | 2 pyclasses | ninguna | mordred → DataFrame; sci-form heurístico sin flag de confianza |
| `crates/python/src/spectroscopy.rs` | 339 | sTDA, fluorescencia, vibracional, IR | 7 pyclasses | ninguna | ORCA/PySCF vibspectrum; **sci-form pierde `thermochemistry` y `label`** |
| `crates/python/src/nmr.rs` | 351 | shifts, GIAO, couplings, spectrum, HOSE | 8 pyclasses | ninguna | nmrshiftdb2; **`hose_codes` devuelve tuplas, no `HoseCodePy` (doc miente)** |
| `crates/python/src/hf_ani_esp.rs` | 160 | HF-3c, ANI, ESP | 3 pyclasses | ninguna | **PySCF siempre Hartree; sci-form no etiqueta unidad** |
| `crates/python/src/mesh.rs` | 78 | orbital mesh | `OrbitalMeshResultPy` | ninguna | mesh f32 vs expectativa f64 |
| `crates/python/src/dynamics.rs` | 243 | MD, NEB, conformer search | 7 pyclasses | ninguna | MDAnalysis/OpenMM: sin Timestep ni acceso a átomos |
| `crates/python/src/stereo.rs` | 88 | estereoquímica | 3 pyclasses | ninguna | RDKit `FindMolChiralCenters`; sci-form no devuelve centros no asignados |
| `crates/python/src/solvation.rs` | 92 | no polar, GB | 2 pyclasses | ninguna | OpenMM `GBSAOBCForce`; nomenclatura de ε divergente |
| `crates/python/src/rings.rs` | 164 | SSSR, ECFP, Tanimoto, Butina, RMSD matrix, filter_diverse | 4 pyclasses | solo lado Rust (F05/F59/F14) | **RDKit devuelve `ExplicitBitVect`; sci-form devuelve índices de bits → no comparables** |
| `crates/python/src/smirks.rs` | 132 | parse/apply SMIRKS | `SmirksTransform`,`SmirksResult` | Rust: 15 tests | **RDKit devuelve `Product` molecules; sci-form devuelve strings SMARTS** |
| `crates/python/src/potentials.rs` | 52 | potenciales analíticos (submódulo `sci_form.potentials`) | 6 fns | ninguna | no hay análogo; **no existe en WASM ni CLI** |
| `crates/python/src/experimental.rs` | 1369 | 8 submódulos tras gate; **~26 fns devuelven constantes hard-coded** | ~45 fns | ninguna | **API falsa duplicada en WASM** (P0-3) |
| `crates/python/src/alpha.rs` | 1011 | `sci_form.alpha`: DFT, ReaxFF, MLFF, OS, CGA, GSM, SDR, reaction-dynamics | ~15 fns | ninguna (gates F77/F82 en Rust, no en CI) | **nombres divergentes vs WASM** (`sdr_embed_py` vs `alpha_sdr_embed`) |
| `crates/python/src/beta.rs` | 446 | `sci_form.beta`: KPM, MBH, RandNLA, Riemannian, CPM | ~6 fns | ninguna | `beta_psd_distance` devuelve float sin wrapper de matriz → no componible con geomstats |

### 2.2 WASM / TypeScript (`crates/wasm/`)

| File | LOC | Qué hace | API expuesta | Testing | Comparación externa |
|---|---|---|---|---|---|
| `crates/wasm/src/lib.rs` | 37 | Árbol de módulos; `init_thread_pool`; `webgpu` y `dynamics_live` tras gate | 22 módulos | ninguna | — |
| `crates/wasm/src/helpers.rs` | 38 | Parsers JSON compartidos | interno | — | — |
| `crates/wasm/src/embed.rs` | 78 | embed (JSON), `embed_coords_typed` (Float64Array), batch NDJSON, parse | 5 fns | script JS roto, no en CI | **único path zero-copy, estilo molstar** |
| `crates/wasm/src/system.rs` | 50 | version, capabilities, plan, compare | 4 fns | ninguna | — |
| `crates/wasm/src/properties.rs` | 128 | charges, sasa, population, bond orders, dipole, dos, rmsd | 7 fns | ninguna | **`compute_rmsd` sin `aligned_coords` → paridad L1 imposible** |
| `crates/wasm/src/forcefield.rs` | 52 | UFF, UFF+heur, MMFF94, pKa | 4 fns | ninguna | **devuelve `{energy,unit}` vs float en Python** |
| `crates/wasm/src/eht.rs` | 1013 | EHT + WebGPU (~600 líneas WGSL) | 13 fns | harness inexistente | **typed arrays OK; mesh JSON no es interleaved (≠ glTF)** |
| `crates/wasm/src/electronic.rs` | 137 | PM3, xTB, GFN1/2, HF-3c, ANI, DFT | 8 fns | ninguna | energía sin etiqueta de unidad |
| `crates/wasm/src/mesh.rs` | 188 | orbital mesh (+accelerated) | 2 fns | ninguna | **el grid acelerado usa siempre base EHT** |
| `crates/wasm/src/esp.rs` | 383 | ESP JSON + typed + acelerado (CPU/GPU/hybrid) | 6 fns | harness inexistente | dims OK |
| `crates/wasm/src/spectroscopy.rs` | 246 | sTDA, fluorescencia, vibracional, IR, NMR | 12 fns | ninguna | **le falta `compute_ir_spectrum_broadened` y `compute_ensemble_j_couplings` que el doc promete** |
| `crates/wasm/src/reactivity.rs` | 83 | frontier, Fukui, ranking, uvvis, graph features, topología | 6 fns | ninguna | — |
| `crates/wasm/src/ml.rs` | 68 | ML properties + descriptores | 2 fns | ninguna | — |
| `crates/wasm/src/rings.rs` | 71 | SSSR, ECFP, Tanimoto, Butina, RMSD matrix | 5 fns | ninguna | **falta `filter_diverse`** |
| `crates/wasm/src/solvation.rs` | 49 | no polar, GB | 2 fns | ninguna | — |
| `crates/wasm/src/stereo.rs` | 20 | estereoquímica | 1 fn | ninguna | — |
| `crates/wasm/src/transport.rs` | 31 | pack arrow, split workers | 3 fns | ninguna | **"pack_batch_arrow" devuelve JSON, no Arrow IPC** |
| `crates/wasm/src/materials.rs` | 54 | unit cell, ensamblaje | 2 fns | ninguna | **`create_unit_cell` sin `lattice` → menos datos que Rust/Python/CLI** |
| `crates/wasm/src/dynamics.rs` | 113 | MD, NEB, conformer search | 4 fns | ninguna | — |
| `crates/wasm/src/dynamics_live.rs` | 194 | Clase `LiveSimulation` con `positions_ptr()` zero-copy | `LiveSimulation` | ninguna | **la parte más fuerte de la superficie WASM** |
| `crates/wasm/src/reaction.rs` | 284 | SMIRKS, NEB por método, reaction dynamics | 8 fns | ninguna | productos como SMARTS (misma divergencia) |
| `crates/wasm/src/experimental.rs` | 1643 | espejo de Python con **~26 fakes** | ~40 fns | ninguna | **P0-3** |
| `crates/wasm/src/alpha.rs` | 1060 | superficie alpha (~30 fns) | `alpha_*` | ninguna | `alpha_modules_info` reporta versión 0.14.4 |
| `crates/wasm/src/beta.rs` | 376 | superficie beta (8 fns) | `beta_*` | ninguna | `beta_modules_info` reporta 0.11.2; `config_json` vs kwargs de Python |
| `crates/wasm/src/webgpu.rs` | 322 | Runtime WebGPU: init, status, dispatch async | `init_webgpu`,`webgpu_status` | harness inexistente | fallback silencioso a CPU |
| `crates/wasm/js/alpha/index.js` | 70 | Shim ESM del subpath `/alpha` | 15 re-exports | solo check de existencia en CI | convención npm ✓ |
| `crates/wasm/js/alpha/index.d.ts` | 204 | Declaraciones TS del subpath alpha | 15 decls | ninguna | **DERIVA de firmas Rust (P2-13)** |
| `crates/wasm/js/beta/index.js` | 46 | Shim ESM `/beta` | 8 re-exports | check de existencia | ✓ |
| `crates/wasm/js/beta/index.d.ts` | 119 | Declaraciones TS `/beta` | 8 decls | ninguna | **consistente con Rust ✓ (único archivo consistente)** |
| `crates/wasm/ts/live_simulation.ts` | 82 | `LiveSimulationView` con Float64Array zero-copy | `LiveSimulationView` | **sin tsconfig, sin type-check en CI** | correcto en principio, nunca compilado |
| `crates/wasm/ts/simulation_loop.ts` | 116 | `SimulationLoop` con rAF | `SimulationLoop` | ninguna | idem |
| `crates/wasm/Cargo.toml` | 64 | Manifest; `default=[]`, 25 features, wgpu 29 (core usa 0.19) | — | CI lo construye | **skew de versión (P2-21)** |
| `crates/wasm/.cargo/config.toml` | 5 | `+atomics,+bulk-memory,+mutable-globals` + `build-std` (requiere nightly; CI usa stable → inerte) | — | — | — |

### 2.3 CLI (`crates/cli/src/`)

| File | LOC | Qué hace | API expuesta | Testing | Comparación externa |
|---|---|---|---|---|---|
| `crates/cli/src/main.rs` | 351 | Dispatch clap de **68 subcomandos**; `ReactionDynamics3d` inlineado (no pasa por `experimental_cmds`) | 68 subcomandos | solo `tests/ci.rs` (lado Rust) + 3 greps en release.yml | obabel/xtb/ase: sin flag global `--json`/`--quiet` |
| `crates/cli/src/args.rs` | 554 | Enum `Commands` completo + flags; gates por feature para alpb/d4/cpm/charts/periodic/gsm/kinetics/reaction-dynamics | 68 subcomandos | ninguna | convención clap ✓ |
| `crates/cli/src/embed_cmds.rs` | 251 | `embed` (json/xyz/sdf), `batch` (chunk 512), `parse`, `info`, `charges`, `uff`, `rmsd` | 7 comandos | release.yml greps | **`cmd_info` anuncia "IR, NMR, UV-Vis" sin subcomandos** (P2-16) |
| `crates/cli/src/calc_cmds.rs` | 298 | `eht`, `sasa`, `population`, `dipole`, `esp`, `dos`, `cell`, `assemble`, `ani`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`, `stereo`, `solvation`, `sssr`, `ecfp`, `tanimoto` | 18 comandos | ninguna | **`cmd_esp` descarta los valores (P2-17)**; `cmd_ani` usa red de test (P2-19); `cmd_hf3c` con CISD off (P2-18) |
| `crates/cli/src/experimental_cmds.rs` | 1072 | 40+ comandos: eeq, alpb, d4, cpm, EDL, kmesh, GSM, NEB, charts, kinetics + **~11 stubs** | 7 handlers **inalcanzables** (P2) | ninguna | **falsos duplicados (P0-3)** |
| `crates/cli/src/format.rs` | 95 | `format_xyz`, `format_sdf` (V2000), `element_symbol` (**13 elementos**), `parse_elems_coords` | interno | ninguna | **`element_symbol` → "X" para metales → SDF inválido (P2-15)** |

---

## 3. TABLA — TESTS (89 archivos)

| Test file | Tests | Módulos cubiertos | LOC | Gates externos |
|---|---|---|---|---|
| `tests/cross_validation.rs` | 38 | alignment, charges_eeq, dipole, dispersion, eht, forcefield, hf, ir, materials, nmr, pm3, rings, scf, xtb | 4043 | F01, F05, F13, F14, F16, F17, F21, F22, F26, F28, F31, F36, F37, F38, F39, F41×2, F42, F43, F44, F44-UKS, F45, F45-UKS, F46, F47, F48, F49, F54, F59, F67×3, F68, F70 — **NO corre en CI** |
| `tests/alpha_reaction_dynamics_3d.rs` | 35 | alpha | 1272 | auto-consistencia |
| `tests/regression/test_gpu_candidates.rs` | 40 | alpha, beta, charges_eeq, dispersion, esp, gpu, hf, scf, solvation_alpb, spectroscopy | 2060 | GPU vs CPU (skip sin GPU) |
| `tests/experimental/test_benchmarks.rs` | 44 | alpha, beta, charges_eeq, dispersion, eht, ir, solvation_alpb | 1428 | — |
| `tests/regression/test_algorithm_validation.rs` | 46 | charges, forcefield, hf, ir, nmr, rings | 1341 | — |
| `tests/regression/test_experimental_comparison.rs` | 54 | eht, hf, nmr, pm3, scf, spectroscopy, xtb | 1581 | rangos laxos |
| `tests/regression/test_phase_c_validation.rs` | 37 | alignment, esp | 718 | cualitativo |
| `tests/regression/test_extended_molecules.rs` | 21 | eht, ml, nmr, pm3, scf, spectroscopy, xtb | 1415 | — |
| `tests/regression/test_new_features.rs` | 38 | forcefield, ir, materials, ml, smarts | 877 | — |
| `tests/benchmarks/test_comprehensive_1k_benchmark.rs` | 38 | eht, materials, xtb | 2219 | — |
| `tests/e3_justification.rs` | 28 | charges_eeq, dispersion, solvation_alpb, spectroscopy | 575 | — |
| `tests/alpha_mlff_consistency.rs` | 10 | ml | 238 | F79 parcial |
| `tests/alpha_reaxff_validation.rs` | 13 | forcefield | 392 | F78 |
| `tests/alpha_scf_eri_accuracy.rs` | 12 | scf | 224 | F80 |
| `tests/alpha_dft_validation.rs` | 11 | dft | 331 | F77 |
| `tests/reaction_battery.rs` | 10 | alpha, eht, pm3, xtb | 807 | — |
| `tests/reaction_dynamics_study.rs` | 17 | alpha | 789 | — |
| `tests/alpha_dynamics_conservation.rs` | 8 | dynamics_live | 234 | F84 |
| `tests/experimental/test_riemannian.rs` | 11 | beta | 446 | F94 |
| `tests/experimental/test_cpm.rs` | 9 | cpm | 285 | F95 |
| `tests/experimental/test_kpm.rs` | 9 | beta | 216 | F91 |
| `tests/experimental/test_mbh.rs` | 8 | beta | 141 | F92 |
| `tests/experimental/test_randnla.rs` | 8 | beta | 308 | F93 |
| `tests/experimental/test_sdr.rs` | 8 | alpha | 171 | F83 |
| `tests/experimental/test_cga.rs` | 20 | alpha | 360 | F81 |
| `tests/experimental/test_gsm.rs` | 8 | alpha | 140 | F82 |
| `tests/experimental/test_eeq.rs` | 8 | charges_eeq | 143 | F22 |
| `tests/experimental/test_d4.rs` | 8 | dispersion | 171 | F46 |
| `tests/experimental/test_alpb.rs` | 8 | solvation_alpb | 140 | — |
| `tests/beta_riemannian_validation.rs` | 4 | riemannian | 201 | F94 |
| `tests/beta_kpm_validation.rs` | 3 | beta | 167 | F91 |
| `tests/beta_mbh_validation.rs` | 2 | beta | 134 | F92 |
| `tests/beta_randnla_validation.rs` | 3 | beta | 157 | F93 |
| `tests/alpha_gsm_validation.rs` | 2 | alpha | 114 | F82 |
| `tests/alpha_cga_validation.rs` | 2 | cga | 140 | F81 |
| `tests/alpha_sdr_validation.rs` | 3 | sdr | 138 | F83 |
| `tests/alpha_edl_validation.rs` | 4 | edl | 263 | F86 |
| `tests/alpha_periodic_validation.rs` | 4 | — | 200 | F88 |
| `tests/alpha_kinetics_validation.rs` | 3 | alpha | 174 | F87 |
| `tests/alpha_dynamics_validation.rs` | 3 | dynamics_live | 182 | F84 |
| `tests/alpha_imd_validation.rs` | 3 | — | 121 | F85 |
| `tests/alpha_render_validation.rs` | 2 | — | 144 | F90 |
| `tests/alpha_cross_method.rs` | 6 | dft, eht, forcefield | 170 | — |
| `tests/alpha_dft_pyscf.rs` | 1 | dft | 104 | F77 |
| `tests/alpha_obara_saika.rs` | 3 | scf | 185 | F80 |
| `tests/ci.rs` | 4 | materials | 191 | smoke — **único integration test en CI** |
| `tests/benchmarks/test_new_features_bench.rs` | 10 | forcefield, ir, materials, ml, smarts | 473 | — |
| `tests/benchmarks/test_10k.rs` | 2 | forcefield | 442 | — |
| `tests/benchmarks/test_chembl_10k.rs` | 1 | — | 419 | — |
| `tests/beta_benchmarks.rs` | 5 | beta, ir | 401 | — |
| `tests/cisd_brute.rs` | 2 | hf | 433 | segunda cuantización independiente |
| `tests/test_smirks_reactions.rs` | 15 | smirks | 594 | — |
| `tests/literature_barriers.rs` | 1 | — | 174 | — |
| `tests/check_atom_counts.rs` | 1 | — | 255 | — |
| `tests/regression/test_benchmark_100.rs` | 1 | gpu, hf, pm3, xtb | 1062 | — |
| `tests/regression/test_gdb20_rmsd.rs` | 1 | forcefield | 664 | — |
| `tests/regression/test_geometry_quality.rs` | 1 | forcefield | 435 | — |
| `tests/regression/test_ensemble_rmsd.rs` | 1 | forcefield, smarts | 556 | — |
| `tests/regression/test_energy_comparison.rs` | 1 | distgeom, forcefield | 302 | — |
| `tests/regression/test_gradient_check.rs` | 1 | distgeom, forcefield | 302 | — |
| `tests/regression/test_gradient_perterm.rs` | 1 | distgeom, forcefield | 449 | — |
| `tests/regression/test_grad_isolate.rs` | 1 | distgeom, forcefield | 295 | — |
| `tests/regression/test_grad_torsion_detail.rs` | 1 | distgeom, forcefield | 271 | — |
| `tests/regression/test_pipeline_trace.rs` | 1 | distgeom, forcefield | 300 | — |
| `tests/regression/test_new_pipeline.rs` | 1 | forcefield | 343 | — |
| `tests/regression/test_bounds_compare.rs` | 2 | distgeom | 229 | **`#[ignore]`d — requiere `/tmp/rdkit_bounds_*.json`** |
| `tests/regression/test_embedding_trace.rs` | 1 | distgeom | 185 | — |
| `tests/regression/test_charges.rs` | 6 | — | 79 | — |
| `tests/regression/test_diverse_molecules.rs` | 1 | — | 249 | — |
| `tests/regression/test_eht.rs` | 9 | eht | 244 | — |
| `tests/regression/test_eht_metal_references.rs` | 1 | eht | 205 | — |
| `tests/regression/test_gdb20.rs` | 1 | — | 240 | — |
| `tests/regression/test_metal_experimental_geometry.rs` | 1 | — | 129 | — |
| `tests/regression/test_parallel_public_api.rs` | 1 | materials, xtb | 264 | — |
| `tests/regression/test_phase_c.rs` | 10 | — | 120 | — |
| `tests/regression/test_roadmap_validation.rs` | 1 | — | 185 | — |
| `tests/regression/test_sasa.rs` | 6 | surface | 220 | — |
| `tests/regression/test_spectroscopy.rs` | 34 | spectroscopy | 1180 | — |
| `tests/regression/test_step_by_step.rs` | 1 | distgeom, forcefield | 460 | — |
| `tests/regression/test_tet_centers.rs` | 1 | — | 180 | — |
| `tests/regression/test_reaction_dft_mesh.rs` | 1 | alpha, eht | 610 | — |
| `tests/analysis/compare_bounds.rs` | 1 | distgeom | 126 | requiere `/tmp/rdkit_bounds.npy` |
| `tests/benchmarks/benchmark_volume_vs_mesh.rs` | 1 | eht | 103 | — |
| `tests/debug/debug_4d_coords.rs` | 3 | distgeom, forcefield | 425 | coords 4D; herramienta de diagnóstico |
| `tests/debug/debug_attempts.rs` | 1 | distgeom, forcefield | 133 | volcado de intentos de embedding |
| `tests/debug/debug_bounds.rs` | 1 | distgeom, forcefield | 153 | inspección de bounds |
| `tests/debug/debug_compare_bounds.rs` | 1 | distgeom | 253 | comparación de bounds (requiere fixture RDKit) |
| `tests/debug/debug_differential.rs` | 1 | distgeom | 89 | test diferencial de embedding |
| `tests/debug/debug_embed.rs` | 1 | distgeom | 162 | traza de embedding |
| `tests/debug/debug_failure_diagnosis.rs` | 1 | forcefield | 264 | diagnóstico de fallos de embedding |
| `tests/debug/debug_ff.rs` | 1 | distgeom, forcefield | 930 | mayor archivo de debug: energía/gradiente FF |
| `tests/debug/debug_ff_energy.rs` | 1 | distgeom, forcefield | 316 | energía FF etapa a etapa |
| `tests/debug/debug_remaining.rs` | 1 | — | 180 | casos restantes de debug |
| `tests/debug/differential_testing.rs` | 2 | distgeom, forcefield | 577 | harness de testing diferencial vs referencia |
| `tests/debug/dump_new_metal_refs.rs` | 1 | eht | 125 | volcado de referencias de metales |
| `tests/debug/trace_failing.rs` | 1 | distgeom, forcefield | 468 | traza del pipeline en fallos |
| `tests/debug.rs` | 0 | — | 26 | agregador de módulos de debug (declara submódulos) |
| `tests/analysis.rs` | 0 | — | 2 | agregador de `tests/analysis/` |
| `tests/benchmarks.rs` | 0 | — | 10 | agregador de `tests/benchmarks/` |
| `tests/regression.rs` | 0 | — | 113 | agregador de `tests/regression/` (gates condicionales) |
| `tests/regression/test_session_algorithms.rs` | 35 | session algorithms | 741 | algoritmos de sesión/regresión de búsqueda de conformeros |
| `tests/alpha_mlff_validation.rs` | 2 | ml | 198 | F79: AEV vs Behler-Parrinello de referencia |
| `tests/integration/test_python_integration.py` | 9 asserts | python | — | **no pytest, no en CI** |
| `tests/integration/test_wasm_integration.js` | asserts | wasm | — | **ruta de import rota → MODULE_NOT_FOUND** |
| `tests/integration/test_cross_language.py` | asserts | cross | — | espera `./target/release/sci-form`; nunca invocado |
| `tests/integration/test_webgpu_harness.js` | asserts | webgpu | — | **harness `crates/wasm/examples/vite-webgpu-hybrid/` no existe en git** |
| `tests/integration/compare_webgpu_reports.js` | asserts | webgpu | — | compara reportes GPU vs CPU; **nunca invocado** (el harness que lo produciría no existe) |
| `tests/integration/test_smirks_reactions.py` | asserts | smirks | — | **no es pytest, no en CI** |

---

## 4. MATRIZ DE PARIDAD DE SUPERFICIES

| Capacidad | Rust | Python | WASM | CLI |
|---|---|---|---|---|
| `embed`/`embed_batch`/`parse`/`version` | ✅ | ✅ | ✅ | ✅ |
| `embed_diverse` | ✅ | ❌ | ❌ | ❌ |
| `to_canonical_smiles` | ✅ | ❌ | ❌ | ❌ |
| `compute_rmsd` | ✅ | ✅ | ⚠️ sin `aligned_coords` | ✅ |
| capabilities / method plan / compare | ✅ | ✅ | ✅ | ❌ |
| `compute_charges(_configured)` | ✅ | ⚠️ sin `_configured` | ⚠️ sin `_configured` | ❌ |
| `compute_sasa`/`population`/`dipole`/`esp`/`dos` | ✅ | ✅ | ✅ | ⚠️ esp sin valores |
| `compute_bond_orders`/frontier/Fukui/ranking/uvvis | ✅ | ✅ | ✅ | ❌ |
| `compute_uff_energy`/`mmff94`/aromatic/pKa | ✅ | ✅ | ✅ | ⚠️ falta mmff94 y pKa |
| `compute_pm3`/`xtb`/GFN1/GFN2/`hf3c`/`ani` | ✅ | ✅ | ✅ | ⚠️ hf3c CISD off |
| `compute_cisd` | ✅ | ❌ | ❌ | ❌ |
| `compute_dos_multimethod` | ✅ | ❌ | ❌ | ❌ |
| `compute_eht_gradient`/`optimize_geometry_eht` | ✅ | ❌ | ❌ | ❌ |
| `compute_aevs_tm` | ✅ | ❌ | ❌ | ❌ |
| `compute_whim`/`rdf`/`getaway` | ✅ | ❌ | ❌ | ❌ |
| `train_random_forest`/`gradient_boosting` | ✅ | ❌ | ❌ | ❌ |
| `assign_peaks` | ✅ | ❌ | ❌ | ❌ |
| `compute_ensemble_j_couplings` | ✅ | ❌ (doc dice que sí) | ❌ (doc dice que sí) | ❌ |
| `compute_ir_spectrum_broadened` | ✅ | ❌ (doc dice que sí) | ❌ (doc dice que sí) | ❌ |
| `compute_vibrational_analysis(_uff)` | ✅ | ⚠️ fusionado | ✅ | ❌ |
| NMR (shifts/couplings/spectrum/GIAO/HOSE) | ✅ | ✅ | ✅ | ❌ |
| `compute_nonpolar_solvation`/`gb_solvation` | ✅ | ✅ | ✅ | ✅ |
| `analyze_stereo` | ✅ | ✅ | ✅ | ✅ |
| SSSR/ECFP/Tanimoto/RMSD matrix | ✅ | ✅ | ✅ | ⚠️ falta rmsd_matrix |
| `butina_cluster` | ✅ | ✅ | ✅ | ❌ |
| `filter_diverse` | ✅ | ✅ | ❌ | ❌ |
| `create_unit_cell`/`assemble_framework` | ✅ | ✅ | ⚠️ sin lattice | ✅ |
| space groups / periodic molecule / hapticity / framework opt | ✅ | ❌ | ❌ | ❌ |
| `parse_smirks`/`apply_smirks` | ✅ | ✅ | ✅ | ❌ |
| transport (pack/split/estimate) | ✅ | ✅ | ✅ | ❌ |
| `potentials::*` | ✅ | ✅ | ❌ | ❌ |
| MD / NEB / conformer search | ✅ | ✅ | ✅ | ⚠️ solo NEB |
| `compute_reaction_dynamics` | ✅ | ❌ | ✅ | ❌ |
| Alpha (DFT, ReaxFF, MLFF, OS, CGA, GSM, SDR, EDL, periodic, kinetics, render, reaction) | ✅ | ✅ | ✅ | ⚠️ EDL inalcanzable |
| Beta (KPM, MBH, RandNLA, Riemannian, CPM) | ✅ | ✅ | ✅ | ⚠️ solo CPM |

**Conclusión de paridad:** 18 capacidades son **Rust-only** aunque `api-surface.md` las presenta como cross-surface. `docs/api/complete-supported-reference.md` marca correctamente 3 de ellas como Rust-only, contradiciendo `docs/api/python.md` y `docs/api/typescript.md`.

---

## 5. AUDITORÍA DE CI

| Job | Qué corre | Qué NO corre |
|---|---|---|
| `test` (ubuntu, full) | `cargo test --lib --release`; `cargo test --lib --release --features "${RUST_EXPERIMENTAL_FEATURES}"`; `cargo test --release --test ci` | **`cross_validation.rs` (38), todos los `alpha_*_validation.rs`, `beta_*_validation.rs`, `regression.rs`, `experimental/*`, `debug/*`, `analysis/*`** |
| `test` (macos/windows, smoke) | solo `cargo test --release --test ci` | todo lo demás |
| `lint` | `cargo fmt --all --check`; `cargo clippy --all-targets --features ...` | **clippy solo del paquete raíz** → `crates/cli` y `crates/python` nunca se lintean |
| `wasm` | `wasm-pack build` | no ejecuta `npm test` (no existe script); el harness WebGPU no existe |
| `docs` | build de VitePress | — |
| `release` | wheel + npm | — |

`RUST_EXPERIMENTAL_FEATURES` = `parallel alpha-dft alpha-reaxff alpha-mlff alpha-obara-saika alpha-cga alpha-gsm alpha-sdr beta-kpm beta-mbh beta-randnla beta-riemannian beta-cpm`.
**Omite:** `alpha-dynamics-live` (15 unit tests nunca corren), `alpha-imd`, `alpha-edl`, `alpha-periodic-linear`, `alpha-kinetics`, `alpha-render-bridge`, `alpha-reaction-dynamics`, `experimental-gpu`.

---

## 6. AUDITORÍA DE DOCS (63 archivos)

| Doc | Estado | Evidencia |
|---|---|---|
| `docs/index.md` | Stale | "(v0.5.0)" vs 0.15.2 |
| `docs/algorithm/*.md` (24) | Accurate salvo 2 | `webgpu-validation.md` documenta harness inexistente |
| `docs/api/cli.md` | Incompleto | 22 de 68 comandos |
| `docs/api/python.md` | Stale | 3 fns inexistentes + shapes erróneos + ~20 fns reales omitidas |
| `docs/api/rust.md` | Accurate | — |
| `docs/api/typescript.md` | Stale | 2 exports inexistentes + ~8 reales omitidos |
| `docs/api/complete-supported-reference.md` | Mayormente accurate, **orphan** | contradice python.md/typescript.md |
| `docs/benchmarks.md` | Plausible, no reproducido | sin job CI que regenere |
| `docs/experimental/*.md` (13) | Accurate | `alpha.md`, `beta.md`, `edl.md` son orphans |
| `docs/guide/*.md` (10) | 2 stale | `python.md` con nombres de función incorrectos en **todos** los bloques; `rust.md`/`typescript.md` con versión 0.9.1 |
| `docs/guide/reaction-viewer.md` | Misleading | el ejemplo usa potenciales analíticos y **nunca importa sci-form-wasm** |
| `docs/validation/cross-library-validation.es.md` (2141 líneas) | **Materialmente inconsistente** | §13.1 declara PASS con valores para F23/F64/F65/F66/F69 que **no tienen test**; `reports/validation/summary.json` dice `value:null, pass:null` para F65/F66/F69 |
| 3 roadmaps + `dynamics_reactivity_roadmap.md` + `nmr_detailed_roadmap.md` + `spectroscopy_roadmap.md` | Orphan | fuera del nav |

---

## 7. PLAN DE ACCIÓN PRIORIZADO

### Sprint 0 — integridad (1–2 días)
1. Quitar los symlinks `pkg`/`pkg-node` de git y añadir al `.gitignore`.
2. P0-3: eliminar o marcar explícitamente como `unimplemented` las ~26 funciones que devuelven constantes.
3. P0-5: arreglar `cmd_ani` (usar la longitud AEV real 4032) o convertirlo en error explícito.
4. P0-8: borrar las 5 filas PASS ficticias de `cross-library-validation.es.md` y de `summary.json`, o implementar los tests.
5. P0-10: extender clippy a `--workspace`.
6. P0-1: envolver `request_adapter` en `catch_unwind` y pre-detectar entorno headless → fallback a CPU determinista.

### Sprint 1 — honestidad numérica (1 semana)
7. P0-2: implementar tablas MMFF94 reales o **marcar el módulo como no funcional** y renombrar la API.
8. P0-6: corregir el índice angular de `ani/aev.rs` (copiar el packing de `ani_tm.rs`) + test de bijectividad; corregir `ani/gradients.rs` para propagar el bloque angular; añadir FD-check de fuerzas.
9. P0-7: unificar unidades en `langevin.rs` (kcal/mol/Å ↔ eV/Å/amu) y añadir test de distribución de temperatura canónica.
10. P0-4: reemplazar `atomicAdd` sobre floats bitcast por reducción de 2 pasadas o fixed-point; añadir test GPU↔CPU real.
11. P2-36: rellenar los 256 slots radiales en `ml/symmetry_functions.rs`; regenerar el oráculo F79.
12. P2-41: implementar `coulomb_energy` real en DFT.
13. P2-28: cross-checks entre implementaciones duplicadas (SSSR, energía FF, RMSD).

### Sprint 2 — cerrar los gates de validación (2–3 semanas)
14. P0-1 de la auditoría de gates: **implementar F69 (XRD vs pymatgen), F64 (space groups vs spglib), F65 (CIF vs pymatgen), F66 (PBC vs ASE)** — las fixtures ya existen, solo falta el test.
15. P1-1/P1-2: gate de paridad ETKDG + bounds vs RDKit (des-ignorar `test_bounds_compare.rs` y commitear las fixtures RDKit).
16. P1-3/P1-4: gate de paridad del FF ETKDG-3D y del BFGS.
17. P1-7: gate de gradientes analíticos vs numéricos para SCF, PM3, xTB, DFT, EHT.
18. P1-5/P1-6: gate de SMARTS (match counts vs RDKit) + hash de la tabla CSD.
19. P1-10 a P1-25: gates de ESP, NPA, población EHT, CIP, GB/ALPB, variable-cell, descriptores ML, NMR, UV-Vis, IR intensities, EEQ2019, ANI, EHT metales, CI-NEB/IRC.

### Sprint 3 — infraestructura de tests (1 semana)
20. P2-1: añadir a CI `cargo test --test cross_validation --test regression` y un job matrix para `alpha_*`/`beta_*`.
21. P2-2: añadir `alpha-dynamics-live` a `RUST_EXPERIMENTAL_FEATURES`.
22. P2-3: migrar los gates de `tests/experimental/test_*` a las features `alpha-*`/`beta-*` actuales.
23. P1-18: añadir tests para `compute_reaction_dynamics` y sus helpers.
24. Convertir `tests/integration/test_python_integration.py` a pytest y añadirlo a CI; arreglar la ruta del test JS; crear `npm test`.
25. P2-24: corregir el `package.json` (script `test`, `main`).

### Sprint 4 — superficie y docs (1 semana)
26. P0-9: limpiar los symlinks y regenerar los paquetes WASM con un único script.
27. P2-22: unificar `post-wasm-build.mjs` y `crates/wasm/build.sh`.
28. P2-13: regenerar los `.d.ts` desde wasm-bindgen en CI (drift check).
29. P2-4 a P2-10: corregir docs (`python.md`, `typescript.md`, `guide/python.md`, `api/cli.md`, `index.md`, versiones), y enlazar los orphans en el nav.
30. P2-14: sincronizar `alpha_modules_info`/`beta_modules_info` con la versión del crate.
31. P2-25: reescribir `experimental_status.rs` con los paths y features actuales.

---

*Fin del informe. Generado por auditoría automatizada + verificación manual de 2026-10-20.*
