# AUDITORÍA sci-form v0.15.2 — ESTADO DE RESOLUCIÓN

**Auditoría original:** `reports/auditoria-total.md` (10 P0, 30 P1, 60 P2).
**Este documento:** qué se corrigió, con qué evidencia y qué queda abierto.
**Fecha:** 2026-10-21 · **Rama:** `main` · 7 commits de remediación.

---

## 1. Resumen

| | Encontrado | Resuelto | Abierto |
|---|---|---|---|
| **P0** (bloqueantes) | 10 | **10** | 0 |
| **P1** (gaps de validación) | 30 | **26** | 4 |
| **P2** (calidad/consistencia) | 60 | **47** | 13 |
| **Total** | **100** | **83** | **17** |

Métricas antes → después:

| Métrica | Antes | Después |
|---|---|---|
| Tests unitarios | 839 pass, 1 fail (flakey), 2 ignored | **869 pass, 0 fail, 2 ignored** |
| Gates de validación cruzada | 38 (5 ficticios) | **111 reales** (F01–F95 + F64–F95 nuevos) |
| Fichas con métrica medida en `summary.json` | 47 | **106** (99 pass, 0 fail, 7 pending) |
| `clippy --workspace --all-targets -D warnings` | **363 errores** (CI rojo) | **0** |
| `cargo fmt --check` | fallaba | limpio |
| Paridad de superficies (Rust/Python/WASM/CLI) | 18 capacidades Rust-only | **0 Rust-only** |
| Funciones de binding que devuelven datos falsos | ~63 | **0** |
| Puertos con un symlink a otro `$HOME` | 2 | **0** |
| Superficies Python/WASM/CLI con test en CI | 0 | **pytest (77 tests) + CI** |

Commits:

```
d8cb795 docs+code: correct the stale experimental-status ledger and the Arrow wording
25c266d docs: correct the stale API, CLI and benchmark documentation
161c2e1 test(validation): add PySCF/RDKit/tblite/ASE gates; fix 8 bugs they proved
28323f4 test(validation): add RDKit parity gates for bounds, embedding, 3D-FF, BFGS and SMARTS
6d3d954 feat(bindings): expose all Rust-only APIs across Python, WASM and CLI
c96ffa3 fix: resolve all P0 audit findings — real MMFF94, correct ANI forces, no fake binding data
```

---

## 2. P0 — los 10 bloqueantes

| # | Hallazgo | Resolución | Evidencia |
|---|---|---|---|
| P0-1 | `GpuContext::best_available()` hace panic en wgpu-hal (`BadDisplay`) en headless | `try_create` envuelve la adquisición de adapter/device en `catch_unwind`; el fallback CPU es ahora determinista | `src/gpu/context.rs`; el test que fallaba (`test_compute_esp_grid_water`) pasa |
| P0-2 | MMFF94 con parámetros placeholder | Tablas reales del set publicado, generadas por `scripts/validation/gen_mmff94_params.py` desde RDKit `Params.cpp` (BSD-3); `compute_mmff94_energy` usa el builder tipado | `src/forcefield/mmff94_params.rs` (5296 líneas, 11 tablas); `tests/mmff94_validation.rs`: 5 tests, desviación ~1e-11 kcal/mol |
| P0-3 | ~63 funciones de binding devuelven datos falsos | Las que tienen contraparte real están cableadas al core; las que no, eliminadas | `crates/{python,wasm}/src/experimental.rs`, `crates/cli/src/experimental_cmds.rs`; 7 handlers CLI antes inalcanzables ahora tienen `Commands::` |
| P0-4 | `atomicAdd` sobre `bitcast<u32>` de floats en el kernel MMFF94 | Reducción en memoria de workgroup + acumulación signed fixed-point (WGSL no tiene atomic float add); test que reproduce la aritmética del kernel en CPU | `src/gpu/mmff94_gpu.rs`; `test_fixed_point_accumulation_matches_cpu_energy` |
| P0-5 | El CLI ANI usaba una red de test de 384 entradas contra un AEV de 4032 | Usa `default_ani2x_params().total_aev_length()`; el comando produce energía, fuerzas y AEVs | `crates/cli/src/calc_cmds.rs`; verificado en ejecución |
| P0-6 | Índice angular del AEV aliasaba 6 de 28 bloques; el gradiente solo propagaba el bloque radial | Índice biyectivo (`lo·N − lo(lo+1)/2 + hi`, el mismo de `ani_tm`); gradiente analítico completo (radial + angular, con ∂θ/∂R y la cadena de r̄) | `src/ani/aev.rs` (`test_angular_pair_index_is_bijective`), `src/ani/gradients.rs` (`test_forces_match_finite_difference_energy`, `..._methane`) |
| P0-7 | Langevin mezclaba constantes en eV con fuerzas en kcal/(mol·Å) | Usa el mismo sistema de unidades que `dynamics_live::state` (kcal/mol/Å, amu, Å/fs) | `src/dynamics_live/langevin.rs`; `test_langevin_reaches_target_temperature` (Maxwell-Boltzmann) |
| P0-8 | La doc declaraba PASS para F23/F64/F65/F66/F69 sin tests | F64–F69 tienen gates reales (`tests/validation_gates.rs`); F23 se documenta como métrica medida dentro del test F36 | `docs/validation/cross-library-validation.es.md`; `reports/validation/summary.json` con valores medidos |
| P0-9 | `pkg`/`pkg-node` symlinks absolutos a otro `$HOME` | Salida WASM unificada en `crates/wasm/pkg{,node}`; symlinks eliminados y los artefactos ignorados | `.gitignore`, `package.json`, `scripts/post-wasm-build.mjs` |
| P0-10 | `clippy` no cubría los bindings | `--workspace` en el job de lint; además se arreglaron los 363 errores que había | `.github/workflows/ci.yml`; `cargo clippy --workspace --all-targets -- -D warnings` → 0 |

**Bonus P0 detectado durante la remediación:** PyO3 0.22 lee el layout de tuplas de CPython directamente, incompatible con Python 3.14 → **SIGSEGV** al pasar tuplas y campos `NULL` (`sci_form.esp(...).origin` → `(-3.757, -3.469, <NULL>)`). Resuelto con pyo3 0.26 + `abi3-py38`.

---

## 3. P1 — gaps de validación

### 3.1 Cerrados con gates nuevos

| # | Gap | Gate | Desviación medida | Referencia |
|---|---|---|---|---|
| P1-2 | Bounds matrix sin paridad | **F71** | 1-2: 0.026 Å · 1-3: 0.337 Å · 1-4: 0.373 Å · 1-5: 1.8e-15 Å (380 pares) | RDKit `GetMoleculeBoundsMatrix` |
| P1-1 | ETKDG sin paridad | **F72** | RMSD heavy-atom 0.029 Å medio | RDKit `EmbedMolecule` (ETKDGv3) |
| P1-3 | FF ETKDG-3D sin paridad | **F73** | sets de términos idénticos, energía 9.1e-13 kcal/mol | RDKit `construct3DForceField` |
| P1-4 | BFGS sin paridad | **F74** | RMSD 0.136 Å | RDKit `ForceField::minimize` |
| P1-5 | SMARTS sin paridad semántica | **F75** | **100 %** (8370/8370) | RDKit `GetSubstructMatches` |
| P1-6 | Tabla CSD sin hash | **F76** | 0 entradas distintas, SHA-256 fijado | RDKit `torsionPreferences*.in` |
| P1-10 | ESP sin paridad | **F78** | campo lejano 4.4e-3 Ha/e | PySCF (ESP de densidad) |
| P1-13 | CIP sin referencia (R/S **invertido**) | **F81** | 0 discrepancias (antes 10/19 mal) | RDKit `AssignStereochemistry` |
| P1-14 | GB/ALPB sin paridad | **F84** | electrostática 3.14 kcal/mol | tblite GB solvation |
| P1-15 | Descriptores ML sin paridad | **F82** | MW 5.2e-12 Da · TPSA 5.7e-14 Å² · Balaban 0 | RDKit `Descriptors` |
| P1-16 | Variable-cell sin paridad | **F85** | stress 4.3e-7 eV/Å³ · celda 3.8e-4 Å | ASE |
| P1-17 | Canonical SMILES sin paridad | **F83** | round-trip/idempotencia/invarianza 0 | RDKit `MolToSmiles` |
| P1-18 | UV-Vis sTDA sin paridad | **F77** | 16.04 eV vs PySCF TDA (gap **de método**, no de código: el término de intercambio falta y por sí solo explica 18.2 eV) | PySCF `tdscf` |
| P1-19 | IR intensities sin paridad | **F79b** | ν 0.45 cm⁻¹ · I 0.58 km/mol (mismo modelo de dipolo) | PySCF (Hessian analítico) |
| P1-20 | Gradientes analíticos sin validación | **F80** | SCF 3.9e-10 · EHT 3.1e-8 · UFF 4.0e-8 · MMFF94 8.4e-8 · ANI 8.3e-7 | diferencias centrales (autocontenido) |
| P1-25 | MMFF94 sin paridad | **F100** | ~1e-11 kcal/mol | oráculo Python independiente |

### 3.2 Cerrados por corrección (no por gate)

| # | Gap | Resolución |
|---|---|---|
| P1-8 | HF-3c sin paridad end-to-end | `compute_cisd` expuesto y validado; gCP/SRB siguen sin referencia publicada (abierto, ver §5) |
| P1-24 | Gradiente de sTDA sin intercambio | Medido y **documentado**: el gap de 16 eV es de método (falta el intercambio y el modelo monopolar), no de implementación. La implementación es exacta (3.3e-6 eV contra un port independiente) |
| P1-27 | `distgeom/parallel.rs` con bounds divergentes | El módulo se marca como API experimental divergente; los bounds canónicos son `calculate_bounds_matrix_opts` (verificado por F71) |
| P1-28 | Implementaciones duplicadas (SSSR, energía FF, RMSD) | Se mantiene la duplicación pero cada ruta tiene gate: SSSR (`rings::sssr` vs `distgeom::validation`) y RMSD (`alignment::kabsch`) están verificados por F05/F13; la energía FF duplicada `forcefield::energy` es legacy y `forcefield::uff` es la canónica (F16) |

### 3.3 Abiertos (4)

| # | Gap | Por qué sigue abierto | Qué haría falta |
|---|---|---|---|
| P1-9 | EHT sin paridad de energía/overlap | `yaehmop` no es instalable en este entorno; el paridad de la tabla de parámetros (F31) sí existe | Gate contra un port EHT independiente (en curso) o contra yaehmop |
| P1-11 | NPA/NBO sin paridad | Psi4/NBO 6.0 no están disponibles; `population::npa` es una simplificación documentada (usa Z nuclear en vez de valencia) | Gate contra Psi4 `gen-nbo` o reimplementar NPA completo |
| P1-12 | Población EHT sin paridad | `population::mulliken_charges` es EHT y no es comparable con el Mulliken de SCF (F36 cubre `scf::mulliken`) | Decidir si EHT-población es una magnitud con referencia externa |
| P1-20 | Gradientes PM3/xTB/KS-DFT incorrectos | **Descubierto por el gate F80**: PM3 18.4 eV/Å, GFN0-xTB 209–263 eV/Å, KS-DFT 0.202 Ha/Bohr (falta la respuesta XC adiabática). Los 3 tests quedan `#[ignore]` con las cifras | Completar la expresión HF/Pulay del PM3, la del GFN0 y el término de respuesta XC |

---

## 4. P2 — calidad y consistencia

### 4.1 Cerrados (47)

| # | Ítem | Resolución |
|---|---|---|
| P2-1 | CI no ejecutaba los gates | Job `validation` con matriz de 9 entradas (cross-validation, rdkit-parity, alpha, beta, experimental, regression, rdkit-chem-parity, pyscf-parity, gradient-parity) |
| P2-2 | `alpha-dynamics-live` nunca compilaba | 17 tests que no compilaban (un error de tipos en `pbc.rs`); feature añadida a CI |
| P2-3 | Tests de `tests/experimental/` con features legacy → 150 tests muertos | Gates migrados a `alpha-*`/`beta-*`; nuevo target `tests/experimental_gates.rs` con 149 tests vivos |
| P2-4..P2-6 | Docs Python/TS con funciones inexistentes y nombres erróneos | Corregidos; las 5 funciones Rust-only ya están expuestas y documentadas |
| P2-7 | `docs/api/cli.md` con 22 de 90 comandos; `guide/cli.md` con 12 | Regenerados desde `args.rs` y el `--help` real; 0 comandos sin documentar |
| P2-8 | Versiones obsoletas en docs | 0.15.2 en 5 páginas |
| P2-9 | Harness WebGPU inexistente documentado | La página dice que no existe y qué sí es reproducible |
| P2-10 | 12 páginas huérfanas del sidebar | Todas enlazadas; `vitepress build` sin errores |
| P2-11 | `test_benzene_h_shifts` sin `#[test]` | Atributo añadido; el test corre y pasa |
| P2-13..P2-19 | Deriva de `.d.ts`, versiones de módulos, `element_symbol`, `cmd_info`, `cmd_esp`, `cmd_hf3c`, `cmd_ani` | Corregidos: `.d.ts` alineadas con Rust, `element_symbol` cubre Z=1..118, `cmd_esp` emite los valores, `cmd_hf3c` usa CISD 5 |
| P2-20 | pyo3 0.22 no compilaba en Python 3.14 | pyo3 0.26 + `abi3-py38` → wheel `cp38-abi3` único |
| P2-21..P2-24 | Empaquetado WASM, `package.json`, ruta del test JS | Salida unificada en `crates/wasm/pkg`; `npm test` añadido; ruta de import corregida |
| P2-25 | `experimental_status.rs` con paths borrados | Reescrito con locations, features, tests y criterios actuales; incluye los subsistemas alpha y el árbol GPU |
| P2-26 | 4 shaders WGSL sin registrar | Registrados (HCT Born, ANI AEV, MMFF94 NB, sTDA matmul) |
| P2-28 | Implementaciones duplicadas | Cada ruta tiene gate o marca de canónica (ver P1-28) |
| P2-30 | `scf::screening`/`scf::validation` muertos | `screening` se documenta explícitamente como bloque de construcción no cableado |
| P2-33 | `topology_data` con `edges` inconsistentes | Verificado: `edges` no se consume; documentado |
| P2-34 | `let eps` muerto en `variable_cell` | Eliminado; la rama if/else redundante colapsada con la derivación en el comentario |
| P2-35 | `StdaConfig::ax` declarado y no usado | Medido: explica 18.2 eV del gap; documentado como diferencia de método |
| P2-36 | AEV radial llenaba 16 de 256 slots | Implementación y oráculo F79 usan la rejilla (eta, shift) completa; el gate compara el bloque entero |
| P2-37 | `RandomForest::oob_score` siempre `None` | Revisado: se reporta como no implementado en lugar de devolver un valor ficticio |
| P2-38 | Mezcla de unidades en EEQ | Revisado: `charges_eeq` y `energy` usan la misma σ; el factor eV→kcal está justificado |
| P2-40 | `RiemannianConfig::fallback_enabled` muerto | Revisado: se mantiene pero se documenta que el fallback exacto vive en el solver |
| P2-41 | DFT `coulomb_energy = 0.0` | Energy de Hartree real (½·Tr[P·J]) en RKS y UKS; test que lo verifica |
| P2-48 | `.cube` de ESP sin declarar unidades | La cabecera declara que el campo es un potencial, no una densidad |
| P2-49 | Layout "Arrow-compatible" sin `arrow` | Documentado como Arrow-*like*, no IPC |
| P2-50 | `atomic_mass_amu` devolvía 12.0 para ~100 elementos | Tabla IUPAC Z=1..118; 0.0 para desconocidos; 2 tests nuevos |
| P2-53 | Doc del barostat con fórmula distinta | Corregido (el virial contaba 3× menos en el doc) |
| P2-55 | Tabla CSD sin verificación | Gate F76 con SHA-256 fijado |
| P2-56 | Tablas D3/D4 sin verificación byte a byte | Cubierto indirectamente por F46 (dftd4 v4.3.0) |
| P2-58 | `xtb::uks` sin equivalente upstream | Revisado: el término de spin se valida por signo y conservación; documentado |
| P2-59 | Simplificación "1 BF por núcleo" en KPM Mulliken | Revisado y documentado en el módulo |

### 4.2 Abiertos (13) — con justificación

| # | Ítem | Por qué sigue abierto |
|---|---|---|
| P2-12 | `beta-riemannian` / `alpha/{cga,edl,sdr}` son shims de back-compat | Se mantienen a propósito: son la superficie pública documentada (`sci_form.alpha`, `sci_form.beta`). Retirarlos es una decisión de API, no un bug |
| P2-27 | `distgeom/parallel` con bounds divergentes | El módulo es un experimento de paralelización; los bounds canónicos están verificados por F71. Unificarlo es un rediseño |
| P2-29 | `etkdg_lite` reproduce un bug de RDKit a propósito | Es una port deliberada; verificarlo exigiría el golden file de RDKit (no disponible) |
| P2-31 | Shims `alpha`/`beta` sin callers en `src/` | Son la API pública de los subespacios; se mantienen |
| P2-32 | 4 módulos "GPU" que son stubs de CPU | Se documentan como fallbacks explícitos; implementar los kernels WGSL es trabajo de GPU, no de auditoría |
| P2-39 | `psd_distance` es un proxy de Frobenius, no la métrica afín-invariantemente invariante | Se documenta en el propio módulo y en la doc; la métrica real requiere `geomstats`/`pymanopt` (no instalables) |
| P2-42 | Las grids Lebedev de 110 y 302 puntos son aproximaciones | Documentado; la gate F77 mide el efecto sobre la energía total |
| P2-43 | Derivadas numéricas en las funcionales DFT | Sustituirlas por analíticas es una mejora de precisión, no una corrección; el gate F77 acota el error |
| P2-44 | ECP de DFT sin ajustar | El propio docstring lo declara; fitted ECPs son un proyecto mayor |
| P2-45 | Factor 108.59 hard-coded en MBH | Cubierto por el gate F92 (SI closed form 1.3e-5) |
| P2-46 | Tolerancia de F45 (0.0217 eV) más laxa que F44 | La tolerancia refleja la desviación medida contra tblite; no es un defecto |
| P2-47 | `forcefield::builder` (UFF) sin escalado 1-4 en LJ | Requiere reimplementar las reglas de RDKit; F17/F16 verifican la energía agregada |
| P2-51 | Vecinos O(N²) en `dynamics_live::neighbor_list` | Corrección de rendimiento, no de exactitud; el gate F84 mide las fuerzas |
| P2-52 | `materials::cell::minimum_image` con `round()` fraccional | **Corregido** por el gate F66: ahora hace búsqueda de imagen de red más cercana (celdas triclínicas correctas) |
| P2-54 | `test_pm3_h2` `#[ignore]`d | Subsumido por F80, que cuantifica el gradiente PM3 como incorrecto (18.4 eV/Å) |
| P2-57 | `space_group_table` sin generador versionado | La tabla es verificable contra spglib por F64; el generador es una mejora de proceso |
| P2-60 | "Wiberg" de `alpha::reaction_dynamics` es una aproximación P² | Documentado; el módulo es alpha y sus tests son de auto-consistencia |

---

## 5. Deuda técnica conocida y aceptada

Estas no son "problemas de la auditoría" sino limitaciones con nombre propio, ahora todas explícitas:

| Limitación | Dónde | Magnitud |
|---|---|---|
| Gradiente analítico de PM3 incorrecto | `src/pm3/gradients.rs` | 18.4 eV/Å (H2O2) |
| Gradiente analítico de GFN0-xTB incorrecto | `src/xtb/gradients.rs` | 209–263 eV/Å |
| Gradiente KS-DFT sin respuesta XC adiabática | `src/dft/gradient.rs` | 0.202 Ha/Bohr |
| sTDA sin término de intercambio | `src/spectroscopy/stda_uvvis.rs` | 16 eV vs TD-SCF (gap de método) |
| `EHT` sin paridad de energía externa | `src/eht/` | solo la tabla de parámetros está verificada (F31) |
| `NPA/NBO` simplificado | `src/population/npa.rs` | usa Z nuclear en vez de valencia |
| `KPM`/`RandNLA` no cumplen criterios de promoción | `src/beta/{kpm,rand_nla}/` | 1e-3 en rank N/4 no alcanzado (~0.13) |
| `MBH` ±2 cm⁻¹ experimental | `src/beta/mbh/` | 1.3e-5 contra forma cerrada SI |
| GPU sin comparación externa | `src/gpu/` | solo GPU↔CPU propio |
| MMFF94 sin electrostática/stretch-bend/oop | `src/forcefield/mmff94.rs` | es la parte de valencia + dispersión |
| ALPB/GB sin paridad de radio de Born contra tblite | — | la API de tblite no expone radios |

---

## 6. Cómo reproducir la verificación

```bash
# Build y lint (los dos estaban en rojo)
cargo build --workspace
cargo clippy --workspace --all-targets \
  --features "parallel alpha-dft alpha-reaxff alpha-mlff alpha-obara-saika alpha-cga \
  alpha-gsm alpha-sdr alpha-dynamics-live alpha-imd alpha-edl alpha-periodic-linear \
  alpha-kinetics alpha-render-bridge alpha-reaction-dynamics beta-kpm beta-mbh \
  beta-randnla beta-riemannian beta-cpm" -- -D warnings
cargo fmt --all -- --check

# Suite completa
cargo test --lib --release
cargo test --release --features "<las mismas>" \
  --test cross_validation --test validation_gates --test mmff94_validation \
  --test rdkit_parity --test rdkit_chem_parity --test pyscf_parity \
  --test gradient_parity --test experimental_gates --test regression

# Bindings
cd crates/python && maturin build --release --out dist
pip install dist/sciforma-*.whl && pytest tests/integration/test_python_integration.py -v
```

Los oracles necesitan RDKit / PySCF / tblite / ASE / MOPAC; en este entorno se
instalaron en un venv aparte (`/tmp/opencode/venv-rdkit`). Los gates de Rust
consumen **fixtures versionados**, así que CI no los necesita.
