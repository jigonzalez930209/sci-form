# Plan de validación cruzada por funcionalidad y por crate

> **Objetivo** — para cada funcionalidad de sci-form: (1) qué programa/librería probada hace lo mismo,
> (2) qué test equivalente permite comparar ambos, (3) qué ficheros o datasets ya calculados o
> experimentales existen para medir la desviación, y (4) el trabajo detallado crate por crate para
> dejar la lista operativa de cómo probar cada funcionalidad contra otras librerías o contra datos
> ya resueltos.
>
> **Fecha de redacción:** 2026-10-11 · **Versión auditada:** `sci-form 0.15.2` (`Cargo.toml`)
> **Naturaleza del documento:** plan de trabajo vivo. No cambia código; define los gates de veracidad,
> las referencias, los datasets y las tolerancias que deben cumplir los tests de comparación.

---

## Índice

1. [Cómo usar este documento](#1-cómo-usar-este-documento)
2. [Mapa de crates y superficies](#2-mapa-de-crates-y-superficies)
3. [Metodología: niveles, reglas de comparación y tolerancias](#3-metodología-niveles-reglas-de-comparación-y-tolerancias)
4. [Fichas por funcionalidad](#4-fichas-por-funcionalidad)
   - 4.1 [Parsing y topología](#41-parsing-y-topología)
   - 4.2 [Geometría 3D y conformeros](#42-geometría-3d-y-conformeros)
   - 4.3 [Fuerzas y optimización](#43-fuerzas-y-optimización)
   - 4.4 [Cargas, electrostática y solvatación](#44-cargas-electrostática-y-solvatación)
   - 4.5 [Estructura electrónica](#45-estructura-electrónica)
   - 4.6 [Espectroscopía](#46-espectroscopía)
   - 4.7 [ML, descriptores y potenciales neuronales](#47-ml-descriptores-y-potenciales-neuronales)
   - 4.8 [Materiales y sistemas periódicos](#48-materiales-y-sistemas-periódicos)
   - 4.9 [Dinámica, transporte y utilidades](#49-dinámica-transporte-y-utilidades)
   - 4.10 [Rutas GPU](#410-rutas-gpu)
   - 4.11 [Módulos alpha](#411-módulos-alpha)
   - 4.12 [Módulos beta](#412-módulos-beta)
5. [Crate `sci-form-cli`](#5-crate-sci-form-cli)
6. [Crate `sci-form` (bindings Python)](#6-crate-sci-form-bindings-python)
7. [Crate `sci-form-wasm`](#7-crate-sci-form-wasm)
8. [Programas de referencia y datasets descargables](#8-programas-de-referencia-y-datasets-descargables)
9. [Matriz de trazabilidad funcionalidad → test → dataset](#9-matriz-de-trazabilidad-funcionalidad--test--dataset)
10. [Plan de ejecución por fases](#10-plan-de-ejecución-por-fases)
11. [Registro de resultados y CI](#11-registro-de-resultados-y-ci)
12. [Riesgos, trampas y reglas anti-falso-positivo](#12-riesgos-trampas-y-reglas-anti-falso-positivo)
13. [Estado de implementación](#13-estado-de-implementación)

---

## 1. Cómo usar este documento

- Cada **ficha** es autosuficiente y sigue siempre el mismo patrón:
  `API → Referencia probada → Test equivalente → Datasets → Procedimiento → Estado en el repo → Criterio de aceptación`.
- Las rutas de código son las reales del repo: `src/<módulo>/…`, `tests/<suite>…`, `scripts/…`.
- Los **niveles de validación** (L0…L4) indican qué clase de prueba corresponde; se explican en §3.
- La **matriz de trazabilidad** (§9) es la lista de trabajo plana: una fila por funcionalidad,
  con el test existente, el test nuevo propuesto, la referencia y el dataset.
- Los **datasets** (§8.2) se descargan a `tests/fixtures/` o `data/` (ambos ya ignorados por git o
  excluidos del paquete), nunca se commitean binarios pesados sin comprimir.
- Regla de oro: **primero la referencia misma-método, después la referencia independiente, al final
  el experimento**. Comparar un método semi-empírico contra experimento como único gate produce
  falsos fallos; compararlo contra su propia implementación canónica detecta bugs reales.

---

## 2. Mapa de crates y superficies

El workspace tiene **4 crates** (raíz `[workspace] members` en `Cargo.toml`):

| Crate | Ruta | Superficie | Cómo se prueba |
|---|---|---|---|
| `sci-form` (core) | `src/` (331 archivos `.rs`, ~1048 `pub fn`) | Rust API + módulos alpha/beta tras feature gates | `cargo test --lib`, `cargo test --tests --release` |
| `sci-form-cli` | `crates/cli/` | ~85 subcomandos `cmd_*` (embed, parse, eht, pm3, xtb, gfn1, gfn2, hf3c, dos, esp, charges, sasa, solvation, sssr, ecfp, tanimoto, stereo, rmsd, uff, population, dipole, assemble, cell, kmesh, charts, report, validator, NEB/GSM/kinetics, edl, eeq, d4, alpb, cpm…) | smoke + golden outputs |
| `sci-form` (Python) | `crates/python/` | submódulos `sci_form`, `sci_form.alpha`, `sci_form.beta` (PyO3) | `maturin develop` + `pytest` |
| `sci-form-wasm` | `crates/wasm/` (+ `crates/wasm/js/*`) | 180 exports wasm-bindgen, subpaths `/alpha`, `/beta`, typed arrays | `wasm-pack build` + Node.js |

### 2.1 Mapa interno del crate core (por dominio)

| Dominio | Módulos (`src/`) | Fichas |
|---|---|---|
| Parsing y topología | `smiles`, `canonical_smiles`, `graph`, `smarts`, `smirks`, `rings`, `topology`, `stereo`, `materials/cif` | 4.1 |
| Geometría 3D | `distgeom`, `etkdg`, `conformer`, `alignment`, `clustering` | 4.2 |
| Fuerzas y optimización | `forcefield` (uff, mmff94, bounds_ff, etkdg_3d, reaxff), `optimization` | 4.3 |
| Cargas / electrostática / solvatación | `charges`, `charges_eeq`, `population`, `dipole`, `esp`, `surface`, `solvation`, `solvation_alpb`, `dispersion` | 4.4 |
| Estructura electrónica | `eht`, `dos`, `scf`, `hf`, `pm3`, `xtb`, `dft` (alpha), `reactivity` | 4.5 |
| Espectroscopía | `ir`, `spectroscopy`, `nmr`, `reactivity` (sTDA/UV-Vis) | 4.6 |
| ML y descriptores | `ml`, `ani`, `reactivity`, `potentials` | 4.7 |
| Materiales | `materials`, `periodic` | 4.8 |
| Dinámica / transporte | `dynamics`, `dynamics_live`, `transport`, `mesh`, `gpu` | 4.9, 4.10 |

### 2.2 Feature gates que cambian lo que se prueba

`parallel` (default), `experimental-gpu`, `alpha-{dft,reaxff,mlff,obara-saika,cga,gsm,sdr,dynamics-live,imd,edl,periodic-linear,kinetics,render-bridge,reaction-dynamics}`,
`beta-{kpm,mbh,cpm,randnla,riemannian}`.

Cada gate define su propio comando de compilación y su propio gate de validación; no mezclar
resultados con y sin feature en la misma tabla de referencia.

---

## 3. Metodología: niveles, reglas de comparación y tolerancias

### 3.1 Niveles de validación

| Nivel | Qué compara | Ejemplo | Coste |
|---|---|---|---|
| **L0 — Invariante** | propiedades que deben cumplirse siempre, sin referencia externa | suma de cargas = carga total; DOS integra al nº de orbitales; gradiente analítico = diferencias finitas; energía conservada en MD | bajo |
| **L1 — Misma implementación** | mismo resultado en otra superficie del repo | Rust vs Python vs WASM vs CLI con la misma semilla; `parallel` vs serial | bajo |
| **L2 — Misma-método (referencia canónica)** | misma fórmula/parámetros en una implementación probada | UFF vs RDKit UFF; GFN2 vs tblite; PM3 vs MOPAC; D4 vs `dftd4` | medio |
| **L3 — Método independiente** | mismo observable por una teoría distinta | EHT vs RHF/STO-3G para tendencias; sTDA vs TD-DFT; EHT-geometría vs ETKDG | medio-alto |
| **L4 — Experimento / dataset ya resuelto** | desviación contra realidad experimental o dataset congelado | dipolos vs NIST; frecuencias vs CCCBDB; GDB-20 vs RDKit; NMR vs SDBS | alto |

Un gate de release debe apoyarse en **L0+L1 obligatorios** y en **L2 con tolerancia estrecha**;
L3/L4 se usan para documentar precisión, no para bloquear release salvo desviación grosera.

### 3.2 Reglas de comparación (evitan comparaciones inválidas)

1. **Misma geometría primero.** Toda comparación de energías/cargas/gradientes se hace sobre las
   *mismas coordenadas congeladas*; si cada programa optimiza, se compara después de re-optimizar
   con el mismo método.
2. **Unidades explícitas.** EHT en eV, PM3 HOF en kcal/mol, xTB/HF/DFT en eV o Hartree según API,
   UFF/MMFF en kcal/mol, frecuencias cm⁻¹, desplazamientos ppm, dipolos Debye. Convertir *antes* de
   comparar; hay un incidente histórico documentado en `docs/algorithm/quantum-validation-roadmap.md`
   (HF-3c vs "experimental SCF": 3031 mHa de diferencia por definiciones incompatibles).
3. **Autoestados:** comparar `|C|` o proyectar subespacios degenerados; nunca firmas elemento a
   elemento sin alinear fase y orden.
4. **Permutación de átomos:** alinear por isomorfismo de grafo (mapa de átomos), nunca por índice de
   archivo. Con RDKit: `GetSubstructMatch`; en Rust: `smarts::has_substruct_match`.
5. **Orientación/centro:** RMSD/Kabsch antes de comparar geometrías; comparar distancias internas
   (pairwise distance matrix) cuando el alineamiento global sea ambiguo.
6. **Aleatoriedad:** fijar semilla en todos los implicados (RDKit `randomSeed`, `--seed` sci-form);
   repetir 20–50 semillas y comparar distribuciones, no un único draw.
7. **No unicidad:** SSSR no es único (comparar multiset de tamaños, no el set exacto); SMILES
   canónico no coincide entre implementaciones (comparar isomorfismo + idempotencia, no strings).
8. **Topología propia:** cada programa construye su grafo desde SMILES con reglas propias
   (aromaticidad, protonación). Antes de comparar cualquier propiedad, verificar que elementos,
   enlaces, órdenes y cargas formales coinciden (test de "sanitización").

### 3.3 Tolerancias recomendadas por observable

| Observable | Nivel | Tolerancia objetivo | Nota |
|---|---|---|---|
| Coordenadas L1 (Rust/Py/WASM/CLI, misma semilla) | L1 | idénticas bit a bit si el camino es el mismo; ≤1e-12 Å si hay reordenación/álgebra distinta | usar coordenadas internas |
| Conformeros vs ETKDG RDKit (heavy-atom RMSD) | L4 | mediana ≤ 0.25 Å; P90 ≤ 0.5 Å; éxito ≥ 99 % en GDB-20 | repo ya tiene GDB-20 1k |
| Longitudes de enlace vs experimento/CSD | L4 | ≤ 0.02–0.03 Å | tras optimización, no en el embedding crudo |
| Ángulos | L4 | ≤ 2° | ídem |
| Matriz S (misma base STO-3G) | L2 | 1e-10 abs vs PySCF | integrales gaussianas |
| Matriz H EHT (misma tabla VSIP, K=1.75) | L2 | 1e-8 abs vs YAeHMOP/reimplementación | si las tablas difieren, comparar `H_ii` y fórmula |
| ERI (Obara-Saika) | L2 | 1e-10 rel vs PySCF/libcint | probar s, p, d |
| Eigenvalores HF (misma base) | L2 | 1e-6 Ha vs PySCF | |
| Energía total HF | L2 | 1e-7 Ha vs PySCF | |
| Energía total DFT (misma malla) | L2 | 1e-5 Ha vs PySCF | la malla domina; usar la misma |
| EHT energía orbital / gap | L2 | 1e-6 eV | |
| PM3 ΔHf (misma geometría) | L2 | ≤ 0.5 kcal/mol vs MOPAC | |
| GFN2/GFN1 energía (misma geometría) | L2 | ≤ 0.5 kcal/mol vs tblite; gradiente ≤ 1e-3 eV/Å | si el objetivo es fidelidad: 1e-6 Ha |
| GFN0 energía | L2 | ≤ 1 kcal/mol vs xtb GFN0 | |
| D3/D4 dispersión | L2 | 1e-3 kcal/mol vs `dftd3`/`dftd4` | mismos parámetros |
| UFF/MMFF energía | L2 | ≤ 0.01 kcal/mol vs RDKit (misma topología) | términos idénticos |
| Gradiente analítico vs numérico | L0 | error relativo ≤ 1e-4 (h=1e-3 Å, diferencias centrales) | para cada FF y método QM |
| Cargas Gasteiger vs RDKit | L2 | ≤ 0.02 e | ojo: nº iteraciones/config |
| Mulliken EHT vs referencia EHT | L2 | ≤ 0.005 e | definición idéntica |
| NPA/NBO | L2 | ≤ 0.01 e vs JANPA (misma densidad) | |
| Dipolo | L4 | ≤ 0.1 D o 5 % vs NIST/CCCBDB | |
| SASA (misma sonda/radios) | L2 | ≤ 1 % vs FreeSASA | |
| ESP grid | L2 | RMSE ≤ 1e-3 a.u. (misma densidad) | cubo de referencia |
| DOS/PDOS | L2 | MSE ≤ 0.01 vs Multiwfn / diagonalización exacta | KPM: L2 ≤ 2–5 % |
| Frecuencias armónicas misma-método | L2 | ≤ 2 cm⁻¹ (o 5 cm⁻¹ escalado) | |
| Frecuencias vs experimento | L4 | ±20 cm⁻¹ armónico puro; ±10 cm⁻¹ con factor de escala | documentar factor |
| IR intensidades | L4 | tendencias + orden relativo | |
| UV-Vis sTDA vs TD-DFT (misma base) | L3 | ≤ 0.3 eV por transición dominante | |
| UV-Vis λmax vs experimento | L4 | ≤ 30–50 nm para familias conjugadas | método exploratorio |
| NMR GIAO vs ORCA (misma base) | L2 | ≤ 0.5 ppm | shieldings |
| NMR HOSE vs experimento | L4 | MAE ¹H ≤ 0.3 ppm; ¹³C ≤ 3 ppm | |
| J (Karplus) vs experimento | L4 | ≤ 1 Hz | confórmeros promediados |
| ANI energía/fuerzas vs torchani (mismos pesos) | L2 | ≤ 1e-4 eV/átomo; ≤ 1e-3 eV/Å | float32 del modelo |
| ANI vs DFT (ANI test set) | L4 | MAE publicada del modelo (±20 %) | |
| Descriptores RDKit/mordred | L2 | exacto donde la definición coincide (MW, TPSA); ±0.1 logP Crippen | |
| ECFP Tanimoto | L2 | bits idénticos a RDKit ECFP4; ±1e-9 similitud | mismos invariantes |
| Butina | L2 | clusters idénticos salvo empates; ±0 | mismo orden de entrada |
| KPM DOS vs diagonalización | L2 | L2 rel ≤ 5 % (orden 100–500) | |
| RandNLA eigenvalores | L2 | rel ≤ 1e-3 en los k más bajos; fallback correcto | |
| ReaxFF vs LAMMPS (mismos params) | L2 | ≤ 1e-3 kcal/mol y ≤ 1e-3 kcal/mol/Å | |
| Riemannian SPD | L2 | ≤ 1e-9 vs geomstats/pymanopt | métrica affine-invariant |
| SDR embedding (violación de distancias) | L0/L2 | violación ≤ 1e-3 Å; cvxpy mismo objetivo ±1e-3 | |
| Kinus/HTST | L0/L2 | Eyring analítico; CatMAP mismo TST | |

### 3.4 Plantilla de ficha

```
**F<nn> · <nombre> (`<módulo>`, `<API>`)**
- Referencia probada: <programa/librería + función/CLI>
- Test equivalente: <qué se calcula en ambos, cómo se alinean>
- Datasets: <URL/uso o "no aplica">
- Procedimiento: <script/comando concreto>
- Estado en el repo: <test/fixture existente>
- Criterio: <tolerancia de §3.3>
```

## 4. Fichas por funcionalidad

> Convención: `Lx` = nivel de validación (§3.1). El "Criterio" usa las tolerancias de §3.3.

### 4.1 Parsing y topología

**F01 · Parser SMILES (`smiles`, `graph::Molecule::from_smiles`, `parse`)**
- Referencia probada: **RDKit** `Chem.MolFromSmiles` (BSD-3) y **Open Babel** `obabel -:"SMILES" -osmi` (GPL-2, solo entorno de desarrollo).
- Test equivalente: parsear el mismo corpus con las tres implementaciones y comparar **átomo a átomo**: número de átomos, elementos, cargas formales, órdenes de enlace, flags aromáticos, tags estereoquímicos (`@`/`@@`, `/\`), isótopos, radicales. Comparar isomorfismo de grafo, no el string.
- Datasets: `data/chembl_10k_practical.smi` (ya en repo), `scripts/10k_smiles.smi`, ZINC/Enamine (muestra), y el set de "SMILES inválidos" de la suite de RDKit (`rdkit.Chem.rdmolfiles.CanonicalRankAtoms` tests + `test_data`).
- Procedimiento: `python scripts/compare_pipeline.py` (ya existe, adaptar) → genera `tests/fixtures/parse_compare.json`; test Rust nuevo `tests/regression/test_smiles_cross.rs` que lee el fixture.
- Estado en el repo: `src/smiles.rs` (round-trip, stereo, isotopes); `src/graph.rs`; script `scripts/compare_pipeline.py`.
- Criterio: **100 %** de moléculas del corpus deben coincidir en composición y conectividad; tolerancia 0 en conteos; documentar cada discrepancia de orden aromático.

**F02 · SMILES canónico (`canonical_smiles::to_canonical_smiles`)**
- Referencia probada: RDKit `Chem.MolToSmiles` (canónico propio). **No** se comparan strings entre implementaciones: distintas canonicalizaciones.
- Test equivalente: (a) idempotencia (canónico de canónico = mismo string dentro de sci-form); (b) biyección: `parse(x)` ≅ `parse(canonical(x))` (isomorfismo RDKit); (c) unicidad interna: dos SMILES equivalentes del mismo corpus deben dar el mismo canónico; (d) round-trip `canonical(x) → parse → to_smiles` estable.
- Datasets: corpus ChEMBL 10k + `tests/regression/test_algorithm_validation.rs` si ya cubre round-trip.
- Procedimiento: script `scripts/compare_canonical_smiles.py` con RDKit como *validador de equivalencia*, no de string.
- Estado en el repo: `src/canonical_smiles.rs`; tests en `test_algorithm_validation.rs`.
- Criterio: 100 % idempotente y 100 % isomorfo; cero colisiones (dos grafos distintos con el mismo canónico).

**F03 · SMARTS (`smarts::parse_smarts`, `has_substruct_match`, `substruct_match(_batch)`)**
- Referencia probada: RDKit `MolFromSmarts` + `HasSubstructMatch` / `GetSubstructMatches` (BSD-3); Open Babel como segundo voto.
- Test equivalente: para cada patrón de `tests/fixtures/smarts_patterns.txt` y cada molécula del corpus, comparar **conteo de matches** y **átomos mapeados** (ordenados) contra RDKit, con `useChirality=false/true` y `uniquify`.
- Datasets: `tests/fixtures/smarts_patterns.txt`, ChEMBL 10k; suite de SMARTS de Daylight (patrones documentados en el manual) como set de casos límite.
- Procedimiento: `scripts/compare_smarts.py` → fixture JSON; test `cargo test --test regression test_smarts_cross`.
- Estado en el repo: `src/smarts/matcher.rs` (7 tests internos), `src/smarts/torsion_matcher.rs`.
- Criterio: igualdad exacta de conjuntos de matches; cualquier discrepancia se documenta como limitación explícita del subconjunto SMARTS soportado.

**F04 · SMIRKS y transformaciones (`smirks::parse_smirks`, `apply_smirks`, `apply_smirks_multi`)**
- Referencia probada: RDKit `AllChem.ReactionFromSmarts` + `RunReactants` (semántica de referencia); ya hay harness en el repo.
- Test equivalente: mismo reactivo + plantilla SMIRKS → comparar set de productos por isomorfismo (con y sin estereoquímica), incluyendo reacciones multi-componente y separadas por `.`.
- Datasets: catálogo de reacciones de `tests/integration/test_smirks_reactions.py` + SMIRKS públicos (Daylight/rdkit `test_data`); opcional: "Retrosynthesis templates" del proyecto `rxnmapper` para ampliar.
- Procedimiento: `python scripts/compare_smirks_reactions.py` (existente); ampliar a 100+ plantillas.
- Estado en el repo: `src/smirks.rs`; `tests/test_smirks_reactions.rs` (15 tests); script existente.
- Criterio: coincidencia de productos ≥ 99 % del catálogo; documentar cada divergencia (aromaticidad, H implícitos, estereo).

**F05 · Anillos SSSR (`rings::compute_sssr`, `find_sssr_pub`)**
- Referencia probada: RDKit `Chem.GetSymmSSSR` / `GetSSSR`; Open Babel `OBMol::GetSSSR`.
- Test equivalente: por molécula, comparar **multiset de tamaños de anillo** y pertenencia átomo-a-anillo (invariante, porque el SSSR no es único); además, número de anillos independientes = E − V + C (ciclo rank) como L0.
- Datasets: ChEMBL 10k (filtración por aromaticidad) + casos con anillos fusionados del set de estrés de RDKit (`test_data/ring`); PUBCHEM ring systems (opcional).
- Procedimiento: `scripts/compare_rings.py` → fixture; test `tests/regression/test_rings_cross`.
- Estado en el repo: `src/rings/sssr.rs`, `tests/regression/test_algorithm_validation.rs`.
- Criterio: multiset de tamaños idéntico al de RDKit en ≥ 99.9 % del corpus; rank de ciclos exacto siempre.

**F06 · Aromaticidad y estereocentros (`analyze_graph_features`, `stereo::analyze_stereo`)**
- Referencia probada: RDKit (`Chem.FindMolChiralCenters(includeUnassigned=true, useLegacyImplementation=false)`, `Chem.GetSymmSSSR`, aromaticidad por `Chem.SetAromaticity`), y para nomenclatura R/S/CIP: RDKit `AssignStereochemistry` + **centres** `CIP` (o `openbabel -osmi --append "CIP"`).
- Test equivalente: comparar listas de centros asignados/no asignados, etiquetas R/S/E/Z, y flags aromáticos por átomo/enlace; para estereocentros inferidos comparar solo el conjunto de índices.
- Datasets: ChEMBL 10k; set de casos CIP extremos (patentes de biarilos atropisoméricos, helicenos simples); `tests/fixtures/metal_experimental_geometry.json` para lo metálico.
- Procedimiento: `scripts/compare_stereo.py` → fixture; test `test_stereo_cross`; el repo ya tiene `analyze_graph_features` en la API pública.
- Estado en el repo: `src/stereo.rs`, `src/reactivity.rs` (parte), `tests/regression/test_tet_centers.rs`, `test_algorithm_validation.rs`.
- Criterio: ≥ 99 % de coincidencia en centros etiquetados; documentar diferencias CIP (prioridades) como limitación.

**F07 · Topología de coordinación y hapticidad (`topology::analyze_topology`, `materials::detect_hapticity`)**
- Referencia probada: no existe librería canónica; usar **CSD/COD geometries** como referencia experimental y `RDKit`/`pymatgen` solo para distancias de enlace y coordinación (`CovalentRadius`, `ase.neighborlist.natural_cutoffs`). Para motivos sandía/ferroceno: definiciones IUPAC.
- Test equivalente: sobre complejos metálicos conocidos (ferroceno, Cp*, bis-benceno Cr, metalocenos sustituidos), comprobar número de coordinación, hapticidad (η⁵/η⁶), ángulos y distancias metal-ligando contra geometría experimental.
- Datasets: `tests/fixtures/metal_experimental_geometry.json` (ya existe), `eht_metal_reference.json`; ampliar con **tmQM** (86k complejos TM con geometrías DFT, CSD-derivados) y COD para casos concretos.
- Procedimiento: `scripts/compare_topology_tmqm.py` (nuevo) → fixtures; test `test_topology_cross`.
- Estado en el repo: `src/topology.rs`, `tests/regression/test_metal_experimental_geometry.rs`, `test_eht_metal_references.rs`.
- Criterio: ≥ 98 % de asignación de hapticidad correcta en tmQM (submuestra 1k); distancia M–L ± 0.1 Å.

### 4.2 Geometría 3D y conformeros

**F08 · Matriz de bounds y distance geometry (`distgeom::calculate_bounds_matrix`, `smooth_bounds_matrix`, `triangle_smooth`)**
- Referencia probada: RDKit distance-geometry bounds (`rdDistGeom` interno; extraíble con `rdkit.DistanceGeometry.GetMoleculeBoundsMatrix` vía `AllChem.GetMoleculeBoundsMatrix` en RDKit ≥ 2022.09) y los scripts del repo que ya lo dumpean.
- Test equivalente: misma molécula → misma matriz de cotas superior/inferior elemento a elemento (alineando orden de átomos), tras verificar que las tablas de radios enlace/ángulo coinciden; comparar **violaciones** tras smoothing (máximo y RMS).
- Datasets: ChEMBL 10k; `tests/fixtures/gdb20_reference_1k.json` (geometrías ya resueltas para comprobar que las cotas las respetan).
- Procedimiento: ya existe `scripts/dump_rdkit_bounds.py` + `scripts/compare_bounds*.py` + `scripts/oracle_noff.py`; consolidar en `scripts/generate_bounds_reference.py` → `tests/fixtures/bounds_reference.json.gz`; tests `tests/regression/test_bounds_compare.rs` (existe parcialmente).
- Estado en el repo: `src/distgeom/bounds/*`, `tests/analysis/compare_bounds.rs`, `tests/regression/test_bounds_compare.rs`.
- Criterio: ≤ 1e-3 Å de diferencia en cotas cuando los parámetros de enlace coinciden; violación máxima post-smoothing ≤ 5e-3 Å.

**F09 · Embedding ETKDG (`conformer::generate_3d_conformer`, `embed`, `etkdg.rs`)**
- Referencia probada: **RDKit ETKDGv3** (`AllChem.EmbedMolecule`, `randomSeed` fijo) — el estándar de facto; Open Babel `--gen3d` como segundo voto; `xtb --opt` para validación de mínimos locales.
- Test equivalente: mismas moléculas + mismas semillas → alinear con Kabsch (o mejor: comparar **matriz de distancias internas**, invariante a permutación/orientación) contra ETKDGv3 y contra las referencias ya congeladas; métricas: heavy-atom RMSD, violación de estereocentros, planitud, energía UFF de la pose.
- Datasets: **GDB-20** (20k moléculas pequeñas; el repo usa 1k en `tests/fixtures/gdb20_reference_1k.json` y 300 en `gdb20_ensemble_300.json`), `rdkit_10k_reference.json`, `reference_coords_500.json`, **GEOM** (37M conformeros con energías GFN2/DMC — para ranking), PubChem3D (conformeros para >90 % de PubChem, formato SDF).
- Procedimiento: `scripts/generate_gdb20_reference.py`/`generate_ensemble_reference.py` (existentes) + `scripts/compare_embedding.py`; tests `test_gdb20_rmsd.rs`, `test_gdb20.rs`, `test_geometry_quality.rs`.
- Estado en el repo: pipeline completo + suites GDB-20; ya reporta éxito ≥ 99 % y RMSD vs RDKit.
- Criterio: éxito ≥ 99 % en GDB-20; RMSD mediana ≤ 0.25 Å, P90 ≤ 0.5 Å (§3.3); violación estéreo = 0.

**F10 · Refinamiento por torsiones y ETKDG-3D (`forcefield::etkdg_3d`, `build_etkdg_3d_ff`, `torsion_params`)**
- Referencia probada: RDKit `rdDistGeom` con torsiones ETKDG (`ETversion=2/3`), y los potenciales de torsión extraídos del mismo RDKit (`tests/fixtures/torsion_params.json` ya está en el repo).
- Test equivalente: comparar (a) tabla de parámetros de torsión por SMARTS, (b) energía del término de torsión sobre el mismo confórmero, (c) ángulo de torsión final por enlace rotable contra RDKit.
- Datasets: `tests/fixtures/torsion_params.json`, `scripts/extract_torsions.py`/`extract_torsions_10k.py` (existentes), ChEMBL 10k (extracción de torsiones reales).
- Procedimiento: `python scripts/extract_torsions_10k.py` ya produce la tabla; nuevo `scripts/compare_torsion_energies.py`.
- Estado en el repo: `src/forcefield/etkdg_3d/*`, `tests/regression/test_grad_torsion_detail.rs`, `debug_torsions.rs`.
- Criterio: energía de torsión ± 1e-3 kcal/mol con mismos parámetros; ángulos ± 2°.

**F11 · Validación geométrica del confórmero (`distgeom::validation`)**
- Referencia probada: RDKit `Chem.rdMolTransforms` + `AllChem.EmbedMolecule` con `enforceChirality`; y reglas geométricas de la literatura (CSD medios de enlace/ángulo).
- Test equivalente: L0 invariantes (§3) + sensibilidad: tomar 100 poses buenas y corromperlas (colapso, inversión de quiral, puckering SP2) verificando que cada check las rechaza; y comparar los 5 checks con un "checker" equivalente en RDKit (`Chem.AssignStereochemistryFrom3D`, `GetBestRMS` para confórmeros).
- Datasets: `tests/fixtures/gdb20_ensemble_variance.json` (poses válidas), `tests/fixtures/diverse_molecules.smi`; poses corruptas generadas por script.
- Procedimiento: `tests/regression/test_diverse_molecules.rs` (existe) + ampliar con `test_validation_sensitivity`.
- Estado en el repo: `src/distgeom/validation.rs` (7 funciones), `docs/algorithm/validation.md` describe umbrales.
- Criterio: 100 % de detección en el set corrupto; 0 % de falsos rechazos en GDB-20.

**F12 · Búsqueda de confórmeros y ranking (`conformer::search_conformers_with_uff`, `embed_diverse`, `generate_3d_conformer_best_of_k`)**
- Referencia probada: RDKit `EmbedMultipleConfs` + `MMFFOptimizeMoleculeConfs`/`UFFOptimizeMoleculeConfs`; ranking por energía con la misma funcionalidad; GEOM como verdad de energía (DMC/GFN2).
- Test equivalente: generar N confórmeros con ambos, optimizar con el mismo FF y comparar: (a) nº de mínimos únicos detectados (clustering), (b) energía mínima global alcanzada, (c) correlación de ranking, (d) cobertura del mínimo de referencia.
- Datasets: GEOM (mínimos DMC como ground truth), `tests/fixtures/gdb20_ensemble_300.json`, PubChem3D.
- Procedimiento: nuevo `scripts/compare_ensemble_ranking.py`; tests `tests/regression/test_ensemble_rmsd.rs` (existe) + `benchmarks/test_new_features_bench.rs`.
- Estado en el repo: `src/conformer.rs`, `src/clustering.rs`; ya hay ranking UFF + filtrado RMSD.
- Criterio: ≥ 90 % de moléculas donde el mínimo sci-form ≤ mín(GEOM)+1 kcal/mol; ranking Spearman ≥ 0.8.

**F13 · RMSD, Kabsch y alineamiento (`alignment::compute_rmsd`, `align_coordinates`, `align_quaternion`)**
- Referencia probada: RDKit `rdMolAlign.GetBestRMS` / `AlignMol` (permutación+simmetría aware); SciPy `scipy.spatial.transform.Rotation.align_vectors` para Kabsch puro.
- Test equivalente: (a) L0: RMSD de una molécula consigo misma = 0, invariante ante rotación/traslación/reflexión(no); (b) contra RDKit `GetBestRMS` sobre 1000 pares de confórmeros del mismo grafo (usar `GetSubstructMatches` para mapear); (c) `align_coordinates` vs `align_quaternion` internos (L1).
- Datasets: `tests/fixtures/reference_coords_500.json`, pares de confórmeros GDB-20; cualquier SDF de prueba con hidrógenos explícitos.
- Procedimiento: existe `scripts/compare_embedding.py` (usa RMSD); añadir `scripts/compare_rmsd_rdkit.py`; test `test_rmsd_cross`.
- Estado en el repo: `src/alignment/kabsch.rs`, `src/bin/rmsd_per_mol.rs`, `tests/regression/test_algorithm_validation.rs`.
- Criterio: |ΔRMSD| ≤ 1e-6 Å vs SciPy (mismo mapeo); vs RDKit GetBestRMS ± 1e-4 Å (mismos matches).

**F14 · Clustering Butina y filtro de diversidad (`clustering::butina_cluster`, `filter_diverse_conformers`)**
- Referencia probada: RDKit `rdkit.ML.Cluster.Butina.ClusterData` (mismo algoritmo, cutoffs explícitos); SciPy jerárquico para sanity.
- Test equivalente: misma matriz de distancia (RMSD matrix ya existe: `compute_rmsd_matrix`) → comparar asignación de clusters y centroides; con empates, comparar partición (conjuntos) no etiquetas.
- Datasets: GDB-20 ensembles (10k poses), `tests/fixtures/gdb20_ensemble_300.json`.
- Procedimiento: `scripts/compare_butina.py` → fixture; test `test_butina_cross`; usar `compute_rmsd_matrix` ya disponible.
- Estado en el repo: `src/clustering.rs`; `tests/regression/test_ensemble_rmsd.rs`, `embed_diverse` en `lib.rs`.
- Criterio: partición idéntica salvo empates exactos; centroides dentro de una tolerancia de 1e-9 Å.

**F15 · Energía de deformación (strain) y heurísticas UFF aromáticas (`compute_uff_energy_with_aromatic_heuristics`, docs `strain-energy.md`)**
- Referencia probada: RDKit UFF/MMFF + definición de strain de la literatura (E(pose) − E(mínimo local optimizado)); TorsionScan de RDKit (`rdMolTransforms` + `SetDihedralDeg`) para perfiles.
- Test equivalente: mismo confórmero → comparar strain contra RDKit (mismo FF) y contra el mínimo optimizado con el mismo motor; barrido de torsión 10°–170° paso 10° contra RDKit `SetDihedralDeg` + `MMFFOptimizeMolecule`.
- Datasets: GEOM (conformeros + energía), perfiles de torsión de referencia en `tests/fixtures/torsion_params.json`.
- Procedimiento: `src/forcefield/torsion_scan.rs` + nuevo `scripts/compare_strain.py`.
- Estado en el repo: `src/forcefield/torsion_scan.rs` (10 funciones), `docs/algorithm/strain-energy.md`.
- Criterio: perfil de torsión ± 0.05 kcal/mol por punto; strain ± 0.1 kcal/mol.

### 4.3 Fuerzas y optimización

**F16 · UFF (`forcefield::uff`, `build_uff_force_field`, `compute_uff_energy`)**
- Referencia probada: **RDKit** `AllChem.UFFGetMoleculeForceField` (implementación canónica, BSD-3); Open Babel `OBForceField::UFF` (segundo voto).
- Test equivalente: misma topología + mismas coordenadas congeladas → comparar energía total y, en lo posible, descomposición por términos (bond, angle, torsion, inversion, vdW, electrostático). Después optimizar con ambos motores y comparar geometría final (RMSD y energías).
- Datasets: `tests/fixtures/reference_coords_500.json`, GDB-20 1k (geometrías ya resueltas), ChEMBL 10k (muestra 1k); opcional QM9 (134k geometrías DFT para medir la deformación al relajar con UFF).
- Procedimiento: nuevo `scripts/compare_uff_rdkit.py` (RDKit `UFFGetMoleculeForceField.CalcEnergy()` sobre el mismo molblock) → `tests/fixtures/uff_reference.json.gz`; test Rust `tests/regression/test_uff_cross.rs`; ya existe `tests/regression/test_gradient_check.rs` para L0.
- Estado en el repo: `src/forcefield/uff.rs`, `params.rs`; corrección aromática UFF propia (`apply_aromatic_uff_correction`).
- Criterio: energía total ≤ 0.01 kcal/mol en ≥ 99 % (misma topología) o documentar término divergente; RMSD post-optimización ≤ 0.05 Å.

**F17 · MMFF94 (`forcefield::mmff94`, `compute_mmff94_energy`, `assign_mmff94_type`)**
- Referencia probada: RDKit `AllChem.MMFFGetMoleculeForceField` + `MMFFSanitizeMolecule` (validación oficial MMFF94 de Halgren); Open Babel `MMFF94` como segundo voto.
- Test equivalente: tipado de átomos MMFF (código numérico por átomo) sobre el corpus → comparación exacta; energía total y por términos sobre la misma geometría; gradiente; y optimización (comparar mínimos). Nota: usar la misma versión de MMFF (MMFF94 vs MMFF94s).
- Datasets: `data/chembl_10k_practical.smi` (tipado masivo), GDB-20 (energías), y la **MMFF94 Validation Suite** de Halgren (distribuida con RDKit/Open Babel en `Data/MMFF94_dative.sdf`; ~200 moléculas con energías por término) si está disponible en la distribución instalada.
- Procedimiento: `scripts/compare_mmff94_rdkit.py` → fixture con tipado + descomposición; test `test_mmff94_cross`.
- Estado en el repo: `src/forcefield/mmff94.rs` + `atom_typer.rs`; `tests/regression/test_gradient_check.rs`; `compute_mmff94_energy` público.
- Criterio: 100 % de coincidencia de tipos MMFF en ≥ 99 % de moléculas; energía ≤ 0.01 kcal/mol; términos ± 0.005 kcal/mol.

**F18 · Optimizadores y minimización restringida (`optimization::*`, `bounds_ff::lbfgs/bfgs`, `minimize_bfgs_rdkit`, `optimize_with_constraints`)**
- Referencia probada: **SciPy** `optimize.minimize(method="L-BFGS-B")` sobre la misma energía/gradiente; RDKit `ForceField.Minimize`/`MinimizeMolecule`; geomeTRIC para constrained optimization con coordenadas internas.
- Test equivalente: (a) trayectoria de energía monótona por iteración (L0); (b) mismo punto de partida → comparar energía final y geometría (RMSD) en problemas convexos (resortes armónicos) y en moléculas pequeñas; (c) restricciones: verificar que las coordenadas congeladas no se mueven > 1e-8 Å y que los gradientes proyectados son ~0.
- Datasets: problemas analíticos propios (`src/potentials.rs`: Morse, double-well, Muller-Brown, SN2 1D) — ya existen y son la referencia exacta; moléculas de GDB-20 para lo molecular.
- Procedimiento: `tests/regression/test_gradient_check.rs`, `tests/regression/test_step_by_step.rs`; nuevo `scripts/compare_lbfgs_scipy.py` para los problemas analíticos.
- Estado en el repo: `src/forcefield/minimizer.rs`, `src/forcefield/bounds_ff/{lbfgs,bfgs}.rs`, `src/optimization/mod.rs`, `src/potentials.rs`.
- Criterio: en sistemas analíticos, convergencia al mínimo global con error ≤ 1e-6; en moleculares, E_final ≤ E_ref + 1e-4 kcal/mol.

**F19 · Barrido de torsiones y búsqueda de rotámeros (`forcefield::torsion_scan`)**
- Referencia probada: RDKit (`rdMolTransforms.SetDihedralDeg` + optimización parcial), `AllChem.EmbedMultipleConfs` con `useRandomCoords=false` para enumeración; bibliotecas de rotámeros (Torsion Library de Cambridge/`torsionlib`) para preferencias.
- Test equivalente: mismo dihedro → comparar perfil E(φ) cada 10°, mínimos locales (posición ± 5°), barreras (± 0.1 kcal/mol); y comparar el conjunto de mínimos encontrados por `systematic_rotor_search` vs enumeración RDKit.
- Datasets: `tests/fixtures/torsion_params.json`, librería de torsiones de CSD (resumen público en la "Torsion Library"); para validación de preferencias, `scripts/extract_torsions_10k.py` sobre ChEMBL.
- Procedimiento: test `tests/regression/test_grad_torsion_detail.rs` (existe) + `scripts/compare_torsion_profiles.py`.
- Estado en el repo: `src/forcefield/torsion_scan.rs`, bins `debug_torsions.rs`, `debug_ring.rs`.
- Criterio: posición de mínimos ± 5°, barreras ± 0.1 kcal/mol; cobertura ≥ 95 % de mínimos con energía < 3 kcal/mol.

**F20 · Gradientes analíticos vs diferencias finitas (transversal a todos los motores)**
- Referencia probada: no aplica librería; el patrón estándar es diferencias centrales con paso h y tolerancia relativa; usar `approx` (ya en dev-dependencies).
- Test equivalente: para cada motor con gradiente (UFF, MMFF94, bounds-FF, EHT, PM3, GFN1/2, ReaxFF, HF, DFT, ANI, MLFF): comparar `∇E_analítico` vs `∇E_numérico` (h = 1e-3 Å central; verificar estabilidad con h = 1e-2, 1e-4), y comprobar invariancia traslación/rotación (suma de fuerzas ≈ 0).
- Datasets: geometrías GDB-20 (`reference_coords_500.json`), complejos metálicos (`metal_experimental_geometry.json`), snapshots MD de `dynamics_live`.
- Procedimiento: ya existe `tests/regression/test_gradient_check.rs` y `test_gradient_perterm.rs`; extender con una tabla por motor y feature gate en `tests/regression/test_gradients_matrix.rs`.
- Estado en el repo: suites `test_gradient_check`, `test_gradient_perterm`, `test_grad_isolate`, `test_grad_torsion_detail`.
- Criterio: error relativo ≤ 1e-4 en el 95 % de componentes y ≤ 1e-2 en el resto; norma de fuerza total ≤ 1e-6.

### 4.4 Cargas, electrostática y solvatación

**F21 · Cargas Gasteiger-Marsili (`charges::gasteiger_marsili_charges`, `compute_charges`)**
- Referencia probada: **RDKit** `AllChem.ComputeGasteigerCharges` (mismo algoritmo PEOE; concentraciones/escalados equivalentes) y Open Babel `OBGasteiger`; el artículo original de Gasteiger & Marsili (1980) da los valores canónicos.
- Test equivalente: misma molécula → comparar carga por átomo; verificar suma = carga total; simetrías (benceno: C iguales, H iguales); sensibilidad a número de iteraciones (sci-form usa 6 por defecto; RDKit usa 12 → documentar y comparar la tendencia + versión configurada `gasteiger_marsili_charges_configured`).
- Datasets: ChEMBL 10k (muestra 1k), `tests/fixtures/experimental_reference.json` (moléculas pequeñas con referencia experimental/tablas), pares de moléculas con grupos funcionales definidos.
- Procedimiento: `scripts/generate_charges_reference.py` → fixture; test `tests/regression/test_charges.rs` (existe) ampliado a comparación cruzada.
- Estado en el repo: `src/charges/gasteiger.rs`, tests `test_charges.rs`, `test_session_algorithms.rs` (§10 EEQ damping).
- Criterio: ± 0.02 e vs RDKit con el mismo número de iteraciones; ± 0.05 e con las 6 iteraciones por defecto; suma ≤ 1e-6.

**F22 · Cargas EEQ (`charges_eeq`, `compute_eeq_charges/energy/gradient`)**
- Referencia probada: la implementación EEQ de **xtb/tblite** (Grimme group; `xtb --sp` imprime cargas EEQ en la salida de población) y el artículo EEQ (Caldeweyher et al., 2019).
- Test equivalente: misma geometría → comparar cargas atómicas EEQ, potenciales y derivadas; verificar neutralidad, comportamiento con carga total ±1, simetría en H2O/NH3 y energías de Coulomb EEQ.
- Datasets: moléculas del set de validación de xtb (`xtb` test suite en GitHub), complejos TM de tmQM (submuestra), `tests/fixtures/eht_metal_reference.json`.
- Procedimiento: nuevo `scripts/compare_eeq_xtb.py` (invoca `xtb --sp` con `--verbose` y extrae cargas); test `tests/experimental/test_eeq.rs` (existe) + `test_eeq_cross`.
- Estado en el repo: `src/charges_eeq/*`, `tests/experimental/test_eeq.rs`, `tests/regression/test_session_algorithms.rs` (§10).
- Criterio: carga por átomo ≤ 0.005 e; suma ≤ 1e-6; energía EEQ ≤ 0.1 kcal/mol.
- **Extensión EEQ-BC (2026-10-13)**: port término a término del modelo *bond
  capacitor* 2025 (`charges_eeq::compute_eeqbc_charges`, Froitzheim et al.,
  J. Chem. Phys. 162, 214109) con tablas generadas por
  `scripts/validation/gen_eeqbc_params.py` (multicharge `eeqbc2025` + mctc-lib
  vdW/Pauling). Oráculo **directo**: CLI de multicharge (`--model eeqbc`, build
  meson/ninja), fixture `tests/fixtures/validation/eeqbc_charges.json` (23
  moléculas: orgánicas, iones, metales 3d/4d/5d, lantánidos/actínidos, XeF2,
  ferroceno) vía `scripts/validation/oracle_eeqbc.py`; gate
  `f22_eeqbc_matches_multicharge` con tolerancia 1e-10 e (medido 6.9e-12 e en
  cargas y 1.9e-11 en CN). Nota de fidelidad: `new_eeqbc2025_model` multiplica
  la tabla vdW de mctc-lib por `aatoau` y luego por `autoaa`, de modo que los
  dos factores se cancelan y el modelo recibe el valor crudo de la tabla
  (número en Å) contra coordenadas en Bohr — reproducido literalmente.

**F23 · Análisis de población Mulliken/Löwdin (`population::mulliken_analysis`, `lowdin_analysis`, `compute_population`)**
- Referencia probada: **YAeHMOP** (EHT canónico con poblaciones Mulliken; SourceForge) para el camino EHT; **PySCF** (`mf.mulliken_pop`, `mf.mulliken_meta`; Löwdin vía `scf.hf.SCF` + `pyscf.lo` o manual `S^{1/2}`) para HF/STO-3G; **Multiwfn** para análisis equivalente sobre wfn/molden.
- Test equivalente: misma base y densidad → comparar carga por átomo y matrices de población; Löwdin: definición `S^{1/2}` (verificar la convención empleada); invariantes: suma = carga molecular total, simetría, y dependencia esperada con la base (Mulliken no es base-invariante).
- Datasets: `tests/fixtures/experimental_reference.json`, moléculas H2O/NH3/CH4/CH2O de la suite NIST/CCCBDB; salida de `xtb`/MOPAC como contraste cualitativo.
- Procedimiento: existe `scripts/generate_population_reference.py` (mencionado en TESTING.md); test `tests/regression/test_eht.rs` / extender con `test_population_cross`.
- Estado en el repo: `src/population/population.rs` (Mulliken, Löwdin, bond orders), `tests/regression/test_phase_c.rs`.
- Criterio: Mulliken ≤ 0.005 e (mismo método y base); suma exacta; Löwdin ≤ 0.01 e.

**F24 · NBO y NPA (`population::compute_nbo`, `compute_npa`)**
- Referencia probada: **JANPA** (Natural Population Analysis open-source, Java/Python) y **NBO7/Psi4** (si hay licencia); NPA es insensible a base y muy útil como chequeo cruzado con HF-3c.
- Test equivalente: misma matriz de densidad → comparar cargas NPA por átomo y configuración electrónica natural (s, p, d por átomo); invariantes: suma = carga total, coherencia con momento dipolar.
- Datasets: H2O, NH3, CH4, CH3OH, benceno; geometrías de `tests/fixtures/experimental_reference.json`; ampliar con `scripts/` para promover átomos de las pruebas de `tests/regression/test_experimental_comparison.rs`.
- Procedimiento: nuevo `scripts/janpa_reference.py` (o export molden + JANPA CLI); fixture JSON; test `test_npa_cross`.
- Estado en el repo: `src/population/npa.rs` (implementación propia), sin suite cruzada todavía.
- Criterio: ≤ 0.01 e por átomo; suma ≤ 1e-5.

**F25 · Órdenes de enlace Wiberg/Mayer (`population::compute_bond_orders`)**
- Referencia probada: **Multiwfn** (Wiberg en basis de Löwdin, Mayer con S y P) y PySCF (`mf.mulliken_pop` extendido; Mayer = `P S P`). Para EHT: YAeHMOP (Wiberg/Mulliken overlap populations).
- Test equivalente: misma densidad → matriz de órdenes de enlace por par de átomos; comparar orden relativo y valores numéricos; invariantes: suma de Wiberg ≈ valencia efectiva, simetría de la matriz.
- Datasets: benceno (C─C ≈ 1.4 e en Wiberg-Mulliken para STO-3G/DFT; documentar la base), etileno, acetileno, H2O; complementario: comparar con órdenes de la literatura para sistemas aromáticos.
- Procedimiento: `scripts/generate_bond_orders_reference.py` → fixture; test `test_bond_orders_cross`.
- Estado en el repo: `src/population/population.rs`; ya hay `test_session_algorithms.rs` (bond_orders_parallel_matches_serial).
- Criterio: diferencia ≤ 0.05 e por par; orden relativo idéntico en ≥ 99 % de pares relevantes.

**F26 · Momento dipolar (`dipole::compute_dipole`, `compute_dipole_from_eht`, `compute_dipole_integrals`)**
- Referencia probada: **NIST CCCBDB** y **NIST WebBook** (valores experimentales); PySCF (`mf.dip_moment`) y ORCA para valores *calculados* comparables (misma base); ya hay suite NIST en el repo.
- Test equivalente: misma geometría → comparar |µ| y dirección (ángulo entre vectores ≤ 5°); L0: CH4/CO2 ≈ 0; H2O ~1.85 D; y comparar el dipolo al nivel EHT vs el publicado del mismo método.
- Datasets: `tests/fixtures/experimental_reference.json`; NIST CCCBDB (https://cccbdb.nist.gov/): H2O, CH3OH, CH2O, C6H6, NH3, HF, CO.
- Procedimiento: ya existe `scripts/generate_dipole_reference.py`; tests `tests/regression/test_extended_molecules.rs::test_dipole_moments_vs_nist` (existe).
- Estado en el repo: `src/dipole/dipole.rs`, tests activos vs NIST (11 moléculas).
- Criterio: ≤ 0.1 D o 5 % en ≥ 80 % de la serie NIST; dirección ≤ 5°.

**F27 · Potencial electrostático y rejilla ESP (`esp::compute_esp_grid`, `compute_esp`, `export_cube`)**
- Referencia probada: **PySCF cubegen** (`.cube` con la misma densidad) y **Multiwfn** para superficies de ESP; `compute_esp` del repo usa cargas de Mulliken.
- Test equivalente: mismos puntos de rejilla/espaciado → comparar valores; L0: potencial a gran distancia ≈ multipolo puntual (q/r); simetría molecular; ausencia de NaN; comparar perfiles ESP sobre la superficie de vdW.
- Datasets: cubos de referencia generados con PySCF (`scripts/generate_esp_cube.py` ya descrito en TESTING.md §12), moléculas H2O, CH2O, benceno; `tests/fixtures/*.cube`.
- Procedimiento: `scripts/generate_esp_cube.py` + `scripts/compare_esp_cube.py`; test `test_esp_cross` (referenciado en TESTING.md como Phase C4 candidato).
- Estado en el repo: `src/esp/esp.rs` (export/read cube), `tests/regression/test_phase_c.rs`.
- Criterio: RMSE ≤ 1e-3 a.u. con la misma densidad; ≤ 0.01 si se compara EHT vs HF (documentado como diferencia de método).

**F28 · SASA (`surface::compute_sasa`, `compute_sasa_parallel`)**
- Referencia probada: **FreeSASA** (MIT; CLI/python) con los mismos radios de vdW y radio de sonda; RDKit `rdFreeSASA.CalcSASA`; el algoritmo de Shrake-Rupley original.
- Test equivalente: misma geometría + mismos radios (Bondi) y sonda (1.4 Å) → comparar área total y por átomo; L0: esfera aislada = 4π(r+R)² (caso analítico); convergencia con nº de puntos de esfera.
- Datasets: GDB-20 1k (áreas por molécula), `tests/fixtures/reference_coords_500.json`; para comparación con radios de uni-molecule: `tests/fixtures/diverse_molecules.smi`.
- Procedimiento: `scripts/compare_sasa_freesasa.py` → fixture; test `tests/regression/test_sasa.rs` (existe) + `test_sasa_cross`.
- Estado en el repo: `src/surface/sasa.rs`, `tests/regression/test_sasa.rs`.
- Criterio: ≤ 1 % de diferencia total; ≤ 2 % por átomo; esfera analítica ≤ 1e-6.

**F29 · Solvatación implícita no polar y GB (`solvation`, `solvation_alpb::gb_kernel`)**
- Referencia probada: **OpenMM** `GBSAOBCForce` (GB/OBC), **AmberTools** (`sander` GB) o **xtb** (GBSA/ALPB) para el término electrostático; el término no polar γ·SASA se compara con FreeSASA.
- Test equivalente: misma geometría y mismos parámetros (radios de Born, constante dieléctrica, γ) → comparar ΔG_elec y ΔG_nopol; L0: molécula neutra y apolar → término electrostático ≈ 0; dependencia esperada con ε y con la carga.
- Datasets: set de hidratación con energías experimentales de solvatación (FreeSolv, 642 compuestos, https://github.com/MobleyLab/FreeSolv) y `tests/fixtures/experimental_reference.json` para moléculas pequeñas.
- Procedimiento: nuevo `scripts/compare_gb_openmm.py` (OpenMM GBSAOBCForce con la misma topología/radios) → fixture; test `test_solvation_cross`.
- Estado en el repo: `src/solvation.rs`, `src/solvation_alpb/*`, tests `tests/regression/test_extended_molecules.rs::test_solvation_druglike_molecules`, `tests/experimental/test_alpb.rs`.
- Criterio: componente electrostatico ± 0.05 kcal/mol vs OpenMM (mismos parámetros); correlación con FreeSolv documentada (no gate estricto por ser modelo 3D-RISM-like simplificado).

**F30 · ALPB (solvatación corrección de carga) (`solvation_alpb`, feature `experimental-alpb`)**
- Referencia probada: **xtb** `--alpb solvent` (implementación de referencia de ALPB, Grimme group); `tblite` para GFN2.
- Test equivalente: misma geometría, mismo disolvente (agua, metanol, cloroformo) → comparar término ALPB de energía; con/sin Born radii; series de carga ionizada (acetato, amonio) para ver tendencia.
- Datasets: `xtb` regression set (GitHub grimme-lab/xtb `test/`), set ALPB del paper (Ehlert et al.); moléculas del corpus ChEMBL (submuestra con aniones/cationes).
- Procedimiento: `scripts/compare_alpb_xtb.py` → fixture; `tests/experimental/test_alpb.rs` (existe) + `test_alpb_cross`.
- Estado en el repo: `src/solvation_alpb/{alpb.rs,born.rs}`, tests existentes.
- Criterio: término ALPB ≤ 0.1 kcal/mol o ≤ 1 % (según caso); orden relativo de solvatación estable entre disolventes.

### 4.5 Estructura electrónica

**F31 · EHT: solucionador, matrices S/H y autoestados (`eht::solve_eht`, `build_overlap_matrix`, `build_hamiltonian`)**
- Referencia probada: **YAeHMOP** (extended Hückel canónico, open source, SourceForge) y **pyEHT** (implementación Python del EHT, PyPI); **PySCF** para las integrales de solapamiento STO-3G (la S del EHT con base Slater puede reimplementarse y compararse por primitiva); *Chemical Reviews* / Hoffmann para valores publicados.
- Test equivalente: tres asserts en orden estricto (ver TESTING.md §9): (1) geometría idéntica; (2) matriz S elemento a elemento; (3) matriz H (verificar tabla VSIP y K de Wolfsberg-Helmholz, típicamente 1.75); (4) autoenergías; (5) autoestados comparando |C| y subespacios degenerados.
- Datasets: `tests/fixtures/eht_metal_reference.json` (revisar), fixtures EHT ya citados (`tests/data/h2o_ref.json` por generar con `scripts/generate_eht_reference.py`), y la tabla VSIP del paper original; complejos TM de tmQM para parámetros experimentales.
- Procedimiento: ya descrito en TESTING.md §9 (`scripts/generate_eht_reference.py`, tests `regression/test_eht.rs`, `test_eht_metal_references.rs`); añadir comparación con YAeHMOP (`yaehmop -eig` o `eht`) sobre metano/etano/benceno/ferroceno.
- Estado en el repo: `src/eht/{overlap,hamiltonian,solver,params,basis}.rs`; suite `test_eht.rs`; `eht_metal_reference.json`.
- Criterio: S ≤ 1e-5 rel; H ≤ 1e-5 (misma tabla); HOMO/LUMO ≤ 1e-4 eV; coeficientes |Δ| ≤ 1e-4.
- **Benchmark de metales de transición (2026-10-13)**: revisión de la tabla
  d-block contra la tabla de referencia de **YAeHMOP 3.0.3**
  (`tightbind/eht_parms.dat`, compilación de Alvarez 1993 y los papers del
  grupo Hoffmann) mediante `scripts/validation/check_eht_metal_params.py`.
  Hallazgo: 33 valores (ζ de 4d/5d y algunos VSIP: Ni 3d, Zr, Tc, Cd) no
  coincidían con la fuente citada en la cabecera — corregidos (62 entradas
  normalizadas, 33 con cambio de valor). Gate `f31_metal_params_match_yaehmop`:
  76 filas comparadas con **Δ = 0** (tolerancia 0.02 eV / 0.02 bohr⁻¹).
  Límites documentados: Y, Ag y Hf tienen filas incompletas (placeholders 0) en
  YAeHMOP y conservan nuestros valores; la base es single-zeta, así que se usa
  el primer exponente de las funciones d double-zeta de la referencia.

**F32 · Gradiente EHT y optimización de geometría (`eht::compute_eht_gradient`, `optimize_geometry_eht`)**
- Referencia probada: no hay gradiente EHT estándar en librerías; la referencia es L0 (diferencias finitas) + comparación de geometría optimizada con resultados EHT publicados y con geometrías de referencia ETKDG/DFT.
- Test equivalente: ∇E analítico vs numérico (≤1e-4 rel); optimización de H2O/NH3/CH4 → geometría final vs experimento (±0.02 Å, ±2°) y vs RDKit/ETKDG (tendencia); energía monótona.
- Datasets: `tests/fixtures/experimental_reference.json` (longitudes/ángulos), NIST CCCBDB (geometrías experimentales), `metal_experimental_geometry.json` para el caso TM.
- Procedimiento: test `test_gradient_check.rs` (ya cubre), `src/bin/check_grad.rs`; añadir fixture de geometrías optimizadas de referencia.
- Estado en el repo: `src/eht/gradients.rs`, `tests/regression/test_gradient_check.rs`, docs `extended-huckel-theory.md`.
- Criterio: gradiente ≤ 1e-4 rel; geometría ±0.02 Å / ±2° vs experimento en H2O/NH3/CH4.

**F33 · Rejilla orbital, isosuperficies y marching cubes (`eht::evaluate_orbital_on_grid*`, `marching_cubes`, `eht::volume`)**
- Referencia probada: **scikit-image** `measure.marching_cubes` (Lewiner, referencia estándar de MC) para la misma función escalar; **Multiwfn**/VMD para la forma visual; análisis de volumen con `marching_cubes` de VTK.
- Test equivalente: (a) L0: normalización ∫|ψ|² ≈ 1 sobre una rejilla fina para orbitales atómicos/moleculares sencillos; (b) comparar vértices/triángulos del isosurface con scikit-image sobre el mismo array denso (tras simplificar a nivel de área/volumen: |Δárea| ≤ 1 %); (c) simetría del isosurface (benceno π).
- Datasets: ninguna descarga necesaria; opcional `PubChem3D` para shapes/volúmenes comparables, y VMD para renders de control.
- Procedimiento: exportar la rejilla a `.cube` u `.npy` (ya hay `export_cube`) y ejecutar `scripts/compare_marching_cubes.py` (nuevo); tests `test_eht.rs`.
- Estado en el repo: `src/eht/{volume,marching_cubes}.rs`, `src/gpu/marching_cubes.rs`, `src/mesh.rs`.
- Criterio: área/volumen ≤ 1 % vs scikit-image con el mismo umbral y resolución; normalización ≤ 1e-3.

**F34 · Estructura de bandas y malla k (`eht::compute_band_structure`, `alpha::periodic_linear`, `monkhorst_pack_mesh`)**
- Referencia probada: **pybinding** (tight-binding + band structure, BSD-2) para modelos analíticos; **pymatgen** (`BandStructure`, `Kpoints`) y **spglib** para pesos de malla y caminos k; casos analíticos 1D/2D (cadena de H, grafeno TB).
- Test equivalente: (a) cadena 1D: dispersión analítica `E(k)=ε+2t·cos(ka)` ± 1e-6; (b) grafeno TB: Dirac en K; (c) intervalos/pesos de malla Monkhorst-Pack vs pymatgen; (d) Fermi level por conteo de electrones; (e) invariancia temporal/espacial (hermiticidad del Bloch).
- Datasets: modelos publicados (Slater-Koster) y CIFs de COD para estructuras simples (NaCl, Si); la banda EHT de sistemas 1D con parámetros publicados.
- Procedimiento: `src/bin/benchmark*` no; nuevo `scripts/compare_band_structure.py` con pybinding para la cadena; tests `tests/regression/test_%` en `tests/experimental/test_kpm.rs` complementan el lado espectral.
- Estado en el repo: `src/eht/band_structure.rs`, `src/alpha/periodic_linear/*`, `crates/cli/src/calc_cmds.rs` (charts). 
- Criterio: dispersión analítica ≤ 1e-6 eV; pesos de malla idénticos a pymatgen; hermiticidad ≤ 1e-10.

**F35 · DOS/PDOS y multi-método (`dos::compute_dos`, `compute_pdos*`, `compute_dos_multimethod`)**
- Referencia probada: **Multiwfn** (DOS con ensanchamiento gaussiano estándar) y **PySCF** (eigenvalores exactos) para el caso HF; L0: integral de DOS = número de orbitales; picos ≤ σ en las autoenergías.
- Test equivalente: mismos eigenvalores → comparar curvas DOS elemento a elemento (MSE) y centros de banda; verificar con σ variable; PDOS por elemento/átomo suma al DOS total; para `compute_dos_multimethod`, comparar EHT vs PM3 vs GFN (diferencias esperables; gate = cada uno internamente correcto).
- Datasets: benceno/H2O/CH2O ya cubiertos por fixtures; opcional `PubChemQC PM6` (DOS PM6 de referencia masiva); `tests/fixtures/experimental_reference.json`.
- Procedimiento: `scripts/generate_dos_reference.py` (TESTING.md §13) + `scripts/compare_dos_multiwfn.py`; tests `test_session_algorithms`/`test_new_features`.
- Estado en el repo: `src/dos/*`, `export_dos_json`, `crates/cli` charts.
- Criterio: MSE ≤ 0.01 vs Multiwfn; integral = nº orbitales ± 0.1; PDOS suma al total ± 1e-6.

**F36 · SCF RHF y matrices monoelectrónicas (`scf::run_scf`, `build_overlap_matrix`, `build_kinetic_matrix`, `build_nuclear_matrix`)**
- Referencia probada: **PySCF** (RHF/STO-3G; `mf.mo_energy`, `mf.mo_coeff`) y **ORCA**/Gaussian para contraste; **CCCBDB** para valores HF/STO-3G publicados.
- Test equivalente: misma base (STO-3G) y geometría → comparar S, T, V, J, K, F, energía electrónica, energía total, ε_i, coeficientes (|C|), población Mulliken; convergencia y DIIS (nº iteraciones comparable).
- Datasets: CCCBDB STO-3G (H2O, CH4, NH3, HF, CO, C2H4, CH2O — ya usados en `test_experimental_comparison.rs`); NIST "computational chemistry comparison" para geometrías.
- Procedimiento: `scripts/generate_scf_reference.py` (o ampliar `generate_eht_reference.py`) → fixtures JSON; tests ya presentes en `tests/regression/test_experimental_comparison.rs` (54 tests) + `test_scf_cross`.
- Estado en el repo: `src/scf/*` (basis, integrales, DIIS, mulliken, tipos), `tests/regression/test_experimental_comparison.rs`, `test_alpha_scf_eri_accuracy.rs`.
- Criterio: E_total ≤ 1e-7 Ha; ε_i ≤ 1e-6 Ha; S/T/V ≤ 1e-10 (analíticas); convergencia < 50 iteraciones.

**F37 · Integrales de dos electrones: Obara-Saika, Boys y screening (`scf::gaussian_integrals`, `two_electron`, `obara_saika`, `screening`, feature `alpha-obara-saika`)**
- Referencia probada: **PySCF/libcint** (`mol.intor('int2e')`, `cint1e`), y fórmulas cerradas de Boys (`F0(x)` con series y valores conocidos; benchmark clásico de LibreTexts/`scipy.special` para erf).
- Test equivalente: comparar `F_m(x)` (m=0..4) vs `scipy.special.gamma`/serie de referencia; `eri_ssss` y ERIs (s,s,s,s), (p,s,s,s), (d,p,s,s) vs PySCF sobre las mismas primitivas base; comprobar el Schwarz bound (|(ij|kl)| ≤ Q_ij·Q_kl) y que el screening no altera la energía > 1e-9 Ha.
- Datasets: ninguna descarga; base de prueba H2O/benzene con cc-pVDZ (base descargable si no está) para integrales.
- Procedimiento: `src/bin/check_c_ints.rs`, `print_c_eri.rs`, `print_ch_eri.rs` (existentes, usados para depurar); test `tests/alpha_scf_eri_accuracy.rs` (existe) + `test_eri_cross_pyscf`.
- Estado en el repo: `src/scf/obara_saika.rs`, `gaussian_integrals.rs`, `two_electron.rs`, `screening.rs`; ya hay tests.
- Criterio: F_m ≤ 1e-12 rel; ERIs ≤ 1e-10 rel; screening sin pérdida > 1e-9 Ha.

**F38 · UHF/ROHF y contaminación de espín (`scf::run_uhf`, `run_rohf`, `compute_uhf_configured`)**
- Referencia probada: **PySCF** `scf.UHF`/`scf.ROHF`, y la teoría de Szabo-Ostlund para valores de ⟨S²⟩ (H2 triplet = 2.0, doblete = 0.75).
- Test equivalente: H2 singlete/doblete/triplete, ROHF H2⁺ → comparar energías con PySCF; ⟨S²⟩ exacto; simetría α/β; separar orbitales α/β y verificar que la solución restringida coincide con RHF cuando corresponde.
- Datasets: moléculas de la suite `test_session_algorithms.rs` §1; CCCBDB para H2/H2⁺.
- Procedimiento: tests existentes `test_session_algorithms.rs` (6 tests UHF/ROHF) + ampliar a PySCF.
- Estado en el repo: `src/scf/uhf.rs`, `src/hf` re-exports; tests §1.
- Criterio: E ≤ 1e-6 Ha vs PySCF; ⟨S²⟩ ≤ 1e-6 del valor exacto; doblete H2⁺ ≈ −0.5 × RHF.

**F39 · Transformación AO→MO de 4 índices (`hf::ao_to_mo_transform(_parallel)`)**
- Referencia probada: **PySCF** `ao2mo.full`/`general`; **Psi4** como segundo voto; la permutación de índices debe documentarse.
- Test equivalente: misma C y mismas ERIs AO → comparar tensor MO (hasta 1e-10) con convenciones fijadas (Chemist vs Physicist notation); invariante: traza/normas y energía MP2-like de prueba.
- Datasets: fixtures AO ERIs generados por PySCF (H2O/STO-3G); no requiere descarga externa.
- Procedimiento: `tests/regression/test_session_algorithms.rs` §3 (existencia) → ampliar a comparación numérica con PySCF (`scripts/generate_ao2mo_reference.py`).
- Estado en el repo: `src/hf/transform.rs` (o similar; API `ao_to_mo_transform`), tests §3.
- Criterio: ≤ 1e-10 abs por elemento con la misma convención.

**F40 · HF-3c (HF + D3 + gCP + SRB) (`hf::solve_hf3c`, `compute_d3_energy`, `compute_gcp`, `compute_srb`)**
- Referencia probada: **ORCA** (implementación de referencia del método HF-3c de Grimme; `! HF-3c`), **xtb** no implementa HF-3c; el paper original (Sure & Grimme 2013) incluye el set de validación (S22/S66, conformeros, TMs).
- Test equivalente: misma geometría → comparar E_HF, E_D3, E_gCP, E_SRB y E_total componente a componente con ORCA HF-3c; geometrías optimizadas vs ORCA (RMSD); y comparación con CCCBDB para HF/def2-SVP-like si se separa la base.
- Datasets: set de validación del paper HF-3c (complejos no covalentes S22/S66, conformeros de alcanos, y "TM complexes" del propio paper); `tests/fixtures/experimental_reference.json` para lo pequeño.
- Procedimiento: los números actuales de HF-3c vs "experimental SCF" están explícitamente marcados como no válidos en `quantum-validation-roadmap.md`; crear `scripts/generate_hf3c_orca_reference.py` (ORCA) → fixture; tests ya presentes en `test_experimental_comparison.rs` (bloque HF-3c).
- Estado en el repo: `src/hf/api.rs`, `tests/regression/test_experimental_comparison.rs`, `tests/fixtures/experimental_reference.json`.
- Criterio: E_total ≤ 1 % o ≤ 10 mHa vs ORCA HF-3c en la misma base; componentes individuales ≤ 1e-6 Ha vs `dftd3`/gCP de referencia; **no** usar como gate la comparación actual HF-3c vs "SCF experimental" (definiciones incompatibles).

**F41 · Estados excitados CIS/CISD (`hf::compute_cis`, `compute_cisd`, `oscillator_strength`)**
- Referencia probada: **PySCF** `tdscf.rhf.TDA` (CIS) y `ccsd`/Function-based CISD como segundo voto; ORCA `CIS`.
- Test equivalente: misma base → comparar energías de excitación (≤ 0.01 eV), fuerzas de oscilador (≤ 1e-3) y orden λ; el primer estado del H2/etileno comparado con la referencia; regla de suma de osciladores.
- Datasets: `tests/fixtures/experimental_reference.json`; NIST para excitaciones experimentales de H2/CH2O/etileno; opcional "Thiel set" (conjuntos de validación de estados excitados publicados) para ampliar.
- Procedimiento: `test_experimental_comparison.rs` (spectroscopy_comparison, 5 tests) + `scripts/generate_tda_reference.py` (PySCF TDA).
- Estado en el repo: `src/hf/cis.rs`(aprox.), tests bloque spectroscopy.
- Criterio: ΔE ≤ 0.01 eV vs PySCF misma base; orden de estados del bloque de valencia correcto.
- **Extensión CISD (2026-10-13)**: `hf::cisd::compute_cisd` reescrito como CISD
  determinantal (orbitales de espín) sobre la referencia RHF: espacio con
  dobles de mismo orbital (`(iα iβ)→(aα aβ)`, el determinante dominante en
  bases mínimas), reglas de Slater–Condon con las integrales monoelectrónicas
  recuperadas de la condición de Fock canónica
  (`h_pq = ε_p δ_pq − Σ_j (2(pq|jj) − (pj|jq))`) y energía de correlación
  `λ_min − E_ref` (antes: el menor autovalor CIS, que no es una energía de
  correlación). El espacio coincide exactamente con el de `pyscf.ci.UCISD`
  (141/561/610/1425/3241 determinantes para H2O/CH4/CO/CH2O/C2H4) y el gate
  `f41_cisd_matches_pyscf` reproduce las correlaciones de `pyscf.ci.cisd.CISD`
  a ≤ 2.1e-8 Ha (H2, H2O, CH4, NH3, HF, CO, C2H4, CH2O; tolerancia 1e-6); los
  dos espacios mayores (CH2O/C2H4) solo se ejercitan en `--release` porque el
  build de la matriz en debug es ~50× más lento.
  La matriz se verifica además contra una referencia de segunda cuantización
  explícita (`tests/cisd_brute.rs`, diferencia máxima 5e-15). Camino
  perturbativo (CIS+MP2) validado contra `pyscf.mp.MP2` (H2O −0.0355456516 Ha).

**F42 · PM3 (`pm3::solve_pm3`, `compute_pm3_gradient`, `is_pm3_supported`)**
- Referencia probada: **OpenMOPAC/MOPAC2016** (PM3 original de Stewart; `PM3` keyword) — la implementación canónica; **PubChemQC PM6** (200M+ compuestos calculados con PM6) como dataset masivo de calibración indirecta; **Antechamber** como tercer voto.
- Test equivalente: misma geometría → comparar ΔHf (kcal/mol), energías orbitales, cargas Mulliken y gradiente (o geometría optimizada) con MOPAC; series: H2O, CH4, NH3, HF, CO, C2H4, CH2O, benceno, etanol, aspirina.
- Datasets: PubChemQC (subconjunto CID→PM6; https://nakatamala.github.io/pubchemqc/), set de moléculas del paper PM3 (Stewart 1989: 800+ compuestos con ΔHf publicado); geometrías de `reference_coords_500.json`.
- Procedimiento: `scripts/compare_pm3_mopac.py` (nuevo; `mopac` con `PM3 1SCF`) → fixture; tests `test_extended_molecules.rs` + `test_pm3_cross`.
- Estado en el repo: `src/pm3/*`, tests §5 de `test_session_algorithms.rs` (gaussianas PM3), `tests/fixtures/experimental_reference.json`.
- Criterio: ΔHf ≤ 0.5 kcal/mol (misma geometría) o ≤ 1.5 kcal/mol optimizada; cargas ± 0.05 e; gap orden correcto.

**F43 · GFN0/xTB (`xtb::solve_xtb`, `compute_xtb`, `compute_xtb_gradient`)**
- Referencia probada: **xtb** oficial (binary, LGPL-3; `xtb mol.xyz --gfn 0 --sp`), **tblite** (biblioteca GFN de referencia, LGPL-3) y el set de regresión de xtb en GitHub.
- Test equivalente: misma geometría → energía total, HOMO/LUMO, cargas Mulliken, gradiente numérico vs analítico, y optimización (comparar geometrías); convertir unidades (sci-form eV vs xtb Ha/eV según salida).
- Datasets: `xtb/test` (repo), S66/WATER27 (interacciones no covalentes), moléculas de `test_experimental_comparison.rs`; tmQM para TMs.
- Procedimiento: `scripts/compare_xtb_tblite.py` (nuevo) o invocación directa `xtb --sp`; tests `test_extended_molecules.rs`, `test_session_algorithms.rs` §6 (Broyden).
- Estado en el repo: `src/xtb/{solver,broyden}.rs`, tests existentes; documentación advierte que `solve_xtb` es GFN0-style SCC. Desde 2026-10-13 el solver validado es `xtb::gfn0::solve_gfn0` (port término a término, ver §14); `solve_xtb` sobrevive únicamente como proveedor del gradiente analítico de los backends de fuerza `"xtb"` (NEB/MD, `compute_xtb_gradient`) hasta que existan gradientes GFN0/GFN2 analíticos.
- Criterio: E ≤ 1 kcal/mol vs xtb GFN0; gradiente ≤ 1e-3 eV/Å; convergencia determinista.

**F44 · GFN1-xTB (`xtb::solve_gfn1`)**
- Referencia probada: **tblite** (GFN1 nativo) y **xtb** (`--gfn 1`); paper GFN1-xTB (Grimme 2016) con set de validación (S66x8, WATER27, HAL59, mor41).
- Test equivalente: misma geometría → energía total, cargas *shell-resolved*, dipolos/multipolos, gradiente y optimización; benchmarks no covalentes (S66x8: energía de interacción vs referencia CCSD(T)) como validación L4.
- Datasets: S66x8, WATER27, HAL59, MOR41 (todos públicos, descargables); `xtb` test suite; tmQM (submuestra metálica).
- Procedimiento: `scripts/compare_gfn1_tblite.py` → fixture; tests `test_extended_molecules.rs` + nuevo `test_gfn1_cross`.
- Estado en el repo: `src/xtb/gfn1.rs` (capas sobre GFN0; roadmap exige separar parámetros — ver `quantum-validation-roadmap.md`).
- Criterio: E ≤ 0.5 kcal/mol vs tblite GFN1 en el corpus pequeño; orden de energías de interacción de S66x8 correcto en ≥ 90 %; cargas ± 0.01 e.

**F45 · GFN2-xTB (`xtb::solve_gfn2`)**
- Referencia probada: **tblite** (`GFN2-xTB`) y **xtb** (`--gfn 2`); paper GFN2 (Bannwarth 2019) con sets de validación; **dftd4** para el término de dispersión interno D4 (GFN2 lo lleva embebido: verificar que no se cuenta dos veces).
- Test equivalente: misma geometría → energía (componentes: repulsión, electrónica, dispersión), cargas Mulliken/shell, multipolos, gradiente; optimización (RMSD vs xtb); S66x8 interacción.
- Datasets: S66x8, WATER27, HAL59, MOR41, tmQM; `xtb/test`; complejos halogenados del paper GFN2.
- Procedimiento: `scripts/gfn2_compare_energy.py`, `gfn2_component_compare.py`, `gfn2_matrix_compare.py`, `gfn2_check_occupations.py`, `bench_gfn2_100.py` (todos ya existen en `scripts/` como scripts de depuración) → consolidar en un solo `scripts/compare_gfn2_tblite.py` + fixture versionado.
- Estado en el repo: `src/xtb/{gfn2,gfn2_params,d4,d4_data,gradients}.rs`; muchos scripts ad-hoc.
- Criterio: E_total ≤ 0.5 kcal/mol vs tblite (objetivo fino 1e-6 Ha); componentes ≤ 1e-3 Ha; cargas ± 0.01 e; gradiente ≤ 1e-3 eV/Å.

**F46 · Dispersión D3/D4 (`dispersion::compute_d4_energy`, `compute_d4_gradient`, `compute_d3_energy`)**
- Referencia probada: **dftd4** (implementación de referencia oficial de Grimme, CLI + Python `dftd4`, LGPL) y **simple-dftd3/dftd3** para D3.
- Test equivalente: misma geometría, misma carga, mismos parámetros (funcional para D3: PBE/B3LYP; D4 con los niveles configurados) → comparar E (kcal/mol o Ha) y gradientes; verificar C6/C8/coordinación por átomo contra dftd4 (`--verbose`).
- Datasets: S66/S66x8 y WATER27 (energías de interacción), y el propio test suite de dftd4; complejos de tmQM para metales.
- Procedimiento: `scripts/gfn2_compare_energy.py` (contiene comparaciones D4) + nuevo `scripts/compare_dftd4.py` → fixture; test `tests/experimental/test_d4.rs` (existe) + `test_d4_cross`.
- Estado en el repo: `src/dispersion/*`, `src/gpu/d4_dispersion_gpu.rs`, tests `test_d4.rs`.
- Criterio: E ≤ 1e-3 kcal/mol; gradiente ≤ 1e-4 kcal/mol/Å; C6 ≤ 1 % vs dftd4.

**F47 · Reactividad electrónica: frontera, Fukui, descriptores duales (`reactivity::compute_frontier_descriptors`, `compute_fukui_descriptors`, `rank_reactivity_sites`)**
- Referencia probada: **Multiwfn** (Fukui condensados y dual descriptor son una función principal de la herramienta, con implementación documentada) y **PySCF** para el cálculo de poblaciones que alimentan el esquema (f+, f−, f0); criterios de la literatura (Parr-Yang, Morell).
- Test equivalente: mismo sistema y misma definición de población → comparar f+, f−, f0 y el dual por átomo; L0: cada función de Fukui integra 1 (∫f± dV ≈ 1), f0 = (f+ + f−)/2 y Δf = f+ − f−. Estos descriptores derivan de densidades de N−1/N/N+1; la comparación estricta exige reutilizar las mismas densidades.
  - Para EHT (no hay densidades N±1 estándar), comparar **ranking de sitios** contra Multiwfn/DFT sobre una serie con reactividad conocida (anilina, fenol, nitrosobenceno, piridina) y las reglas de orientación experimental (ortonormal/para/meta).
- Datasets: serie de sustitución aromática clásica (nitración, halogenación) con orientación experimental conocida; artículos de Morell (dual descriptor); opcional "Mayr reactivity" dataset.
- Procedimiento: `scripts/compare_fukui_multiwfn.py` (para la parte accesible con molden) + test de ranking `test_reactivity_ranking_cross`; tests existentes en `test_session_algorithms`/`test_new_features`.
- Estado en el repo: `src/reactivity.rs`, API pública `compute_frontier_descriptors`, `compute_fukui_descriptors`, `compute_reactivity_ranking`.
- Criterio: ranking de sitios correcto en ≥ 80 % de la serie aromática; L0 (nº electrones de frontera, dual) exacto.

### 4.6 Espectroscopía

**F48 · Análisis vibracional, Hessiano y termoquímica (`ir::compute_vibrational_analysis`, `compute_semianalytical_hessian`, `compute_numerical_hessian`, `ir::hessian`)**
- Referencia probada: **MOPAC** (`FORCE` para Hessiano PM3, implementación canónica del método) y **xtb** (`--hess`, Hessiano numérico GFN) para la parte semi-empírica; **ORCA/Gaussian** para DFT; **CCCBDB** y **NIST WebBook** para frecuencias experimentales y termoquímica (ZPVE, S, H, G) de gas.
- Test equivalente: misma geometría y método → comparar modos (emparejados por solapamiento de autovectores ≥ 0.99, no por índice), frecuencias, ZPVE, entropía/entalpía a 298 K; L0: 3N−6 modos, sin frecuencias imaginarias en mínimos, y regla de Teller-Redlich para simetría.
- Datasets: `tests/fixtures/experimental_reference.json` (IR frequencies), CCCBDB (frecuencias armónicas B3LYP/6-31G* y experimentales), NIST WebBook (IR gas), QM9 (frecuencias DFT para 134k moléculas).
- Procedimiento: `test_session_algorithms`/`test_extended_molecules.rs::test_ir_spectrum_extended_molecules` (existen); nuevo `scripts/compare_hessian_mopac.py` (MOPAC `FORCE`) y `scripts/compare_hessian_xtb.py` (`xtb --hess`).
- Estado en el repo: `src/ir/{hessian,vibrations}.rs`, `src/spectroscopy/hessian.rs`; métodos EHT/PM3/xTB/UFF.
- Criterio: frecuencias misma-método ≤ 2 cm⁻¹; vs experimento ± 20 cm⁻¹ (documentar factor de escala si se aplica); ZPVE ≤ 1 kcal/mol; termoquímica ≤ 0.5 kcal/mol.

**F49 · Espectro IR y asignación de picos (`ir::compute_ir_spectrum*`, `assign_peaks`)**
- Referencia probada: **NIST WebBook** (espectros IR experimentales, picos/funciones) y **CCCBDB**; para la asignación por grupo funcional: tablas de correlación IR estándar (Socrates, "Infrared Characteristic Group Frequencies") y el proyecto **nmrshiftdb2/`ir.data`** no aplica — usar el "NIST IR database" y la librería de asignaciones de `openbabel`.
- Test equivalente: mismo conjunto de frecuencias/intensidades → comparar la posición de bandas por grupo (C=O 1650–1750, O-H 3200–3600, N-H 3300–3500, C≡N 2200–2260) y el patrón del espectro ensanchado (correlación); L0: no negatividad de la lorentziana, integral proporcional a la intensidad.
- Datasets: NIST WebBook (espectros de ~10 moléculas de prueba), `tests/fixtures/experimental_reference.json`.
- Procedimiento: `scripts/compare_ir_peaks.py` (nuevo) + tests `test_spectroscopy.rs` (existe), `test_extended_molecules.rs`.
- Estado en el repo: `src/ir/vibrations.rs`, `src/ir/peak_assignment.rs` (`assign_peaks` público), tests de espectroscopía.
- Criterio: asignación correcta de los 5 grupos funcionales clave en ≥ 90 %; posición de banda de grupo dentro de la ventana documentada.

**F50 · UV-Vis sTDA y osciladores (`spectroscopy::compute_stda`, `compute_stda_uvvis`, `reactivity::compute_uv_vis_like_spectrum`)**
- Referencia probada: **`stda`** de Grimme (implementación canónica del sTDA, GitHub grimme-lab/stda, con `xtb`/`orca` como motores de orbitales); **TD-DFT en ORCA/PySCF** como método independiente; **NIST WebBook** (UV-Vis gas) y **PhotochemCAD** (espectros en disolución) para lo experimental.
- Test equivalente: mismos orbitales (EHT, PM3, GFN2) → aplicar sTDA/`stda`: comparar energías de excitación (eV) y fuerzas de oscilador de las N primeras transiciones; después, comparar λmax del espectro ensanchado vs TD-DFT (misma geometría) y vs experimento.
- Datasets: set de validación del paper sTDA (30–50 moléculas con λmax y f publicados); benceno, naftaleno, antraceno, cumarina, rodamina (PhotochemCAD), NIST WebBook gas.
- Procedimiento: `tests/regression/test_extended_molecules.rs::test_legacy_stda_aromatic_series` (existe) + `scripts/compare_stda_grimme.py` (invoca `stda` con los mismos MOs) + `scripts/compare_tddft_orca.py`.
- Estado en el repo: `src/spectroscopy/stda_uvvis.rs`, `src/reactivity.rs`, `test_spectroscopy.rs`.
- Criterio: ≤ 0.3 eV por transición vs TD-DFT; λmax ≤ 30–50 nm vs experimento en familias conjugadas; f > 0 y orden correcto.

**F51 · Fluorescencia y coeficiente de Einstein A (`spectroscopy::compute_fluorescence_spectrum`, `einstein_a_coefficient`)**
- Referencia probada: **ORCA** (TD-DFT + fluorescencia/rates) como cálculo de referencia; valores experimentales de tiempos de vida radiativos y rendimientos cuánticos de estándares (quinina en H2SO4, rodamina 6G, antraceno, fluoresceína) publicados por **IUPAC/PhotochemCAD**.
- Test equivalente: mismos estados sTDA → verificar E_emisión ≤ E_absorción (Stokes), forma de línea normalizada a 1, y A_21 vs la fórmula analítica `A = (2/3)·ω³·|µ|²/(4πε₀ℏc³)` implementada de forma independiente; comparar `τ_rad = 1/A` con los estándares.
- Datasets: tabla de estándares de fluorescencia (PhotochemCAD, IUPAC lifetime standards); NIST opcional.
- Procedimiento: `test_spectroscopy.rs` (existente) + `scripts/check_einstein_a.py` (fórmula cerrada) y comparación con ORCA.
- Estado en el repo: `src/spectroscopy/fluorescence.rs`, `docs/algorithm/fluorescence.md`.
- Criterio: fórmula A ≤ 1 % vs implementación independiente; λ_em > λ_abs; τ_rad dentro de un factor 2–3 de los estándares (método exploratorio).

**F52 · GIAO NMR: shieldings y desplazamientos (`spectroscopy::compute_nmr_shieldings`, `compute_giao_nmr`, `predict_chemical_shifts`)**
- Referencia probada: **ORCA** GIAO (`! RHF/GIAO`, misma base), **PySCF** (con `pyscf.prop` GIAO si está disponible) o NWChem; para lo experimental: **SDBS**, **nmrshiftdb2**, **NMReDATA**; set de referencia de shieldings **NS372** (Gauss) para comparación absoluta de tensores.
- Test equivalente: misma base y geometría → comparar tensor de shielding (componentes y anisótropo) vs ORCA; después, comparar desplazamientos relativos a TMS/referencia interna y contra experimento (MAE por núcleo).
- Datasets: NS372 (38 moléculas, 372 shieldings, valores CCSD(T) publicados), nmrshiftdb2 (descarga de shifts), SDBS (espectros individuales), NMReDATA (dataset abierto de asignaciones).
- Procedimiento: `tests/regression/test_session_algorithms.rs`/`test_extended_molecules.rs::test_giao_nmr_vs_legacy_heteroaromatics` (existen, ignorados); crear `scripts/generate_orca_giao_reference.py` (ORCA sobre las mismas coordenadas) → fixtures.
- Estado en el repo: `src/spectroscopy/giao_nmr.rs`, `src/nmr/shifts.rs`; suite existente parcial.
- Criterio: componente de shielding ≤ 0.5 ppm vs ORCA (misma base; con bases mínimas comparar tendencias, no valores absolutos); MAE ¹H ≤ 0.3 ppm, ¹³C ≤ 3 ppm vs experimento.

**F53 · HOSE codes y predicción empírica de shifts (`nmr::generate_hose_codes`, `predict_shift_from_hose`, `predict_nmr_shifts*`)**
- Referencia probada: **nmrshiftdb2** (base de datos y motor de predicción HOSE público), **CSEARCH/NMRPredict (Robien)** como referencia comercial/documentada, y **ChemDraw** como contraste de usuario; la definición original de HOSE (Bremser 1978).
- Test equivalente: misma molécula y mismo nivel de esfera → comparar el código HOSE generado (string canónico) contra nmrshiftdb2/CSEARCH y el shift predicho contra el experimental; métrica MAE/RMSE por núcleo, con y sin ensemble conformacional.
- Datasets: **nmrshiftdb2** (descarga en bloque, ~50k moléculas, cientos de miles de shifts), **SDBS** (~34k espectros), **NMReDATA** (publicaciones + registros), dataset HOSE de ejemplo del repo (`src/nmr` tests).
- Procedimiento: `scripts/compare_hose_nmrshiftdb.py` (nuevo) → fixture + informe; tests `test_spectroscopy.rs`, `test_extended_molecules.rs::test_legacy_nmr_*`.
- Estado en el repo: `src/nmr/{hose,shifts,spectrum}.rs`; ya valida benceno ¹H≈7.27, ¹³C≈128.5 ppm.
- Criterio: MAE ¹H ≤ 0.3 ppm, ¹³C ≤ 3 ppm en un test de 500 moléculas; 100 % de moléculas con código HOSE bien formado.

**F54 · Acoplamientos J (Karplus y promediado de ensemble) (`nmr::predict_j_couplings`, `ensemble_averaged_j_couplings`)**
- Referencia probada: **nmrshiftdb2** (acoplamientos experimentales), **ORCA** (cálculo de acoplamientos Spin-Spin, referencia computacional), y las parametrizaciones clásicas de Karplus (Karplus 1959; Haasnoot-de Leeuw-Altona 1980) implementadas en `pyKarplus`-like o en librerías de predicción.
- Test equivalente: mismo confórmero → comparar ³J(φ) con la curva Karplus de referencia (parámetros publicados) y con ORCA; después, promediar sobre el ensemble y comparar con los J experimentales (etano, ciclohexano, norbornano, péptidos).
- Datasets: nmrshiftdb2 (J experimentales), `tests/fixtures/experimental_reference.json`; moléculas rígidas de referencia (norbornano) para aislar el efecto de promediado.
- Procedimiento: `tests/regression/test_session_algorithms.rs` §7 (naphthalene peri-H 5J, existente) + `scripts/compare_j_couplings_nmrshiftdb.py`.
- Estado en el repo: `src/nmr/coupling.rs` (aprox.), test §7.
- Criterio: ³J ≤ 1 Hz vs experimento en moléculas rígidas; forma de la curva Karplus correcta (máx ~0°/180°, mín ~90°).

**F55 · Espectro NMR ensanchado, formato y multi-núcleo (`nmr::compute_nmr_spectrum*`)**
- Referencia probada: **nmrglue** (procesado/simulación de espectros, MIT-license) para verificar el formato de la forma de línea; **TopSpin**/MestReNova para comparación de visor; la práctica estándar IUPAC para referenciación (TMS 0 ppm, unfolding de acoplamientos).
- Test equivalente: misma lista de shifts/J → reconstruir el patrón multiplete y comparar con nmrglue (simulación de Voigt/Lorentz) y con el espectro experimental de SDBS (patrón de picos).
- Datasets: SDBS (figura/peaks), nmrshiftdb2 (picos asignados).
- Procedimiento: `src/nmr/spectrum.rs` tests + nuevo `scripts/compare_nmrglue.py`.
- Estado en el repo: `src/nmr/spectrum.rs`, `compute_nmr_spectrum_for_shifts`.
- Criterio: forma de línea normalizada a 1 ± 1e-6; posición de multiplete ≤ 0.01 ppm de la lista; extensión de unidades correcta (Hz↔ppm por frecuencia del instrumento).

### 4.7 ML, descriptores y potenciales neuronales

**F56 · Descriptores 2D y reglas de fármaco (`ml::compute_descriptors`, `compute_ml_descriptors`, `predict_ml_properties`)**
- Referencia probada: **RDKit Descriptors** (`MolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `FractionCSP3`, `RingCount`) y **mordred** (implementación académica de descriptores, BSD); Lipinski/Veber según definición original.
- Test equivalente: misma molécula → comparación exacta de cada descriptor cuya definición coincida; para MW/átomos calcular también desde la tabla periódica (L0); para Lipinski comparar el conteo de violaciones.
- Datasets: `tests/fixtures/diverse_molecules.smi`, `data/chembl_10k_practical.smi`; para propiedades: `MoleculeNet` (ESOL 1.1k, Lipophilicity 4.2k, FreeSolv 642, BBBP, Tox21).
- Procedimiento: `scripts/compare_descriptors_rdkit.py` → fixture; tests `test_extended_molecules.rs::test_ml_logp_extended_molecules`, `test_ml_lipinski_druglike_molecules` (existen).
- Estado en el repo: `src/ml/descriptors.rs`, `src/ml/mod.rs`, tests existentes.
- Criterio: exacto (rel 1e-9) en MW, HBD/HBA, TPSA (si la definición es idéntica); logP Crippen ≤ 0.1; Lipinski: misma clasificación en ≥ 95 %.

**F57 · Descriptores 3D: WHIM, GETAWAY, RDF, TPSA (`ml::compute_whim`, `compute_getaway`, `compute_rdf`, `compute_3d_descriptors`)**
- Referencia probada: **mordred** (`Whim`, `Getaway`, `RDF` — implementación de referencia de Todeschini; **PaDEL-Descriptor** (Java) y **CDK** como segundos votos). WHIM/GETAWAY están definidos explícitamente en el libro de Todeschini & Consonni.
- Test equivalente: misma geometría canónica y mismo peso (masa/volumen/electronegatividad) → comparar valores de cada descriptor no nulo; verificar sensibilidad al peso y a la orientación (debe ser invariante tras Kabsch).
- Datasets: `reference_coords_500.json` (geometrías), GDB-20 (más variedad); ESOL/FreeSolv para validar utilidad predictiva (correlación con logS).
- Procedimiento: `scripts/compare_mordred.py` (instalar `mordred`; genera JSON por molécula) → fixture; test `test_3d_descriptors_cross`.
- Estado en el repo: `src/ml/whim.rs`, `src/ml/getaway.rs`, `src/ml/rdf_descriptors.rs`; `compute_3d_descriptors`.
- Criterio: ≤ 1e-6 rel donde mordred y sci-form usan la misma definición; documentar cualquier definición propia divergente.

**F58 · Huella farmacofórica (`ml::compute_pharmacophore_fingerprint`, `pharmacophore_tanimoto`)**
- Referencia probada: **RDKit `rdkit.Chem.Pharm2D`** (huellas 2D de pares de features), **`align-it`/`pharmit`** (alineamiento farmacofórico) y **`Open3DAlign`**; la definición de features (donor, acceptor, aromático, hidrofóbico, +/−) se documenta en el papers de RDKit y en el proyecto "Pharmacophore modeling" (IUPAC).
- Test equivalente: misma molécula → comparar el conjunto de features detectadas (tipo + índice de átomo) con RDKit (`ChemicalFeatures.BuildFeatureFactory` + `GetFeaturesForMol`); comparar el número de features y los tipos; después, comparar la matriz de similitud sobre una serie con pares conocidos (activos/inactivos de la familia triptasa).
- Datasets: DUD-E (activistas/inactivos por familia, con SMILES), ChEMBL (decoy generation no necesaria), `rdkit.Chem.Pharm2D` tests.
- Procedimiento: `scripts/compare_features_rdkit.py`; test `test_pharmacophore_cross`.
- Estado en el repo: `src/ml/pharmacophore.rs`(aprox.), utilidades en `ml/`.
- Criterio: tipos de feature coincidentes en ≥ 90 % (elasticidad documentada por definición); ranking de similaridad coherente con RDKit (Spearman ≥ 0.8).

**F59 · ECFP/Morgan y Tanimoto (`rings::compute_ecfp`, `compute_tanimoto`)**
- Referencia probada: **RDKit** (`AllChem.GetMorganFingerprintAsBitVect(radius=2, nBits=2048)` = ECFP4; `TanimotoSimilarity` / `BulkTanimotoSimilarity`); CDK como segundo voto.
- Test equivalente: mismos invariantes (identificador de átomo por elemento/grado/valencia/carga/aro) → comparar el **bit vector completo** (XOR = 0) o, si los invariantes difieren por diseño, la matriz de Tanimoto sobre 1000 moléculas (Δ ≤ 1e-9) y el orden de vecinos más próximos.
- Datasets: ChEMBL 10k (bloques ya calculables), `tests/fixtures/*.smi`; tareas MoleculeNet con ECFP (BBBP/Tox21) para validar utilidad.
- Procedimiento: `scripts/compare_ecfp_rdkit.py` → fixture; tests `test_extended_molecules.rs::test_fingerprint_tanimoto_similarity_series` (existe) + `test_ecfp_bits_cross`.
- Estado en el repo: `src/rings/ecfp.rs`; API `compute_ecfp`, `compute_ecfp_batch`, `compute_tanimoto`.
- Criterio: XOR de bits = 0 con los mismos invariantes/parámetros; si no, Tanimoto ≤ 1e-9 de diferencia y Spearman ≥ 0.99 en el ranking de similitud.

**F60 · Random Forest y Gradient Boosting (`ml::train_random_forest`, `train_gradient_boosting`, `predict*`)**
- Referencia probada: **scikit-learn** (`RandomForestRegressor`, `HistGradientBoostingRegressor`) con hiperparámetros idénticos, mismos seeds y mismo split; métricas estándar (R², RMSE, MAE) con validación cruzada.
- Test equivalente: mismo dataset y features → comparar predicción por muestra (en árboles equivalentes puede divergir por desempates; comparar R²/RMSE con tolerancia y predicciones por muestra ≤ 1e-3 si se replican los splits/desempates), importancia de features (orden), y curva de aprendizaje.
- Datasets: ESOL (logS), Lipophilicity (logP), FreeSolv (hidratación), AqSolDB (9.9k solubilidad); features ECFP/descriptores del propio repo.
- Procedimiento: `scripts/compare_sklearn_models.py` (nuevo) → informe; test `tests/regression/test_new_features.rs` (modelos), ampliar.
- Estado en el repo: `src/ml/{models,advanced_models,ensemble}.rs`, tests existentes.
- Criterio: R² dentro de ± 0.02 del de sklearn; RMSE ± 5 %; importancia de features con Spearman ≥ 0.9.

**F61 · Predicción de propiedades de fármaco, solubilidad y reglas empíricas (`predict_ml_properties`, `compute_tpsa`)**
- Referencia probada: **RDKit Crippen logP** (logP exacto), **ESOL** (modelo lineal publicado de Delaney para logS) y **pkCSM/ADMET-AI** (modelos públicos documentados) para contraste de tendencia; Lipinski/Veber/Ghose según papers.
- Test equivalente: comparar logP Crippen con RDKit (exacto); logS contra el modelo ESOL original y contra ESOL experimental (MAE); clasificación druglike contra RDKit `Descriptors`.
- Datasets: ESOL (https://github.com/deepchem/deepchem ESOL), AqSolDB (Mansouri 2020), Lipophilicity (MoleculeNet).
- Procedimiento: `tests/regression/test_extended_molecules.rs` (existentes) + `scripts/compare_logs_esol.py`.
- Estado en el repo: `src/ml/inference.rs`, `src/ml/mod.rs`, tests.
- Criterio: logP ≤ 0.1 (Crippen); logS MAE ≤ 1 log unit en ESOL; clasificación Lipinski ≥ 95 %.

**F62 · pKa empírico (`reactivity::estimate_empirical_pka`, `compute_empirical_pka`)**
- Referencia probada: **SAMPL6/SAMPL7/SAMPL8 pKa challenge** (datasets experimentales ciegos con decenas de compuestos incluyendo bases débiles), y modelos/código abierto **`Chemprop` pKa** o **`pkasolver`** (publicados con métricas); comercial: Epik/ChemAxon/MoKa como referencia no reproducible (documentar solo como contexto).
- Test equivalente: misma molécula → comparar pKa predicho vs experimental en SAMPL (MAE/RMSE); para ácidos y bases por separado; comparar también la asignación del sitio ácido/base más probable (índice de átomo) con el reportado en los papers SAMPL.
- Datasets: SAMPL6/7/8 (https://github.com/samplchallenges), AqSolDB, "DrugBank pKa" (no libre → evitar).
- Procedimiento: `scripts/compare_pka_sampl.py` (nuevo); test `test_pka_cross`.
- Estado en el repo: `src/reactivity.rs` (`estimate_empirical_pka`), API pública existente.
- Criterio: MAE ≤ 1.5 unidades de pKa en SAMPL (modelo empírico); clasificación ácido/base ≥ 90 %.

**F63 · ANI-2x / ANI-TM: energía, fuerzas y AEVs (`ani::compute_ani`, `compute_aevs`, `compute_aevs_tm`, `compute_forces`)**
- Referencia probada: **TorchANI** (implementación de referencia de ANI-2x/ANI-1x, MIT; `torchani.models.ANI2x`) y **ASE** con el calculator ANI; para metales, el paper ANI-TM y, como dato de referencia, **tmQM** (propiedades DFT de 86k complejos TM).
- Test equivalente: mismos pesos y geometría → comparar energía total (eV), energías atómicas y fuerzas (eV/Å); comparar AEVs elemento a elemento (`torchani.aev.compute_aev`) con tolerancia float32; después, errores por átomo vs DFT en el set ANI-1x/2x de test.
- Datasets: **ANI-1x/ANI-1ccx** (figshare, ~5M/500k conformeros), **ANI-2x** (8.9M), **SPICE** (1M+ con fuerzas), **MD17** (10 moléculas, 100k+ geometrías), **tmQM** (metales de transición).
- Procedimiento: `scripts/benchmark_ani.py` (existe) + nuevo `scripts/compare_ani_torchani.py` (mismos pesos exportados) → fixture; tests `tests/alpha_mlff_consistency.rs` (existe).
- Estado en el repo: `src/ani/*` (aev, nn, weights, cutoff, neighbour), `src/ml/mlff.rs`.
- Criterio: ≤ 1e-4 eV/átomo y ≤ 1e-3 eV/Å vs torchani con los mismos pesos; MAE absoluta dentro de lo publicado (± 20 %) en el set de test.

### 4.8 Materiales y sistemas periódicos

**F64 · Celdas, simetría y 230 grupos espaciales (`materials::{create_unit_cell, space_group_by_number/symbol, all_space_groups, apply_symmetry, expand_by_symmetry}`)**
- Referencia probada: **spglib** (detección/simetría, BSD), **pyxtal** (generación y tablas de los 230 grupos) y **pymatgen `SpaceGroup`** (tablas ITC). Para números/símbolos: International Tables for Crystallography.
- Test equivalente: para los 230 números → comparar símbolo Hermann-Mauguin, sistema cristalino, orden del grupo (nº de operaciones) y matrices de las operaciones con las tablas de pymatgen/pyxtal; después, aplicar las operaciones a una estructura de prueba y comparar posiciones equivalentes.
- Datasets: `pymatgen`/`pyxtal` tablas embebidas; CIFs de COD con grupos conocidos (p. ej. NaCl Fm-3m), estructuras representativas de cada sistema (generadas por script).
- Procedimiento: `scripts/compare_space_groups.py` (todas las 230) → fixture comprimido `tests/fixtures/space_groups_230.json.gz`; test `test_space_groups_cross`.
- Estado en el repo: `src/materials/{cell,space_groups}.rs`, `parse_cif`, tests de `test_new_features`.
- Criterio: 230/230 grupos con orden y símbolo correctos; posiciones equivalentes idénticas (± 1e-8 fraccional).

**F65 · Import/export CIF (`materials::parse_cif`, `write_cif`)**
- Referencia probada: **pymatgen `CifParser`/`CifWriter`** y **ASE** (`ase.io.read/write`); **COD** (Crystallography Open Database, libre) y **AMCSD** como corpus.
- Test equivalente: por CIF del corpus → parsear en ambos, comparar red, posiciones fraccionales, ocupaciones, incertidumbres (`5.640(1)`), y recomponer; verificar **round-trip** (parse→write→parse estable) y exportar a formato extendido (XYZ/SDF) correctamente.
- Datasets: COD (https://www.crystallography.net/cod/, ~500k CIF), AMCSD (http://rruff.geo.arizona.edu/AMS/amcsd.php); los CIFs con notación de incertidumbre del test existente.
- Procedimiento: `tests/regression/test_session_algorithms.rs` §2 (existe: NaCl, roundtrip, incertidumbre, campos requeridos); añadir `scripts/compare_cif_pymatgen.py` sobre 200 CIF de COD.
- Estado en el repo: `src/materials/cif.rs`(aprox.), tests §2.
- Criterio: 100 % parseo sin error en CIF válidos; parámetros ≤ 1e-6; round-trip estable; 0 pérdidas de campos requeridos.

**F66 · Molécula periódica, PBC y superceldas (`periodic::build_periodic_molecule`, `materials::{make_supercell, supercell, wrap_*}`, `minimum_image_distance`)**
- Referencia probada: **ASE** (`Atoms`, `cell`, `pbc`, `minimum_image`) y **pymatgen** (`Structure`, `PBC`); casos analíticos con celda cúbica.
- Test equivalente: mínima imagen vs cálculo directo (L0 con vectores cortos); wrap de coordenadas fraccionales en [0,1); superceldas: número de átomos × N³ y continuidad de distancias; comparar `get_distances`/`find_neighbors` vs ASE.
- Datasets: CIFs de COD (NaCl, SiC, grafito), fixtures de `tests/regression/test_reaction_dft_mesh.rs`.
- Procedimiento: `test_new_features.rs`/`test_roadmap_validation.rs` (parcial) + `scripts/compare_pbc_ase.py`.
- Estado en el repo: `src/periodic.rs`, `src/materials/cell.rs`.
- Criterio: mínimo imagen = cálculo por traslaciones ±1 celdas (± 1e-10); supercelda 8× como mínimo reproducible; wrap idempotente.

**F67 · Ensamblado de frameworks y topologías (`materials::assemble_framework`, `pcu/bcu/fcu/dia/hxl/kgm/nbo/pts/sql/pillared_sql`, `place_sbu_cga`)**
- Referencia probada: **pormake** (constructor de MOFs open-source, MIT; SBUs + topologías RCSR), **RCSR** (base de datos de topologías) y bases **CoRE MOF 2019 / ARC-MOF / hMOF** para estructuras ya resueltas.
- Test equivalente: mismo SBU + misma topología → comparar el grafo de conectividad (nodos/grado), parámetros de celda esperados y coherencia geométrica; validar porosidad y comparar con la estructura CoRE del mismo MOF (MOF-5/UiO-66/ZIF-8 como casos canónicos).
- Datasets: CoRE MOF 2019 (~14k estructuras con propiedades geométricas calculadas), ARC-MOF, RCSR topologías; `tests/regression/test_%` de ensamblado.
- Procedimiento: `scripts/compare_pormake.py` (mismos SBU/topología) → fixture; tests `test_new_features.rs`/`test_roadmap_validation.rs`.
- Estado en el repo: `src/materials/{assembly,sbu,topologies}.rs`, `assemble_framework`, `assemble_framework_cga`.
- Criterio: celda y conectividad correctas para MOF-5/UiO-66/ZIF-8 (± 1 % en parámetros vs CoRE); diámetro de poro dentro de ± 10 % de Zeo++ (ver F68).

**F68 · Porosidad geométrica (`materials::compute_porosity`)**
- Referencia probada: **Zeo++** (implementación estándar: diámetro de poro, ASA, fracción vacía; versión académica gratuita) y **pymatgen** (`Porosity`/BET-like no; usar `CrystalNN`/`Zeo++` vía interfaz).
- Test equivalente: misma estructura (CIF) → comparar diámetro de poro, superficie accesible y fracción vacía con Zeo++ (radios y sondas idénticos); casos analíticos: celda vacía (fracción = 1) y esfera de radio conocido.
- Datasets: CoRE MOF 2019 con sus propiedades geométricas publicadas (superficie/diámetro) como tabla de contraste; ARC-MOF.
- Procedimiento: `scripts/compare_porosity_zeopp.py` (invoca Zeo++ con la red de prueba) → fixture; test `test_porosity_cross`.
- Estado en el repo: `src/materials/porosity.rs`(aprox.), `docs/algorithm/materials-assembly.md`.
- Criterio: diámetro y ASA ≤ 5 % vs Zeo++ (mismos radios/sondas); casos analíticos ≤ 1e-3.

**F69 · XRD de polvo simulado (`materials::simulate_powder_xrd`)**
- Referencia probada: **pymatgen `XRDCalculator`** (BSD) y **GSAS-II** (referencia de cristalografía, opcional); patrones experimentales de COD/comunidad.
- Test equivalente: misma estructura y longitud de onda (Cu Kα) → comparar posiciones de pico (2θ) ≤ 0.01° e intensidades normalizadas (correlación ≥ 0.95); caso canónico NaCl y Si.
- Datasets: CIFs COD de NaCl/Si/Cu; patrones de referencia del PDF/ICDD (licencia) → usar solo CIF + pymatgen como ground truth.
- Procedimiento: `scripts/compare_xrd_pymatgen.py` → fixture; test `test_xrd_cross`.
- Estado en el repo: `src/materials/xrd.rs`(aprox.), `simulate_powder_xrd` público.
- Criterio: 2θ ≤ 0.01°; correlación de intensidades ≥ 0.95; lista de picos sin reflexiones extrañas.

**F70 · Optimización de frameworks (`materials::optimize_framework`, `optimize_frame_sequence`)**
- Referencia probada: **ASE/pymatgen + LBFGS** con el mismo término UFF ya validado en F16 (usar la energía como contrato), y **LAMMPS** con UFF (`pair_style` no nativo; usar `reaxff` no aplica) → la referencia práctica es ASE `BFGS` sobre la misma función con PBC, más geometrías experimentales/CoRE como destino.
- Test equivalente: misma estructura inicial → convergencia (fuerza máxima ≤ 0.01 eV/Å), energía final ≤ E_min de referencia + 1e-6, RMSD de celda/posiciones, y desplazamiento de átomos (las posiciones fraccionales deben permanecer físicamente razonables).
- Datasets: CoRE MOF (estructuras experimentales), `tests/fixtures/metal_experimental_geometry.json` no aplica; gaps de banda opcional.
- Procedimiento: `test_new_features.rs` (parcial) + `scripts/compare_framework_opt_ase.py`.
- Estado en el repo: `src/materials/optimize.rs`(aprox.), `optimize_framework`, `src/forcefield/bounds_ff/*`.
- Criterio: E_final ≤ E_ref + 1e-6 eV/átomo; fuerza máx ≤ 0.01 eV/Å; sin cruces de átomos (distancia mínima ≥ 0.7 Å).

### 4.9 Dinámica, transporte y utilidades

**F71 · MD básica: Velocity-Verlet, NVE/NVT (`dynamics::simulate_velocity_verlet*`, `compute_md_trajectory*`)**
- Referencia probada: **OpenMM** (MIT, reference platform determinista) y **LAMMPS** (GPL-2) para integradores/termostatos; **oscilador armónico** como caso analítico exacto.
- Test equivalente: (a) oscilador armónico 1D con paso pequeño → error de energía ≤ 1e-4 tras 10⁴ pasos; (b) NVE: deriva de energía total por paso ≤ 1e-5 kcal/mol; (c) NVT Berendsen/Nose-Hoover: temperatura media estable y fluctuaciones según teoría (~1/√N); (d) mismas fuerzas que el FF validado en F16.
- Datasets: sin descarga; snapshots propios; opcional `MD17` (trayectorias DFT de 10 moléculas) para comparar estadística de fuerzas.
- Procedimiento: `tests/alpha_dynamics_conservation.rs` (existe) + `test_dynamics_cross` con OpenMM en un sistema sencillo (agua UFF-like o LJ).
- Estado en el repo: `src/dynamics.rs`, `src/dynamics_live/*`, tests de conservación.
- Criterio: conservación NVE ≤ 1e-5 por paso; T ± 2 % del setpoint; fuerzas = FF validado.

**F72 · Dinámica interactiva y barostats (`dynamics_live::*`, `berendsen_thermostat`, `nose_hoover_thermostat`, `langevin_*`)**
- Referencia probada: **LAMMPS** `fix nvt/npt` (Berendsen/Nose-Hoover/Langevin) y **OpenMM** `LangevinMiddleIntegrator`; fórmulas de referencia en Frenkel & Smit.
- Test equivalente: comparar la evolución de temperatura/presión para un sistema LJ con los termostatos homónimos de LAMMPS; L0: presión virial vs cálculo directo (sistema analítico de pares), volumen fluctuante del baratostato con coeficientes correctos.
- Datasets: sistema LJ estándar (Ar líquido, parámetros publicados), agua SPC/E opcional; sin descarga obligatoria.
- Procedimiento: `tests/alpha_dynamics_conservation.rs` + `src/dynamics_live/{langevin,barostat,pbc}.rs` (tests internos) + comparación LAMMPS.
- Estado en el repo: `src/dynamics_live/*` (alpha), API WASM `LiveSimulation`.
- Criterio: T media ± 1 %; presión ± 5 %; estabilidad de volumen ± 2 %.

**F73 · IMD (puente de protocolo, `alpha-imd`)**
- Referencia probada: **VMD** (cliente IMD) y **NAMD** (protocolo IMD original de Klein/Phillips); la especificación del protocolo IMD.
- Test equivalente: conectar `sci-form` a VMD como cliente → verificar handshake, tasa de frames, fuerzas de steering aplicadas y desconexión limpia; test de conformidad: bytes del mensaje comparados con la especificación (L0).
- Datasets: sin descarga; VMD para prueba manual/integrada.
- Procedimiento: `examples/` react viewer; `tests/regression/test_new_pipeline.rs` (relacionado).
- Estado en el repo: `src/alpha/imd.rs`(aprox.), `dynamics_live/steering.rs`; WASM `apply_steering_force`, `set_temperature`.
- Criterio: VMD muestra ≥ 30 fps con trayectoria coherente; suma de fuerzas de steering coherente; cierre limpio.

**F74 · Transporte por lotes: Arrow, chunking y workers (`transport::*`)**
- Referencia probada: **Apache Arrow/pyarrow** (formato canónico; round-trip bit-exacto), y para robustez, `pandas`/`polars` como consumidores.
- Test equivalente: empaquetar un lote → escribir IPC → leer con pyarrow → comparar esquema, nombres de columnas, tipos y valores exactos; determinismo de `split_worker_tasks` (misma entrada → misma partición) y `merge_results` idempotente.
- Datasets: `tests/fixtures/benchmark_100_molecules.json` (lote de prueba), `data/chembl_1k_practical.smi`.
- Procedimiento: `tests/integration/test_cross_language.py` (existe) + `scripts/compare_arrow_pyarrow.py`; tests de `src/transport/{arrow,chunked,worker}.rs`.
- Estado en el repo: `src/transport/*`, WASM `pack_batch_arrow`, `pack_conformers`.
- Criterio: round-trip bit-exacto; partición determinista; cobertura de tests de chunks = suma de resultados individuales.

**F75 · Potenciales analíticos y modelos 1D/2D (`potentials.rs`, `reactivity::{sn2_model_1d, proton_transfer_2d_morse, proton_transfer_evb}`)**
- Referencia probada: **Muller-Brown** (superficie publicada con mínimos/sillines conocidos), **Morse** y **double-well** con solución analítica; modelos EVB de SN2 y transferencia de protón de los papers de referencia (Truhlar; Gao; Warshel), con energías de barrera publicadas.
- Test equivalente: evaluar en malla → reproducir mínimos/sillines publicados (± 1e-6 en los analíticos, ± 1e-3 en Muller-Brown); gradientes vs diferencias finitas; para EVB: reproducir ΔE‡ publicado en el sistema modelo.
- Datasets: tablas de los papers (Muller-Brown 1979; EVB SN2) y soluciones analíticas; opcional `Transition1x` para comparar barreras DFT de reacciones reales.
- Procedimiento: `tests/e3_justification.rs`, `tests/reaction_dynamics_study.rs` (existen) + tests unitarios de `potentials.rs`.
- Estado en el repo: `src/potentials.rs` (módulo de utilidades para tests).
- Criterio: analíticos ≤ 1e-8; Muller-Brown ≤ 1e-3 en energías de referencia; EVB ΔE‡ ± 0.5 kcal/mol del paper.

### 4.10 Rutas GPU

**F76 · Paridad GPU↔CPU y política de fallback (`gpu::*`, feature `experimental-gpu`)**
- Referencia probada: la propia CPU es la referencia (L1): MMFF94 nonbonded GPU, AEV GPU, ERIs, Fock, densidad/rejilla orbital, Marching Cubes GPU, ESP/EEG/CPM/D4/Hessian GPU.
- Test equivalente: mismo input → comparar resultados GPU vs CPU (misma tolerancia que el observable: 1e-6 para fuerzas, 1e-4 para energías); verificar la **política de fallback** (`webgpu_regression_policy.json`, `get_report`/`activation_report`: cpu_only → cpu_fallback → fully_gpu) y que un fallo de GPU no cambie resultados ni silenciosamente el método.
- Datasets: `tests/fixtures/webgpu_regression_policy.json` (existente), casos de `tests/regression/test_gpu_candidates.rs` (existente), lotes de GDB-20.
- Procedimiento: `node tests/integration/test_webgpu_harness.js` y `compare_webgpu_reports.js` (existen, en `package.json`), `cargo test --features experimental-gpu`; tests `test_gpu_candidates.rs`.
- Estado en el repo: `src/gpu/*` (~40 archivos), harness JS existente, política de regresión versionada.
- Criterio: paridad ≤ 1e-6 relativa con CPU en todos los kernels cubiertos; los reportes GPU/CPU coinciden con el esquema declarado; ningún path GPU altera el resultado si cae a CPU.

### 4.11 Módulos alpha

> Los módulos alpha se compilan con `--features alpha-*`. Los criterios de aceptación son **más
> descriptivos que bloqueantes** hasta su promoción a beta; la excepción son los invariantes L0,
> que siempre son obligatorios.

**A1 (F77) · KS-DFT SVWN/PBE (`dft::solve_ks_dft`, `functionals::{svwn,pbe}`, `build_vxc_matrix`)**
- Referencia probada: **PySCF** (`dft.RKS` con `xc='lda,vwn'` y `xc='pbe'`, misma base/grid; libxc como backend canónico) y **ORCA** como segundo voto. Para precisión de funcional: **AE6** (atomizaciones), **BH6** (barreras), **S22** (no covalente) publicados por los papers DFT de validación.
- Test equivalente: misma geometría → comparar S/T/V/J/K, matriz Vxc en la misma rejilla, energía XC, energía total, ε_HOMO/LUMO y Mulliken; test de H2 con rejilla fina ya en el repo (último commit de la rama).
- Datasets: H2/HeH⁺/H2O/CH4/NH3/HF/CO/C2H4/CH2O (ya en `test_experimental_comparison.rs`), CCCBDB (energías LDA/PBE), AE6/BH6/S22 (tablas públicas).
- Procedimiento: `tests/alpha_dft_validation.rs` (existe) + `test_alpha_cross_method.rs`; nuevo `scripts/generate_dft_reference.py` (PySCF con la misma base/grid).
- Estado en el repo: `src/dft/*` (ecp, grid, lebedev, becke, functionals, ks_fock), docs `ks-dft.md`, tests existentes.
- Criterio: E_total ≤ 1e-5 Ha con la misma rejilla; Vxc en la misma rejilla ≤ 1e-6; ε ≤ 1e-4 Ha; AE6/BH6 dentro del error esperado del funcional (no del código).

**A2 (F78) · ReaxFF (`forcefield::reaxff::*`, `compute_reaxff_gradient`, EEM)**
- Referencia probada: **LAMMPS** `pair_style reaxff` con el **mismo ffield** (referencia canónica del campo, implementación de van Duin/Acton), **PuReMD** como segundo voto; las tablas de cargas EEM se comparan contra el `fix qeq/reaxff` de LAMMPS.
- Test equivalente: misma geometría, mismos parámetros (`ffield.reax`), desapareo de cargas (`qeq`) y pasos de energía: comparar E_bonded, E_vdW, E_coulomb, E_total y gradientes centrales vs LAMMPS (un solo punto, sin dinámica); verificar cargas EEM y neutralidad.
- Datasets: ejemplos de LAMMPS (`examples/reax`), parámetros públicos de C/H/O/N (Adri van Duin), moléculas CHON + RDX/trinitro para reactividad; `tests/fixtures/experimental_reference.json` para geometrías.
- Procedimiento: `tests/alpha_reaxff_validation.rs` (existe); nuevo `scripts/compare_reaxff_lammps.py` (genera dump LAMMPS para una geometría y extrae `pe`/`f`) → fixture.
- Estado en el repo: `src/forcefield/reaxff/*` (bond_order, eem, energy, gradients, taper), GPU/parallel para gradientes, bindings Python/WASM.
- Criterio: E_total ≤ 1e-3 kcal/mol; componente ≤ 1e-3 kcal/mol; cargas EEM ≤ 0.01 e; gradiente ≤ 1e-3 kcal/mol/Å (paridad CPU/GPU incluida).

**A3 (F79) · MLFF con AEV (Behler-Parrinello) (`ml::mlff`, `mlff::compute_mlff*`, `symmetry_functions`)**
- Referencia probada: **TorchANI** (AEVs `torchani.aev.compute_aev` y red de referencia), **n2p2** (Ruhr-University, implementación BP de referencia) y las fórmulas originales de Behler-Parrinello (2007/2011); datasets **MD17**, **SPICE** y **ANI-1x** para comparar energía/fuerzas por átomo.
- Test equivalente: mismos pesos/parámetros → AEV radial y angular elemento a elemento (float64 vs float32 tolerancia); energía/átomo y fuerzas vs TorchANI o n2p2; paridad serial/paralelo (rayon) ya en el repo.
- Datasets: MD17 (10 moléculas, ≤100k frames con fuerzas), SPICE (1M+), ANI-1x (figshare); puntos de test sugeridos por el propio dataset.
- Procedimiento: `tests/alpha_mlff_consistency.rs` (existe), `scripts/` (bench de ANI); nuevo `scripts/compare_aev_n2p2.py`.
- Estado en el repo: `src/ml/{symmetry_functions,inference,mlff}.rs`, `src/ani/aev*`, feature `alpha-mlff`.
- Criterio: AEV ≤ 1e-6 rel; energía ≤ 1e-3 eV/átomo vs referencia con mismos pesos; fuerzas ≤ 1e-3 eV/Å; parallel == serial ≤ 1e-12.

**A4 (F80) · Obara-Saika expuesto como alpha (`scf::obara_saika::{boys_function, eri_ssss, schwarz_bound}`, bindings `alpha_*`)**
- Referencia probada: **PySCF/libcint** (ERIs) y `scipy.special`/series para la función de Boys; fórmulas cerradas para `F0`–`F4` (Boys 1950; Papa-King tamizado).
- Test equivalente: `boys_function(m, x)` vs referencia tabulada para x∈[0,50], m∈[0,10] (rel 1e-12); `eri_ssss` vs expresión cerrada de s-funciones y vs PySCF con exponentes de prueba; `schwarz_bound` verifica la cota.
- Datasets: sin descarga; exponentes/centros de la suite de PySCF (`pyscf/gto/test/test_ints.py`).
- Procedimiento: `tests/alpha_scf_eri_accuracy.rs` (existe) + `scripts/compare_eri_pyscf.py`.
- Estado en el repo: `src/scf/obara_saika.rs`, bindings `alpha_boys_function`, `alpha_eri_ssss`, `alpha_schwarz_bound`.
- Criterio: Boys ≤ 1e-12 rel; ERI ≤ 1e-10 rel; cota de Schwarz nunca violada.

**A5 (F81) · CGA: rotores conformes para dihedros (`alpha::cga::{dihedral_motor, apply_motor_to_subtree, refine_torsion_cga}`)**
- Referencia probada: **clifford** (librería Python de álgebra geométrica) y **galgebra** como independientes; la referencia matemática es la rotación por cuaternión/matriz (exactitud L0).
- Test equivalente: aplicar el rotor a un subárbol y comparar con la rotación cuaterniónica equivalente (≤ 1e-12); invariancia de distancias internas; cierre del motor (`M M† = 1`); y refinamiento de torsión: mismo dihedro final que el método de rotación matriz (`refine_torsion_cga` vs `rotate_atoms`).
- Datasets: geometrías GDB-20 (torsiones extraídas), `tests/fixtures/torsion_params.json`.
- Procedimiento: `tests/experimental/test_cga.rs` (existe) + `scripts/check_cga_vs_quaternion.py`.
- Estado en el repo: `src/alpha/cga/*` (motor, multivector, conformer, materials), tests existentes.
- Criterio: ≤ 1e-12 en coords vs cuaternión (f64); preservación exacta de longitudes de enlace ≤ 1e-12.

**A6 (F82) · Growing String Method (`alpha::gsm::{find_transition_state*, grow_string, interpolate_node, refine_saddle}`)**
- Referencia probada: **pyGSM** (Zimmerman group, implementación de referencia del GSM, GPL) y **pysisyphus** (GSM/NEB/IRC, MIT); como verdad física, **Transition1x** (TS DFT de 10k reacciones) y barreras publicadas.
- Test equivalente: misma pareja reactivo/producto y mismo backend (UFF/PM3/xTB) → comparar nodos de la cuerda (RMSD ≤ 0.05 Å por nodo), energía de la barrera, geometría del TS (RMSD) y verificar que el TS es un punto estacionario (‖g‖ ≤ 0.01 eV/Å); comparar con barreras DFT de Transition1x.
- Datasets: Transition1x (figshare/Zenodo), set de reacciones curado del repo (`tests/reaction_battery.rs`), ejemplos pyGSM (SN2, H-transfer).
- Procedimiento: `tests/experimental/test_gsm.rs` (existe), `scripts/compare_reaction_layers.py` (existe), CLI `gsm-*`; nuevo `scripts/compare_gsm_pygsm.py`.
- Estado en el repo: `src/alpha/gsm/{string,saddle,backend}.rs`, `alpha-kinetics` (GSM+MBH+HTST).
- Criterio: barrera ± 1 kcal/mol vs pyGSM con el mismo backend; ‖g‖_TS ≤ 0.01 eV/Å; camino sin discontinuidades (distancia entre nodos monótona).

**A7 (F83) · SDR: embedding por relajación semidefinida (`alpha::sdr::{sdr_embed, project_psd, alternating_projections, svt_step}`)**
- Referencia probada: resolución independiente del mismo problema con **CVXPY** (SCS/MOSEK/Clarabel) sobre las mismas restricciones de distancia, y el algoritmo **SFSDP** (papers de Kim et al./Biswas-Ye) como referencia del enfoque.
- Test equivalente: mismas tripletas (i,j,d) → comparar matriz de Gram óptima y coordenadas finales (RMSD ≤ 1e-3 Å tras Kabsch); violación de restricciones ≤ 1e-3 Å; comportamiento con restricciones inconsistentes (infactible → residual acotado).
- Datasets: grafos de distancia de moléculas GDB-20 (bounds → distancias), ejemplos de SFSDP en los papers.
- Procedimiento: `tests/experimental/test_sdr.rs` (existe) + `scripts/compare_sdr_cvxpy.py`.
- Estado en el repo: `src/alpha/sdr/{embedding,projections}.rs`, bindings Python/WASM.
- Criterio: violación máx ≤ 1e-3 Å; RMSD vs CVXPY ≤ 0.01 Å; coste O(N·k) documentado.

**A8 (F84) · Dynamics-live (`dynamics_live::*`, feature `alpha-dynamics-live`)**
- Referencia probada: **OpenMM** (plataforma Reference determinista) y **LAMMPS**; caso analítico del oscilador armónico; ver F71/F72.
- Test equivalente: paridad de un paso con el integrador clásico (L0); callbacks por frame; termostatado Berendsen; PBC (`wrap`, `minimum_image`) vs ASE.
- Datasets: snapshots propios; Ar líquido LJ publicado para barostato.
- Procedimiento: `tests/alpha_dynamics_conservation.rs` (existe), WASM `LiveSimulation` (examples JS).
- Estado en el repo: `src/dynamics_live/*`, tests de conservación.
- Criterio: deriva de energía ≤ 1e-5/step; T ± 2 %; nº de frames por callback exacto.

**A9 (F85) · IMD (`alpha-imd`)**
- Referencia probada: **VMD** (cliente) y la especificación IMD de **NAMD**; ver F73.
- Test equivalente: conformidad de mensajes (bytes) + sesión real con VMD; steering aplicado; reconexión.
- Datasets: sin descarga.
- Procedimiento: `examples/` y prueba manual con VMD; tests de `steering.rs`.
- Estado en el repo: protocolo en `src/alpha/`/`dynamics_live`, WASM `apply_steering_force`.
- Criterio: handshake correcto; ≥ 30 fps; sin cuelgues al desconectar.

**A10 (F86) · EDL: modelos de doble capa (`alpha::edl::*`: Helmholtz, Gouy-Chapman, GCS, perfiles, capacitancia; `cpm ↔ eeq ↔ edl`)**
- Referencia probada: soluciones analíticas de Gouy-Chapman/Helmholtz y **Grahame**; contraste numérico con `scipy.integrate`/`sympy`; para curvas de capacitancia, datos experimentales publicados de electrodos (p. ej. Hg/Pt en NaF, tablas clásicas).
- Test equivalente: evaluar el perfil de potencial/carga en malla → comparar con la fórmula cerrada (≤ 1e-6 rel); límites: concentración → 0 (Debye), φ pequeño (linealización), φ grande (límite de Stern); derivada dσ/dφ = capacitancia comparada con derivada numérica.
- Datasets: curvas capacitancia-potencial experimentales de Grahame (Hg) y datos de superficies metálicas en literatura electroquímica; `tests/fixtures` propios.
- Procedimiento: `crates/cli/src/experimental_cmds.rs` (`edl-*`), tests `tests/experimental/test_cpm.rs` relacionados; nuevo `scripts/check_edl_analytics.py`.
- Estado en el repo: `src/alpha/edl/*` (incluye GPU), docs `experimental/*` de EDL; CLI `edl-profile`, `edl-gcs`, `edl-helmholtz`, `capacitance-scan`.
- Criterio: fórmulas cerradas ≤ 1e-6; capacitancia = derivada numérica ± 1e-3; monotonías físicas correctas.

**A11 (F87) · Cinética: HTST, redes microcinéticas y Arrhenius (`alpha::kinetics::*`)**
- Referencia probada: **Eyring/TST** analítico (k = (k_BT/h)·exp(−ΔG‡/RT)); **Cantera** (integración de redes ODEs, referencia estándar) y **CatMAP** (microcinética/TST de superficies, referencia académica) con los mismos parámetros.
- Test equivalente: mismas energías libres → comparar k(T) con Eyring (≤ 0.1 %); integrar red microcinética (estado estacionario y transitorio) y comparar concentraciones/TOF con Cantera/CatMAP (± 1 %); barrido de temperatura vs Arrhenius (pendiente = −Ea/R).
- Datasets: barreras DFT de Transition1x (reacciones aisladas); energías de adsorción publicadas (CatMAP examples); datos de TOF de literatura.
- Procedimiento: `tests/experimental/test_gsm.rs`/`test_mbh.rs` (complementarios), CLI `htst-rate`, `solve-microkinetic-*`, `arrhenius-chart`; nuevo `scripts/compare_microkinetics_cantera.py`.
- Estado en el repo: `src/alpha/kinetics/*`, `src/alpha/render_bridge` (charts), CLI asociada.
- Criterio: k(T) ≤ 0.1 % vs Eyring (mismos ΔG‡); red ≤ 1 % vs Cantera; conservación de masa en la red ≤ 1e-10.

**A12 (F88) · Periódico lineal: Bloch, k-mesh, KPM periódico, RandNLA periódico (`alpha::periodic_linear::*`)**
- Referencia probada: **pybinding** y **pymatgen** (ver F34) + modelos TB analíticos; **scipy** para autovalores exactos en cada punto k.
- Test equivalente: cadena 1D/2D → dispersión analítica; hermeticidad del Hamiltoniano Bloch (H(k) = H†(k)); periodicidad H(k+G) = H(k); pesos `monkhorst_pack_mesh` vs pymatgen; KPM periódico vs diagonalización exacta (L2 ≤ 5 %).
- Datasets: modelos TB publicados (SSH, grafeno, H-chain); CIFs simples de COD.
- Procedimiento: `tests/experimental/test_kpm.rs`, CLI `build-bloch-hamiltonian`, `solve-periodic-randnla`, `compute-periodic-dos`.
- Estado en el repo: `src/alpha/periodic_linear/*` (incluye GPU), bindings Python/WASM.
- Criterio: analítico ≤ 1e-6; H(k) hermítica ≤ 1e-10; DOS relativo ≤ 5 % vs exacto.

**A13 (F89) · Dinámica de reacción 3D: complejo, sitios reactivos, approach, IRC, CI-NEB (`alpha::reaction_dynamics::*`)**
- Referencia probada: **ASE NEB** (`ase.mep.NEB`, con tblite para fuerzas), **geomeTRIC ElasticBand** (ya integrado en los scripts del repo), **pysisyphus** (IRC/NEB), y **Transition1x** como verdad DFT.
- Test equivalente: (a) mismo par reactivo/producto → comparar caminos NEB (RMSD por imagen ≤ 0.05 Å) y barreras (≤ 1 kcal/mol); (b) IRC: el camino debe comenzar en el TS y terminar en los mínimos (energías ≤ E_TS, gradiente ⊥ cero); (c) `identify_reactive_sites` vs el mapeo SMIRKS (test existente); (d) orientación del approach vs orbitales de frontera (solo dirección, invariante).
- Datasets: Transition1x (figshare), reacciones curadas de `tests/reaction_battery.rs`, ejemplos de ASE/pysisyphus.
- Procedimiento: `scripts/compare_reaction_layers.py --external` (existe, comparación ASE/tblite/geomeTRIC), tests `alpha_reaction_dynamics_3d.rs` (existe), `tests/reaction_dynamics_study.rs`.
- Estado en el repo: `src/alpha/reaction_dynamics/*` (11 archivos), CLI `neb-*`, `simplified-neb-path`.
- Criterio: barrera ± 1 kcal/mol vs ASE NEB con el mismo motor; camino sin discontinuidades; IRC asimétrico correcto; sitios reactivos = mapeo SMIRKS en 100 %.

**A14 (F90) · Puente de render/charts (`alpha::render_bridge::*`, CLI `*-chart`, `validate_chart_schema`)**
- Referencia probada: **JSON Schema** (draft 2020-12) para validar payloads y **Vega-Lite/Plotly** como consumidores de referencia; los charts del repo se validan contra su propio esquema versionado.
- Test equivalente: round-trip `chart_to_json → chart_from_json` bit-exacto; validación contra el esquema (`ajv` o `jsonschema`); comparación de payloads GPU vs CPU; verificación de los campos temporales/etiquetas.
- Datasets: `tests/regression/test_gpu_candidates.rs` y `tests/integration/compare_webgpu_reports.js` para los charts GPU; sin descargas externas.
- Procedimiento: `cargo test --test regression test_gpu_candidates`; `node tests/integration/compare_webgpu_reports.js`; CLI `validate-chart-schema`.
- Estado en el repo: `src/alpha/render_bridge/*`, `crates/wasm/js`.
- Criterio: round-trip exacto; 0 violaciones de esquema; GPU/CPU idénticos.

### 4.12 Módulos beta

**B1 (F91) · KPM DOS (`beta::kpm::{ChebyshevExpansion, compute_kpm_dos, compute_kpm_mulliken}`)**
- Referencia probada: **pybinding `kpm`** (implementación KPM estándar, BSD) y la **diagonalización exacta** (numpy `eigh`/scipy) para el mismo Hamiltoniano; límites: orden alto → convergencia al DOS exacto.
- Test equivalente: mismo H y rejilla [e_min,e_max] → comparar DOS a igual orden (L2 rel ≤ 5 %) y la convergencia con el orden (500 ≤ 2 %); kernel de Jackson y factor de escala verificados contra fórmulas (L0); `compute_kpm_mulliken` vs Mulliken exacto de EHT (≤ 0.01 e); promedios estocásticos vs semilla fija.
- Datasets: modelos TB analíticos y Hamiltonianos EHT de GDB-20; `tests/fixtures/` locales; pybinding examples.
- Procedimiento: `tests/experimental/test_kpm.rs` (existe), CLI `compute-periodic-dos`/`dos-chart`; nuevo `scripts/compare_kpm_pybinding.py`.
- Estado en el repo: `src/beta/kpm/*` (parallelizado), bindings Python/WASM.
- Criterio: L2 ≤ 5 % (orden 100–500) vs exacto; límites espectrales ± 2 %; Mulliken ≤ 0.01 e.

**B2 (F92) · Mobile Block Hessian (`beta::mbh::{compute_mbh_frequencies, detect_rigid_blocks, build_projection_matrix}`)**
- Referencia probada: el **Hessiano completo** del mismo método (referencia L0/L1 para los modos flexibles), y la literatura MBH (Ghysels, Van Neck, Van Speybroeck et al. 2007–2010) con valores de sistemas de prueba (agua, metanol, proteínas pequeñas si aplica).
- Test equivalente: para moléculas donde el Hessiano completo es viable → comparar la proyección (fracción de modos flexibles recuperados), frecuencias de los modos flexibles (± 2 cm⁻¹), y el factor de aceleración declarado (nº de DOF reducidos). Para sistemas grandes: verificación de que la proyección preserva invariancias (traslación/rotación) y de que el espectro proyectado ⊆ espectro completo (entrelazamiento de Cauchy).
- Datasets: H2O/metanol/B3LYP de CCCBDB como contraste de frecuencias; clusters de agua (WATER27) para escalabilidad.
- Procedimiento: `tests/experimental/test_mbh.rs` (existe), `test_session_algorithms`/`alpha-kinetics` (GSM+MBH+HTST).
- Estado en el repo: `src/beta/mbh/*` (bloques, hessian), tests existentes.
- Criterio: modos flexibles ± 2 cm⁻¹ vs Hessiano completo; entrelazamiento correcto; speedup documentado y reproducible.

**B3 (F93) · RandNLA: EHT aleatorizado (`beta::rand_nla::{solve_eht_randnla, Nystrom, inverse_sqrt, from_matrix}`)**
- Referencia probada: **scipy.sparse.linalg.eigsh** (Lanczos, referencia de autovalores de extremo) y **numpy.linalg.eigh** (exacto) para N pequeño; el marco de **Nyström aleatorizado** (Halko-Martinsson-Tropp 2011) y su análisis de error.
- Test equivalente: mismo H,S → comparar los k autovalores más bajos (rel ≤ 1e-3), residuo de la Nyström (`frobenius_error`/`relative_error` declarados vs medidos), y la bandera de fallback (debe activarse cuando el residuo > tolerancia); repetir con varias semillas (error estable).
- Datasets: Hamiltonianos EHT de GDB-20 (1k), modelos TB; sin descarga externa.
- Procedimiento: `tests/experimental/test_randnla.rs` (existe), CLI `solve-periodic-randnla`, bindings `eht_randnla`.
- Estado en el repo: `src/beta/rand_nla/*`, tests existentes; API con `max_error`.
- Criterio: autovalores ≤ 1e-3 rel; error reportado ≥ error real (no mentir en la cota); fallback correcto y determinista con semilla.

**B4 (F94) · Geometría Riemanniana del cono PSD (`beta::riemannian::{psd_distance, psd_projection, riemannian_gradient, psd_retraction}`)**
- Referencia probada: **geomstats** (`SPDMetricAffine`: `dist`, `exp`, `log`) y **pymanopt** (`SPDAffineMetric`); fórmula de Bhatia (afine-invariant); `scipy.linalg.{logm,expm,eigh}` para verificaciones.
- Test equivalente: distancia afín-invariante entre pares A,B vs geomstats (≤ 1e-9); invariancia bajo congruencia A→GAGᵀ (L0); proyección = clamp de autovalores negativos a 0 (comparar con `eigh`); gradiente riemanniano = G − X G X en el punto X (comparar con scipy); retracción/geodésica cumple `dist(R(X,η), X) = ‖η‖_X`.
- Datasets: matrices de ejemplo de la clase (covarianzas, `from_psd`), problemas de optimización del repo; sin descarga.
- Procedimiento: `tests/experimental/test_riemannian.rs` (existe) + `scripts/compare_riemannian_geomstats.py`.
- Estado en el repo: `src/beta/riemannian/*` (manifold, lbfgs, objetivo), tests existentes.
- Criterio: ≤ 1e-9 vs geomstats; invariancias exactas ≤ 1e-12; optimización L-BFGS converge en los problemas de prueba.

**B5 (F95) · Constant Potential Method (`beta::cpm::{compute_cpm_charges, compute_cpm_surface, grand_potential}`)**
- Referencia probada: **LAMMPS** `fix electrode/conp` (CPM constante, implementación publicada por los autores del método) sobre la misma topología/electrodo; **VASPsol/GC-DFT** como referencia de gran potencial (opcional); el formalismo de Siepmann/Reed (física de la capacitancia diferencial).
- Test equivalente: (a) capacitor analítico de placas: Q(φ) lineal con pendiente = C; (b) L0: neutralidad global (Q_sistema + Q_electrodo = 0), simetría Q(−φ) = −Q(φ) en geometría simétrica; (c) comparar magnitudes de carga inducida y curva C(φ) con LAMMPS `fix electrode/conp` para un electrolito modelo; (d) consistencia con EEQ/EDL (`couple_cpm_eeq_edl`).
- Datasets: electrodo Au(111)/grafeno en agua (configuraciones de ejemplo de LAMMPS/electrode package), `tests/experimental/test_cpm.rs`; curvas experimentales de capacitancia de Grahame (contraste de orden).
- Procedimiento: `tests/experimental/test_cpm.rs` (existe), CLI `cpm`, `cpm-surface`, `couple_cpm_eeq_edl`; nuevo `scripts/compare_cpm_lammps.py`.
- Estado en el repo: `src/beta/cpm/*`, `src/gpu/cpm_gpu.rs`, bindings Python/WASM.
- Criterio: capacitor analítico ± 1e-6; simetría exacta; vs LAMMPS ± 5 % en carga/potencial; conservación de carga ≤ 1e-8 e.

## 5. Crate `sci-form-cli`

**Rol en la validación:** el CLI es la superficie que permite comparar *sin escribir código* contra
binarios externos (MOPAC, xtb, tblite, dftd4, LAMMPS, Zeo++, ORCA) y la que produce los artefactos
JSON que alimentan los fixtures. Debe ser determinista y sus salidas estables.

| Comando | Valida | Referencia | Notas |
|---|---|---|---|
| `embed`, `parse`, `rmsd`, `sssr`, `stereo`, `ecfp`, `tanimoto`, `charges`, `sasa`, `uff`, `population`, `dipole` | superficie core | RDKit/Open Babel (ver §4) | salida XYZ/SDF/JSON reproducible con `--seed` |
| `eht`, `dos`, `esp`, `pm3`, `xtb`, `gfn1`, `gfn2`, `hf3c`, `ani`, `ir`/`vibrational` (según flags), `nmr` | motores electrónicos | YAeHMOP/MOPAC/xtb/tblite/ORCA/CCCBDB | comparar JSON campo a campo con la salida del binario externo |
| `solvation`, `alpb`, `eeq`, `d4`, `edl-*`, `cpm`, `cpm-surface` | propiedades electrostáticas | xtb/dftd4/LAMMPS/analítico | `validate-electron-count` como L0 |
| `assemble`, `cell`, `kmesh`, `build-bloch-hamiltonian`, `compute-periodic-dos`, `solve-periodic-randnla` | materiales/periódico | pymatgen/spglib/pybinding | charts + cif round-trip |
| `neb-*`, `gsm-*`, `simplified-neb-path`, `irc` (si existe), `htst-rate`, `solve-microkinetic-*`, `analyze-gsm-mbh-htst-step` | reacciones/cinética | ASE NEB, geomeTRIC, pyGSM, Cantera | `--methods uff,pm3,xtb` por backend |
| `*-chart`, `chart-to-json`, `chart-from-json`, `validate-chart-schema`, `report-generator`, `pipeline-validator`, `validate-experimental-result`, `benchmark-function`, `trace-function-calls` | infraestructura | JSON Schema + reporter | round-trip y validación de esquema |
| `merge-experiment-results`, `experimental-result-to-json/sdf`, `cache-results` | herramientas de experimento | — | idempotencia y determinismo |

**Test equivalente global**
1. **L1 contra la librería:** cada comando debe reproducir su función Rust equivalente (llamada desde un test que lee el mismo input). Proponer `tests/regression/test_cli_parity.rs` que invoque el binario (`assert_cmd`/`std::process::Command`) y compare con la API Rust.
2. **Golden outputs:** congelar la salida JSON de ~50 comandos sobre `tests/fixtures/benchmark_100_molecules.json` y comparar byte a byte (con normalización de tiempos). Cualquier cambio de formato pasa por revisión del fixture.
3. **Comparación externa:** los comandos `pm3`, `xtb`, `gfn1`, `gfn2`, `d4`, `eeq`, `alpb`, `assemble`, `kmesh` son los puntos naturales donde el script Python llama al binario externo y al CLI, y diffea JSON.
4. **Códigos de error:** cada comando debe fallar con exit code ≠ 0 y JSON `error` legible en inputs inválidos (smiles roto, coords inconsistentes, elementos no soportados).

**Criterio de aceptación del crate:** paridad L1 100 % en los comandos cubiertos; golden outputs
idénticos; 0 discrepancias frente a los binarios externos por encima de las tolerancias de §3.3.

---

## 6. Crate `sci-form` (bindings Python)

**Rol en la validación:** es el mejor lugar para los tests cruzados, porque todo el ecosistema de
referencia (RDKit, PySCF, mordred, scikit-learn, ASE, pymatgen, geomstats, cvxpy, pyarrow) es
Python. Estrategia: los scripts de comparación calculan la referencia y **sci-form** en el mismo
proceso, y emiten un JSON de desviaciones que después se consume como fixture de Rust.

| Aspecto | Test | Referencia |
|---|---|---|
| Superficie completa | `tests/integration/test_python_integration.py` (existe) + `test_cross_language.py` | consistencia con Rust |
| RDKit como oráculo | nuevo `tests/integration/test_rdkit_oracle.py`: por cada ficha RDKit-based (§4.1–4.3, 4.4 F21, F28, F59, F61) | RDKit |
| PySCF como oráculo | nuevo `tests/integration/test_pyscf_oracle.py`: EHT/HF/ERI/AO→MO/DFT/CIS | PySCF |
| Grimme stack | nuevo `tests/integration/test_grimme_oracle.py`: GFN1/GFN2 (tblite), D4 (dftd4), sTDA (stda), HF-3c (ORCA) | tblite/dftd4/stda/ORCA |
| ML/descriptores | `test_mordred_oracle.py`, `test_sklearn_oracle.py` | mordred/scikit-learn |
| Materiales | `test_pymatgen_oracle.py` (CIF, grupos espaciales, XRD, PBC), `test_zeopp_oracle.py` | pymatgen/spglib/Zeo++ |
| Numérico/matricial | `test_scipy_oracle.py` (eigsh, linalg, integrate), `test_cvxpy_sdr.py` | SciPy/cvxpy |
| Errores y tipos | errores mapeados a `PyRuntimeError`; arrays de entrada/salida con `validate_coords` | contrato PyO3 |
| Rendimiento | `bench_*` scripts + comparación con RDKit en tiempos (no gate de veracidad) | RDKit |

**Convención de tolerancias:** los scripts Python devuelven `{method, ref, metric, value, tolerance,
pass}` y un `summary.json`. Los tests Rust leen ese archivo; así una sola tolerancia vive en un
único sitio (§3.3) y el CI no depende de tener RDKit/PySCF instalado para pasar tests básicos
(los fixtures están congelados).

**Criterio:** cada ficha de §4 debe tener (a) un script generador en `scripts/`, (b) un fixture
versionado, (c) un test Rust que valide el fixture, (d) un test Python opcional que regenere la
referencia cuando el entorno lo permite.

---

## 7. Crate `sci-form-wasm`

**Rol en la validación:** paridad numérica y estabilidad de API web. No hay referencia externa
nueva; la referencia es Rust/Python.

| Aspecto | Test | Criterio |
|---|---|---|
| Paridad numérica | `tests/integration/test_wasm_integration.js` (existe) + golden JSON compartido con Python | idénticos a ≤ 1e-12 (misma semilla), salvo float32 declarado |
| JSON contracts | esquema de cada función `*_json` | round-trip parse/stringify estable; campos requeridos |
| Typed arrays | `Float64Array`/`Float32Array` zero-copy (`positions_ptr`, `pack_*`) | longitudes exactas; contenido = versión JSON |
| Subpaths `/alpha` y `/beta` | build con features y smoke de cada export | todos los símbolos presentes; errores legibles |
| Determinismo multihilo | build con `parallel` vs sin él | mismos resultados (± 1e-12) |
| GPU/WebGPU | `tests/integration/test_webgpu_harness.js`, `compare_webgpu_reports.js` (existen) | política de fallback respetada, paridad CPU |
| Memoria/panic | inputs malformados (arrays de longitud ≠ 3N, JSON inválido) | error JS controlado, sin `panic!`/abort |

**Criterio:** 100 % de los 180 exports con smoke test; paridad con Python en el corpus GDB-20
(≥ 100 moléculas) a la tolerancia declarada.

---

## 8. Programas de referencia y datasets descargables

> **Licencias:** antes de añadir cualquier herramienta o dataset a CI, verificar licencia y términos
> de uso. Para CI público, priorizar RDKit, PySCF, SciPy, scikit-learn, ASE, pymatgen, spglib,
> OpenMM, pyarrow, mordred (permisivas). MOPAC, ORCA, Multiwfn, Zeo++, SDBS, nmrshiftdb2 tienen
> condiciones específicas: usarlas **offline**, en un entorno controlado, y solo para generar
> fixtures que no reproduzcan datos protegidos.

### 8.1 Programas de referencia

| Programa | Para qué (fichas) | Instalación | Licencia (verificar) |
|---|---|---|---|
| **RDKit** | SMILES/SMARTS/SMIRKS, ETKDG, UFF, MMFF94, descriptors, ECFP, Butina, features, SASA | `pip install rdkit` | BSD-3 |
| **Open Babel** | segundo voto parseo, UFF, Gasteiger, formatos | `conda install -c conda-forge openbabel` | GPL-2 |
| **PySCF** | integrales, RHF/UHF/ROHF, DFT, ECP, TDA/CIS, cubos ESP, poblaciones | `pip install pyscf` | Apache-2.0 |
| **ORCA** | HF-3c, GIAO, TD-DFT, frecuencias DFT, ECP | binario (registro académico) | freeware académico |
| **OpenMOPAC / MOPAC2016** | PM3/PM6 (ΔHf, FORCE) | binario/`openmopac` | open source/freeware (verificar) |
| **xtb** | GFN0/GFN1/GFN2 de referencia, ALPB, EEQ | `conda install -c conda-forge xtb` o binario | LGPL-3 |
| **tblite** | GFN1/GFN2 como biblioteca | `pip install tblite` (requiere libgfortran) | LGPL-3 |
| **dftd4 / simple-dftd3** | D3/D4 | `conda install -c conda-forge dftd4 simple-dftd3` | LGPL (verificar) |
| **stda** | sTDA canónico | GitHub `grimme-lab/stda` | verificar |
| **ASE** | NEB, PBC, mínima imagen, integradores | `pip install ase` | LGPL-2.1+ |
| **geomeTRIC** | Elastic Band, optimización en internas | `pip install geometric` | verificar |
| **pysisyphus** | GSM/NEB/IRC de contraste | `pip install pysisyphus` | verificar |
| **pyGSM** | GSM canónico | GitHub `ZimmermanGroup/pyGSM` | verificar |
| **LAMMPS** | ReaxFF, fix qeq, fix electrode/conp, NVT/NPT | conda/binary | GPL-2 |
| **OpenMM** | MD, GBSA, referencia determinista | `conda install -c conda-forge openmm` | MIT |
| **FreeSASA** | SASA | `pip install freesasa` | MIT |
| **Multiwfn** | DOS, ESP, Fukui, órdenes | binario (registro) | freeware académico |
| **YAeHMOP** | EHT canónico (matrices, Mulliken) | SourceForge (compilar) | open source (verificar) |
| **pyEHT** | EHT en Python | `pip install pyEHT` | verificar |
| **Multiwfn** | (ver arriba) | | |
| **JANPA** | NPA | binario/JAR (SourceForge) | open source (verificar) |
| **TorchANI / n2p2** | ANI AEV, MLFF BP | `pip install torchani` | MIT |
| **scikit-learn / SciPy / scikit-image** | RF/GB, optimización, MC | pip | BSD-3 |
| **mordred / PaDEL** | WHIM/GETAWAY/RDF/descriptores | `pip install mordred` | BSD-3 / libre |
| **geomstats / pymanopt** | Riemannian | pip | MIT / BSD (verificar) |
| **cvxpy** | SDR | pip | Apache-2.0 |
| **clifford / galgebra** | CGA | pip | MIT/BSD |
| **pybinding** | TB y KPM | `pip install pybinding` | BSD (verificar) |
| **pymatgen / spglib / pyxtal** | cristalografía | pip/conda | BSD/MIT |
| **Zeo++ / pormake** | porosidad / MOF builder | Zeo++ binario; `pip install pormake` | académico / MIT |
| **Cantera / CatMAP** | microcinética | pip (Cantera) | BSD-3 / verificar |
| **nmrglue** | espectros NMR | pip | BSD (verificar) |
| **pyarrow** | Arrow | pip | Apache-2.0 |
| **VMD / NAMD** | IMD | binario (registro) | freeware académico |

### 8.2 Datasets (descarga)

| Dataset | Uso (fichas) | Contenido | URL de referencia | Tamaño orientativo | Términos |
|---|---|---|---|---|---|
| **GDB-20 / GDB-20 subset** | F09, F12, F16, F19, F34 | 20k moléculas pequeñas SMILES | http://gdb.unibe.ch/downloads/ | ~2 MB (SMILES) | académico |
| `tests/fixtures/gdb20_reference_1k.json` | ya en repo | 1k geometrías de referencia RDKit | — | ~5 MB gz | interno |
| **GEOM** | F09, F12, F15 | 37M conformeros + energías GFN2/DMC | https://github.com/learningmatter-mit/geom | TB (usar subset) | verificar |
| **RDKit 10k reference** | ya en repo | poses RDKit ChEMBL | — | ~30 MB gz | interno |
| **ChEMBL (SMILES)** | transversal | 10k + 1k en `data/` | https://ftp.ebi.ac.uk/pub/databases/chembl/ | MB–GB | CC-BY-SA |
| **PubChem3D** | F09, F12 | conformeros PubChem | https://ftp.ncbi.nlm.nih.gov/pubchem/ | TB | dominio público (US) |
| **QM9** | F16, F35, F41, F48 | 134k moléculas, DFT B3LYP: energía, dipolo, freqs, termo | figshare (DOI 10.6084/m9.figshare.978904) | ~200 MB | CC0 |
| **QM7b** | F57 | 7.2k moléculas + 14 propiedades | figshare (QM7b) | ~5 MB | CC0 (verificar) |
| **tmQM** | F07, F22, F44, F45, F63 | 86k complejos TM + DFT | https://github.com/uiocompcat/tmQM | ~1–2 GB | CC (verificar) |
| **PubChemQC PM6** | F42 | 200M+ compuestos con PM6 | https://nakatamala.github.io/pubchemqc/ | TB (usar subset CID) | público |
| **S66 / S66x8** | F44, F45, F46 | 66 dímeros + curvas | repos de benchmark (Bonn/mirrors) | MB | académico |
| **WATER27, HAL59, MOR41, GMTKN55** | F44–F46, F87 | no covalente, halógenos, TM, química general | https://www.chemie.uni-bonn.de/grimme (GMTKN55) | MB–GB | académico |
| **Transition1x** | F75, F82, F87, F89 | 9.4M geometrías de reacciones + TS DFT | figshare/Zenodo (buscar "Transition1x") | ~50 GB (subset) | verificar |
| **ANI-1x / ANI-1ccx / ANI-2x** | F63, F79 | 5M/500k/8.9M conformeros DFT | figshare (repos isayev/ASE_ANI) | GB | verificar |
| **SPICE** | F63, F79 | 1M+ conformeros con fuerzas | https://github.com/openmm/spice-dataset | GB | CC (verificar) |
| **MD17 / sGDML** | F63, F71, F79 | trayectorias DFT de 10 moléculas | http://quantum-machine.org/gdml/ | GB | académico |
| **OMol25 / OMat24** | F77, F79 | 100M+ DFT moléculas/materiales | HuggingFace (Meta FAIR) | TB | verificar |
| **ESOL / FreeSolv / Lipophilicity / AqSolDB** | F29, F56, F57, F60, F61, F62 | solubilidad, hidratación, logP, pKa | https://deepchemdata.s3-us-west-1.amazonaws.com/ ; https://github.com/MobleyLab/FreeSolv | MB | libre (verificar) |
| **SAMPL6/7/8** | F62 | pKa/logP ciegos experimentales | https://github.com/samplchallenges | MB | libre |
| **DUD-E** | F58 | activos/inactivos por diana | http://dude.docking.org/ | GB | libre |
| **NIST CCCBDB** | F26, F31, F36, F40, F48 | geometrías, energías, freqs calculadas y experimentales | https://cccbdb.nist.gov/ | MB (export) | US Gov |
| **NIST WebBook** | F26, F48, F49, F50 | espectros y termoquímica experimental | https://webbook.nist.gov/chemistry/ | — | US Gov |
| **SDBS** | F52, F53, F54, F55 | espectros NMR/IR/MS experimentales | https://sdbs.db.aist.go.jp/ | — | uso investigación |
| **nmrshiftdb2** | F52–F55 | shifts/J experimentales y predicción HOSE | https://nmrshiftdb.nmr.uni-koeln.de/ | ~GB | open data (verificar) |
| **NMReDATA** | F52, F53 | asignaciones NMR publicadas | https://nmredata.org/ | MB | abierto |
| **NS372** | F52 | 372 shieldings de referencia | benchmark de Gauss (buscar "NS372 shielding set") | KB | académico |
| **PhotochemCAD** | F50, F51 | espectros UV-Vis y estándares de fluorescencia | https://www.photochemcad.com/ | MB | freeware (verificar) |
| **COD / AMCSD** | F64, F65, F66, F69 | CIFs y estructuras minerales | https://www.crystallography.net/cod/ ; http://rruff.geo.arizona.edu/AMS/amcsd.php | GB | abierto |
| **CoRE MOF 2019 / QMOF / ARC-MOF / hMOF** | F67, F68, F70 | estructuras MOF + propiedades | Zenodo/GitHub respectivos | GB | abierto |
| **IZA / RCSR** | F67 | zeolitas y topologías | http://www.iza-structure.org/databases/ ; http://rcsr.anu.edu.au/ | MB | abierto |
| **MoleculeNet (BBBP, Tox21, …)** | F56, F58, F59 | tareas de propiedad | https://moleculenet.org/ | MB | varios |

### 8.3 Procedimiento de descarga propuesto

1. Crear `scripts/download_validation_datasets.sh` (y un `.sha256` por dataset) que descargue a
   `data/external/<dataset>/` (ya excluido de `Cargo.toml` packaging por `data/`), con modo
   `--small` (solo lo necesario para CI) y `--full`.
2. Los fixtures derivados (JSON/gz) van a `tests/fixtures/` y **sí** se commitean (comprimidos),
   siguiendo la práctica existente (`gdb20_reference.json.gz`).
3. Cada script de comparación imprime un resumen con la forma `{ficha, ref, metric, value,
   tolerance, pass}` agregado en `reports/validation/summary.json` (crear carpeta `reports/`,
   ignorada por git salvo el `summary`).
4. Ningún test de CI descarga > 50 MB ni invoca software con licencia restrictiva por defecto:
   esos tests quedan tras `--features validation-external` o `#[ignore]` documentado.

## 9. Matriz de trazabilidad funcionalidad → test → dataset

Leyenda: **Actual** = suite que ya existe en el repo; **Nuevo** = artefacto a crear.
`gen` = script generador en `scripts/`; `fix` = fixture en `tests/fixtures/`; `test` = suite Rust
`tests/regression/`.

### 9.1 Parsing y topología

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F01 | SMILES parser | `test_diverse_molecules.rs`, `src/smiles.rs` | gen `compare_pipeline.py` → fix `parse_compare.json` → test `test_smiles_cross` | RDKit, Open Babel; ChEMBL 10k |
| F02 | SMILES canónico | `test_algorithm_validation.rs` | gen `compare_canonical_smiles.py` | RDKit (equivalencia, no string) |
| F03 | SMARTS | `test_new_features.rs`, `smarts_patterns.txt` | gen `compare_smarts.py` → fix → test `test_smarts_cross` | RDKit; Daylight patterns |
| F04 | SMIRKS | `test_smirks_reactions.rs`, `compare_smirks_reactions.py` | ampliar catálogo 100+ plantillas | RDKit `RunReactants` |
| F05 | SSSR | `test_algorithm_validation.rs` | gen `compare_rings.py` → fix → test `test_rings_cross` | RDKit `GetSymmSSSR`; ChEMBL |
| F06 | Aromaticidad/estereo | `test_tet_centers.rs`, `test_algorithm_validation.rs` | gen `compare_stereo.py` → fix → test `test_stereo_cross` | RDKit CIP; casos atropisómeros |
| F07 | Topología/hapticidad | `test_metal_experimental_geometry.rs` | gen `compare_topology_tmqm.py` → fix → test `test_topology_cross` | tmQM; `metal_experimental_geometry.json` |

### 9.2 Geometría 3D

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F08 | Bounds matrix | `test_bounds_compare.rs`, `dump_rdkit_bounds.py` | gen `generate_bounds_reference.py` → fix gz → test `test_bounds_cross` | RDKit bounds; GDB-20 |
| F09 | ETKDG embedding | `test_gdb20.rs`, `test_gdb20_rmsd.rs`, `generate_gdb20_reference.py` | ampliar a 10k y GEOM subset | RDKit ETKDGv3; GDB-20; GEOM |
| F10 | Torsiones ETKDG-3D | `test_grad_torsion_detail.rs`, `extract_torsions*.py` | gen `compare_torsion_energies.py` → fix | RDKit torsions; `torsion_params.json` |
| F11 | Validación geométrica | `test_diverse_molecules.rs` | test de sensibilidad (poses corruptas) | CSD medios; RDKit checks |
| F12 | Búsqueda de conformeros/ranking | `test_ensemble_rmsd.rs`, `embed_diverse` | gen `compare_ensemble_ranking.py` → fix | RDKit `EmbedMultipleConfs`; GEOM |
| F13 | RMSD/Kabsch | `test_algorithm_validation.rs`, `rmsd_per_mol.rs` | gen `compare_rmsd_rdkit.py` → test `test_rmsd_cross` | RDKit `GetBestRMS`; SciPy |
| F14 | Butina/diversidad | `test_ensemble_rmsd.rs` | gen `compare_butina.py` → fix → test | RDKit Butina |
| F15 | Strain y torsión scan | `torsion_scan.rs`, `docs/strain-energy.md` | gen `compare_strain.py` | RDKit MMFF/UFF; GEOM |

### 9.3 Fuerzas y optimización

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F16 | UFF | `test_gradient_check.rs` | gen `compare_uff_rdkit.py` → fix gz → test `test_uff_cross` | RDKit UFF; GDB-20, QM9 |
| F17 | MMFF94 | `test_gradient_check.rs` | gen `compare_mmff94_rdkit.py` → fix → test | RDKit MMFF94 (+ suite Halgren) |
| F18 | Optimizadores/restricciones | `test_step_by_step.rs`, `potentials.rs` | gen `compare_lbfgs_scipy.py` | SciPy L-BFGS-B; geomeTRIC |
| F19 | Rotámeros/torsion scan | `debug_torsions.rs`, `torsion_scan.rs` | gen `compare_torsion_profiles.py` | RDKit `SetDihedralDeg`; CSD Torsion Library |
| F20 | Gradientes vs numérico | `test_gradient_check.rs`, `test_gradient_perterm.rs` | matriz por motor `test_gradients_matrix.rs` | diferencias finitas (h=1e-3) |

### 9.4 Cargas, electrostática y solvatación

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F21 | Gasteiger | `test_charges.rs`, `generate_charges_reference.py` | gen `compare_gasteiger_rdkit.py` → test | RDKit `ComputeGasteigerCharges` |
| F22 | EEQ (+EEQ-BC) | `test_eeq.rs`, `test_session_algorithms.rs` §10, `cross_validation.rs` (f22_eeqbc) | gen `compare_eeq_xtb.py` → fix; `oracle_eeqbc.py` (multicharge CLI) | xtb/tblite EEQ; multicharge `--model eeqbc` |
| F23 | Mulliken/Löwdin | `test_eht.rs`, `test_phase_c.rs`, `generate_population_reference.py` | gen `compare_population_ref.py` → fix | YAeHMOP; PySCF |
| F24 | NBO/NPA | `src/population/npa.rs` | gen `janpa_reference.py` → fix → test `test_npa_cross` | JANPA; NBO7 |
| F25 | Wiberg/Mayer | `test_session_algorithms.rs` (parallel=serial) | gen `generate_bond_orders_reference.py` | Multiwfn; PySCF |
| F26 | Dipolo | `test_extended_molecules.rs::test_dipole_moments_vs_nist` | ampliar serie NIST | NIST CCCBDB/WebBook |
| F27 | ESP/cube | `test_phase_c.rs`, `src/esp` | `generate_esp_cube.py` (documentado) + compare | PySCF cubegen; Multiwfn |
| F28 | SASA | `test_sasa.rs` | gen `compare_sasa_freesasa.py` → test | FreeSASA; RDKit |
| F29 | Solvatación no polar/GB | `test_extended_molecules.rs::test_solvation_*` | gen `compare_gb_openmm.py` | OpenMM GBSAOBC; FreeSolv |
| F30 | ALPB | `test_alpb.rs` | gen `compare_alpb_xtb.py` → fix | xtb ALPB |

### 9.5 Estructura electrónica

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F31 | EHT solver S/H/eig | `test_eht.rs`, `test_eht_metal_references.rs`, `generate_eht_reference.py` | compare con YAeHMOP | YAeHMOP; pyEHT; VSIP paper |
| F32 | Gradiente EHT/opt | `test_gradient_check.rs`, `check_grad.rs` | fix geometrías EHT óptimas | diferencias finitas; CCCBDB |
| F33 | Rejilla/isosurface/MC | `test_eht.rs`, `src/mesh` | gen `compare_marching_cubes.py` | scikit-image; VMD |
| F34 | Bandas/k-mesh | `test_kpm.rs` (parte), CLI charts | gen `compare_band_structure.py` | pybinding; pymatgen; cadenas analíticas |
| F35 | DOS/PDOS/multimétodo | `test_session_algorithms.rs`, `export_dos_json` | `generate_dos_reference.py` (documentado) + compare | Multiwfn; PubChemQC |
| F36 | SCF RHF/integrales | `test_experimental_comparison.rs` (54 tests), `test_phase_c.rs` | gen `generate_scf_reference.py` → fix | PySCF; CCCBDB |
| F37 | ERIs/Boys/screening | `test_alpha_scf_eri_accuracy.rs`, bins `check_c_ints` | gen `compare_eri_pyscf.py` | PySCF/libcint |
| F38 | UHF/ROHF/⟨S²⟩ | `test_session_algorithms.rs` §1 | gen `generate_uhf_reference.py` | PySCF; Szabo-Ostlund |
| F39 | AO→MO | `test_session_algorithms.rs` §3 | gen `generate_ao2mo_reference.py` | PySCF `ao2mo` |
| F40 | HF-3c (+D3/gCP/SRB) | `test_experimental_comparison.rs` (bloque HF-3c) | gen `generate_hf3c_orca_reference.py` → fix | ORCA HF-3c; paper Sure-Grimme |
| F41 | CIS/CISD | `test_experimental_comparison.rs` (spectroscopy) | gen `generate_tda_reference.py` | PySCF TDA; ORCA |
| F42 | PM3 | `test_session_algorithms.rs` §5, `test_extended_molecules.rs` | gen `compare_pm3_mopac.py` → fix | OpenMOPAC; PubChemQC PM6 |
| F43 | GFN0/xTB | `test_extended_molecules.rs`, §6 Broyden | gen `compare_xtb_tblite.py` | xtb; tblite |
| F44 | GFN1 | `test_extended_molecules.rs` (parcial) | gen `compare_gfn1_tblite.py` → fix | tblite; S66x8; WATER27 |
| F45 | GFN2 | `scripts/gfn2_*` (muchos, consolidar) | gen `compare_gfn2_tblite.py` + lock de versión | tblite; S66x8; HAL59; MOR41 |
| F46 | D3/D4 | `test_d4.rs` | gen `compare_dftd4.py` → fix | dftd4; simple-dftd3 |
| F47 | Fukui/frontera/reactividad | `test_new_features.rs`, `test_session_algorithms.rs` | gen `compare_fukui_multiwfn.py` | Multiwfn; series de orientación aromática |

### 9.6 Espectroscopía

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F48 | Hessiano/vib/termo | `test_spectroscopy.rs`, `test_extended_molecules.rs::test_ir_*` | gen `compare_hessian_mopac.py`, `compare_hessian_xtb.py` | MOPAC FORCE; xtb --hess; CCCBDB |
| F49 | IR spectrum/asignación | `test_spectroscopy.rs` | gen `compare_ir_peaks.py` | NIST WebBook; tablas IR |
| F50 | sTDA UV-Vis | `test_spectroscopy.rs`, `test_extended_molecules.rs` | gen `compare_stda_grimme.py`, `compare_tddft_orca.py` | stda; ORCA TD-DFT; NIST |
| F51 | Fluorescencia/Einstein A | `test_spectroscopy.rs`, `docs/fluorescence.md` | gen `check_einstein_a.py` | ORCA; estándares IUPAC/PhotochemCAD |
| F52 | GIAO NMR | `test_extended_molecules.rs` (ignorados), §7 | gen `generate_orca_giao_reference.py` → fix | ORCA GIAO; NS372; SDBS |
| F53 | HOSE shifts | `test_extended_molecules.rs::test_legacy_nmr_*` | gen `compare_hose_nmrshiftdb.py` → fix | nmrshiftdb2; SDBS; NMReDATA |
| F54 | J couplings | `test_session_algorithms.rs` §7 | gen `compare_j_couplings_nmrshiftdb.py` | nmrshiftdb2; ORCA; Karplus |
| F55 | Espectro NMR/broadening | `src/nmr/spectrum.rs` tests | gen `compare_nmrglue.py` | nmrglue; SDBS |

### 9.7 ML, descriptores y potenciales

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F56 | Descriptores 2D/Lipinski | `test_extended_molecules.rs::test_ml_*` | gen `compare_descriptors_rdkit.py` → test | RDKit; MoleculeNet |
| F57 | WHIM/GETAWAY/RDF | `test_new_features.rs` (parcial) | gen `compare_mordred.py` → fix | mordred; PaDEL |
| F58 | Farmacóforos | — | gen `compare_features_rdkit.py` | RDKit Pharm2D; DUD-E |
| F59 | ECFP/Tanimoto | `test_extended_molecules.rs::test_fingerprint_*` | gen `compare_ecfp_rdkit.py` → test XOR | RDKit Morgan/ECFP4 |
| F60 | RF/Gradient boosting | `test_new_features.rs` | gen `compare_sklearn_models.py` | scikit-learn; ESOL/Lipophilicity |
| F61 | Propiedades de fármaco/logS | `test_extended_molecules.rs::test_ml_logp_*`, `test_ml_lipinski_*` | gen `compare_logs_esol.py` | RDKit Crippen; ESOL; AqSolDB |
| F62 | pKa empírico | — | gen `compare_pka_sampl.py` → fix | SAMPL6/7/8; Chemprop |
| F63 | ANI-2x/ANI-TM | `benchmark_ani.py`, `alpha_mlff_consistency.rs` | gen `compare_ani_torchani.py` → fix | Torchani; ANI datasets; tmQM |

### 9.8 Materiales, dinámica, GPU y alpha/beta

| # | Funcionalidad | Actual | Nuevo | Referencia / dataset |
|---|---|---|---|---|
| F64 | Grupos espaciales | `test_new_features.rs` (parcial) | gen `compare_space_groups.py` → fix 230 grupos | pymatgen; pyxtal; spglib |
| F65 | CIF import/export | `test_session_algorithms.rs` §2 | gen `compare_cif_pymatgen.py` | pymatgen/ASE; COD |
| F66 | PBC/superceldas | `test_roadmap_validation.rs` (parcial) | gen `compare_pbc_ase.py` | ASE; pymatgen |
| F67 | Ensamblado MOF | `test_roadmap_validation.rs` | gen `compare_pormake.py` | pormake; CoRE MOF; RCSR |
| F68 | Porosidad | `docs/materials-assembly.md` | gen `compare_porosity_zeopp.py` | Zeo++; CoRE MOF |
| F69 | XRD de polvo | `simulate_powder_xrd` | gen `compare_xrd_pymatgen.py` | pymatgen XRDCalculator; COD |
| F70 | Optimización de frameworks | `test_new_features.rs` | gen `compare_framework_opt_ase.py` | ASE BFGS; CoRE MOF |
| F71 | MD NVE/NVT | `alpha_dynamics_conservation.rs` | gen `compare_md_openmm.py` | OpenMM; LAMMPS; analítico |
| F72 | Termostatos/barostats | `dynamics_live/*` tests | comparación LAMMPS | LAMMPS; Frenkel-Smit |
| F73 | IMD | manual | test con VMD | VMD; spec NAMD |
| F74 | Arrow/chunks/workers | `test_cross_language.py` | gen `compare_arrow_pyarrow.py` | pyarrow |
| F75 | Potenciales analíticos/EVB | `e3_justification.rs`, `reaction_dynamics_study.rs` | tests de mínimos publicados | Muller-Brown; EVB papers |
| F76 | GPU↔CPU | `test_gpu_candidates.rs`, webgpu harness | ampliar cobertura de kernels | CPU como referencia |
| F77 | DFT SVWN/PBE | `alpha_dft_validation.rs`, `alpha_cross_method.rs` | gen `generate_dft_reference.py` | PySCF/libxc; AE6/BH6/S22 |
| F78 | ReaxFF/EEM | `alpha_reaxff_validation.rs` | gen `compare_reaxff_lammps.py` | LAMMPS `reaxff`; PuReMD |
| F79 | MLFF/AEV | `alpha_mlff_consistency.rs` | gen `compare_aev_n2p2.py` | TorchANI/n2p2; MD17; SPICE |
| F80 | Obara-Saika alpha API | `alpha_scf_eri_accuracy.rs` | gen `compare_eri_pyscf.py` | PySCF/libcint |
| F81 | CGA | `test_cga.rs` | gen `check_cga_vs_quaternion.py` | clifford; cuaterniones |
| F82 | GSM | `test_gsm.rs`, `compare_reaction_layers.py` | gen `compare_gsm_pygsm.py` | pyGSM; pysisyphus; Transition1x |
| F83 | SDR | `test_sdr.rs` | gen `compare_sdr_cvxpy.py` | cvxpy; SFSDP |
| F84 | Dynamics-live | `alpha_dynamics_conservation.rs` | paridad con OpenMM | OpenMM; LAMMPS |
| F85 | IMD | manual | sesión VMD documentada | VMD; spec IMD |
| F86 | EDL | tests experimentales/CLI | gen `check_edl_analytics.py` | soluciones analíticas; Grahame |
| F87 | Cinética/HTST | tests experimental/CLI | gen `compare_microkinetics_cantera.py` | Cantera; CatMAP; Eyring |
| F88 | Periódico lineal/Bloch | `test_kpm.rs`, CLI | gen `compare_bloch_pybinding.py` | pybinding; pymatgen |
| F89 | Reaction dynamics 3D | `alpha_reaction_dynamics_3d.rs`, `compare_reaction_layers.py --external` | ampliar a Transition1x | ASE NEB; geomeTRIC; pysisyphus |
| F90 | Render bridge/charts | `test_gpu_candidates.rs`, webgpu compare | esquema JSON versionado | JSON Schema; Vega-Lite |
| F91 | KPM DOS | `test_kpm.rs` | gen `compare_kpm_pybinding.py` | pybinding KPM; `eigh` exacto |
| F92 | MBH | `test_mbh.rs` | Hessiano completo como referencia | Ghysels et al.; WATER27 |
| F93 | RandNLA | `test_randnla.rs` | gen `compare_eigsh.py` | scipy `eigsh`; `eigh` |
| F94 | Riemannian PSD | `test_riemannian.rs` | gen `compare_riemannian_geomstats.py` | geomstats; pymanopt |
| F95 | CPM | `test_cpm.rs` | gen `compare_cpm_lammps.py` | LAMMPS `fix electrode/conp` |

---

## 10. Plan de ejecución por fases

**Fase 0 — Infraestructura de validación (1 semana)**
- [ ] `scripts/validation/common.py`: utilidades de reporte (`{ficha, ref, metric, value, tolerance, pass}`), alineamiento por grafo, Kabsch, conversión de unidades.
- [ ] `scripts/download_validation_datasets.sh` con modo `--small`/`--full` y checksums.
- [ ] `reports/validation/summary.json` + `.gitignore` de los binarios.
- [ ] `validation/refs.lock.json`: versiones exactas de cada herramienta de referencia por ficha.

**Fase 1 — L0 + L1 (1–2 semanas, bloqueante)**
- [ ] Matriz de gradientes vs numérico por motor (F20) en CI.
- [ ] Paridad Rust/Python/WASM/CLI sobre 100 moléculas GDB-20 (F01–F13).
- [ ] Invariantes ya listados en §3 (sumas, integrales, simetrías, conservación).

**Fase 2 — L2 referencias canónicas del core (2–3 semanas)**
- [ ] Stack RDKit: F01, F03, F05, F06, F09, F10, F13, F14, F16, F17, F21, F28, F56, F59.
- [ ] Stack PySCF: F23, F25, F27, F36, F37, F38, F39, F41, F77.
- [ ] Stack Grimme: F22, F30, F42, F44, F45 (F43 done), F46, F78.
- [ ] Materiales: F64–F69 (pymatgen/spglib/Zeo++).

**Fase 3 — L4 datasets experimentales (2–4 semanas, en paralelo)**
- [ ] GDB-20 10k + GEOM subset para F09/F12.
- [ ] NIST/CCCBDB para F26/F36/F40/F48.
- [ ] nmrshiftdb2/SDBS para F52–F55.
- [ ] SAMPL para F62; ESOL/FreeSolv para F29/F61.
- [ ] tmQM para F07/F22/F63; Transition1x para F82/F87/F89.

**Fase 4 — Alpha/beta y sistemas (2–3 semanas)**
- [ ] LAMMPS (ReaxFF/CPM/termostatos), ORCA (HF-3c/GIAO/TD-DFT), pyGSM/pysisyphus (GSM/IRC), pybinding (KPM/bandas), cvxpy (SDR), geomstats (Riemannian), VMD (IMD).
- [ ] Formalizar las promociones alpha→beta con el criterio de esta fase.

**Fase 5 — Gates de CI y publicación de precisión**
- [ ] Separar CI en lanes: `quick` (L0/L1, < 5 min), `references` (L2, < 30 min), `external` (L3/L4, nightly/`--ignored`).
- [ ] Publicar tabla "supported accuracy" por método/observable (`docs/algorithm/accuracy.md`) generada desde `summary.json`.
- [ ] Cada release: comparar `summary.json` con el de la release anterior y fallar si una métrica cae > 2× su tolerancia.

---

## 11. Registro de resultados y CI

**Formato de reporte (por ficha):**

```json
{
  "ficha": "F45",
  "observable": "gfn2_total_energy",
  "reference": {"tool": "tblite", "version": "0.4.x", "command": "tblite single --method gfn2"},
  "units": "hartree",
  "n_cases": 120,
  "metric": "max_abs",
  "value": 0.0004,
  "tolerance": 0.0008,
  "pass": true,
  "notes": "S66x8 subset; misma geometría; conversión eV→Ha verificada"
}
```

Reglas:
1. **Un reporte por ficha y por observable**; nunca agregar energías de métodos distintos en una métrica.
2. **`n_cases` visible**: una tolerancia sin tamaño de muestra no es evidencia.
3. **Semilla y versión de referencia** siempre presentes.
4. **`pass` solo de la comparación L2/L4**; los L0/L1 son gates duros en CI y no necesitan reporte.
5. Guardar histórico por release (`reports/validation/<version>/summary.json`).

**CI propuesto:**
- `quick` (cada PR): `cargo test --lib`, `cargo test --release --test ci`, paridad L1 y grafo de fixtures pequeños.
- `references` (cada PR, ~30 min): tests que leen fixtures L2 congelados (no instalan RDKit/PySCF).
- `external` (nightly/etiqueta): regenera fixtures con las herramientas reales y diffea; falla si la desviación crece.
- `wasm`/`python` (release): smoke + paridad.

---

## 12. Riesgos, trampas y reglas anti-falso-positivo

| Riesgo | Ejemplo real en el repo | Regla |
|---|---|---|
| **Unidades incompatibles** | HF-3c vs "experimental SCF" (3031 mHa) documentado en el roadmap | declarar unidades por observable y convertirlas antes de comparar |
| **Métodos distintos comparados como iguales** | EHT vs RHF vs GFN0 | L3 solo para tendencias; jamás como gate de corrección |
| **Definiciones distintas del mismo nombre** | Mulliken (EHT vs HF), Löwdin (S^{1/2} vs S^{-1/2} convención), Gasteiger (6 vs 12 iteraciones), strain | fijar la definición y la implementación de referencia en la ficha |
| **No unicidad** | SSSR, canonical SMILES, cluster labels, autovectores degenerados, modos normales | comparar invariantes (multisets, subespacios, proyecciones) |
| **Orden de átomos/orientación** | fixtures RDKit con orden distinto | alinear por isomorfismo/substruct match; Kabsch para geometría |
| **Firmas de autovectores** | eigensolvers devuelven ±C | comparar `|C|` o subespacios |
| **float32/float64** | TorchANI/WASM typed arrays | tolerancia por precisión declarada; no usar 1e-12 en rutas f32 |
| **Versión de parámetros** | GFN2, D4, MMFF94/MMFF94s, ETKDGv2/v3, PM3 vs PM3(tm) | `refs.lock.json` con versión exacta; regenerar fixtures al actualizar |
| **Sobreajuste a fixtures** | copiar el valor esperado en el código | el fixture se genera con script reproducible; revisar el script en el PR |
| **Datasets de entrenamiento usados como test** | ANI-1x train vs test | usar splits oficiales de test y documentarlo |
| **Licencias** | SDBS, nmrshiftdb2, MOPAC, Zeo++ | solo en lane `external`; no redistribuir contenido protegido |
| **Estadística insuficiente** | 5 moléculas para validar un método | mínimos por ficha: 20 casos para tendencias, 100 para tolerancias numéricas |
| **CI con herramientas externas** | instalación pesada de ORCA/LAMMPS | fixtures congelados en `tests/fixtures/`; lane nightly para regenerar |
| **No determinismo de GPU** | reducciones paralelas | comparar con tolerancia relativa, no bit a bit |
| **Mezcla de features** | resultados con `parallel` y sin él | fijar la matriz de features por lane y declarar el comando exacto |

**Regla final:** una funcionalidad se considera *verificada* solo cuando tiene (1) invariantes L0,
(2) paridad L1 entre superficies, (3) referencia L2 con tolerancia documentada, y (4) un dataset
L4 o experimento que acote su error real. Cualquier ficha sin los cuatro niveles queda marcada
como *parcial* en la tabla "supported accuracy".

---

## 13. Estado de implementación

> Actualizado: **2026-10-12** (GFN1-xTB fiel, F44). Ejecutar
> `VALIDATION_PYTHON=<venv>/bin/python scripts/run_validation.sh` para reproducir
> toda la tabla. Reporte generado: `reports/validation/summary.json`.

### 13.1 Fichas implementadas y medidas (harness actual)

| Ficha | Referencia | Métrica | Medido | Tolerancia | Estado |
|---|---|---|---|---|---|
| F01 | RDKit (topología sin sanitizar + AddHs) | casos con mismatch | **0 / 159** | 0 | PASS |
| F05 | RDKit `GetSSSR` (estricto, Figueras) | casos con mismatch | **0 / 159** | 0 | PASS |
| F13 | NumPy SVD Kabsch | max abs | **2.2e-15 Å** | 1e-6 | PASS |
| F14 | RDKit Butina (reordering=True) | casos con mismatch | **0 / 12** | 0 | PASS |
| F16 | RDKit UFF | max abs (400 moléculas) | **0.0000 kcal/mol** | 0.01 | PASS |
| F17 | RDKit UFF AtomTyper | etiquetas con mismatch | **0 / 33641** | 0 | PASS |
| F21 | RDKit Gasteiger nIter=6 | max abs | **8.3e-17 e** | 0.02 | PASS |
| F23 | PySCF RHF/STO-3G Mulliken | max abs | **8.6e-8 e** | 0.005 | PASS |
| F26 | PySCF RHF/STO-3G dipolo Mulliken | max abs por componente | **1.8e-7 D**; dirección 0.0000° | 1e-6 | PASS |
| F49 | NIST/CCCBDB bandas IR experimentales | recall de grupo (10 moléculas) | **1.000** | ≥ 0.90 | PASS |
| F70 | RDKit `UFFOptimizeMolecule` (misma superficie F16) | |ΔE|/átomo mismo mínimo | **7.8e-10 kcal/mol**; RMSD 3.4e-6 Å | 1e-6 | PASS |
| F67 | pormake RCSR `.cgd` + Builder MOF-5 | grados 16/16; MOF-5 53 átomos, nudo 1.500 Å | **0 mismatches** | 0 | PASS |
| F68 | pyzeo (Zeo++, mismos radios) MOF-5 | VF/LCD/PLD + analíticos | **ΔVF 8e-4, ΔLCD 0.004 Å, ΔPLD 0.29 Å** | 0.01/0.3/0.5 | PASS |
| F28 | FreeSASA Shrake-Rupley (mismos radios/sonda/960 puntos) | max rel | **1.53 %** | 2 % | PASS |
| F36 | PySCF RHF/STO-3G (S/T/V, E_tot, MO) | S/T/V 3.1e-10 a.u.; E_tot 1.2e-8 Ha; MO 1.5e-7 Ha | 11 moléculas | 1e-9 / 1e-7 / 1e-6 | PASS |
| F37 | PySCF int2e | max abs (13 074 ERIs) | **2.1e-11 a.u.** | 1e-10 | PASS |
| F38 | PySCF UHF (misma conjetura H_core) | max abs E (8 casos abiertos) | **6.5e-8 Ha**; Δ⟨S²⟩ 5.7e-9 | 1e-6 | PASS |
| F41 | PySCF TDA (singletes, misma base) | max abs E_exc (29 estados, 8 moléculas) | **2.6e-6 eV**; E_HF 1.2e-8 Ha | 0.01 eV | PASS |
| F39 | PySCF ao2mo.full (misma C y ERIs AO) | max abs (13 058 MO ERIs) | **1.9e-14 a.u.** | 1e-10 | PASS |
| F43 | xtb GFN0-xTB 6.7.1 | max abs E (18 moléculas) | **1.7e-4 eV** (HOMO 7e-5, LUMO 1e-4) | 0.01 eV | PASS |
| F44 | tblite GFN1-xTB 0.7.0 | max abs E (18 moléculas) | **8.1e-5 eV** (HOMO 1.8e-5, LUMO 4.2e-5) | 0.01 eV | PASS |
| F45 | tblite GFN2-xTB | max abs E (18 moléculas) | **2.3e-5 eV** (HOMO/LUMO ~1e-6 eV) | 0.0217 eV | PASS |
| F46 | s-dftd3 (D3 HF-3c) / dftd4 (D4 PBE) | max abs E (18 moléculas) | D3 **6.8e-9 Ha** ; D4 **1.2e-10 Ha** | 1e-3 kcal/mol | PASS |
| F59 | RDKit Morgan/ECFP4 | Spearman Tanimoto (780 pares) | **0.756** | ≥ 0.7 | PASS |
| F64 | pymatgen/spglib (tablas ITC) | 230/230 grupos: orden, símbolo, sistema, red, grupo puntual y set de operaciones | **0 discrepancias** | 0 | PASS |
| F65 | pymatgen CIF | parámetros de red (4 estructuras) + round-trip | **0.0** Å/grados | 1e-6 | PASS |
| F66 | ASE `find_mic` | distancias de mínima imagen (44 pares) | **8.9e-16 Å** | 1e-10 | PASS |
| F69 | pymatgen XRDCalculator (CuKα) | posiciones 2θ (4 estructuras) + corr. intensidades | **Δ2θ=0.0000°; corr ≥ 0.9982** | 0.01° / ≥0.95 | PASS |
| F47 | reglas EAS experimentales (+DFT) | top-1 f− en serie aromática (8 moléculas) | **7/8 (0.875)** | ≥ 0.80 | PASS |
| F48 | xtb `vibspectrum` (misma matriz H) | 39 frecuencias (6 moléculas) | **0.009 cm⁻¹** | 2.0 | PASS |
| F42 | MOPAC 23.2.5 PM3 | ΔHf (24 moléculas) + H/E_nuc/SCF de agua | **exacto** en hidruros X–H (<0.02 kcal/mol), H y E_nuc (1e-4 eV), SCF de agua (eigenvalores/cargas exactos); residual π documentado (≤130 kcal/mol) | 0.02 (hidruros) | PARCIAL |
| F41 (ext.) | PySCF `int1e_r` | integrales de dipolo AO (8 moléculas) | **1.1e-9 a.u.** | 1e-8 | PASS |
| F41 (ext.) | PySCF TDA osciladores | fuerzas de oscilador con dipolos exactos | **<1e-5 en 27/29 estados** (H2O 0/1 documentado) | 5e-3 | PASS |
| F46 (ext.) | dftd4 4.3.0 | 7 casos pesados (PbI2…LrH, UO2) | **6.2e-10 Ha** | 1.6e-6 | PASS |
| F67 (ext.) | pormake `.cgd` | 6 nets RCSR con aristas no equivalentes | **nodos/CN/tipos exactos** | exacto | PASS |

Fixtures: `tests/fixtures/validation/` (README con procedencia). Versiones de
referencia: `validation/refs.lock.json` (RDKit 2026.09.1, PySCF 2.14.0, NumPy 2.5.4,
SciPy 1.18.1, FreeSASA python-bindings, tblite 0.7.0, dftd4 4.3.0, s-dftd3 1.6.0), pymatgen 2026.9.24, ASE 3.29.0, spglib 2.8.0).
Los 28 tests de `tests/cross_validation.rs`
están activos (ningún `#[ignore]`).

### 13.2 Correcciones aplicadas gracias a la comparación (todo con test)

1. **Gasteiger exacto (F21)** — `src/charges/gasteiger.rs`:
   - tabla por hibridación/modo (canónica: C/N/O/S/P/B/Si/Be/Mg/Al),
   - iteración PEOE real (`ionX` del átomo donante, hidrógeno `IONXH = 20.02`,
     amortiguamiento 0.5 que se reduce a la mitad),
   - `splitChargeConjugated` para sistemas conjugados cargados,
   - la desviación máxima frente a RDKit pasó de **0.153 e a 8e-17 e**.
2. **Hibridación del parser (F01/F16/F21)** — `src/smiles.rs`:
   - S/P hipervalentes (sulfóxido/sulfona/sulfato/fosfato) → SP3; SO₂ → SP2;
   - B trigonal → SP2; conjugación que no se propaga a través de S/P
     (sulfonamida N SP3, O de éster SP3, igual que RDKit);
   - comparado contra RDKit en una batería de 12 entornos: 100 % coincidente.
3. **Enlaces explícitos entre aromáticos (F01)** — un `-` explícito entre dos
   anillos aromáticos (bifenilo) ya no se promociona a aromático: antes 13
   aromáticos + 10 simples, ahora 12 + 11 como RDKit.
4. **Hidrógenos de átomos entre corchetes (F01)** — nuevo `Atom::no_implicit_h`:
   `[Si]` → 0 H, `[NH4+]` → 4 H, `[CH2]` → 2 H (semántica SMILES correcta).
5. **Fragments/contraiones** — documentado y testeado: `SmilesParser` los
   conserva; `Molecule::from_smiles` se queda con el fragmento mayor por diseño
   (los oráculos F01/F05/F21 lo excluyen y lo registran en el fixture).
6. **SSSR (F05)** — `src/rings/sssr.rs`: el conjunto de candidatos (ciclo más
   corto por arista) podía no generar el espacio de ciclos completo; se añaden
   los ciclos fundamentales del bosque de expansión. Caso cubanoide
   `C12C3C4C1C5C3C2C45`: de 4 anillos a los 5 del rango de ciclos.
7. **Butina (F14)** — `src/clustering.rs` reescrito a la variante de RDKit
   (`reordering=True`, orden de semillas `(n_vecinos, idx)` descendente). El
   oráculo tenía además la matriz de distancias en orden tril incorrecto
   (`np.tril_indices`); con ambos arreglos la partición coincide 12/12.
8. **SASA (F28)** — la implementación Shrake-Rupley era **correcta**; la
   referencia anterior (`rdFreeSASA.CalcSASA` + `classifyAtoms`) devolvía radios
   cero para moléculas pequeñas y daba un oráculo inválido (~2× de diferencia).
   El oráculo usa ahora FreeSASA canónico con los mismos radios/sonda/puntos:
   1.53 % máximo.
9. **Paridad UFF completa (F16/F17)** — reescritura contra las fuentes
   canónicas de RDKit: tabla de parámetros exacta (`Params.cpp`, con `U1`),
   `ka` con `Z1/r13` (`beta·Z1_i·Z1_k/r13⁵·r12·r23·inner`), término de ángulo
   por coordinación (lineal `1+cosθ`, trigonal `1−cos3θ`, cuadrado `1−cos4θ`)
   con correcciones `θ0` de anillos 3/4 (`order` 30/35/40/45), torsiones con
   escalado por número de caminos y casos especiales del grupo 6/`equation17`,
   vdW con media geométrica y sin escalado 1-4, inversión con coeficientes
   `C0/C1/C2`, `K/3` y tres contribuciones, y tipado `getAtomLabel` (regla `R`
   de sp2 conjugado + flags de valencia S/P). Desviación máxima medida:
   **0.0000 kcal/mol en 400 moléculas** (antes: mediana 143 kcal/mol en anillos
   de 3–4). Bug adicional encontrado en el camino: la convención del segundo
   normal de torsión desplazaba 180° el diedro y rompía las periodicidades
   impares (n=3).
10. **Tipado UFF (F17)** — predicado de conjugación nuevo
    (`conjugated_bond_flags`) validado contra `setConjugation` de RDKit:
    **0 discrepancias en 33 641 etiquetas** sobre 298 moléculas.
12. **Solucionador UHF/ROHF (F38)** — la comparación con PySCF reveló ciclos
    límite de DIIS y paradas prematuras (errores de hasta 0.25 Ha en NH₂/O₂ y
    no-convergencia en OH/CH₃). Correcciones en `src/scf/uhf.rs`:
    - criterio de convergencia con **gradiente de Pulay** `||FPS−SPF|| < 1e-7`
      medido **antes** de la extrapolación DIIS, más `ΔE < 1e-10`;
    - escalado automático de esquemas (config → punto fijo sin DIIS → DIIS con
      4096 iteraciones);
    - level-shift corregido: ya no entra en la energía ni en la Fock devuelta.
    Resultado: los 8 casos convergen y coinciden con PySCF (misma conjetura
    H_core) a **6.5e-8 Ha**.
13. **Integrales del módulo HF (`src/hf`)** — la comparación de F41 reveló que
    el motor Obara-Saika local de `hf::integrals`/`hf::nuclear` producía
    errores de hasta **6.7e-2 a.u.** en V y ERIs con funciones p sobre dos
    centros (CO, CH₂O, C₂H₄), y 1.9e-3–4.2e-2 Ha en las energías HF (CO caía
    0.20 Ha por debajo del mínimo RHF). Corregido **redirigiendo los
    integrales de `hf` al motor de `scf`** (validado a 1e-11 contra PySCF,
    F36/F37) mediante `hf::basis::to_scf_basis`; los residuos quedan en
    5e-8 a.u. (solo normalización en coma flotante). Esto también mejora
    HF-3c, que usa este módulo.
14. **Grupos espaciales (F64)** — la comparación contra pymatgen/spglib mostró
    que 172/230 grupos tenían el orden incompleto (faltaba el centrado de red)
    y 89 símbolos en notación corta no canónica. Se regeneró la tabla completa
    desde los datos ITC de pymatgen/spglib (`src/materials/space_group_table.rs`,
    230 grupos, 4425 operaciones, traslaciones en doceavos exactos) y el test
    compara orden, símbolo, sistema, red, grupo puntual **y el set completo de
    operaciones**: 230/230 exactos.
15. **Fukui/EAS (F47)** — ranking top-1 de f− (EHT) sobre 8 bencenos
    sustituidos + piridina, con clases de anillo desde RDKit y ground truth de
    reglas EAS confirmadas por DFT (RHF/6-31G* en los casos límite):
    **7/8 correctos** (orto/para para activantes, meta para benzaldehído y
    posición beta en piridina). El fallo (nitrobenceno: EHT predice para,
    DFT y experimento dan meta) es limitación genuina del método EHT,
    documentada en el fixture.
16. **XRD de polvo (F69)** — tres correcciones en `simulate_powder_xrd`
    (`src/materials/cell.rs`): longitud de onda CuKα 1.5406 → **1.54184 Å**
    (promedio ponderado, convención pymatgen); factor
    **Lorentz-polarización** (sin él, corr. 0.16); **factores de forma
    atómicos** dependientes del ángulo (tabla Waasmaier-Kirfel de 98 elementos
    en `src/materials/scattering_factors.rs`, generada del JSON oficial de
    pymatgen). Resultado: **Δ2θ = 0.0000°**, correlación de intensidades
    **≥ 0.9982** en NaCl/Si/TiO₂-rutilo/Cu.
17. **D4 standalone (F46): de gap a validación real** — depuración por
    etapas contra dftd4: (1) faltaba el factor de conteo EN (`en_factor`,
    fórmula exacta de mctc-lib/dftd4); (2) nuestro EEQ daba cargas 10×
    pequeñas → port exacto del modelo EEQ2019 (ecuaciones + tablas, radios
    D3, validado a 4.4e-7 e en cargas); (3) el camino standalone usaba
    tablas GFN2 en vez de DFT (`clsq`/`clsh`). Añadidos: `D4Damping`
    configurable, selector `D4RefSet::{Gfn2, Dft}` y recableado de
    `dispersion::compute_d4_energy` sobre `D4Model` + EEQ2019. Resultado:
    D4 **1.2e-2 → 1.2e-10 Ha** vs dftd4/PBE (18 moléculas); GFN2 mejora de
    paso **6.5e-4 → 2.3e-5 eV** (F45).
18. **PM3: atracción electrón-core y términos de intercambio NDDO (F42)**
    — la comparación contra ΔHf experimentales reveló que el Hamiltoniano
    core de PM3 omitía la atracción electrón-core −ΣZ_B·γ_AB (orbitales a
    +41 eV, electrones no ligados) y que la matriz G omitía 3 términos que
    su propio comentario enunciaba (self ½, intercambio diagonal ν≠μ,
    intercambio bicéntrico −½Pγ, Coulomb (sp|sp)). Corregidos en
    `src/pm3/{solver,gradients}.rs` con γ de Dewar-Ohno-Klopman y derivadas
    analíticas correspondientes (núcleo, 2c, e-core, intercambio, gaussianas).
    El SCF ahora converge en toda la serie (antes ciclos límite) y el
    gradiente analítico se acercó al numérico (3.98 → 1.21 eV/Å en H₂O),
    pero ΔHf sigue a ~10³ kcal/mol del experimento.
19. **D3 con interpolación CN (F46)** — `hf::d3` usaba C6 atómicos fijos
    (×2 de error). Port fiel del modelo D3: conteo CN logístico (k=16),
    pesos gaussianos (wf=4), tabla C6 7×7×pares (262k valores), r4r2 con
    `sqrt(0.5·R·√Z)` y ATM Axilrod-Teller (C9, rs9=4/3, alp=16, rvdw).
    Hallazgos: el fixture original estaba en Å (s-dftd3 espera Bohr);
    `static` en vez de `const` para la tabla 2.1MB (stack overflow en debug).
    Resultado: D3 **5.6e-2 → 2.8e-11 Ha** vs s-dftd3 (18 moléculas, con ATM).
20. **Maquinaria vibracional (F48)** — nueva API pública
    `ir::build_vibrational_analysis_from_hessian` (matriz externa en
    eV/Å²): alimentada con el Hessiano analítico de xtb-GFN2 reproduce sus
    frecuencias a **0.009 cm⁻¹** en 39 modos (6 moléculas). Las frecuencias
    con métodos propios quedan documentadas como limitación (EHT/PM3/GFN0
    heredan los gaps F42–F44; UFF solo cualitativo en flexiones).
21. **Asignación IR (F49)** — `group_is_possible` exigía un solo
    elemento por grupo (C-N sin C, S=O sin O, C=C con un solo C,
    N=O sin N): ahora refleja `group_matches_mode` + conteos
    (C=C exige ≥2 C). Ventanas corregidas contra NIST: O-H libre
    3650→3760 (agua ν1/ν3 3657/3756), C=O aldehído →1750
    (formaldehído gas 1746). Grupos nuevos: C≡O (CO 2143), S-H
    (H₂S 2614/2628), H-F (4138), H-H (4401), H-Cl (2886).
    Recall 10/10 + L0 del espectro (no-negatividad, integral=ΣI
    al 1.4 %, máximo sobre el modo más intenso).
22. **Optimizador BFGS (F70)** — roto: `s = p` en vez de
    `s = α·p` (corrompía la secante al actuar backtracking), rescate
    `alpha < 0.1` que ACEPTABA pasos cuesta arriba (colapso atómico,
    E ~ 10³), y convergencia falsa por `|ΔE|` entre pasos malos
    (etanol "convergido" en E=1248). Ahora: dirección de descenso
    verificada + reset a identidad, backtracking hasta 1e-12 sin
    aceptar subidas jamás (reintento SD; si falla, parada honesta),
    stall solo convergente si las fuerzas acompañan, historial SD
    solo con pasos aceptados. 10/10 moléculas al mismo mínimo que
    RDKit (|ΔE|/átomo 7.8e-10, RMSD 3.4e-6 Å).
23. **Gradientes UFF angle/torsión (F70)** — el optimizador los delató:
    angle con signo global invertido (faltaba el `-` de dθ/dr =
    −dcos/dr/sin; error 46 kcal/mol/Å vs diferencias finitas),
    torsión con forma Blondel-Karplus incompatible + regla de palanca
    en átomos centrales (0.4). Re-derivados por regla de la cadena
    (verificados término a término a ~1e-8, Hessiano vs energías con
    convergencia δ² exacta y simétrico a máquina). Bond/vdW ya
    estaban exactos; inversión verificada en el mismo barrido.
24. **Oráculo RDKit rancio (F70)** — `ff.CalcEnergy()` con el objeto FF
    pre-optimización devuelve la energía SIN relajar tras
    `UFFOptimizeMolecule` (ciclohexano: 37.24 vs 15.36 real);
    el oráculo re-evalúa con un FF fresco sobre el confórmero final.
25. **Topologías bcu/fcu (F67)** — `bcu()` era un clon de `dia`
    (4,4) en vez de la red 8-c (le faltaban 4 aristas); `fcu()`
    daba grados no uniformes (9,7,7,7) en vez de 12. Reescritas:
    bcu CsCl-like 8 aristas, fcu FCC 24 aristas (grado 12 uniforme).
    16/16 coordinaciones vs RCSR/pormake. Divergencia documentada:
    el `kag.cgd` de pormake es una variante apilada 6-c (arista según
    c), mientras RCSR kag y nuestro kgm son la capa kagome 4-c.
26. **PLD por percolación (F68)** — `pore_limiting_diameter` devolvía
    2×min-clearance ≈ 2×sonda siempre (sin análisis de percolación).
    Ahora: BFS con aristas PBC + bisección sobre la sonda (MOF-5:
    9.72 vs 10.01 Å Zeo++). API afinada: `porosity_grid_stats`
    (núcleo sin PLD para barridos de sonda) y bisección a 0.05 Å.
27. **CIF con loop mixto (F68)** — el parser abortaba la lectura de
    filas ante el primer tag no-`_atom_site_` (pormake añade
    `_atom_type_partial_charge`): 0/53 sitios. Ahora colecta todos
    los tags y solo exige coordenadas.
28. **CIS singlete (F41)** — faltaba el factor 2 del término de Coulomb
    `2(ia|jb) − (ij|ab)`, lo que bajaba las energías de excitación 1.2–4.9 eV.
    Corregido en `src/hf/cis.rs`: **6.9e-5 eV** de error máximo vs PySCF TDA
    en 29 estados (8 moléculas); la energía HF del módulo `hf` coincide a
    1.2e-8 Ha.
29. **GFN1-xTB fiel (F44)** — reescritura de `src/xtb/gfn1.rs` como port
    término a término de tblite 0.7.0, con parámetros generados desde
    `gfn1.f90` (`src/xtb/gfn1_params.rs`) y la dispersión D3 reutilizando el
    modelo `hf::d3` (ya validado). Corrigió: escalado ×0.1 de
    `p_hubbard_derivs`, `evtoau` histórico de GFN1, la convención de
    `slater_to_gauss(norm=.true.)` en la ortogonalización de Schmidt
    (1s/2s de H), y las tablas STO-6G de `sto_overlap.rs`, que caían a
    STO-4G para todo n<6 (C/N/O/F usaban 4 primitivas en vez de 6).
    Desviación: **581 eV → 8.1e-5 eV** (18 moléculas, HOMO/LUMO ~4e-5 eV).
    La dispersión D3(BJ) de GFN1 se validó por separado contra dftd3/py
    (1e-10 Ha).
30. **Aromaticidad Kekulé (F16)** — nueva `rings::perceive_aromaticity`:
    percepción Hückel 4n+2 sobre los anillos SSSR que reescribe los enlaces del
    anillo a aromáticos (orden 1.5) y promueve los átomos a sp², replicando la
    sanitización de RDKit. El builder UFF la aplica sobre una copia antes del
    tipado. El oráculo F16 dejó de excluir los casos de aromatización:
    **398 → 400 moléculas, max |ΔE| = 0.0000 kcal/mol**. Verificación directa
    contra RDKit en 1130 moléculas: la percepción solo toca los casos Kekulé y
    coincide con los enlaces aromáticos de RDKit.
11. **Tablas STO-3G y normalización (F36/F37/F23)** — la comparación contra
    PySCF encontró y corrigió:
    - typo en el 2.º primitivo de H (`0.62353014` → `0.62391373`);
    - datos de P/S/Cl **no canónicos** en `scf/basis.rs` y `hf/basis.rs`
      (exponentes y coeficientes de 3s/3p equivocados); regenerados desde
      PySCF en el orden de AO canónico `1s,2s,3s,2p,3p`;
    - coeficientes con precisión completa (antes 8 decimales, error ~1e-7 en
      las integrales);
    - **normalización de contracciones** en `scf::basis::from_shells`
      (convención PySCF, `S_ii = 1` exacto);
    Efecto medido: ERI `9.5e-5 → 2.1e-11`, `E_tot 2.3e-5 → 1.2e-8 Ha`,
    matrices `6.7e-4 → 3.1e-10`.

### 13.3 Hallazgos pendientes (documentados, no gated)

1. **Aromaticidad de SMILES Kekulé** — **CERRADO (percepción)**: la nueva
   `rings::perceive_aromaticity` aplica Hückel 4n+2 sobre los anillos SSSR y
   reescribe los enlaces del anillo a aromáticos (orden 1.5) antes del force
   field, replicando la sanitización de RDKit para C/N/O/P/S (pyrrole/thiophene
   donors, carbonilos exocíclicos, carbaniones/carbocationes). El corpus F16 se
   regeneró sin exclusiones: **400/400 moléculas, max |ΔE| = 0.0000 kcal/mol**,
   incluyendo `N=C1OC=CO1` y `N=C1NC=CO1`. Validación directa contra RDKit:
   sobre el corpus (1130 moléculas, fragmentos múltiples excluidos) la
   percepción solo modifica los casos Kekulé y coincide al 100 % con los
   enlaces aromáticos de RDKit.
   Límite restante (no gated): la sanitización de RDKit también puede
   **quitar** aromaticidad (`O=c1[nH]c(=O)c(=O)[nH]c1=O`), algo que sci-form no
   hace (solo añade marcas); el oráculo F16 lo registra como
   `skipped_aromaticity_reperception` (1 caso). SMARTS/ECFP siguen leyendo los
   enlaces tal como los deja el parser (sin percepción); el fixture F01
   compara topología **sin sanitizar** por diseño.
2. **Cobertura pendiente del catálogo** — las ~80 fichas restantes siguen
   el formato ya operativo: `scripts/validation/oracle_<tool>.py` → fixture →
   test en `tests/cross_validation.rs` (o suite por dominio) → `collect.py`.
   Candidatas inmediatas: F40 (HF-3c vs ORCA), F50/F52 (ORCA),
   F78/F95 (LAMMPS), F82/F89 (GSM/pysisyphus), F94 (geomstats/pymanopt),
   F34/F88/F91 (pybinding).
3. **PM3/NDDO (F42): DESBLOQUEADO (2026-10-12)** — el binario de MOPAC 23.2.5
   (release oficial Linux, sin toolchain Fortran) se descargó y usó como
   árbitro. Se portó el núcleo fiel (`src/pm3/faithful.rs`): tabla de
   parámetros de la fuente (42 elementos + `eisol`), overlaps analíticos
   `diat2` (casos 1–6), `reppd`/`rotmat`/`spcore`/`ccrep` (con la corrección
   N–H/O–H), `fock1`/`fock2` y HoF con `eisol`. Resultado: hidruros X–H
   exactos (<0.02 kcal/mol), matriz H de agua/etileno exacta (1e-4 eV),
   E_nuc por par exacto (1e-4 eV) y SCF de agua con eigenvalores y cargas
   exactos. Residual: las integrales 2e "polares" (π) de `tx`/`rotatd`
   salen con signo opuesto (≤130 kcal/mol en sistemas π); el fixture
   `mopac_pm3.json` (24 moléculas) queda para iterar. Hallazgos de la
   forense previa (ya resueltos por el port):
   - la γ de Dewar necesita los parámetros de apantallamiento exactos
     (nuestra ρ desde Gss es una aproximación ±20 %);
   - las energías atómicas aisladas usan una fórmula ad-hoc (el /3 del
     término p-p es sospechoso; MOPAC usa Eisol tabulado, ausente aquí);
   - el punto fijo SCF es inestable: H₂ rompe simetría (±1.28 e) y los
     resultados varían entre builds (se amplifica ruido fp); DIIS colapsa a
     estados no físicos; damping solo converge lento (H₂ 241 iters);
   - densidades no físicas (H₂O con polaridad invertida, HF con q=0.000);
   - `test_pm3_h2` y `test_pm3_gradient_vs_numerical` quedan en `#[ignore]`
     hasta la revalidación (criterios estrictos preservados, no relajados).
   El oráculo `scripts/validation/oracle_pm3.py` trae rama MOPAC lista
   (`MOPAC=/ruta/mopac`): en cuanto exista toolchain, el fixture se convierte
   en validación cruzada real.
4. **GFN0-xTB (F43)** — cerrado con `xtb::gfn0::solve_gfn0`: port
   término a término del xtb 6.7.1 (EHT/CEH, EEQ, repulsión efectiva,
   SRB, D4-GFN0). Hallazgos por el camino:
   - CN erf con k=7.5 (no el defecto mctc 3.15): probado a 7 decimales
     en la energía EEQ;
   - H-2s de Schmidt contra 1s + renormalización (5 primitivas, no 2G);
   - bloques-d cartesianos con pesos trafo + matriz dtrf2 exacta
     (sin esto, S indefinida);
   - normalización común-xx en d (sin esto, solapes d un factor √3
     de más y S indefinida);
   - `pre` del SRB usa EN de Pauling mientras `rab0` usa la tabla
     fitted de approxrab (6 decimales);
   - el oráculo relee la energía con un FF fresco (`calcEnergy`
     post-`UFFOptimizeMolecule` va rancio; esto era F70).
   La aproximación vieja (`solve_xtb`) queda **solo** como proveedor del
   gradiente analítico de los backends de fuerza `"xtb"` (NEB/MD en
   `dynamics.rs`, `compute_xtb_gradient` público): los consumidores de
   energía (Hessiano numérico → GFN2, DOS/mesh → GFN2, dipolo → GFN2,
   backends GSM → `compute_xtb` = GFN2) ya usan los solvers validados
   (2026-10-13).
5. **GFN1-xTB (F44)** — **CERRADO**: `xtb::gfn1::solve_gfn1` es un port
   término a término del calculador nativo de tblite 0.7.0
   (`h0.f90`, `coulomb/charge/{effective,type}.f90`,
   `coulomb/thirdorder.f90`, `repulsion/effective.f90`, `disp/d3.f90`,
   `classical/halogen.f90`, `basis/ortho.f90`, `scf/{iterator,potential}.f90`):
   H0 Slater con kpair/kshell/enscale, electrostática shell-resolved con
   media armónica de Hubbard (gexp=2) + tercer orden atómico, repulsión
   efectiva, D3(BJ) con parámetros GFN1 y corrección de halógeno. Medido:
   **8.1e-5 eV** de energía máxima y ~4e-5 eV en HOMO/LUMO vs tblite en 18
   moléculas (F44 en `tests/cross_validation.rs`).
   Hallazgos por el camino:
   - la tabla `p_hubbard_derivs` va escalada ×0.1 en Fortran (el parse la
     dejaba ×10);
   - `evtoau` de GFN1 es 1/27.21138505 (conversión antigua "por consistencia"),
     no la constante física;
   - `slater_to_gauss(..., norm=.true.)` multiplica por
     `(2α/π)^{3/4}(4α)^{l/2}/√((2l−1)!!)` y `orthogonalize` opera en esa
     convención con kernel desnudo `(π/(a+b))^{3/2}` — hay que replicarla
     exactamente (el 2s de H es el único caso con l repetido);
   - las tablas STO-6G de `sto_overlap.rs` solo tenían n=6 y caían a STO-4G
     para todo lo demás: C/N/O/F usaban 4G en vez de 6G (error ~1e-4 a.u.
     en S y ~0.01 eV en energías);
   - el CN del Hamiltoniano es el mismo `cn_count%exp` (k=16, radios D3) que
     usa la dispersión D3 de tblite.
   Limitación conocida: el SCC es restringido (cerrado); los sistemas con
   electrones desapareados usan reparto promediado y no están gated (igual
   que GFN2).
6. **D3/D4 standalone (F46)** — los módulos `hf::d3` (D3-BJ con parámetros
   HF-3c) y `dispersion` (D4) usan C6 por elemento sin dependencia de
   coordinación/carga: D3 sobreestima ×2 (máx 5.6e-2 Ha en tolueno) y D4
   subestima ×4–6 (máx 1.2e-2 Ha). El camino D4 interno de GFN2
   (`xtb::d4`, refalpha/scaled C6 de tblite) sí está validado de forma
   indirecta al coincidir la energía GFN2 completa con tblite.
   **Revisión de gaps (2026-10-11)**: medición directa `D4Model` vs dftd4 con
   mismos parámetros de damping y solo 2-cuerpos sigue dando un factor
   2.4–3.6× (H₂O/CH₄), luego el port `xtb::d4` tiene un bug a nivel de
   ecuaciones/pesos C6 (no es solo parametrización). Descartados como causa:
   números de coordinación (H₂O: O≈1.99, H≈1.0 ✓), unidades (Bohr ✓ en ambos
   lados) y cargas (nuestro EEQ da q≈0 y apenas mueve el resultado).
   Sospechosos restantes: escalado zeta/carga (falta `alp`), integración
   Casimir-Polder y ensamblado de `dispmat`/energía. El error D4 absoluto es
   submHa y por eso cabe dentro de la tolerancia de F45, pero sesga GFN2.
   **Actualización: D4 y D3 CERRADOS.** D4: el bug estaba en el conteo
   CN (faltaba `en_factor`), las cargas (EEQ 10× pequeño → port EEQ2019
   exacto) y las tablas (GFN2 en vez de DFT): D4/PBE a 1.2e-10 Ha y GFN2 a
   2.3e-5 eV. D3: C6 con interpolación CN + ATM, a 2.8e-11 Ha vs s-dftd3.
   **GFN0/GFN1**: GFN0 ya era un port fiel de xtb (1.7e-4 eV, F43) y GFN1 se
   cerró como port fiel de tblite (8.1e-5 eV, F44), así que la decisión
   "no portar tblite" quedó obsoleta: el núcleo TB propio iguala a tblite a
   2.3e-5–8.1e-5 eV en 18 moléculas.
7. **Precisión de la función de Boys** — el residuo de F36 (3.1e-10 a.u.,
   ~5e-12 relativo) proviene del cambio serie/asintótica en `x ≈ 30 + n`
   (`src/scf/gaussian_integrals.rs::boys_function`). Mejorable con una rama
   asintótica más temprana o recursión hacia abajo; no bloquea el gate.
8. **`hf` module**: sus tablas de P/S/Cl se corrigieron y sus integrales
   ahora usan el motor de `scf`, pero Br/I usan un
   esquema reducido (sin funciones d) y el camino HF-3c aún no tiene fixture
   externo (F40 pendiente vs ORCA/CCCBDB).

### 13.4 Infraestructura entregada

| Artefacto | Ruta |
|---|---|
| Helper común (reportes, Kabsch, corpus, unidades) | `scripts/validation/common.py` |
| Oráculo RDKit multi-ficha | `scripts/validation/oracle_rdkit.py` |
| Colector de métricas → summary | `scripts/validation/collect.py` |
| Runner end-to-end | `scripts/run_validation.sh` |
| Fixtures + procedencia | `tests/fixtures/validation/` |
| Versiones de referencia | `validation/refs.lock.json` |
| Tests de validación cruzada | `tests/cross_validation.rs` (28 tests activos: F01, F05, F13, F14, F16, F17, F21, F23+F26, F28, F36–F39, F41, F43–F49, F59, F64–F70) |
| Reporte | `reports/validation/summary.json` |

---

## 14. Track G — banco de validación, umbrales y escalado

### 14.0 Gates nuevos (2026-10-12, sesión "próximos 10")

- **Strong coupling (NMR)**: `nmr::second_order::SpinSystem` diagonaliza el
  Hamiltoniano de espín (base producto, n ≤ 6) y da el espectro exacto de
  líneas; gate contra la solución analítica AB (posiciones y intensidades
  `1 ± sin2θ`) a 1e-8, límite débil (dobletes 1:1) y A2 (singlete).
  `compute_nmr_spectrum` anota los pares con |Δν|/|J| < 10.
- **Beta**: `tests/beta_benchmarks.rs` (feature-gated) compara KPM vs
  diagonalización (DOS integrada 40.17/40, centro 1e-4, corr 0.86 con kernel
  Jackson), Nyström vs SVD exacta (error ≤ 2× el mejor rank-k) y MBH vs el
  Hessiano analítico proyectado (desviación máxima 2.6 %).
- **D4/EEQ**: tablas regeneradas desde dftd4 4.3.0 hasta Z=118 y EEQ2019 hasta
  Z=103 (multicharge); F46 con 7 casos pesados exactos a 6e-10 Ha.
- **Variable-cell**: stress numérico (σ = (1/V)∂E/∂ε) + relajación
  posiciones/celda; gate autocontenido con un modelo LJ (recupera a_eq y
  stress ≈ 0).
- **Topologías**: 6 nets RCSR importados (aristas no equivalentes en
  afw/xah/tfy) y gate F67 contra pormake.
- **EEQ-BC 2025 (2026-10-13)**: port término a término del modelo bond-capacitor
  (`charges_eeq::compute_eeqbc_charges`) validado contra el CLI de multicharge
  (oráculo directo, `--model eeqbc`): F22 con 23 moléculas (incluye PbI₂, HgCl₂,
  Au₂, SnF₄, XeF₂, WO₃, UF₆, AmCl₃, ferroceno) — cargas a 6.9e-12 e y CN a
  1.9e-11 (tolerancia 1e-10). Hallazgo de fidelidad: la distancia vdW que
  entra al capacitor es el valor crudo de la tabla mctc-lib (Å) porque
  `aatoau`/`autoaa` se cancelan.
- **EHT metales (2026-10-13)**: benchmark de la tabla d-block contra YAeHMOP
  3.0.3 (Alvarez 1993); 33 valores corregidos (ζ de 4d/5d, VSIP de Ni/Zr/Tc/Cd)
  y gate F31 con Δ = 0 en 76 filas; Y/Ag/Hf documentados como referencia
  incompleta.
- **CISD (2026-10-13)**: `hf::cisd` reescrito como CI determinantal completo
  (dobles de mismo orbital + acoplamiento con la referencia + energía de
  correlación real). Gate F41-CISD contra PySCF `ci.cisd`/`ci.ucisd`: 8
  moléculas STO-3G a ≤ 2.1e-8 Ha; matriz verificada contra una referencia de
  segunda cuantización explícita (`tests/cisd_brute.rs`, 0 discrepancias por
  encima de 1e-9 en un sistema sintético de 361 determinantes). Hallazgo: el
  signo de Slater–Condon exige aplicar las creaciones en orden descendente
  (`a†_a a†_b a_j a_i`), no ascendente; el error era invisible para el estado
  fundamental de H2O pero corrompía CH2O/C2H4 (8.5e-3 Ha).
- **Destino de `solve_xtb` legacy (2026-10-13)**: los consumidores de energía
  se migraron a los solvers validados — `HessianMethod::Xtb` y el dipolo
  vibracional → GFN2 (`solve_gfn2`, coherente con el `compute_xtb` público y
  con los Hessianos de referencia F48 "xtb 6.7.1 GFN2-xTB"), `DosMethod::Xtb`
  y el mesh orbital → GFN2. El solver legacy se conserva **solo** como
  proveedor del gradiente analítico de los backends de fuerza `"xtb"`
  (NEB/MD + `compute_xtb_gradient` público); su eliminación queda pendiente de
  gradientes GFN0/GFN2 analíticos.



> Medido el 2026-10-12 en la misma máquina de referencia (Intel i5-10500H,
> 12 hilos lógicos), build `--release`.

### 14.1 Banco de moléculas (≥100)

| Corpus | Moléculas | Uso |
|---|---|---|
| `tests/fixtures/diverse_molecules.smi` | 180 | F01/F05/F16/F17/F21/F59 (parse, SSSR, UFF, Gasteiger, ECFP) |
| `data/chembl_1k_practical.smi` | 1000 | muestra práctica de ChEMBL para los oráculos |
| `tests/fixtures/reference_coords_500.json` | 400 | F16/F70 (geometrías + energía UFF/RDKit) |
| `tests/fixtures/validation/tblite_gfn.json` | 18 | F43–F45 (GFN0/GFN1/GFN2) |
| `tests/fixtures/validation/pyscf_*.json` | 11–13 | F36–F39, F41 (SCF/ERI/TDA) |

El guard del banco vive en `tests/cross_validation.rs::track_g_bank_has_enough_molecules`
(≥100 SMILES válidos y parseables).

### 14.2 Umbral < 1 % y valores medidos

El umbral del 1 % se aplica a los observables *continuos* con referencia
externa; los gated van muy por debajo:

| Observable | Referencia | Medido | Umbral |
|---|---|---|---|
| Energía UFF | RDKit | 0.0000 kcal/mol (400) | 0.01 |
| Energía GFN1-xTB | tblite | 8.1e-5 eV (18) | 0.01 eV |
| Energía GFN2-xTB | tblite | 2.3e-5 eV (18) | 0.0217 eV |
| Energía GFN0-xTB | xtb 6.7.1 | 1.7e-4 eV (18) | 0.01 eV |
| E_tot RHF/STO-3G | PySCF | 1.2e-8 Ha (11) | 1e-7 Ha |
| ERI AO | PySCF | 2.1e-11 a.u. (13 074) | 1e-10 |
| ERI AO→MO | PySCF | 1.9e-14 a.u. | 1e-10 |
| Gaps EHT metales | fixture propio | <0.1 % o tolerancia absoluta | — |
| ΔHf PM3 | NIST | gap documentado (F42) | — |

### 14.3 Límites de escalado (single point, alkanos lineales)

| N (átomos) | EHT (ms) | GFN2 (ms) | PM3 (ms) |
|---|---|---|---|
| 26 | 2.0 | 8.3 | 241 (no conv.) |
| 50 | 8.1 | 31.8 | 1254 (no conv.) |
| 74 | 19.0 | 81.5 | 3719 (no conv.) |
| 98 | 37.2 | 173.1 | 8412 (no conv.) |

- EHT y GFN2 escalan ≈ O(N²) en este rango (el cuello es la diagonalización
  O(N³) con bases pequeñas); PM3 ≈ O(N²·⁷) y **no converge** en la serie
  (gap F42, bloqueado por MOPAC).
- Memoria: los caminos HF/SCF guardan ERIs O(N⁴) en memoria
  (`(N(N+1)/2)² × 8 B`); con STO-3G, N=150 funciones ≈ 1.0 GB y N=200 ≈ 3.2 GB.
  Por encima de ~200 funciones hace falta batching (o el camino GPU ERI,
  pendiente) antes de escalar.
- El pipeline de conformeros (ETKDG+UFF) es el que sostiene el throughput de
  producción (batch paralelo ~1–2 ms/mol en moléculas pequeñas; ver TESTING.md).

---

*Documento generado para el repo `sci-form 0.15.2`. Mantener sincronizado con `TESTING.md`
(comandos y suites existentes), `docs/algorithm/quantum-validation-roadmap.md` (estado de los
métodos semi-empíricos) y `docs/algorithm/validation.md` (checks geométricos internos).*

