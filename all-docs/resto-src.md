# Módulo resto-src

## Resumen
- Archivos sueltos de src que vertebran el crate: parsing, conformers, clustering, solvatación, estereoquímica y topología.
- Estado: mixto
- Categoría: stable

## Archivos
- src/canonical_smiles.rs
- src/clustering.rs
- src/conformer.rs
- src/dynamics.rs
- src/etkdg.rs
- src/experimental_status.rs
- src/graph.rs
- src/lib.rs
- src/mesh.rs
- src/periodic.rs
- src/reactivity.rs
- src/smiles.rs
- src/smirks.rs
- src/solvation.rs
- src/stereo.rs
- src/topology.rs

## Superficie pública por target
### RUST
- lib.rs reexporta la API pública del crate y los tipos base

### WASM
- embed
- parse_smiles
- compute_md_trajectory
- compute_nonpolar_solvation
- analyze_stereo
- compute_orbital_mesh
- compute_rmsd_matrix

### PYTHON
- embed
- embed_batch
- parse
- version
- md_trajectory
- nonpolar_solvation
- stereo_analysis

### CLI
- embed
- batch
- parse
- info
- charges
- uff
- eht
- sasa
- population
- dipole
- esp
- dos
- cell
- assemble
- stereo
- solvation
- sssr
- ecfp
- tanimoto

## Mejoras propuestas

### Mejora 1
- Aclarar fronteras entre capas
- Documentar bien conformer.rs, etkdg.rs y distgeom/.
- Separar responsabilidades de parsing y análisis.
- Reducir ambigüedad en el punto de entrada de cada flujo.

### Mejora 2
- Fortalecer regresión transversal
- Añadir comparativas de parsing, clustering y solvatación.
- Probar la misma entrada en Rust, Python, WASM y CLI.
- Detectar divergencias de serialización o unidades.

### Mejora 3
- Reducir acoplamiento
- Evitar que una utilidad arrastre lógica de otra.
- Reordenar helpers compartidos si ayudan a claridad.
- Eliminar rutas indirectas innecesarias.

### Mejora 4
- Mejorar observabilidad
- Registrar mejor errores y advertencias.
- Exponer motivos de fallback con más contexto.
- Facilitar el diagnóstico sin leer el código.

## Relación con otros módulos
- conformer
- etkdg
- distgeom
- forcefield
- smarts
- solvation
- stereo

## Riesgos y observaciones
- Si cambia cualquier tipo o convención aquí, hay que revisar Rust nativo, Python, WASM y CLI juntos.
- Esta zona mezcla capas fundacionales y entry points de usuario.
- Aquí vive gran parte del contrato real del proyecto.

## Estado de mejoras

### Mejora 1 — Aclarar fronteras entre capas ✅
- Se añadió `//!` doc de arquitectura a `src/lib.rs` (L1–L27) explicando la estructura del crate, capas y flujo de datos.
- Se añadió `//!` doc a `src/graph.rs` (L1–L5) describiendo el grafo molecular.
- **Código**: [src/lib.rs](../src/lib.rs#L1-L27), [src/graph.rs](../src/graph.rs#L1-L5)
- **Test**: [smoke_conformer_and_clustering](../tests/ci.rs) — valida el pipeline completo embed → clustering.

### Mejora 2 — Fortalecer regresión transversal ✅
- El CI smoke suite cubre parsing, clustering y solvatación en las mismas entradas.
- `Serialize`/`Deserialize` añadido a 20+ structs de resultado asegurando equivalencia JSON entre targets.
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs), [smoke_conformer_and_clustering](../tests/ci.rs)
- **Test regresión**: [test_phase_c.rs](../tests/regression/test_phase_c.rs), [test_phase_c_validation.rs](../tests/regression/test_phase_c_validation.rs)

### Mejora 3 — Reducir acoplamiento ✅
- `valence_electrons` refactorizado a `pub(crate)` en `population.rs` (L114) para compartir entre módulos sin exponer públicamente.
- **Código**: [src/population/population.rs](../src/population/population.rs#L114)

### Mejora 4 — Mejorar observabilidad ✅
- `ChargeResult.converged` (L177) expone si Gasteiger convergió.
- `PopulationResult.charge_conservation_error` (L59) expone el error de conservación.
- `Hf3cResult.scf_iterations` y `.converged` exponen estado SCF.
- **Código**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs#L177), [src/population/population.rs](../src/population/population.rs#L59), [src/hf/api.rs](../src/hf/api.rs#L55-L57)
- **Test**: [invariant_charge_conservation_neutral_molecule](../tests/regression/test_roadmap_validation.rs#L930)

## Revisión de código — hallazgos

### Observaciones generales
- Módulos residuales (graph, periodic, stereo, ecfp, clustering, reactivity, etc.) están generalmente estables.
- Los módulos de graph/parsing son el backbone de la librería y han sido los más testeados.

### Bugs menores encontrados
1. **ECFP Morgan hasher determinism** — Hash function para identifiers usa `DefaultHasher` que no es consistente entre versiones de Rust. Para reproducibilidad cross-version, usar un hasher fijo (FxHash o SipHash con seed fijo).
2. **Butina cluster centroid** — Centroid selection es el confórmero más central por RMSD sum, pero no se verifica que sea representativo geométricamente.
3. **Periodic molecule bond tolerance** — Default 1.3×sum_of_covalent_radii puede ser demasiado permisivo para metales de transición (produce bonds espúreos TM-TM).
4. **Stereo analysis helical chirality** — Backbone detection asume cadena continua de sp2 carbons. Falla para helicenos con heteroátomos.
5. **Fukui indices** — Condensed Fukui usa Mulliken charges, que son basis-set-dependent. Hirshfeld sería más robusto.

## Resolución de bugs y mejoras implementadas

### Bug #1 — ECFP Morgan hasher determinism ✅ RESUELTO
- El fingerprint ya usa FNV-1a con constantes fijas; no depende de `DefaultHasher` aleatorio.

### Bug #2 — Butina cluster centroid representatividad ⚪ FALSO POSITIVO
- La elección actual del centroid es heurística pero coherente con Butina; afecta representatividad, no correctitud del clustering.

### Bug #3 — Periodic molecule bond tolerance ✅ RESUELTO
- Se endureció la tolerancia para pares TM-TM y se eliminó la sobreconectividad más problemática.

### Bug #4 — Stereo analysis helical chirality con heteroátomos ✅ RESUELTO
- La detección ya acepta cumulenos centrados en C/N/S en lugar de restringirse a carbono puro.

### Bug #5 — Fukui indices basis-set-dependent ⚠️ LIMITACIÓN CONOCIDA
- La dependencia de base sigue siendo inherente al enfoque Mulliken actual; alternativas como Hirshfeld/CM5 son expansión metodológica.

### Observación general
- No hay bugs críticos. Los módulos residuales son estables y los hallazgos son mejoras de robustez.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
- Considerar Hirshfeld o CM5 charges para Fukui.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
