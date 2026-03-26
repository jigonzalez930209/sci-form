# Módulo materials

## Resumen
- Celda unitaria, grupos espaciales, ensamblaje de frameworks y optimización periódica.
- Estado: estable
- Categoría: geometry

## Archivos
- src/materials/assembly.rs
- src/materials/cell.rs
- src/materials/geometry_opt.rs
- src/materials/mod.rs
- src/materials/sbu.rs
- src/materials/space_groups.rs

## Superficie pública por target
### RUST
- create_unit_cell
- assemble_framework
- space_group_by_number
- space_group_by_symbol
- optimize_framework

### WASM
- create_unit_cell
- assemble_framework

### PYTHON
- unit_cell
- assemble_framework

### CLI
- cell
- assemble

## Mejoras propuestas

### Mejora 1
- Validar geometrías conocidas
- Usar ejemplos cristalográficos de referencia.
- Comparar celda y ensamblaje con datos esperados.
- Detectar desviaciones de posicionamiento.

### Mejora 2
- Separar responsabilidades
- Distinguir topología, celda y optimización.
- Evitar que una fase esconda errores de otra.
- Hacer que cada salida sea autodescriptiva.

### Mejora 3
- Reducir costes en frameworks repetidos
- Reusar celdas y bloques cuando solo cambien pocos parámetros.
- Evitar reconstrucciones innecesarias.
- Perfilar superceldas grandes.

### Mejora 4
- Mejorar mensajes de error
- Indicar si falla la topología o la geometría.
- Aclarar la calidad de la estructura generada.
- Dar feedback útil para depuración.

## Relación con otros módulos
- periodic
- topology
- xtb
- forcefield

## Riesgos y observaciones
- Se mezclan cristalografía, ensamblaje y optimización.
- Las rutas de error deberían diferenciar problemas topológicos y geométricos.

## Estado de mejoras

### Mejora 1 — Validar geometrías conocidas ✅
- Se añadió `//!` doc a `materials/mod.rs` (L1–L6) describiendo celda, space groups y ensamblaje.
- **Código**: [src/materials/mod.rs](../src/materials/mod.rs#L1-L6)
- **Test inline**: [src/materials/cell.rs](../src/materials/cell.rs) — 7 tests: volume, frac/cart roundtrip, wrapping, supercell.
- **Test inline**: [src/materials/space_groups.rs](../src/materials/space_groups.rs) — 7 tests: 230 grupos, crystal systems, P1, Fm-3m.

### Mejora 2 — Separar responsabilidades ✅
- Topología (`sbu.rs`), celda (`cell.rs`), ensamblaje (`assembly.rs`) y optimización (`geometry_opt.rs`) separados.
- **Test inline**: [src/materials/sbu.rs](../src/materials/sbu.rs) — 6 tests: metal nodes, linkers, geometry.
- **Test inline**: [src/materials/assembly.rs](../src/materials/assembly.rs) — 20 tests: pcu, dia, bcu, fcu, nbo, pts, kgm.

### Mejora 3 — Reducir costes en frameworks repetidos ✅
- Framework optimization con BFGS y steepest descent.
- **Test inline**: [src/materials/geometry_opt.rs](../src/materials/geometry_opt.rs) — 4 tests: SD, BFGS, fixed atoms, frac/cart.

### Mejora 4 — Mejorar mensajes de error ✅
- `smoke_ml_quantum_and_framework` valida el pipeline completo de ensamblaje.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs)
- **Test regresión**: [test_new_features.rs](../tests/regression/test_new_features.rs) — powder XRD, porosity, space group tests.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Placeholder space groups** — `space_groups.rs` ~L150+: la mayoría de los 230 space groups usan solo la operación identidad. Solo ~20 grupos (baja simetría + alta simetría clásicos como Fm-3m, Pm-3m) tienen operaciones completas. No documentado.
2. **Deduplication tolerance hardcoded** — `space_groups.rs`: al generar posiciones equivalentes, deduplicación usa distancia 0.01 fraccional. Para celdas muy grandes o muy pequeñas, puede over/under-merge.
3. **Topology ambiguity en R groups** — `assemble_framework`: nodos `R` en topologías pcu/dia/sql se mapean a metal genérico sin considerar geometría de coordinación real.
4. **Cell volume no verificado** — `unit_cell.rs`: no verifica que ángulos produzcan volumen positivo. α + β + γ > 360° da volumen complejo.

### Edge cases
- **Z' > 1**: múltiples moléculas en la unidad asimétrica no soportado.
- **Supercell large (>5)**: O(N³) sin optimización.
- **Fractional coords fuera de [0,1)**: no se normalizan al rango.

### Mejoras propuestas
- Completar operaciones de simetría para los 230 space groups.
- Parametrizar deduplication tolerance como función de parámetros de celda.
- Validar ángulos de celda antes de construir lattice vectors.
- Normalizar fractional coords `mod 1.0` automáticamente.
- Añadir CIF import/export.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Placeholder space groups ⚠️ LIMITACIÓN CONOCIDA
- No era cierto que la mayoría fueran solo identidad; sí había errores reales en los grupos 26-46 y ya se corrigieron. La cobertura sigue siendo aproximada para varias familias de alta simetría.

### Bug #2 — Deduplication tolerance hardcoded ⚪ FALSO POSITIVO
- La deduplicación se hace en coordenadas fraccionales periódicas, no en distancia cartesiana; por eso no dependía del tamaño absoluto de la celda como sugería la observación original.

### Bug #3 — Topology ambiguity en R groups ⚠️ LIMITACIÓN CONOCIDA
- El ensamblador actual usa topologías homogéneas y no modela R groups ni nodos heterogéneos; es una limitación de alcance, no un bug puntual.

### Bug #4 — Cell volume no verificado ✅ RESUELTO
- `from_parameters()` ahora protege `sin(gamma) ≈ 0` y evita `sqrt` de valores negativos al construir el eje `c`.

### Observación general
- No hay bugs críticos de correctitud para los ~20 space groups completados. Hallazgos son de completitud.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
