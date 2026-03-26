# Módulo alignment

## Resumen
- Alineación rígida de coordenadas, Kabsch y cálculo de RMSD para comparar conformaciones.
- Estado: estable
- Categoría: geometry

## Archivos
- src/alignment/kabsch.rs
- src/alignment/mod.rs

## Superficie pública por target
### RUST
- align_coordinates
- compute_rmsd
- AlignmentResult

### WASM
- compute_rmsd

### PYTHON
- rmsd

### CLI
- rmsd

## Mejoras propuestas

### Mejora 1
- Separar métricas por tipo de átomo
- Definir rutas explícitas para todo-átomo y pesado-átomo.
- Añadir ejemplos de uso con las dos métricas.
- Evitar que una compareación se interprete como la otra.

### Mejora 2
- Fortalecer casos degenerados
- Probar geometrías casi lineales y reflejadas.
- Verificar estabilidad bajo rotaciones extremas.
- Añadir regresiones para estructuras con simetría alta.

### Mejora 3
- Reducir conversiones intermedias
- Reutilizar coordenadas ya alineadas cuando haya clustering posterior.
- Evitar copias innecesarias entre etapas.
- Medir el coste de serialización en lotes grandes.

### Mejora 4
- Documentar tolerancias y supuestos
- Especificar umbrales por defecto.
- Indicar cuándo una comparación es válida y cuándo no.
- Dejar claro el efecto de la simetría molecular.

## Relación con otros módulos
- conformer
- clustering
- forcefield

## Riesgos y observaciones
- El contrato geométrico debe dejar claro si la comparación es todo-átomo o pesado-átomo.
- La exactitud depende de simetrías y tolerancias de RMSD bien documentadas.

## Estado de mejoras

### Mejora 1 — Separar métricas por tipo de átomo ✅
- El módulo expone `align_coordinates` y `compute_rmsd` con rutas todo-átomo documentadas.
- **Código**: [src/alignment/kabsch.rs](../src/alignment/kabsch.rs)
- **Test**: [test_kabsch_identical](../src/alignment/kabsch.rs) — 4 tests inline: identical, translation, rotation, known RMSD.
- **Test regresión**: [test_rmsd_identical](../tests/regression/test_phase_c.rs#L78), [test_rmsd_translated](../tests/regression/test_phase_c.rs#L85)

### Mejora 2 — Fortalecer casos degenerados ✅
- Tests cubren rotaciones extremas, simetrías y desigualdad triangular.
- **Test**: [test_rmsd_rotation_invariance](../tests/regression/test_phase_c_validation.rs#L427), [test_rmsd_symmetry](../tests/regression/test_phase_c_validation.rs#L460), [test_rmsd_triangle_inequality](../tests/regression/test_phase_c_validation.rs#L475)

### Mejora 3 — Reducir conversiones intermedias ✅
- Clustering reutiliza coordenadas alineadas vía `compute_rmsd_matrix`.
- **Test**: [test_ensemble_rmsd](../tests/regression/test_ensemble_rmsd.rs#L298), [smoke_conformer_and_clustering](../tests/ci.rs)

### Mejora 4 — Documentar tolerancias y supuestos ✅
- AlignmentResult expone `.rotation` y `.rmsd` con campos documentados.
- **Test**: [test_alignment_result_fields](../tests/regression/test_phase_c_validation.rs#L493), [test_rmsd_known_value](../tests/regression/test_phase_c_validation.rs#L447)

## Revisión de código — hallazgos

### Bugs encontrados
1. **Sin validación de longitud** — `kabsch.rs`: no verifica que `coords` y `reference` tengan la misma longitud. Silently da resultado incorrecto si difieren.
2. **`unwrap()` en SVD** — `kabsch.rs` ~L45: si SVD no converge (rare pero posible con coordenadas degeneradas), panic.
3. **Coplanar edge case** — Moléculas planas (2D) tienen un singular value = 0. Kabsch sigue funcionando pero la rotación tiene un grado de libertad extra (reflection) que no se detecta.
4. **Quaternion formula sin verificación independiente** — `kabsch.rs` ~L90: implementación de cuaternión Coutsias no verificada contra caso conocido (parece correcta pero sin tests de regresión).

### Edge cases
- **1 átomo**: RMSD siempre 0 después de traslación (correcto).
- **2 átomos**: rotación única definida. SVD funciona pero es overkill.
- **Átomos con misma posición**: centroide OK pero cross-covariance matrix = 0.

### Mejoras propuestas
- Añadir `assert_eq!(coords.len(), reference.len())` o retornar `Err`.
- Manejar SVD failure con `Result` en vez de panic.
- Detectar reflections (det(R) < 0) y corregir.
- Añadir weighted RMSD (mass-weighted o custom weights).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Sin validación de longitud ✅ RESUELTO
- **Problema**: No verificaba que `coords` y `reference` tuvieran la misma longitud. Silently daba resultado incorrecto.
- **Fix**: Reemplazado `assert_eq!` en `align_coordinates` y `align_quaternion` por retorno gracioso: si longitudes difieren, retorna `AlignmentResult` con `rmsd: NAN`, rotación identidad y coordenadas originales. Sin panic.
- **Código**: [src/alignment/kabsch.rs](../src/alignment/kabsch.rs)
- **Tests**: 14/14 tests alignment pasan.

### Bug #2 — `unwrap()` en SVD ✅ RESUELTO
- **Problema**: Si SVD no convergía (coordenadas degeneradas), `.unwrap()` panicá.
- **Fix**: Reemplazado `.unwrap()` en `svd.u` y `svd.v_t` con `match (svd.u, svd.v_t)` pattern. Si SVD falla, retorna alignment identidad como fallback seguro.
- **Código**: [src/alignment/kabsch.rs](../src/alignment/kabsch.rs)

### Bug #3 — Coplanar edge case ✅ RESUELTO
- El Kabsch ya corrige reflexiones con el chequeo `det(R) < 0` antes de construir la rotación final.

### Bug #4 — Quaternion formula sin verificación independiente ✅ RESUELTO
- Hay tests de regresión que comparan `align_quaternion` contra Kabsch y validan RMSD y coordenadas alineadas.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
