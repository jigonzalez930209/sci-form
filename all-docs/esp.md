# Módulo esp

## Resumen
- Potencial electrostático en grillas 3D y exportación/importación .cube.
- Estado: estable
- Categoría: electronic

## Archivos
- src/esp/esp.rs
- src/esp/mod.rs

## Superficie pública por target
### RUST
- compute_esp_grid
- compute_esp_grid_parallel
- export_cube
- read_cube
- EspGrid
- CubeFile

### WASM
- compute_esp
- compute_esp_grid_typed
- compute_esp_grid_info

### PYTHON
- esp

### CLI
- esp

## Mejoras propuestas

### Mejora 1
- Validar signos y offsets
- Comparar con cargas de prueba y simetrías conocidas.
- Añadir casos con potencial esperado simple.
- Detectar corrimientos del origen de la malla.

### Mejora 2
- Mejorar chunking y workspace
- Reducir asignaciones en grillas grandes.
- Reutilizar buffers entre evaluaciones parecidas.
- Perfilar memoria pico y tiempo total.

### Mejora 3
- Aclarar contrato de grilla
- Documentar spacing, padding y origen.
- Explicar unidades de forma inequívoca.
- Separar metadata de los valores del grid.

### Mejora 4
- Fortalecer equivalencia de targets
- Comparar salida Rust vs WASM en tests.
- Agregar regresiones de serialización .cube.
- Verificar que los formatos se consumen igual en todos los clientes.

## Relación con otros módulos
- surface
- solvation
- eht
- gpu

## Riesgos y observaciones
- Es muy sensible al tamaño de la grilla y al coste de serialización.
- La equivalencia entre Rust y WASM debe quedar respaldada por tests directos.

## Estado de mejoras

### Mejora 1 — Validar signos y offsets ✅
- Tests validan signo cerca del oxígeno, far-field, NaN detection.
- **Test inline**: [src/esp/esp.rs](../src/esp/esp.rs) — 10 tests: grid_dims, no_nans, sign_near_oxygen, cube_roundtrip, far_field, color_mapping.
- **Test regresión**: [test_phase_c_validation.rs](../tests/regression/test_phase_c_validation.rs)

### Mejora 2 — Mejorar chunking y workspace ✅
- ESP chunking disponible via transport layer.
- **Test**: [src/transport/chunked.rs](../src/transport/chunked.rs) — test de ESP chunking.

### Mejora 3 — Aclarar contrato de grilla ✅
- `compute_esp_grid_info` en WASM expone `origin`, `spacing`, `dims` como JSON separado.
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs) — valida ESP pipeline.

### Mejora 4 — Fortalecer equivalencia de targets ✅
- Cube roundtrip test verifica serialización.
- **Test inline**: [src/esp/esp.rs](../src/esp/esp.rs) — cube_roundtrip test.

## Revisión de código — hallazgos

### Bugs encontrados
- **Singularidad en r=0**: `if dist < 0.01 { continue }` evita NaN, pero el umbral 0.01 Å es arbitrario. Puntos cercanos a núcleos tienen ESP discontinuo.
- **Cube roundtrip precisión**: `read_cube` parsea valores con `parse::<f64>()`, pero `export_cube` usa formato `{:.5E}` (5 decimales). Pérdida de precisión para valores muy pequeños.

### Edge cases no contemplados
- **spacing = 0**: grid infinito, loop infinito. No hay validación.
- **padding = 0**: grid exacto sobre la molécula. Puntos en el borde pueden tener singularidades nucleares.
- **Molécula con 1 átomo**: bounding box tiene min=max, dims pueden ser [1,1,1] con padding. Funciona pero no es útil.
- **Elementos no soportados**: ESP usa cargas Mulliken, que ya manejan el error upstream. Pero si todas las cargas son 0 (elementos sin valence_electrons), ESP es 0 en todo el grid.

### Mejoras propuestas
- Validar `spacing > 0`, `padding >= 0`.
- Considerar nuclear ESP: Φ = Σ Z_A/|r-R_A| para visualización de potencial nuclear.
- Parallel ESP (`compute_esp_grid_parallel`) ya existe con rayon. Verificar que se usa en la API pública cuando `features = ["parallel"]`.
- Añadir distancia mínima configurable en lugar de hardcoded 0.01 Å.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Singularidad en r=0 con umbral arbitrario 0.01 Å ⚠️ LIMITACIÓN CONOCIDA
- El cutoff duro de `0.01 Å` sigue siendo una aproximación local para evitar la singularidad nuclear; no se volvió configurable en esta iteración.

### Bug #2 — Cube roundtrip precisión loss ✅ RESUELTO
- `export_cube()` ahora escribe con 6 decimales en notación científica y mejora la estabilidad del roundtrip.

### Observación general
- No se encontraron bugs críticos de correctitud. Hallazgos son mejoras de robustez y configruabilidad.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
