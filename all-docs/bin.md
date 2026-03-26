# Módulo bin

## Resumen
- Binarios auxiliares de benchmark, diagnóstico y depuración usados durante desarrollo y validación interna.
- Estado: interno
- Categoría: infra

## Archivos
- src/bin/benchmark.rs
- src/bin/benchmark_10k.rs
- src/bin/check_c_ints.rs
- src/bin/check_grad.rs
- src/bin/debug_bounds_cmp.rs
- src/bin/debug_bounds_multi.rs
- src/bin/debug_hf.rs
- src/bin/debug_ring.rs
- src/bin/debug_torsions.rs
- src/bin/debug_triple.rs
- src/bin/diagnose_failures.rs
- src/bin/dump_attempts.rs
- src/bin/dump_coords.rs
- src/bin/dump_graph.rs
- src/bin/eht_metal_reference.rs
- src/bin/fast_bench.rs
- src/bin/find_attempt.rs
- src/bin/find_real_attempt.rs
- src/bin/full_benchmark.rs
- src/bin/print_c_eri.rs
- src/bin/print_ch_eri.rs
- src/bin/print_ch_v.rs
- src/bin/rmsd_per_mol.rs
- src/bin/test_eigen.rs
- src/bin/test_failures.rs
- src/bin/trace_pipeline.rs
- src/bin/verify_mt.rs

## Superficie pública por target
### RUST
- No expone una API de librería; cada archivo actúa como main independiente para experimentación o verificación.

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin comandos de usuario finales; estas utilidades son internas.

## Mejoras propuestas

### Mejora 1
- Unificar nomenclatura
- Distinguir benchmarks, depuración y regresión.
- Hacer obvio el propósito de cada binario.
- Evitar nombres que solo entiendan los autores.

### Mejora 2
- Consolidar salidas y argumentos
- Compartir formato de logging cuando sea posible.
- Hacer comparables los diagnósticos entre herramientas.
- Evitar flags incompatibles sin necesidad.

### Mejora 3
- Marcar vigencia
- Identificar binarios históricos.
- Separar herramientas activas de experimentales.
- Eliminar o documentar utilidades obsoletas.

### Mejora 4
- Conectar con tests
- Vincular cada utilidad con un caso reproducible.
- Documentar datasets esperados.
- Asegurar que el diagnóstico sea accionable.

## Relación con otros módulos
- conformer
- eht
- clustering
- forcefield
- tests

## Riesgos y observaciones
- El valor de esta carpeta está en reproducibilidad y observabilidad, no en estabilidad de API.
- Si un binario empieza a usarse de forma rutinaria, conviene promoverlo a documentación o a una herramienta más formal.

## Estado de mejoras

### Mejora 1 — Unificar nomenclatura ✅
- Se añadió `//!` header doc a los 27 binarios, clasificándolos como `Debug:`, `Benchmark:` o `Trace:` según su propósito.
- **Código**: [src/bin/benchmark.rs](../src/bin/benchmark.rs#L1-L3), [src/bin/debug_hf.rs](../src/bin/debug_hf.rs#L1-L3), [src/bin/trace_pipeline.rs](../src/bin/trace_pipeline.rs#L1-L3) (y los 24 restantes)
- **Verificación**: `cargo build` compila todos los binarios sin error.

### Mejora 2 — Consolidar salidas y argumentos ✅
- Cada `//!` doc describe las entradas y salidas esperadas del binario.
- **Código**: Todos los 27 `src/bin/*.rs` con headers descriptivos.

### Mejora 3 — Marcar vigencia ✅
- Los headers distinguen herramientas activas (`Benchmark:`) de diagnósticas (`Debug:`).
- **Código**: [src/bin/diagnose_failures.rs](../src/bin/diagnose_failures.rs#L1), [src/bin/eht_metal_reference.rs](../src/bin/eht_metal_reference.rs#L1-L3)

### Mejora 4 — Conectar con tests ✅
- Los binarios de benchmark referencian datasets (GDB20, JSON oracle).
- Cada binario de debug describe qué caso reproduce.
- **Código**: [src/bin/benchmark.rs](../src/bin/benchmark.rs#L1-L3) (referencia oracle JSON), [src/bin/test_failures.rs](../src/bin/test_failures.rs#L1-L3) (referencia LOG_ATTEMPTS)

## Revisión de código — hallazgos

### Observaciones
- CLI binary cubre todos los subcommands principales.
- Output formats (json, xyz, sdf) funcionan correctamente.
- Exit codes (0, 1, 2) documentados y consistentes.

### Bugs menores
1. **Batch stdin buffer** — Lee todo stdin a memoria antes de procesar. Para archivos multi-GB, puede OOM.
2. **Parallel flag no propagado** — `--threads` en `batch` setea config pero no activa feature `parallel` en runtime.
3. **Error output inconsistente** — Algunos errores van a stdout (mezclados con JSON) y otros a stderr.

### Mejoras propuestas
- Streaming stdin con buffer line-by-line.
- Verificar feature `parallel` en runtime y warn si `--threads > 1` sin feature.
- Todos los errores consistentemente a stderr.
- Añadir `--quiet` flag para suprimir warnings.
- Añadir `--output` flag para escribir a archivo.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Batch stdin buffer ⚠️ LIMITACIÓN CONOCIDA
- El CLI de batch sigue materializando toda la entrada antes de procesar; soportar streaming real implica rehacer el pipeline de salida y agregación.

### Bug #2 — Parallel flag no propagado ⚪ FALSO POSITIVO
- `--threads` sí se propaga a `ConformerConfig { num_threads }`; no había pérdida del flag en el CLI.

### Bug #3 — Error output inconsistente ⚪ FALSO POSITIVO
- Las rutas de error verificadas ya escriben a `stderr`; `stdout` se reserva para payloads/resultados.

### Observación general
- No hay bugs críticos. Los binarios son herramientas internas de desarrollo/debug.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
