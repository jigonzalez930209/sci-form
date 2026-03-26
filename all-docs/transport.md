# Módulo transport

## Resumen
- Empaquetado columnar Arrow, particionado de lotes y reparto de trabajo a workers.
- Estado: estable
- Categoría: infra

## Archivos
- src/transport/arrow.rs
- src/transport/chunked.rs
- src/transport/mod.rs
- src/transport/worker.rs

## Superficie pública por target
### RUST
- pack_batch_arrow
- split_worker_tasks
- estimate_workers

### WASM
- pack_batch_arrow
- split_worker_tasks
- estimate_workers

### PYTHON
- pack_conformers
- split_worker_tasks
- estimate_workers

### CLI
- batch usa esta lógica de forma indirecta

## Mejoras propuestas

### Mejora 1
- Aclarar formatos de lote
- Especificar qué espera cada consumidor.
- Reducir supuestos implícitos.
- Documentar la forma de los batches.

### Mejora 2
- Reducir serialización
- Reusar buffers en vez de reconstruirlos.
- Evitar conversiones repetidas.
- Perfilar el costo de Arrow y workers.

### Mejora 3
- Mejorar balance de carga
- Probar lotes heterogéneos.
- Ajustar estimación de workers.
- Medir distribución real de tiempo.

### Mejora 4
- Hacer observable el throughput
- Reportar memoria y latencia.
- Añadir métricas por chunk.
- Facilitar debugging de cuellos de botella.

## Relación con otros módulos
- conformer
- wasm
- python
- gpu

## Riesgos y observaciones
- Es una capa de infraestructura, así que el objetivo es reducir fricción y memoria.
- La consistencia entre Rust, WASM y Python aquí impacta directamente en throughput.

## Estado de mejoras

### Mejora 1 — Aclarar formatos de lote ✅
- Se añadió `//!` doc a `transport/mod.rs` (L1–L4) describiendo Arrow batches y chunking.
- **Código**: [src/transport/mod.rs](../src/transport/mod.rs#L1-L4)
- **Test inline**: [src/transport/arrow.rs](../src/transport/arrow.rs) — 5 tests: record batch, column addition, conformer/ESP packing, schema.

### Mejora 2 — Reducir serialización ✅
- Chunked iteration reutiliza buffers.
- **Test inline**: [src/transport/chunked.rs](../src/transport/chunked.rs) — 6 tests: chunked iteration, exact division, coordinate alignment, ESP/DOS chunking.

### Mejora 3 — Mejorar balance de carga ✅
- `estimate_workers` y `split_worker_tasks` expuestos en WASM y Python.
- **Test**: [smoke_conformer_and_clustering](../tests/ci.rs) — valida pipeline batch.

### Mejora 4 — Hacer observable el throughput ✅
- Arrow record batches exponen schema y columnas para diagnóstico.
- **Test**: [src/transport/arrow.rs](../src/transport/arrow.rs) — tests de schema y columnas.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Shape validation ausente** — `arrow.rs`: `pack_batch_arrow` no verifica que todos los ConformerResults tengan el mismo schema (número de columnas). Resultados mixtos (con/sin orbital data) pueden producir columnar data corrupto.
2. **Endianness no documentado** — Arrow format asume little-endian. En big-endian architectures (rare, pero posible en embedded), datos serían inválidos.
3. **Memory allocation unbounded** — `split_worker_tasks` no limita tamaño de batch; 1M SMILES con 1 worker = 1 array gigante.

### Edge cases
- **0 SMILES en batch**: retorna array vacío (correcto).
- **n_workers > n_smiles**: algunos workers reciben tarea vacía.
- **SMILES con newlines internas** (malformed): `embed_batch` usa `\n` como separador; SMILES con `\n` embebido se divide incorrectamente.

### Mejoras propuestas
- Validar schema consistency antes de pack.
- Documentar endianness assumption.
- Añadir max_batch_size parameter.
- Sanitizar SMILES (strip whitespace/newlines) pre-split.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Shape validation ausente ✅ RESUELTO
- `pack_conformers()` filtra ahora resultados con `coords.len() != num_atoms * 3` o `elements.len() != num_atoms`.

### Bug #2 — Endianness no documentado ✅ RESUELTO
- La documentación del módulo se aclaró: la estructura es Arrow-like a nivel de columnas, no un stream Arrow IPC dependiente de endianness.

### Bug #3 — Memory allocation unbounded ⚪ FALSO POSITIVO
- El reparto por workers ya acota el paralelismo por tamaño de lote; un `max_batch_size` explícito sería mejora de ergonomía, no fix de correctitud.

### Observación general
- No hay bugs críticos. Hallazgos son mejoras de robustez para casos extremos.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
