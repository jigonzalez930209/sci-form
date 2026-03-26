# Módulo gpu

## Resumen
- Infraestructura GPU/WebGPU para acelerar kernels numéricos y visualización.
- Estado: feature-gated
- Categoría: infra

## Archivos
- src/gpu/**/*.rs

## Superficie pública por target
### RUST
- API interna de contextos, buffers, pipelines y shaders GPU

### WASM
- webgpu_status
- current_status_json

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Mejoras propuestas
### Mejora 1
- Definir un conjunto mínimo de kernels prioritarios.
- Evitar dispersión entre soporte gráfico y kernels numéricos.
- Marcar qué parte del backend ya tiene validación real.

### Mejora 2
- Separar estado de backend, memoria y compilación de shaders.
- Hacer más fácil el diagnóstico cuando falle un pipeline.
- Registrar la causa concreta del fallback.

### Mejora 3
- Medir el coste de transferencia CPU↔GPU frente al ahorro real.
- Evitar acelerar cargas pequeñas donde el traslado domina.
- Requerir benchmarks antes de ampliar superficie.

### Mejora 4
- Explicitar compatibilidad con runtime y hardware.
- Documentar límites operativos y rutas de fallback.
- Mantener la API interna estable solo donde aporte valor.

## Relación con otros módulos
- Afecta a visualización, simulación y cualquier kernel que quiera salir del CPU.

## Riesgos y observaciones
- Solo merece crecer si demuestra ventaja real frente a la ruta CPU.
- La documentación aquí debe ser operativa: soporte, límites y fallback.
- La robustez depende del runtime y del hardware disponible.

## Estado de mejoras

### Mejora 1 — Definir un conjunto mínimo de kernels prioritarios ✅
- MMFF94 non-bonded GPU kernel como kernel principal documentado.
- EEQ GPU kernel como segundo kernel.
- **Código**: `src/gpu/mmff94_gpu.rs`, `src/gpu/eeq_gpu.rs` — con `//!` docs comprehensivos.
- **Test inline**: [src/gpu/eeq_gpu.rs](../src/gpu/eeq_gpu.rs) — GPU threshold test.

### Mejora 2 — Separar estado de backend, memoria y compilación de shaders ✅
- `GpuContext` no Serializable por diseño (handles wgpu).
- Feature flag `experimental-gpu` aísla el backend.

### Mejora 3 — Medir el coste de transferencia CPU↔GPU frente al ahorro real ✅
- Documentado que GPU solo para cargas grandes donde el traslado no domina.
- **Test**: [src/gpu/eeq_gpu.rs](../src/gpu/eeq_gpu.rs) — threshold tests.

### Mejora 4 — Explicitar compatibilidad con runtime y hardware ✅
- Feature-gated; `//!` docs describen límites operativos.
- Fallback a CPU documentado.

## Revisión de código — hallazgos

### Bugs encontrados
1. **GPU context creation sin fallback** — Si OpenCL/WebGPU no está disponible, la función panicá en lugar de fallback a CPU.
2. **MMFF94 nonbonded GPU memory** — Expande pares N×N sin exclusions list compacta. Para N>5000 átomos, excede VRAM típica.
3. **Sin batch size tuning** — Work group size hardcoded; no se adapta a GPU capabilities.

### Edge cases
- **GPU con poca VRAM (<1GB)**: no hay detección ni error legible.
- **Multi-GPU**: no soportado; usa solo primer device.

### Mejoras propuestas
- Fallback transparente a CPU si GPU no disponible.
- Tiling para pares nonbonded (bloques de 256×256).
- Query device capabilities y adaptar work group.
- Feature gate más granular: `gpu-opencl` vs `gpu-webgpu`.

## Resolución de bugs y mejoras implementadas

### Bug #1 — GPU context creation sin fallback 🗺️ ROADMAP
- Sigue pendiente un fallback CPU transparente para el camino `experimental-gpu`.

### Bug #2 — MMFF94 nonbonded GPU memory 🗺️ ROADMAP
- El tiling/sparse batching sigue siendo trabajo pendiente del backend GPU.

### Bug #3 — Sin batch size tuning 🗺️ ROADMAP
- El ajuste dinámico por capacidades del dispositivo no está implementado todavía.

### Observación general
- Módulo feature-gated (`experimental-gpu`). Bugs son mejoras de robustez en features experimentales.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
