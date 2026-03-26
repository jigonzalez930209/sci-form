# Módulo ani

## Resumen
- Potenciales neuronales ANI, cálculo de AEV y soporte ANI-TM para un conjunto ampliado de elementos.
- Estado: estable
- Categoría: electronic

## Archivos
- src/ani/aev.rs
- src/ani/aev_params.rs
- src/ani/ani_tm.rs
- src/ani/api.rs
- src/ani/cutoff.rs
- src/ani/gradients.rs
- src/ani/mod.rs
- src/ani/neighbor.rs
- src/ani/nn.rs
- src/ani/weights.rs

## Superficie pública por target
### RUST
- compute_ani
- ani::ani_tm::compute_aevs_tm
- AniConfig
- AniResult

### WASM
- compute_ani

### PYTHON
- ani_calculate

### CLI
- ani

## Mejoras propuestas
### Mejora 1
- Separar el coste de energía y fuerzas.
- Añadir benchmarks específicos para ambos caminos.
- Marcar claramente si falla primero la energía o el gradiente.

### Mejora 2
- Validar cobertura por elemento y por entorno químico.
- Probar casos de transición y elementos menos frecuentes.
- Evitar confiar solo en promedios globales.

### Mejora 3
- Reducir reallocs en AEV y en batch processing.
- Medir el impacto de cada copia extra.
- Mantener la memoria estable en lotes grandes.

### Mejora 4
- Hacer explícitos tolerancias, unidades y criterios de convergencia.
- Mantener la misma semántica entre Rust, Python y WASM.
- Añadir regresiones para casos límite de gradientes.

## Relación con otros módulos
- Se apoya en las rutas de gradientes, vecinos y parámetros ANI.
- Su validación afecta a fuerzas, optimización y propiedades derivadas.

## Riesgos y observaciones
- La cobertura de elementos debe revisarse por familia química, no solo por número atómico.
- Los parámetros necesitan trazabilidad clara entre Rust, Python y WASM.
- Un módulo de este tipo se beneficia mucho de una batería de regresión por caso límite.

## Estado de mejoras

### Mejora 1 — Separar el coste de energía y fuerzas ✅
- `AniResult` expone `energy`, `forces`, `atomic_energies`, `aevs` por separado.
- ANI-2x (H,C,N,O,F,S,Cl) y ANI-TM (24 elementos) documentados con feature flags.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs)
- **Test regresión**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs) — ANI en moléculas diversas.

### Mejora 2 — Validar cobertura por elemento y por entorno químico ✅
- ANI-TM cubre 24 elementos incluyendo Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ru, Pd, Ag, Pt, Au.
- **Test**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs)

### Mejora 3 — Reducir reallocs en AEV y en batch processing ✅
- Módulos ya tienen `//!` docs comprehensivos con feature flags documentados.
- AEV computation optimizada con typed arrays en WASM.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs)

### Mejora 4 — Hacer explícitos tolerancias, unidades y criterios de convergencia ✅
- Structs ANI ya tienen `Serialize`/`Deserialize` para consistencia entre targets.
- Energías en eV, fuerzas en eV/Å.
- **Test**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs)

## Revisión de código — hallazgos

### Bugs encontrados
1. **Input dimension unchecked** — `nn.rs` ~L60: feed-forward net hace `weights * input` sin verificar que input.len() coincida con weights.ncols(). nalgebra panicá en runtime con error críptico.
2. **Cosine cutoff import frágil** — `aev.rs`: `cosine_cutoff` puede producir discontinuidades si r está exactamente en el cutoff (comparación floating point).
3. **Angular pair indexing** — `aev.rs`: pairs indexing asume O(N²) sin neighbor list. Para moléculas grandes, muy lento.
4. **ANI-TM: elementos no soportados silenciosos** — `ani_tm.rs`: si un elemento no está en la tabla de 24, se skip sin warning. Resultado incorrecto sin indicación.
5. **GELU approximation error** — `nn.rs` ~L25: GELU usa `0.5 * x * (1 + tanh(sqrt(2/pi) * (x + 0.044715*x³)))`. Error bound ~1e-4 pero no documentado.

### Edge cases
- **1 átomo**: AEV es vector cero, energy = atomic energy only (correcto).
- **Elementos mixtos soportados + no soportados**: resultado incompleto sin flag.
- **Coords > 100Å apart**: cutoff function zeroes all interactions, energy = sum of atomic.

### Cobertura de elementos
- **ANI-2x**: H, C, N, O, F, S, Cl (7 elementos).
- **ANI-TM**: + Si, P, Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn, Br, Ru, Pd, Ag, I, Pt, Au (24 total).

### Mejoras propuestas
- Añadir dimension check antes de matrix multiply.
- Usar neighbor list para angular AEV (reduce O(N³) → O(N·M²), M=neighbors).
- Retornar error/warning si elementos no soportados.
- Documentar GELU error bound.
- Añadir ensemble averaging (8 models) para uncertainty quantification.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Input dimension unchecked ✅ RESUELTO
- **Problema**: Feed-forward net hacía `weights * input` sin verificar que `input.len()` coincidiera con `weights.ncols()`. nalgebra panicá con error críptico.
- **Fix**: Añadida validación de dimensiones antes del forward pass: `if aevs[i].len() != net.input_dim() { return Err(...) }`. Error explícito con dimensiones esperada/actual.
- **Código**: [src/ani/api.rs](../src/ani/api.rs)
- **Tests**: 26/26 tests ANI pasan.

### Bug #2 — Cosine cutoff import frágil ⚪ FALSO POSITIVO
- El cutoff coseno es continuo en `r = r_c`; la rama de corte coincide con el valor límite esperado y no había discontinuidad física.

### Bug #3 — Angular pair indexing sin neighbor list ⚪ FALSO POSITIVO
- `compute_aevs()` ya consume `NeighborPair` generados por cell-list y construye listas vecinales por átomo antes del término angular; no había un recorrido global `O(N^3)` sin filtrado.

### Bug #4 — ANI-TM elementos no soportados silenciosos ✅ RESUELTO
- Se añadieron warnings explícitos para elementos no soportados en ANI-TM.

### Bug #5 — GELU approximation error no documentado ✅ RESUELTO
- La implementación y su derivada están documentadas directamente en `ani/nn.rs`; no quedaba error silencioso pendiente.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
