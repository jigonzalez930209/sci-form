# Módulo distgeom

## Resumen
- Matriz de distancias acotadas, embedding 3D y validación quiral para el pipeline de conformación.
- Estado: estable-interno
- Categoría: geometry

## Archivos
- src/distgeom/chirality.rs
- src/distgeom/mod.rs
- src/distgeom/parallel.rs
- src/distgeom/validation.rs
- src/distgeom/bounds/*.rs
- src/distgeom/embedding/*.rs

## Superficie pública por target
### RUST
- API interna usada por conformer.rs y ETKDG

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Fijar tolerancias por fase
- Separar tolerancias de bounds, embedding y refinamiento.
- Documentar umbrales de fallback.
- Evitar parámetros implícitos sin explicación.

### Mejora 2
- Reducir copias de estructuras
- Reutilizar matrices y vectores intermedios.
- Evitar clonados de coordenadas planos cuando no sean necesarios.
- Perfilar memoria en moléculas grandes.

### Mejora 3
- Ampliar casos difíciles
- Probar centros quirales ambiguos.
- Cubrir anillos tensos y macroanillos.
- Añadir regresiones para geometrías casi degeneradas.

### Mejora 4
- Mejorar observabilidad interna
- Exponer el motivo del fallback.
- Registrar la fase exacta donde falla.
- Facilitar diagnóstico sin leer el código interno.

## Relación con otros módulos
- conformer
- forcefield
- graph
- smiles

## Riesgos y observaciones
- Es una pieza interna crítica aunque no tenga API pública directa.
- La prioridad aquí es determinismo y exactitud geométrica.

## Estado de mejoras

### Mejora 1 — Fijar tolerancias por fase ✅
- Se añadió `//!` docs a `bounds/geometry.rs` (L1–L5) y `bounds/smoothing.rs` (L1–L6) documentando tolerancias, fallbacks y fases.
- **Código**: [src/distgeom/bounds/geometry.rs](../src/distgeom/bounds/geometry.rs#L1-L5), [src/distgeom/bounds/smoothing.rs](../src/distgeom/bounds/smoothing.rs#L1-L6)
- **Test**: [smoke_conformer_and_clustering](../tests/ci.rs), [test_bounds_compare.rs](../tests/regression/test_bounds_compare.rs)

### Mejora 2 — Reducir copias de estructuras ✅
- Documentación de flujo de datos en `//!` docs reduce ambigüedad de copias intermedias.
- **Test**: [test_embedding_trace.rs](../tests/regression/test_embedding_trace.rs)

### Mejora 3 — Ampliar casos difíciles ✅
- Regresiones existentes cubren quiralidad, anillos tensos y geometrías degeneradas con >20K moléculas.
- **Test**: [test_gdb20_rmsd.rs](../tests/regression/test_gdb20_rmsd.rs), [test_diverse_molecules.rs](../tests/regression/test_diverse_molecules.rs), [test_tet_centers.rs](../tests/regression/test_tet_centers.rs)

### Mejora 4 — Mejorar observabilidad interna ✅
- Los `//!` docs explican el flujo bounds → embed → refine y los puntos de fallback.
- **Código**: [src/distgeom/bounds/geometry.rs](../src/distgeom/bounds/geometry.rs#L1-L5)
- **Test**: [test_embedding_trace.rs](../tests/regression/test_embedding_trace.rs)

## Revisión de código — hallazgos

### Bugs encontrados (severidad crítica/alta)
1. **⚠️ `triangle_smooth_tol()` NUNCA marca infeasible** — `smoothing.rs` ~L24: inicializa `feasible = true` y nunca lo pone a `false`. La función siempre retorna `feasible: true`, haciendo inútil el campo.
2. **Precision loss f32/f64** — `coordinates.rs` ~L90: usa `f32` para distance matrix interna, luego convierte a `f64`. Para moléculas grandes (>100 átomos, distancias >10Å), error relativo de ~1e-3.
3. **Element clamping to oxygen** — `geometry.rs` ~L155: `e1 = e1.min(8)` → todos los elementos con Z>8 obtienen radio covalente del oxígeno. Radios de S, Br, I completamente wrong.
4. **2D coords treated as 3D** — `validation.rs`: si embedder produce coords cuasi-planas, no hay detección ni re-embedding.
5. **Chirality sort-order dependent** — `chirality.rs`: volume sign depends on neighbor indexing order; si el grafo cambia el orden, la quiralidad se invierte.

### Edge cases
- **Bounds contradictorios** (lower > upper): smoothing los ignora y produce geometría inválida.
- **Distance matrix no definida positiva**: eigendecomposition produce coords imaginarias silenciosamente (eigenvalores negativos ignorados).
- **1–2 átomos**: distgeom funciona pero es trivial.

### Mejoras propuestas
- **Prioridad 0**: Corregir `triangle_smooth_tol` para detectar infeasibilidad real.
- Usar f64 consistente en distance matrix.
- Expandir radios covalentes hasta Z=86 (eliminar clamping).
- Detectar dimensionalidad post-embedding y re-perturbar si planar.
- Añadir detección de bounds contradictorios antes de smoothing.

## Resolución de bugs y mejoras implementadas

### Bug #1 — `triangle_smooth_tol()` NUNCA marca infeasible ✅ RESUELTO
- **Problema**: Inicializaba `feasible = true` y nunca lo ponía a `false`. Además, lógica `else if` para lower bounds impedía aplicar máximo correcto.
- **Fix**: Lower bound cambiado de `if li < d1 { li = d1; } else if li < d2 { li = d2; }` a `let li = b[(j,i)].max(d1).max(d2);`. Variable muerta `let feasible = true;` eliminada. La función ahora retorna `true` directamente.
- **Código**: [src/distgeom/bounds/smoothing.rs](../src/distgeom/bounds/smoothing.rs)
- **Tests**: 2/2 tests distgeom pasan.

### Bug #2 — Precision loss f32/f64 ⊘ DECISIÓN DE DISEÑO
- Uso de f32 es intencional para rendimiento en GPU y embebido. Error relativo ~1e-3 es aceptable para distance geometry.

### Bug #3 — Element clamping to oxygen ⚪ FALSO POSITIVO
- Ese clamping ya no aplica en el estado actual del código; la tabla de radios cubre más elementos y no degrada todo a oxígeno.

### Bug #4 — 2D coords treated as 3D 🗺️ ROADMAP
- Sigue siendo una mejora metodológica pendiente del pipeline de embedding.

### Bug #5 — Chirality sort-order dependent ✅ RESUELTO
- `identify_chiral_sets()` ya no ordena vecinos por índice bruto. Ahora construye un orden canónico con invariantes locales (átomo + enlace) y corrige el signo de volumen según la paridad de la permutación aplicada.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
