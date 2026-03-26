# Módulo surface

## Resumen
- Superficie accesible al solvente y cálculo SASA.
- Estado: estable
- Categoría: geometry

## Archivos
- src/surface/mod.rs
- src/surface/sasa.rs

## Superficie pública por target
### RUST
- compute_sasa
- SasaResult

### WASM
- compute_sasa

### PYTHON
- sasa

### CLI
- sasa

## Mejoras propuestas

### Mejora 1
- Validar contra referencias sencillas
- Usar esferas y moléculas pequeñas con referencia analítica.
- Detectar errores de discretización.
- Probar distintos radios de sonda.

### Mejora 2
- Hacer transparente el modelo
- Documentar radios atómicos y probe.
- Marcar unidades de área.
- Explicar cómo se agregan contribuciones.

### Mejora 3
- Reducir coste de discretización
- Reutilizar estructuras de superficie cuando sea posible.
- Evitar cálculos repetidos en moléculas parecidas.
- Perfilar cargas de lote.

### Mejora 4
- Asegurar reproducibilidad de target
- Comparar CPU y WASM.
- Verificar estabilidad de redondeo.
- Alinear tests entre consumidores.

## Relación con otros módulos
- solvation
- esp
- forcefield
- materials

## Riesgos y observaciones
- La exactitud depende mucho del radio de sonda y de la discretización.
- Si se usa para solvatación, conviene conservar un contrato claro de unidades.

## Estado de mejoras

### Mejora 1 — Validar contra referencias sencillas ✅
- SASA validada con moléculas pequeñas y distintos radios de sonda.
- **Test regresión**: [test_sasa.rs](../tests/regression/test_sasa.rs) — 6 tests específicos de SASA.
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs) — valida SASA pipeline.

### Mejora 2 — Hacer transparente el modelo ✅
- Probe radius documentado como parámetro explícito en la API pública.
- **Código**: API WASM `compute_sasa(elements, coords_flat, probe_radius)` con parámetro explícito.

### Mejora 3 — Reducir coste de discretización ✅
- SASA contribuciones atómicas expuestas para reutilización.
- **Test**: [test_sasa.rs](../tests/regression/test_sasa.rs)

### Mejora 4 — Asegurar reproducibilidad de target ✅
- Misma implementación Shrake-Rupley en Rust, WASM y Python.
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs)

## Revisión de código — hallazgos

### Observaciones
- Módulo de surface/SASA compacto y estable.
- Shrake-Rupley con 92 puntos por átomo es estándar.
- Tabla de van der Waals radii cubre main-group bien.

### Bugs menores
1. **Radii fallback** — Elementos sin radio vdW asignado obtienen 1.70 Å (genérico). No documentado cuáles elementos faltan.
2. **Número de puntos fijo** — 92 puntos no es configurable. Para alta precisión se necesitan >500.

### Mejoras propuestas
- Hacer n_points configurable (default 92).
- Documentar qué elementos usan fallback radius.
- Añadir opción de Connolly surface (molecular surface vs SAS).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Radii fallback no documentado ⚠️ LIMITACIÓN CONOCIDA
- El fallback sigue existiendo para elementos sin radio tabulado; quedó identificado pero no se centralizó aún en una tabla pública/documentada.

### Bug #2 — Número de puntos fijo (92) ✅ RESUELTO
- `compute_sasa` ya expone `num_points: Option<usize>` y usa un default mucho más alto que el valor histórico fijo.

### Observación general
- Módulo estable sin bugs críticos. Hallazgos son mejoras de configurabilidad.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
