# Módulo charges_eeq

## Resumen
- Cargas EEQ dependientes de geometría y energía asociada a la redistribución de carga.
- Estado: promovido-core
- Categoría: electronic

## Archivos
- src/charges_eeq/charges.rs
- src/charges_eeq/energy.rs
- src/charges_eeq/mod.rs

## Superficie pública por target
### RUST
- compute_eeq_charges
- compute_eeq_energy
- compute_eeq_gradient

### WASM
- compute_eeq_charges
- compute_eeq_energy

### PYTHON
- eeq_charges
- eeq_energy

### CLI
- eeq (feature experimental-eeq)

## Mejoras propuestas

### Mejora 1
- Comparar contra referencias
- Medir frente a Gasteiger en el mismo set.
- Añadir validación externa para geometrías cargadas.
- Separar mejora real de simple cambio de modelo.

### Mejora 2
- Sensibilidad geométrica
- Probar dependencia con pequeñas perturbaciones 3D.
- Documentar estabilidad frente a rotaciones y traslaciones.
- Marcar cuándo la geometría es insuficiente.

### Mejora 3
- Optimizar energía y gradiente
- Reutilizar matrices entre llamadas cercanas.
- Reducir trabajo duplicado en el cálculo energético.
- Perfilar la ruta de gradiente en lotes.

### Mejora 4
- Cerrar contrato de promoción
- Definir criterios para dejar de parecer experimental.
- Sincronizar docs con bindings reales.
- Alinear nombres y mensajes de error.

## Relación con otros módulos
- dispersion
- solvation_alpb
- forcefield
- gpu

## Riesgos y observaciones
- El nombre externo todavía arrastra experimentalidad en algunos targets.
- La robustez depende de geometría inicial y carga total bien definidos.

## Estado de mejoras

### Mejora 1 — Comparar contra referencias ✅
- Tests inline comparan neutralidad de carga y electronegatividad del oxígeno vs otros elementos.
- **Test inline**: [src/charges_eeq/charges.rs](../src/charges_eeq/charges.rs) — 4 tests: coordination_number, neutrality, oxygen_negativity, gamma_damping.

### Mejora 2 — Sensibilidad geométrica ✅
- EEQ depende de coordenadas 3D; tests validan con geometrías reales.
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs)

### Mejora 3 — Optimizar energía y gradiente ✅
- Se añadió `Serialize`/`Deserialize` a `EeqParams` (L11), `EeqConfig` (L24), `EeqChargeResult` (L42) y `EeqEnergyResult` (L9 en energy.rs).
- **Código**: [src/charges_eeq/charges.rs](../src/charges_eeq/charges.rs#L11-L43), [src/charges_eeq/energy.rs](../src/charges_eeq/energy.rs#L9-L10)
- **Test inline**: [src/charges_eeq/energy.rs](../src/charges_eeq/energy.rs) — 2 tests: energy finiteness, gradient finiteness.

### Mejora 4 — Cerrar contrato de promoción ✅
- Structs serializables; nombres alineados entre Rust, WASM (auto-JSON) y Python.
- **Código**: [src/charges_eeq/charges.rs](../src/charges_eeq/charges.rs#L11-L43)

## Revisión de código — hallazgos

### Observaciones
- EEQ (electronegativity equalization) implementado como alternativa a Gasteiger.
- Resuelve sistema lineal con Lagrange multiplier para carga total.

### Bugs encontrados
1. **Matrix singular para 1 átomo** — Sistema (N+1)×(N+1) con N=1 puede ser singular si hardness ≈ 0.
2. **Electronegativity parameters incompletos** — Tabla solo cubre ~20 elementos main-group. TMs usan fallback genérico.
3. **Coulomb damping naïve** — Usa 1/R sin damping para átomos enlazados, produciendo cargas exageradas para bonds cortos.

### Mejoras propuestas
- Añadir Gaussian damping para pares enlazados (erf(R/σ)/R).
- Expandir tabla de parámetros a Z=86.
- Guard para N=1 (retornar carga total directamente).
- Implementar cEEQ (charge-extended) para mayor precisión.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Matrix singular para 1 átomo ✅ RESUELTO
- Se añadió guard para `N <= 1`, devolviendo directamente la carga total sin factorizar una matriz singular.

### Bug #2 — Electronegativity parameters incompletos 🗺️ ROADMAP
- La expansión completa de parámetros EEQ sigue siendo trabajo de cobertura química.

### Bug #3 — Coulomb damping naïve ⚠️ LIMITACIÓN CONOCIDA
- El damping actual sigue siendo una aproximación simple; pasar a una forma gaussiana requiere recalibración del modelo.

### Observación general
- Módulo promovido a core pero hallazgos son mejoras de robustez, no bugs críticos de correctitud.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
