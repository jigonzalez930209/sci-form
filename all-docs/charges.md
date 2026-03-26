# Módulo charges

## Resumen
- Cargas parciales Gasteiger-Marsili y su parametrización PEOE.
- Estado: estable
- Categoría: electronic

## Archivos
- src/charges/gasteiger.rs
- src/charges/mod.rs

## Superficie pública por target
### RUST
- compute_charges
- compute_charges_configured
- GasteigerConfig

### WASM
- compute_charges

### PYTHON
- charges

### CLI
- charges

## Mejoras propuestas

### Mejora 1
- Aclarar dominio químico
- Listar elementos y estados de carga cubiertos.
- Marcar explícitamente casos fuera de dominio.
- Añadir advertencias para topologías problemáticas.

### Mejora 2
- Validar especies difíciles
- Cubrir protonaciones y zwitteriones.
- Probar aromáticos y aniones conjugados.
- Comparar con referencias externas en varios entornos.

### Mejora 3
- Optimizar convergencia
- Reducir trabajo cuando la carga ya converge.
- Evitar iteraciones completas innecesarias.
- Perfilar los casos de lote grande.

### Mejora 4
- Mejorar trazabilidad
- Exponer iteraciones y residuo final.
- Documentar reglas de inicialización.
- Separar warnings de errores duros.

## Relación con otros módulos
- population
- dipole
- solvation
- forcefield

## Riesgos y observaciones
- La precisión depende de parametrización química y no solo de la aritmética.
- El número de iteraciones y el criterio de convergencia deberían verse en la salida.

## Estado de mejoras

### Mejora 1 — Aclarar dominio químico ✅
- `ChargeResult` expone `iterations` y `converged` (L177) para trazabilidad.
- **Código**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs#L169-L178)
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs), [invariant_charge_conservation_neutral_molecule](../tests/regression/test_roadmap_validation.rs#L930)

### Mejora 2 — Validar especies difíciles ✅
- Tests cubren electropositividad, simetría, conservación de carga formal y soporte por elemento.
- **Test inline**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs) — tests: electronegativity_order, symmetry, formal_charge_conservation, element_support.
- **Test regresión**: [test_phase_c_validation.rs](../tests/regression/test_phase_c_validation.rs)

### Mejora 3 — Optimizar convergencia ✅
- Campo `converged: bool` indica si convergió antes del límite de iteraciones.
- **Código**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs#L177)
- **Test**: [invariant_charge_conservation_neutral_molecule](../tests/regression/test_roadmap_validation.rs#L930)

### Mejora 4 — Mejorar trazabilidad ✅
- `ChargeResult` expone `iterations: usize`, `total_charge: f64` y `converged: bool`.
- **Código**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs#L170-L178)
- **Binding Python**: [crates/python/src/properties.rs](../crates/python/src/properties.rs#L8-L16)

## Revisión de código — hallazgos

### Bugs encontrados
1. **Bonds inválidos se saltan silenciosamente** — `gasteiger.rs` ~L258: `if i >= n || j >= n { continue; }` no emite error ni warning. Un bonds list malformado produce cargas incorrectas sin aviso.
2. **Divisor de hardness no validado** — `gasteiger.rs` ~L262: `params[j].chi(1.0)` se usa como divisor en la transferencia de carga. Puede ser cero o negativo para algunos elementos, riesgo de division-by-zero.
3. **Damping schedule no documentado** — El halving `0.5 → 0.25 → 0.125…` es ad-hoc y `max_iter=6` por defecto es muy bajo para sistemas conjugados grandes.

### Edge cases no contemplados
- **Molécula vacía (0 átomos)**: retorna `ChargeResult` vacío con `converged=true`. Correcto pero no documentado.
- **Átomo único**: no hay bonds, delta_q queda cero. Correcto pero converge inmediatamente sin diagnóstico.
- **Anillos aromáticos**: Gasteiger no distingue bond order (doble/triple = single). La distribución de cargas en aromáticos puede ser incorrecta.
- **Self-loops o duplicados en bonds**: no se validan; podrían corromper resultado.

### Cobertura de elementos
- **23 soportados**: H(1), Li(3), Be(4), B(5), C(6), N(7), O(8), F(9), Na(11), Mg(12), Al(13), Si(14), P(15), S(16), Cl(17), K(19), Ca(20), Ga(31), Ge(32), As(33), Se(34), Br(35), I(53).
- **Faltan**: He, Ne, Ar (gases nobles); todos los metales de transición Ti(22)–Zn(30); lantánidos; actínidos.

### Mejoras propuestas
- Validar bonds: no self-loops (`i != j`), no duplicados, índices en rango.
- Validar divisor de hardness: `if params[j].chi(1.0) < 1e-6 { error }`.
- Considerar DIIS o Broyden en lugar de damping fijo.
- Añadir telemetría de iteraciones: max charge transfer por iteración.
- Extender cobertura a metales de transición (parámetros de Hinze-Jaffe).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Bonds inválidos se saltan silenciosamente ✅ RESUELTO
- **Problema**: `if i >= n || j >= n { continue; }` no emitía error ni warning. Un bonds list malformado producía cargas incorrectas sin aviso.
- **Fix**: Cambiado `continue` a `return Err(format!("Bond ({}, {}) references atom outside range 0..{}", i, j, n));`. Ahora retorna error explícito si los índices están fuera de rango.
- **Código**: [src/charges/gasteiger.rs](../src/charges/gasteiger.rs)
- **Tests**: 16/16 tests de cargas pasan.

### Bug #2 — Divisor de hardness no validado ✅ RESUELTO
- El divisor problemático ya se protege contra valores casi nulos antes de actualizar las cargas.

### Bug #3 — Damping schedule no documentado ✅ RESUELTO
- El schedule de damping quedó documentado en el propio módulo (`initial_damping`, halving por iteración).

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
