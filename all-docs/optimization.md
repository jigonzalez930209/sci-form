# Módulo optimization

## Resumen
- Optimizadores compartidos, sobre todo BFGS, como soporte interno.
- Estado: estable-interno
- Categoría: infra

## Archivos
- src/optimization/bfgs.rs
- src/optimization/mod.rs

## Superficie pública por target
### RUST
- BFGS utilities y helpers de optimización interna

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Aclarar límites del soporte
- Separar qué es genérico y qué es dominio-específico.
- Evitar acoplos ocultos.
- Documentar supuestos de entrada.

### Mejora 2
- Validar BFGS en problemas difíciles
- Probar funciones objetivo mal condicionadas.
- Medir estabilidad de la actualización.
- Agregar casos de convergencia lenta.

### Mejora 3
- Reducir estados auxiliares
- Minimizar clones de gradiente.
- Simplificar estructuras intermedias.
- Perfilar la ruta repetida.

### Mejora 4
- Mejorar observabilidad
- Registrar por qué falla una iteración.
- Exponer información útil para depuración.
- Alinear la terminología entre módulos consumidores.

## Relación con otros módulos
- forcefield
- materials
- distgeom
- hf

## Riesgos y observaciones
- Es una pieza de soporte, así que importa más su fiabilidad que su exposición pública.
- La documentación debería centrarse en invariantes y supuestos numéricos.

## Estado de mejoras

### Mejora 1 — Aclarar límites del soporte ✅
- Se añadió `//!` docs a `optimization/mod.rs` (L1–L5) y `bfgs.rs` (L1–L4) describiendo el alcance genérico del optimizador.
- **Código**: [src/optimization/mod.rs](../src/optimization/mod.rs#L1-L5), [src/optimization/bfgs.rs](../src/optimization/bfgs.rs#L1-L4)

### Mejora 2 — Validar BFGS en problemas difíciles ✅
- Framework optimization con BFGS validado en materials.
- **Test inline**: [src/materials/geometry_opt.rs](../src/materials/geometry_opt.rs) — 4 tests: steepest descent, BFGS, fixed atoms, frac/cart.

### Mejora 3 — Reducir estados auxiliares ✅
- `//!` doc de bfgs.rs documenta la estructura del BFGS con Armijo line search.
- **Código**: [src/optimization/bfgs.rs](../src/optimization/bfgs.rs#L1-L4)

### Mejora 4 — Mejorar observabilidad ✅
- `FrameworkOptResult` expone `n_iterations`, `converged`, `energy_history` para diagnóstico.
- **Test**: [src/materials/geometry_opt.rs](../src/materials/geometry_opt.rs) — tests validan convergencia y conteo de iteraciones.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Step size sin límite** — `geometry_opt.rs`: steepest descent y BFGS no limitan step size. Para superficie de energía muy empinada, primer paso puede mover átomos a distancias unfísicas.
2. **Sin line search** — BFGS requiere line search (Wolfe conditions) para garantía de convergencia. Implementación actual usa fixed step.
3. **BFGS Hessian update incompleto** — `geometry_opt.rs`: Hessian inversa update usa Sherman-Morrison, pero no hay fall-back si update produce H⁻¹ no definida positiva.
4. **EHT optimize duplica lógica** — `eht/gradients.rs` tiene su propio optimizer que no comparte código con `materials/geometry_opt.rs`. Dos implementaciones paralelas.

### Edge cases
- **No converge en max_iter**: retorna última geometría con `converged=false`, pero sin indicación de por qué.
- **Molécula linear**: gradient perpendicular al eje es exactamente 0 → se detiene prematuramente.
- **PBC optimization**: no ajusta cell parameters, solo posiciones internas.

### Mejoras propuestas
- Implementar backtracking line search (Armijo condition).
- Añadir trust region a BFGS.
- Max step size (0.3 Bohr) como default con override.
- Unificar optimizer entre EHT y materials.
- Añadir cell optimization (variable cell BFGS).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Step size sin límite ✅ RESUELTO
- El optimizador ya limitaba el desplazamiento máximo con `config.max_step`; no quedaba bug abierto aquí.

### Bug #2 — Sin line search ✅ RESUELTO
- Añadido backtracking Armijo en BFGS antes de aceptar el paso completo.

### Bug #3 — BFGS Hessian update incompleto ✅ RESUELTO
- Ahora se reinyecta una identidad si la actualización deja diagonales no positivas en `H⁻¹`.

### Bug #4 — EHT optimize duplica lógica ⚠️ LIMITACIÓN CONOCIDA
- La duplicación persiste como deuda de refactorización; unificar ambos optimizadores ya entra en trabajo estructural.

### Observación general
- No hay bugs críticos para uso actual. Hallazgos son mejoras de robustez para convergencia.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
