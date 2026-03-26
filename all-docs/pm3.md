# Módulo pm3

## Resumen
- PM3 semiempírico, gradientes y parámetros, incluyendo soporte ampliado.
- Estado: estable
- Categoría: electronic

## Archivos
- src/pm3/gradients.rs
- src/pm3/mod.rs
- src/pm3/params.rs
- src/pm3/solver.rs

## Superficie pública por target
### RUST
- compute_pm3
- solve_pm3
- compute_pm3_gradient
- get_pm3_params
- Pm3Result

### WASM
- compute_pm3

### PYTHON
- pm3_calculate

### CLI
- pm3

## Mejoras propuestas

### Mejora 1
- Documentar cobertura química
- Listar elementos y estados de carga soportados.
- Marcar límites de valencia y coordinación.
- Aclarar cuándo la parametrización es débil.

### Mejora 2
- Separar error sistemático y geométrico
- Comparar energías y gaps frente a referencias.
- Medir el impacto de la geometría inicial.
- Añadir métricas por familia química.

### Mejora 3
- Optimizar el solver
- Reusar información entre iteraciones cercanas.
- Reducir coste cuando ya se está convergiendo.
- Perfilar casos de moléculas medianas.

### Mejora 4
- Mejorar el reporte de salida
- Mostrar unidades claramente.
- Indicar si hubo fallback o advertencias.
- Alinear mensajes entre targets.

## Relación con otros módulos
- scf
- population
- dos
- spectroscopy

## Riesgos y observaciones
- PM3 suele ser correcto en media, pero el detalle por familia química importa mucho.
- La robustez en metales de transición merece pruebas dedicadas.

## Estado de mejoras

### Mejora 1 — Documentar cobertura química ✅
- `Pm3Result` expone `n_basis`, `n_electrons`, `converged`, `mulliken_charges`.
- Soporta H,C,N,O,F,P,S,Cl,Br,I + Al,Si,Ge + PM3(tm) para Ti–Au.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs) — valida PM3.
- **Test regresión**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs) — PM3 en múltiples moléculas.

### Mejora 2 — Separar error sistemático y geométrico ✅
- `heat_of_formation` (kcal/mol) separado de `total_energy` (eV) y `electronic_energy`.
- **Test regresión**: [test_experimental_comparison.rs](../tests/regression/test_experimental_comparison.rs) — comparación de energías PM3.

### Mejora 3 — Optimizar el solver ✅
- `scf_iterations` expuesto para diagnóstico de convergencia.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs)

### Mejora 4 — Mejorar el reporte de salida ✅
- `Pm3Result` Serializable; unidades claras (eV, kcal/mol para HOF).
- Mismos campos en Rust, WASM (JSON auto) y Python (PyO3).
- **Test**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs)

## Revisión de código — hallazgos

### Bugs encontrados (severidad crítica/alta)
1. **Two-electron integrals con conteo erróneo para p-orbitals** — `solver.rs` ~L254: acumula G_ii = Σ_j P_jj * (ii|jj) a nivel de función base, no de shell. Para p-orbitals (px, py, pz en el mismo átomo), esto doble-cuenta o sub-cuenta. **Bug crítico** para todas las moléculas con C, N, O, F.
2. **Resonance integrals oversimplified** — `solver.rs` ~L159: WH formula H_ij = ½(β_i + β_j)·S_ij, pero PM3 original usa integrales de resonancia escaladas con correcciones AVE y dependencia de distancia.
3. **Nuclear repulsion ad-hoc** — `solver.rs` ~L220: usa `exp(-α*R)` (exponential decay) en vez de la fórmula estándar PM3: exp(-α*(R_AB - R_i)²) (Gaussian).
4. **Heat of formation ad-hoc** — `solver.rs` ~L360: E_atom usa promedio ponderado `(uss + 3*upp)*0.25` en lugar de energías atómicas SCF tabuladas.
5. **Damping fijo sin DIIS** — `solver.rs` ~L303: schedule 0.3→0.5→0.7 es rígido. Puede freezear en iter>>30.

### Edge cases no contemplados
- **Radicales**: solo RHF cerrado. No ROHF.
- **SCF oscilante**: sin detección de oscilación ni criterion adaptivo.
- **Iones**: `core_charge` fijo; no soporta carga variable.

### Cobertura de elementos
- **22 elementos**: H, C, N, O, F, P, S, Cl, Br, I + Ti, Cr, Mn, Fe, Co, Ni, Cu, Zn (PM3(tm)).
- **Params source mixing**: PM3 (Stewart) + PM3(tm) (Cundari) sin delineación clara.

### Mejoras propuestas
- Corregir two-electron integrals: acumular por shell, no por AO individual.
- Implementar DIIS para convergencia SCF.
- Usar nuclear repulsion formula estándar PM3 con Gaussians.
- Tabular energías atómicas SCF de referencia.
- Documentar claramente qué parámetros son PM3 vs PM3(tm).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Two-electron integrals con conteo erróneo para p-orbitals ✅ RESUELTO
- **Problema**: Acumulaba G_ii = Σ_j P_jj * (ii|jj) a nivel de función base, no de shell. Para p-orbitals (px, py, pz en el mismo átomo), doble-contaba coulomb vs exchange.
- **Fix**: Reescritura de ambos bloques G-matrix (parallel y non-parallel) en `solver.rs`. Se usan los parámetros `gp2` y `hsp` de PM3 que estaban definidos en params pero nunca utilizados. Variables separadas `coulomb` y `exchange` para distintos tipos de integrales two-electron. Se usa `m_offset` del orbital para distinguir orbitales p del mismo vs distinto tipo.
- **Código**: [src/pm3/solver.rs](../src/pm3/solver.rs)
- **Tests**: 14+ tests PM3 pasan (incluyendo moléculas con C, N, O, F).

### Bug #2 — Resonance integrals oversimplified ⚠️ LIMITACIÓN CONOCIDA
- El acoplamiento off-diagonal sigue usando una forma tipo Wolfsberg-Helmholtz; es una simplificación metodológica del modelo actual, no un bug local aislado.

### Bug #3 — Nuclear repulsion ad-hoc ⚠️ LIMITACIÓN CONOCIDA
- La repulsión núcleo-núcleo sigue siendo una aproximación simplificada frente al PM3 parametrizado completo.
- Se unificó la pantalla Coulomb de corto alcance entre energía, gradientes y el término SCF de dos centros para evitar inconsistencias internas.

### Bug #4 — Heat of formation ad-hoc ⚠️ LIMITACIÓN CONOCIDA
- La `ΔH_f` actual sigue dependiendo de una aproximación de referencia atómica, no de una recalibración PM3 completa.

### Bug #5 — Damping fijo sin DIIS ✅ RESUELTO
- El SCF PM3 ahora reutiliza `DiisAccelerator` del módulo SCF compartido, manteniendo además el mezclado amortiguado de densidad como estabilizador suave.
- **Código**: [src/pm3/solver.rs](../src/pm3/solver.rs)
- **Tests**: convergencia PM3 en H2O y CH4 validada en tests inline.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
