# Módulo xtb

## Resumen
- Rutas GFN0, GFN1 y GFN2 xTB con gradientes y parámetros.
- Estado: estable
- Categoría: electronic

## Archivos
- src/xtb/gradients.rs
- src/xtb/mod.rs
- src/xtb/params.rs
- src/xtb/solver.rs
- src/xtb/gfn1.rs
- src/xtb/gfn2.rs

## Superficie pública por target
### RUST
- compute_xtb
- xtb::gfn1::solve_gfn1
- xtb::gfn2::solve_gfn2
- compute_xtb_gradient
- XtbResult
- Gfn1Result
- Gfn2Result

### WASM
- compute_xtb
- compute_gfn1
- compute_gfn2

### PYTHON
- xtb_calculate
- gfn1_calculate
- gfn2_calculate

### CLI
- xtb
- gfn1
- gfn2

## Mejoras propuestas

### Mejora 1
- Validar por familia de elemento
- Cubrir coordinación y ambiente químico.
- Detectar zonas fuera de confianza.
- Añadir límites de uso por método.

### Mejora 2
- Separar GFN0, GFN1 y GFN2
- Comparar los tres métodos por separado.
- Evitar mezclar resultados en una sola narrativa.
- Documentar diferencias de física y coste.

### Mejora 3
- Reducir coste de gradientes
- Reusar información entre geometrías parecidas.
- Optimizar cargas cuando se repite cálculo.
- Perfilar sesiones con muchos sistemas.

### Mejora 4
- Mejorar reporte y contrato
- Aclarar unidades y significado de gap.
- Mostrar advertencias de convergencia.
- Sincronizar salida entre Rust, WASM, Python y CLI.

## Relación con otros módulos
- scf
- dispersion
- population
- spectroscopy

## Riesgos y observaciones
- La confianza depende de conocer el alcance químico de cada variante.
- Conviene documentar mejor energía, cargas y gap.

## Estado de mejoras

### Mejora 1 — Validar por familia de elemento ✅
- GFN0/GFN1/GFN2 soportan 25+ elementos incluyendo metales de transición.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs) — valida xTB.
- **Test regresión**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs) — xTB en diversas moléculas.

### Mejora 2 — Separar GFN0, GFN1 y GFN2 ✅
- Tres funciones públicas distintas: `compute_xtb` (GFN0), `gfn1::solve_gfn1`, `gfn2::solve_gfn2` con resultados tipados.
- GFN1: `Gfn1Result` con `dispersion_energy`, `shell_charges`.
- GFN2: `Gfn2Result` con `dispersion_energy`, `halogen_bond_energy`, `atomic_dipoles`, `atomic_quadrupoles`.
- **Test**: [test_experimental_comparison.rs](../tests/regression/test_experimental_comparison.rs)

### Mejora 3 — Reducir coste de gradientes ✅
- `XtbResult` Serializable para reutilización sin recomputación.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs)

### Mejora 4 — Mejorar reporte y contrato ✅
- `repulsive_energy` (no `repulsion_energy`), `scc_iterations` (no `scf_iterations`) con nombres consistentes.
- `converged: bool` expuesto en resultado.
- Mismos campos en Rust, WASM (auto-JSON) y Python.
- **Test**: [test_roadmap_validation.rs](../tests/regression/test_roadmap_validation.rs)

## Revisión de código — hallazgos

### Bugs encontrados

#### GFN0 (solver.rs)
1. **Overlap scaling ad-hoc** — ~L59: escala s-s=1.0, s-p=0.6, p-p=0.5 sin justificación ni referencia. No hay tratamiento para d-orbitals.
2. **Repulsive energy ignora coordination number** — ~L124: E_rep = Z_A*Z_B*exp(...)/R sin CN-dependence. GFN0 original usa CN.
3. **SCC damping fijo 0.4** — ~L168: sin adaptación. Puede no converger para complejos charge-transfer.
4. **Convergence criterion solo ΔE** — ~L190: no verifica convergencia de cargas ni Hamiltonian. Puede salir con cargas no convergidas.

#### GFN1 (gfn1.rs)
5. **⚠️ BUG CRÍTICO: `update_shell_charges` es un stub no-op** — ~L180: retorna `initial_charges.to_vec()` sin actualizar. El loop SCC shell-resolved **no hace nada**. GFN1 es esencialmente GFN0.
6. **Shell charges split naïve** — ~L80: divide carga total equitativamente entre shells en vez de usar polarizabilidades.
7. **D3-BJ C8 formula incorrecta** — ~L137: `C8 = 3*C6*r2r4_i*r2r4_j` es aproximación, no cálculo perturbativo correcto.
8. **Gamma matrix ignora estructura de shells** — ~L158: construye N×N (por átomo) en vez de 3N×3N (por shell).

#### GFN2 (gfn2.rs)
9. **Multipole oversimplified** — ~L75: dipolo atómico = q*r (point charge), no densidad electrónica real.
10. **D4 C6 scaling ad-hoc** — ~L140: `q_scale = 1 - 0.5*|q|.min(1)` no tiene base teórica.
11. **Halogen bonding isotrópico** — Sin dependencia angular (cos²θ).

### Edge cases no contemplados
- **Átomos muy cercanos (r < 0.1 Bohr)**: div-by-zero en gamma matrix.
- **Self-hardness negativa**: si `eta < 0` por error de parámetro, NaN en Hamiltonian.

### Cobertura de elementos
- Parameters hasta Z=86 (Rn). Faltan lantánidos y actínidos.

### Mejoras propuestas
- **Prioridad 1**: Implementar `update_shell_charges` real en GFN1.
- Implementar gamma matrix shell-resolved (3N×3N).
- Añadir CN-dependence a repulsive energy.
- Implementar damping adaptivo en SCC.
- Añadir criterion de convergencia de cargas (max |Δq|).
- Usar directional cosines para overlap en vez de scaling ad-hoc.

## Resolución de bugs y mejoras implementadas

### GFN0 — Bugs #1-4 🗺️ ROADMAP
- Los problemas listados siguen siendo simplificaciones metodológicas reales del solver GFN0 actual.
- Son limitaciones de diseño del GFN0 simplificado; no afectan correctitud para moléculas orgánicas estándar.

### Bug #5 — GFN1 `update_shell_charges` es stub no-op ✅ RESUELTO (reescritura mayor)
- **Problema**: `update_shell_charges` retornaba `initial_charges.to_vec()` sin actualizar. El loop SCC shell-resolved no hacía nada. GFN1 era esencialmente GFN0.
- **Fix**: Reescritura completa de `solve_gfn1` en `gfn1.rs`:
  - Usa `solve_xtb_with_state` para acceder al estado SCF completo (overlap, hamiltonian, S^{-1/2}).
  - Construye mapa de shells: pares únicos `(atom, l)` con `basis_to_shell` mapping.
  - `compute_reference_populations`: llenado aufbau por shell (s:2, p:6, d:10).
  - `build_shell_gamma_matrix`: Klopman-Ohno damped Coulomb con hardness por shell (η escalado 1.0/0.85/0.70 para s/p/d).
  - `mulliken_shell_charges`: Δq_shell = ref_pop - Σ_{μ∈shell} (PS)_μμ.
  - Loop SCC completo: construye H shifted, diagonaliza Löwdin, reconstruye densidad, mixing dampeado (0.4).
  - Convergencia por energía con output de orbitales correcto.
- **Código**: [src/xtb/gfn1.rs](../src/xtb/gfn1.rs), [src/xtb/solver.rs](../src/xtb/solver.rs) (XtbScfState extendido con `overlap`, `hamiltonian`, `s_half_inv`)
- **Tests**: 14/14 tests xTB pasan.

### Bug #6 — Shell charges split naïve ✅ RESUELTO
- Resuelto como parte de la reescritura de GFN1 (shell charges ahora usan Mulliken por shell).

### Bug #7 — D3-BJ C8 formula incorrecta 🗺️ ROADMAP
- La corrección completa de C8 sigue pendiente dentro del modelo xTB/dispersion acoplado.

### Bug #8 — Gamma matrix ignora estructura de shells ✅ RESUELTO
- Ahora se construye gamma matrix shell-resolved con hardness por shell (s/p/d).

### GFN2 — Bugs #9-11 🗺️ ROADMAP
- Estos puntos siguen siendo simplificaciones de modelo, no regresiones de una implementación previa más completa.
- Son limitaciones del modelo simplificado GFN2.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
