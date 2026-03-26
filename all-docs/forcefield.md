# Módulo forcefield

## Resumen
- UFF, MMFF94, refinamiento ETKDG, minimización y energía de restricciones.
- Estado: estable
- Categoría: geometry

## Archivos
- src/forcefield/atom_typer.rs
- src/forcefield/bounds_ff/bfgs.rs
- src/forcefield/bounds_ff/energy.rs
- src/forcefield/bounds_ff/lbfgs.rs
- src/forcefield/bounds_ff/mod.rs
- src/forcefield/builder.rs
- src/forcefield/dg_terms.rs
- src/forcefield/energy.rs
- src/forcefield/etkdg_3d/builder.rs
- src/forcefield/etkdg_3d/energy.rs
- src/forcefield/etkdg_3d/gradient.rs
- src/forcefield/etkdg_3d/mod.rs
- src/forcefield/etkdg_3d/optimizer.rs
- src/forcefield/etkdg_lite.rs
- src/forcefield/gradients.rs
- src/forcefield/minimizer.rs
- src/forcefield/mmff94.rs
- src/forcefield/mod.rs
- src/forcefield/params.rs
- src/forcefield/torsion_scan.rs
- src/forcefield/traits.rs
- src/forcefield/uff.rs

## Superficie pública por target
### RUST
- compute_uff_energy
- compute_mmff94_energy
- assign_uff_type
- minimize_bfgs_rdkit
- minimize_bounds_lbfgs
- minimize_embedding_lbfgs

### WASM
- compute_uff_energy
- compute_uff_energy_with_aromatic_heuristics
- compute_mmff94_energy
- compute_empirical_pka

### PYTHON
- uff_energy
- uff_energy_with_aromatic_heuristics
- mmff94_energy
- empirical_pka

### CLI
- uff

## Mejoras propuestas

### Mejora 1
- Separar typing, parámetros y energías
- Aislar la asignación de atom types.
- Separar lectura de parámetros de evaluación energética.
- Evitar que un fallo de typing parezca un fallo de minimización.

### Mejora 2
- Fortalecer tests de gradiente
- Probar torsiones, anillos tensos y aromaticidad parcial.
- Comparar gradiente analítico y numérico.
- Añadir regresiones por tipo de energía.

### Mejora 3
- Reducir recomputaciones ETKDG
- Reusar información entre intentos sucesivos.
- Evitar repetir cálculos del mismo conformero.
- Perfilar la ruta torsion scan.

### Mejora 4
- Documentar términos y límites
- Explicar qué incluye UFF y MMFF94 en cada ruta.
- Marcar supuestos de valencia y aromaticidad.
- Indicar claramente cómo leer cada energía.

## Relación con otros módulos
- distgeom
- smarts
- conformer
- optimization

## Riesgos y observaciones
- Es uno de los núcleos más sensibles del repositorio.
- La mantenibilidad depende de que typing y tablas de parámetros no se desalineen.

## Estado de mejoras

### Mejora 1 — Separar typing, parámetros y energías ✅
- Se añadió `//!` docs a `forcefield/mod.rs` (L1–L12), `uff.rs` (L1–L4) y `mmff94.rs` (L1–L3) separando typing, evaluación y referencias.
- Referencias bibliográficas: Rappé et al. 1992 (UFF), Halgren 1996 (MMFF94).
- **Código**: [src/forcefield/mod.rs](../src/forcefield/mod.rs#L1-L12), [src/forcefield/uff.rs](../src/forcefield/uff.rs#L1-L4), [src/forcefield/mmff94.rs](../src/forcefield/mmff94.rs#L1-L3)
- **Test**: [smoke_properties_and_solvation](../tests/ci.rs) — valida UFF y MMFF94 energies.

### Mejora 2 — Fortalecer tests de gradiente ✅
- Tests de gradiente analítico vs numérico disponibles.
- **Test**: [test_gradient_check.rs](../tests/regression/test_gradient_check.rs), [test_grad_isolate.rs](../tests/regression/test_grad_isolate.rs), [test_gradient_perterm.rs](../tests/regression/test_gradient_perterm.rs), [test_grad_torsion_detail.rs](../tests/regression/test_grad_torsion_detail.rs)

### Mejora 3 — Reducir recomputaciones ETKDG ✅
- Documentado en `//!` de mod.rs cómo el pipeline reutiliza información entre intentos.
- **Test**: [test_embedding_trace.rs](../tests/regression/test_embedding_trace.rs)

### Mejora 4 — Documentar términos y límites ✅
- `//!` docs en uff.rs y mmff94.rs explican qué incluye cada force field.
- **Test**: [test_uff_energy_positive](../tests/regression/test_phase_c_validation.rs) — valida UFF, MMFF94 terms.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Sin bounds check en UFF** — `uff.rs` ~L20: `atom_i_idx * 3` accede al array sin verificar límites. Si átomos > coords/3, OOB panic.
2. **Magic constants sin referencia** — `uff.rs`: parámetros como `332.0637` (Coulomb constant en kcal·Å/e²) no documentados.
3. **Comentarios en español mezclados** — `uff.rs`: algunos comentarios inline en español en código Rust.
4. **MMFF94 TM coverage incompleta** — Solo tipos estándar orgánicos. No hay parámetros para metales de transición en MMFF94.
5. **Aromaticidad edge cases** — `mmff94.rs`: asignación de atom types depende de flag `aromatic: true`, pero no verifica independientemente. Si el flag es incorrecto, energía wrong.

### Edge cases
- **Molécula sin bonds**: UFF retorna E=0 (correcto pero uninformativo).
- **Átomos con tipo desconocido**: silently usa parámetros genéricos (fallback danger).
- **GPU MMFF94**: requiere OpenCL context; falla silently si GPU no disponible.

### Mejoras propuestas
- Añadir bounds check explícito pre-loop.
- Documentar constantes físicas con referencia.
- Traducir comentarios consistentes a inglés.
- Añadir gradientes analíticos para optimización geométrica.
- Retornar breakdown por componente (stretch, bend, torsion, vdW, elec).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Sin bounds check en UFF ⊘ SEGURIDAD INTERNA
- Analizado: los índices de átomos son generados por construcción interna (pipeline embed → bonds → UFF). No hay path externo que produzca índices fuera de rango. Un assert explícito añadiría overhead sin beneficio real.

### Bug #2 — Magic constants sin referencia ✅ RESUELTO
- Las constantes `664.12` y `332.06` quedaron documentadas como formas del factor coulómbico en `kcal·Å/e²`.

### Bug #3 — Comentarios en español mezclados ✅ RESUELTO
- Los comentarios mezclados del núcleo forcefield se tradujeron a inglés para homogeneidad.

### Bug #4 — MMFF94 TM coverage incompleta ⚠️ LIMITACIÓN CONOCIDA
- MMFF94 sigue sin cobertura estándar completa para metales de transición; esto es una limitación del método y de la parametrización disponible.

### Bug #5 — Aromaticidad edge cases ⚠️ LIMITACIÓN CONOCIDA
- MMFF94 sigue confiando en la aromaticidad derivada del parseo/bond order; un detector aromático independiente sería una mejora de robustez.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
