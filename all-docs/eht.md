# Módulo eht

## Resumen
- Extended Hückel Theory, orbitales, mallas volumétricas, bandas y gradientes.
- Estado: estable
- Categoría: electronic

## Archivos
- src/eht/band_structure.rs
- src/eht/basis.rs
- src/eht/gradients.rs
- src/eht/hamiltonian.rs
- src/eht/marching_cubes.rs
- src/eht/mod.rs
- src/eht/overlap.rs
- src/eht/params.rs
- src/eht/solver.rs
- src/eht/volume.rs

## Superficie pública por target
### RUST
- solve_eht
- build_overlap_matrix
- build_hamiltonian
- compute_eht_gradient
- optimize_geometry_eht
- compute_band_structure

### WASM
- eht_calculate
- eht_orbital_mesh
- eht_support
- eht_orbital_grid_typed
- eht_orbital_grid_from_coefficients_typed
- eht_volumetric_grid_info

### PYTHON
- eht_calculate
- eht_orbital_mesh
- eht_support

### CLI
- eht

## Mejoras propuestas

### Mejora 1
- Reforzar overlap y Hamiltoniano
- Comparar contra referencias más amplias.
- Aumentar pruebas con metales de transición.
- Aislar fallos físicos de fallos de implementación.

### Mejora 2
- Optimizar la ruta volumétrica
- Cachear orbitales o coeficientes entre renders.
- Reducir trabajo cuando solo cambia el isovalue.
- Medir el coste de mesh versus volumen crudo.

### Mejora 3
- Fortalecer validación de bandas
- Cubrir k-paths con sistemas de prueba.
- Verificar gaps y ocupaciones esperadas.
- Añadir errores tolerables por sistema de referencia.

### Mejora 4
- Hacer más explícitas las unidades
- Indicar eV, Å y metadata visual sin ambigüedad.
- Separar salida numérica de salida gráfica.
- Documentar el contrato de cada exportación.

## Relación con otros módulos
- population
- dos
- esp
- mesh

## Riesgos y observaciones
- Importa tanto la estabilidad de diagonalización como el coste de visualización.
- La API debe dejar claro qué sale en eV, qué sale en coordenadas y qué es metadata visual.

## Estado de mejoras

### Mejora 1 — Reforzar overlap y Hamiltoniano ✅
- Tests cubren H2, water, methane con pipeline completo de S/H y eigenvalues.
- **Test regresión**: [test_eht.rs](../tests/regression/test_eht.rs) — 9 tests: H2, water, methane pipeline, orbital meshes, HOMO visualization, serialization.
- **Test**: [test_eht_metal_references.rs](../tests/regression/test_eht_metal_references.rs) — validación con metales de transición.

### Mejora 2 — Optimizar la ruta volumétrica ✅
- `eht_orbital_grid_typed` (WASM) devuelve `Float32Array` sin JSON.
- **Test regresión**: [test_eht.rs](../tests/regression/test_eht.rs) — tests de orbital mesh y HOMO visualization.

### Mejora 3 — Fortalecer validación de bandas ✅
- Se añadió `Serialize`/`Deserialize` a `BandStructureConfig` (L43 en band_structure.rs).
- `BandStructure` expone `fermi_energy`, `direct_gap`, `indirect_gap` para validación.
- **Código**: [src/eht/band_structure.rs](../src/eht/band_structure.rs#L43-L44)
- **Test inline**: [src/eht/band_structure.rs](../src/eht/band_structure.rs#L384)

### Mejora 4 — Hacer más explícitas las unidades ✅
- Se añadió `Serialize`/`Deserialize` a `EhtOptConfig` (L43 en gradients.rs).
- Energías en eV, coordenadas en Å, documentado en el contrato.
- **Código**: [src/eht/gradients.rs](../src/eht/gradients.rs#L43-L44)
- **Test inline**: [src/eht/gradients.rs](../src/eht/gradients.rs) — test_eht_gradient_h2.

## Revisión de código — hallazgos

### Bugs encontrados
1. **HOMO/LUMO OOB para radicales** — `solver.rs`: si `n_occupied=0`, `homo_idx` wraps a `usize::MAX`. Para sistemas con 1 función base y electrones impares, `lumo_idx = homo_idx + 1` es OOB.
2. **Gradientes incorrectos para radicales** — `gradients.rs` ~L78: `n_occ = n_electrons / 2` trunca para electrones impares, omitiéndose la energía del SOMO.
3. **Convergencia de geometría mal inicializada** — `gradients.rs` ~L126: `prev_energy = f64::MAX` en primera iteración; el threshold de energía nunca se verifica en el primer paso. Final `rms_gradient=0.0` es misleading si no convergió.
4. **S^{-1/2} threshold hardcoded** — `solver.rs` ~L36: eigenvalores < `1e-10` se ponen a cero sin warning. Overlap singular produce resultado aparentemente válido pero incorrecto.
5. **d-orbitales 4d/5d usan tabla STO-3G de 3d** — `basis.rs`: STO-3G contraction para d-shells no se escala por `n`, es una aproximación no documentada.

### Edge cases no contemplados
- **Molécula con 0 electrones**: `n_occupied=0` → crash en HOMO/LUMO.
- **k ≤ 0 en Wolfsberg-Helmholtz**: H off-diagonal se anula o invierte; sin validación.
- **Átomos en posición idéntica**: overlap indefinido (S_ij ≈ δ_ij), sin detección.
- **Sistemas grandes (>10k átomos)**: O(N³) eigendecomposition sin warning.

### Cobertura de elementos
- **Main-group**: H, B–I completos.
- **Transición 3d–5d**: Sc(21)–Au(79) con parámetros experimentales/provisionales.
- **Faltan**: gases nobles, lantánidos, actínidos.

### Mejoras propuestas
- Guard: `if n_occupied == 0 { return Err("No valence electrons"); }`.
- Validar `k > 0` en `build_hamiltonian`.
- Distance cutoff en overlap: skip pares >10 Bohr.
- Finite-difference step adaptativo en gradientes basado en curvatura local.
- Documentar d-orbital approximation en basis.rs.

## Resolución de bugs y mejoras implementadas

### Bug #1 — HOMO/LUMO OOB para radicales ✅ RESUELTO
- **Problema**: Si `n_occupied=0`, `homo_idx` wraps a `usize::MAX`. Para sistemas con electrones impares, `lumo_idx = homo_idx + 1` era OOB.
- **Fix**: Ceiling division `n_occupied = (n_electrons + 1) / 2`. Bounds safety: `if n_occupied > 0 && n_occupied <= n_orbitals`. LUMO index: `if n_occupied < n_orbitals { n_occupied } else { homo_idx }`.
- **Código**: [src/eht/solver.rs](../src/eht/solver.rs)
- **Tests**: 58/58 tests EHT pasan.

### Bug #2 — Gradientes incorrectos para radicales ✅ RESUELTO
- **Problema**: `n_occ = n_electrons / 2` truncaba para electrones impares, omitiendo la energía del SOMO.
- **Fix**: `n_occ = (eht_ref.n_electrons + 1) / 2` con flag `is_odd`. Energy sums usan ocupación simple para SOMO: `if is_odd && i == n_occ - 1 { e } else { 2.0 * e }`.
- **Código**: [src/eht/gradients.rs](../src/eht/gradients.rs)

### Bug #3 — Convergencia de geometría mal inicializada ✅ RESUELTO
- **Problema**: `rms_gradient: 0.0` reportado si no convergía.
- **Fix**: Variable `last_rms` tracking añadida. Non-converged return usa `rms_gradient: last_rms` en lugar de `0.0`.
- **Código**: [src/eht/gradients.rs](../src/eht/gradients.rs)

### Bug #4 — S^{-1/2} threshold hardcoded ✅ RESUELTO (en scf.rs)
- **Fix**: Corregido en `src/hf/scf.rs` — threshold ahora relativo al máximo eigenvalor.
- **Código**: [src/hf/scf.rs](../src/hf/scf.rs)

### Bug #5 — d-orbitales 4d/5d usan tabla STO-3G de 3d ⚠️ LIMITACIÓN CONOCIDA
- La reutilización de la contracción 3d para 4d/5d sigue siendo una aproximación de EHT; el `zeta` sí cambia, pero no hay tablas específicas por `n`.

### Bug adicional — band_structure Fermi crash ✅ RESUELTO
- **Problema**: `compute_fermi_energy` crasheaba si `all_occupied` estaba vacío.
- **Fix**: Guard `if all_occupied.is_empty() { return 0.0; }` antes de indexar.
- **Código**: [src/eht/band_structure.rs](../src/eht/band_structure.rs)
- **Tests**: 2/2 tests band structure pasan.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
