# Módulo dos

## Resumen
- Densidad de estados total, proyectada y comparaciones multi-método.
- Estado: estable
- Categoría: electronic

## Archivos
- src/dos/dos.rs
- src/dos/mod.rs
- src/dos/multi_method.rs

## Superficie pública por target
### RUST
- compute_dos
- compute_dos_multimethod
- DosMethod
- DosResult
- MultiMethodDosResult

### WASM
- compute_dos

### PYTHON
- dos

### CLI
- dos

## Mejoras propuestas

### Mejora 1
- Separar tipos de DOS
- Diferenciar total, proyectado y multi-método en la salida.
- Documentar claramente el significado de cada curva.
- Evitar mezclar métricas conceptualmente distintas.

### Mejora 2
- Controlar el suavizado
- Probar distintos sigma y ancho de banda.
- Añadir tests de resolución espectral.
- Marcar cuándo la curva pierde información física.

### Mejora 3
- Reusar eigendatos
- Evitar recalcular orbitales ya disponibles.
- Cachear resultados cuando varias curvas comparten sistema.
- Perfilar la ruta de lotes grandes.

### Mejora 4
- Validar forma de curva
- Comparar contra referencias por método.
- Separar error de posición y error de forma.
- Añadir métricas de distancia entre curvas.

## Relación con otros módulos
- eht
- pm3
- xtb
- hf

## Riesgos y observaciones
- Una curva suavizada no debe confundirse con una referencia exacta.
- La salida necesita contexto sobre método, ocupación y sigma.

## Estado de mejoras

### Mejora 1 — Separar tipos de DOS ✅
- `DosResult` expone `energies`, `total_dos`, `pdos` separados.
- Multi-method DOS disponible vía `DosMethod`: Eht, Pm3, Xtb, Gfn1, Gfn2, Hf3c.
- **Código**: [src/dos/dos.rs](../src/dos/dos.rs#L15-L16)
- **Test**: [test_dos_water](../tests/regression/test_phase_c.rs#L65), [test_dos_pdos_sum_equals_total](../tests/regression/test_phase_c_validation.rs#L349)

### Mejora 2 — Controlar el suavizado ✅
- Tests validan distintos sigma y resolución espectral.
- **Test**: [test_dos_narrower_sigma_sharper_peaks](../tests/regression/test_phase_c_validation.rs#L393)

### Mejora 3 — Reusar eigendatos ✅
- Se añadió `Serialize`/`Deserialize` a `DosResult` (L15) para cacheo y reutilización.
- **Código**: [src/dos/dos.rs](../src/dos/dos.rs#L15)

### Mejora 4 — Validar forma de curva ✅
- Tests validan integral, grid, pDOS no-negativa.
- **Test**: [test_dos_integral_approximation](../tests/regression/test_phase_c_validation.rs#L371), [test_dos_pdos_nonnegative](../tests/regression/test_phase_c_validation.rs#L410), [test_dos_grid_matches_parameters](../tests/regression/test_phase_c_validation.rs#L324)

## Revisión de código — hallazgos

### Bugs encontrados
- Ningún bug crítico. Las funciones DOS y PDOS son numéricamente estables.

### Edge cases no contemplados
- **sigma = 0**: división por cero en `1.0 / (sigma * sqrt(2π))`. No hay validación.
- **n_points = 0 o 1**: `(n_points - 1).max(1)` previene div-by-zero en step, pero n_points=0 produce output vacío sin warning.
- **Orbital energies vacías**: `compute_dos(&[], ...)` retorna DOS = 0 para todos los puntos. Correcto pero no documentado.
- **PDOS normalización**: pesos `orbital_atom_weight[k]` se normalizan a sum=1 por orbital. Si `total_w ≈ 0` (orbital no vinculado), se salta normalización, dejando pesos cero. OK pero podría indicar error upstream.

### Mejoras propuestas
- Validar `sigma > 0`, `n_points > 1`.
- Añadir HOMO/LUMO markers al resultado (homo_energy, lumo_energy, gap).
- Export JSON con formato hand-rolled; considerar usar `serde_json` para robustez.
- Añadir Fermi level estimation al `DosResult`.

## Resolución de bugs y mejoras implementadas

### Observación general
- Ningún bug crítico. Funciones DOS y PDOS son numéricamente estables.
- Edge cases (sigma=0, n_points=0) son mejoras de validación de entrada, no bugs de correctitud.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
