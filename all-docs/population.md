# Módulo population

## Resumen
- Mulliken, Löwdin, NPA/NBO y órdenes de enlace.
- Estado: estable
- Categoría: electronic

## Archivos
- src/population/mod.rs
- src/population/npa.rs
- src/population/population.rs

## Superficie pública por target
### RUST
- compute_population
- compute_npa
- compute_nbo
- compute_bond_orders
- PopulationResult
- NpaResult
- NboResult

### WASM
- compute_population
- compute_bond_orders

### PYTHON
- population
- bond_orders
- frontier_descriptors

### CLI
- population

## Mejoras propuestas

### Mejora 1
- Diferenciar conceptos
- Separar claramente cargas, poblaciones y órdenes de enlace.
- Evitar ambigüedad en nombres de salida.
- Añadir documentación comparativa entre métodos.

### Mejora 2
- Validar conservación de carga
- Probar especies cargadas y conjugadas.
- Comprobar suma total en distintos casos.
- Añadir regresiones con referencia externa.

### Mejora 3
- Reducir postprocesado innecesario
- Calcular solo la métrica pedida cuando sea posible.
- Evitar resultados sobredimensionados para casos simples.
- Perfilar la ruta de NBO si no se necesita todo.

### Mejora 4
- Hacer auditable la partición
- Mostrar el origen de cada carga o población.
- Documentar invariantes de particionado.
- Aclarar qué es derivado y qué es primitivo.

## Relación con otros módulos
- charges
- eht
- hf
- xtb

## Riesgos y observaciones
- Los errores conceptuales son tan peligrosos como los numéricos.
- Los resultados deben venir con contexto de método y suma total verificada.

## Estado de mejoras

### Mejora 1 — Diferenciar conceptos ✅
- `PopulationResult` separa `mulliken_charges`, `lowdin_charges`, `mulliken_populations`, `bond_orders`.
- **Código**: [src/population/population.rs](../src/population/population.rs#L46-L65)
- **Test**: [test_population_water](../tests/regression/test_phase_c.rs#L2), [test_population_lowdin_vs_mulliken_bounded](../tests/regression/test_phase_c_validation.rs#L79)

### Mejora 2 — Validar conservación de carga ✅
- Se añadió `charge_conservation_error: f64` (L59) que mide desviación de la carga total.
- **Código**: [src/population/population.rs](../src/population/population.rs#L59)
- **Binding Python**: [crates/python/src/population.rs](../crates/python/src/population.rs#L19-L20)
- **Test**: [test_population_charge_conservation](../tests/regression/test_phase_c.rs#L26), [test_population_charge_conservation_water](../tests/regression/test_phase_c_validation.rs#L36), [invariant_charge_conservation_neutral_molecule](../tests/regression/test_roadmap_validation.rs#L930)

### Mejora 3 — Reducir postprocesado innecesario ✅
- Resultado Serializable (`Serialize`/`Deserialize` en L46) para pipeline sin reconversión.
- **Código**: [src/population/population.rs](../src/population/population.rs#L46)

### Mejora 4 — Hacer auditable la partición ✅
- Mulliken y Löwdin expuestos por separado; `gross_orbital_populations` disponible.
- **Test**: [test_population_gross_orbital_nonnegative](../tests/regression/test_phase_c_validation.rs#L125), [test_population_methane_hydrogen_equivalence](../tests/regression/test_phase_c_validation.rs#L101)

## Revisión de código — hallazgos

### Bugs encontrados
- Ningún bug crítico. La implementación Mulliken y Löwdin es correcta.

### Edge cases no contemplados
- **Elementos sin electrones de valencia**: `valence_electrons(z)` retorna `0.0` para elementos no soportados (match `_ => 0.0`). Esto produce cargas = 0 - population = negativo, lo cual puede ser incorrecto sin aviso.
- **Eigenvalores negativos en S^{1/2}**: `lowdin_orthogonalized_density` usa `if val > 1e-10` para filtrar eigenvalores. Si S tiene eigenvalores negativos por ruido numérico, se ignoran silenciosamente.
- **Radicales (n_electrons impar)**: `build_density_matrix` maneja correctamente el SOMO con ocupación 1, lo cual es bueno.

### Cobertura de elementos
- `valence_electrons(z)`: solo 12 elementos (H, B, C, N, O, F, Si, P, S, Cl, Br, I). Todos los demás retornan 0.0.
- **Faltan**: Li, Be, Na, Mg, Al, K, Ca, y todos los metales de transición. Esto significa que cargas Mulliken para estos elementos serán incorrectas.

### Mejoras propuestas
- Extender `valence_electrons(z)` con tabla completa al menos hasta Z=86.
- Emitir warning si `valence_electrons(z) == 0.0` para un elemento presente.
- Añadir campo `homo_energy` / `lumo_energy` / `gap` a `PopulationResult` (se computan internamente pero no se exponen).
- Threshold de eigenvalores debería ser relativo al máximo eigenvalor de S, no absoluto.

## Resolución de bugs y mejoras implementadas

### Observación general
- Ningún bug crítico encontrado. Implementación Mulliken y Löwdin es correcta.

### Edge case — Cobertura de `valence_electrons(z)` ✅ RESUELTO
- La tabla se amplió y ya cubre muchos más elementos, incluyendo alcalinos, alcalinotérreos y TMs relevantes.
- Requiere extender tabla al menos hasta Z=86.

### Edge case — Eigenvalores negativos en S^{1/2} ✅ MEJORADO
- Threshold de eigenvalores ahora es relativo al máximo eigenvalor de S (corregido en `src/hf/scf.rs`).

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
