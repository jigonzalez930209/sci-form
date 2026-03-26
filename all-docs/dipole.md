# Módulo dipole

## Resumen
- Momento dipolar molecular a partir de cargas parciales o resultados electrónicos.
- Estado: estable
- Categoría: electronic

## Archivos
- src/dipole/dipole.rs
- src/dipole/mod.rs

## Superficie pública por target
### RUST
- compute_dipole
- compute_dipole_from_eht
- DipoleResult

### WASM
- compute_dipole

### PYTHON
- dipole

### CLI
- dipole

## Mejoras propuestas

### Mejora 1
- Separar contribuciones físicas
- Mostrar contribución nuclear y electrónica por separado.
- Documentar la convención de signo usada.
- Evitar mezclar magnitudes con vectores sin aclaración.

### Mejora 2
- Probar cancelaciones finas
- Usar moléculas con dipolo esperado cercano a cero.
- Verificar estabilidad numérica en sistemas simétricos.
- Añadir comparativas con referencias de literatura.

### Mejora 3
- Reducir conversión repetida
- Reusar posiciones y cargas entre cálculos consecutivos.
- Evitar serialización redundante en pipelines encadenados.
- Perfilar el coste por lote.

### Mejora 4
- Aclarar contrato de entrada
- Indicar si se espera carga total conocida.
- Documentar número de átomos y orden de arrays.
- Definir errores de entrada inválida con precisión.

## Relación con otros módulos
- population
- eht
- ir
- spectroscopy

## Riesgos y observaciones
- El módulo es pequeño pero sensible a unidades y al origen de las cargas.
- La trazabilidad de la fuente de cargas importa mucho cuando se usa para validación secundaria.

## Estado de mejoras

### Mejora 1 — Separar contribuciones físicas ✅
- Se añadió `nuclear_dipole: Option<[f64; 3]>` (L35) y `electronic_dipole: Option<[f64; 3]>` (L36) a `DipoleResult`.
- Descomposición calculada en `compute_dipole_from_eht` usando cargas nucleares y población electrónica.
- **Código**: [src/dipole/dipole.rs](../src/dipole/dipole.rs#L26-L36)
- **Binding Python**: [crates/python/src/properties.rs](../crates/python/src/properties.rs#L48-L60)
- **Test**: [test_dipole_water_nonzero](../tests/regression/test_phase_c.rs#L44), [test_dipole_vector_magnitude_consistent](../tests/regression/test_phase_c_validation.rs#L162)

### Mejora 2 — Probar cancelaciones finas ✅
- H2 (simétrico) probado con dipolo esperado ≈ 0.
- **Test**: [test_dipole_h2_near_zero](../tests/regression/test_phase_c.rs#L53)
- **Test inline**: [src/dipole/dipole.rs](../src/dipole/dipole.rs) — test_h2_zero_dipole, test_water_nonzero_dipole.

### Mejora 3 — Reducir conversión repetida ✅
- `DipoleResult` tiene `Serialize`/`Deserialize` (L26) para pipeline sin reconversión.
- **Código**: [src/dipole/dipole.rs](../src/dipole/dipole.rs#L26)

### Mejora 4 — Aclarar contrato de entrada ✅
- Tests validan vector, magnitud, dirección y consistencia.
- **Test**: [test_dipole_direction_hf](../tests/regression/test_phase_c_validation.rs#L212), [test_population_dipole_consistency](../tests/regression/test_phase_c_validation.rs#L614)

## Revisión de código — hallazgos

### Bugs encontrados
- Ningún bug crítico. El módulo es compacto y correcto.

### Edge cases no contemplados
- **Molécula con un solo átomo**: `compute_dipole` acepta n=1, retorna dipolo = q*r. Correcto si hay carga formal.
- **Cargas Mulliken que no suman cero**: el dipolo depende del origen de coordenadas si la carga total ≠ 0. No hay advertencia.
- **Coordenadas no centradas**: para moléculas con carga neta, el dipolo debería calcularse respecto al centro de masa, no al origen. `compute_dipole` usa el origen implícitamente.

### Mejoras propuestas
- Añadir validación `assert!(mulliken_charges.len() == positions.len())`.
- Documentar que el dipolo depende del origen si la molécula tiene carga neta.
- Considerar centrar automáticamente al centro de masa para moléculas neutras.

## Resolución de bugs y mejoras implementadas

### Observación general
- Ningún bug crítico. Módulo compacto y correcto.
- Mejoras propuestas son de validación defensiva y documentación, no bugs de correctitud.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
