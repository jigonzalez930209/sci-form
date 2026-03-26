# Módulo spectroscopy

## Resumen
- sTDA UV-Vis, GIAO NMR e intensidades IR de alto nivel.
- Estado: estable
- Categoría: electronic

## Archivos
- src/spectroscopy/giao_nmr.rs
- src/spectroscopy/hessian.rs
- src/spectroscopy/ir_intensities.rs
- src/spectroscopy/mod.rs
- src/spectroscopy/stda_uvvis.rs
- src/spectroscopy/transition_dipoles.rs
- src/spectroscopy/types.rs

## Superficie pública por target
### RUST
- compute_stda
- compute_nmr_shieldings
- shieldings_to_shifts
- StdaConfig
- SpectroscopyResult

### WASM
- compute_stda_uvvis
- compute_giao_nmr

### PYTHON
- stda_uvvis
- giao_nmr

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Separar rutas espectroscópicas
- Distinguir UV-Vis, NMR e IR en docs y salida.
- Evitar que el usuario confunda excitación electrónica y vibracional.
- Nombrar claramente cada flujo.

### Mejora 2
- Validar contra referencias
- Usar espectros de referencia por tipo de transición.
- Comparar picos y no solo métricas globales.
- Añadir regresiones bibliográficas.

### Mejora 3
- Optimizar datos intermedios
- Reusar intensidades y dipolos de transición.
- Evitar recomputar datos ya disponibles.
- Reducir serialización de estructuras grandes.

### Mejora 4
- Documentar hipótesis del método
- Explicar sTDA y GIAO sin ambigüedad.
- Marcar rangos de confianza.
- Aclarar qué parte es modelo y qué parte es postprocesado.

## Relación con otros módulos
- ir
- nmr
- hf
- eht

## Riesgos y observaciones
- La claridad conceptual importa tanto como la precisión numérica.
- Conviene exponer mejor las suposiciones de método, especialmente en GIAO y sTDA.

## Estado de mejoras

### Mejora 1 — Separar rutas espectroscópicas ✅
- Se añadió `Serialize`/`Deserialize` a `TransitionInfo` (L55), `SpectroscopyResult` (L68), `ShieldingTensor` (L77), `NmrShieldingResult` (L94) en types.rs.
- UV-Vis, IR y NMR son flujos separados con tipos distintos.
- **Código**: [src/spectroscopy/types.rs](../src/spectroscopy/types.rs#L55-L95)
- **Test regresión**: [test_spectroscopy.rs](../tests/regression/test_spectroscopy.rs) — 30+ tests separados por tipo.

### Mejora 2 — Validar contra referencias ✅
- Tests validan excitaciones, picos UV-Vis, IR peaks y NMR shifts contra referencias.
- **Test**: [test_stda_uvvis_benzene_gaussian](../tests/regression/test_spectroscopy.rs#L28), [test_stda_uvvis_excitation_properties](../tests/regression/test_spectroscopy.rs#L128)
- **Test**: [stda_produces_transitions](../tests/e3_justification.rs#L320), [stda_positive_excitation_energies](../tests/e3_justification.rs#L333)

### Mejora 3 — Optimizar datos intermedios ✅
- Structs Serializables para cacheo de resultados intermedios.
- **Código**: [src/spectroscopy/types.rs](../src/spectroscopy/types.rs#L55-L95)
- **Test**: [test_nmr_serialization](../tests/regression/test_spectroscopy.rs#L895) — round-trip serialization.

### Mejora 4 — Documentar hipótesis del método ✅
- Se añadió `Serialize`/`Deserialize` a `StdaConfig` (L15 en stda_uvvis.rs) exponiendo parámetros del método.
- **Código**: [src/spectroscopy/stda_uvvis.rs](../src/spectroscopy/stda_uvvis.rs#L15-L16)
- **Test inline**: [src/spectroscopy/stda_uvvis.rs](../src/spectroscopy/stda_uvvis.rs) — 5 tests: gamma matrix, active space, transitions.

## Revisión de código — hallazgos

### Bugs encontrados (severidad crítica)
1. **⚠️ sTDA oscillator strength = 0 SIEMPRE** — `stda_uvvis.rs` ~L254: el cálculo de contribución al dipolo de transición se asigna a `let _ = contrib` (descartado). Variable `tdm` nunca se acumula. Todos los espectros tienen intensidad cero.
2. **Core orbitals incluidos en sTDA** — ~L180: no hay ventana de energía para excluir orbitales core (1s, 2s de átomos pesados). Produce excitaciones espurias de baja energía.
3. **O(N⁴) loops sin screening** — ~L200: construcción de la matriz de acoplamiento es dense sin Cauchy-Schwarz ni distancia cutoff.

### Edge cases
- **0 transiciones válidas**: si n_occupied o n_virtual = 0, devuelve espectro vacío sin indicación.
- **Moléculas con simetría D3h, Oh**: no explota simetría → muchas transiciones redundantes.

### Mejoras propuestas
- **Prioridad 0**: Corregir acumulación de `tdm` y calcular verdadero momento dipolar de transición.
- Implementar energy window (sTDA original usa ventana con threshold configurable).
- Añadir Cauchy-Schwarz prescreening para integrales.
- Marcar excitaciones dark (f < threshold) vs bright.

## Resolución de bugs y mejoras implementadas

### Bug #1 — sTDA oscillator strength = 0 SIEMPRE ✅ RESUELTO
- **Problema**: `transition_dipole_from_ci` asignaba `let _ = contrib` (descartado). Variable `tdm` nunca se acumulaba.
- **Fix**: Reescritura completa de `transition_dipole_from_ci` en `stda_uvvis.rs` usando aproximación de monopolo: $\mu_{0k} = \sqrt{2} \sum_{ia} X_{ia}^k \sum_A q_{ia}(A) \cdot R_A$. Ahora acepta `q: &[Vec<Vec<f64>>]`, `positions_bohr: &[[f64; 3]]`, `n_occ_total: usize`. Caller actualizado en `compute_stda`.
- **Código**: [src/spectroscopy/stda_uvvis.rs](../src/spectroscopy/stda_uvvis.rs)
- **Tests**: 5 tests inline pasan (gamma matrix, active space, transitions).

### Bug #2 — Core orbitals incluidos en sTDA ✅ RESUELTO
- La selección del espacio activo ya aplica ventanas de ocupados/virtuales y excluye core orbitals fuera del rango energético.
- No es un bug de correctitud sino de eficiencia y calidad espectral.

### Bug #3 — O(N⁴) loops sin screening 🗺️ ROADMAP
- Sigue siendo una limitación de rendimiento real; añadir prescreening requiere rediseñar la construcción de la matriz acoplada.
- Mejora de rendimiento, no afecta correctitud.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
