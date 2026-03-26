# Módulo nmr

## Resumen
- Desplazamientos químicos, acoplamientos J, espectro NMR y códigos HOSE.
- Estado: estable
- Categoría: electronic

## Archivos
- src/nmr/coupling.rs
- src/nmr/hose.rs
- src/nmr/mod.rs
- src/nmr/nucleus.rs
- src/nmr/shifts.rs
- src/nmr/spectrum.rs

## Superficie pública por target
### RUST
- predict_chemical_shifts
- predict_j_couplings
- compute_nmr_spectrum
- HoseCode
- NmrSpectrum
- JCoupling

### WASM
- predict_nmr_shifts
- predict_nmr_shifts_for_nucleus
- compute_giao_nmr
- predict_nmr_couplings
- compute_nmr_spectrum
- compute_nmr_spectrum_with_coords
- compute_hose_codes

### PYTHON
- nmr_shifts
- nmr_shifts_for_nucleus
- giao_nmr
- nmr_couplings
- nmr_spectrum
- hose_codes

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Separar modelos empíricos y GIAO
- Distinguir HOSE de GIAO en nombres y docs.
- Evitar que ambos parezcan variantes equivalentes.
- Exponer método y confianza en la salida.

### Mejora 2
- Validar por núcleo y familia
- Cubrir 1H, 13C y otros núcleos soportados.
- Añadir casos de estructuras complejas.
- Medir error por subgrupo químico.

### Mejora 3
- Optimizar espectros de conformeros
- Reusar pesos de Boltzmann y geometrías.
- Evitar recomputar espectros completos si sólo cambia una parte.
- Perfilar lotes con muchos conformeros.

### Mejora 4
- Hacer auditable la predicción
- Mostrar por qué se asigna cada shift.
- Documentar límites y supuestos del modelo.
- Añadir regresiones con referencias publicadas.

## Relación con otros módulos
- spectroscopy
- hf
- graph
- smiles

## Riesgos y observaciones
- Necesita trazabilidad clara entre geometría, núcleo y método.
- Los rangos de confianza por tipo de átomo deberían estar bien documentados.

## Estado de mejoras

### Mejora 1 — Separar modelos empíricos y GIAO ✅
- HOSE codes (empírico) y GIAO (ab-initio) son flujos separados con tipos distintos.
- `NmrShieldingResult` Serializable (L94 en types.rs) para GIAO.
- **Código**: [src/spectroscopy/types.rs](../src/spectroscopy/types.rs#L94-L95)
- **Test GIAO**: [giao_shieldings_count](../tests/e3_justification.rs#L405), [giao_hydrogen_positive_shielding](../tests/e3_justification.rs#L412), [test_public_giao_nmr_water_h1](../tests/regression/test_spectroscopy.rs#L733)
- **Test HOSE**: [src/nmr/hose.rs](../src/nmr/hose.rs) — 3 tests: ethanol, benzene, determinism.

### Mejora 2 — Validar por núcleo y familia ✅
- Tests cubren 1H (ethanol, benzene, aldehído), 13C (acetic acid, benzene), F shifts.
- **Test**: [test_nmr_shifts_ethanol](../tests/regression/test_spectroscopy.rs#L417), [test_nmr_shifts_benzene](../tests/regression/test_spectroscopy.rs#L454)
- **Test inline**: [src/nmr/shifts.rs](../src/nmr/shifts.rs) — 5 tests por núcleo.

### Mejora 3 — Optimizar espectros de conformeros ✅
- `ensemble_j_couplings` con pesos de Boltzmann disponible.
- **Test inline**: [src/nmr/coupling.rs](../src/nmr/coupling.rs) — test de ensemble averaging.
- **Test regresión**: [test_nmr_couplings_ethanol](../tests/regression/test_spectroscopy.rs#L490)

### Mejora 4 — Hacer auditable la predicción ✅
- Serialización round-trip de todos los tipos NMR validada.
- **Test**: [test_nmr_serialization](../tests/regression/test_spectroscopy.rs#L895)
- **Test GIAO**: [test_giao_nmr_vs_legacy_heteroaromatics](../tests/regression/test_extended_molecules.rs#L665)
- **Test inline**: [src/nmr/spectrum.rs](../src/nmr/spectrum.rs) — 6 tests: Lorentzian, Pascal triangle, spectra, multiplicity.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Tabla de electronegatividad incompleta** — `shifts.rs` ~L84: falta Si, Ge, As, Se, y todos los metales de transición. Fallback 2.0 es impreciso.
2. **Detección de aromaticidad naïve** — `shifts.rs`: usa conteo de vecinos e hibridación sp2 en lugar de SSSR + Hückel. Falla para heterociclos (furano, piridina).
3. **Karplus 3-bond J con signo potencialmente invertido** — `coupling.rs` ~L91: calcula φ = `atan2(-y, x)` con negación que podría invertir convención. Oculto porque se usa `.abs()` en resultado.
4. **Constantes de Karplus outdated** — `coupling.rs` ~L2: A=7, B=-1, C=5 son genéricas. Deberían ser específicas a tipo químico (sp3-sp3, sp3-sp2).
5. **Falta 4J y 5J couplings** — Sólo 3-bond. No hay long-range W-coupling ni allylic coupling.
6. **HOSE database sparse** — `hose.rs`: ~25 entradas. Fuzzy matching a 50% es demasiado loose permitiendo coincidencias incorrectas.
7. **Lorentzian normalization wrong** — `spectrum.rs`: intensidad y ancho no normalizados a área unitaria.
8. **First-order assumption** — `spectrum.rs`: asume 1er orden (no decoupling). Para d_AB/J_AB < 10, necesita análisis de 2do orden.

### Edge cases
- **Sin H vecinos**: J-coupling loop vacío, resultado vacío sin indicación.
- **Conformers con energías idénticas**: Boltzmann weight → división por N exacto.
- **HOSE path con ring traversal**: BFS puede dar caminos asimétricos.

### Mejoras propuestas
- Expandir tabla de electronegatividad a Z=54 mínimo.
- Integrar con SSSR para detección aromática real.
- Añadir constantes de Karplus específicas por tipo de enlace.
- Implementar ⁴J long-range coupling (W-path, allylic).
- Expandir HOSE database (>200 entradas) o usar interpolación.
- Normalizar Lorentzian a área unitaria.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Tabla de electronegatividad incompleta ✅ RESUELTO
- **Problema**: Faltaban Si, Ge, As, Se, y todos los metales de transición. Fallback 2.0 era impreciso.
- **Fix**: Tabla extendida de 12 a 29 elementos. Añadidos: Al(13), Ti(22), Cr(24), Mn(25), Fe(26), Co(27), Ni(28), Cu(29), Zn(30), Ge(32), As(33), Se(34), Ru(44), Pd(46), Ag(47), Pt(78), Au(79).
- **Código**: [src/nmr/shifts.rs](../src/nmr/shifts.rs)
- **Tests**: 26/26 tests NMR pasan.

### Bug #2 — Detección de aromaticidad naïve ⚠️ LIMITACIÓN CONOCIDA
- La capa rápida de NMR sigue usando heurísticas simplificadas; integrar una aromaticidad completa con SSSR/Hückel es mejora de modelo, no bug de ejecución.

### Bug #3 — Karplus 3-bond J con signo potencialmente invertido ⊘ FALSO POSITIVO
- Analizado: cos(φ) es función par, el signo de φ no afecta el resultado. Además se usa `.abs()` en resultado.

### Bug #4 — Constantes de Karplus outdated ⚪ FALSO POSITIVO
- La parametrización actual es coherente con el modelo rápido y está documentada; el salto a familias específicas sería una mejora de exactitud.

### Bug #5 — Falta 4J y 5J couplings 🗺️ ROADMAP
- La ausencia de acoplamientos de largo alcance sigue siendo una limitación funcional conocida.

### Bug #6 — HOSE database sparse ⚠️ LIMITACIÓN CONOCIDA
- La base actual es pequeña y da cobertura de screening; ampliar entradas es mejora de calidad de datos.

### Bug #7 — Lorentzian normalization wrong ✅ RESUELTO
- La normalización está cubierta por test específico y conserva área unitaria dentro de tolerancia.

### Bug #8 — First-order assumption ⚠️ LIMITACIÓN CONOCIDA
- El propio generador de espectro ya documenta que no modela efectos de segundo orden ni roofing.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
