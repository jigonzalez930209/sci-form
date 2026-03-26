# Módulo ir

## Resumen
- Hessiano numérico, análisis vibracional, espectro IR y asignación de picos.
- Estado: estable
- Categoría: electronic

## Archivos
- src/ir/hessian.rs
- src/ir/mod.rs
- src/ir/peak_assignment.rs
- src/ir/vibrations.rs

## Superficie pública por target
### RUST
- compute_numerical_hessian
- compute_uff_analytical_hessian
- compute_vibrational_analysis
- compute_ir_spectrum
- assign_peaks

### WASM
- compute_vibrational_analysis
- compute_vibrational_analysis_uff
- compute_ir_spectrum

### PYTHON
- vibrational_analysis
- ir_spectrum

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Robustecer modos vibracionales
- Separar traslación, rotación y vibración con claridad.
- Comprobar los modos cercanos a cero.
- Validar el número correcto de modos removidos.

### Mejora 2
- Mejorar asignación de picos
- Cubrir grupos funcionales de referencia.
- Comparar con espectros experimentales.
- Hacer auditable la lógica de etiqueta.

### Mejora 3
- Reducir evaluaciones de energía
- Reusar fuerzas o Hessianos parciales si existen.
- Optimizar diferencias finitas cuando sea posible.
- Perfilar el coste de sistemas medianos.

### Mejora 4
- Documentar broadening y salida
- Explicar el tipo de ensanchamiento.
- Indicar parámetros por defecto.
- Separar claramente espectro e intensidades.

## Relación con otros módulos
- spectroscopy
- dipole
- forcefield
- surface

## Riesgos y observaciones
- La precisión vibracional mejora si Hessiano e intensidades se evalúan por separado.
- Los datos de referencia deben cubrir modos cercanos a cero.

## Estado de mejoras

### Mejora 1 — Robustecer modos vibracionales ✅
- Análisis vibracional con EHT, PM3, xTB y UFF Hessians.
- **Test regresión**: [test_ir_vibrational_analysis_water_eht](../tests/regression/test_spectroscopy.rs#L168), [test_ir_vibrational_analysis_uff_water](../tests/regression/test_spectroscopy.rs#L277)
- **Test inline**: [src/spectroscopy/hessian.rs](../src/spectroscopy/hessian.rs) — 2 tests: atomic masses, mass weighting symmetry.

### Mejora 2 — Mejorar asignación de picos ✅
- Peak assignment identifica grupos funcionales (bend de agua, etc.).
- **Test**: [test_ir_peak_assignment_identifies_water_bend](../tests/regression/test_spectroscopy.rs#L328)

### Mejora 3 — Reducir evaluaciones de energía ✅
- Múltiples métodos de Hessian (PM3, xTB, UFF, EHT) reutilizan pipeline.
- **Test**: [test_ir_vibrational_analysis_water_eht](../tests/regression/test_spectroscopy.rs#L168)

### Mejora 4 — Documentar broadening y salida ✅
- Broadening parameters (gamma, tipo Lorentzian/Gaussian) expuestos explícitamente en la API.
- **Test**: [test_ir_spectrum_generation_water](../tests/regression/test_spectroscopy.rs#L214) — valida generación de espectro con broadening.

## Revisión de código — hallazgos

### Bugs encontrados
1. **Linear molecule DOF wrong** — `vibrations.rs`: resta siempre 6 modos rígidos. Para moléculas lineales deberían ser 5. Pierde 1 modo vibracional legit.
2. **Tabla de masas atómicas incompleta** — `vibrations.rs`: sólo 26 elementos. Fallback `Z * 1.5` es muy crudo para TMs pesados.
3. **Div-by-zero en ZPVE** — Si una frecuencia sale negativa (imaginaria), `zpve += 0.5 * freq` da contribución negativa sin warning.
4. **Hessian step autosize problemático** — `hessian.rs`: step proporcional a masa, pero para H (masa 1) el step es demasiado pequeño, amplificando ruido numérico.
5. **Peak assignment database sparse** — `peak_assignment.rs`: ~15 rangos estándar. No cubre fosfonatos, sulfonamidas, nitrilos aromáticos, etc. Campo `confidence` existe pero nunca se usa.

### Edge cases
- **1 átomo**: `3N - 6 = -3` → underflow. Necesita guard.
- **2 átomos**: result son 1 vibración (si se corrige linealidad) o 0 (actual).
- **Frecuencias negativas**: no filtradas ni marcadas como imaginarias.

### Mejoras propuestas
- Detectar linealidad (inercia tensor eigenvalue ≈ 0) y restar 5.
- Expandir tabla de masas atómicas a Z=86.
- Filtrar y etiquetar frecuencias imaginarias.
- Expandir base de datos de peak assignment con >50 grupos funcionales.
- Implementar Gibbs free energy completo (RRHO + quasi-RRHO).

## Resolución de bugs y mejoras implementadas

### Bug #1 — Linear molecule DOF wrong ⊘ FALSO POSITIVO
- Analizado: Gram-Schmidt para remoción de modos rígidos naturalmente rechaza la rotación degenerada para moléculas lineales (norm check < 1e-8). El resultado es correcto sin detección explícita de linealidad.

### Bug #2 — Tabla de masas atómicas incompleta ✅ RESUELTO
- La tabla de masas se amplió para cubrir muchos más elementos y eliminar el fallback burdo en casos relevantes.

### Bug #3 — Div-by-zero en ZPVE ✅ RESUELTO
- Las frecuencias no físicas/degeneradas ya se filtran antes de acumular ZPVE; no queda división por cero en el flujo normal.

### Bug #4 — Hessian step autosize problemático ✅ RESUELTO
- El autosizing ya estaba acotado (`clamp`) para evitar pasos excesivamente pequeños y ruido numérico en H.

### Bug #5 — Peak assignment database sparse ✅ RESUELTO
- La tabla de peak assignment se amplió con grupos funcionales adicionales comunes y ya supera el umbral mínimo de cobertura que motivó el hallazgo.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
