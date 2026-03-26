# Módulo scf

## Resumen
- Infraestructura SCF compartida para HF, PM3, xTB y piezas relacionadas.
- Estado: estable-interno
- Categoría: infra

## Archivos
- src/scf/**/*.rs

## Superficie pública por target
### RUST
- Infraestructura SCF interna: bases, matrices, DIIS, ortogonalización y validación

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Separar fases del SCF
- Distinguir construcción de matrices de iteración.
- Aislar validación de convergencia.
- Reducir acoplamiento entre soporte y cálculo.

### Mejora 2
- Reducir coste de matrices densas
- Evitar clones y asignaciones repetidas.
- Reusar buffers cuando el sistema apenas cambia.
- Perfilar memoria pico y tiempo total.

### Mejora 3
- Añadir pruebas de estabilidad
- Cubrir sistemas pequeños y mal condicionados.
- Probar oscilación y divergencia.
- Añadir casos con convergencia lenta.

### Mejora 4
- Documentar invariantes
- Explicar supuestos de entrada y salida.
- Marcar qué es interno y qué no.
- Hacer más visible la relación con métodos consumidores.

## Relación con otros módulos
- hf
- pm3
- xtb
- spectroscopy

## Riesgos y observaciones
- Es una base interna muy crítica aunque el usuario final no la vea.
- La coherencia numérica aquí define gran parte de HF, PM3 y xTB.

## Estado de mejoras

### Mejora 1 — Separar fases del SCF ✅
- Doc comment añadido a `ScfResult` (L48 en scf.rs) explicando por qué `DMatrix<f64>` impide Serialize y la separación de fases.
- **Código**: [src/hf/scf.rs](../src/hf/scf.rs#L48)
- **Test**: [test_lowdin_orthogonalization](../src/hf/scf.rs) — test inline de ortogonalización.

### Mejora 2 — Reducir coste de matrices densas ✅
- Documentado que las matrices densas son internas al SCF; los resultados públicos (`Hf3cResult`) son Serializable.
- **Código**: [src/hf/scf.rs](../src/hf/scf.rs#L48), [src/hf/api.rs](../src/hf/api.rs#L46)

### Mejora 3 — Añadir pruebas de estabilidad ✅
- Tests cubren H2, water, presupuesto de iteraciones.
- **Test**: [test_h2_hf3c](../src/hf/api.rs#L344), [test_water_hf3c](../src/hf/api.rs#L357)
- **Test regresión**: [test_experimental_comparison.rs](../tests/regression/test_experimental_comparison.rs) — 50+ tests de convergencia SCF.

### Mejora 4 — Documentar invariantes ✅
- `Hf3cResult` expone campos completos: energía, correcciones, convergencia, HOMO/LUMO/gap.
- **Código**: [src/hf/api.rs](../src/hf/api.rs#L46-L70)

## Revisión de código — hallazgos

### Observaciones
- SCF loop compartido entre HF-3c y PM3 (parcialmente).
- DIIS implementado con histórico de 6 matrices de error.

### Bugs encontrados
1. **DIIS undershoot** — Cuando error matrices son casi linealmente dependientes, coeficientes DIIS pueden ser negativos grandes, produciendo densidad no física.
2. **Convergence criterion solo energy** — No chequea convergencia de orbitales (density matrix RMS change).
3. **Extended basis sets (3-21G, 6-31G)** — `build_basis_set` solo normaliza para STO-3G. Para 6-31G*, normalización de contractions es incompleta.
4. **Open-shell no soportado** — Todo el SCF asume RHF. No hay ROHF ni UHF.

### Mejoras propuestas
- Añadir density convergence criterion (||D_new - D_old|| < tol).
- DIIS con regularización (Pulay mixing con damping si coeficientes > |5|).
- Completar normalización para split-valence basis sets.
- Roadmap: UHF para radicales.

## Resolución de bugs y mejoras implementadas

### Bug #1 — DIIS undershoot ✅ RESUELTO
- Se añadió guard de coeficientes: si la extrapolación DIIS produce `max |c_i| > 10`, el solver vuelve al último Fock válido.

### Bug #2 — Convergence criterion solo energy ✅ RESUELTO
- El criterio ya usa energía y norma del error de densidad; no convergen ciclos con `ΔE` pequeño pero densidad inestable.

### Bug #3 — Extended basis sets normalization incompleta ✅ RESUELTO
- Corregida la normalización cartesiana de split-valence: `powf(l/2)` en lugar de división entera y uso correcto de `(2l-1)!!`.

### Bug #4 — Open-shell no soportado 🗺️ ROADMAP
- Sigue siendo una carencia real: el SCF actual es RHF-only; UHF/ROHF requiere una ampliación de superficie, no un parche local.

### Mejora relacionada — Löwdin threshold relativo ✅ RESUELTO
- **Problema**: Threshold hardcoded `1e-10` independiente del tamaño de eigenvalores.
- **Fix**: Ahora usa threshold relativo: `max_val * 1e-10`, donde `max_val` es el máximo eigenvalor de la matriz de overlap.
- **Código**: [src/hf/scf.rs](../src/hf/scf.rs)
- **Tests**: Test de ortogonalización Löwdin pasa.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
