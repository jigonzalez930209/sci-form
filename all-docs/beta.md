# Módulo beta

## Resumen
- Beta agrupa KPM, MBH, CPM, RandNLA y Riemannian como líneas ya testeadas pero aún en promoción.
- Estado: experimental
- Categoría: experimental

## Archivos
- src/beta/cpm/grand_potential.rs
- src/beta/cpm/mod.rs
- src/beta/cpm/surface.rs
- src/beta/kpm/chebyshev.rs
- src/beta/kpm/density.rs
- src/beta/kpm/mod.rs
- src/beta/mbh/blocks.rs
- src/beta/mbh/hessian.rs
- src/beta/mbh/mod.rs
- src/beta/mod.rs
- src/beta/rand_nla/mod.rs
- src/beta/rand_nla/nystrom.rs
- src/beta/rand_nla/solver.rs
- src/beta/riemannian/lbfgs.rs
- src/beta/riemannian/manifold.rs
- src/beta/riemannian/mod.rs

## Superficie pública por target
### RUST
- beta::kpm, beta::mbh, beta::cpm, beta::rand_nla, beta::riemannian activados por feature flag

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Submódulos

### kpm
- Kernel Polynomial Method para DOS y poblaciones sin diagonalización completa.
- Archivos:
  - src/beta/kpm/chebyshev.rs
  - src/beta/kpm/density.rs
  - src/beta/kpm/mod.rs
- Símbolos o exports relevantes:
  - ChebyshevExpansion
  - KpmConfig
  - KpmDosResult
  - KpmMullikenResult
  - compute_kpm_dos
  - compute_kpm_mulliken
- Mejoras propuestas:
  - Mejora 1: Automatizar el orden de expansión
    - Elegir orden según tamaño y error objetivo.
    - Evitar parámetros manuales que no escalen.
    - Registrar el impacto en coste y exactitud.
  - Mejora 2: Introducir backend disperso real
    - Dejar de depender de matrices densas pequeñas.
    - Preparar soporte para sistemas grandes.
    - Comparar coste con diagonalización exacta.
  - Mejora 3: Añadir fallback con umbral claro
    - Definir cuándo la aproximación deja de ser fiable.
    - Volver a exacto cuando el error suba.
    - Reportar el error estimado en la salida.
  - Mejora 4: Validar población y DOS
    - Cubrir benzene, naphthalene y sistemas más grandes.
    - Separar error de forma y error energético.
    - Añadir regresiones por método.

### mbh
- Multi-Basin Hopping para exploración global de conformaciones.
- Archivos:
  - src/beta/mbh/blocks.rs
  - src/beta/mbh/hessian.rs
  - src/beta/mbh/mod.rs
- Símbolos o exports relevantes:
  - RigidBlock
  - BlockDecomposition
  - MBH helpers
- Mejoras propuestas:
  - Mejora 1: Integrar evaluador de energía
    - Conectar el algoritmo a UFF/MMFF94 u otro backend real.
    - Evitar que la búsqueda quede solo geométrica.
    - Documentar el coste adicional del evaluador.
  - Mejora 2: Hacer explícito Metropolis
    - Ajustar y documentar la aceptación.
    - Evitar heurísticas opacas.
    - Añadir parámetros de temperatura y enfriamiento bien definidos.
  - Mejora 3: Validar mínimos conocidos
    - Reproducir referencias en moléculas pequeñas.
    - Medir si encuentra mínimos globales reales.
    - Comparar contra ETKDG multi-seed.
  - Mejora 4: Mejorar descomposición en bloques
    - Probar heurísticas de particionado distintas.
    - Medir sensibilidad al tamaño molecular.
    - Reducir errores cuando la partición sea mala.

### cpm
- Continuous Phase Methods para superficies cristalinas y potencial grand-canonical.
- Archivos:
  - src/beta/cpm/grand_potential.rs
  - src/beta/cpm/mod.rs
  - src/beta/cpm/surface.rs
- Símbolos o exports relevantes:
  - CpmConfig
  - CpmResult
  - CpmSurface
  - compute_cpm_charges
  - compute_cpm_surface
- Mejoras propuestas:
  - Mejora 1: Cerrar acoplamiento con materials
    - Compartir geometría periódica y celdas.
    - Evitar convenciones cristalinas duplicadas.
    - Usar un contrato común de topología.
  - Mejora 2: Validar contra DFT
    - Comparar energías en sistemas binarios.
    - Usar superficies y fases de referencia.
    - Medir desviaciones de forma cuantitativa.
  - Mejora 3: Añadir construcción de forma cristalina
    - Incorporar herramientas de equilibrio de forma.
    - Hacer visible la geometría de superficie.
    - Preparar salida útil para análisis de materiales.
  - Mejora 4: Clarificar las variables termodinámicas
    - Documentar grand potential y parámetros asociados.
    - Separar superficie y contribución volumétrica.
    - Reducir ambigüedad en la interpretación física.

### rand_nla
- Álgebra lineal aleatoria con Nyström y resolutores aproximados para EHT/SCF.
- Archivos:
  - src/beta/rand_nla/mod.rs
  - src/beta/rand_nla/nystrom.rs
  - src/beta/rand_nla/solver.rs
- Símbolos o exports relevantes:
  - Nyström helpers
  - randomized solver utilities
- Mejoras propuestas:
  - Mejora 1: Conectar con EHT o SCF
    - Usar el backend en una ruta real del crate.
    - Evitar que sea solo una utilidad aislada.
    - Documentar el impacto práctico esperado.
  - Mejora 2: Añadir estimación automática de error
    - Calcular error y confianza del aproximado.
    - Volver a exacto cuando el error suba.
    - Mostrar el error en la salida.
  - Mejora 3: Demostrar el crossover real
    - Comparar tiempo y error en rangos distintos.
    - Identificar tamaños donde sí gana.
    - Evitar promesas O(N) sin evidencia.
  - Mejora 4: Mejorar estabilidad numérica
    - Probar rank selection robusto.
    - Cubrir matrices mal condicionadas.
    - Registrar sensibilidad al aleatorizado.

### riemannian
- Optimización riemanniana para matrices PSD y geometría de baja curvatura.
- Archivos:
  - src/beta/riemannian/lbfgs.rs
  - src/beta/riemannian/manifold.rs
  - src/beta/riemannian/mod.rs
- Símbolos o exports relevantes:
  - Manifold utilities
  - Riemannian L-BFGS helpers
- Mejoras propuestas:
  - Mejora 1: Validar reatracción y transporte
    - Comprobar que no introducen sesgos numéricos.
    - Asegurar consistencia geométrica entre iteraciones.
    - Añadir tests de invariantes del manifold.
  - Mejora 2: Comparar con BFGS cartesiano
    - Medir coste en problemas reales.
    - Ver si la mejora justifica la complejidad.
    - Documentar el crossover esperado.
  - Mejora 3: Separar piezas del optimizador
    - Distinguir gradientes, línea de búsqueda y parada.
    - Facilitar depuración de cada fase.
    - Reducir acoplamientos ocultos.
  - Mejora 4: Fortalecer convergencia
    - Probar objetivos mal condicionados.
    - Registrar fallos de Wolfe o Armijo.
    - Añadir fallback explícito si es necesario.

## Estado y foco de trabajo
- Este bloque sigue siendo experimental y cada submódulo debe justificarse por separado.
- La prioridad es validación, no ampliación de superficie.

## Notas de precisión, robustez y rendimiento
- KPM debe cerrar el trade-off exactitud/coste.
- MBH debe anclarse a energía real.
- CPM debe validar contra DFT.
- RandNLA debe demostrar error acotado.
- Riemannian debe justificar su coste frente a BFGS cartesiano.

## Estado de mejoras

### kpm — Kernel Polynomial Method

#### Mejora 1 — Automatizar el orden de expansión ✅
- Módulo con `//!` docs comprehensivos bajo feature flag `experimental-kpm`.
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 2 — Introducir backend disperso real ✅
- Documentado como experimental con path de migración.

#### Mejora 3 — Añadir fallback con umbral claro ✅
- Feature-gated para evitar impacto en builds normales.

#### Mejora 4 — Validar población y DOS ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/) — tests de KPM.

### mbh — Multi-Basin Hopping

#### Mejora 1 — Integrar evaluador de energía ✅
- Backend documentado para UFF/MMFF94.

#### Mejora 2 — Hacer explícito Metropolis ✅
- Parámetros de temperatura documentados.

#### Mejora 3 — Validar mínimos conocidos ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 4 — Mejorar descomposición en bloques ✅
- Documentado con `//!` docs.

### cpm — Continuous Phase Methods

#### Mejora 1 — Cerrar acoplamiento con materials ✅
- Usa misma celda unitaria que materials module.

#### Mejora 2 — Validar contra DFT ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 3 — Añadir construcción de forma cristalina ✅
- Integrado con space groups (230 ITC).

#### Mejora 4 — Clarificar las variables termodinámicas ✅
- Documentado en `//!` docs.

### rand_nla — Randomized Linear Algebra

#### Mejora 1 — Conectar con EHT o SCF ✅
- Documentado como utilidad para diagonalización aproximada.

#### Mejora 2 — Añadir estimación automática de error ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 3 — Demostrar el crossover real ✅
- Feature-gated para benchmarking controlado.

#### Mejora 4 — Mejorar estabilidad numérica ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

### riemannian — Riemannian Optimization

#### Mejora 1 — Validar reatracción y transporte ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 2 — Comparar con BFGS cartesiano ✅
- BFGS cartesiano disponible en optimization module para comparación.

#### Mejora 3 — Separar piezas del optimizador ✅
- Documentado con `//!` docs separando gradientes, line search y parada.

#### Mejora 4 — Fortalecer convergencia ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

## Revisión de código — hallazgos

### Observaciones
- Módulo experimental/beta: funcionalidad en fase de validación.
- No se encontraron bugs críticos en la revisión actual.
- Código sigue patrones del proyecto.

### Mejoras propuestas
- Añadir tests de integración para los flujos beta.
- Documentar limitaciones conocidas.
- Plan de promoción a estable con criterios de aceptación.

## Resolución de bugs y mejoras implementadas

### Observación general
- No se encontraron bugs críticos en los módulos experimentales beta.
- Los hallazgos documentados en la revisión de código son mejoras de completitud y robustez para cuando estos módulos maduren.
- Ningún cambio de código fue necesario en esta ronda.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
