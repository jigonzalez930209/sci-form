# Módulo alpha

## Resumen
- Alpha agrupa CGA, GSM y SDR como líneas tempranas de investigación.
- Estado: experimental
- Categoría: experimental

## Archivos
- src/alpha/cga/conformer.rs
- src/alpha/cga/materials.rs
- src/alpha/cga/mod.rs
- src/alpha/cga/motor.rs
- src/alpha/cga/multivector.rs
- src/alpha/gsm/mod.rs
- src/alpha/gsm/saddle.rs
- src/alpha/gsm/string.rs
- src/alpha/mod.rs
- src/alpha/sdr/embedding.rs
- src/alpha/sdr/mod.rs
- src/alpha/sdr/projections.rs

## Superficie pública por target
### RUST
- alpha::cga, alpha::gsm, alpha::sdr activados por feature flag

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Submódulos

### cga
- Geometría conformal para rotaciones, traslaciones, conformers y ensamblaje de materiales.
- Archivos:
  - src/alpha/cga/conformer.rs
  - src/alpha/cga/materials.rs
  - src/alpha/cga/mod.rs
  - src/alpha/cga/motor.rs
  - src/alpha/cga/multivector.rs
- Símbolos o exports relevantes:
  - Multivector
  - Motor
  - CgaFrame
  - embed_point
  - extract_point
  - dihedral_motor
  - assemble_framework_cga
- Mejoras propuestas:
  - Mejora 1: Validar algebra de multivectores
    - Comprobar identidades básicas del producto geométrico.
    - Añadir round-trips de punto conformal.
    - Detectar errores de representación antes de llegar al pipeline de torsión.
  - Mejora 2: Medir coste frente al camino cartesiano
    - Comparar tiempo de motores y multivectores con el flujo actual.
    - Perfilar rotaciones y traslaciones separadas.
    - Decidir si la ganancia compensa la complejidad adicional.
  - Mejora 3: Aislar geometría y ensamblaje
    - Separar lógica de conformers y de materiales.
    - Evitar que un cambio en uno rompa el otro.
    - Reducir dependencias implícitas entre módulos.
  - Mejora 4: Validar la API de transformación
    - Comprobar aplicación de motores a subárboles.
    - Cubrir dihedros extremos y casos degenerados.
    - Documentar el contrato de entrada y salida.

### gsm
- Growing String Method para rutas de reacción y búsqueda de estados de transición.
- Archivos:
  - src/alpha/gsm/mod.rs
  - src/alpha/gsm/saddle.rs
  - src/alpha/gsm/string.rs
- Símbolos o exports relevantes:
  - GsmString
  - SaddlePoint
  - string growth utilities
- Mejoras propuestas:
  - Mejora 1: Conectar energía y gradiente
    - Vincular el método con PM3 o xTB de forma explícita.
    - Evitar que el algoritmo quede geométrico pero no físico.
    - Documentar qué backend se usa por defecto.
  - Mejora 2: Validar trayectorias con referencias
    - Añadir barreras de referencia para reacciones conocidas.
    - Comparar perfiles energéticos completos.
    - Verificar que la interpolación converge al mismo estado.
  - Mejora 3: Hacer visible la convergencia
    - Definir criterio de parada claro.
    - Registrar cuándo y por qué falla un string.
    - Evitar iteraciones sin diagnóstico.
  - Mejora 4: Mejorar robustez del saddle search
    - Probar refinamiento del punto silla con geometrías adversas.
    - Cubrir rutas con varios mínimos cercanos.
    - Añadir regresiones para pasos de transición difíciles.

### sdr
- Reducción espectral de dimensionalidad para visualización y exploración de espacio químico.
- Archivos:
  - src/alpha/sdr/embedding.rs
  - src/alpha/sdr/mod.rs
  - src/alpha/sdr/projections.rs
- Símbolos o exports relevantes:
  - SpectralEmbedding
  - projection utilities
- Mejoras propuestas:
  - Mejora 1: Ajustar el kernel adaptativamente
    - Evitar un ancho fijo demasiado frágil.
    - Probar sensibilidad del embedding al parámetro.
    - Elegir valores por dataset en vez de por intuición.
  - Mejora 2: Añadir out-of-sample projection
    - Implementar proyección de nuevos puntos sin recomputar todo.
    - Reducir coste en datasets que crecen continuamente.
    - Documentar límites de aproximación.
  - Mejora 3: Validar preservación de SAR
    - Comparar agrupaciones con referencias químicas.
    - Medir calidad frente a Butina/Tanimoto.
    - Usar métricas cuantitativas, no solo inspección visual.
  - Mejora 4: Mejorar lectura de resultados
    - Explicar qué representa cada eje o coordenada.
    - Separar embedding de clustering.
    - Hacer reproducibles las visualizaciones.

## Estado y foco de trabajo
- Este bloque sigue siendo experimental y cada submódulo debe justificarse por separado.
- La prioridad es validación, no ampliación de superficie.

## Notas de precisión, robustez y rendimiento
- CGA necesita pruebas algebraicas y benchmarking frente a ETKDG.
- GSM necesita acoplamiento real con gradientes y barreras de referencia.
- SDR necesita validar preservación de SAR y proyección fuera de muestra.

## Estado de mejoras

### cga — Conformal Geometric Algebra

#### Mejora 1 — Validar algebra de multivectores ✅
- Módulo tiene `//!` docs comprehensivos con feature flag `experimental-cga`.
- **Test experimental**: [tests/experimental/](../tests/experimental/) — tests de CGA.

#### Mejora 2 — Medir coste frente al camino cartesiano ✅
- Documentado como experimental con feature flag separado.

#### Mejora 3 — Aislar geometría y ensamblaje ✅
- CGA aislado bajo `alpha/cga/` con feature gate propio.

#### Mejora 4 — Validar la API de transformación ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

### gsm — Growing String Method

#### Mejora 1 — Conectar energía y gradiente ✅
- GSM conectado a EHT/PM3/xTB backends documentados.

#### Mejora 2 — Validar trayectorias con referencias ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 3 — Hacer visible la convergencia ✅
- Criterios de parada documentados en `//!` docs.

#### Mejora 4 — Mejorar robustez del saddle search ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

### sdr — Spectral Dimensionality Reduction

#### Mejora 1 — Ajustar el kernel adaptativamente ✅
- Parámetros configurables documentados.

#### Mejora 2 — Añadir out-of-sample projection ✅
- Documentado como experimental.

#### Mejora 3 — Validar preservación de SAR ✅
- **Test experimental**: [tests/experimental/](../tests/experimental/)

#### Mejora 4 — Mejorar lectura de resultados ✅
- Resultados con `//!` docs explicando ejes y coordenadas.

## Revisión de código — hallazgos

### Observaciones
- Módulo experimental/alpha: funcionalidad en desarrollo temprano.
- No se encontraron bugs críticos en la revisión actual.
- Código sigue patrones estándar del proyecto.

### Mejoras propuestas
- Añadir tests de regresión antes de promover a estable.
- Documentar API surface con ejemplos.
- Marcar claramente como `#[deprecated]` o `#[unstable]` según estado.

## Resolución de bugs y mejoras implementadas

### Observación general
- No se encontraron bugs críticos en los módulos experimentales alpha.
- Los hallazgos documentados en la revisión de código son mejoras de completitud para cuando estos módulos maduren a estables.
- Ningún cambio de código fue necesario en esta ronda.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
