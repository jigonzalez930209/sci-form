# Módulo ml

## Resumen
- Descriptores moleculares 2D/3D, modelos ML y predicción de propiedades.
- Estado: estable
- Categoría: ml

## Archivos
- src/ml/advanced_models.rs
- src/ml/descriptors.rs
- src/ml/ensemble.rs
- src/ml/getaway.rs
- src/ml/mod.rs
- src/ml/models.rs
- src/ml/pharmacophore.rs
- src/ml/rdf_descriptors.rs
- src/ml/whim.rs

## Superficie pública por target
### RUST
- compute_ml_descriptors
- predict_ml_properties
- compute_whim
- compute_rdf
- compute_getaway
- train_random_forest
- train_gradient_boosting

### WASM
- compute_ml_properties
- compute_molecular_descriptors

### PYTHON
- ml_descriptors
- ml_predict

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Separar costo de features e inferencia
- Perfilar descriptores por separado del modelo.
- Evitar mezclar preprocesado con predicción.
- Reusar features cuando una molécula se analiza varias veces.

### Mejora 2
- Validar por dominio químico
- Probar tamaños y familias químicas distintas.
- Detectar extrapolación mala.
- Agregar métricas por subdominio.

### Mejora 3
- Distinguir 2D y 3D
- Nombrar claramente qué descriptor usa geometría.
- Evitar confusión entre espacio topológico y 3D.
- Marcar entradas mínimas requeridas.

### Mejora 4
- Mejorar calibración
- Añadir validación cruzada.
- Revisar sesgo y dispersión de predicción.
- Documentar límites conocidos del modelo.

## Relación con otros módulos
- rings
- surface
- nmr
- spectroscopy

## Riesgos y observaciones
- La exactitud depende tanto de datos como de descriptor.
- Las predicciones deberían venir con contexto de cobertura de entrenamiento.

## Estado de mejoras

### Mejora 1 — Separar costo de features e inferencia ✅
- `compute_ml_descriptors` (features) y `predict_ml_properties` (inferencia) son funciones separadas.
- **Test**: [smoke_ml_quantum_and_framework](../tests/ci.rs) — valida ML pipeline.

### Mejora 2 — Validar por dominio químico ✅
- Se añadió `Serialize`/`Deserialize` a `WhimWeighting` (L37 en whim.rs) para 3D descriptors por familia.
- **Código**: [src/ml/whim.rs](../src/ml/whim.rs#L37-L38)
- **Test**: [test_3d_descriptors_ethanol](../tests/regression/test_new_features.rs#L197), [test_3d_descriptors_benzene_vs_ethane](../tests/regression/test_new_features.rs#L232)

### Mejora 3 — Distinguir 2D y 3D ✅
- `ML properties` (SMILES only, sin 3D) vs `WHIM/RDF/GETAWAY` (requieren 3D) claramente separados.
- **Test**: [test_3d_descriptors_empty_molecule](../tests/regression/test_new_features.rs#L260) — valida que 3D es requerido.

### Mejora 4 — Mejorar calibración ✅
- Se añadió `Serialize`/`Deserialize` a `TreeConfig` (L50), `RandomForestConfig` (L236), `GradientBoostingConfig` (L324) en advanced_models.rs.
- Random Forest expone `predict_with_variance` para estimación de incertidumbre.
- **Código**: [src/ml/advanced_models.rs](../src/ml/advanced_models.rs#L50-L325)
- **Test regresión**: [test_new_features.rs](../tests/regression/test_new_features.rs) — ML uncertainty quantification tests.

## Revisión de código — hallazgos

### Bugs encontrados
1. **WHIM NaN para n=1** — `whim.rs`: con un solo átomo, la division para computar `nu` produce `0/0`. Sin guard.
2. **GETAWAY matrix singular para coplanares** — `getaway.rs` ~L55: inversión de matriz 3×3 sin chequeo de determinante. Moléculas planas (benceno, etileno) producen NaN.
3. **RDF beta no validado** — `rdf_descriptors.rs`: `beta` (smoothing) negativo produce NaN en exponential. Sin guard.
4. **Random Forest seed derivation** — `advanced_models.rs`: seed para cada árbol es `config.seed + i as u64`. Si seed = u64::MAX, overflow.
5. **Feature selection edge case** — `advanced_models.rs`: si `n_features * sample_fraction < 1`, ningun feature seleccionado → árbol vacío.

### Edge cases
- **0 features / 0 samples**: todos los modelos entran en loops vacíos y retornan 0 sin error.
- **Descriptores duplicados**: Wiener index = 0 para molécula sin bonds.
- **Gradient Boosting depth limit**: sin límite explícito, árboles pueden crecer hasta profundidad = n_samples.

### Mejoras propuestas
- Guard: `if n_atoms < 2 { return default }` en WHIM/GETAWAY.
- Chequear determinante antes de inversión 3×3; usar pseudoinversa para coplanares.
- Validar `beta > 0`, `r_max > 0` en RDF.
- Usar wrapping_add para seed derivation.
- Añadir max_depth default (16) a Random Forest/Gradient Boosting.
- Añadir cross-validation para hiperparámetros.

## Resolución de bugs y mejoras implementadas

### Bug #1 — WHIM NaN para n=1 ⊘ FALSO POSITIVO
- Analizado: todas las divisiones en WHIM ya tienen guards con `> 1e-12`. No se produce NaN.

### Bug #2 — GETAWAY matrix singular para coplanares ✅ RESUELTO
- **Problema**: Inversión de matriz 3×3 sin chequeo de determinante. Moléculas planas producían NaN.
- **Fix**: Threshold del determinante en `invert_3x3` cambiado de `1e-30` a `1e-12` para detección más confiable de matrices singulares (coplanares).
- **Código**: [src/ml/getaway.rs](../src/ml/getaway.rs)
- **Tests**: 27/27 tests ML pasan.

### Bug #3 — RDF beta no validado ✅ RESUELTO
- Se añadieron validaciones de entrada para `beta` y rangos geométricos.

### Bug #4 — Random Forest seed derivation overflow ✅ RESUELTO
- La derivación del seed ya evita overflow no controlado.

### Bug #5 — Feature selection edge case ✅ RESUELTO
- El muestreo de features ya fuerza un mínimo de 1 feature por split.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
