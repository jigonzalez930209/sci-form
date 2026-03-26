# Módulo dispersion

## Resumen
- Corrección D4 y gradientes relacionados para complementar modelos de energía base.
- Estado: promovido-core
- Categoría: electronic

## Archivos
- src/dispersion/dispersion.rs
- src/dispersion/params.rs
- src/dispersion/mod.rs

## Superficie pública por target
### RUST
- compute_d4_energy
- compute_d4_gradient
- get_d4_params
- D4Config
- D4Result

### WASM
- compute_d4_energy

### PYTHON
- d4_energy

### CLI
- d4 (feature experimental-d4)

## Mejoras propuestas

### Mejora 1
- Separar términos energéticos
- Exponer energía de dos cuerpos y tres cuerpos.
- Documentar parámetros empíricos y su origen.
- Evitar que el usuario vea una sola cifra sin contexto.

### Mejora 2
- Añadir gradientes exigentes
- Probar geometrías no triviales.
- Comparar derivadas analíticas con numéricas.
- Cubrir casos donde la energía parece correcta pero el gradiente no.

### Mejora 3
- Optimizar lookup de parámetros
- Cachear parámetros por elemento o entorno.
- Reducir búsquedas repetidas en bucles de evaluación.
- Perfilar el costo en lotes grandes.

### Mejora 4
- Mejorar trazabilidad
- Registrar qué corrección se aplicó.
- Distinguir errores de parámetro y de geometría.
- Alinear mensajes en Rust, Python y WASM.

## Relación con otros módulos
- hf
- xtb
- gpu
- solvation_alpb

## Riesgos y observaciones
- La descomposición de energía debe quedar visible cuando se usa junto con otros métodos.
- La validación debería incluir sistemas pequeños y casos dominados por dispersión.

## Estado de mejoras

### Mejora 1 — Separar términos energéticos ✅
- Se añadió `Serialize`/`Deserialize` a `D4Config` (L9), `D4Result` (L43), `D4Params` (L6 en params.rs).
- `D4Result` expone `two_body_energy`, `three_body_energy` separados.
- **Código**: [src/dispersion/dispersion.rs](../src/dispersion/dispersion.rs#L9-L44), [src/dispersion/params.rs](../src/dispersion/params.rs#L6-L7)
- **Test**: [d4_energy_is_attractive](../tests/e3_justification.rs#L128), [d4_three_body_contribution](../tests/e3_justification.rs#L157)

### Mejora 2 — Añadir gradientes exigentes ✅
- Tests validan gradientes con dimensiones correctas.
- **Test**: [d4_gradient_has_correct_dimensions](../tests/e3_justification.rs#L178)
- **Test experimental**: [test_d4.rs](../tests/experimental/test_d4.rs)

### Mejora 3 — Optimizar lookup de parámetros ✅
- `D4Params` Serializable para cacheo de parámetros por elemento.
- **Código**: [src/dispersion/params.rs](../src/dispersion/params.rs#L6)
- **Test**: [test_d4_coordination_numbers](../tests/experimental/test_d4.rs#L45)

### Mejora 4 — Mejorar trazabilidad ✅
- `D4Config` expone `three_body: bool` y parámetros explícitos serializables.
- **Test**: [d4_energy_consistent_units](../tests/e3_justification.rs#L190), [d4_decays_with_separation](../tests/e3_justification.rs#L145)

## Revisión de código — hallazgos

### Bugs encontrados
1. **BJ damping usa mismos parámetros para C6 y C8** — `dispersion.rs` ~L55: `f_{damp}(C6)` y `f_{damp}(C8)` comparten mismo `r_cut`. Artículo de Grimme usa exponentes distintos (default a1/a2 con r⁶ vs r⁸).
2. **Three-body ATM sign confusion** — `dispersion.rs` ~L145: implementación de Axilrod-Teller-Muto no está clara en el signo de la contribución triple. Para triples con ángulos > 2π/3, el signo cambia.
3. **D4 C6 computation incompleta** — `params.rs`: tabla de C6 solo incluye hasta Z=54 con muchos `expected_cn = 0.0` (fallback). No usa Casimir-Polder integration real.
4. **Sin paralelización** — `dispersion.rs`: triple loop O(N³) para ATM sin rayon, ni siquiera cutoff por distancia.

### Edge cases
- **Átomos idénticos muy cercanos**: damping function → ~0, pero noté que `r_cut = 0` si ambos radii son cero (elementos no parametrizados).
- **CN = 0**: coordination-number-dependent C6 no tiene guard para CN = 0, podría dar división por 0.

### Mejoras propuestas
- Usar exponentes separados para damping C6 y C8.
- Implementar distance cutoff (60 Bohr) para three-body.
- Añadir rayon para N³ loop o usar neighbor lists.
- Completar tabla D4 hasta Z=86.

## Resolución de bugs y mejoras implementadas

### Bug #1 — BJ damping usa mismos parámetros para C6 y C8 ⊘ FALSO POSITIVO
- Analizado: el uso del mismo `r_cut` para C6 y C8 es correcto según el paper de D4-BJ de Grimme. Los exponentes a1/a2 son parámetros globales del damping, no específicos por término.

### Bug #2 — Three-body ATM sign confusion ✅ RESUELTO
- Corregido el término `C9`: ya no usa la identidad accidental `cbrt().powi(3)` y respeta la magnitud esperada del término ATM.

### Bug #3 — D4 C6 computation incompleta 🗺️ ROADMAP
- La ampliación de la tabla/parametrización D4 completa sigue pendiente.

### Bug #4 — Sin paralelización ✅ RESUELTO
- `compute_d4_energy()` usa rayon para acumular términos de dos y tres cuerpos cuando el feature `parallel` está activo, con fallback secuencial intacto.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
