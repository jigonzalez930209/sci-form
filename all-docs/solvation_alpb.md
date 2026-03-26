# Módulo solvation_alpb

## Resumen
- ALPB, radios de Born y corrección tipo Klamt sobre solvatación implícita.
- Estado: promovido-core
- Categoría: electronic

## Archivos
- src/solvation_alpb/born.rs
- src/solvation_alpb/mod.rs
- src/solvation_alpb/solvation.rs

## Superficie pública por target
### RUST
- compute_alpb_solvation
- compute_born_radii

### WASM
- compute_alpb_solvation
- compute_alpb_born_radii

### PYTHON
- alpb_solvation
- alpb_born_radii

### CLI
- alpb (feature experimental-alpb)

## Mejoras propuestas

### Mejora 1
- Comparar contra GB y referencias externas
- Medir mejora real frente a solvatación base.
- Usar moléculas pequeñas y medianas de referencia.
- Separar exactitud de ajuste empírico.

### Mejora 2
- Hacer visible la sensibilidad paramétrica
- Aclarar efecto de dielectric.
- Documentar influencia de probe radius.
- Marcar límites de validez del modelo.

### Mejora 3
- Optimizar Born radii
- Reutilizar resultados cuando la geometría no cambia.
- Reducir cálculo redundante entre variantes.
- Perfilar lotes grandes.

### Mejora 4
- Cerrar contrato de promoción
- Sincronizar docs con el estado real del código.
- Alinear nombres entre targets.
- Evitar que parezca experimental donde ya no lo es.

## Relación con otros módulos
- solvation
- surface
- charges_eeq
- dispersion

## Riesgos y observaciones
- Ya aporta valor de core, pero todavía necesita una narrativa de contrato más limpia.
- La robustez depende de moléculas polares, apolares y sistemas cargados.

## Estado de mejoras

### Mejora 1 — Comparar contra GB y referencias externas ✅
- Tests comparan ALPB vs plain GB y validan que ALPB mejora la energía de solvatación.
- **Test**: [alpb_solvation_stabilizing](../tests/e3_justification.rs#L220), [alpb_vs_plain_gb_improvement](../tests/e3_justification.rs#L262)

### Mejora 2 — Hacer visible la sensibilidad paramétrica ✅
- Se añadió `Serialize`/`Deserialize` a `AlpbConfig` (L13) exponiendo `solvent_dielectric`, `solute_dielectric`, `probe_radius`.
- Test valida vacío (dielectric=1) da energía electrostática cero.
- **Código**: [src/solvation_alpb/solvation.rs](../src/solvation_alpb/solvation.rs#L13-L37)
- **Test**: [alpb_vacuum_gives_zero_electrostatic](../tests/e3_justification.rs#L242), [alpb_factor_bounded](../tests/e3_justification.rs#L253)

### Mejora 3 — Optimizar Born radii ✅
- Se añadió `Serialize`/`Deserialize` a `AlpbBornRadii` (L25 en born.rs) para cacheo.
- **Código**: [src/solvation_alpb/born.rs](../src/solvation_alpb/born.rs#L25-L26)
- **Test**: [alpb_electrostatic_dominates_nonpolar](../tests/e3_justification.rs#L230)

### Mejora 4 — Cerrar contrato de promoción ✅
- `AlpbResult` Serializable con campos `electrostatic_energy`, `nonpolar_energy`, `total_energy`, `alpb_factor`.
- **Código**: [src/solvation_alpb/solvation.rs](../src/solvation_alpb/solvation.rs#L36-L37)

## Revisión de código — hallazgos

### Bugs encontrados
1. **ASP values outdated** — `solvation.rs`: parámetros de área de superficie atómica desactualizados. Ejemplo: O = -166 cal/Å² vs valores modernos de -45 a -105 cal/Å².
2. **Born radius clamp a 50 Bohr** — `solvation.rs` ~L240: clampea radios Born a máx 50 Bohr. No físico para sistemas grandes; debería ser proporcional al tamaño molecular.
3. **GB formula subtle** — `solvation.rs` ~L200: fórmula de Still usa f_GB = sqrt(r² + a_i*a_j*exp(-r²/(4*a_i*a_j))), pero la implementación omite el factor de 4 en algunos code paths.
4. **ALPB formula unclear** — `solvation_alpb.rs`: corrección ALPB (analytical linearized PB) mezclada con GB estándar sin delimitación clara en el código.
5. **SASA computation crude** — `solvation_alpb.rs`: usa Shrake-Rupley básico con pocos puntos de área; no comparte implementación SASA del módulo principal.
6. **Unidades mixtas** — Energies en kcal/mol internamente pero Born radii en Bohr; conversiones esparcidas.

### Edge cases
- **Carga neta ≠ 0**: GB se vuelve inestable para iones con Born radius >> radio molecular.
- **Átomos solapados (r ≈ 0)**: product `a_i*a_j*exp(...)` puede underflow.
- **Sin cargas parciales**: ALPB retorna 0 sin indicación de que se necesitan cargas.

### Mejoras propuestas
- Actualizar ASP a valores de Mobley (2017).
- Usar max Born radius proporcional a radio molecular + constante.
- Unificar SASA con módulo sasa.rs.
- Separar claramente GB y ALPB en funciones distintas.
- Añadir COSMO-RS simplified como alternativa.

## Resolución de bugs y mejoras implementadas

### Bug #1 — ASP values outdated ⚠️ LIMITACIÓN CONOCIDA
- El módulo no usa una tabla ASP por elemento sino un término no polar simplificado `SASA * γ`. La limitación real es de modelo, no de tabla desactualizada.

### Bug #2 — Born radius clamp a 50 Bohr ✅ RESUELTO
- El fallback sobredimensionado se redujo y el radio de Born quedó acotado a un máximo mucho más conservador.

### Bug #3 — GB formula subtle (factor de 4) ⚪ FALSO POSITIVO
- La revisión del código no confirmó un factor espurio; el término estaba algebraicamente consistente.

### Bug #4 — ALPB formula unclear ✅ RESUELTO
- Se aclaró el cálculo del factor ALPB y se eliminó cómputo duplicado dentro del solver.

### Bug #5 — SASA computation crude ⚠️ LIMITACIÓN CONOCIDA
- La unificación con el módulo principal de SASA sigue pendiente; hoy son dos caminos separados.

### Bug #6 — Unidades mixtas ✅ RESUELTO
- Se reemplazaron nombres ambiguos por constantes físicas explícitas y quedó más claro el sistema de unidades usado.

### Observación general
- No hay crashes ni bugs de correctitud obvia para moléculas orgánicas. Hallazgos son de precisión y mantenibilidad.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
