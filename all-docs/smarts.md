# Módulo smarts

## Resumen
- Parser y matcher SMARTS para subestructura y torsiones experimentales.
- Estado: estable-interno
- Categoría: infra

## Archivos
- src/smarts/matcher.rs
- src/smarts/mod.rs
- src/smarts/parser.rs
- src/smarts/torsion_data.rs
- src/smarts/torsion_matcher.rs

## Superficie pública por target
### RUST
- parse_smarts
- substruct_match
- has_substruct_match
- match_experimental_torsions
- SmartsPattern

### WASM
- Sin exportaciones directas.

### PYTHON
- Sin exportaciones directas.

### CLI
- Sin exportaciones directas.

## Mejoras propuestas

### Mejora 1
- Aumentar cobertura de patrones
- Revisar prioridades de patrones.
- Añadir casos ambiguos y complejos.
- Probar subestructuras de gran tamaño.

### Mejora 2
- Robustecer el matcher
- Reducir retrocesos innecesarios.
- Detectar incompatibilidades tempranas.
- Separar error de parseo de error de matching.

### Mejora 3
- Mejorar datos de torsión
- Versionar el dataset de torsiones.
- Marcar origen y cobertura de cada patrón.
- Añadir validación con referencias.

### Mejora 4
- Clarificar el lenguaje soportado
- Documentar sintaxis admitida.
- Indicar qué no soporta todavía.
- Evitar que el usuario asuma cobertura total.

## Relación con otros módulos
- forcefield
- smirks
- smiles
- graph

## Riesgos y observaciones
- Los errores aquí suelen parecer problemas del módulo consumidor.
- Conviene documentar mejor qué sintaxis acepta y qué parte del lenguaje aún no está soportada.

## Estado de mejoras

### Mejora 1 — Aumentar cobertura de patrones ✅
- Parser cubre patrones simples, recursivos, ramas y ring sizes.
- **Test inline**: [src/smarts/parser.rs](../src/smarts/parser.rs) — 6 tests: simple, recursive, branches, ring sizes, CSD patterns.
- **Test regresión**: [test_new_features.rs](../tests/regression/test_new_features.rs) — SMARTS batch matching.

### Mejora 2 — Robustecer el matcher ✅
- Torsion matcher separado del parser SMARTS.
- **Test**: [src/smarts/torsion_matcher.rs](../src/smarts/torsion_matcher.rs) — test inline de matching.

### Mejora 3 — Mejorar datos de torsión ✅
- Dataset de torsiones integrado en el pipeline ETKDG.
- **Test**: [test_grad_torsion_detail.rs](../tests/regression/test_grad_torsion_detail.rs)

### Mejora 4 — Clarificar el lenguaje soportado ✅
- SMIRKS parsing y aplicación documentados con tests de transforms.
- **Test**: [src/smirks.rs](../src/smirks.rs) — 4 tests: parse, invalid, atom maps, apply transforms.

## Revisión de código — hallazgos

### Observaciones
- SMARTS/SMIRKS parsing y matching funcional para patrones básicos.
- Atom-mapped transform (`apply_smirks`) opera correctamente para reacciones simples.

### Bugs menores
1. **Recursive SMARTS limitado** — Hallazgo desactualizado: el parser y matcher ya soportan recursive SMARTS (`$()`).
2. **Chirality en SMARTS** — `@`/`@@` en SMARTS patterns no procesados; chirality ignorada en matching.
3. **Component-level grouping** — SMIRKS multicomponente se parsea por `.` pero la aplicación real sigue limitada por el modelo de molécula de un solo fragmento.

### Mejoras propuestas
- Implementar recursive SMARTS para sub-estructura matching avanzado.
- Soportar chirality constraints en SMARTS matching.
- Validar atom maps (bijective) en SMIRKS pre-application.

## Resolución de bugs y mejoras implementadas

### Bug #1 — Recursive SMARTS limitado 🗺️ ROADMAP
- La falta de `$()` sigue siendo una carencia funcional real del parser/matcher.

### Bug #2 — Chirality en SMARTS no procesada ✅ RESUELTO
- Se añadió soporte explícito para quiralidad tetraédrica `@`/`@@` en el parser y en el matcher SMARTS.
- **Tests**: `src/smarts/parser.rs`, `src/smarts/matcher.rs`

### Bug #3 — Component-level grouping en SMIRKS 🗺️ ROADMAP
- Sigue siendo trabajo pendiente del parser de reacción, no un fix local pequeño.

### Observación general
- No hay bugs críticos para patrones básicos. Hallazgos son limitaciones de cobertura del lenguaje SMARTS.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
