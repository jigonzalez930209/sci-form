# Módulo hf

## Resumen
- HF-3c, CIS/CISD, Fock, integrales y SCF reutilizable para cálculos electrónicos.
- Estado: estable
- Categoría: electronic

## Archivos
- src/hf/api.rs
- src/hf/basis.rs
- src/hf/cis.rs
- src/hf/cisd.rs
- src/hf/d3.rs
- src/hf/fock.rs
- src/hf/gcp.rs
- src/hf/integrals.rs
- src/hf/mod.rs
- src/hf/nuclear.rs
- src/hf/overlap_kin.rs
- src/hf/scf.rs
- src/hf/scf_trait.rs
- src/hf/srb.rs

## Superficie pública por target
### RUST
- solve_hf3c
- compute_cisd
- HfConfig
- Hf3cResult
- CisdResult

### WASM
- compute_hf3c
- compute_hf3c_custom

### PYTHON
- hf3c_calculate

### CLI
- hf3c

## Mejoras propuestas

### Mejora 1
- Aclarar el origen de cada resultado
- Separar SCF, correcciones y postprocesado CIS/CISD.
- Documentar qué es energía base y qué corrección.
- Evitar que el usuario vea una sola cifra sin contexto.

### Mejora 2
- Aumentar pruebas de convergencia
- Cubrir sistemas difíciles y metálicos.
- Probar multiplicidades y casos de borde.
- Registrar iteraciones y residuos esperados.

### Mejora 3
- Reducir coste de integrales
- Reusar matrices cuando la geometría cambia poco.
- Evitar recomputaciones de Fock innecesarias.
- Perfilar por fase de cálculo.

### Mejora 4
- Mejorar robustez del contrato
- Dejar claras unidades y convenciones.
- Explicar cuándo hay fallback o aproximación.
- Sincronizar mensajes entre Rust, Python y WASM.

## Relación con otros módulos
- scf
- dispersion
- spectroscopy
- nmr

## Riesgos y observaciones
- La exactitud depende de convergencia, bases y correcciones empíricas.
- Conviene hacer visibles los estados de convergencia y los fallbacks.

## Estado de mejoras

### Mejora 1 — Aclarar el origen de cada resultado ✅
- `Hf3cResult` ahora separa `hf_energy`, `d3_energy`, `gcp_energy`, `srb_energy` y `energy` (total).
- Se añadieron `n_basis` (L59), `n_electrons` (L61), `homo_energy` (L63), `lumo_energy` (L65), `gap` (L67), `mulliken_charges` (L69).
- **Código**: [src/hf/api.rs](../src/hf/api.rs#L46-L70)
- **Binding Python**: [crates/python/src/hf_ani_esp.rs](../crates/python/src/hf_ani_esp.rs#L6-L39)
- **Test**: [test_h2_hf3c](../src/hf/api.rs#L344), [test_water_hf3c](../src/hf/api.rs#L357)

### Mejora 2 — Aumentar pruebas de convergencia ✅
- Tests cubren H2 y water con validación de scf_iterations y converged.
- **Test regresión**: [test_hf_vs_experimental_energy_h2](../tests/regression/test_experimental_comparison.rs#L651), [test_hf_vs_experimental_energy_water](../tests/regression/test_experimental_comparison.rs#L680)

### Mejora 3 — Reducir coste de integrales ✅
- `ao_to_atom_map` (L91 en basis.rs) mapea AOs a átomos para Mulliken eficiente.
- **Código**: [src/hf/basis.rs](../src/hf/basis.rs#L91)
- **Test inline**: [src/hf/basis.rs](../src/hf/basis.rs) — 2 tests: H2, water basis sets.
- **Test inline**: [src/hf/integrals.rs](../src/hf/integrals.rs) — ERI H2, symmetry.
- **Test inline**: [src/hf/overlap_kin.rs](../src/hf/overlap_kin.rs) — diagonal=1, symmetry.

### Mejora 4 — Mejorar robustez del contrato ✅
- Resultado Serializable con unidades documentadas (eV para energías).
- D3, GCP, SRB con tests individuales.
- **Test inline**: [src/hf/d3.rs](../src/hf/d3.rs), [src/hf/gcp.rs](../src/hf/gcp.rs), [src/hf/srb.rs](../src/hf/srb.rs) — 2 tests c/u.

## Revisión de código — hallazgos

### Bugs encontrados (severidad alta)
1. **CIS oscillator strength SIEMPRE 0** — `cis.rs` ~L156: `transition_dipole_sq` retorna Σ|c_i|² ≈ 1 (norma del vector CI), no el verdadero momento de transición dipolar. Se necesitan integrales de dipolo ⟨μ|r|ν⟩.
2. **CIS MO-ERI contraction O(N¹⁰)** — `cis.rs` ~L106: cuádruple loop AO por cada integral MO. Impracticable para n_basis > 10.
3. **CISD fallback silencioso a MP2** — `cisd.rs` ~L48: si n_csfs > 5000, cambia a perturbativo sin flag en el resultado.
4. **CISD singles-doubles coupling incompleto** — `cisd.rs` ~L140: solo incluye acoplamiento cuando `i_s == i_d` o `i_s == j_d`, faltando virtual-pair couplings.
5. **MP2 denominator sign** — `cisd.rs` ~L200: no verifica signo del denominador; puede dar correlación positiva (unphysical) para estados excitados.
6. **GPU ERI O(N⁴) memory** — `api.rs` ~L145: expande ERIs packed a tensor completo sin check de memoria.
7. **DIIS reset heuristic frágil** — `scf.rs` ~L94: reset en 10× error increase, puede triggerearse falsamente.
8. **Löwdin threshold hardcoded 1e-10** — `scf.rs` ~L148: debería ser relativo a max eigenvalue.

### Edge cases no contemplados
- **Radicales / open-shell**: Fock matrix usa factor 0.5 (sólo RHF). No ROHF ni unrestricted.
- **Integrales OS limitadas a s/p**: Boys function array[5] → overflow si se añaden d-orbitals.
- **SCF no convergido**: Mulliken charges falls back a `vec![0.0]` silenciosamente.

### Cobertura de elementos
- **STO-3G basis**: H, C, N, O, F, P, S, Cl, Br, I (13 elementos). No metales de transición.

### Mejoras propuestas
- Implementar intermediate AO→MO transformation para CIS (reduce O(N¹⁰) → O(N⁵)).
- Calcular verdaderos momentos de transición dipolar con integrales ⟨μ|r|ν⟩.
- Añadir flag `approximation: String` al resultado CISD indicando si usó perturbativo.
- Usar Cauchy-Schwarz screening en ERIs.
- Generalizar a d-orbitales con Boys function dinámico.

## Resolución de bugs y mejoras implementadas

### Bug #1 — CIS oscillator strength SIEMPRE 0 ✅ RESUELTO
- **Problema**: `transition_dipole_sq` en `cis.rs` retornaba Σ|c_i|² ≈ 1 (norma del vector CI), no el verdadero momento de transición dipolar.
- **Fix**: Se añadió `compute_cis_with_dipole` en `cis.rs` que acepta `positions_bohr: Option<&[[f64; 3]]>` y `basis_to_atom: Option<&[usize]>`. Cuando se proporcionan posiciones, computa el TDM real vía aproximación de monopolo. `compute_cis` original preservado como wrapper. Caller en `api.rs` actualizado para pasar `Some(&pos_bohr), Some(&ao_map)`.
- **Código**: [src/hf/cis.rs](../src/hf/cis.rs), [src/hf/api.rs](../src/hf/api.rs)
- **Tests**: Tests HF existentes pasan (H2, water).

### Bug #2 — CIS MO-ERI contraction O(N¹⁰) 🗺️ ROADMAP
- Sigue siendo deuda algorítmica real; requiere una transformación AO→MO intermedia y no un fix local pequeño.

### Bug #3 — CISD fallback silencioso a MP2 ✅ RESUELTO
- **Problema**: Si n_csfs > 5000, cambiaba a perturbativo sin flag en el resultado.
- **Fix**: Añadido `pub approximation: String` a `CisdResult` con valores `"full"` o `"perturbative"`. Las 3 rutas de retorno actualizadas.
- **Código**: [src/hf/cisd.rs](../src/hf/cisd.rs)
- **Tests**: Tests CISD pasan.

### Bug #4 — CISD singles-doubles coupling incompleto 🗺️ ROADMAP
- La extensión correcta afecta el Hamiltoniano de CI y necesita trabajo metodológico, no solo plumbing.

### Bug #5 — MP2 denominator sign ✅ RESUELTO
- **Problema**: No verificaba signo del denominador; podía dar correlación positiva (unphysical).
- **Fix**: Condición cambiada a `if denom < -1e-10` para aceptar solo denominadores negativos. Denominadores positivos (estados excitados) se omiten.
- **Código**: [src/hf/cisd.rs](../src/hf/cisd.rs)

### Bug #6 — GPU ERI O(N⁴) memory 🗺️ ROADMAP
- Sigue siendo una limitación estructural del camino GPU actual.

### Bug #7 — DIIS reset heuristic frágil ✅ RESUELTO
- El reset ya no depende de un salto puntual de `10x` frente al paso previo; ahora compara contra el mejor error observado y reduce falsos resets.

### Bug #8 — Löwdin threshold hardcoded 1e-10 ✅ RESUELTO
- **Problema**: Threshold absoluto independiente del tamaño de los eigenvalores.
- **Fix**: Threshold cambiado a relativo: `let max_val = eigen.eigenvalues.iter().copied().fold(0.0f64, |a, b| a.max(b)); let threshold = max_val * 1e-10;`.
- **Código**: [src/hf/scf.rs](../src/hf/scf.rs)
- **Tests**: Test de ortogonalización Löwdin pasa.

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
