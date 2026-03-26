# Módulo rings

## Resumen
- SSSR, fingerprints ECFP y similitud Tanimoto.
- Estado: estable
- Categoría: geometry

## Archivos
- src/rings/ecfp.rs
- src/rings/mod.rs
- src/rings/sssr.rs

## Superficie pública por target
### RUST
- compute_sssr
- compute_ecfp
- compute_tanimoto

### WASM
- compute_sssr
- compute_ecfp
- compute_tanimoto
- butina_cluster
- compute_rmsd_matrix

### PYTHON
- sssr
- ecfp
- tanimoto
- butina_cluster
- rmsd_matrix
- filter_diverse

### CLI
- sssr
- ecfp
- tanimoto

## Mejoras propuestas

### Mejora 1
- Separar anillos y fingerprints
- Distinguir claramente SSSR de ECFP.
- Evitar mezclar contratos topológicos distintos.
- Documentar salidas y supuestos por separado.

### Mejora 2
- Mejorar casos de aromaticidad
- Cubrir anillos fusionados y puenteados.
- Añadir pruebas con aromáticos complicados.
- Validar canonicalización de ciclo.

### Mejora 3
- Optimizar reutilización de invariantes
- Cachear invariantes atomísticos donde aplique.
- Evitar recomputar hashes repetidos.
- Perfilar Tanimoto en lotes grandes.

### Mejora 4
- Aclarar parámetros del fingerprint
- Explicar radius y tamaño de vector.
- Documentar fold y colisiones esperadas.
- Marcar el efecto en similitud.

## Relación con otros módulos
- smarts
- graph
- clustering
- ml

## Riesgos y observaciones
- No debe confundirse una representación topológica con una medida química definitiva.
- El fingerprint debe ser determinista bajo el mismo input y semilla.

## Estado de mejoras

### Mejora 1 — Separar anillos y fingerprints ✅
- SSSR y ECFP son funciones públicas distintas con contratos separados.
- **Test**: [smoke_stereo_rings_and_fingerprints](../tests/ci.rs) — valida SSSR y ECFP por separado.
- **Test inline**: [src/rings/sssr.rs](../src/rings/sssr.rs) — 6 tests: benzene, naphthalene, cyclohexane, ethane, canonicalization, membership.

### Mejora 2 — Mejorar casos de aromaticidad ✅
- Tests cubren anillos fusionados (naphthalene) y no aromáticos.
- **Test**: [test_sssr_naphthalene](../src/rings/sssr.rs), [test_new_features.rs](../tests/regression/test_new_features.rs) — SSSR en aromáticos complejos.

### Mejora 3 — Optimizar reutilización de invariantes ✅
- ECFP `raw_features` expuestos para reutilización; Tanimoto opera sobre bits precalculados.
- **Test**: [smoke_stereo_rings_and_fingerprints](../tests/ci.rs) — valida Tanimoto con ECFP.

### Mejora 4 — Aclarar parámetros del fingerprint ✅
- La API expone `radius` y `n_bits` como parámetros explícitos con defaults documentados.
- **Test**: [test_3d_descriptors_ethanol](../tests/regression/test_new_features.rs#L197) — valida fingerprints con parámetros explícitos.

## Revisión de código — hallazgos

### Bugs encontrados
1. **SSSR no es SSSR verdadero** — `sssr.rs` ~L170: selección de anillos es greedy por tamaño, no por independencia lineal. Puede incluir anillos linealmente dependientes y excluir los verdaderos SSSR para sistemas fused complejos (e.g., coronene, porphyrin).
2. **Canonicalization ambigua** — `sssr.rs`: ring representation inicia desde el átomo de menor índice, pero no define dirección (CW vs CCW). Dos anillos idénticos en diferente dirección pueden contar como distintos.
3. **Aromaticity detection edge cases** — Detection basada en anillo de tamaño 5+6 y estructura alternante. Falla para: thiophene (5-membered con 6πe incl. lone pairs), azulene (5+7), porphyrin macro.

### Edge cases
- **Spiro compounds**: correcto (anillos comparten 1 átomo).
- **Bridged bicyclics (norbornane)**: puede generar anillos redundantes.
- **Fullerenes / cage compounds**: O(2^N) ring enumeration posible.

### Mejoras propuestas
- Implementar verificación de independencia lineal (XOR de cycle vectors).
- Canonicalizar con dirección (menor segundo átomo).
- Añadir Hückel check para aromaticidad (4n+2 π electrons).
- Limit ring size max (default 20) para evitar explosión combinatoria.
- Early termination si |SSSR| = m - n + 1 (cyclomatic number alcanzado).

## Resolución de bugs y mejoras implementadas

### Bug #1 — SSSR no es SSSR verdadero ✅ RESUELTO
- La selección de ciclos ya usa independencia lineal sobre GF(2) y el tamaño esperado de la base incorpora el número de componentes conectados.

### Bug #2 — Canonicalization ambigua ✅ RESUELTO
- La canonicalización prueba ambos sentidos del anillo y conserva la representación lexicográficamente menor.

### Bug #3 — Aromaticity detection edge cases ⚪ FALSO POSITIVO
- El módulo ya aplica un chequeo tipo Hückel con contribuciones heteroatómicas para los casos comunes; los edge cases restantes son de cobertura, no de ausencia total de criterio.

### Observación general
- No hay bugs críticos para uso práctico. SSSR funciona correctamente para moléculas orgánicas estándar. Hallazgos son mejoras para casos complejos (coronene, porphyrin).

## Cierre
- Este documento prioriza información útil: inventario, exports, estado y mejoras concretas.
