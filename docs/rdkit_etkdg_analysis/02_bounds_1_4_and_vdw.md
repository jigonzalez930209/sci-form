# Análisis Algorítmico RDKit: Bounds Matrix (Distancias 1-4 y VDW)

Para el rango topológico 1-4 y mayor ($N \ge 4$), la contribución conformacional depende drásticamente de ángulos de torsión (Diedros) e interrupciones estéricas (VDW).

## 1. Distancias 1-4 (Torsiones y Cadenas)

Para átomos separados de topología A-B-C-D, RDKit evalúa cada caso dependiendo del enlace principal (B-C).
El error principal de los motores genéricos radica en obviar que un diedro dicta la distancia 1-4.

### Reglas:
1. Distancias Cis ($D_{\text{cis}}$) calculadas usando Ley de Cosenos espaciales forzando Diedro = $0^\circ$.
2. Distancias Trans ($D_{\text{trans}}$) calculadas usando Ley de Cosenos espaciales forzando Diedro = $180^\circ$ ($\pi$).
3. **Casos Predefinidos (Estereoquímica Exacta):**
   - Enlaces Dobles (ej. `C=C` o amidas estables): RDKit evalúa la configuración `STEREOZ/STEREOCIS` obligando a que los límites sean estrechos alrededor del valor cis: Límite $[D_{\text{cis}} - 0.06, D_{\text{cis}} + 0.06]$.
   - `STEREOE/STEREOTRANS` fija el Límite en $[D_{\text{trans}} - 0.06, D_{\text{trans}} + 0.06]$.
4. **Rango Libre Genérico:** Para rotos simples sin preferencias topológicas, RDKit establece Límite Inferior en $D_{\text{cis}}$ y Superior en $D_{\text{trans}}$ (o viceversa), añadiendo un `GEN_DIST_TOL = 0.06` en los bordes.
5. **Casos Disulfuro (`S-S`):** Se asume forzosamente un diedro exacto de $90^\circ$ ($\pi/2$) para la distancia `C-S-S-C`.
6. **Amidas y Ésteres:** Evalúa profundamente patrones `O=C-N` o `O=C-O` (macro-anillos y cadenas). Dependiendo de las flags (ej. `forceTransAmides`), marca cis o trans restrictivamente y bloquea variaciones locas de torsión dentro de un margen estricto $\pm 0.06$ Å.

## 2. Interacciones de Van Der Waals ($\ge$ 1-5 topológicamente)

Los átomos que caen a una distancia de enlaces $\ge 4$ o de otras ramas utilizan los radios de VDW para garantizar que el sistema molecular **no colapse**. RDKit tiene reglas blandas para contactos 1-5 respecto al resto.

### Algoritmo (`setLowerBoundVDW`):
Sea $vw_1$ el radio de VDW del átomo 1 y $vw_2$ el radio de VDW del átomo 2.
- **Topología de 5 enlaces (1-5 Distancias) ($d=4.0$ en topological matrix):** La fuerza de repulsión se debilita. Límite interior escalado a `VDW_SCALE_15 = 0.7`. Límite Inferior = $0.7 \cdot (vw_1 + vw_2)$.
- **Topología de 6 enlaces (1-6 Distancias) ($d=5.0$):** Escalado intermedio. Límite Inferior = $0.85 \cdot (vw_1 + vw_2)$.
- **Otras Interacciones (Lejanas):** Se impone la impenetrabilidad sólida: Límite Inferior = $vw_1 + vw_2$.

Ciertas variaciones de ETKDG (versiones 2 y 3) alteran estas torsiones introduciendo perfiles probabilísticos derivados del PDB (CSD - Cambridge Structural Database).
