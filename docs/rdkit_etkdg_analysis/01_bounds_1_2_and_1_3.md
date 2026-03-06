# Análisis Algorítmico RDKit: Bounds Matrix (Distancias 1-2 y 1-3)

Para alcanzar un error RMSD $< 0.1$ Å respecto a RDKit durante la generación de conformadores 3D mediante ETKDG, es estricto replicar las lógicas aplicadas en **`BoundsMatrixBuilder.cpp`**.

## 1. Distancias 1-2 (Enlaces Directos)

RDKit no usa longitudes de enlace genéricas, sino parámetros del Force Field Universal (UFF) o aproximaciones empíricas.

### Algoritmo:
1. Extraer el orden del enlace ($1.0, 1.5, 2.0, 3.0$).
2. Solicitar `calcBondRestLength(bOrder, atomParams[begId], atomParams[endId])` desde el UFF. Esto cruza propiedades atómicas (radio covalente, electronegatividad) para un valor $r_{\text{eq}}$ exacto.
3. **Aplicar Tolerancias (`DIST12_DELTA = 0.01`):**
   - Límite Superior: $r_{\text{eq}} + 0.01$
   - Límite Inferior: $r_{\text{eq}} - 0.01$
4. **Excepción Empírica (Squish):** Si el enlace está en un anillo conjugado de tamaño 5 y contiene heteroátomos grandes (Número Atómico $> 10$), RDKit amplía la tolerancia agregando $\pm 0.2$ Å para darle más flexibilidad durante el "Triangle Smoothing".
5. **Fallback (Faltan Parámetros):** Si el UFF falla, aproxima la longitud como el promedio de los radios de Van der Waals $bl = (vdw_1 + vdw_2) / 2$, estableciendo tolerancia amplificada $[0.5 \cdot bl, 1.5 \cdot bl]$.

## 2. Distancias 1-3 (Ángulos de Enlace)

Se implementan usando la **Ley de los Cosenos**:
$$ d_{13} = \sqrt{ r_{12}^2 + r_{23}^2 - 2 \cdot r_{12} \cdot r_{23} \cdot \cos(\theta) } $$

### Algoritmo de Hibridación de Ángulos $\theta$:
1. **Átomos en Anillo (`_setRingAngle`):**
   - Si es **`SP2`** y tamaño de anillo $\le 8$, o anillo de tamaño $3, 4$: El ángulo estricto es el de un polígono regular: $\theta = \pi \cdot (1 - \frac{2}{ringSize})$.
   - Si es **`SP3`** en anillo tamaño 5: $\theta = 104^\circ$. En otros tamaños: $109.5^\circ$.
   - **`SP3D`**: $105^\circ$, **`SP3D2`**: $90^\circ$. Resto empírico: $120^\circ$.
2. **Átomos Lineales / Ramificados:**
   - **`SP`**: $180^\circ$ ($\pi$).
   - **`SP2`**: $120^\circ$ ($2\pi / 3$).
   - **`SP3`**: $109.5^\circ$.
   - Coordinación Octaédrica **`SP3D2`**: $135^\circ$ o especialidades estereoquímicas.
3. **Tolerancias (`DIST13_TOL = 0.04`):**
   - Si cualquiera del par 1-2-3 es un átomo grande **`SP2`** (ej. Fósforo/Azufre aromático), se dobla la tolerancia a $\pm 0.08$ Å por cada átomo conflictivo.

Cualquier mínima desviación en $\theta$ o $r_{12}$ introduce una cascada de varianza en el smoothing subsiguiente. RDKit prioriza aplicar la topología del anillo como primera regla global.
