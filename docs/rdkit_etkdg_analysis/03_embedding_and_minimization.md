# Análisis Algorítmico RDKit: Triángulo Smoothing y Embebido ETKDG

Luego de inyectar las distancias físicas locales (1-2, 1-3, 1-4, VDW), la Matriz de Distancias original se encuentra fragmentada y con información faltante entre átomos extremos. La conversión a geometría 3D requiere la resolución del Sistema Métrico.

## 1. Triangle Inequality Smoothing (Floyd-Warshall)

RDKit en **`TriangleSmooth.cpp`** es crítico para resolver distancias $N^2$.
Ningún átomo A, B, C debe violar la de-igualdad $D(A, C) \le D(A, B) + D(B, C)$ o $D(A, C) \ge |D(A, B) - D(B, C)|$.
- RDKit empíricamente itera cada trío, actualizando asimétricamente los límites superior e inferior iteración tras iteración hasta estabilizar la matriz. Sin esta estabilización, el SVD implosionará produciendo complejos auto-valores negativos irreales, elevando el RMSD a $>$ 2.0 Å.

## 2. Embedding y Coord Generation

Se implementa en **`Embedder.cpp`** y **`DistGeomUtils.cpp`**.
### Selección Estocástica:
1. Una vez los Límites superior e inferior son válidos globalmente, se selecciona una matriz de de distancias "Randomizada" (DistMat):
   - $d_{ij} = random(Lower_{ij}, Upper_{ij})$.
2. Con la `DistMat`, se calcula la Metrica Geocéntrica a través del Teorema de Cayley-Menger.
3. Se sacan los 3 Vector Propios Principales (3 eigenvalues mayores de la Matriz). Dichos vectores actúan como coordenadas 3D para la configuración atómica en crudo.

## 3. Minimización por Campo de Fuerza Local

Una vez se poseen las Coordenadas SVD:
- ETKDG ("Experimental Torsion Knowledge Distance Geometry") somete rápidamente esas coordenadas a un ForceField.
- ETKDGv3 usa un `basinThresh=5.0` para L-BFGS. Las coordenadas se optimizan usando Tolerancia Fuerza = `1e-3`. RDKit evalúa el sistema en un modelo $O(N)$ optimizado o con paralelismo OpenMP, permitiendo correr miles de moléculas sin explotar el cálculo.
- **Chequeo Topológico Posterior:** Las coordenadas pasan por escáner `finalChiralChecks` y de volumen Tetraédrico. Si un Centro Quiral perdió estabilidad o la energía sobrepasa $0.05$ por átomo, el modelo falla o reminimiza forzando en un espacio vectorial de 4-Dimensiones.

El uso exhaustivo de `ETKDGDetails` para incluir torsiones impropias asegura que configuraciones inestables no terminen subiendo el RMSD final de las pruebas. Para conseguir `< 0.1` Å en sci-form, es necesario pasar de "Steepest descent" a "L-BFGS" e imitar las torsiones MMFF94 del punto 3.
