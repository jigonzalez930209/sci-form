# Roadmap: VeloSci 3D Conformer Engine (RDKit -> Rust)

**Objetivo Principal:** Traducir y aislar el pipeline de generación de conformaciones 3D de RDKit (C++) a Rust idiomático, eliminando dependencias pesadas (Boost, Eigen) para compilar un módulo `.wasm` ultra-ligero y rápido para la web.

**Pila Tecnológica:**
* **Lenguaje:** Rust (Target: `wasm32-unknown-unknown`)
* **Álgebra Lineal:** `nalgebra` o `ndarray` (Reemplazo de Eigen)
* **Teoría de Grafos:** `petgraph` (Reemplazo de Boost.Graph)
* **Testing Oracle:** Python + RDKit nativo

---

## Fase 1: Infraestructura y Estructuras de Datos Base (Mes 1)
El objetivo de esta fase no es hacer química, sino preparar el terreno matemático y la arquitectura de pruebas.

- [ ] **1.1 Setup del Entorno y Repositorio:**
  - Inicializar el proyecto con `cargo new --lib velosci-embed`.
  - Configurar `wasm-bindgen` para las pruebas iniciales de interoperabilidad.
- [ ] **1.2 Pipeline de Differential Testing:**
  - Crear un script en Python (`oracle.py`) que tome un dataset de 1000 SMILES de prueba, genere las coordenadas 3D con RDKit (`AllChem.EmbedMolecule`) y las guarde en formato JSON o XYZ.
  - Escribir la suite de tests en Rust que calcule el RMSD entre las coordenadas de Rust y las de C++.
- [ ] **1.3 Replicar `ROMol` (El Grafo Molecular):**
  - Implementar la estructura central usando `petgraph`.
  - Definir `Atom` (con propiedades: número atómico, hibridación, carga) y `Bond` (orden de enlace, aromaticidad).
  - *Opcional/Atajo:* Evaluar usar crates de Rust existentes como `purr` para el parsing de SMILES a Grafo 2D y no reinventar esa rueda.

## Fase 2: Geometría de Distancias y Matriz de Límites (Meses 2-3)
Aquí comienza la traducción pura de RDKit. Nos enfocaremos en el *namespace* `DistGeom`.

- [ ] **2.1 Implementar la *Bounds Matrix*:**
  - Escribir la lógica para calcular las distancias mínimas y máximas permitidas entre átomos basándose en las topologías de enlace (1-2, 1-3, 1-4).
  - Programar la suavización de la matriz usando el algoritmo del camino más corto de Floyd-Warshall o Dijkstra (traducido de RDKit/Boost).
- [ ] **2.2 Embebimiento Métrico (Metric Geometry):**
  - Traducir la selección de distancias aleatorias dentro de los límites de la matriz.
  - Implementar la matriz de distancias cuadráticas y el centrado de la matriz.
- [ ] **2.3 Diagonalización y Proyección 3D:**
  - Usar `nalgebra` para calcular los autovalores y autovectores (Eigenvalue decomposition) de la matriz de distancias.
  - Proyectar los 3 autovectores principales para obtener las coordenadas $(x, y, z)$ iniciales en bruto.

## Fase 3: Refinamiento Químico y Estereoquímica (Mes 4)
Las coordenadas de la Fase 2 serán físicamente inestables. Necesitamos aplicar reglas empíricas.

- [ ] **3.1 Percepción de Estereoquímica:**
  - Implementar la lógica para detectar centros quirales (R/S) e isómeros (E/Z) a partir del grafo 2D.
  - Ajustar la matriz de límites para forzar los volúmenes quirales correctos.
- [ ] **3.2 Implementar ETKDG (Opcional pero altamente recomendado):**
  - Traducir el algoritmo *Experimental Torsion-Angle Preference* de RDKit.
  - Agregar las penalizaciones de energía basadas en bibliotecas de torsión experimentales (esto es lo que hace que los anillos se plieguen correctamente).

## Fase 4: Minimización de Energía (Force Fields) (Meses 5-6)
Las coordenadas crudas necesitan relajarse para alcanzar el mínimo de energía local.

- [ ] **4.1 Estructura del Campo de Fuerza (UFF):**
  - Traducir los parámetros del *Universal Force Field* (UFF) desde RDKit.
  - Implementar las funciones de energía para enlaces (stretching), ángulos (bending) y torsiones (dihedrals).
  - Implementar las fuerzas no enlazadas (Van der Waals / Lennard-Jones).
- [ ] **4.2 Motor de Minimización:**
  - Escribir un algoritmo iterativo (Steepest Descent o L-BFGS, idealmente L-BFGS por eficiencia) usando `nalgebra` para calcular los gradientes de energía y mover las posiciones atómicas iterativamente.

## Fase 5: Optimización WASM y API Pública (Mes 6+)
Hacer que la librería sea usable en la web para el ecosistema VeloSci.

- [ ] **5.1 Bindings de WASM:**
  - Exponer una API limpia en `wasm-bindgen` (ej. `generate_3d_conformer(smiles: &str) -> String`).
- [ ] **5.2 Benchmarking:**
  - Medir el tiempo de ejecución en el navegador contra RDKit-JS (MinimalLib).
  - Optimizar cuellos de botella en la asignación de memoria (Memory Profiling).
- [ ] **5.3 Reducción de Tamaño:**
  - Configurar `Cargo.toml` con `lto = true`, `opt-level = 'z'`, y stripping para asegurar que el binario `.wasm` pese menos de 1MB.

---
> **Nota de Desarrollo:** No intentes traducir archivo por archivo linealmente. C++ orientado a objetos con herencia múltiple no se traduce bien a Rust. Extrae la **lógica matemática** de las funciones en C++ e impleméntala usando la composición basada en traits de Rust.