# Roadmap: VeloSci 3D Conformer Engine (RDKit -> Rust)

**Objetivo Principal:** Traducir y aislar el pipeline de generación de conformaciones 3D de RDKit (C++) a Rust idiomático, eliminando dependencias pesadas (Boost, Eigen) para compilar un módulo `.wasm` ultra-ligero y rápido para la web. Una vez estable, el motor evolucionará a una plataforma de simulación de materiales porosos y electroquímica nativa en la web (Visión Post-Core).

**Pila Tecnológica:**
* **Lenguaje:** Rust (Target: `wasm32-unknown-unknown`)
* **Álgebra Lineal:** `nalgebra` o `ndarray` (Reemplazo de Eigen)
* **Teoría de Grafos:** `petgraph` (Reemplazo de Boost.Graph)
* **Testing Oracle:** Python + RDKit nativo

---

## 🏗 Estructura de Datos Base Propuesta (sci-form)
Para soportar componentes avanzados (como periodicidad para MOFs, análisis electrostático e integración de memoria "zero-copy" con Arrow), la estructura en Rust debe diseñarse con visión a futuro desde el inicio:

```rust
use nalgebra::Vector3;

pub struct Atom {
    pub element: u8,            // Número atómico
    pub position: Vector3<f32>, // Coordenadas (x, y, z)
    pub charge: f32,            // Carga parcial calculada (Fase 6)
    pub formal_charge: i8,      
    pub hybridization: Hybridization,
}

pub struct Bond {
    pub start: usize,           // Índice del átomo A
    pub end: usize,             // Índice del átomo B
    pub order: BondOrder,       // Single, Double, Triple, Aromatic
}

pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub metadata: MoleculeMeta, // Nombre, SMILES original, etc.
    // Para la Fase 7 (Materiales Complejos):
    pub unit_cell: Option<UnitCell>, 
}
```

> **Nota de Desarrollo:** No intentes traducir archivo por archivo linealmente. C++ orientado a objetos con herencia múltiple no se traduce bien a Rust. Extrae la **lógica matemática** de las funciones en C++ e impleméntala usando la composición basada en traits de Rust.

---

## ⚙️ Módulos Core: Generador de Conformaciones 3D

### Fase 1: Infraestructura y Estructuras de Datos Base (Mes 1)
El objetivo de esta fase no es hacer química, sino preparar el terreno matemático y la arquitectura de pruebas.

- [ ] **1.1 Setup del Entorno y Repositorio:**
  - Inicializar el proyecto con `cargo new --lib sci-form`.
  - Configurar `wasm-bindgen` para las pruebas iniciales de interoperabilidad.
- [ ] **1.2 Pipeline de Differential Testing:**
  - Crear un script en Python (`oracle.py`) que tome un dataset de 1000 SMILES de prueba, genere coordenadas 3D locales con RDKit (`AllChem.EmbedMolecule`) y las guarde en JSON o XYZ.
  - Escribir la suite de tests en Rust que calcule el RMSD entre las coordenadas generadas por el WASM en Rust y las de C++ como fuente de verdad.
- [ ] **1.3 Replicar `ROMol` (El Grafo Molecular):**
  - Implementar la estructura central usando `petgraph` basada en la propuesta de `Atom`, `Bond` y `Molecule`.
  - *Opcional/Atajo:* Evaluar usar crates de Rust existentes como `purr` para el parsing de SMILES a Grafo 2D y no reinventar esa rueda.

### Fase 2: Geometría de Distancias y Matriz de Límites (Meses 2-3)
Aquí comienza la traducción pura de RDKit (*namespace* `DistGeom`).

- [ ] **2.1 Implementar la *Bounds Matrix*:**
  - Escribir la lógica para calcular las distancias mínimas y máximas permitidas entre átomos basándose en las topologías de enlace (1-2, 1-3, 1-4).
  - Programar la suavización de la matriz usando el algoritmo del camino más corto de Floyd-Warshall o Dijkstra.
- [ ] **2.2 Embebimiento Métrico (Metric Geometry):**
  - Traducir la selección de distancias aleatorias dentro de los límites de la matriz.
  - Implementar la matriz de distancias cuadráticas y el centrado de la matriz.
- [ ] **2.3 Diagonalización y Proyección 3D:**
  - Usar `nalgebra` para calcular los autovalores y autovectores (Eigenvalue decomposition) de la matriz de distancias.
  - Proyectar los 3 autovectores principales para obtener las coordenadas $(x, y, z)$ iniciales en bruto.

### Fase 3: Refinamiento Químico y Estereoquímica (Mes 4)
Las coordenadas de la Fase 2 serán físicamente inestables. Se aplican reglas empíricas para acercarlas a la realidad química.

- [ ] **3.1 Percepción de Estereoquímica:**
  - Implementar la lógica para detectar centros quirales (R/S) e isómeros (E/Z) a partir del grafo 2D.
  - Ajustar la matriz de límites para forzar los volúmenes quirales correctos.
- [ ] **3.2 Implementar ETKDG:**
  - Traducir el algoritmo *Experimental Torsion-Angle Preference* de RDKit.
  - Agregar penalizaciones de energía basadas en bibliotecas de torsión experimentales para garantizar pliegues de anillos correctos.

### Fase 4: Minimización de Energía (Force Fields) (Meses 5-6)
Las coordenadas crudas necesitan relajarse para alcanzar el mínimo de energía local.

- [ ] **4.1 Estructura del Campo de Fuerza (UFF):**
  - Traducir los parámetros del *Universal Force Field* (UFF) desde RDKit.
  - Implementar las funciones de energía para enlaces (stretching), ángulos (bending) y torsiones (dihedrals).
  - Implementar las fuerzas no enlazadas (Van der Waals / Lennard-Jones).
- [ ] **4.2 Motor de Minimización:**
  - Escribir un algoritmo iterativo (preferentemente L-BFGS por eficiencia) usando `nalgebra` para calcular gradientes de energía y relajar la molécula iterativamente.

### Fase 5: Optimización WASM y API Pública (Mes 6+)
Hacer que la librería sea usable en la web para el ecosistema VeloSci.

- [ ] **5.1 Bindings de WASM:**
  - Exponer una API limpia en `wasm-bindgen` (ej. `generate_3d_conformer(smiles: &str) -> String`).
- [ ] **5.2 Benchmarking:**
  - Medir tiempo de ejecución del código WASM en el navegador contra RDKit-JS (MinimalLib).
  - Optimizar cuellos de botella en asignación de memoria.
- [ ] **5.3 Reducción de Tamaño:**
  - Ajustar `Cargo.toml` (`lto = true`, `opt-level = 'z'`, stripping) buscando un WASM ultraligero (< 1MB).

---

## 🚀 Módulos "Killer" Post-Core: Visión de Futuro (sci-form)
Estas características transformarán este motor en una plataforma innovadora de simulación inorgánica web que superará a las herramientas web convencionales enfocadas solo en industria farmacéutica.

### Fase 6: Análisis Electrostático y Superficies (Meses 7-8)
**Objetivo:** Dotar a la molécula de "personalidad" química y modelar propiedades de interacción. El motor debe conocer la distribución electrónica para simulaciones reales de reactividad.
- [ ] **6.1 Implementación Nativa de Cargas Parciales (Gasteiger-Marsili):**
  - Desarrollar un solver iterativo para la ecualización de electronegatividad.
  - Asignar cargas parciales $q_i$ útiles en posteriores simulaciones (ej. Espectroscopía de Impedancia Electroquímica o reactividad de superficies).
- [ ] **6.2 Módulo SASA (Solvent Accessible Surface Area):**
  - Integrar un algoritmo tipo Shrake-Rupley acelerado en Rust en el paso de generación 3D.
  - Calcular Área de Superficie Accesible al Disolvente, volumen libre y tamaño de cavidades (vital para capacidades de filtrado y almacenamiento en materiales porosos).

### Fase 7: Ensamblaje de Materiales Complejos (Meses 9-10)
**Objetivo:** Romper la barrera de las moléculas orgánicas aisladas y crear un módulo de Ensamblaje Topológico para modelar inorgánica y polímeros periódicos (MOFs) en el navegador.
- [ ] **7.1 Arquitectura SBU (Secondary Building Units) y Nodos:**
  - Definir la lógica y algoritmia para centros metálicos y coordinación de redes de manera estérica (ej. acople geométrico de ligandos orgánicos con iones compuestos).
  - Permitir acople dinámico sin requerir pesados programas de escritorio.
- [ ] **7.2 Soporte de Periodicidad (Celdas Unitarias):**
  - Añadir estructuras y algoritmos de vectores de celda con condiciones de contorno periódicas (PBC).
  - Viabilizar cálculos y visualización en WASM de redes cristalinas infinitas por instanciación.

### Fase 8: Optimización de Datos para la Web (Mantenimiento Continuo)
**Objetivo:** Optimización de latencia en browser en interops JS <-> WASM. Eliminar cuellos de botella serializados para una experiencia interactiva instantánea.
- [ ] **8.1 Transferencia de Memoria "Zero-Copy" con Apache Arrow:**
  - Implementar soporte nativo del layout de memoria Arrow para compartir RAM directo con el Frontend de VeloSci. 
  - Evadir completamente la ineficiente traducción a JSON (`JSON.parse()`) al intercambiar matrices de 50,000+ coordenadas atómicas.
- [ ] **8.2 Multithreading con Web Workers:**
  - Paralelizar cálculos pesados (generación múltiple de conformadores, iteraciones del Force Field) usando librerías de concurrencia de Rust (ej. `rayon`) adaptadas a WASM (vía SharedArrayBuffer).
