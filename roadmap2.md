Apéndice al ROADMAP.md: Visión de Futuro (Post-Core)Una vez que el motor de embebido 3D sea estable, sci-form evolucionará de un clon de utilidades a una plataforma de simulación de materiales porosos y electroquímica nativa en la web.
Fase 6: Análisis Electrostático y Superficies (Meses 7-8)Objetivo: Dotar a la molécula de "personalidad" química más allá de su forma.
  [ ] 6.1 Implementación de Gasteiger-Marsili:Crear un solver iterativo para la ecualización de electronegatividad.Asignar cargas parciales $q_i$ a cada átomo en el struct.
  [ ] 6.2 Módulo SASA (Solvent Accessible Surface Area):Implementar el algoritmo de Shrake-Rupley acelerado por Rust.Calcular volumen de poros y diámetros de ventana (crucial para MOFs).
Fase 7: Ensamblaje de Materiales Complejos (Meses 9-10)Objetivo: Romper la barrera de las "moléculas pequeñas" para entrar en el diseño de redes.
  [ ] 7.1 Arquitectura SBU (Secondary Building Units):Crear lógica para definir centros metálicos y puntos de coordinación.Implementar el "acoplamiento geométrico" de ligandos orgánicos a nodos metálicos.
  [ ] 7.2 Soporte de Periodicidad (Celdas Unitarias):Añadir soporte para vectores de celda y condiciones de contorno periódicas (PBC).Permitir la visualización de cristales infinitos mediante instanciación en WASM.
Fase 8: Optimización de Datos para la Web (Mantenimiento Continuo)Objetivo: Hacer que la comunicación Rust <-> JS sea invisible y ultrarrápida.
  [ ] 8.1 Integración con Apache Arrow:Implementar el layout de memoria de Arrow para las coordenadas y cargas.Eliminar el overhead de JSON.parse() en el frontend de VeloSci.
  [ ] 8.2 Multithreading con Web Workers:Paralelizar la generación de múltiples conformadores usando rayon compilado a WASM (requiere SharedArrayBuffer).Propuesta de Estructura de Datos Base (sci-form)Para soportar todo lo anterior, el struct principal en Rust debe ser algo así. ¿Qué te parece esta arquitectura para empezar?Rustuse nalgebra::Vector3;

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
    // Para la Fase 7:
    pub unit_cell: Option<UnitCell>, 
}


Para que esta librería pase de ser "un generador 3D" a una herramienta verdaderamente innovadora en el ecosistema científico, aquí tienes 4 módulos "killer" que deberíamos considerar en el diseño base:1. Cálculo Nativo de Cargas Parciales (Esencial para Electroquímica)Tener las coordenadas $(x, y, z)$ es solo la mitad de la física. Para hacer simulaciones moleculares reales, el motor necesita saber cómo se distribuyen los electrones.La idea: Implementar un estimador ultrarrápido como el algoritmo de ecualización de electronegatividad de Gasteiger-Marsili.Por qué es clave: Si en algún momento el software necesita escalar para evaluar propiedades de conductividad, reactividad de superficies o integrarse con módulos de simulación de Espectroscopía de Impedancia Electroquímica, las posiciones geométricas no sirven de nada sin la carga electrostática de cada átomo.2. Ensamblaje de Retículos y Nodos (Soporte para MOFs)Las librerías tradicionales como RDKit o OpenBabel fueron diseñadas casi exclusivamente para la industria farmacéutica (moléculas orgánicas pequeñas y aisladas). Odian la química inorgánica y las estructuras periódicas infinitas.La idea: Crear un módulo de "Ensamblaje Topológico". En lugar de solo leer un SMILES, la librería podría tomar ligandos orgánicos pre-calculados e intentar acoplarlos geométricamente alrededor de un centro metálico (imagina coordinar ligandos con iones complejos derivados de nitrato de cobre o similares).Por qué es clave: Te daría la capacidad de generar y visualizar unidades de construcción secundarias (SBUs) para modelar Metal-Organic Frameworks (MOFs) directamente en el navegador, algo que actualmente requiere software de escritorio pesadísimo.3. Transferencia de Memoria "Zero-Copy" (Apache Arrow)Este es un problema puramente de ingeniería de software, pero vital para WASM.La idea: Cuando tu motor en Rust genere una matriz con 50,000 coordenadas para una macromolécula, serializar eso a un string JSON para enviarlo a la interfaz en JavaScript hundirá el rendimiento del navegador. Podríamos integrar soporte nativo para Apache Arrow en memoria.Por qué es clave: Permite que Rust y JavaScript compartan exactamente el mismo bloque de memoria RAM para leer las matrices de coordenadas de los átomos sin tener que copiarlas ni traducirlas, logrando un renderizado 3D en la pantalla virtualmente instantáneo.4. Cálculo de Superficie y Volumen de Poros (SASA)Ya que tenemos la molécula en 3D, el siguiente paso analítico más demandado es saber cuánto espacio ocupa y qué tan accesible es.La idea: Integrar un algoritmo para calcular el Área de Superficie Accesible al Disolvente (SASA) y el volumen libre.Por qué es clave: Volviendo a los materiales porosos, la capacidad de un material para almacenar gases o filtrar moléculas depende enteramente del tamaño de sus cavidades. Tener esta métrica calculada directamente en el paso de generación 3D ahorraría horas de post-procesamiento.Si logramos empaquetar la generación 3D, la asignación de cargas y una interfaz de memoria ultrarrápida en un solo binario .wasm de un par de megabytes, estaríamos creando algo que no existe actualmente en el panorama del código abierto.