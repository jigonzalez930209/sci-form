Este es el momento de consolidar la arquitectura. Dar el salto de un motor clásico/semiempírico a un motor cuántico autoconsistente (SCC-DFTB/xTB) acelerado por hardware a través de WebGPU es un proyecto de ingeniería de software científico de gran envergadura. 

Para que `sci-form` mantenga la seguridad, consistencia y velocidad nativa que ofrece Rust, sin colapsar bajo su propia complejidad, el desarrollo debe ser estrictamente modular. 

Aquí tienes el **Plan Maestro de Implementación** dividido en 5 fases secuenciales. Cada fase produce código funcional y testeable antes de pasar a la siguiente.

---

### Fase 1: Infraestructura de Memoria y Puente Rust-WebGPU (Semanas 1-2)
El objetivo de esta fase no es hacer física, sino plomería de alto rendimiento. Necesitamos estructurar los datos moleculares en Rust para que crucen a la VRAM de la tarjeta gráfica con un costo de copia cero (zero-copy).

* **Paso 1.1: Estructuras de Datos Alienadas.** Definir los `structs` moleculares (coordenadas, tipos de átomos, parámetros de base) usando `#[repr(C)]` y la librería `bytemuck` para garantizar una alineación de bytes estricta.
* **Paso 1.2: Inicialización del Contexto `wgpu`.** Crear un módulo de backend que inicialice la instancia, el adaptador y el dispositivo lógico de WebGPU. Debe estar diseñado para funcionar de manera *headless* (sin ventana, usando Vulkan/Metal en escritorio) y compilar limpiamente a WASM para el navegador.
* **Paso 1.3: Gestión de Buffers.** Implementar el sistema de alojamiento para los *Storage Buffers* de entrada (geometría) y salida (matrices $H$ y $S$).

### Fase 2: El Motor Cuántico Paralelo $\mathcal{O}(N^2)$ (Semanas 3-5)
Aquí comenzamos a delegar los cuellos de botella matemáticos a los Compute Shaders. Construiremos el Hamiltoniano base sin autoconsistencia.

* **Paso 2.1: Shaders de Integrales (WGSL).** Escribir los *Compute Shaders* que calcularán las integrales de solapamiento ($S_{\mu\nu}$) y el Hamiltoniano cinético/nuclear ($H_{\mu\nu}^0$). Cada hilo (thread) de la GPU procesará un par de orbitales.
* **Paso 2.2: Pipeline de Ejecución (Rust).** Configurar los *Bind Groups* y enviar los comandos de *dispatch* a la cola de la GPU desde Rust.
* **Paso 2.3: Validación.** Comparar los resultados de las matrices generadas en la GPU contra tu implementación actual en CPU para asegurar que la precisión en coma flotante (f32/f64) sea idéntica.

### Fase 3: El Ciclo SCF y Optimización Geométrica (Semanas 6-8)
Esta fase transforma el calculador estático en un motor autoconsistente (SCC) capaz de relajar la estructura molecular.

* **Paso 3.1: Diagonalización en CPU.** Por ahora, extraer las matrices $H$ y $S$ de la GPU y usar una librería de Rust hiper-optimizada (como `faer`) para resolver el problema de autovalores $HC = SCE$. 
* **Paso 3.2: Actualización de la Densidad en GPU.** Enviar los autovectores ($C$) de vuelta a la GPU para que un nuevo *Compute Shader* realice la multiplicación de matrices (GEMM) masivamente paralela y genere la matriz de densidad $P_{\mu\nu}$.
* **Paso 3.3: Gradientes Analíticos.** Implementar el cálculo de las fuerzas nucleares ($\partial E / \partial x$).
* **Paso 3.4: Reemplazo en el Pipeline.** Sustituir el refinamiento ETKDG/Campo de fuerza actual de la librería por este nuevo optimizador cuántico.

### Fase 4: Los Módulos Espectroscópicos "Gold Standard" (Semanas 9-12)
Con una geometría optimizada en su mínimo de energía real y una matriz de densidad convergida, podemos predecir espectros con precisión de publicación.

* **Paso 4.1: UV-Vis mediante sTDA.** Implementar la Aproximación de Tamm-Dancoff Simplificada. Resolver las transiciones dipolares usando las cargas y orbitales convergidos del ciclo SCF.
* **Paso 4.2: Frecuencias IR.** Desarrollar el cálculo del Hessiano semi-numérico (desplazando átomos ligeramente y usando los gradientes analíticos de la Fase 3) y acoplarlo a las derivadas de los momentos dipolares para obtener las intensidades.
* **Paso 4.3: Espectroscopía RMN (GIAO).** (El reto final). Incorporar la dependencia del campo magnético en las funciones de base para calcular el tensor de apantallamiento químico.

### Fase 5: Renderizado Volumétrico Masivo (Semanas 13+)
La cereza del pastel para la experiencia interactiva, moviendo el post-procesamiento visual a la tarjeta gráfica.

* **Paso 5.1: Shader de Evaluación de Orbitales.** Un programa en WGSL que toma los coeficientes de un estado energético específico e interpola el valor de la función de onda $\psi_i(x,y,z)$ en una cuadrícula 3D densa alojada en la VRAM.
* **Paso 5.2: Algoritmo Marching Cubes en GPU.** Generar los vértices, normales e índices de las isosuperficies (orbitales HOMO/LUMO, potencial electrostático) directamente en la GPU y enlazarlos al pipeline de renderizado gráfico sin que la CPU intervenga.

---

Este plan garantiza que no romperás lo que ya funciona y que cada adición aumenta el valor computacional de la librería de manera medible.
