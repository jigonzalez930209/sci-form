# Portabilidad del backend GPU en sci-form

Este documento resume qué tan portable es el backend GPU actual de `sci-form`, en qué sistemas puede funcionar y cuáles son sus límites reales.

## Resumen corto

Sí, el backend está planteado para ser multiplataforma y multi-vendor.

No depende de CUDA, ROCm ni de una marca concreta de GPU. La implementación usa `wgpu`, que abstrae el backend nativo del sistema operativo y selecciona un adaptador compatible automáticamente.

Eso significa que, en principio, puede correr en:

- macOS con GPUs compatibles con Metal
- Windows con GPUs/drivers compatibles con Direct3D 12 o Vulkan
- Linux con GPUs/drivers compatibles con Vulkan

Y puede hacerlo sobre hardware de:

- NVIDIA
- AMD
- Intel
- Apple Silicon

Pero no significa que funcione en cualquier tarjeta sin excepción. La compatibilidad real depende de cuatro cosas:

- que se compile con la feature `experimental-gpu`
- que el sistema operativo exponga un backend soportado por `wgpu`
- que el driver de la GPU esté bien instalado y sea suficientemente moderno
- que el adaptador soporte los límites de cómputo y buffers que necesitan estos kernels

## Cómo funciona hoy la selección del backend

El contexto GPU de `sci-form` usa `GpuContext::best_available()` y trata de crear un dispositivo `wgpu` nativo.

El flujo actual es:

1. Si `experimental-gpu` no está habilitada, no se intenta usar GPU.
2. Si la feature está habilitada, `wgpu` intenta pedir un adaptador de alto rendimiento.
3. Si encuentra un adaptador y puede crear el dispositivo, se usa GPU.
4. Si no encuentra adaptador o falla la inicialización, el código cae a CPU.

En otras palabras: el comportamiento no es “GPU o error”, sino “GPU si está disponible; CPU si no”.

## Por sistema operativo

## macOS

En macOS, lo normal es que `wgpu` use Metal.

Eso incluye normalmente:

- Macs con Apple Silicon
- Macs Intel con GPUs compatibles con Metal

No hace falta que la GPU sea NVIDIA, AMD o Intel en particular. Lo importante es que el sistema tenga soporte Metal operativo.

## Windows

En Windows, lo habitual es usar Direct3D 12, y en algunos equipos también puede entrar Vulkan según el entorno y los drivers.

En la práctica, esto puede funcionar en:

- NVIDIA
- AMD
- Intel

siempre que el driver soporte bien el backend que `wgpu` termine usando.

## Linux

En Linux, el camino importante para este proyecto es Vulkan.

Otra vez, no depende de la marca de la GPU, pero sí del estado del stack gráfico:

- driver Mesa o propietario
- soporte Vulkan disponible
- permisos/dispositivo visibles desde el entorno donde corre el proceso

## Lo que sí significa “multiplataforma” aquí

Significa que el código shader y el pipeline de cómputo están escritos contra una abstracción portable, no contra una API propietaria de un solo proveedor.

Eso permite reutilizar los mismos kernels en equipos distintos sin mantener una versión CUDA para NVIDIA, otra Metal para Apple y otra DirectX para Windows.

## Lo que no significa

No significa que cualquier GPU vieja, cualquier VM gráfica o cualquier driver roto vaya a funcionar.

Tampoco significa que todos los adaptadores den el mismo rendimiento.

Por ejemplo:

- una iGPU de Intel puede ejecutar el kernel pero ser más lenta que una GPU discreta
- una Apple GPU puede ejecutar bien el kernel, pero con un perfil de rendimiento distinto al de una NVIDIA
- un driver viejo puede aceptar el adaptador pero fallar al crear buffers o compilar shaders complejos

## Estado actual en sci-form

El backend GPU actual es nativo y está pensado para escritorio. En el código actual, la inicialización GPU está condicionada por:

- feature `experimental-gpu`
- objetivo nativo, no `wasm32`

Eso quiere decir que este camino no está habilitado hoy para WebAssembly en navegador con la misma infraestructura.

## Estado actual en WASM / navegador

Ahora existe una ruta específica para WASM con WebGPU, separada del backend nativo síncrono.

Esto es importante: en navegador no conviene forzar el mismo modelo del runtime nativo porque WebGPU en WASM es asíncrono por naturaleza.

Por eso el soporte WebGPU en WASM se plantea como una capa específica, manteniendo el backend nativo existente casi intacto.

### Qué soporta esta primera base WASM

- inicialización asíncrona de WebGPU
- consulta de estado del runtime WebGPU
- ejecución acelerada en WASM para rutas experimentales Tier 3
- selector de modo `cpu`, `gpu`, `hybrid` y `auto`

En esta iteración, la ruta WASM acelerada se conectó sobre todo a cálculos donde el patrón O(N²) ya estaba muy claro:

- EEQ
- D4
- CPM

### Qué significa “híbrido” en WASM

En navegador, “híbrido” no significa exactamente lo mismo que en un binario nativo de escritorio.

Aquí significa aprovechar que una parte del algoritmo corre bien en GPU y otra parte sigue siendo mejor en CPU, o que la CPU puede ir adelantando preparación/postproceso mientras la GPU ejecuta el kernel dominante.

Ejemplos:

- EEQ: la GPU arma la matriz de interacción amortiguada y la CPU ensambla y resuelve el sistema lineal
- D4: la GPU ataca el término dominante O(N²) y la CPU sigue encargándose de partes auxiliares o correctivas
- CPM: la GPU construye la matriz Coulomb y la CPU conserva la iteración de carga

Es un híbrido cooperativo, no una duplicación ciega del cálculo completo en los dos lados.

## Restricciones reales del soporte WASM

El soporte WASM con WebGPU tiene varias limitaciones prácticas que hay que dejar explícitas:

- requiere navegador con WebGPU disponible
- normalmente requiere contexto seguro, es decir HTTPS o `localhost`
- la API expuesta es asíncrona
- el rendimiento depende mucho del navegador, del driver y del sistema operativo
- algunos kernels del backend nativo todavía no están expuestos en la capa WASM WebGPU

## CPU paralela en WASM

También existe la posibilidad de acelerar CPU en WASM con `wasm-bindgen-rayon`, pero ahí hay otra restricción importante.

Para que el modo paralelo CPU funcione en navegador, no basta compilar la feature `parallel`. Además hacen falta capacidades del entorno web para threads/atomics.

En términos prácticos, eso implica:

- `SharedArrayBuffer` disponible
- aislamiento de origen adecuado en el navegador
- soporte de `atomics` y `bulk-memory` en el target/configuración WASM

Si eso no está configurado, el build paralelo WASM falla o el runtime no puede usar el pool de hilos.

## Ruta reproducible para navegador híbrido completo

La ruta reproducible recomendada en este repositorio ya no es invocar `cargo check` desde el root para validar threads en navegador, sino usar el empaquetado canónico de `crates/wasm`.

Paso 1: construir el paquete web con la ruta oficial.

```bash
cd crates/wasm
./build.sh --web-only --web-features "parallel experimental-gpu"
```

Ese script fuerza el target `web`, conserva `snippets` para `wasm-bindgen-rayon` y añade los target features WASM necesarios cuando `parallel` está activo.

Paso 2: servir el frontend con cabeceras COOP y COEP.

Cabeceras mínimas:

- `Cross-Origin-Opener-Policy: same-origin`
- `Cross-Origin-Embedder-Policy: require-corp`

Paso 3: inicializar los dos runtimes en el navegador.

```ts
await init();
await initThreadPool(navigator.hardwareConcurrency ?? 4);
await init_webgpu();
```

Con eso queda cerrado el camino híbrido de punta a punta:

- CPU paralela en navegador vía `wasm-bindgen-rayon`
- GPU vía WebGPU
- modo `hybrid` para repartir trabajo entre ambos en kernels volumétricos y rutas experimentales

Hay un ejemplo mínimo ya preparado en [crates/wasm/examples/vite-webgpu-hybrid](../../crates/wasm/examples/vite-webgpu-hybrid/README.md).

## Resumen práctico actualizado

La situación actual queda así:

- escritorio nativo: backend GPU principal, síncrono, con fallback a CPU
- WASM: backend WebGPU específico, asíncrono, con fallback a CPU
- híbrido en WASM: cooperativo y limitado por la naturaleza del navegador
- CPU paralela en WASM: posible, pero depende de atomics, bulk-memory y configuración del entorno web

## Consecuencia práctica para usuarios

La respuesta correcta a “¿se puede correr en cualquier adaptador gráfico?” es:

Se puede correr en una gran variedad de adaptadores y marcas, pero no hay garantía universal por modelo. La compatibilidad real depende del backend gráfico del sistema, del driver y de los límites del adaptador.

La respuesta correcta a “¿sirve en NVIDIA, Intel o AMD?” es:

Sí, en principio sí. El diseño no está atado a una marca y puede funcionar con NVIDIA, Intel, AMD y Apple, siempre que el sistema exponga un backend soportado por `wgpu` y la inicialización del dispositivo tenga éxito.

## Comportamiento de seguridad

Si la GPU no está disponible, `sci-form` no debería quedar inutilizable por eso. El diseño actual hace fallback a CPU cuando:

- falta la feature `experimental-gpu`
- no se encuentra adaptador
- falla la creación del dispositivo GPU

Eso reduce el riesgo operativo: la GPU acelera, pero no es obligatoria para que el código siga funcionando.

## Recomendación operativa

Si quieres usar esta ruta en producción o en validación seria, conviene asumir esta regla:

Validar por combinación de sistema operativo + driver + familia de GPU, no solo por marca.

Una matriz razonable sería:

- macOS + Apple Silicon
- Windows + NVIDIA
- Windows + Intel
- Linux + AMD o NVIDIA con Vulkan

Con eso se verifica la portabilidad real del backend, no solo la intención multiplataforma del código.