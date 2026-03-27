# TODO consolidado desde all-docs

Fecha de revision: 2026-03-26

## Alcance y criterio

Este documento resume la revision completa de `all-docs/`.

- Solo se listan bugs no resueltos, limitaciones conocidas, roadmaps y mejoras propuestas que siguen pendientes.
- No se repiten items marcados en los propios documentos como `RESUELTO`, `FALSO POSITIVO` o `DECISION DE DISENO`.
- Cuando un archivo no tiene deuda funcional clara, queda marcado como cerrado o con mejoras menores de validacion/documentacion.

## Prioridades sugeridas

1. `scf.md`: soporte open-shell real (`UHF`/`ROHF`), porque hoy el solver sigue siendo `RHF-only`.
2. `hf.md`: bajar coste y deuda metodologica en `CIS`/`CISD` (`AO->MO`, coupling incompleto, memoria GPU).
3. `pm3.md`: completar parametrizacion PM3/PM3(tm) mas alla del solver ya corregido.
4. `xtb.md`: cerrar deuda de `GFN0` y `GFN2`, mas correccion completa de `D3-BJ C8`.
5. `nmr.md`: aromaticidad real, `4J/5J`, ampliacion de HOSE database y constantes de Karplus por entorno.
6. `materials.md`: completar operaciones para los `230 space groups` y anadir `CIF import/export`.
7. `gpu.md`: fallback CPU transparente, tiling/sparse batching y tuning por capacidades del dispositivo.
8. `smarts.md`: `recursive SMARTS` (`$()`) y `SMIRKS` multicomponente.
9. `charges_eeq.md`: ampliar parametros a `Z=86` y mejorar damping.
10. `spectroscopy.md`: prescreening real para evitar el cuello de botella `O(N^4)`.

## Estado por archivo

- `alignment.md`: sin bugs abiertos relevantes; documento funcionalmente cerrado.
- `alpha.md`: sigue experimental; pendiente validacion cuantitativa por submodulo, benchmarking frente a rutas estables y criterios claros de promocion.
- `ani.md`: sin bugs abiertos relevantes; queda como modulo estable con deuda baja.
- `beta.md`: sigue experimental; faltan tests de integracion mas representativos, validacion cuantitativa por submodulo y plan de promocion a estable.
- `bin.md`: sin bugs criticos, pero queda una limitacion conocida en el CLI batch: la entrada por `stdin` se materializa completa antes de procesar; streaming real sigue pendiente.
- `charges.md`: sin bugs abiertos criticos; quedan mejoras de cobertura por elemento, posible DIIS/Broyden y mejor telemetria de iteraciones.
- `charges_eeq.md`: pendientes la expansion de parametros hasta `Z=86`, un damping gaussiano menos naive, `cEEQ` y limpiar restos de experimentalidad en nombres/superficies externas.
- `dipole.md`: sin bugs abiertos; solo mejoras menores de validacion y documentacion.
- `dispersion.md`: pendiente completar la parametrizacion/tablas `D4 C6` hasta cobertura amplia y mantener visible la descomposicion energetica cuando se combine con otros metodos.
- `distgeom.md`: pendiente detectar embeddings 2D o demasiado planos y re-perturbar; tambien conviene reforzar deteccion temprana de bounds contradictorios.
- `dos.md`: sin bugs abiertos relevantes; documento funcionalmente cerrado.
- `eht.md`: sigue pendiente la limitacion de usar tabla `3d STO-3G` para orbitales `4d/5d`, junto con la revision de parametros provisionales en TMs mas alla del camino actual.
- `esp.md`: pendiente hacer configurable la distancia minima dura (`0.01 A`) y endurecer validacion de `spacing`/`padding`.
- `forcefield.md`: `MMFF94` sigue sin cobertura completa para metales de transicion; tambien queda pendiente un detector aromatico independiente y mejoras como gradientes analiticos y breakdown por terminos.
- `gpu.md`: pendientes un fallback CPU transparente, tiling/sparse batching para `MMFF94`, tuning dinamico de batch/workgroup segun el device y mensajes utiles para hardware con poca VRAM.
- `hf.md`: pendientes una transformacion intermedia `AO->MO` para `CIS`, completar `CISD singles-doubles coupling`, reducir la ruta `GPU ERI O(N^4)` en memoria y mejorar momentos de transicion mas alla de la aproximacion de monopolo.
- `ir.md`: sin bugs criticos abiertos; quedan tareas de validacion mas amplia contra espectros experimentales y posibles mejoras de broadening/asignacion.
- `materials.md`: pendientes completar operaciones de simetria para los `230 space groups`, `CIF import/export`, normalizacion automatica de coordenadas fraccionales y expansion mas alla de topologias homogeneas.
- `ml.md`: sin bugs abiertos relevantes; deuda baja y mayormente de evolucion de calidad, no de correctitud.
- `nmr.md`: pendientes aromaticidad real con `SSSR/Huckel`, `4J/5J` y otros long-range couplings, ampliacion fuerte de `HOSE database`, constantes de Karplus mas especificas y salida mas alla del supuesto `first-order`.
- `optimization.md`: sin bugs criticos abiertos; sigue pendiente unificar optimizadores entre rutas (`EHT`/`materials`), anadir `variable-cell BFGS` y seguir endureciendo la estrategia de paso/convergencia.
- `pm3.md`: pendiente completar la parte metodologica del modelo (`resonance integrals`, `nuclear repulsion`, referencias atomicas para `heat of formation`) y dejar bien delimitado `PM3` vs `PM3(tm)`.
- `population.md`: sin bug abierto critico; siguen pendientes ampliar `valence_electrons` hacia `Z=86`, emitir warning cuando falte cobertura y exponer `homo/lumo/gap` en el resultado publico.
- `resto-src.md`: la limitacion principal pendiente es que los indices de Fukui siguen apoyandose en cargas Mulliken y por tanto son basis-set-dependent; si se quiere mas robustez habria que evaluar `Hirshfeld`/`CM5`.
- `rings.md`: sin bugs abiertos relevantes; solo quedan mejoras marginales para casos fused muy complejos.
- `roadmap.md`: no tiene bugs propios; la tarea pendiente es mantener el orden de cierre y la coherencia con la superficie publica real del proyecto.
- `scf.md`: pendiente soporte open-shell real (`UHF`/`ROHF`); hoy el solver compartido sigue siendo `RHF-only`.
- `smarts.md`: pendientes `recursive SMARTS` (`$()`) y `component-level grouping` en `SMIRKS`; tambien conviene endurecer validacion de atom maps antes de aplicar transformaciones.
- `solvation_alpb.md`: pendientes la unificacion con el modulo principal de `SASA`, una separacion mas clara entre `GB` y `ALPB`, y futuras alternativas mas precisas si se amplia el modelo.
- `spectroscopy.md`: pendiente anadir prescreening real (`Cauchy-Schwarz` o similar) para reducir el coste `O(N^4)` en `sTDA`; como mejora secundaria, clasificar excitaciones `dark`/`bright` de forma mas explicita.
- `surface.md`: pendiente documentar/publicar de forma clara la tabla de `fallback radii` y valorar una opcion tipo `Connolly surface`.
- `transport.md`: sin bugs abiertos relevantes; documento funcionalmente cerrado.
- `xtb.md`: persisten simplificaciones de modelo en `GFN0` y `GFN2`; sigue pendiente la correccion completa de `D3-BJ C8`, mas refinamientos en gamma shell-resolved, `SCC damping`, criterio de convergencia de cargas y repulsion dependiente de coordinacion.

## Bloques que hoy se pueden considerar cerrados o de baja prioridad

- Cerrados funcionalmente: `alignment.md`, `ani.md`, `dipole.md`, `dos.md`, `ml.md`, `rings.md`, `transport.md`.
- Con limitacion puntual pero sin bug critico: `bin.md` por el batch via `stdin` sin streaming real.
- Con deuda baja o de validacion, pero sin bug abierto critico: `charges.md`, `ir.md`, `optimization.md`, `population.md`, `roadmap.md`, `surface.md`.
- Experimentales sin bug critico, pero todavia no listos para promocion automatica: `alpha.md`, `beta.md`, `gpu.md`.

## Resumen ejecutivo

La mayor deuda pendiente no esta en bugs de correctitud basica ya detectados, sino en cuatro frentes:

1. Metodos electronicos y semiempiricos aun simplificados (`SCF`, `HF`, `PM3`, `xTB`, `EHT`).
2. Cobertura quimica incompleta en algunos modulos (`EEQ`, `MMFF94`, `materials`, `population`).
3. Lenguajes/modelos aun parciales (`SMARTS/SMIRKS`, `NMR`, `spectroscopy`).
4. Features experimentales que requieren validacion y endurecimiento antes de promoverse (`alpha`, `beta`, `gpu`).

Si hubiera que ordenar la siguiente ronda de trabajo, empezaria por `scf.md`, `hf.md`, `pm3.md`, `xtb.md`, `nmr.md`, `materials.md` y `gpu.md`.