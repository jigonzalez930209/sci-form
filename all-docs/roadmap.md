# Roadmap de desarrollo de all-docs

## Objetivo

Este documento ordena el trabajo de `all-docs/` desde la base más estable hasta los módulos más experimentales. La idea no es listar carpetas en orden alfabético, sino definir una secuencia de cierre que respete dependencias reales, reduzca retrabajo y permita que cada archivo de módulo concentre sus mejoras propuestas sin duplicación.

## Regla general de ejecución

1. Primero se cierran los documentos que describen el contrato transversal del crate.
2. Después se avanza por la geometría y la construcción molecular, porque esas capas alimentan casi todo lo demás.
3. Luego se completan cargas, superficie, solvatación y propiedades locales.
4. Más tarde se cierran los métodos electrónicos y semiempíricos.
5. La espectroscopía y la respuesta molecular van después de la base electrónica.
6. Materiales, ML e infraestructura quedan en una segunda mitad del recorrido.
7. Los módulos experimentales y de aceleración se dejan para el final, cuando la base estable ya está bien definida.

## Fase 0. Contrato transversal y documentación interna

Esta fase se empieza aquí porque define cómo leer el resto del árbol.

### [resto-src.md](resto-src.md)
- Cerrar primero este archivo, porque resume la superficie pública real del crate y conecta Rust, Python, WASM y CLI.
- Las mejoras deben centrarse en fronteras de capa, regresión transversal, acoplamiento y observabilidad.
- Este documento sirve de referencia para detectar si una mejora en un módulo cambia el contrato global.

### [bin.md](bin.md)
- Cerrar después de `resto-src.md` para separar herramientas internas de API pública.
- Las mejoras deben enfocarse en nomenclatura, salidas consistentes, vigencia de utilidades y conexión con pruebas.
- Si un binario sigue siendo útil, este es el sitio para decidir si debe documentarse, depurarse o retirarse.

## Fase 1. Geometría, topología y preparación molecular

Aquí se resuelven las piezas que sostienen conformación, validación geométrica y manejo estructural.

### [distgeom.md](distgeom.md)
- Prioridad alta porque el bounds matrix y la geometría inicial alimentan conformer, forcefield y alignment.
- Las mejoras deben orientarse a exactitud de distancias, estabilidad numérica, manejo de restricciones y coste en moléculas grandes.

### [alignment.md](alignment.md)
- Va justo después de `distgeom.md` porque depende de coordenadas bien formadas.
- Las mejoras deben reforzar RMSD, Kabsch, rotaciones, comparación de conformeros y manejo de simetrías.

### [forcefield.md](forcefield.md)
- Debe cerrarse en esta fase porque da el mínimo energético que estabiliza el resto del pipeline geométrico.
- Las mejoras deben centrarse en typing, parámetros, gradientes, reducción de recomputación y documentación de límites.

### [rings.md](rings.md)
- Conviene cerrarlo aquí para que la percepción de anillos esté disponible antes de tocar módulos de aromaticidad, force field o smarts.
- Las mejoras deben validar SSSR, tamaños de anillo, casos fusionados y rendimiento en grafos complejos.

### [smarts.md](smarts.md)
- Debe cerrarse junto con rings porque la lógica de patrones y transformaciones se apoya en el grafo químico ya consistente.
- Las mejoras deben concentrarse en parser, matching, transformaciones seguras y diagnóstico de fallos.

### [transport.md](transport.md)
- Se deja en esta fase porque la serialización y el reparto de lotes dependen de que la base molecular ya esté estable.
- Las mejoras deben cubrir packing, chunking, balanceo, compatibilidad de formatos y costes de copia.

### [optimization.md](optimization.md)
- Puede cerrarse aquí o en la fase 2 si se quiere tratar junto con force field, pero su dependencia natural es geométrica.
- Las mejoras deben enfocarse en criterios de parada, robustez, reporte de energía y coherencia con la geometría inicial.

## Fase 2. Cargas, superficie, dipolos y propiedades locales

Cuando la geometría ya es fiable, se puede cerrar la capa de propiedades derivadas.

### [charges.md](charges.md)
- Empieza después de la geometría porque necesita coordenadas coherentes y convenciones de carga claras.
- Las mejoras deben tratar cálculo, estabilidad, cobertura por elemento y trazabilidad de resultados.

### [charges_eeq.md](charges_eeq.md)
- Debe seguir a `charges.md` porque es una variante especializada de la misma familia de problemas.
- Las mejoras deben compararse contra Gasteiger y otras rutas de carga para medir desviaciones reales.

### [dipole.md](dipole.md)
- Se cierra aquí porque el dipolo depende de cargas, posiciones y convenciones de unidades.
- Las mejoras deben validar vector, magnitud, referencia y sensibilidad al origen.

### [esp.md](esp.md)
- Conviene cerrarlo junto con dipolo porque el potencial electrostático comparte la misma base de cargas y geometría.
- Las mejoras deben atender grillas, padding, densidad de muestreo, coste y lectura del mapa.

### [population.md](population.md)
- Va en esta fase porque la población electrónica se apoya en matrices y convenciones ya estabilizadas.
- Las mejoras deben distinguir Mulliken, Löwdin, HOMO/LUMO, consistencia de cargas y comparaciones con referencias.

### [solvation_alpb.md](solvation_alpb.md)
- Debe cerrarse después de cargas y superficie, porque la solvatación necesita ambas capas bien definidas.
- Las mejoras deben separar energía polar y no polar, validar parámetros y marcar supuestos del modelo.

### [surface.md](surface.md)
- Va en el mismo bloque porque la superficie molecular soporta SASA, cavidades y contribuciones a solvatación.
- Las mejoras deben reducir ambigüedad geométrica, mejorar trazabilidad y medir el coste por átomo.

### [dispersion.md](dispersion.md)
- Debe cerrarse en esta fase porque la dispersión entra como corrección energética transversal.
- Las mejoras deben validar términos, escalado, compatibilidad con forcefield y sensibilidad al entorno.

## Fase 3. Estructura electrónica y métodos semiempíricos

Esta fase se apoya en la geometría ya validada y en las propiedades locales ya coherentes.

### [scf.md](scf.md)
- Es una base importante para métodos posteriores, así que conviene cerrarlo antes de PM3, xTB y EHT.
- Las mejoras deben reforzar convergencia, precisión de matrices, estabilidad y reutilización de resultados.

### [hf.md](hf.md)
- Debe ir después de `scf.md` porque comparte el mismo lenguaje de convergencia y matrices electrónicas.
- Las mejoras deben tocar integrales, gradientes, criterio de parada, coste y validez de ocupaciones.

### [pm3.md](pm3.md)
- Conviene cerrarlo una vez estable la base SCF/HF, porque PM3 depende de una implementación electrónica bien controlada.
- Las mejoras deben centrarse en cobertura de elementos, parámetros, energía, cargas y comparación contra referencias.

### [xtb.md](xtb.md)
- Debe venir después de PM3 por dependencia conceptual y por la necesidad de una base estable de energía y gradiente.
- Las mejoras deben cubrir robustez SCC, dispersión, contribuciones repulsivas y validación por familia química.

### [eht.md](eht.md)
- Se cierra en este bloque porque conecta estructura electrónica con orbitales, densidad y visualización.
- Las mejoras deben atender matrices S/H, eigenproblema generalizado, metadata orbital y coste de volúmenes.

### [dos.md](dos.md)
- Cierra esta fase porque depende de orbitales, energías o métodos electrónicos ya disponibles.
- Las mejoras deben centrarse en muestreo energético, smearing, resolución, método seleccionado y coste para moléculas grandes.

## Fase 4. Espectroscopía y respuesta molecular

Aquí ya se trabaja con resultados electrónicos y vibracionales suficientemente estables.

### [spectroscopy.md](spectroscopy.md)
- Este archivo debe cerrar primero en la fase espectroscópica porque actúa como paraguas de la familia.
- Las mejoras deben ordenar qué métodos son núcleo, cuáles son derivados y cómo se conectan con IR, UV-Vis y NMR.

### [ir.md](ir.md)
- Va después de spectroscopy porque necesita vibraciones, análisis térmico y tratamiento de intensidades.
- Las mejoras deben validar frecuencias, broadening, asignación de picos y reproducibilidad.

### [nmr.md](nmr.md)
- Debe cerrarse al final del bloque espectroscópico porque depende de geometría, cargas, acoplamientos y, a veces, ensembles.
- Las mejoras deben reforzar HOSE, shifts, J-couplings, promedio por conformeros y comparación con datos de referencia.

## Fase 5. Materiales y predicción

Este bloque ya puede cerrarse cuando la base geométrica y electrónica esté suficientemente madura.

### [materials.md](materials.md)
- Debe venir antes de cualquier expansión experimental de materiales, porque define celda, grupos espaciales y ensamblaje.
- Las mejoras deben centrarse en validación cristalográfica, separación de responsabilidades, reducción de coste y mensajes de error.

### [ml.md](ml.md)
- Se cierra después de que descriptors, cargas y geometría estén estables, porque depende de una representación consistente.
- Las mejoras deben atacar calidad de features, validación cruzada, interpretabilidad, cobertura y estabilidad de entrenamiento.

## Fase 6. Módulos experimentales y aceleración

Esta fase se deja para el final porque necesita una base estable para justificar complejidad adicional.

### [ani.md](ani.md)
- Va aquí porque es un método de energía y fuerzas que debe compararse con la base electrónica ya consolidada.
- Las mejoras deben cerrar energía vs fuerzas, cobertura por entorno químico, eficiencia de AEV y coherencia de tolerancias.

### [alpha.md](alpha.md)
- Se deja hacia el final porque agrupa líneas tempranas de investigación que no deben bloquear la base estable.
- Las mejoras deben cerrarse submódulo por submódulo: `cga`, `gsm` y `sdr`.
- Orden sugerido dentro de alpha: primero `cga`, luego `gsm`, al final `sdr`, porque así se validan primero las piezas geométricas más cercanas al resto del crate.

### [beta.md](beta.md)
- También se deja para el final porque agrupa técnicas numéricas avanzadas que conviene justificar con métricas claras.
- Las mejoras deben cerrarse submódulo por submódulo: `kpm`, `mbh`, `cpm`, `rand_nla` y `riemannian`.
- Orden sugerido dentro de beta: primero `kpm`, luego `mbh` y `cpm`, después `rand_nla`, y al final `riemannian`, porque ese recorrido va de aproximaciones más directas a optimización más especializada.

### [gpu.md](gpu.md)
- Debe cerrarse al final porque solo tiene sentido cuando ya existen kernels o rutas que merezcan aceleración real.
- Las mejoras deben definir prioridad de kernels, diagnóstico del backend, coste CPU↔GPU y límites de hardware.

## Fase 7. Cierre de consistencia y revisión final

Esta fase no añade nuevos módulos; valida que todo lo anterior quedó homogéneo.

### Lista de chequeo final
- Confirmar que cada archivo de módulo tenga resumen, archivos, superficie por target, mejoras y observaciones sin duplicación.
- Confirmar que los módulos estables no dependan de observaciones experimentales para justificar su API.
- Confirmar que `alpha.md`, `beta.md` y `gpu.md` no introduzcan requisitos que contradigan el resto del árbol.
- Confirmar que los documentos de base (`resto-src.md`, `bin.md`, `roadmap.md`) sigan siendo consistentes con el contrato público real.

## Orden recomendado de lectura y cierre

Si el objetivo es trabajar archivo por archivo sin perder dependencias, este es el orden recomendado:

1. `resto-src.md`
2. `bin.md`
3. `distgeom.md`
4. `alignment.md`
5. `forcefield.md`
6. `rings.md`
7. `smarts.md`
8. `transport.md`
9. `optimization.md`
10. `charges.md`
11. `charges_eeq.md`
12. `dipole.md`
13. `esp.md`
14. `population.md`
15. `solvation_alpb.md`
16. `surface.md`
17. `dispersion.md`
18. `scf.md`
19. `hf.md`
20. `pm3.md`
21. `xtb.md`
22. `eht.md`
23. `dos.md`
24. `spectroscopy.md`
25. `ir.md`
26. `nmr.md`
27. `materials.md`
28. `ml.md`
29. `ani.md`
30. `alpha.md`
31. `beta.md`
32. `gpu.md`

## Criterio de terminación

El trabajo de `all-docs/` termina cuando:

- Cada archivo describe su módulo sin bloques repetidos.
- Cada archivo deja claro qué mejora toca primero y cuál depende de la base anterior.
- `alpha.md` y `beta.md` mantienen su enfoque por submódulo, sin mezclar prioridades.
- Los documentos de cierre (`resto-src.md`, `bin.md`, `roadmap.md`) ya no contradicen la superficie pública real del proyecto.
- El orden de ejecución queda suficientemente claro para que una siguiente ronda de mejoras pueda avanzar sin rehacer la fase anterior.

## Resolución de bugs y mejoras implementadas

Este documento es de planificación, no de módulo. No tiene hallazgos de código asociados.

### Resumen de la ronda de resolución
- **17 bugs corregidos** en código fuente (de 25 encontrados en la revisión profunda).
- **8 clasificados** como falso positivo o decisión de diseño.
- **565/565 tests pasan** tras las correcciones.
- **Todos los 33 archivos all-docs** actualizados con sección "Resolución de bugs y mejoras implementadas".
