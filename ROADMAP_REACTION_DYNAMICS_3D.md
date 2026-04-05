# Roadmap: Simulación Realista de Dinámicas de Reacción en 3D

> **Objetivo**: Corregir la orientación solo en eje X del sistema actual de reacciones y evolucionar hacia una simulación de dinámicas de reacción que represente fielmente la realidad experimental: movimiento en los 3 ejes, orientaciones moleculares guiadas por orbitales, repulsiones electrostáticas, efectos estéricos, y caminos de mínima energía calculados frame-a-frame.

---

## Diagnóstico del Estado Actual

### Problemas Identificados

#### 1. Orientación forzada en eje X (BUG CRÍTICO)

**En Rust** (`dynamics.rs`, `build_reaction_complex()`):
- Molécula 0 se rota para que su átomo reactivo apunte a `[1, 0, 0]` (+X)
- Molécula 1 se rota para que su átomo reactivo apunte a `[-1, 0, 0]` (-X)
- Las moléculas extra (≥2) se apilan más en +X: `extra += 4.0`
- Resultado: **colisiones frontales siempre a lo largo del eje X**

**En JS** (`useReactionCalculation.ts`, V-path fallback):
- `rotateToAlign()` usa `dir = [1,0,0]` y `dir = [-1,0,0]` explícitamente
- Las moléculas siempre se mueven en ±X

**Consecuencia**: Reacciones como ataques π laterales (Diels-Alder, adiciones electrofílicas), ataques nucleofílicos a ángulos específicos (Bürgi-Dunitz ~107°), o eliminaciones E2 anti-periplanares **no pueden representarse correctamente**.

#### 2. NEB simplificado (no CI-NEB)

El NEB actual es una versión simplificada:
- **Sin climbing image** → el TS no está en el punto silla real
- **Sin proyección tangencial** → las fuerzas de resorte interfieren con las fuerzas químicas reales
- Spring-coupled relaxation mezcla fuerzas de spring con gradientes → imágenes cortan esquinas en la PES

#### 3. Complejo reactivo sin optimización

`build_reaction_complex()` no minimiza la estructura del complejo inicial:
- Los fragmentos se colocan geométricamente pero no se relajan
- La energía de interacción inicial puede ser muy alta o artificial

#### 4. Sin per-frame quantum refinement de la trayectoria

El NEB genera la geometría, pero no hay:
- Optimización de geometría constraint en cada frame
- Recálculo de enlace/ruptura basado en distancias y órdenes de enlace Wiberg/Mayer
- Detección dinámica de formación/ruptura de enlaces

#### 5. Approach/Departure simplista

`slide_molecules()` mueve radialmente respecto al COM global, pero:
- No considera la orientación de orbitales HOMO/LUMO
- No considera estérica (volumen de van der Waals)
- No hay angular sampling

---

## Fundamentos Teóricos para la Implementación

### Teoría del Estado de Transición (TST)
- La velocidad de reacción depende de la energía del complejo activado (punto silla)
- $k = \kappa \frac{k_B T}{h} e^{-\Delta G^\ddagger / RT}$ (ecuación de Eyring)
- El TS es un **máximo en una dirección** (coordenada de reacción) y **mínimo en todas las demás**

### Coordenada Intrínseca de Reacción (IRC)
- Curva paramétrica que conecta dos mínimos pasando por el punto silla
- Sigue el camino de máximo descenso en masa ponderada desde el TS
- Es la proyección 1D de la PES multidimensional que realmente importa

### Superficie de Energía Potencial (PES)
- Sistema de N átomos: 3N-6 grados de libertad internos
- La reacción transcurre en este espacio multidimensional, no solo en X
- Los puntos estacionarios (∂E/∂q_i = 0 para todo i) definen: mínimos (reactivos/productos), puntos silla (TS)

### Ángulo de Ataque Bürgi-Dunitz
- Ataques nucleofílicos al C=O ocurren a ~107° respecto al plano del carbonilo
- **No se puede representar con approach solo en eje X**

### Reglas de Woodward-Hoffmann / FMO
- Las reacciones pericíclicas dependen de la simetría de orbitales frontera
- La orientación relativa de HOMO/LUMO determina si la reacción es suprafacial/antarafacial
- Requiere conocer la distribución espacial de orbitales para determinar el ángulo de approach

---

## Arquitectura de la Solución

### Principios Rectores

1. **Product-Guided Orientation** (ya parcialmente implementado) — usar la geometría del producto para guiar la orientación del reactivo
2. **Orbital-Guided Approach** (NUEVO) — el ángulo de approach debe ser dictado por la densidad orbital HOMO/LUMO
3. **Full 3D NEB** (ya funciona en coordenadas internas) — pero necesita mejoras
4. **Per-frame Energy/Property Calculation** (parcialmente implementado) — necesita ser completo
5. **Climbing-Image NEB** (NUEVO) — para localizar el TS real

---

## Fases de Implementación

### Fase 0: Corrección Urgente — Eliminar Forzado a Eje X
**Prioridad: CRÍTICA** | **Complejidad: Media** | **Archivos: `dynamics.rs`, frontend `useReactionCalculation.ts`**

#### 0.1 Rust: Reemplazar orientación ±X por orientación basada en producto

**Estado actual** en `build_reaction_complex()`:
```rust
// Orienta mol 0 reactive atom → +X
rotate_to_align(&mut mols[0].0, [rx/rl, ry/rl, rz/rl], [1.0, 0.0, 0.0]);
// Orienta mol 1 reactive atom → -X
rotate_to_align(&mut mols[1].0, [rx/rl, ry/rl, rz/rl], [-1.0, 0.0, 0.0]);
```

**Solución propuesta**: Reemplazar con orientación basada en la dirección inter-fragmento derivada de la geometría del producto. `build_product_guided_reactant_complex()` ya hace algo parecido con Kabsch, pero `build_reaction_complex()` (que se usa como fallback y para moléculas extra) sigue forzando a X.

**Tareas:**
- [ ] Deprecar `build_reaction_complex()` como assembly primario
- [ ] Hacer que `build_product_guided_reactant_complex()` sea el único camino
- [ ] Añadir rotación Kabsch para >2 fragmentos (actualmente solo funciona bien para 2)
- [ ] Para moléculas extra (≥2), posicionar según sus átomos equivalentes en el producto, no apiladas en +X

#### 0.2 Frontend: Eliminar V-path X-axis fallback

**JS V-path** en `useReactionCalculation.ts`:
```typescript
rotateToAlign(a, closestIdx_a, [1, 0, 0]);
rotateToAlign(b, closestIdx_b, [-1, 0, 0]);
```

**Solución**: El V-path JS no debería existir como camino principal. Si no hay productos (single-reactant mode), usar el Rust NEB igualmente con un product guess, o eliminarlo para evitar confusión.

**Tareas:**
- [ ] Eliminar `rotateToAlign` con vectores ±X hardcodeados
- [ ] Si se mantiene V-path como fallback rápido, usar la dirección COM→COM entre fragmentos como eje de approach, no X
- [ ] Mejor: enviar siempre al backend Rust para que use product-guided orientation

#### 0.3 Test de regresión

- [ ] Test: SN2 Br⁻ + CH₃Cl → debe mostrar approach desde el lado opuesto al Cl, no solo en X
- [ ] Test: Diels-Alder → el dieno debe acercarse al dienófilo de manera suprafacial
- [ ] Test: H₂ + C₂H₂ → H₂ debe acercarse lateralmente al triple enlace

---

### Fase 1: Climbing-Image NEB (CI-NEB) con Tangent Projection
**Prioridad: Alta** | **Complejidad: Alta** | **Archivo: `dynamics.rs` (nueva función)**

El NEB actual mezcla fuerzas spring con gradientes → imágenes no encuentran el TS real.

#### 1.1 Implementar tangent vector (mejorado de Henkelman-Jónsson)

```
τ_i = R_{i+1} - R_{i-1}  (tangente simple)
// Mejorado: usar pesos basados en energía
τ_i = τ^+ si E_{i+1} > E_{i-1}
      τ^- si E_{i-1} > E_{i+1}
      interpolación si E_i es máximo/mínimo local
```

#### 1.2 Proyectar fuerzas perpendiculares al tangente

```
F_i^⊥ = F_i^real - (F_i^real · τ_i_hat) τ_i_hat   // fuerza real perpendicular
F_i^spring = k (|R_{i+1} - R_i| - |R_i - R_{i-1}|) τ_i_hat  // spring en tangente
F_i^NEB = F_i^⊥ + F_i^spring   // imágenes intermedias
```

#### 1.3 Climbing image

La imagen de máxima energía no usa spring forces sino que **sube** a lo largo del tangente:
```
F_CI = F^real - 2(F^real · τ_hat) τ_hat
```

Esto hace que la climbing image converja al **punto silla real** de la PES.

**Tareas:**
- [ ] `fn compute_ci_neb_path(...)` — nueva función en `dynamics.rs`
- [ ] Tangente mejorada con pesos de energía
- [ ] Proyección de fuerzas (⊥spring, ∥real)
- [ ] Climbing image (la de máxima energía usa F_CI)
- [ ] Convergencia: criterio basado en `max |F_i^⊥|` < threshold
- [ ] Exponer vía WASM: `compute_ci_neb_path_with_method()`
- [ ] Mantener el NEB simplificado como fallback rápido

---

### Fase 2: Orientación Guiada por Orbitales (Orbital-Guided Approach)
**Prioridad: Alta** | **Complejidad: Alta** | **Archivo: nuevo `src/dynamics_orbitals.rs`**

Usar la distribución espacial del HOMO/LUMO para determinar la dirección de approach óptima.

#### 2.1 Orbital direction extraction

Dado un cálculo EHT/PM3/xTB, extraer la dirección principal del HOMO (nucleófilo) y LUMO (electrófilo):

```rust
fn orbital_principal_direction(
    elements: &[u8],
    positions: &[[f64; 3]],
    mo_coefficients: &[f64],  // del EHT/PM3
    mo_index: usize,  // HOMO o LUMO
) -> [f64; 3]
```

**Método**: 
1. Calcular contribuciones AO → MO para cada átomo
2. Ponderar posiciones atómicas por |c_i|² de las contribuciones orbitales
3. La dirección principal es el eigenvector dominante de la matriz de covarianza ponderada
4. Esto da el "lóbulo principal" del orbital

#### 2.2 FMO-guided approach direction

Para reacción nucleófilo + electrófilo:
- Dirección de approach = línea que conecta centro de HOMO(nucleófilo) → centro de LUMO(electrófilo)
- Esto automáticamente produce ángulos Bürgi-Dunitz para C=O, approach lateral para π systems, etc.

#### 2.3 Integration con build complex

```rust
fn build_orbital_guided_complex(
    r_confs: &[ConformerResult],
    reactive_dist: f64,
    homo_direction: [f64; 3],  // del nucleófilo
    lumo_direction: [f64; 3],  // del electrófilo
) -> (Vec<f64>, Vec<u8>)
```

**Tareas:**
- [ ] Extraer dirección principal de orbital MO
- [ ] Identificar HOMO/LUMO automáticamente
- [ ] Calcular dirección de approach basada en FMO
- [ ] Rotar fragmentos para alinear HOMO→LUMO
- [ ] Fallback a product-guided si orbital info no disponible

---

### Fase 3: Optimización de Geometría del Complejo Reactivo
**Prioridad: Media** | **Complejidad: Media** | **Archivo: `dynamics.rs`**

El complejo reactivo debe optimizarse parcialmente antes del NEB para evitar energías artificiales.

#### 3.1 Constrained optimization del complejo inicial

```rust
fn optimize_reactive_complex(
    smiles: &str,
    coords: &[f64],
    elements: &[u8],
    frozen_atoms: &[usize],  // átomos a no mover
    max_steps: usize,
    method: NebBackend,
) -> Result<Vec<f64>, String>
```

- Usar L-BFGS o steepest descent con restricciones
- Congelar el ángulo/distancia entre los centros reactivos
- Relajar todo lo demás (conformación interna de cada fragmento)

#### 3.2 Integration con el pipeline

Antes de NEB:
1. Embed fragmentos → 3D coords
2. Posicionar con product-guided/orbital-guided
3. **NUEVO**: Optimizar complejo con constraints → complex relajado
4. NEB entre complex relajado y product

**Tareas:**
- [ ] Implementar constrained geometry optimization
- [ ] Freeze mask para átomos
- [ ] L-BFGS con projected gradients
- [ ] Integrar en `compute_reaction_dynamics()` pipeline

---

### Fase 4: Per-Frame Property Calculation Completo
**Prioridad: Media** | **Complejidad: Media** | **Archivos: `dynamics.rs`, WASM `reaction.rs`, frontend**

Actualmente el frontend hace un segundo pass con `computeReactionProfile()`, pero es parcial.

#### 4.1 Rust-side per-frame properties

Para cada frame del NEB, calcular en Rust (no solo energía):

```rust
pub struct ReactionFrameProperties {
    pub energy_kcal_mol: f64,
    pub mulliken_charges: Vec<f64>,
    pub wiberg_bond_orders: Vec<(usize, usize, f64)>,
    pub homo_energy: f64,
    pub lumo_energy: f64,
    pub gap: f64,
    pub dipole: [f64; 3],
    pub dipole_magnitude: f64,
}
```

#### 4.2 Bond formation/breaking detection

```rust
fn detect_bond_changes(
    frame_bonds: &[(usize, usize, f64)],  // Wiberg orders
    ref_bonds: &[(usize, usize, f64)],    // reactant bonds
    threshold: f64,  // 0.3 Wiberg order
) -> BondChangeReport
```

- Compare bond orders between frames to detect:
  - Bond breaking: order decreases below threshold
  - Bond forming: order increases above threshold
  - Bond order change: significant change without break/form

#### 4.3 Frontend streaming

- Cambiar de "run NEB → then property pass" a "NEB returns frames with properties"
- Eliminar el segundo pass en el frontend
- Cada frame ya trae charges, bond orders, HOMO/LUMO

**Tareas:**
- [ ] Struct `ReactionFrameProperties` en Rust
- [ ] Calcular properties para cada NEB image en el pipeline
- [ ] Bond change detection basado en Wiberg/Mayer orders
- [ ] Integrar en `ReactionDynamicsFrame`
- [ ] Frontend: eliminar segundo pass, usar propiedades embebidas

---

### Fase 5: IRC de Máximo Descenso (Steepest Descent from TS)
**Prioridad: Media** | **Complejidad: Alta** | **Archivo: nuevo módulo `src/irc.rs`**

Una vez localizado el TS con CI-NEB, trazar el IRC real.

#### 5.1 Mass-weighted steepest descent

```rust
pub fn compute_irc(
    elements: &[u8],
    ts_coords: &[f64],
    ts_gradient: &[f64],
    method: NebBackend,
    step_size: f64,
    max_steps: usize,
    direction: IrcDirection,  // Forward | Backward | Both
) -> Result<IrcResult, String>
```

**Algoritmo**:
1. Desde el TS, dar un paso en la dirección del eigenvector negativo de la Hessiana
2. Seguir el gradiente en coordenadas mass-weighted
3. Integrar con Euler o Gonzalez-Schlegel second-order:
   ```
   q_{n+1} = q_n - Δs * ∇E(q_n) / |∇E(q_n)|  (mass-weighted)
   ```
4. Terminar cuando la energía suba o el gradiente sea < threshold

#### 5.2 Hessiana numérica en el TS

- Calcular la Hessiana por finite differences en el TS
- Diagonalizar → el eigenvector con eigenvalue negativo es la dirección de reacción
- Los otros eigenvectors son modos de vibración del TS

**Tareas:**
- [ ] Módulo `src/irc.rs`
- [ ] Mass-weighted steepest descent integration
- [ ] Numerical Hessian at TS
- [ ] Dirección de reacción desde eigenvector negativo
- [ ] Gonzalez-Schlegel corrector (opcional, para precisión)
- [ ] WASM binding

---

### Fase 6: Electrostatics-Aware Approach
**Prioridad: Media** | **Complejidad: Media**

Durante el approach, la orientación debería ajustarse por fuerzas electrostáticas.

#### 6.1 Electrostatic steering

```rust
fn compute_approach_with_electrostatics(
    mol_a: &FragmentData,  // con charges parciales
    mol_b: &FragmentData,
    n_frames: usize,
    far_dist: f64,
    reactive_dist: f64,
) -> Vec<Vec<f64>>  // frames con orientación adaptiva
```

**Método**:
1. Calcular cargas parciales (Gasteiger) de cada fragmento aislado
2. En cada step del approach:
   - Calcular torque electrostático: $\vec{\tau} = \sum_i q_i \vec{r}_i \times \vec{E}(\vec{r}_i)$
   - Rotar fragmentos según el torque (las cargas opuestas se atraen)
   - Avanzar un paso hacia el otro fragmento
3. Resultado: las moléculas rotan naturalmente para minimizar la energía electrostática durante el approach

#### 6.2 Van der Waals pre-screening

- Antes de approach, calcular el volumen vdW de cada fragmento
- Evitar colisiones estéricas: si el approach en una dirección causa overlap vdW, buscar orientación alternativa
- Usar repulsión Lennard-Jones suavizada para guiar

**Tareas:**
- [ ] Torque electrostático para rotación durante approach
- [ ] VdW pre-screening para evitar colisiones estéricas
- [ ] Combinar con orbital-guided approach (Fase 2)
- [ ] Angular sampling: probar múltiples orientaciones, elegir la de menor energía

---

### Fase 7: Multi-Angular Sampling (Monte Carlo sobre orientaciones)
**Prioridad: Baja** | **Complejidad: Alta**

Para reacciones donde el ángulo de approach no es obvio, probar múltiples orientaciones.

#### 7.1 Rotational sampling

```rust
fn sample_approach_orientations(
    mol_a: &FragmentData,
    mol_b: &FragmentData,
    n_samples: usize,
    method: NebBackend,
) -> Vec<(f64, QuaternionRotation)>  // (energy, rotation)
```

1. Generar `n_samples` orientaciones aleatorias uniformes en SO(3) (uniform quaternion sampling)
2. Para cada orientación, colocar fragmentos cara a cara, evaluar energía
3. Tomar las K mejores orientaciones
4. Para cada una, correr un NEB corto
5. Elegir el path de menor barrera de activación

**Tareas:**
- [ ] Uniform quaternion sampling on SO(3)
- [ ] Short NEB for each candidate orientation
- [ ] Select lowest-barrier path
- [ ] Parallel evaluation (feature "parallel")

---

### Fase 8: Visualización Mejorada en Frontend
**Prioridad: Media** | **Complejidad: Media** | **Archivos: frontend `3d-orbitals-simulations`**

#### 8.1 Bond formation/breaking visual feedback

- Cambiar el estilo de visualización de enlaces que están formándose (---) o rompiéndose (···)
- Color de enlace basado en bond order: fuerte=sólido, débil=transparente
- Usar los `wiberg_bond_orders` del per-frame properties

#### 8.2 Charge evolution overlay

- Mostrar cargas parciales cambiando de color por frame
- Mapa de calor sobre átomos: rojo=negativo, azul=positivo
- Animación de redistribución de carga durante la reacción

#### 8.3 Energy profile synchronized

- Marker en el diagrama de energía que sigue al frame actual
- Highlight del TS con color diferente
- Mostrar ΔG‡ y ΔG° directamente en el gráfico

#### 8.4 Orbital evolution

- Para frames clave (approach, TS, departure), mostrar los orbitales HOMO/LUMO
- Isosurface del orbital que cambia durante la reacción
- Esto ya está parcialmente en `ReactionDialog` con `FrameOrbitalData`

---

## Priorización y Dependencias

```
Fase 0 (X-axis fix) ─────────────────────── INMEDIATA
  │
  ├── Fase 1 (CI-NEB) ──────────────────── ALTA PRIORIDAD
  │     │
  │     └── Fase 5 (IRC from TS) ───────── MEDIA
  │
  ├── Fase 2 (Orbital-guided) ──────────── ALTA PRIORIDAD
  │     │
  │     └── Fase 6 (Electrostatics) ────── MEDIA
  │           │
  │           └── Fase 7 (Multi-angular)── BAJA
  │
  ├── Fase 3 (Complex optimization) ────── MEDIA
  │
  ├── Fase 4 (Per-frame properties) ────── MEDIA
  │
  └── Fase 8 (Frontend visualization) ──── MEDIA (paralelo)
```

---

## Capacidades Existentes que se Reutilizan

| Capacidad | Módulo | Estado |
|-----------|--------|--------|
| EHT orbitals + MO coefficients | `src/eht/` | ✅ Funcional |
| PM3 SCF + gradients | `src/pm3/` | ✅ Funcional |
| GFN0/1/2-xTB | `src/xtb/` | ✅ Funcional |
| Wiberg/Mayer bond orders | `compute_bond_orders` | ✅ Funcional |
| Mulliken/Löwdin charges | `compute_population` | ✅ Funcional |
| Gasteiger charges | `compute_charges` | ✅ Funcional |
| Kabsch alignment | `src/alignment/` | ✅ Funcional |
| EHT gradients + geometry opt | `src/eht/gradients.rs` | ✅ Funcional |
| NEB multi-method | `dynamics.rs` | ✅ Simplificado |
| MD Velocity Verlet/Nosé-Hoover | `dynamics.rs` | ✅ Funcional |
| Framework geometry opt (BFGS) | `src/materials/geometry_opt.rs` | ✅ Reutilizable |
| ECFP fingerprints | `src/ecfp/` | ✅ Para similarity checks |
| Dipole moment | `compute_dipole` | ✅ Funcional |
| ESP grid | `compute_esp` | ✅ Para visualization |
| Vibrational analysis | `src/ir/` | ✅ Para Hessiana en TS |
| Product-guided complex | `dynamics.rs` | ✅ Parcial |
| SMIRKS transforms | `src/smirks/` | ✅ Para prediction |

---

## Métricas de Validación

### Para cada fase, validar contra:

1. **SN2 (CH₃Cl + Br⁻)**: Walden inversion, approach desde lado opuesto al leaving group, ángulo ~180°
2. **Diels-Alder (butadieno + etileno)**: Approach suprafacial [4+2], ambas moléculas coplanares
3. **Esterificación (ácido acético + metanol)**: Bürgi-Dunitz angle ~107° al C=O
4. **H₂ + C₂H₂ → C₂H₄**: Approach lateral al triple enlace (en plano π)
5. **Sustitución electrofílica aromática**: Approach perpendicular al plano del anillo

### Métricas cuantitativas:
- **Barrera de activación**: ΔE‡ dentro del 20% del valor experimental/NIST
- **Energía de reacción**: ΔE_rxn con signo correcto
- **Ángulo de approach**: Dentro de ±15° del valor experimental para reacciones conocidas
- **Conservación de energía**: Drift <1% en MD trajectories

---

## Notas Técnicas

### Coordenadas mass-weighted vs cartesianas
- IRC usa coordenadas mass-weighted: $q_i = \sqrt{m_i} \cdot x_i$
- NEB puede funcionar en cartesianas pero converge mejor en mass-weighted
- La Hessiana en mass-weighted da frecuencias vibracionales directamente

### WASM constraints
- Todo debe cruzar la frontera WASM como JSON strings o typed arrays
- No hay threads en WASM → el NEB debe ser secuencial en el hilo principal o usar Web Workers vía `split_worker_tasks`
- Para CI-NEB con muchas imágenes, considerar progress callbacks

### Performance
- CI-NEB con PM3: ~0.5-2s para 30 imágenes (orgánicos pequeños)
- CI-NEB con GFN2: ~5-15s para 30 imágenes
- CI-NEB con HF-3c: ~30-120s para 30 imágenes
- Orbital-guided approach: ~200ms extra (un cálculo EHT)
- Per-frame properties: ~1-3s total para 60 frames con PM3

### Feature flags sugeridas
```toml
[features]
ci-neb = []           # Climbing-image NEB
irc = []              # IRC steepest descent
orbital-approach = []  # Orbital-guided complex assembly
```
