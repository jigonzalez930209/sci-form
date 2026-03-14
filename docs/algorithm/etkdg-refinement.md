# ETKDG Refinement

The **Experimental Torsion Knowledge Distance Geometry** (ETKDG) refinement is what distinguishes this algorithm from plain distance geometry. It uses experimentally observed torsion angle preferences from the **Cambridge Structural Database** (CSD) to guide conformer geometry toward chemically realistic conformations.

## The Key Insight

Plain distance geometry produces geometrically valid structures, but the torsion angles around rotatable bonds are essentially random. The ETKDG approach adds a **torsion preference force field** that biases the conformer toward experimentally observed dihedral angles.

<img src="/svg/etkdg-pipeline.svg" alt="etkdg-pipeline" class="svg-diagram" />

## CSD Torsion Pattern Library

sci-form includes **837 SMARTS patterns** with associated Fourier coefficients, derived from the Cambridge Structural Database analysis by Guba et al.

### Pattern Categories

| Category | Count | Description |
|----------|-------|-------------|
| v2 patterns | 365 | General torsion patterns |
| Macrocycle patterns | 472 | Patterns specific to macrocyclic systems |

### Fourier Representation

Each pattern encodes the preferred torsion angle distribution as a 6-term Fourier series:

$$V(\phi) = \sum_{k=1}^{6} V_k (1 + s_k \cos(k\phi))$$

where:
- $V_k$ is the amplitude for the $k$-th Fourier component
- $s_k \in \{-1, +1\}$ is the sign
- $\phi$ is the dihedral angle

The coefficients $V_k$ are derived by fitting to the observed torsion angle histograms from crystallographic data.

### Pattern Matching Priority

For each rotatable bond, patterns are matched in order:

1. **CSD patterns** — first-match-wins among the 837 SMARTS
2. **Basic knowledge** — fallback rules for common chemical environments

::: tip
Pattern matching uses a first-match-wins strategy. More specific patterns are listed before general ones, so a pattern for "amide C-N" will match before a generic "any C-N" pattern.
:::

## Basic Knowledge Torsion Rules

When no CSD pattern matches a rotatable bond, these rules provide reasonable defaults:

| Environment | Rule | $V_k$ |
|-------------|------|-------|
| Ring bond (4-member) | Flat | $V_2 = 100.0$ |
| Ring bond (5-member) | Flat | $V_2 = 100.0$ |
| Ring bond (6-member) | Flat | $V_2 = 100.0$ |
| Double bond | Planar | $V_2 = 100.0$ |
| Amide C-N | Planar preference | $V_2 = 7.0$ |
| Ester C-O | Planar preference | $V_2 = 7.0$ |
| Aromatic-X | Semi-planar | $V_2 = 5.0$ |
| SP3-SP3 | Staggered | $V_3 = 7.0$ |
| Ether/Amine | Soft staggered | $V_3 = 2.5$ |
| Biaryl | Semi-planar | $V_2 = 5.0$ |

## ETKDG 3D Force Field Components

The complete ETKDG 3D force field combines torsion preferences with structural constraints:

### 1. Torsion Contributions (from CSD or basic knowledge)

For each matched torsion:
$$E_{\text{tors}} = V(\phi) = \sum_{k=1}^{6} V_k(1 + s_k \cos(k\phi))$$

Computed via Chebyshev recurrence for efficiency — only one $\cos$ evaluation per torsion angle.

### 2. UFF Inversions (Out-of-Plane)

For SP2 centers with 3 heavy neighbors:

$$E_{\text{inv}} = K \cdot (1 - \sin Y)$$

where $Y$ is the Wilson out-of-plane angle. The energy is zero when the atom is perfectly planar ($Y = 90°$) and increases as it deviates.

Three permutations are evaluated per improper center, cycling through the neighbor triple:

$$E_{\text{inv}}^{\text{total}} = \frac{K_{\text{base}} \cdot 10}{3} \sum_{p=1}^{3} (1 - \sin Y_p)$$

### 3. Distance Constraints

Maintain bond lengths and angles via flat-bottom potentials:

$$E_{\text{dist}} = \frac{k}{2} \max(0, |d - d_0| - \epsilon)^2$$

| Bond type | $k$ | $\epsilon$ |
|-----------|-----|------------|
| 1-2 (bonds) | 100 | 0.01 Å |
| 1-3 (improper) | 100 | varies |
| Long-range | 10 | varies |

### 4. Linear Angle Constraints

For SP atoms (triple bonds, allenes), maintain 180° angle:

$$E_{\text{angle}} = k \cdot (\theta - 180°)^2$$

## Optimization

The ETKDG 3D force field is minimized with a single BFGS pass:

- **Maximum iterations**: 300
- **No restarts** — this is a refinement step, not a global optimization
- **Early skip**: if initial energy < $10^{-5}$, skip entirely
- **Result**: final 3D coordinates ready for validation

## Example: Butane Torsion

For butane (CCCC), the central C-C bond matches a CSD pattern with staggered preference ($V_3 = 7.0$):

<img src="/svg/etkdg-torsion.svg" alt="etkdg-torsion" class="svg-diagram" />

The torsion force field naturally drives the dihedral toward the anti (180°) or gauche (±60°) conformations, which match experimental observation.
