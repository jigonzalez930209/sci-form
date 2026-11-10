# GCMC adsorption in porous materials

The first v2.x materials workflow adds **grand-canonical Monte Carlo (GCMC)**
adsorption on a rigid periodic framework. It answers the question the existing
materials module could not answer: how much gas does a framework take up at a
given pressure and temperature?

> **Status.** Implemented in `feat/wave-a` for rigid Lennard-Jones frameworks
> with optional charges, Ewald electrostatics and multi-site guests. Framework
> flexibility and bonded MM terms remain follow-ups. See the
> [v2.x roadmap](/v2-roadmap).

## Grand-canonical ensemble

In the grand-canonical ensemble the number of adsorbed molecules fluctuates.
The probability of a configuration with `N` guests is proportional to

$$
P(N) \propto \frac{z^N}{N!}\,
e^{-\beta U(N)},
$$

where `z` is the activity, `U` the potential energy and `β = 1/(kBT)`. The
first implementation uses the ideal-gas relation

$$
zV = \frac{PV}{k_BT}
$$

for the insertion/deletion factors, with pressure in Pa and volume in m³.

## Insertion and deletion

For insertion,

$$
P_{\mathrm{ins}} =
\min\left[
1,\,
e^{-\beta\Delta U}
\frac{PV}{k_BT(N+1)}
\right].
$$

For deletion,

$$
P_{\mathrm{del}} =
\min\left[
1,\,
e^{-\beta\Delta U}
\frac{Nk_BT}{PV}
\right].
$$

The implementation also supports random translations. Every move is accepted
with the Metropolis rule, and the framework is kept rigid.

## Energy model

The host–guest and guest–guest interactions use the Lennard-Jones potential

$$
U_{ij}(r) =
4\epsilon_{ij}
\left[
\left(\frac{\sigma_{ij}}{r}\right)^{12}
-
\left(\frac{\sigma_{ij}}{r}\right)^{6}
\right],
$$

with Lorentz–Berthelot mixing

$$
\sigma_{ij} = \frac{\sigma_i+\sigma_j}{2},
\qquad
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j}.
$$

Periodic boundaries use the existing `UnitCell::minimum_image_distance`, which
searches the neighbouring lattice images for skewed cells. A spherical cutoff
is applied to both host–guest and guest–guest pairs.

## API

```rust
use sci_form::gcmc::{
    GcmcConfig, GcmcFramework, GuestModel, LennardJonesAtom,
};
use sci_form::compute_gcmc_isotherm;

let framework = GcmcFramework {
    cell: sci_form::materials::UnitCell::cubic(20.0),
    atoms: vec![LennardJonesAtom {
        position: [10.0, 10.0, 10.0],
        epsilon_kcal_mol: 1.0,
        sigma_angstrom: 3.0,
    }],
};

let guest = GuestModel {
    epsilon_kcal_mol: 0.2,
    sigma_angstrom: 3.0,
};

let config = GcmcConfig {
    temperature_k: 298.15,
    pressures_bar: vec![0.1, 1.0, 10.0],
    ..GcmcConfig::default()
};

let isotherm = compute_gcmc_isotherm(&framework, &guest, &config)?;
for point in isotherm.points {
    println!("{:.2} bar -> {:.3} molecules/cell", point.pressure_bar, point.mean_loading);
}
```

### Input

| Type | Meaning |
|---|---|
| `GcmcFramework` | Unit cell plus host atoms with `epsilon`, `sigma` and optional `charge` |
| `GuestModel` | Single-site or multi-site guest, with optional charges |
| `GcmcConfig` | Temperature, pressure points, move counts, seed, cutoff |

### Output

| Type | Meaning |
|---|---|
| `GuestMolecule` | Centre plus rotation matrix for a multi-site guest |
| `AdsorptionPoint` | Pressure, mean loading, standard error, sample count |
| `AdsorptionIsotherm` | Temperature, guest model, points and notes |

Loading is **molecules per unit cell**. Divide by the cell volume to obtain a
number density, or convert to mmol/g once the framework mass is known.

## Charges and Ewald

`LennardJonesAtom` and `GuestSite` carry a point charge. When
`GcmcConfig::electrostatics` is `Ewald`, the total electrostatic energy is

$$
E_{\mathrm{el}} =
\frac{1}{2}\sum_{i\ne j} \frac{q_iq_j\,\operatorname{erfc}(\alpha r_{ij})}{r_{ij}}
+
\frac{2\pi}{V}\sum_{\mathbf{k}\ne0}
\frac{e^{-k^2/(4\alpha^2)}}{k^2}|S(\mathbf{k})|^2
-
\frac{\alpha}{\sqrt{\pi}}\sum_i q_i^2
-
\frac{\pi}{2V\alpha^2}\left(\sum_i q_i\right)^2.
$$

The last term is the neutralizing background for a non-neutral cell. The
reciprocal sum uses the reciprocal lattice of the existing `UnitCell`. The
Ewald parameters are `alpha_angstrom_inv`, `k_max` and
`real_cutoff_angstrom`.

## Multi-site guests

A guest can define several interaction sites:

```rust
let co2 = GuestModel {
    sites: vec![
        GuestSite {
            offset_angstrom: [-1.16, 0.0, 0.0],
            epsilon_kcal_mol: 0.2,
            sigma_angstrom: 3.0,
            charge: -0.35,
        },
        GuestSite {
            offset_angstrom: [0.0, 0.0, 0.0],
            epsilon_kcal_mol: 0.1,
            sigma_angstrom: 2.8,
            charge: 0.70,
        },
        GuestSite {
            offset_angstrom: [1.16, 0.0, 0.0],
            epsilon_kcal_mol: 0.2,
            sigma_angstrom: 3.0,
            charge: -0.35,
        },
    ],
    epsilon_kcal_mol: 0.0,
    sigma_angstrom: 1.0,
    charge: 0.0,
};
```

Each insertion samples a random centre and rotation. Translation moves can also
rotate with probability `rotation_probability`. Host-guest and guest-guest
interactions are summed over every site pair with Lorentz-Berthelot mixing.

## Validation

| Gate | Observable | Result |
|---|---|---|
| Lennard-Jones minimum | `U(2^(1/6)σ) = -ε` | < 1e-9 kcal/mol |
| Ideal gas | Mean loading vs `PV/(kBT)` | within 35 % on a 30 Å box |
| Attractive host | Loading vs empty framework | host ≥ empty |
| Ewald pair | Two opposite charges, 20 Å cell | matches direct Coulomb within 1 % |
| Multi-site | CO2-style three-site model | total charge and site count |
| Invalid input | Negative epsilon, zero cell volume | rejected |

The ideal-gas test is a statistical check of the insertion/deletion factors. It
is intentionally loose because the run is short; the point is to catch a wrong
activity factor, not to replace a production equilibration.

## Bindings

| Surface | Entry point |
|---|---|
| CLI | `sci-form gcmc '<framework-json>' '<guest-json>' '<config-json>'` |
| WASM | `compute_gcmc_isotherm(framework_json, guest_json, config_json)` |
| Python | pending |

The CLI prints the full JSON isotherm. The JSON types are the serialized
`GcmcFramework`, `GuestModel` and `GcmcConfig`.

## Limitations

- The framework is rigid. Framework flexibility and swelling are not modelled.
- The framework is rigid; charges are fixed and no polarization is included.
- Ewald uses a neutralizing background; a non-neutral framework is allowed but
  the physical interpretation of the background must be documented by the user.
- Multi-site guests are rigid; internal degrees of freedom and bonded terms are
  not sampled.
- The isotherm uses a spherical cutoff and no tail correction. The cutoff must be
  large enough for the tested density; the default is 12 Å.
- Ewald parameters are not auto-tuned; the user must check convergence.
- The standard error is a simple block-independent estimate. Production work
  should use block averaging or a correlated-time estimator.
- The framework atom parameters must be supplied explicitly. Automatic UFF
  typing for MOF nodes is a follow-up.
- The implementation has no Widom insertion or volume move.

## Next steps

1. Add Ewald parameter auto-tuning and convergence diagnostics.
2. Add internal degrees of freedom for multi-site guests.
3. Add automatic UFF/MMFF94 typing for framework atoms.
4. Add block averaging and tail corrections.
5. Add Python, WASM and CLI wrappers.
6. Benchmark CO2/N2 isotherms against RASPA and experiment.
