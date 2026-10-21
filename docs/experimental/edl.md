# EDL: Electrical Double Layer (core)

**Module:** `sci_form::edl` (core since 2026-10-14; `alpha::edl` re-export)
**Feature flag:** none (core)

---

## Overview

Continuum models of the electrochemical double layer at charged interfaces:
Helmholtz compact layer, Gouy-Chapman diffuse layer, and the combined
Gouy-Chapman-Stern (GCS) model, plus coupling adapters to CPM (`cpm`),
EEQ charges (`charges_eeq`), and ALPB solvation (`solvation_alpb`).

---

## Theory

### Gouy-Chapman (1:1 electrolyte)

$$\phi(x) = \frac{2kT}{ze}\ln\frac{1+\gamma e^{-\kappa x}}{1-\gamma e^{-\kappa x}}, \quad \gamma = \tanh\frac{ze\phi_0}{4kT}$$

with Debye parameter $\kappa^{-1} = \sqrt{\varepsilon_0\varepsilon_r kT / 2c_0 z^2 e^2 N_A}$.

### Grahame equation

$$\sigma = \sqrt{8c_0\varepsilon_0\varepsilon_r kT}\,\sinh\frac{ze\phi_d}{2kT}$$

inverted for the Stern split $\phi_d = (2kT/ze)\,\mathrm{asinh}(\sigma/\sqrt{\cdots})$.

### Stern series capacitance

$$C^{-1} = C_H^{-1} + C_\mathit{diff}^{-1}, \quad C_H = \varepsilon_0\varepsilon_c/d, \quad C_\mathit{diff} = \varepsilon_0\varepsilon_r\kappa\cosh\frac{ze\phi_d}{2kT}$$

---

## API

- `compute_helmholtz_profile(σ, config)` — compact capacitor.
- `compute_gouy_chapman_profile(φ₀, config)` — diffuse layer.
- `compute_gcs_profile(σ, config)` — Stern split + series capacitance.
- `compute_edl_profile(x, config)` — dispatches on `EdlModel`.
- `scan_edl_capacitance(σmin, σmax, n, config)` — C(σ) curve (even, minimum at σ=0).
- Couples: `compute_alpb_stern_correction`, `compute_eeq_edl_profile`,
  `compute_cpm_helmholtz_scan`, `couple_cpm_eeq_edl`.

---

## Validation (F86, 2026-10-14)

External gate `tests/alpha_edl_validation.rs` (no feature needed) vs
`scripts/validation/oracle_edl.py` (independent numpy implementation of the
textbook closed forms in SI): 3 temperatures × 3 ionic strengths × 6
potentials/charges — profiles (φ, E, ρ), Stern split, all three
capacitances, worst 2.0e-14 (tolerance 1e-9). Fixture
`tests/fixtures/validation/edl_analytic.json`.

---

## Tests

```bash
cargo test --test alpha_edl_validation
```

4 tests in `tests/alpha_edl_validation.rs` (gate F86) comparing the Gouy-Chapman, Grahame-inverted Stern split, Helmholtz compact capacitor and Stern series capacitance against the analytic closed forms in the fixture `tests/fixtures/validation/edl_analytic.json` (oracle `scripts/validation/oracle_edl.py`).
