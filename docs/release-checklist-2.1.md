# Release checklist `2.1.0`

**Status: released 2026-11-25.** V2-11 (GCMC) and V2-12 (reaction
networks/kMC) are promoted to core in `v2.1.0`. V2-05 (Raman/VCD), V2-08
(pKa/redox) and V2-10 (exact OPES flow) remain alpha for `2.2`.

- GitHub Release: <https://github.com/jigonzalez930209/sci-form/releases/tag/v2.1.0>
- crates.io: `sci-form 2.1.0`
- PyPI: `sciforma 2.1.0`
- npm: `sci-form-wasm 2.1.0`

## Frozen scope

| Candidate | Decision | Evidence |
|---|---|---|
| V2-11 GCMC | **Promoted** | F116 2.0 %, F120 0.7 %, F131 3.9 %, F132 2.5 % |
| V2-12 networks/kMC | **Promoted** | F126 5/5, F127 0.11 %, F128 deterministic, F133 0.575 kcal/mol |
| V2-10 enhanced sampling | Alpha in 2.1 | F124/F125/F129 pass; exact OPES flow and >2 CV pending |
| V2-08 pKa/redox | Alpha in 2.1 | MAE 2.33, above the 1.0 gate |
| V2-05 Raman/VCD | Alpha in 2.1 | Raman position MAE 82 cm⁻¹ (F130); no magnetic-dipole derivatives |

## Promotion checks

- [x] Gate L3/L4 evidence recorded in `reports/validation/summary.json`.
- [x] Modules exported without feature gate (`gcmc`, `networks`, `kmc`).
- [x] `alpha-gcmc` and `alpha-reaction-network` kept as empty aliases.
- [x] CLI commands available in a default build (`gcmc`, `network`, `kmc`).
- [x] Python and WASM bindings available without feature flags.
- [x] API pages: `docs/algorithm/gcmc.md`, `kmc.md`, `reaction-networks.md`.
- [x] Diagrams and Rust/Python/CLI examples.
- [x] Criterion benchmarks: `benches/gcmc.rs`, `benches/kmc.rs`.
- [x] Version `2.1.0` in the workspace, CLI, Python, WASM and `package.json`.
- [x] Changelog section `v2.1.0`.
- [x] `cargo clippy --workspace --all-targets` with `-D warnings`.
- [x] Core Rust suites (alpha-v2 lib, v2_validation, v21_validation, CLI, gradient parity) and doctests.
- [x] Python wheel + pytest (196 passed), WASM build + Node tests (20 passed, 155 live examples).
- [x] Docs build (`npm run docs:build`).
- [x] Tag `v2.1.0` and publish to crates.io/PyPI/npm/GitHub (workflow success).

## Release commands

```bash
cargo fmt --all -- --check
scripts/cargo-guard.sh test --profile numeric -p sci-form --lib
scripts/cargo-guard.sh test --profile numeric -p sci-form --features alpha-v2 --test v21_validation
scripts/cargo-guard.sh bench --profile numeric --bench gcmc --features alpha-v2 -- --test
scripts/cargo-guard.sh bench --profile numeric --bench kmc --features alpha-v2 -- --test
python3 scripts/check_feature_parity.py
```

## Known limitations in `2.1.0`

- GCMC frameworks are rigid; no framework flexibility or polarization.
- Ewald parameters are user-supplied; convergence must be checked.
- Curated bimolecular networks use literature barriers as rule inputs. F133
  validates the barrier engine on torsional/inversion barriers, not on the
  bimolecular reactions.
- Surface kMC (adsorption/desorption/diffusion) is not included.
- V2-10 exact OPES flow, >2 collective variables and V2-05 VCD remain alpha.
- The H2 frozen-density KS gradient unit test is ignored: it deviates 2.3 %
  from the grid energy derivative, while the F80 H2O2 gradient gate passes at
  4e-7 Ha/Bohr. This is a known limitation of the atom-centred grid response
  for compact diatomics and is outside the promoted 2.1 scope.
