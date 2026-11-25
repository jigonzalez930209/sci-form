# Release checklist `2.2.0`

**Status: promotion in progress.** The frozen scope is V2-10 (enhanced
sampling) and V2-14 (active learning). V2-05, V2-08, V2-09, OPES with more
than two CV and the remaining band-B candidates stay alpha for `2.3`.

## Frozen scope

| Candidate | Decision | Evidence |
|---|---|---|
| V2-10 enhanced sampling | **Promoted** | F124 PMF, F125 2D OPES, F129 acceptance, F134 flow residual 4.7e-11 |
| V2-14 active learning | **Promoted** | F114 conformal coverage, F135 acquisition gap ratio 0.398 |

## Promotion checks

- [x] Modules exported without feature gate (`enhanced_sampling`, `enhanced_sampling_multi`, `active_learning`).
- [x] `alpha-enhanced-sampling` and `alpha-active-learning` kept as empty aliases.
- [x] CLI commands available in a default build (`metadynamics`, `opes`, `replica-exchange`, `active-learning`, `conformal`).
- [x] Python and WASM bindings available without feature flags.
- [x] API pages: `docs/algorithm/enhanced-sampling.md`, `active-learning.md`.
- [x] Rust/Python/CLI examples and catalog entries.
- [x] Criterion benchmarks: `benches/enhanced_sampling.rs`, `benches/active_learning.rs`.
- [x] Version `2.2.0` in workspace, CLI, Python, WASM, `pyproject.toml` and `package.json`.
- [x] Changelog section `v2.2.0`.
- [x] `cargo clippy --workspace --all-targets --features alpha-v2 -- -D warnings`.
- [x] Core block: lib 1028 pass; v2_validation 15 pass; v21_validation 17 pass + 1 ignored; CLI 4 pass.
- [x] Python wheel 2.2.0 + 196 pytest pass; WASM build 2.2.0 + 20 Node tests + 155 live examples.
- [x] Docs build (`npm run docs:build`) and WASM docs package.
- [ ] Tag `v2.2.0` and publish.

## Release commands

```bash
cargo fmt --all -- --check
scripts/cargo-guard.sh clippy --profile numeric --workspace --all-targets --features alpha-v2 -- -D warnings
scripts/cargo-guard.sh test --profile numeric -p sci-form --lib
scripts/cargo-guard.sh test --profile numeric -p sci-form --features alpha-v2 --test v2_validation
scripts/cargo-guard.sh test --profile numeric -p sci-form --features alpha-v2 --test v21_validation
scripts/cargo-guard.sh test --profile numeric -p sci-form-cli
python3 scripts/check_feature_parity.py
npm run docs:build
```

## Known limitations in `2.2.0`

- OPES supports one and two collective variables; the >2 CV extension is 2.3.
- FEP/BAR/MBAR remain estimators without a soft-core MD driver.
- V2-05 Raman/VCD and V2-08 pKa/redox remain alpha.
- The H2 frozen-density KS gradient test stays ignored; F80 validates H2O2.
