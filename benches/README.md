# Benchmarks

Criterion micro-benchmarks for the hot paths of the library:

- `parse` — SMILES → molecule graph.
- `embed` — full ETKDG-style conformer generation (ethanol, aspirin, cholesterol).
- `force_field` — UFF and MMFF94 energy evaluation on the same aspirin geometry.
- `quantum` — EHT, PM3 and GFN2-xTB on water, methane and benzene.

Run them with:

```bash
cargo bench --bench core
cargo bench --bench core -- embed
```

Criterion writes the HTML report to `target/criterion/report/index.html`. The
inputs are deterministic (`seed = 42`) so two runs on the same machine are
comparable. These benchmarks are **not** part of `scripts/release-check.sh`:
the release gate runs correctness tests, not timings.
