# v1.0.0 release checklist

Everything that must be true before `sci-form` 1.0.0 is published to
crates.io, PyPI and npm.

Legend: `[x]` done, `[ ]` pending.

---

## 1. The code

- [x] `cargo fmt --all -- --check` is clean.
- [x] `cargo build --workspace --all-targets --features "<all>"` succeeds with
      every feature flag on.
- [x] `cargo clippy --workspace --all-targets --features "<all>" -- -D warnings`
      is clean (0 warnings).
- [x] `cargo test --lib --release` → 881 passed, 0 failed, 2 ignored.
- [x] `cargo test --doc --release` → 15 passed, 0 failed, 1 ignored.
- [x] All 9 integration test targets pass (78 tests with every feature on;
      60 of them are compiled without the alpha/beta features).
- [x] `cargo doc --no-deps --all-features` → 0 warnings.
- [x] `cargo clippy -- -W missing_docs` → 0 undocumented public items
      (was 1038 of 2121).
- [x] No panicking site in the library code is reachable from user input.
- [x] The force-field energy and gradient are bit-identical regardless of the
      rayon thread count (verified in a subprocess with 1, 4 and 16 threads).
- [x] Every public enum is `#[non_exhaustive]` (57 of 60; 3 deliberately
      closed, documented in the migration guide).
- [x] Every returned energy, length, frequency and charge states its unit.

## 2. The bindings

- [x] Python: `cargo build -p sci-form-python` succeeds; the wheel builds with
      `maturin build --release`.
- [x] Python: the pytest suite passes (`tests/integration/test_python_integration.py`).
- [x] Python: `sci_form.__all__` matches the registered functions.
- [x] Python: exactly one `cpm_charges` is registered (behind the
      `experimental-cpm` / `beta-cpm` feature flag) and it is the μ =
      −4.44 eV one; the β copy is `beta_compute_cpm_charges`.
- [x] WASM: `crates/wasm/build.sh` produces a working `pkg/` and `pkg-node/`
      with the `alpha`/`beta` subpaths and the `snippets/` directory.
- [x] WASM: the published package is dual. `pkg/cjs/` holds the CommonJS
      half, the root manifest points its `require` condition at it, and
      `npm pack` + `npm install` + `require('sci-form-wasm')` is what the CI
      and release gates exercise — not the raw build directory.
- [x] WASM: every error path emits valid JSON.
- [x] CLI: `sci-form --help` and every subcommand's `--help` render.
- [x] CLI: an optional argument is absent, not a sentinel value.

## 3. The validation gates

- [x] The 127 validation fixtures are versioned in the repository, and
      `reports/validation/summary.json` records 120 passing and 7 pending
      (F44, F46, F71, F82, F83, F84, F85).
- [x] The CI validation matrix runs 9 hermetic gate jobs (no RDKit, PySCF,
      ASE, tblite or MOPAC needed at test time).
- [x] `reports/validation/summary.json` is regenerated and committed.
- [x] The 12 open findings of the first audit are documented in
      `reports/auditoria-resolucion.md` with the reason each one is still open.
- [x] The 88 findings of the second and third audits are dispositioned in
      `reports/v1.0-blockers.md`: 7 are accepted with a written justification
      and the rest are closed, each with a regression test or a documented
      change. The five closed in the last round (B-H7, C-H2, C-H3, C-H4 and
      A-C13) carry a `Estado:` line in their own entry.

## 4. The documentation

- [x] The crate-level rustdoc states what the library is AND IS NOT, the
      coordinate/element/sign conventions, the unit table, the capability
      lookup table and the feature-flag table.
- [x] Every public item has a doc comment; `# Errors` and `# Panics` sections
      where they apply.
- [x] `docs/api/{rust,python,typescript,cli}.md` are regenerated against the
      current code.
- [x] `docs/api/complete-supported-reference.md` covers every capability on
      every surface, with 0 orphans in either direction.
- [x] `docs/guide/migration-1.0.md` lists every breaking change.
- [x] `CHANGELOG.md` has a `v1.0.0` section.
- [x] The VitePress site builds.

## 5. Packaging

- [x] The version is `1.0.0` in all four manifests, `pyproject.toml`,
      `package.json` and `pkg-node/package.json`.
- [x] The root `Cargo.toml` uses an allow-list `include`; the published crate
      is 340 files / 8.0 MiB.
- [x] The 27 developer binaries are not in the published crate.
- [x] `rust-version = "1.87"` is declared in all four manifests.
- [x] `license`, `repository`, `homepage`, `keywords`, `categories` are set in
      every manifest.
- [x] `cargo package` succeeds for the root crate. The three binding crates
      depend on `sci-form = "1.0.0"`, so their `cargo package` only succeeds
      **after** the root crate is on crates.io — run it again once the publish
      step is done.
- [x] `cargo deny check advisories bans licenses` passes.
- [x] `cargo audit` reports no unfixed advisory other than RUSTSEC-2024-0436
      (`paste` via `simba` → `nalgebra`), which has no safe fix.

## 6. CI and release

- [x] `ci.yml` runs fmt, clippy (`--workspace`), the unit tests, the 9-target
      integration suite, the validation matrix, the Python wheel, the WASM
      build and the audit.
- [x] `docs.yml` rebuilds the site on a change in `src/**` or `crates/**` as
      well as in `docs/**`.
- [x] `release.yml` compares the git tag with the version of all four
      packages before publishing anything.
- [x] `release.yml` publishes the crates.io crate, the PyPI wheel and the npm
      package.

---

## Publish procedure

```bash
# 1. Make sure the tree is green and the version is right.
cargo test --workspace --all-targets --release
grep -n '^version' Cargo.toml crates/*/Cargo.toml     # all 1.0.0
git status --short                                     # clean

# 2. Tag. The release workflow refuses to publish a tag that does not match
#    the manifest versions, so a mismatch fails fast and safely.
git tag -a v1.0.0 -m "sci-form 1.0.0"
git push origin main --tags

# 3. Watch the workflow.
#    It publishes, in order: crates.io (sci-form, sci-form-cli), the PyPI
#    wheel (sciforma) and the npm package (sci-form-wasm).
```

## Post-release

- [ ] Tag the release notes from `CHANGELOG.md`.
- [ ] Verify `cargo install sci-form-cli` and `pip install sciforma` and
      `npm install sci-form-wasm` from a clean environment.
- [ ] Verify `https://docs.rs/sci-form/1.0.0` builds without warnings.
- [ ] Open the 1.1 milestone for the typed error enum (the error type is
      still `String`).

## Known limitations carried into 1.0

These are documented, not hidden. See `reports/auditoria-resolucion.md` and
`reports/v1.0-blockers.md` for the full list.

- The error type is `String`, not a typed enum (planned for 1.1).
- The ERI tensor is a full M⁴ `Vec<f64>` — 12.8 GB at M = 200. No packed,
  Cholesky or density-fitting path.
- The conformer minimiser rebuilds the force field on every MD step.
- The ETKDG retry loop is 10 × N attempts of an O(N²) search.
- The PM3, GFN0-xTB and KS-DFT analytic gradients are quantitatively wrong
  (18.4 eV/Å, 209–263 eV/Å and 0.202 Ha/Bohr). The affected tests are
  `#[ignore]`d with the measured numbers in them.
- The 3D reaction path (`alpha-reaction-dynamics`) is **not a validated
  reaction path**: the CI-NEB frequently does not converge and the image
  energies are computed with the reactant topology. The frames are always
  finite, and the module says so in its notes and in
  `docs/experimental/alpha.md` (A14).
- NPA/NBO has no parity reference (it needs Psi4/NBO 6.0).
