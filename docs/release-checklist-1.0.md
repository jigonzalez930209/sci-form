# v1.0.0 release checklist

> **Historical document.** This checklist records the 1.0.0 release gate. The
> 2.0 follow-up (2026-11-08) closed MMFF94, ReaxFF analytical gradients, typed
> root-API errors, ERI Cholesky and the RandNLA N>500 path; see
> `reports/1.1-experimental-2026-11-07.es.md` and `CHANGELOG.md`. Items below
> that mention `String` errors, M⁴-only ERI storage or central-difference
> ReaxFF gradients describe the 1.0.0 state, not the current one.

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
- [x] `cargo test --lib --release` → 890 passed, 0 failed, 0 ignored
      (1022 passed with the CI feature set).
- [x] `cargo test --doc --release` → 17 passed, 0 failed, 1 ignored.
- [x] `cargo test --workspace --exclude sci-form-python --all-targets --release
      --features "<CI feature set>"` → 1962 passed, 0 failed, 78 ignored
      (52 are `--ignored` benchmarks, 17 are the `reaction_dynamics_study` sweep (GFN2, >10 min); the rest are heavy SCF runs and two that
      need an RDKit-generated file). The `sci-form-python` crate is excluded
      because a PyO3 extension module cannot link in `cargo test`; its tests
      are the pytest suites below.
- [x] `tests/gradient_parity.rs` (F80): 10 tests, none ignored. PM3 1.1e-5 eV/Å,
      xTB 4.8e-5 eV/Å, KS-DFT 3.1e-7 Ha/Bohr, UKS (OH) 1.0e-6 Ha/Bohr, UFF/MMFF94/EHT/ANI ≤ 1e-6.
- [x] `cargo doc --no-deps --all-features` → 0 warnings, and
      `cargo doc --workspace` works too (the Python library and the CLI binary
      are marked `doc = false` so they do not overwrite the root crate's
      document).
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
- [x] Python: the pytest suites pass — `tests/integration/test_python_integration.py`
      (99) and `crates/python/tests/test_bindings.py` (92); both run in CI.
- [x] Python / CLI / WASM cross-language: `tests/integration/test_cross_language.py`
      70 of 70 (run manually: it needs the release CLI, a wheel and the WASM
      `pkg-node`).
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
- [x] WASM: `tests/integration/test_wasm_integration.js` 19 of 19 (runs in CI).
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

- [x] **No tests run on GitHub** (limited minutes). `ci.yml` only *builds* — the
      library with every feature, the CLI, the Python wheel and the WASM
      packages — with `Swatinem/rust-cache`, and skips documentation-only
      changes.
- [x] **All tests run locally**: `scripts/release-check.sh` (versions, fmt,
      clippy, the Rust suite, doctests, the CLI guide, the CLI packages, pytest,
      WASM, the 155 live examples, cross-language, the docs build). It writes
      `.release-check.ok` with the commit it verified.
- [x] `scripts/release.sh <version>` creates and pushes the tag **only** if the
      stamp matches the current commit and the tree is clean.
- [x] `release.yml` (on the tag) compares the tag with every manifest, then
      builds: the CLI for Linux x86_64/aarch64, macOS x86_64/arm64 and Windows
      (`.tar.gz`/`.zip`, plus `.deb` and `.rpm` on Linux), the wheels + sdist,
      and the npm package; publishes crates.io, PyPI and npm; creates the GitHub
      Release with the packages and `SHA256SUMS`; and finally deploys the
      documentation.
- [x] The documentation is deployed **after** the builds (`docs.yml`, called by
      `release.yml`), on the WASM package that was just compiled, so its live
      examples run the released code. `docs.yml` can also be run by hand.
- [x] Every public WASM function has a live example (`docs/examples/catalog.json`,
      155 examples) and `scripts/test-live-examples.mjs` runs them all against
      the compiled package.

---

## Publish procedure

```bash
# 1. Commit everything, then run the whole gate on this machine.
git status --short                      # clean
scripts/release-check.sh                # add --heavy for the #[ignore]d tests

# 2. Tag and push — the script refuses unless the gate passed for this commit.
scripts/release.sh 1.0.0

# 3. Watch the workflow (builds only; no tests run on GitHub).
#    It publishes crates.io (sci-form, sci-form-cli), PyPI (sciforma), npm
#    (sci-form-wasm), creates the GitHub Release with the CLI packages
#    (.tar.gz, .zip, .deb, .rpm) and deploys the documentation.
```

## Post-release

These four items could only be verified once the tag was pushed and the
registries had the release. Verified on 2026-11-07:

- [x] Tag the release notes from `CHANGELOG.md`. The GitHub Release `v1.0.0`
      carries a release body, and the changelog entry is now dated
      `2026-11-03` with the final-build fixes folded into it.
- [x] Verify `cargo install sci-form-cli` and `pip install sciforma` and
      `npm install sci-form-wasm` from a clean environment. Verified against
      the published 1.0.0 artefacts: the CLI installs with
      `cargo install --locked sci-form-cli --version 1.0.0`, the Python wheel
      reports `1.0.0`, and the npm package embeds ethanol with 9 atoms.
- [x] Verify `https://docs.rs/sci-form/1.0.0` builds without warnings. Build
      #5141712 succeeded; the only log warning is docs.rs' own
      "target filter specified, but no targets matched" no-op, not a rustdoc
      warning.
- [x] Open the 1.1 milestone for the typed error enum (the error type is
      still `String`). Opened as milestone #1:
      <https://github.com/jigonzalez930209/sci-form/milestone/1>.

## Known limitations carried into 1.0

These are documented, not hidden. See `reports/auditoria-resolucion.md` and
`reports/v1.0-blockers.md` for the full list.

- The error type is `String`, not a typed enum (planned for 1.1).
- The ERI tensor is a full M⁴ `Vec<f64>` — 12.8 GB at M = 200. No packed,
  Cholesky or density-fitting path.
- The conformer minimiser rebuilds the force field on every MD step.
- The ETKDG retry loop is 10 × N attempts of an O(N²) search.
- The PM3 and xTB (GFN2) gradients are central differences of the reported
  energy (`6·N` SCF runs per gradient), not analytic expressions. They agree
  with the energy to ~1e-4 eV/Å (gates F80). An analytic NDDO/GFN2 gradient is
  future work. The KS-DFT gradient is analytic and was never wrong; its F80
  gate had a unit error (Ha/Å vs Ha/Bohr) and is now enforced.
- `"xtb"` is GFN2 on every surface (NEB, MD, reaction dynamics, GSM,
  `compute_xtb`); `"gfn0"` is only a legacy alias of it. The faithful GFN0 is
  `xtb::solve_gfn0` (Rust only).
- The 3D reaction path (`alpha-reaction-dynamics`) exists **only in the Rust
  core**; the Python, WASM and CLI bindings do not expose it. It is **not a validated
  reaction path**: the CI-NEB frequently does not converge and the image
  energies are computed with the reactant topology. The frames are always
  finite, and the module says so in its notes and in
  `docs/experimental/alpha.md` (A14).
- NPA/NBO has no parity reference (it needs Psi4/NBO 6.0).
