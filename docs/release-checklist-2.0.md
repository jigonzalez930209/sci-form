# v2.0.0 release checklist

> **Released 2026-11-09.** This is the current release gate. The
> [1.0.0 checklist](./release-checklist-1.0.md) is kept as the historical record
> of the first stable release.

`sci-form` 2.0 completes the MMFF94 energy expression, adds analytical ReaxFF
gradients, moves the validated experimental modules into the core crate, and
replaces the `String` error type of the root Rust API with `SciFormError`.
Version 2.0.0 is published on crates.io, PyPI and npm, and the GitHub Release
carries the CLI packages, wheels, sdist and `SHA256SUMS`.

Legend: `[x]` done, `[ ]` pending.

---

## 1. API and migration

- [x] The 65 root-level functions in `src/lib.rs` return
      `Result<_, SciFormError>` instead of `Result<_, String>`.
- [x] `SciFormError` is `#[non_exhaustive]` with six variants: `Parse`,
      `Validation`, `Geometry`, `Numerical`, `Unsupported` and `Other`.
- [x] `SciFormError` implements `Display`, `std::error::Error`,
      `Deref<Target = str>`, `AsRef<str>`, `From<String>` and `From<&str>`, so
      existing `?`-based code and `error.to_string()` call sites keep working.
- [x] The compatibility aliases for the promoted modules still compile; the
      old `alpha::*` / `beta::*` paths resolve to the same implementation.
- [x] MBH, DFT, MLFF and reaction dynamics remain gated behind their research
      feature flags.
- [x] `docs/guide/migration-2.0.md` documents every rename, the new error type
      and the four remaining research flags.
- [x] The Python, WASM and CLI packages report version `2.0.0`; old Python and
      TypeScript names remain available as deprecation aliases.

## 2. Science and coverage

- [x] **MMFF94 is complete.** Stretch-bend, out-of-plane and buffered Coulomb
      electrostatics with PBCI/BCI charges are implemented; the seven-molecule
      RDKit parity check reaches a maximum energy deviation of
      `2.8e-11 kcal/mol`, and the partial charges match exactly.
- [x] **ReaxFF analytical gradients.** Bond-order scaling, bond/angle,
      over/under-coordination, van der Waals, Coulomb and EEM charge response
      (`dq = -A^-1 dA q`) are differentiated and gated against central
      differences on H2, water, methane and CH3O.
- [x] **ReaxFF ffield parser** covers the general, atoms, bonds, off-diagonal,
      angles, torsions and hydrogen-bond sections.
- [x] **RandNLA N>500.** `solve_eht_randnla_lowrank` combines a Nyström sketch
      with re-orthogonalized Lanczos iterations; at N=600, k=160 the four
      lowest states reach `9.6e-4` relative error in about 80 ms versus 143 ms
      for the exact path.
- [x] **Cholesky ERI.** `hf::integrals::CholeskyEris` reconstructs the packed
      tensor from a thresholded pivoted Cholesky factor; benzene STO-3G has
      rank 237/666, a 0.714 ratio and a maximum error of `9.9e-7 Ha`.
- [x] **PM3 and xTB gradient cost.** PM3 uses the stationary-density
      (Pulay-free) form, so a gradient no longer runs one SCF per displacement;
      xTB warm-starts every displaced SCC from the reference charges and
      multipoles.
- [x] **NMR HOSE.** `HoseDatabase` reads JSON/CSV/TSV exports from NMRShiftDB2
      or SDBS and overrides 1H/13C shifts on the longest-prefix match. No
      database data is bundled.
- [x] **Promotion sweep.** GSM, kinetics/HTST, periodic-linear adapters, live
      MD, RandNLA, KPM, CGA, SDR, EDL and the render-bridge domain builders are
      core modules with compatibility aliases for the historical feature flags.

## 3. Validation gate

- [x] The external campaign was re-run against pinned references:
      RDKit 2026.09.1, PySCF 2.14.0, tblite 0.7.0, dftd4 4.3.0, xtb 6.7.1,
      multicharge 0.6.0, MOPAC 23.2.5, LAMMPS 29Aug2024_update2, ASE 3.29.0,
      pymatgen 2026.9.24, pormake 0.2.3 and FreeSASA.
- [x] `reports/validation/summary.json` records **129 rows: 122 PASS, 7 method
      or provenance gaps and 0 FAIL**.
- [x] **36/36** oracle scripts and **117/117** Rust validation tests pass.
- [x] Python 3.14.4 is the oracle runtime.
- [x] The PySCF AO→MO fixture canonicalises degenerate MO subspaces before
      comparison.
- [x] `tests/gpu_parity.rs` adds eleven GPU-vs-CPU gates (ERI, Fock, gamma,
      orbital and density grids, ESP, EEQ, CPM, HTST, Bloch and EDL). The suite
      skips cleanly without an adapter; `SCIFORM_REQUIRE_GPU=1` turns that skip
      into a failure.
- [x] The ERI GPU kernel rejects p/d/f bases instead of silently returning an
      `(ss|ss)`-only result.
- [x] Every row, tolerance and owning test is published in the
      [full validation matrix](/validation/full-matrix); the narrative report is
      [Scientific validation](/validation/scientific-validation).

## 4. Documentation

- [x] The VitePress site builds with VitePress 1.6.4 and deploys to GitHub
      Pages.
- [x] The documentation audit covers **71 pages and 830 MathJax expressions**
      at 1280, 1024, 768 and 414 px: **0 MathJax errors, 0 raw LaTeX
      delimiters and 0 page overflows**.
- [x] All prose in `docs/` is English. The historical Spanish validation
      notebook lives outside `docs/` for provenance.
- [x] `algorithm/ks-dft.md` and `algorithm/fluorescence.md` use the MathJax
      delimiters the renderer supports (`$...$` / `$$...$$`) and render 13 and
      15 expressions respectively.
- [x] **56 animated SVG diagrams** are transparent, theme-aware, in motion and
      bounded by their `viewBox`; `scripts/docs/check_animated_svgs.py` enforces
      all four properties.
- [x] Mermaid flowcharts render through the local `Mermaid.vue` component with
      light and dark themes; no external CDN is required at runtime.
- [x] Wide tables are bounded to `75vh` with sticky headers and internal
      vertical/horizontal scrolling; display equations scale to the content
      width on small screens.
- [x] `docs/guide/migration-2.0.md` and the v2 sections of `CHANGELOG.md` are
      current.

## 5. Dependencies, packaging and CI

- [x] nalgebra 0.35 / simba 0.10 removed the unmaintained `paste` crate;
      RUSTSEC-2024-0436 is closed and `deny.toml` has no advisory ignores.
- [x] Workspace MSRV is `1.89`.
- [x] Native `src/gpu/` and the WASM runtime share wgpu 29; the duplicated
      naga/ash/wgpu-core trees are gone.
- [x] Criterion infrastructure lives in `benches/core.rs` with the `core`
      benchmark target.
- [x] Runnable doctests rose from 15 to 28 (plus one ignored).
- [x] `scripts/release-check.sh` is green and stamps `.release-check.ok` with
      the verified commit.
- [x] `scripts/release.sh 2.0.0` refuses to tag unless the stamp matches the
      current commit and the working tree is clean.
- [x] The GitHub Release `v2.0.0` contains 14 assets: Linux `.tar.gz` for
      x86_64/aarch64, macOS `.tar.gz` for x86_64/arm64, the Windows `.zip`,
      `.deb` and `.rpm` packages for both Linux architectures, three Python
      wheels, the sdist and `SHA256SUMS`.
- [x] crates.io (`sci-form`, `sci-form-cli`), PyPI (`sciforma`) and npm
      (`sci-form-wasm`) serve `2.0.0`.
- [x] CI run `42802896874` and documentation deployment `42803485203` completed
      successfully for the released `main`.

## 6. Known limitations carried into 2.0

These are documented, not hidden. The full list is in
[Scientific validation](/validation/scientific-validation) and
`reports/validation/summary.json`.

- Seven validation rows are method or provenance gaps rather than passing
  gates; none is a failing parity check.
- No NMRShiftDB2 or SDBS data is bundled with `HoseDatabase`; the CC-BY-SA
  review remains external.
- Symbolic NDDO and GFN2 derivatives remain research work. PM3 no longer runs
  an SCF per displacement, and xTB warm-starts each displaced SCC, but the
  gradients are still numerical.
- MBH, DFT, MLFF and reaction dynamics remain gated research tracks.
- Heavy `#[ignore]` benchmarks and the GFN2 reaction-dynamics sweep are run
  locally with `scripts/release-check.sh --heavy`, not on every CI push.

---

## Publish procedure

```bash
# 1. Commit everything, then run the whole gate on this machine.
git status --short                      # clean
scripts/release-check.sh                # add --heavy for the #[ignore]d tests

# 2. Tag and push; the script refuses unless the gate passed for this commit.
scripts/release.sh 2.0.0

# 3. Watch the release workflow: it compares the tag with every manifest,
#    builds the CLI packages, wheels and npm package, publishes crates.io,
#    PyPI and npm, creates the GitHub Release with SHA256SUMS, and deploys
#    the documentation on the freshly compiled WASM package.
```

## Post-release

- [x] The GitHub Release `v2.0.0` is published with the release body and all
      14 assets.
- [x] `cargo install sci-form-cli --version 2.0.0`, `pip install sciforma` and
      `npm install sci-form-wasm` resolve the 2.0.0 artefacts.
- [x] The documentation site deploys from the rewritten `main` and renders the
      equations, tables and diagrams described above.
