# Changelog

All notable changes to this repository are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows semantic versioning for published packages.

## [0.10.6] - 2026-03-26

### Added
- Expanded the supported scientific surface across Rust, Python, WASM, and CLI for
  PM3/PM3(tm), GFN0/GFN1/GFN2-xTB, HF-3c, ANI, spectroscopy, solvation,
  stereochemistry, rings/fingerprints, clustering, periodic systems, and materials.
- Added a comprehensive VitePress documentation site covering algorithms, guides, API
  references, experimental modules, benchmarks, and WebGPU validation workflows.
- Added browser-side WebGPU validation and benchmark flows, including the documented
  `webgpu:benchmark:*` scripts for the WASM hybrid example.

### Changed
- Updated published package versions across the workspace to `0.10.6`.
- Refreshed the README and docs landing pages so the project description matches the
  current public surface instead of the earlier ETKDG/EHT-only summary.
- Expanded local validation and release documentation for Rust, Python, WASM, CLI,
  docs publishing, and cross-ecosystem releases.

### Infrastructure
- Added GitHub Actions coverage for CI, docs publication, and automated multi-package
  release workflows.

