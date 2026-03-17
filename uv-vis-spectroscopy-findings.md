# sci-form Spectroscopy Search - Complete Findings

## UV-Vis Spectroscopy - FOUND ✅

Complete UV-Vis implementation exists. All code is in **exploratory/qualitative mode** (not calibrated TDDFT).

### Core Functions:

**src/reactivity.rs** (lines 366-428):
- `compute_uv_vis_like_spectrum()` — Main function builds spectrum from EHT occupied→virtual transitions
- Uses **Gaussian broadening** (not Lorentzian)
- Computes intensity from coefficient overlap approximation 
- Up to 24 peaks tracked, sorted by intensity
- Outputs energy grid + intensities + peak list

**src/lib.rs** (lines 808-823):
- Public wrapper `compute_uv_vis_spectrum()` 
- Calls EHT solver, then `compute_uv_vis_like_spectrum()`

### Data Structures (src/reactivity.rs lines 129-152):
- `UvVisPeak` — energy_ev, wavelength_nm, intensity, from_mo, to_mo
- `UvVisSpectrum` — energies_ev, intensities, peaks, sigma, notes

### Bindings:
- Python: `uvvis_spectrum(elements, coords, sigma=0.25, e_min=0.5, e_max=8.0, n_points=600)` (crates/python/src/lib.rs:1260-1280)
- WASM: `compute_uv_vis_spectrum()` returns JSON (crates/wasm/src/lib.rs:568-591)

---

## Hessian / Vibrational / IR - NOT FOUND ❌
- Mentions in BFGS optimization code only (src/optimization/bfgs.rs, src/forcefield/bounds_ff/bfgs.rs)
- These are for geometry optimization, not vibrational analysis

## NMR - PLANNED NOT IMPLEMENTED ❌
- See docs/nmr_detailed_roadmap.md for full 4-phase strategy:
  - Phase 1: HOSE code generation
  - Phase 2: Static database via phf crate
  - Phase 3: 3D J-couplings via Karplus equation
  - Phase 4: Lorentzian broadening & API

## Lorentzian Broadening - PLANNED ONLY ❌
- Mentioned in docs/spectroscopy_roadmap.md line 45
- Mentioned in docs/nmr_detailed_roadmap.md lines 56-58 with formula
- Current code uses Gaussian only (see fn gaussian() in reactivity.rs)

## Summary Table
| Feature | Status | Location |
|---------|--------|----------|
| UV-Vis spectrum | ✅ Implemented | src/reactivity.rs, src/lib.rs |
| Gaussian broadening | ✅ Implemented | src/reactivity.rs:fn gaussian() |
| Hessian/Vibration | ❌ Not found | — |
| IR spectrum | ❌ Not found | — |
| NMR shifts | ❌ Planned | docs/nmr_detailed_roadmap.md |
| Lorentzian broadening | ❌ Planned | docs/*.md roadmaps |
