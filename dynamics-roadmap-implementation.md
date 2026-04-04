> **Historical document** — This session log describes the initial dynamics module implementation. Alpha module architecture is now documented in `docs/experimental/alpha.md`.

k# Dynamics Roadmap Implementation Session

## Key patterns observed:
- Alpha modules go under `src/alpha/` with feature flags `alpha-*` 
- But for these new modules, we're creating them as top-level modules with `alpha-` feature gates
- Actually, looking at the roadmap, these are NEW modules at `src/dynamics/`, `src/dft/`, `src/forcefield/reaxff/`, `src/scf/obara_saika.rs` etc.
- The user wants ALL new modules gated as ALPHA until development/testing complete
- We CAN use existing stable core code

## Feature flag pattern:
```toml
alpha-dynamics-live = []
alpha-dft = []
alpha-reaxff = []
alpha-obara-saika = []
alpha-mlff = []
alpha-imd = []
```

## Implementation plan:
1. Task 1: DynamicAtom/MolecularSystem + LiveSimulation WASM (src/dynamics/ refactor)
2. Task 2: Obara-Saika integrals + SCF hardening (src/scf/)
3. Task 3: DFT (src/dft/ new module)
4. Task 4: PBC + ensembles (src/dynamics/)
5. Task 5: ReaxFF (src/forcefield/reaxff/)
6. Task 6: MLFF inference (src/ml/)
7. Task 7: IMD browser loop (crates/wasm/)

## Decisions:
- All new code behind feature flags `alpha-dynamics-live`, `alpha-dft`, etc.
- Reuse existing infrastructure: forcefield traits, SCF loop, dynamics engine
- No breaking changes to existing stable code
