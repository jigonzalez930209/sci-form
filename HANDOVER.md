# Handover — `feat/wave-c`: gates v2 y filas de validación

> Estado al 2026-11-19. Rama: `feat/wave-c` (no fusionada en `main`).
> Último commit de trabajo: ver `git log`. Borra este archivo al fusionar.

## 1. Contexto

`feat/wave-c` contiene los 15 candidatos v2.x (V2-01…V2-15) de
`docs/v2-roadmap.md`. `feat/wave-a` y `feat/wave-b` son ancestros suyos: no hay
que fusionarlas por separado. Las ramas de marzo se borraron del remoto.

## 2. Qué está hecho (y subido)

| Commit | Qué |
|---|---|
| `6d784a4` fix(dft) | PBE0 contaba el intercambio dos veces; el PBE de la malla se escala por `1 − a` (`pbe_with_exchange_scale`) |
| `3139dd9` feat(v2) | Gates `alpha-*` por candidato + paraguas `alpha-v2`; bindings CLI `pbe0`, `nmr-all`, `nmr-isotopes`; WASM `compute_pbe0_energy`; workflows con `alpha-v2`; roadmap §5.1 en una tabla; línea **Feature** en cada página de algoritmo |
| `5e68494` fix(raman) | Polarizabilidad UCHF sin el factor 2 (α era la mitad) |
| `5b68aef` fix(rism) | Clausura KH ramificaba en `h − c` en vez de `d = −βu + h − c` |
| `d5b6597` test(validation) | Oráculos `scripts/validation/oracle_pyscf_v2.py` y `oracle_pymbar_v2.py` + fixtures `tests/fixtures/validation/v2_pyscf.json`, `v2_pymbar.json` |

Mapa de features (core, Python, WASM y CLI usan los mismos nombres):

| Feature | Módulos |
|---|---|
| `alpha-dft-hybrid` (activa `alpha-dft`) | `DftMethod::Pbe0` |
| `alpha-mp2` / `alpha-qmmm` / `alpha-epr` / `alpha-rism` / `alpha-gcmc` | `mp2` / `qmmm` / `epr` / `rism` / `gcmc` |
| `alpha-raman` | `raman`, `vcd_roa` |
| `alpha-core-spectroscopy` | `xps` |
| `alpha-redox-pka` | `electrochemistry` |
| `alpha-free-energy` | `free_energy` |
| `alpha-enhanced-sampling` | `enhanced_sampling`, `enhanced_sampling_multi` |
| `alpha-reaction-network` | `networks`, `kmc` |
| `alpha-mlip` (`mlip-onnx` lo activa) | `mlip` (+ `mlip_onnx`) |
| `alpha-active-learning` | `active_learning` |
| `alpha-nmr-all` | `nmr_all` |

Verificado: `cargo check` de core + bindings con y sin `alpha-v2` limpio;
tests de los 17 módulos + `dft::{ks_fock,vxc_matrix,functionals}` → 95 ok;
`raman/vcd_roa/rism` tras los fixes → 12 ok; test de catálogo de la CLI ok.
**No** se ejecutó: `every_example_runs` de la CLI, tests Python (maturin),
build `wasm32`.

## 3. Estado de la validación v2 (2026-11-19, sesión 2)

### 3.1 Hecho — `tests/v2_validation.rs` (F101–F115)

`scripts/cargo-guard.sh test --profile numeric --features alpha-v2 --test v2_validation`
→ 15/15 tests; `collect.py` → 23 filas v2: **20 PASS, 3 gap**, 0 fail.
`collect.py` copia ahora `status/level/units/reference/notes` y cuenta los gap
aparte. Narrativa en `docs/validation/scientific-validation.md` (sección
"v2.x candidates"); matriz regenerada (156 filas).

Bugs encontrados por los gates y corregidos (cada uno con test de regresión):

| Módulo | Defecto |
|---|---|
| `scf` (RHF/UHF) | Guess H_core → estados excitados (NH2 ²A₁, O2 +0.255 Ha); ahora guess SAD. F38 compara con `uhf_default` de PySCF |
| `scf` (RHF/UHF/FON) | La energía usaba la Fock extrapolada por DIIS (y con level shift en RHF); ahora E y F son las de la densidad devuelta (F41 osc. 4.1e-3 → 9.2e-7) |
| `dft::lebedev` | Grid 26 (`Medium`, el default) con pesos arista/cubo intercambiados; grid 110 (`Fine`) no era Lebedev (PBE −0.25 Ha). Ahora Lebedev–Laikov exacto hasta grado 7/17 |
| `gcmc` | El borrado aceptaba ΔU ≤ 0 sin el factor N/(βfV): gas ideal a ~0.6 PV/kT |
| `rism` | Clausura `c = g−1−h` en vez de `c = g−1−t`; ahora itera sobre t = h − c |
| `enhanced_sampling(_multi)` | FES WT-metaD/OPES era `−V`; ahora `−γ/(γ−1)·V` |
| `free_energy` | FEP rechazaba muestras solo del estado 0 |
| `mlip` | `CallbackMlip` ignoraba `energy_unit: ev` |

Gaps registrados (valor medido, `pass: null`):

| Fila | Valor | Causa |
|---|---|---|
| F103 embedding electrostático (L3) | 4.2 kcal/mol | sin polarización MM→QM |
| F105 CPHF (L3) | 21 % | sin respuesta acoplada |
| F108 pKa (L4) | MAE 6.6 | single-point sobre geometrías ETKDG sin optimizar; con 8 confórmeros cambia hasta 8 unidades |

### 3.2 Decisiones abiertas

- **pKa**: optimizar geometría (QM o FF) antes del single-point.
- Grid 302 (`VeryFine`) sigue siendo una espiral de Fibonacci (aproximada,
  pero F77/F101 dan < 1e-4 Ha); sustituir por Lebedev 302 real.
- `MetaDynamicsConfig::default().friction_ps_inv = 1.0` deja al caminante
  caliente con la deposición por defecto (documentado; F110 usa 10 ps⁻¹).

### 3.3 Pendientes menores

- Enhanced sampling no tiene comando CLI (anotado en roadmap §5.1).
- La CLI de release se compila sin features (`release.yml`), como ya pasa con
  `alpha-edl`: los comandos v2 no están en ese binario.
- Física pendiente por candidato: ver columna "Remaining before promotion" de
  roadmap §5.1 (RI-MP2, ΔSCF/XANES, energía libre RISM, embedding
  electrostático, tensor g, MD soft-core, OPES exacto, MBH/atom mapping…).

## 4. Datasets para los gates L4 (encontrados, aún no descargados)

| Candidato | Dataset | Fuente |
|---|---|---|
| V2-06 XPS | CORE65 (65 BE 1s, 32 moléculas; exp. de Jolly et al. 1984) | Golze et al. JPCL 2020, Table S1 / NOMAD "CORE65 benchmarkset" — https://arxiv.org/abs/1911.08428 |
| V2-07 RISM, V2-09 FEP | FreeSolv (ΔG_hyd exp.) | https://github.com/jchodera/FreeSolv |
| V2-08 pKa | SAMPL6 pKa (24 moléculas) | https://github.com/samplchallenges/SAMPL6 |
| V2-08 redox | ROP313 (313 E° exp.) | Neugebauer, Grimme et al. 2020 |
| V2-15 RMN | NS372 (372 apantallamientos CCSD(T)) | Schattenberg & Kaupp JCTC 2021; estructuras en DepositOnce (TU Berlin) |
| V2-13 MLIP | rMD17 (CC0) y SPICE (CC-BY 4.0) | https://figshare.com/articles/dataset/Revised_MD17_dataset_rMD17_/12672038 · https://doi.org/10.5281/zenodo.7338495 (verificar) |
| V2-11 GCMC | NIST/ARPA-E ISODB | https://github.com/NIST-ISODB/isodb-library |
| V2-12 barreras | BH76 (GMTKN55) | https://chemie.uni-bonn.de/grimme/de/software/gmtkn/gmtkn55 |
| V2-04 EPR | Benchmarks QTP (56 radicales) / Hermosilla 2015 | https://impact.ornl.gov/en/publications/benchmarking-isotropic-hyperfine-coupling-constants-using-qtp-dft/ |

**No meter valores experimentales de memoria.** Los de CORE65 por molécula no
se pudieron verificar desde una fuente abierta.

## 5. Entorno

```bash
# oráculos (documentado en scripts/run_validation.sh)
uv venv .venv-validation -p 3.12
uv pip install -p .venv-validation/bin/python pyscf pymbar numpy scipy
.venv-validation/bin/python -I scripts/validation/oracle_pyscf_v2.py   # ~150 s
.venv-validation/bin/python -I scripts/validation/oracle_pymbar_v2.py
```

- `.venv-validation/` está en `.gitignore`.
- `pyscf-properties` (GIAO/EPR en PySCF) se instala desde git; en la sesión
  anterior esa instalación fue denegada por permisos → por eso F115 es L1.
- Preferencia del usuario: **no correr toda la suite**; solo los tests de los
  módulos tocados.
- **Guard de RAM**: compilar/testear siempre con `scripts/cargo-guard.sh`
  (scope systemd con MemoryHigh/MemoryMax, 8 jobs) y `--profile numeric`
  (release sin LTO fat). `cargo test --release` con LTO tumbó la sesión por OOM.
  Ej.: `scripts/cargo-guard.sh test --profile numeric -p sci-form --lib --features alpha-v2 -- <modulo>::`
