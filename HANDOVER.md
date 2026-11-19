# Handover — `feat/wave-c`: gates v2 y filas de validación

> Estado al 2026-11-19. Rama: `feat/wave-c` (no fusionada en `main`).
> Último commit de trabajo: `d5b6597`. Borra este archivo al fusionar.

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

## 3. Qué falta (siguiente paso inmediato)

### 3.1 `tests/v2_validation.rs` — NO existe todavía

Fichas **F101–F115 ↔ V2-01…V2-15** (F100 era el último ID usado;
`collect.py` solo recoge `target/validation/F*.json`). Declarar el test en
`Cargo.toml` con `required-features = ["alpha-v2"]`.

**Regla de registro (decidida):** `collect.py` devuelve 1 si alguna fila tiene
`pass: false`. Por eso:
- Gates L0–L2 que el prototipo debe cumplir → `assert!` + `pass: true`.
- Gates L3/L4 que el prototipo aún no cumple → registrar el valor medido con
  `pass: null`, `status: "gap"`. Nunca inventar un pass.
- Hay que **extender `scripts/validation/collect.py`** para copiar campos
  opcionales de la fila (`status`, `level`, `units`, `reference`, `notes`);
  hoy solo copia ficha/observable/metric/value/tolerance/pass/n_cases/source.
- Helper `record(...)` propio en el test (copiar el de `tests/pyscf_parity.rs`
  pero con `pass: Option<bool>` y los campos extra).

Plan por ficha (oráculo → API sci-form):

| Ficha | Gate | Nivel | Cómo |
|---|---|---|---|
| F101 | PBE y PBE0 total vs PySCF (`v2_pyscf.json["dft"]`) | L2 → probablemente **gap** | `dft::ks_fock::solve_ks_dft` con `DftMethod::Pbe`/`Pbe0`, base `Minimal`. La fila PBE aísla si el error viene del motor PBE ("PW92 simplificado") o del híbrido. PySCF agua: PBE −75.225377, PBE0 −75.245543; sci-form PBE0 agua dio −75.0743 |
| F102 | MP2 correlación vs PySCF (`["mp2"]`) | L2, tol 1e-6 Ha | `compute_mp2_energy(elements, positions)`; todos los electrones |
| F103 | QM/MM mecánico vs referencia PySCF (`["qmmm"]["mechanical_total_kcal_mol"]`) | L2, tol ~1e-3 kcal/mol | `compute_qmmm_energy` (UHF, `include_mm_mm: true`). Gap L3: `electrostatic_total_kcal_mol` (−4.2 kcal/mol de polarización ausente) |
| F104 | A_iso vs densidad de espín UHF de PySCF (`["epr"]`, H, CH3, NH2, OH) | L2 | `compute_epr_hyperfine`. Gap L4: átomo H exacto 1420.405751 MHz (STO-3G da ~1761) |
| F105 | α isotrópica no acoplada vs PySCF (`["polarizability"]`) | L2 | `raman::compute_polarizability`. Gap L3: `coupled_isotropic_ang3` (CPHF campo finito) |
| F106 | BE Koopmans vs PySCF (`["xps"]`, emparejar por `atom_index`) | L2 | `compute_xps_binding_energies` |
| F107 | RISM límite ideal: χ = 0 ⇒ h = c ⇒ g = KH(−βu) exacto | L0, tol 1e-8 | `solve_3d_rism` con pozo atractivo **y** repulsivo, χ = 0 |
| F108 | pKa relativo a ácido acético (4.76) | L4 → gap si MAE > 1 | `compute_pka`, UHF, 1 conf. Exp (CRC Handbook): fórmico 3.75, propanoico 4.87, cloroacético 2.87, benzoico 4.20 |
| F109 | FEP/BAR/MBAR vs pymbar (`v2_pymbar.json`) | L2, tol 1e-6 | `free_energy::estimate_free_energy`, `n_states: 3` (BAR: 2 estados, muestras de 0 y 1). L0 vs `exact_delta_f` con tolerancia estadística (MBAR pymbar: [0, 0.3412, 0.7032] vs exacto [0, 0.3466, 0.6931]) |
| F110 | Barrera WT-metaD en `DoubleWell::new(3.0, 1.0)` vs exacta | L0, tol 1 kcal/mol | `run_metadynamics`; CV = coordenada ⇒ F(x) = U(x). Ver config de ejemplo en tests de `enhanced_sampling.rs` (~l.1577) |
| F111 | GCMC gas ideal ⟨N⟩ = PV/kT | L0, 5 % | `GcmcFramework::empty_cubic`, `GuestModel::single_site(0,3,0)`; ver test `ideal_gas_loading_follows_pv_over_kbt` |
| F112 | kMC A→B: ⟨N_A(t=1/k)⟩/N₀ vs e⁻¹ | L0, 5 % | `run_kinetic_monte_carlo` con `max_time_s = 1/k`, ~200 semillas; k de `kinetics::evaluate_htst_rate(..).forward_rate_s_inv` |
| F113 | Contrato MLIP: callback en eV → kcal/mol y F = −∇E por diferencias finitas | L0 | `mlip::CallbackMlip` + `MlipEnergyUnit::Ev`; fila ONNX solo con `cfg(feature = "mlip-onnx")` |
| F114 | Cobertura conformal empírica 0.90 ± 0.02 en datos sintéticos | L0 | `active_learning::ConformalCalibrator::fit`, intervalo = pred ± `radius` |
| F115 | Driver all-nuclide == GIAO directo por núcleo | L1 | `nmr_all::compute_all_nmr` vs `compute_giao_nmr_configured` |

Luego: `cargo test --features alpha-v2 --test v2_validation` y
`python3 scripts/validation/collect.py` para escribir las filas en
`reports/validation/summary.json`; añadir la narrativa a
`docs/validation/scientific-validation.md` (lo pide el roadmap §8).

### 3.2 Pendientes menores

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
  módulos tocados (`cargo test -p sci-form --lib --features alpha-v2 -- <modulo>::`).
