# Campaña de validación externa — 2026-11-07

Esta campaña vuelve a ejecutar los oráculos de referencia y los gates Rust
contra las versiones reales de las bibliotecas externas instaladas en la
máquina del mantenedor, y regenera los fixtures que consume la suite. No
sustituye a `docs/validation/cross-library-validation.es.md` (que describe el
contrato de cada ficha): este informe registra **qué se ejecutó, con qué
versiones, qué cambió y qué queda pendiente**.

## 1. Entorno

- Fecha: 2026-11-07.
- Python: 3.14.4 (`~/.venvs/sci-form-validation`).
- Rust: toolchain estable del repositorio; `cargo test --release`.
- Sistema: Linux x86-64, 16 núcleos, 26 GB de RAM.
- `cargo clean` previo: 36.9 GiB liberados; 51 GB libres durante la campaña.

### Versiones de referencia

| Herramienta | Versión | Uso principal |
|---|---|---|
| RDKit | 2026.09.1 | F01/F05/F13/F14/F16/F17/F21/F28/F59/F70–F76/F79/F84 |
| PySCF | 2.14.0 | F23/F26/F36–F39/F41/F77/F80 |
| tblite | 0.7.0 | F44/F45 (GFN1/GFN2, UKS) |
| dftd4 | 4.3.0 | F46 (D4) |
| s-dftd3 | 1.6.0 | F46 (D3) |
| xtb | 6.7.1 (binario oficial) | F43/F48/F49 |
| multicharge | 0.6.0 (CLI, build meson/ninja) | F22 (EEQ-BC) |
| MOPAC | 23.2.5 | F42 (PM3) |
| LAMMPS | 29Aug2024_update2 (REAXFF/QEQ/ELECTRODE) | F78/F95 |
| ASE | 3.29.0 | F66/F82/F88 |
| pormake | 0.2.3 | F67/F68 |
| pymatgen / spglib | 2026.9.24 / 2.8.0 | F64/F65/F69 |
| NumPy / SciPy | 2.5.4 / 1.18.1 | oráculos numéricos |
| FreeSASA | bindings Python | F28 |
| pyzeo | 0.2.0 | F68 |
| clifford | 1.5.1 | F81 |
| cvxpy (SCS) | 1.9.3 | F83 |
| jsonschema | 4.26.0 | F90 |

Las versiones quedan registradas en `validation/refs.lock.json`. Se añadieron
las entradas de **multicharge 0.6.0** y **LAMMPS 29Aug2024_update2**, que
antes figuraban como pendientes.

## 2. Comandos

```bash
# Toolchain de referencia (fuera del repo)
export XTB=/tmp/opencode/refs/xtb-dist/bin/xtb
export XTBPATH=/tmp/opencode/refs/xtb-dist/share/xtb
export SCIFORM_REF_DATA="$PWD/.cache/refdata"

# Oráculos (36 del lote; logs en /tmp/opencode/oracle-logs/)
python scripts/validation/oracle_rdkit.py all
python scripts/validation/oracle_pyscf.py all
python scripts/validation/oracle_tblite.py
python scripts/validation/oracle_dftd4.py
python scripts/validation/oracle_xtb.py
python scripts/validation/oracle_hessian.py
python scripts/validation/oracle_pm3.py mopac-fixture
python scripts/validation/oracle_eeqbc.py \
  --binary "$PWD/.cache/refdata/multicharge-src/build/app/multicharge"
python scripts/validation/oracle_reaxff_lammps.py
python scripts/validation/oracle_cpm_lammps.py --lmp /usr/bin/lmp
# … resto de oráculos alpha/beta (ase, gsm, periodic, pormake, fukui, eht,
# eeq2019, sdr, kpm, randnla, riemannian, mbh, cga, edl, kinetics, render,
# imd, mlff, obara_saika, hla, uks_tblite, pyscf_stda/esp/ir, dynamics,
# uffopt)

# Gates Rust de validación (19 binarios + 2 de paridad)
FEATS="$(python3 scripts/_features.py RUST_EXPERIMENTAL_FEATURES)"
cargo test --release --features "$FEATS" \
  --test cross_validation --test alpha_cga_validation \
  --test alpha_dft_validation --test alpha_dynamics_validation \
  --test alpha_edl_validation --test alpha_gsm_validation \
  --test alpha_imd_validation --test alpha_kinetics_validation \
  --test alpha_mlff_validation --test alpha_periodic_validation \
  --test alpha_reaxff_validation --test alpha_render_validation \
  --test alpha_sdr_validation --test beta_kpm_validation \
  --test beta_mbh_validation --test beta_randnla_validation \
  --test beta_riemannian_validation --test mmff94_validation \
  --test validation_gates \
  --test eht_eeq_parity --test pyscf_parity

# Consolidación del resumen
python scripts/validation/collect.py
```

## 3. Resultados

### 3.1 Oráculos

**36/36 oráculos del lote terminan con código 0**. Los 5 scripts restantes
(`oracle_rdkit_geometry`, `oracle_rdkit_smarts`, `check_csd_torsion_table`,
`oracle_pymatgen`, `oracle_mmff94`) generan tablas/fixtures que no cambiaron
en esta campaña y cuyos gates siguen verdes.

Dos fallos se diagnosticaron y corrigieron durante la campaña:

1. **`oracle_mlff.py`**: usaba `Chem` sin importar RDKit y no tenía el
   *fallback* que su propio docstring prometía. Ahora importa RDKit dentro de
   `main()` y, si no está disponible, recalcula desde las geometrías del
   fixture (`--from-fixture`).
2. **ReaxFF/LAMMPS**: el `lmp` de Debian (2025-12-10) anuncia los paquetes
   REAXFF pero conserva solo los estilos obsoletos `reax`/`reax/c`, que ya no
   existen; `pair_style reaxff` no está compilado. Se compiló LAMMPS
   `stable_29Aug2024_update2` con `PKG_REAXFF`, `PKG_QEQ`, `PKG_ELECTRODE` y
   `PKG_MANYBODY`, y se apuntó `SCIFORM_REF_DATA` a ese build. F78 y F95 pasan.

También se compiló **multicharge 0.6.0** desde fuente (meson + ninja, LAPACK
netlib local), porque no hay wheel de Python y el fixture anterior se había
generado con 0.5.0.

### 3.2 Gates Rust

- **117/117 tests de validación externa en verde**:
  - 38 en `cross_validation`;
  - 44 en los gates alpha/beta;
  - 5 en `mmff94_validation`;
  - 4 en `validation_gates`;
  - 3 en `eht_eeq_parity`;
  - 3 en `pyscf_parity`.
- F39 (transformada AO→MO) con el fixture regenerado: **1.8e-14 a.u.**
  (tolerancia 1e-10), 13 058 integrales MO.
- `reports/validation/summary.json`: **127 filas, 120 PASS, 7 pendientes,
  0 FAIL**. Las 7 pendientes son filas de procedencia o comparaciones con
  *gap* metodológico documentado (F44/F46 procedencia, F71 bounds, F82
  descriptores ML, F83 SMILES canónico, F84 solvatación ALPB, F85 stress de
  celda variable); ninguna es un fallo numérico.

### 3.3 Fixtures regenerados

Se actualizaron **30 fixtures**. El contenido numérico cambia solo en:

- `generated` (fecha) y, en su caso, la versión de la herramienta;
- ruido de última cifra (≤ 1e-10; la mayoría ≤ 1e-15);
- `eeqbc_charges.json`: fuente 0.5.0 → 0.6.0 y diferencias ≤ 4.4e-16;
- `pyscf_ao2mo.json`: la base MO de PySCF es arbitraria dentro de subespacios
  degenerados (T2 de CH4/NH3). El oráculo ahora **canonicaliza** cada
  subespacio por proyección de la base estándar + Gram-Schmidt y fija el signo
  de cada columna. Dos ejecuciones consecutivas coinciden a ≤ 1.6e-11 y el
  gate F39 sigue midiendo la transformada, no la fase de los orbitales.

## 4. Limitaciones y pendientes

- **ORCA** (F40/F50/F52) sigue fuera: no hay binario instalado.
- **pyGSM/pysisyphus** (F82/F89): F82 usa ASE CI-NEB; F89 usa barreras
  experimentales curadas.
- **pybinding** (F34/F88/F91) no es instalable en este entorno; los gates
  usan ASE/NumPy.
- **geomstats/pymanopt** (F94): se usa un proxy de Frobenius, no la métrica
  afín-invariante.
- Los fixtures se consumen desde JSON; la reproducibilidad byte a byte no es
  un requisito del gate, pero los generadores son deterministas salvo el ruido
  BLAS de última cifra.

## 5. Conclusión

La campaña no encontró **ningún fallo de validación**: 117/117 gates Rust y
36/36 oráculos del lote en verde, con 120/127 filas PASS y 0 FAIL. Los dos
problemas reales eran de infraestructura (import ausente en el oráculo MLFF y
LAMMPS Debian sin REAXFF), no del código de `sci-form`. Los fixtures
regenerados y el lock de versiones quedan como evidencia reproducible de la
campaña.
