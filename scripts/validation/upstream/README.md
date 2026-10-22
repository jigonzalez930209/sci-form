# Vendored upstream inputs

Third-party parameter sources that the generators in `scripts/validation/`
parse.  They are vendored so that the generated Rust tables can be **verified
offline** (`python3 scripts/validation/gen_<x>.py --check` in CI) instead of
depending on a network fetch at run time.

## `mmff94/Params.cpp`

* **Upstream**: <https://raw.githubusercontent.com/rdkit/rdkit/master/Code/ForceField/MMFF/Params.cpp>
  (RDKit, BSD-3-Clause).
* **SHA-256**: `9e6c21cf735e1fefe9083d04c00b61d2f2eb105297a3c8d1e8ecb3f40fa16700`
  — the same digest that `src/forcefield/mmff94_params.rs` records in its
  header, so the provenance of every number stays traceable.
* **Data licence**: the MMFF94 parameter set itself is
  Copyright (c) Merck and Co., Inc., 1994, 1995, 1996 (the published set of
  Halgren, *J. Comput. Chem.* **1996**, *17*, 490-519).

Regenerate/verify with:

```bash
python3 scripts/validation/gen_mmff94_params.py --check
```

The generator still accepts `--params-cpp PATH` to point at another copy, and
falls back to downloading the file when the vendored copy is missing.

## The other generators

`gen_d4_refs.py`, `gen_eeq2019_params.py`, `gen_eeqbc_params.py`,
`gen_gfn0_params.py`, `gen_pm3_params.py` and `gen_topologies.py` need
upstream source trees (dftd4, grimme-lab/multicharge, grimme-lab/mctc-lib,
xtb, MOPAC, pormake) that are **not** vendored here — they are whole
third-party repositories with their own licences.

For those, point `SCIFORM_REF_DATA` at a directory that holds the checkouts:

```bash
export SCIFORM_REF_DATA=/path/to/refdata
python3 scripts/validation/gen_d4_refs.py --check
```

When the input is missing, the `--check` mode exits with code **2** ("skipped")
instead of failing, so the CI gate reports the gap without hiding a real
drift (exit code **1**) in the tables it can verify.
