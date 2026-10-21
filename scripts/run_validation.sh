#!/usr/bin/env bash
# Cross-library validation runner.
#
# See docs/validation/cross-library-validation.es.md for the fichas, tolerances
# and datasets. Requires the oracle toolchain in $PYTHON:
#
#   python3 -m venv .venv-validation
#   .venv-validation/bin/pip install rdkit freesasa pyscf pyarrow scikit-learn tblite dftd4 dftd3 pymatgen ase spglib
#   VALIDATION_PYTHON=.venv-validation/bin/python scripts/run_validation.sh
#
set -euo pipefail
cd "$(dirname "$0")/.."

PYTHON="${VALIDATION_PYTHON:-python3}"
FAST="${VALIDATION_FAST:-0}"

step() { printf '\n== %s ==\n' "$1"; }

step "[1/5] build CLI"
cargo build --release -p sci-form-cli

step "[2/5] generate RDKit fixtures + reports (${PYTHON})"
"${PYTHON}" scripts/validation/oracle_rdkit.py all
"${PYTHON}" scripts/validation/oracle_uffopt.py

step "[3/6] generate PySCF fixtures (optional)"
if "${PYTHON}" -c "import pyscf" >/dev/null 2>&1; then
  "${PYTHON}" scripts/validation/oracle_pyscf.py all
  "${PYTHON}" scripts/validation/oracle_pyscf_stda.py
else
  echo "pyscf not available in ${PYTHON}; skipping (F36/F37/F23/F38/F39 use the checked-in fixtures)"
fi

step "[4/7] generate materials fixtures (optional)"
if "${PYTHON}" -c "import pymatgen" >/dev/null 2>&1; then
  "${PYTHON}" scripts/validation/oracle_pymatgen.py all
else
  echo "pymatgen not available in ${PYTHON}; skipping (F64/F65/F66/F69 use the checked-in fixtures)"
fi
if "${PYTHON}" -c "import pormake, pyzeo" >/dev/null 2>&1; then
  "${PYTHON}" scripts/validation/oracle_pormake.py
else
  echo "pormake/pyzeo not available in ${PYTHON}; skipping (F67/F68 use the checked-in fixture)"
fi

step "[5/7] generate tblite/dftd4/xtb fixtures (optional)"
if "${PYTHON}" -c "import tblite" >/dev/null 2>&1; then
  "${PYTHON}" scripts/validation/oracle_tblite.py
else
  echo "tblite not available in ${PYTHON}; skipping (F44/F45 use the checked-in fixture)"
fi
if "${PYTHON}" -c "import dftd4" >/dev/null 2>&1; then
  "${PYTHON}" scripts/validation/oracle_dftd4.py
else
  echo "dftd4 not available in ${PYTHON}; skipping (F46 uses the checked-in fixture)"
fi
if command -v xtb >/dev/null 2>&1 || [[ -n "${XTB:-}" ]]; then
  "${PYTHON}" scripts/validation/oracle_xtb.py
  "${PYTHON}" scripts/validation/oracle_hessian.py
else
  echo "xtb binary not found (set XTB=/path/to/xtb); skipping (F43/F48 use the checked-in fixtures)"
fi

step "[6/7] generate Fukui/PM3/EHT/EEQ2019 fixtures"
"${PYTHON}" scripts/validation/oracle_fukui.py
"${PYTHON}" scripts/validation/oracle_pm3.py
# numpy/scipy only (no external tool needed)
"${PYTHON}" scripts/validation/oracle_eht.py
"${PYTHON}" scripts/validation/oracle_eeq2019.py

step "[7/7] run cross-validation tests"
if [[ "${FAST}" == "1" ]]; then
  CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 \
    cargo test --test cross_validation -- --nocapture
else
  cargo test --test cross_validation -- --nocapture
fi

step "[8/8] collect summary"
"${PYTHON}" scripts/validation/collect.py
