"""Cross-library validation infrastructure for sci-form.

This package contains the *oracle* scripts that generate reference data from
established implementations (RDKit, PySCF, ASE, ...) and compare them against
sci-form. Design rules (see docs/validation/cross-library-validation.es.md):

1. Oracles only *generate* fixtures and reports; they never modify src/.
2. Scient-form results are obtained through the CLI or the Rust test suite.
3. Every metric lives in a report of the shape:
       {ficha, observable, reference, units, n_cases, metric, value,
        tolerance, pass, notes}
4. Fixtures are stored as JSON under tests/fixtures/validation/.
5. Summary reports are appended to reports/validation/summary.json.

Usage:
    python scripts/validation/oracle_rdkit.py parse
    python scripts/validation/oracle_pyscf.py scf
"""
