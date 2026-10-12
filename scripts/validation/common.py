"""Shared helpers for cross-library validation oracles.

Dependencies: numpy (required).  RDKit/PySCF/... are only imported by the
individual oracle scripts that need them.
"""

from __future__ import annotations

import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Paths and constants
# ---------------------------------------------------------------------------

REPO = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO / "tests" / "fixtures" / "validation"
REPORT_DIR = REPO / "reports" / "validation"
SUMMARY_PATH = REPORT_DIR / "summary.json"
DATA_DIR = REPO / "data"

# Unit conversion factors
HARTREE2KCAL = 627.5094740631
HARTREE2EV = 27.211386245988
EV2KCAL = 23.060548867
BOHR2ANG = 0.529177210903
ANG2BOHR = 1.0 / BOHR2ANG
DEBYE_AU = 2.541746473
# e·Angstrom -> Debye conversion. NOTE: this intentionally matches
# sci-form's `EANG_TO_DEBYE = 4.80321` (not the CODATA value) so that
# point-charge dipoles agree to machine precision in F26.
DEBYE_EANG = 4.80321

# Tolerance defaults (mirrors docs/validation/cross-library-validation.es.md §3.3)
TOL = {
    "smiles_parse": 0.0,
    "sssr_multiset": 0.0,
    "smarts_matches": 0.0,
    "stereo_label": 0.0,
    "gasteiger_charge": 0.02,
    "lowdin_charge": 0.01,
    "mulliken_charge": 0.005,
    "ecfp_bits": 0.0,
    "tanimoto": 1e-9,
    "sasa_rel": 0.01,
    "uff_energy": 0.01,
    "mmff94_energy": 0.01,
    "rmsd": 1e-6,
    "embedding_rmsd_median": 0.25,
    "embedding_rmsd_p90": 0.5,
    "scf_energy": 1e-7,
    "scf_mo_energy": 1e-6,
    "overlap_element": 1e-10,
    "eri_element": 1e-10,
    "dos_mse": 0.01,
    "kpm_rel": 0.05,
    "dipole_debye": 0.1,
}


def cli_path() -> Path:
    """Return the sci-form CLI binary path (release build preferred)."""
    env = os.environ.get("SCI_FORM_CLI")
    if env:
        return Path(env)
    for candidate in (
        REPO / "target" / "release" / "sci-form",
        REPO / "target" / "debug" / "sci-form",
    ):
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "sci-form CLI not found; run `cargo build --release -p sci-form-cli` "
        "or set SCI_FORM_CLI"
    )


def run_cli(args: Sequence[str], timeout: float = 300.0) -> dict:
    """Run the sci-form CLI and parse its JSON stdout.

    Raises RuntimeError with the captured stderr/stdout on failure.
    """
    cmd = [str(cli_path()), *args]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if proc.returncode != 0:
        raise RuntimeError(
            f"CLI failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout: {proc.stdout[:500]}\nstderr: {proc.stderr[:500]}"
        )
    out = proc.stdout.strip()
    if not out:
        return {}
    try:
        return json.loads(out)
    except json.JSONDecodeError as exc:
        raise RuntimeError(f"CLI did not return JSON: {out[:300]}") from exc


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def fixture_path(name: str) -> Path:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    return FIXTURE_DIR / name


def write_fixture(name: str, payload: Any) -> Path:
    path = fixture_path(name)
    with path.open("w") as fh:
        json.dump(payload, fh, indent=1)
    return path


def load_fixture(name: str) -> Any:
    with fixture_path(name).open() as fh:
        return json.load(fh)


def read_smiles_corpus(limit: int | None = None, extra: Iterable[str] = ()) -> list[str]:
    """Corpus for graph-level comparisons: diverse molecules + ChEMBL sample."""
    smiles: list[str] = []
    diverse = REPO / "tests" / "fixtures" / "diverse_molecules.smi"
    for path in (diverse, DATA_DIR / "chembl_1k_practical.smi"):
        if not path.exists():
            continue
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            smiles.append(line.split()[0])
    smiles.extend(extra)
    # Deduplicate, preserve order
    seen: set[str] = set()
    unique = []
    for smi in smiles:
        if smi not in seen:
            seen.add(smi)
            unique.append(smi)
    if limit is not None:
        unique = unique[:limit]
    return unique


# ---------------------------------------------------------------------------
# Numeric helpers
# ---------------------------------------------------------------------------


def kabsch_rmsd(p: np.ndarray, q: np.ndarray) -> float:
    """Kabsch-aligned RMSD between two (N,3) coordinate arrays."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    if p.shape != q.shape:
        raise ValueError(f"shape mismatch: {p.shape} vs {q.shape}")
    pc = p - p.mean(axis=0)
    qc = q - q.mean(axis=0)
    h = pc.T @ qc
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    dmat = np.diag([1.0, 1.0, d])
    rot = vt.T @ dmat @ u.T
    aligned = pc @ rot.T
    diff = aligned - qc
    return float(np.sqrt((diff * diff).sum() / len(p)))


def max_abs_diff(a: Sequence[float], b: Sequence[float]) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    return float(np.max(np.abs(a - b))) if a.size else 0.0


def rmse(a: Sequence[float], b: Sequence[float]) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.sqrt(np.mean((a - b) ** 2)))


def rel_error(value: float, reference: float) -> float:
    if reference == 0.0:
        return abs(value - reference)
    return abs(value - reference) / abs(reference)


def percentile(values: Sequence[float], p: float) -> float:
    if not values:
        return float("nan")
    return float(np.percentile(np.asarray(values, dtype=float), p))


# ---------------------------------------------------------------------------
# Reports
# ---------------------------------------------------------------------------


@dataclass
class Report:
    """One validation metric (see docs §11 for the schema)."""

    ficha: str
    observable: str
    reference_tool: str
    units: str
    tolerance: float
    level: str = "L2"
    n_cases: int = 0
    metric: str = "max_abs"
    value: float = float("nan")
    notes: str = ""
    cases: list[dict] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    @property
    def passed(self) -> bool:
        try:
            return bool(self.value <= self.tolerance)
        except TypeError:
            return False

    def to_dict(self) -> dict:
        return {
            "ficha": self.ficha,
            "observable": self.observable,
            "reference": {"tool": self.reference_tool},
            "units": self.units,
            "level": self.level,
            "n_cases": self.n_cases,
            "metric": self.metric,
            "value": None if math.isnan(self.value) else self.value,
            "tolerance": self.tolerance,
            "pass": self.passed,
            "notes": self.notes,
            **({"extra": self.extra} if self.extra else {}),
        }

    def emit(self, verbose: bool = True) -> "Report":
        REPORT_DIR.mkdir(parents=True, exist_ok=True)
        summary: list[dict]
        if SUMMARY_PATH.exists():
            try:
                summary = json.loads(SUMMARY_PATH.read_text())
            except json.JSONDecodeError:
                summary = []
        else:
            summary = []
        payload = self.to_dict()
        payload["timestamp"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        summary = [
            row
            for row in summary
            if not (row.get("ficha") == self.ficha and row.get("observable") == self.observable)
        ]
        summary.append(payload)
        SUMMARY_PATH.write_text(json.dumps(summary, indent=1))
        if verbose:
            status = "PASS" if self.passed else "FAIL"
            print(
                f"[{status}] {self.ficha} {self.observable}: "
                f"{self.metric}={self.value:.6g} tol={self.tolerance:g} "
                f"(n={self.n_cases}, ref={self.reference_tool})"
            )
        return self


def emit_fixture_report(
    name: str,
    payload: Any,
    report: Report,
    verbose: bool = True,
) -> Path:
    """Write a fixture and emit its report."""
    path = write_fixture(name, payload)
    report.emit(verbose=verbose)
    return path
