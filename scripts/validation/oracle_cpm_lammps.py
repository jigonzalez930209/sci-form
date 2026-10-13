#!/usr/bin/env python3
"""Generate `tests/fixtures/validation/cpm_lammps.json` from LAMMPS ELECTRODE.

Runs `fix electrode/conp` (constant-potential method, Siepmann/Reed) on a
graphene bilayer capacitor (two 3x3-cell flakes, AB stacked at 3.35 A, 40 A
cubic box, p p p) at several symmetric voltages and records per-atom and
per-plate charges.

Findings baked into the companion gate (`tests/experimental/test_cpm.rs`):
- LAMMPS needs two coupled plates (`couple`); a lone electrode group stays
  neutral (cell neutrality, `symm on` default).
- `run 0` leaves charges at zero; the fix solves during stepping (`run 1`
  with nve) — data come from the stepped frame.
- eta = 1.979 1/A (graphite value from the LAMMPS electrode example).

Our CPM (`beta::cpm`, per-atom mu + direct LU solve) is run in-test; the
fixture holds the LAMMPS reference + geometry only.

Usage:
    python scripts/validation/oracle_cpm_lammps.py [--out FILE] [--lmp LMP]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "tests/fixtures/validation/cpm_lammps.json"

A_CC = 2.46  # graphene lattice constant (A)
SEPARATION = 3.35  # AB interlayer distance (A)
BOX = 40.0
ETA = 1.979  # Gaussian reciprocal width (1/A, graphite)
VOLTAGES = [0.5, 1.0, 2.0]
NCELL = 3  # 3x3 unit cells per plate (18 atoms)

INPUT_TEMPLATE = """units real
atom_style charge
boundary p p p
pair_style coul/long 12.0
kspace_style pppm/electrode 1e-5
kspace_modify gewald 0.3
atom_modify map yes
read_data {data}
pair_coeff * *
group bot id 1:{nplate}
group top id {nplate1}:{ntot}
variable qq atom q
compute qbot bot reduce sum v_qq
compute qtop top reduce sum v_qq
thermo 1
thermo_style custom step c_qbot c_qtop
dump 1 all custom 1 {dump} id v_qq
fix conp bot electrode/conp {vm:.3f} {eta} couple top {vp:.3f} symm on
fix nve all nve
run 1
"""


def bilayer(ncell: int = NCELL) -> np.ndarray:
    """Two AB-stacked honeycomb patches (1-based plate split at nplate)."""
    d1 = np.array([A_CC, 0.0])
    d2 = np.array([A_CC / 2, A_CC * np.sqrt(3) / 2])
    basis = [np.array([0.0, 0.0]), np.array([A_CC / 2, A_CC / (2 * np.sqrt(3))])]
    pts = np.array([i * d1 + j * d2 + b for i in range(ncell) for j in range(ncell) for b in basis])
    shift = np.array([A_CC / 2, A_CC / (2 * np.sqrt(3))])
    bot = np.column_stack([pts, np.zeros(len(pts))])
    top = np.column_stack([pts + shift, np.full(len(pts), SEPARATION)])
    both = np.vstack([bot, top])
    off = np.array(
        [(BOX - both[:, 0].max()) / 2, (BOX - both[:, 1].max()) / 2, (BOX - SEPARATION) / 2]
    )
    return both + off


def write_data(path: Path, pos: np.ndarray) -> None:
    with open(path, "w") as f:
        f.write(
            "LAMMPS bilayer\n\n%d atoms\n1 atom types\n\n"
            % len(pos)
            + f"0.0 {BOX} xlo xhi\n0.0 {BOX} ylo yhi\n0.0 {BOX} zlo zhi\n\n"
            + "Masses\n\n1 12.011\n\nAtoms # charge\n\n"
        )
        for k, (x, y, z) in enumerate(pos):
            f.write(f"{k + 1} 1 0.0 {x:.4f} {y:.4f} {z:.4f}\n")


def run_lammps(lmp: str, workdir: Path, data: Path, v: float) -> tuple[float, float, list[float]]:
    inp = INPUT_TEMPLATE.format(
        data=data, nplate=len_charges // 2, nplate1=len_charges // 2 + 1,
        ntot=len_charges, dump=workdir / f"out_{v}.dump", vm=-v, vp=v, eta=ETA,
    )
    (workdir / f"in_{v}").write_text(inp)
    r = subprocess.run(
        ["mpirun", "-np", "1", lmp, "-in", str(workdir / f"in_{v}")],
        capture_output=True, text=True, cwd=workdir,
    )
    if r.returncode != 0:
        raise RuntimeError(f"LAMMPS failed at V={v}:\n{r.stderr[-2000:]}")
    m = re.search(r"^\s*\d+\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$", r.stdout, re.M)
    if not m:
        raise RuntimeError(f"no thermo totals at V={v}")
    qbot, qtop = float(m.group(1)), float(m.group(2))
    charges: dict[int, float] = {}
    for line in open(workdir / f"out_{v}.dump"):
        p = line.split()
        try:
            charges[int(p[0])] = float(p[1])
        except (ValueError, IndexError):
            pass
    q = [charges[i] for i in range(1, len_charges + 1)]
    return qbot, qtop, q


len_charges = 0


def main() -> None:
    global len_charges
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--lmp", default="lmp")
    args = ap.parse_args()

    pos = bilayer()
    len_charges = len(pos)
    nplate = len(pos) // 2

    with tempfile.TemporaryDirectory(prefix="cpm_oracle_") as tmp:
        workdir = Path(tmp)
        data = workdir / "bilayer.data"
        write_data(data, pos)
        volts = []
        for v in VOLTAGES:
            qbot, qtop, q = run_lammps(args.lmp, workdir, data, v)
            print(f"V=±{v}: Qbot={qbot:.5f} Qtop={qtop:.5f}", flush=True)
            volts.append(
                {"voltage": v, "lammps_qbot": qbot, "lammps_qtop": qtop, "lammps_charges": q}
            )

    payload = {
        "reference": "LAMMPS fix electrode/conp (10Dec2025), graphene bilayer 3x3 AB @3.35A, eta=1.979, box=40A ppp",
        "generated": "scripts/validation/oracle_cpm_lammps.py",
        "elements": [6] * len(pos),
        "coords": [c for p in pos for c in (float(p[0]), float(p[1]), float(p[2]))],
        "nplate": nplate,
        "voltages": volts,
    }
    args.out.write_text(json.dumps(payload, indent=1) + "\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
