#!/usr/bin/env python3
"""Generate `tests/fixtures/validation/reaxff_lammps.json` from LAMMPS REAXFF.

Runs the source-built LAMMPS (`pair_style reaxff` + `fix qeq/reaxff` +
`fix reaxff/bonds`, `ffield.reax.cho`) on single-bond stretch series
(CH4 C-H pull, H2O O-H pull) and records total energies (kcal/mol) and the
pulled-bond bond order.

Our side (`forcefield::reaxff` with the transcribed ffield pair table) runs
in-test; the fixture holds the LAMMPS reference + geometries only. Energies
are compared RELATIVELY (E(r) - E(eq)): absolute zeros differ (LAMMPS
includes atomic heat increments). Bond orders compare absolutely.

Usage:
    python scripts/validation/oracle_reaxff_lammps.py [--out FILE] [--lmp LMP]
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "tests/fixtures/validation/reaxff_lammps.json"
# The source-built LAMMPS and its REAXFF force field are not vendored:
# point SCIFORM_REF_DATA at the checkout, or pass --ffield / --lmp.
_REF_DATA = Path(os.environ.get("SCIFORM_REF_DATA", REPO / ".cache" / "refdata"))
DEFAULT_FFIELD = _REF_DATA / "lammps-src" / "potentials" / "ffield.reax.cho"
DEFAULT_LMP = _REF_DATA / "lammps-src" / "build" / "lmp"

INPUT_TEMPLATE = """units real
atom_style charge
boundary p p p
pair_style reaxff NULL
atom_modify map yes
read_data {data}
pair_coeff * * {ffield} {symbols}
fix qeq all qeq/reaxff 1 0.0 10.0 1e-6 reaxff
fix bonds all reaxff/bonds 1 {bonds}
thermo 1
thermo_style custom step pe etotal
run 0
"""

MASSES = {1: 1.0080, 6: 12.0000, 8: 15.9990}


def ch4_series() -> tuple[list[int], list[np.ndarray], int, int]:
    """CH4, pull H[1] along +x. Returns (elements, geometries, i, j)."""
    a = 1.09
    base = np.array(
        [
            [0.0, 0.0, 0.0],
            [a, 0.0, 0.0],
            [-a / 3, a * 0.9428, 0.0],
            [-a / 3, -a * 0.4714, a * 0.8165],
            [-a / 3, -a * 0.4714, -a * 0.8165],
        ]
    )
    elements = [6, 1, 1, 1, 1]
    dists = [0.8, 0.9, 1.0, 1.09, 1.2, 1.4, 1.6, 2.0, 2.5]
    geos = []
    for d in dists:
        g = base.copy()
        g[1] = np.array([d, 0.0, 0.0])
        geos.append((d, g + 10.0))
    return elements, geos, 0, 1


def h2o_series() -> tuple[list[int], list[np.ndarray], int, int]:
    """H2O, pull H[1]. Returns (elements, geometries, i, j)."""
    ang = np.deg2rad(104.5)
    r = 0.96
    base = np.array(
        [
            [0.0, 0.0, 0.0],
            [r, 0.0, 0.0],
            [r * np.cos(ang), r * np.sin(ang), 0.0],
        ]
    )
    elements = [8, 1, 1]
    dists = [0.8, 0.9, 0.96, 1.1, 1.3, 1.6, 2.0]
    geos = []
    for d in dists:
        g = base.copy()
        g[1] = np.array([d, 0.0, 0.0])
        geos.append((d, g + 10.0))
    return elements, geos, 0, 1


def write_data(path: Path, elements: list[int], pos: np.ndarray, types: list[str]) -> None:
    type_id = {s: k + 1 for k, s in enumerate(types)}
    with open(path, "w") as f:
        f.write(f"LAMMPS reaxff\n\n{len(elements)} atoms\n{len(types)} atom types\n\n")
        f.write("0.0 20.0 xlo xhi\n0.0 20.0 ylo yhi\n0.0 20.0 zlo zhi\n\nMasses\n\n")
        for s in types:
            z = {"C": 6, "H": 1, "O": 8}[s]
            f.write(f"{type_id[s]} {MASSES[z]}\n")
        f.write("\nAtoms # charge\n\n")
        z_of = {"C": 6, "H": 1, "O": 8}
        sym_of = {v: k for k, v in z_of.items()}
        for k, (z, (x, y, zc)) in enumerate(zip(elements, pos)):
            f.write(f"{k + 1} {type_id[sym_of[z]]} 0.0 {x:.4f} {y:.4f} {zc:.4f}\n")


def run_one(lmp: str, workdir: Path, elements, pos, types, tag) -> tuple[float, float]:
    data = workdir / f"{tag}.data"
    write_data(data, elements, pos, types)
    inp = INPUT_TEMPLATE.format(
        data=data, ffield=DEFAULT_FFIELD, symbols=" ".join(types), bonds=workdir / f"{tag}.bonds"
    )
    (workdir / f"{tag}.in").write_text(inp)
    r = subprocess.run([lmp, "-in", str(workdir / f"{tag}.in")], capture_output=True, text=True, cwd=workdir)
    if r.returncode != 0:
        raise RuntimeError(f"LAMMPS failed ({tag}):\n{r.stderr[-2000:]}")
    m = re.search(r"^\s*0\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*$", r.stdout, re.M)
    if not m:
        raise RuntimeError(f"no thermo energy ({tag})")
    energy = float(m.group(1))
    # pulled bond order: partners are atoms 1 and 2 (C-H / O-H)
    bo = None
    for line in open(workdir / f"{tag}.bonds"):
        p = line.split()
        if len(p) < 6 or not p[0].isdigit():
            continue
        idx, nb = int(p[0]), int(p[2])
        ids = [int(x) for x in p[3 : 3 + nb]]
        mol = int(p[3 + nb])
        bos = [float(x) for x in p[3 + nb + 1 : 3 + nb + 1 + nb]]
        if idx == 2 and 1 in ids:
            bo = bos[ids.index(1)]
        if idx == 1 and 2 in ids:
            bo = max(bo or 0.0, bos[ids.index(2)])
    if bo is None:
        # Below the coarse bond-order cutoff (0.3): bond considered broken.
        bo = 0.0
    return energy, bo


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--lmp", default=str(DEFAULT_LMP))
    args = ap.parse_args()

    systems = []
    for name, gen, types in [("CH4_CH_pull", ch4_series, ["C", "H"]), ("H2O_OH_pull", h2o_series, ["O", "H"])]:
        elements, geos, i, j = gen()
        with tempfile.TemporaryDirectory(prefix="reaxff_oracle_") as tmp:
            workdir = Path(tmp)
            points = []
            for d, g in geos:
                e, bo = run_one(args.lmp, workdir, elements, g, types, f"{name}_{d}")
                print(f"{name} r={d}: E={e:.3f} BO={bo:.4f}", flush=True)
                points.append(
                    {
                        "distance": d,
                        "coords": [float(c) for row in g for c in row],
                        "lammps_energy_kcal": e,
                        "lammps_bo": bo,
                    }
                )
        systems.append({"name": name, "elements": elements, "bond": [i, j], "points": points})

    payload = {
        "reference": "LAMMPS REAXFF (29Aug2024-U2 source build) ffield.reax.cho + qeq/reaxff; energies relative in gate",
        "generated": "scripts/validation/oracle_reaxff_lammps.py",
        "systems": systems,
    }
    args.out.write_text(json.dumps(payload, indent=1) + "\n")
    print("wrote", args.out)


if __name__ == "__main__":
    main()
