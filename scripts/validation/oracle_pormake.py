#!/usr/bin/env python3
"""pormake + pyzeo (Zeo++ engine) oracle for sci-form validation (F67/F68).

1. Topology nets: parse pormake's RCSR `.cgd` files (NODE coordination
   numbers) for every net with a sci-form counterpart.
2. MOF-5 end-to-end: pormake Builder on (pcu, N16 = Zn4O node,
   E14 = phenylene edge) -> primitive cell, composition, bond count,
   CIF, and the building-block geometries (X connection points removed,
   canonical frames) that the Rust test feeds to its own assembler.
3. pyzeo reference porosities on the pormake MOF-5 CIF with sci-form's
   own vdW radii (same radii/probe on both sides, per the plan):
   accessible-volume fraction (probes 0 and 1.2 A) and the
   largest-included / free-sphere diameters (Di/Df).

Generates `tests/fixtures/validation/pormake_mof.json`.

Usage:
    python scripts/validation/oracle_pormake.py
"""

from __future__ import annotations

import re
import tempfile
import time
from collections import Counter
from pathlib import Path

import numpy as np

from common import Report, emit_fixture_report

TOPO_DIR = (
    Path(__import__("pormake").__file__).resolve().parent
    / "database"
    / "topologies"
)
BB_DIR = (
    Path(__import__("pormake").__file__).resolve().parent / "database" / "bbs"
)

# sci-form net -> pormake .cgd file (None = no counterpart in pormake's DB).
NET_MAP = {
    "pcu": "pcu",
    "dia": "dia",
    "sql": None,  # no sql.cgd ships with pormake
    "bcu": "bcu",
    "fcu": "fcu",
    "nbo": "nbo",
    "pts": "pts",
    "kgm": "kag",  # RCSR kagome; NOTE below
    "hxl": None,  # RCSR hex.cgd is the 8-c net, not hxl (6-c)
    "pillared_sql": None,  # derived net, no RCSR entry
    # Extended set: multiple node types and nets with inequivalent edges
    # (validation F67 extension).
    "acs": "acs",
    "afw": "afw",
    "sod": "sod",
    "pts": "pts",
    "xah": "xah",
    "tfy": "tfy",
}

# sci-form's vdW table (src/materials/cell.rs) for the pyzeo .rad file.
VDW = {1: 1.20, 6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.80, 16: 1.80,
       17: 1.75, 22: 2.00, 26: 2.00, 29: 1.40, 30: 1.39, 35: 1.85, 53: 1.98}
SYM2Z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
         "Cl": 17, "Ti": 22, "Fe": 26, "Cu": 29, "Zn": 30, "Br": 35, "I": 53}

# MOF-5 conventional cubic cell, Fm-3m (Rowsell et al., JACS 2004).
MOF5_A_CONV_LIT = 25.832


def parse_cgd(path: Path) -> dict:
    nodes, edges, centres = [], [], []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("#"):
            line = line.lstrip("#").strip()
        if line.startswith("NODE"):
            t = line.split()
            nodes.append({"cn": int(t[2]), "frac": [float(x) for x in t[3:6]]})
        elif line.startswith("EDGE_CENTER"):
            centres.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("EDGE"):
            t = [float(x) for x in line.split()[1:]]
            edges.append(t)
    return {
        "nodes": nodes,
        "n_edges": len(edges),
        "n_edge_centres": len(centres),
        "n_edge_types": len(centres),
    }


def read_bb(name: str):
    """(elements, coords) with X connection points removed."""
    lines = (BB_DIR / f"{name}.xyz").read_text().splitlines()
    n = int(lines[0])
    els, xyz, xpts = [], [], []
    for line in lines[2: 2 + n]:
        parts = line.split()
        if parts[0] == "X":
            xpts.append([float(x) for x in parts[1:4]])
        else:
            els.append(SYM2Z[parts[0]])
            xyz.append([float(x) for x in parts[1:4]])
    return els, np.array(xyz), np.array(xpts)


def canonical_node() -> dict:
    els, xyz, _ = read_bb("N16")
    xyz = xyz - xyz.mean(axis=0)
    return {"elements": els, "positions": xyz.tolist()}


def canonical_edge() -> dict:
    els, xyz, xpts = read_bb("E14")
    # Align the X-X axis with +x, center at origin (matches the assembler's
    # "linker along x-axis" convention).
    axis = xpts[1] - xpts[0]
    axis = axis / np.linalg.norm(axis)
    target = np.array([1.0, 0.0, 0.0])
    v = np.cross(axis, target)
    c = float(np.dot(axis, target))
    if np.linalg.norm(v) < 1e-12:
        rot = np.eye(3) if c > 0 else -np.eye(3)
    else:
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rot = np.eye(3) + vx + vx @ vx / (1 + c)
    xyz = (rot @ (xyz - xyz.mean(axis=0)).T).T
    return {"elements": els, "positions": xyz.tolist()}


def build_mof5() -> dict:
    from pormake.builder import Builder
    from pormake.database import Database

    db = Database()
    topo = db.get_topo("pcu")
    fw = Builder().build_by_type(
        topo, {0: db.get_bb("N16")}, {(0, 0): db.get_bb("E14")}
    )
    atoms = fw.atoms
    cell = [float(x) for x in atoms.cell.lengths()]
    symbols = atoms.get_chemical_symbols()
    comp = dict(Counter(symbols))
    tmp = Path(tempfile.gettempdir()) / "pormake_mof5.cif"
    fw.write_cif(str(tmp))
    return {
        "cell_lengths": cell,
        "cell_angles": [float(x) for x in atoms.cell.angles()],
        "n_atoms": len(atoms),
        "composition": comp,
        "n_bonds": int(len(fw.bonds)),
        "cif": tmp.read_text(),
    }


def zeo_reference(cif_text: str) -> dict:
    import pyzeo

    tmp = Path(tempfile.gettempdir()) / "pormake_mof5.cif"
    tmp.write_text(cif_text)
    rad = Path(tempfile.gettempdir()) / "sc_radii.rad"
    rad.write_text(
        "".join(
            f"{sym} {r}\n"
            for sym, r in sorted(
                [(s, v) for s, v in
                 [("H", 1.20), ("C", 1.70), ("N", 1.55), ("O", 1.52),
                  ("Zn", 1.39)]]
            )
        )
    )
    atmnet = pyzeo.AtomNetwork.read_from_CIF(str(tmp), rad_file=str(rad))
    out = {}
    for probe, mc in ((0.0, 200000), (1.2, 200000)):
        raw = pyzeo.volume(atmnet, probe, probe, mc)
        text = raw.decode() if isinstance(raw, bytes) else str(raw)
        m = re.search(r"AV_Volume_fraction:\s*([0-9.]+)", text)
        out[f"vf_probe_{probe}"] = float(m.group(1))
    # High-accuracy network: MOF-5's exact symmetries create degenerate
    # Voronoi vertices that crash raw voro++ ("Unable to add edge").
    pyzeo.high_accuracy_atomnet(atmnet)
    atmnet.perform_voronoi_decomposition()
    dia = Path(tempfile.gettempdir()) / "pormake_mof5.dia"
    atmnet.calculate_free_sphere_parameters(str(dia))
    nums = [float(x) for x in dia.read_text().split()[-3:]]
    out["lcd_diameter"] = nums[0]
    out["pld_diameter"] = nums[1]
    return out


def main() -> int:
    topologies = {}
    for ours, cgd in NET_MAP.items():
        if cgd is None:
            topologies[ours] = {"cgd": None}
            continue
        topologies[ours] = {"cgd": cgd, **parse_cgd(TOPO_DIR / f"{cgd}.cgd")}
        info = topologies[ours]
        print(
            f"  net {ours:13} <- {cgd:5}: nodes={len(info['nodes'])} "
            f"cn={[n['cn'] for n in info['nodes']]}",
            flush=True,
        )

    mof5 = build_mof5()
    print(
        f"  MOF-5 pcu+N16+E14: a={mof5['cell_lengths'][0]:.4f} "
        f"n={mof5['n_atoms']} comp={mof5['composition']} bonds={mof5['n_bonds']}",
        flush=True,
    )
    a_conv = 2 * mof5["cell_lengths"][0]
    print(
        f"  conventional 2a={a_conv:.4f} vs lit {MOF5_A_CONV_LIT} "
        f"({(a_conv / MOF5_A_CONV_LIT - 1) * 100:+.2f}%)",
        flush=True,
    )
    zeo = zeo_reference(mof5["cif"])
    print(f"  pyzeo: {zeo}", flush=True)

    payload = {
        "reference": "pormake 0.2.3 (RCSR nets + Builder) + pyzeo 0.2.0 (Zeo++ engine)",
        "topologies": topologies,
        "mof5": {**mof5, "a_conv_lit": MOF5_A_CONV_LIT,
                 "a_conv_dev": a_conv / MOF5_A_CONV_LIT - 1},
        "node_bb": canonical_node(),
        "edge_bb": canonical_edge(),
        "pyzeo": zeo,
        "radii": "sci-form vdw table (H1.20 C1.70 N1.55 O1.52 Zn1.39)",
        "generated": time.strftime("%Y-%m-%d"),
    }
    report = Report(
        ficha="F67",
        observable="mof5_assembly_vs_pormake",
        reference_tool="pormake Builder (pcu/N16/E14)",
        units="atoms/equation",
        tolerance=0.0,
        level="L2",
        n_cases=1,
        metric="max_abs",
        value=float("nan"),
        notes=(
            "Same BB geometries + same cell: atom count, composition and "
            "porosity must agree. Topology degrees vs RCSR .cgd "
            "coordinations gated in tests/cross_validation.rs::f67_topology. "
            "Conventional-cell literature anchor (25.832 A) reported, not "
            "gated (assembler takes the cell as input; pormake scaling is "
            "+2.55% off experiment)."
        ),
    )
    emit_fixture_report("pormake_mof.json", payload, report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
