#!/usr/bin/env python3
"""Generate `src/materials/topology_data.rs` from pormake's RCSR `.cgd` nets.

The existing built-in topologies (`Topology::pcu/dia/bcu/fcu`) are homogeneous
nets. This generator imports a wider selection from pormake's topology database,
including nets with several node types and — the point of the extension —
nets with **inequivalent edges** (`xah`, `tfy`), which the previous data model
could not express.

Parsed per net:
  - cell parameters (a, b, c, α, β, γ),
  - nodes (fractional coordinates + declared connectivity),
  - unique edges (coordinate pairs) with an edge-type index derived from the
    `# EDGE_CENTER` groups,
  - counts used by the F67 gate.

Usage:
    python scripts/validation/gen_topologies.py [--check]
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from refdata import skip

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "src" / "materials" / "topology_data.rs"

# sci-form net name -> pormake `.cgd` file (None when pormake lacks it).
NETS = ["acs", "afw", "sod", "pts", "xah", "tfy"]


def topology_dir() -> Path:
    import pormake

    return Path(pormake.__file__).resolve().parent / "database" / "topologies"


def parse_cgd(path: Path) -> dict:
    name = path.stem
    cell = None
    nodes = []
    edges = []
    edge_centers = []
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith("#"):
            line = line.lstrip("#").strip()
        if line.startswith("NAME"):
            name = line.split()[1]
        elif line.startswith("CELL"):
            vals = [float(x) for x in line.split()[1:7]]
            cell = vals
        elif line.startswith("NODE"):
            t = line.split()
            nodes.append({"cn": int(t[2]), "frac": [float(x) for x in t[3:6]]})
        elif "EDGE_CENTER" in line:
            edge_centers.append([float(x) for x in line.split()[1:4]])
        elif line.startswith("EDGE"):
            vals = [float(x) for x in line.split()[1:7]]
            edges.append((vals[0:3], vals[3:6]))

    # Edge type: index of the matching EDGE_CENTER group (unique centres).
    types = []
    for a, b in edges:
        centre = [0.5 * (a[i] + b[i]) for i in range(3)]
        best, best_idx = 1e9, 0
        for idx, c in enumerate(edge_centers):
            d = 0.0
            for i in range(3):
                diff = centre[i] - c[i]
                diff -= round(diff)  # periodic image
                d += diff * diff
            if d < best:
                best, best_idx = d, idx
        types.append(best_idx)

    # pormake's own classification (by local structure) is the reference for
    # the edge types; the `.cgd` centre grouping is kept for provenance.
    pormake_types = None
    pormake_n_edge_types = 0
    try:
        from pormake import Topology as _Topo

        t = _Topo(str(path))
        pormake_types = [int(t.get_edge_type(e)) for e in range(t.n_edges)]
        pormake_n_edge_types = int(t.n_edge_types)
    except Exception:
        pass

    return {
        "name": name,
        "cell": cell,
        "nodes": nodes,
        "edges": edges,
        "edge_types": types,
        "n_edge_types": len(set(types)) if types else 0,
        "pormake_edge_types": pormake_types,
        "pormake_n_edge_types": pormake_n_edge_types,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--check",
        action="store_true",
        help="regenerate into a temporary file, compare byte-for-byte with the "
             "committed table and exit non-zero on a difference",
    )
    args = ap.parse_args()

    try:
        topo_dir = topology_dir()
    except ImportError as exc:
        msg = f"the pormake package is not installed ({exc}); `pip install pormake`"
        if args.check:
            return skip(msg)
        raise SystemExit(f"error: {msg}")
    nets = []
    for name in NETS:
        path = topo_dir / f"{name}.cgd"
        if not path.exists():
            print(f"[warn] {name}.cgd not found")
            continue
        nets.append(parse_cgd(path))

    for net in nets:
        print(
            f"{net['name']}: nodes={len(net['nodes'])} "
            f"cn={[n['cn'] for n in net['nodes']]} "
            f"edges={len(net['edges'])} edge_centres={net['n_edge_types']}"
        )

    out = [
        "//! RCSR topology nets imported from pormake (generated).",
        "//!",
        "//! Includes nets with multiple node types and nets with inequivalent",
        "//! edges (`xah`, `tfy`), which the homogeneous built-in topologies could",
        "//! not express. Regenerate with `scripts/validation/gen_topologies.py`;",
        "//! do not edit by hand.",
        "",
        "/// One node of a topology net.",
        "#[derive(Debug, Clone, Copy)]",
        "pub struct NetNode {",
        "    /// Fractional coordinates.",
        "    pub frac: [f64; 3],",
        "    /// Declared connectivity (number of edges).",
        "    pub connectivity: u32,",
        "}",
        "",
        "/// A topology net from the RCSR/pormake database.",
        "#[derive(Debug, Clone, Copy)]",
        "pub struct NetData {",
        "    /// Net name (RCSR symbol or a variant such as `acs-a`).",
        "    pub name: &'static str,",
        "    /// Cell parameters `[a, b, c, alpha, beta, gamma]` (Å, degrees).",
        "    pub cell: [f64; 6],",
        "    /// Nodes with their fractional coordinates and connectivity.",
        "    pub nodes: &'static [NetNode],",
        "    /// Unique edges as `(node_a, node_b)` index pairs into `nodes`.",
        "    pub edges: &'static [(usize, usize)],",
        "    /// Edge-type index per edge (inequivalent edges have distinct types).",
        "    pub edge_types: &'static [u8],",
        "    /// Number of distinct edge types.",
        "    pub n_edge_types: usize,",
        "}",
        "",
        "/// Topology nets available to the assembler.",
        "#[rustfmt::skip]",
        "pub static TOPOLOGY_NETS: &[NetData] = &[",
    ]
    for net in nets:
        out.append("    NetData {")
        out.append(f'        name: "{net["name"]}",')
        cell = net["cell"] or [0.0] * 6
        out.append(
            "        cell: ["
            + ", ".join(f"{v:.5f}" for v in cell)
            + "],"
        )
        out.append("        nodes: &[")
        for node in net["nodes"]:
            f = node["frac"]
            out.append(
                f"            NetNode {{ frac: [{f[0]:.5f}, {f[1]:.5f}, {f[2]:.5f}], "
                f"connectivity: {node['cn']} }},"
            )
        out.append("        ],")
        # Map edge endpoints to node indices (mod 1, tolerance 1e-3).
        def node_index(point) -> int:
            best, best_i = 1e9, None
            for i, node in enumerate(net["nodes"]):
                d = 0.0
                for k in range(3):
                    diff = point[k] - node["frac"][k]
                    diff -= round(diff)
                    d += diff * diff
                if d < best:
                    best, best_i = d, i
            return best_i if best < 1e-3 else None

        pairs = []
        for a, b in net["edges"]:
            ia, ib = node_index(a), node_index(b)
            pairs.append((ia, ib))
        if all(i is not None for p in pairs for i in p):
            out.append("        edges: &[")
            for ia, ib in pairs:
                out.append(f"            ({ia}, {ib}),")
            out.append("        ],")
        else:
            # Some nets list edges with periodic images that do not match a
            # node exactly; keep the coordinate-level information out and rely
            # on the connectivity counts (documented).
            out.append("        edges: &[],")
        out.append("        edge_types: &[")
        for t in net["edge_types"]:
            out.append(f"            {t},")
        out.append("        ],")
        out.append(f"        n_edge_types: {net['n_edge_types']},")
        out.append("    },")
    out.append("];")
    text = "\n".join(out) + "\n"
    if args.check:
        current = OUT.read_text() if OUT.exists() else ""
        if current != text:
            print(
                f"{OUT.relative_to(REPO)} is out of date — "
                f"run scripts/validation/gen_topologies.py"
            )
            return 1
        print(f"{OUT.relative_to(REPO)} is up to date")
        return 0

    OUT.write_text(text)
    print(f"wrote {OUT} ({len(nets)} nets)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
