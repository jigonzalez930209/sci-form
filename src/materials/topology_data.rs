//! RCSR topology nets imported from pormake (generated).
//!
//! Includes nets with multiple node types and nets with inequivalent
//! edges (`xah`, `tfy`), which the homogeneous built-in topologies could
//! not express. Regenerate with `scripts/validation/gen_topologies.py`;
//! do not edit by hand.

/// One node of a topology net.
#[derive(Debug, Clone, Copy)]
pub struct NetNode {
    /// Fractional coordinates.
    pub frac: [f64; 3],
    /// Declared connectivity (number of edges).
    pub connectivity: u32,
}

/// A topology net from the RCSR/pormake database.
#[derive(Debug, Clone, Copy)]
pub struct NetData {
    /// Net name (RCSR symbol or a variant such as `acs-a`).
    pub name: &'static str,
    /// Cell parameters `[a, b, c, alpha, beta, gamma]` (Å, degrees).
    pub cell: [f64; 6],
    /// Nodes with their fractional coordinates and connectivity.
    pub nodes: &'static [NetNode],
    /// Unique edges as `(node_a, node_b)` index pairs into `nodes`.
    pub edges: &'static [(usize, usize)],
    /// Edge-type index per edge (inequivalent edges have distinct types).
    pub edge_types: &'static [u8],
    /// Number of distinct edge types.
    pub n_edge_types: usize,
}

/// Topology nets available to the assembler.
#[rustfmt::skip]
pub static TOPOLOGY_NETS: &[NetData] = &[
    NetData {
        name: "acs",
        cell: [1.41420, 1.41420, 1.15471, 90.00000, 90.00000, 120.00000],
        nodes: &[
            NetNode { frac: [0.66667, 0.33333, 0.25000], connectivity: 6 },
        ],
        edges: &[],
        edge_types: &[
            0,
        ],
        n_edge_types: 1,
    },
    NetData {
        name: "afw",
        cell: [2.88749, 2.88749, 2.23492, 90.00000, 90.00000, 120.00000],
        nodes: &[
            NetNode { frac: [0.00000, 0.19995, 0.00000], connectivity: 4 },
        ],
        edges: &[],
        edge_types: &[
            0,
            1,
        ],
        n_edge_types: 2,
    },
    NetData {
        name: "sod",
        cell: [2.82843, 2.82843, 2.82843, 90.00000, 90.00000, 90.00000],
        nodes: &[
            NetNode { frac: [0.00000, 0.25000, 0.50000], connectivity: 4 },
        ],
        edges: &[],
        edge_types: &[
            0,
        ],
        n_edge_types: 1,
    },
    NetData {
        name: "pts",
        cell: [1.63299, 1.63299, 2.30940, 90.00000, 90.00000, 90.00000],
        nodes: &[
            NetNode { frac: [0.50000, 0.00000, 0.00000], connectivity: 4 },
            NetNode { frac: [0.00000, 0.00000, 0.25000], connectivity: 4 },
        ],
        edges: &[
            (1, 0),
        ],
        edge_types: &[
            0,
        ],
        n_edge_types: 1,
    },
    NetData {
        name: "xah",
        cell: [1.80117, 1.80117, 1.86935, 90.00000, 90.00000, 90.00000],
        nodes: &[
            NetNode { frac: [0.00000, 0.00000, 0.26747], connectivity: 5 },
            NetNode { frac: [0.00000, 0.50000, 0.50000], connectivity: 4 },
        ],
        edges: &[],
        edge_types: &[
            0,
            1,
        ],
        n_edge_types: 2,
    },
    NetData {
        name: "tfy",
        cell: [9.20325, 9.20325, 9.20325, 90.00000, 90.00000, 90.00000],
        nodes: &[
            NetNode { frac: [0.00000, 0.00000, 0.25000], connectivity: 4 },
            NetNode { frac: [0.00000, 0.03237, 0.14627], connectivity: 3 },
            NetNode { frac: [0.05960, 0.05960, 0.05960], connectivity: 3 },
        ],
        edges: &[
            (1, 2),
            (0, 1),
        ],
        edge_types: &[
            0,
            1,
        ],
        n_edge_types: 2,
    },
];
