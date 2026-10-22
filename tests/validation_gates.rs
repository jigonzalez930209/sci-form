//! Cross-library validation gates for fichas F64, F65, F66 and F69.
//!
//! These four gates were documented as PASS in
//! `docs/validation/cross-library-validation.es.md` (§4.8 and the results
//! table) but no test consumed the corresponding reference fixtures. This
//! file implements them against the committed fixtures — no RDKit /
//! pymatgen / ASE installation is required, only the JSON references:
//!
//! | Gate | Fixture | Reference | Criterion |
//! |---|---|---|---|
//! | F64 | `space_groups_230.json` | pymatgen/spglib ITC tables | 230/230 groups: order, symbol, system, lattice, point group, operation set, equivalent positions (±1e-8 fractional) |
//! | F65 | `pymatgen_cif.json` | pymatgen `CifWriter`/`CifParser` | cell parameters ≤ 1e-6 Å/deg, site count, elements, stable round-trip |
//! | F66 | `pymatgen_pbc.json` | ASE `find_mic` | minimum-image distances ≤ 1e-10 Å + wrap/lattice-shift invariants |
//! | F69 | `pymatgen_xrd.json` | pymatgen `XRDCalculator` (Cu-Kα) | 2θ ≤ 0.01°, intensity correlation ≥ 0.95, no extra/missing peaks |
//!
//! The tolerances are the ones recorded in the fixture provenance rows of
//! `reports/validation/summary.json` (F65: 1e-6, F66: 1e-10, F69: 0.01°)
//! and in the validation plan (§4.8).

use std::path::PathBuf;

use serde::Deserialize;

use sci_form::materials::cell::CellParameters;
use sci_form::materials::cif::{parse_cif, write_cif};
use sci_form::materials::space_groups::{
    crystal_system_for_number, space_group_by_number, space_group_by_symbol, CrystalSystem,
    LatticeType,
};
use sci_form::materials::{simulate_powder_xrd, UnitCell};

// ═══════════════════════════════════════════════════════════════════════════
// Fixture plumbing (mirrors tests/cross_validation.rs)
// ═══════════════════════════════════════════════════════════════════════════

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("validation")
        .join(name)
}

/// Load a JSON fixture (`read_text_fixture` also resolves a `.gz` sibling).
fn load<T: for<'de> Deserialize<'de>>(name: &str) -> T {
    let path = fixture_path(name);
    let data = sci_form::fixture_io::read_text_fixture(&path).unwrap_or_else(|err| {
        panic!(
            "missing fixture {} ({err}); run `scripts/validation/oracle_pymatgen.py` \
             (or the matching oracle) to regenerate it",
            path.display()
        )
    });
    serde_json::from_str(&data).expect("invalid fixture JSON")
}

/// Write the measured metric to `target/validation/{ficha}.json` so
/// `scripts/validation/collect.py` can merge it into
/// `reports/validation/summary.json` after the test run.
fn record(
    ficha: &str,
    observable: &str,
    metric: &str,
    value: f64,
    tolerance: f64,
    passed: bool,
    n_cases: usize,
) {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target")
        .join("validation");
    let _ = std::fs::create_dir_all(&dir);
    let payload = serde_json::json!({
        "ficha": ficha,
        "observable": observable,
        "metric": metric,
        "value": value,
        "tolerance": tolerance,
        "pass": passed,
        "n_cases": n_cases,
        "source": "tests/validation_gates.rs",
    });
    let path = dir.join(format!("{ficha}.json"));
    let mut rows: Vec<serde_json::Value> = std::fs::read_to_string(&path)
        .ok()
        .and_then(|text| serde_json::from_str::<serde_json::Value>(&text).ok())
        .map(|value| match value {
            serde_json::Value::Array(items) => items,
            other => vec![other],
        })
        .unwrap_or_default();
    rows.retain(|row| row.get("observable").and_then(|v| v.as_str()) != Some(observable));
    rows.push(payload);
    let _ = std::fs::write(&path, serde_json::Value::Array(rows).to_string());
}

/// Pearson correlation (constant series count as perfectly correlated).
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    let ma = a.iter().sum::<f64>() / n as f64;
    let mb = b.iter().sum::<f64>() / n as f64;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for i in 0..n {
        let x = a[i] - ma;
        let y = b[i] - mb;
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da == 0.0 || db == 0.0 {
        return 1.0;
    }
    num / (da.sqrt() * db.sqrt())
}

/// Element symbols for Z = 1..=54 (used to map the reference element
/// labels of the fixtures onto atomic numbers).
const ELEMENT_SYMBOLS: [&str; 54] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe",
];

fn atomic_number(symbol: &str) -> u8 {
    ELEMENT_SYMBOLS
        .iter()
        .position(|s| *s == symbol)
        .map(|idx| (idx + 1) as u8)
        .unwrap_or_else(|| panic!("unknown element symbol {symbol:?}"))
}

// ═══════════════════════════════════════════════════════════════════════════
// F64 — 230 space groups vs pymatgen/spglib ITC tables
// ═══════════════════════════════════════════════════════════════════════════

#[derive(Deserialize)]
struct SpaceGroupFixture {
    reference: String,
    groups: Vec<SpaceGroupCase>,
}

#[derive(Deserialize)]
struct SpaceGroupCase {
    number: u16,
    symbol: String,
    crystal_system: String,
    lattice: String,
    point_group: String,
    operations: Vec<ReferenceOp>,
}

#[derive(Deserialize)]
struct ReferenceOp {
    rotation: [[i32; 3]; 3],
    translation: [f64; 3],
}

fn crystal_system_from_name(name: &str) -> Option<CrystalSystem> {
    match name {
        "Triclinic" => Some(CrystalSystem::Triclinic),
        "Monoclinic" => Some(CrystalSystem::Monoclinic),
        "Orthorhombic" => Some(CrystalSystem::Orthorhombic),
        "Tetragonal" => Some(CrystalSystem::Tetragonal),
        "Trigonal" => Some(CrystalSystem::Trigonal),
        "Hexagonal" => Some(CrystalSystem::Hexagonal),
        "Cubic" => Some(CrystalSystem::Cubic),
        _ => None,
    }
}

fn lattice_type_from_name(name: &str) -> Option<LatticeType> {
    match name {
        "P" => Some(LatticeType::P),
        "I" => Some(LatticeType::I),
        "F" => Some(LatticeType::F),
        "C" => Some(LatticeType::C),
        "R" => Some(LatticeType::R),
        "A" => Some(LatticeType::A),
        "B" => Some(LatticeType::B),
        _ => None,
    }
}

/// Apply a reference symmetry operation to fractional coordinates and wrap
/// the result into [0, 1).
fn apply_reference_op(op: &ReferenceOp, frac: [f64; 3]) -> [f64; 3] {
    let mut out = [0.0f64; 3];
    for (i, out_i) in out.iter_mut().enumerate() {
        let v = op.rotation[i][0] as f64 * frac[0]
            + op.rotation[i][1] as f64 * frac[1]
            + op.rotation[i][2] as f64 * frac[2]
            + op.translation[i];
        *out_i = v - v.floor();
    }
    out
}

/// Periodic (minimum-image) distance between two fractional positions.
fn frac_periodic_distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3)
        .map(|k| {
            let d = (a[k] - b[k]).abs();
            d.min(1.0 - d)
        })
        .fold(0.0f64, f64::max)
}

/// F64 — the 230 ITC space groups must match the reference tables exactly:
/// Hermann-Mauguin symbol, crystal system, lattice centering, point group,
/// group order, the full set of symmetry operations and the equivalent
/// positions they generate (validation plan: "0 discrepancias", ±1e-8
/// fractional for the positions).
#[test]
fn f64_space_groups_match_spglib() {
    let fixture: SpaceGroupFixture = load("space_groups_230.json");
    assert_eq!(
        fixture.groups.len(),
        230,
        "the ITC fixture must carry all 230 space groups (reference: {})",
        fixture.reference
    );

    // A generic probe site: a general position for every group, so the
    // orbit size equals the group order.
    let probe = [0.1234f64, 0.2345, 0.3456];

    let mut max_translation_dev = 0.0f64;
    let mut max_position_dev = 0.0f64;
    let mut n_ops = 0usize;

    for (idx, case) in fixture.groups.iter().enumerate() {
        assert_eq!(
            case.number as usize,
            idx + 1,
            "the fixture groups must be ordered by ITC number"
        );

        let sg = space_group_by_number(case.number)
            .unwrap_or_else(|| panic!("space group #{} is missing from the table", case.number));

        // ── identity data ──
        assert_eq!(sg.number, case.number);
        assert_eq!(
            sg.symbol, case.symbol,
            "Hermann-Mauguin symbol mismatch for group #{}",
            case.number
        );
        let reference_system = crystal_system_from_name(&case.crystal_system)
            .unwrap_or_else(|| panic!("unknown crystal system {:?}", case.crystal_system));
        assert_eq!(
            sg.crystal_system, reference_system,
            "crystal system mismatch for group #{}",
            case.number
        );
        assert_eq!(
            crystal_system_for_number(case.number),
            Some(reference_system),
            "crystal_system_for_number disagrees for group #{}",
            case.number
        );
        let reference_lattice = lattice_type_from_name(&case.lattice)
            .unwrap_or_else(|| panic!("unknown lattice type {:?}", case.lattice));
        assert_eq!(
            sg.lattice_type, reference_lattice,
            "lattice centering mismatch for group #{}",
            case.number
        );
        assert_eq!(
            sg.point_group, case.point_group,
            "point group mismatch for group #{}",
            case.number
        );

        // ── group order (multiplicity) ──
        assert_eq!(
            sg.multiplicity(),
            case.operations.len(),
            "group order mismatch for group #{}",
            case.number
        );
        n_ops += case.operations.len();

        // ── symbol lookup round-trip ──
        let by_symbol = space_group_by_symbol(&case.symbol).unwrap_or_else(|| {
            panic!(
                "symbol {:?} of group #{} is not found by space_group_by_symbol",
                case.symbol, case.number
            )
        });
        assert_eq!(
            by_symbol.number, case.number,
            "space_group_by_symbol({:?}) resolved to group #{}",
            case.symbol, by_symbol.number
        );

        // ── symmetry operation set ──
        let mut used = vec![false; sg.operations.len()];
        for op in &case.operations {
            let mut best: Option<(usize, f64)> = None;
            for (op_idx, rust_op) in sg.operations.iter().enumerate() {
                if used[op_idx] {
                    continue;
                }
                let same_rotation = (0..3)
                    .all(|i| (0..3).all(|j| rust_op.rotation[i][j] as i32 == op.rotation[i][j]));
                if !same_rotation {
                    continue;
                }
                let dev = (0..3)
                    .map(|k| (rust_op.translation[k] - op.translation[k]).abs())
                    .fold(0.0f64, f64::max);
                match best {
                    Some((_, best_dev)) if best_dev <= dev => {}
                    _ => best = Some((op_idx, dev)),
                }
            }
            let (op_idx, dev) = best.unwrap_or_else(|| {
                panic!(
                    "group #{}: no Rust operation matches the reference rotation {:?}",
                    case.number, op.rotation
                )
            });
            used[op_idx] = true;
            max_translation_dev = max_translation_dev.max(dev);
            assert!(
                dev <= 1e-12,
                "group #{}: translation mismatch {:.3e} (ours {:?}, ref {:?})",
                case.number,
                dev,
                sg.operations[op_idx].translation,
                op.translation
            );
        }

        // ── equivalent positions ──
        let ours = sg.generate_equivalent_positions(&probe);
        let reference: Vec<[f64; 3]> = case
            .operations
            .iter()
            .map(|op| apply_reference_op(op, probe))
            .collect();
        assert_eq!(
            ours.len(),
            reference.len(),
            "group #{}: orbit size {} != reference {}",
            case.number,
            ours.len(),
            reference.len()
        );
        for p in &ours {
            let dev = reference
                .iter()
                .map(|q| frac_periodic_distance(*p, *q))
                .fold(f64::INFINITY, f64::min);
            max_position_dev = max_position_dev.max(dev);
            assert!(
                dev <= 1e-8,
                "group #{}: equivalent position {p:?} is {dev:.3e} away from the reference orbit",
                case.number
            );
        }
    }

    println!(
        "F64: 230/230 groups match ({} operations, max translation dev {:.3e}, max position dev {:.3e})",
        n_ops, max_translation_dev, max_position_dev
    );
    record(
        "F64",
        "space_groups_vs_spglib",
        "max_abs_translation",
        max_translation_dev,
        1e-12,
        max_translation_dev <= 1e-12 && max_position_dev <= 1e-8,
        fixture.groups.len(),
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// F65 — CIF parsing vs pymatgen
// ═══════════════════════════════════════════════════════════════════════════

#[derive(Deserialize)]
struct CifFixture {
    reference: String,
    cases: Vec<CifCase>,
}

#[derive(Deserialize)]
struct CifCase {
    name: String,
    cif: String,
    lattice_params: [f64; 6],
    n_sites: usize,
    elements: Vec<String>,
}

/// Value of a single-value CIF tag, read independently of the Rust parser.
fn cif_tag_value<'a>(cif: &'a str, tag: &str) -> Option<&'a str> {
    for line in cif.lines() {
        let line = line.trim();
        if let Some(rest) = line.strip_prefix(tag) {
            let rest = rest.trim();
            if rest.is_empty() {
                continue;
            }
            return Some(rest.trim_matches('\'').trim_matches('"').trim());
        }
    }
    None
}

/// F65 — CIFs written by pymatgen must parse to the same cell, the same
/// number of sites and the same elements, and survive a parse → write →
/// parse round-trip (validation plan: "parámetros ≤ 1e-6; round-trip
/// estable").
#[test]
fn f65_cif_matches_pymatgen() {
    let fixture: CifFixture = load("pymatgen_cif.json");
    assert!(!fixture.cases.is_empty());

    let mut max_param_dev = 0.0f64;
    let mut max_roundtrip_dev = 0.0f64;

    for case in &fixture.cases {
        let parsed = parse_cif(&case.cif)
            .unwrap_or_else(|err| panic!("{}: parse_cif failed: {err}", case.name));

        // ── cell parameters ──
        let ours = [
            parsed.cell_params.a,
            parsed.cell_params.b,
            parsed.cell_params.c,
            parsed.cell_params.alpha,
            parsed.cell_params.beta,
            parsed.cell_params.gamma,
        ];
        for (k, (o, r)) in ours.iter().zip(case.lattice_params.iter()).enumerate() {
            let dev = (o - r).abs();
            max_param_dev = max_param_dev.max(dev);
            assert!(
                dev <= 1e-6,
                "{}: cell parameter #{k} deviates by {dev:.3e} (ours {o}, ref {r})",
                case.name
            );
        }

        // The lattice vectors built from those parameters must reproduce
        // them (cell ⇄ parameters consistency).
        let rebuilt = parsed.cell.parameters();
        for (o, r) in [
            (rebuilt.a, parsed.cell_params.a),
            (rebuilt.b, parsed.cell_params.b),
            (rebuilt.c, parsed.cell_params.c),
            (rebuilt.alpha, parsed.cell_params.alpha),
            (rebuilt.beta, parsed.cell_params.beta),
            (rebuilt.gamma, parsed.cell_params.gamma),
        ] {
            let dev = (o - r).abs();
            assert!(
                dev <= 1e-6,
                "{}: cell vectors do not reproduce the cell parameters ({dev:.3e})",
                case.name
            );
        }

        // ── sites ──
        assert_eq!(
            parsed.atom_sites.len(),
            case.n_sites,
            "{}: parsed {} sites, reference has {}",
            case.name,
            parsed.atom_sites.len(),
            case.n_sites
        );
        let mut ours_elements: Vec<String> = parsed
            .atom_sites
            .iter()
            .map(|site| site.type_symbol.clone())
            .collect();
        ours_elements.sort();
        ours_elements.dedup();
        assert_eq!(
            ours_elements, case.elements,
            "{}: element set mismatch",
            case.name
        );
        for site in &parsed.atom_sites {
            assert!(
                (site.occupancy - 1.0).abs() <= 1e-12,
                "{}: unexpected occupancy {} for site {}",
                case.name,
                site.occupancy,
                site.label
            );
            assert!(
                site.frac_x >= 0.0
                    && site.frac_x < 1.0
                    && site.frac_y >= 0.0
                    && site.frac_y < 1.0
                    && site.frac_z >= 0.0
                    && site.frac_z < 1.0,
                "{}: fractional coordinates not wrapped for site {}",
                case.name,
                site.label
            );
        }

        // ── space-group tags ──
        let ref_hm = cif_tag_value(&case.cif, "_symmetry_space_group_name_H-M")
            .unwrap_or_else(|| panic!("{}: the reference CIF has no H-M tag", case.name));
        assert_eq!(
            parsed.space_group_hm.as_deref(),
            Some(ref_hm),
            "{}: space-group symbol mismatch",
            case.name
        );
        let ref_number = cif_tag_value(&case.cif, "_symmetry_Int_Tables_number")
            .and_then(|v| v.parse::<u16>().ok())
            .unwrap_or_else(|| panic!("{}: the reference CIF has no ITC number", case.name));
        assert_eq!(
            parsed.space_group_number,
            Some(ref_number),
            "{}: space-group number mismatch",
            case.name
        );

        // ── round-trip: parse → write → parse must be stable ──
        // Tolerances follow the writer's fixed precision (6 decimals for
        // lengths and fractional coordinates, 4 decimals for the angles).
        let reparsed = parse_cif(&write_cif(&parsed))
            .unwrap_or_else(|err| panic!("{}: round-trip parse failed: {err}", case.name));
        for (o, r) in [
            (reparsed.cell_params.a, parsed.cell_params.a),
            (reparsed.cell_params.b, parsed.cell_params.b),
            (reparsed.cell_params.c, parsed.cell_params.c),
        ] {
            let dev = (o - r).abs();
            max_roundtrip_dev = max_roundtrip_dev.max(dev);
            assert!(
                dev <= 1e-6,
                "{}: round-trip cell length drifted by {dev:.3e}",
                case.name
            );
        }
        for (o, r) in [
            (reparsed.cell_params.alpha, parsed.cell_params.alpha),
            (reparsed.cell_params.beta, parsed.cell_params.beta),
            (reparsed.cell_params.gamma, parsed.cell_params.gamma),
        ] {
            let dev = (o - r).abs();
            max_roundtrip_dev = max_roundtrip_dev.max(dev);
            assert!(
                dev <= 1e-4,
                "{}: round-trip cell angle drifted by {dev:.3e}",
                case.name
            );
        }
        assert_eq!(
            reparsed.atom_sites.len(),
            parsed.atom_sites.len(),
            "{}: round-trip changed the site count",
            case.name
        );
        for (rt, orig) in reparsed.atom_sites.iter().zip(parsed.atom_sites.iter()) {
            assert_eq!(
                rt.type_symbol, orig.type_symbol,
                "{}: round-trip site mismatch",
                case.name
            );
            for (a, b) in [
                (rt.frac_x, orig.frac_x),
                (rt.frac_y, orig.frac_y),
                (rt.frac_z, orig.frac_z),
            ] {
                let dev = (a - b).abs();
                max_roundtrip_dev = max_roundtrip_dev.max(dev);
                assert!(
                    dev <= 1e-6,
                    "{}: round-trip fractional coordinate drifted by {dev:.3e}",
                    case.name
                );
            }
        }
    }

    println!(
        "F65: {} CIFs match {} (max cell-parameter dev {:.3e}, max round-trip dev {:.3e})",
        fixture.cases.len(),
        fixture.reference,
        max_param_dev,
        max_roundtrip_dev
    );
    record(
        "F65",
        "cif_roundtrip_vs_pymatgen",
        "max_abs_lattice_param",
        max_param_dev,
        1e-6,
        max_param_dev <= 1e-6 && max_roundtrip_dev <= 1e-6,
        fixture.cases.len(),
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// F66 — minimum-image distances vs ASE `find_mic`
// ═══════════════════════════════════════════════════════════════════════════

#[derive(Deserialize)]
struct PbcFixture {
    reference: String,
    cases: Vec<PbcCase>,
}

#[derive(Deserialize)]
struct PbcCase {
    name: String,
    lattice: LatticeParams,
    frac_coords: Vec<[f64; 3]>,
    pairs: Vec<MicPair>,
}

#[derive(Deserialize)]
struct MicPair {
    i: usize,
    j: usize,
    mic_distance: f64,
}

#[derive(Deserialize, Clone, Copy)]
struct LatticeParams {
    a: f64,
    b: f64,
    c: f64,
    alpha: f64,
    beta: f64,
    gamma: f64,
}

impl LatticeParams {
    fn cell_parameters(&self) -> CellParameters {
        CellParameters {
            a: self.a,
            b: self.b,
            c: self.c,
            alpha: self.alpha,
            beta: self.beta,
            gamma: self.gamma,
        }
    }
}

/// F66 — the minimum-image distance must match ASE's `find_mic` for every
/// pair of the reference corpus (≤ 1e-10 Å), and the implementation must
/// respect the PBC invariants: wrapping a point into the cell or translating
/// it by a lattice vector cannot change the distance.
#[test]
fn f66_minimum_image_matches_pymatgen() {
    let fixture: PbcFixture = load("pymatgen_pbc.json");
    assert!(!fixture.cases.is_empty());

    let mut max_dev = 0.0f64;
    let mut n_pairs = 0usize;

    for case in &fixture.cases {
        let cell = UnitCell::from_parameters(&case.lattice.cell_parameters());
        let cart: Vec<[f64; 3]> = case
            .frac_coords
            .iter()
            .map(|frac| cell.frac_to_cart(*frac))
            .collect();

        for pair in &case.pairs {
            let (i, j) = (pair.i, pair.j);
            assert!(
                i < cart.len() && j < cart.len(),
                "{}: pair ({i}, {j}) out of range",
                case.name
            );
            let (pi, pj) = (cart[i], cart[j]);
            let ours = cell.minimum_image_distance(pi, pj);
            let dev = (ours - pair.mic_distance).abs();
            max_dev = max_dev.max(dev);
            assert!(
                dev <= 1e-10,
                "{}: pair ({i}, {j}) minimum-image distance {:.12} vs ASE {:.12} (dev {dev:.3e})",
                case.name,
                ours,
                pair.mic_distance
            );
            n_pairs += 1;

            // ── PBC invariants ──
            // (a) symmetry
            let mirrored = cell.minimum_image_distance(pj, pi);
            assert!(
                (mirrored - ours).abs() <= 1e-10,
                "{}: pair ({i}, {j}) minimum image is not symmetric",
                case.name
            );
            // (b) wrapping into the primary cell cannot change the distance
            let wrapped = cell.wrap_cart(pj);
            let after_wrap = cell.minimum_image_distance(pi, wrapped);
            assert!(
                (after_wrap - ours).abs() <= 1e-10,
                "{}: pair ({i}, {j}) distance changed after wrapping ({ours} → {after_wrap})",
                case.name
            );
            // (c) translating by a lattice vector cannot change the distance
            for shift in [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [-1.0, 2.0, -3.0],
            ] {
                let moved = [
                    pj[0]
                        + shift[0] * cell.lattice[0][0]
                        + shift[1] * cell.lattice[1][0]
                        + shift[2] * cell.lattice[2][0],
                    pj[1]
                        + shift[0] * cell.lattice[0][1]
                        + shift[1] * cell.lattice[1][1]
                        + shift[2] * cell.lattice[2][1],
                    pj[2]
                        + shift[0] * cell.lattice[0][2]
                        + shift[1] * cell.lattice[1][2]
                        + shift[2] * cell.lattice[2][2],
                ];
                let shifted = cell.minimum_image_distance(pi, moved);
                assert!(
                    (shifted - ours).abs() <= 1e-10,
                    "{}: pair ({i}, {j}) distance changed after a lattice translation ({ours} → {shifted})",
                    case.name
                );
            }
            // (d) the minimum image is never longer than the direct distance
            let direct = {
                let d = [pj[0] - pi[0], pj[1] - pi[1], pj[2] - pi[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            };
            assert!(
                ours <= direct + 1e-12,
                "{}: pair ({i}, {j}) minimum image {ours} exceeds the direct distance {direct}",
                case.name
            );
        }
    }

    // ── skewed cell: the fractional shortcut is not the closest image ──
    // a = b = 1 Å, γ = 60°: the closest image of the fractional point
    // (0.5, 0.5, 0) is (0.5, -0.5, 0), i.e. 0.5 Å away — not the 0.866 Å
    // the "round each fractional component" shortcut returns.
    let skew = UnitCell::from_parameters(&CellParameters {
        a: 1.0,
        b: 1.0,
        c: 2.0,
        alpha: 90.0,
        beta: 90.0,
        gamma: 60.0,
    });
    let p0 = skew.frac_to_cart([0.0, 0.0, 0.0]);
    let p1 = skew.frac_to_cart([0.5, 0.5, 0.0]);
    let skew_distance = skew.minimum_image_distance(p0, p1);
    assert!(
        (skew_distance - 0.5).abs() <= 1e-12,
        "skewed-cell minimum image should be 0.5 Å, got {skew_distance:.12}"
    );

    assert_eq!(
        n_pairs, 44,
        "the reference corpus must cover all 44 atom pairs"
    );

    println!(
        "F66: {n_pairs} pairs match {} (max dev {max_dev:.3e} Å); skewed-cell MIC {skew_distance:.6} Å",
        fixture.reference
    );
    record(
        "F66",
        "minimum_image_distance_vs_ase",
        "max_abs",
        max_dev,
        1e-10,
        max_dev <= 1e-10,
        fixture.cases.len(),
    );
}

// ═══════════════════════════════════════════════════════════════════════════
// F69 — powder XRD vs pymatgen XRDCalculator (Cu-Kα)
// ═══════════════════════════════════════════════════════════════════════════

#[derive(Deserialize)]
struct XrdFixture {
    reference: String,
    wavelength: f64,
    cases: Vec<XrdCase>,
}

#[derive(Deserialize)]
struct XrdCase {
    name: String,
    lattice: LatticeParams,
    elements: Vec<String>,
    frac_coords: Vec<[f64; 3]>,
    two_theta: Vec<f64>,
    intensity: Vec<f64>,
}

/// F69 — the simulated powder pattern must reproduce the pymatgen pattern:
/// same peaks (no missing and no extra reflections), positions within
/// 0.01° and intensities with a correlation ≥ 0.95.
#[test]
fn f69_powder_xrd_matches_pymatgen() {
    let fixture: XrdFixture = load("pymatgen_xrd.json");
    assert!(!fixture.cases.is_empty());
    assert!(
        (fixture.wavelength - 1.54184).abs() <= 1e-9,
        "unexpected reference wavelength {}",
        fixture.wavelength
    );

    let mut max_two_theta_dev = 0.0f64;
    let mut max_intensity_dev = 0.0f64;
    let mut min_correlation = 1.0f64;

    for case in &fixture.cases {
        let cell = UnitCell::from_parameters(&case.lattice.cell_parameters());
        let elements: Vec<u8> = case.elements.iter().map(|s| atomic_number(s)).collect();
        assert_eq!(elements.len(), case.frac_coords.len());

        let ours = simulate_powder_xrd(&cell, &elements, &case.frac_coords, 90.0)
            .expect("the powder XRD simulation must succeed");
        assert_eq!(
            ours.two_theta.len(),
            case.two_theta.len(),
            "{}: {} peaks simulated, reference has {}",
            case.name,
            ours.two_theta.len(),
            case.two_theta.len()
        );

        for (idx, ref_two_theta) in case.two_theta.iter().enumerate() {
            let dev = (ours.two_theta[idx] - ref_two_theta).abs();
            max_two_theta_dev = max_two_theta_dev.max(dev);
            assert!(
                dev <= 0.01,
                "{}: peak {idx} at {:.5}° vs reference {:.5}° (dev {dev:.3e}°)",
                case.name,
                ours.two_theta[idx],
                ref_two_theta
            );
        }

        // The intensities must be reproduced, not merely correlated: the
        // reference stores full f64 precision, so a 1e-6 window on the
        // 0–100 scale is still a machine-precision-level agreement.
        let mut case_max_intensity_dev = 0.0f64;
        for (idx, ref_intensity) in case.intensity.iter().enumerate() {
            let dev = (ours.intensities[idx] - ref_intensity).abs();
            case_max_intensity_dev = case_max_intensity_dev.max(dev);
        }
        max_intensity_dev = max_intensity_dev.max(case_max_intensity_dev);
        assert!(
            case_max_intensity_dev <= 1e-6,
            "{}: intensity deviation {case_max_intensity_dev:.3e} exceeds 1e-6",
            case.name
        );

        let correlation = pearson(&ours.intensities, &case.intensity);
        min_correlation = min_correlation.min(correlation);
        assert!(
            correlation >= 0.95,
            "{}: intensity correlation {correlation:.4} below the 0.95 gate",
            case.name
        );
    }

    println!(
        "F69: {} patterns match {} (max Δ2θ {:.3e}°, max ΔI {:.3e}, min correlation {:.6})",
        fixture.cases.len(),
        fixture.reference,
        max_two_theta_dev,
        max_intensity_dev,
        min_correlation
    );
    record(
        "F69",
        "powder_xrd_vs_pymatgen",
        "max_abs_two_theta",
        max_two_theta_dev,
        0.01,
        max_two_theta_dev <= 0.01 && min_correlation >= 0.95,
        fixture.cases.len(),
    );
}
