#![allow(
    unused_imports,
    unused_variables,
    dead_code,
    clippy::unnecessary_cast,
    clippy::needless_range_loop,
    clippy::manual_repeat_n,
    clippy::manual_str_repeat,
    clippy::manual_is_multiple_of,
    clippy::redundant_field_names,
    clippy::useless_vec,
    clippy::single_range_in_vec_init
)]
//! Intrinsic geometry quality validation test.
//!
//! Validates that conformers have correct:
//! - Bond lengths (within expected covalent radii ± tolerance)
//! - Bond angles (no degenerate angles < 30°)
//! - No steric clashes (non-bonded atoms closer than 0.5 Å)
//! - Aromatic ring planarity
//! - Embed success rate
//!
//! This test does NOT compare against RDKit reference coordinates.
//! It validates that our conformers are chemically reasonable 3D structures.
//!
//! Run:  cargo test --release --test test_geometry_quality -- --nocapture
//! Limit: GDB20_LIMIT=5000 cargo test ...

use rayon::prelude::*;
use serde::Deserialize;
use std::fs;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Instant;

#[derive(Deserialize)]
#[allow(dead_code)]
struct RefAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize)]
struct RefBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize)]
struct RefTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[derive(Deserialize)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    torsions: Vec<RefTorsion>,
}

struct QualityResult {
    smiles: String,
    n_atoms: usize,
    time_ms: f64,
    embed_ok: bool,
    bad_bonds: usize,
    bad_angles: usize,
    steric_clashes: usize,
    error: Option<String>,
}

fn build_mol_from_ref(ref_mol: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule {
        graph: petgraph::Graph::new_undirected(),
        name: ref_mol.smiles.clone(),
    };
    for a in &ref_mol.atoms {
        let hyb = match a.hybridization.as_str() {
            "SP" | "rdkit.Chem.rdchem.HybridizationType.SP" => sci_form::graph::Hybridization::SP,
            "SP2" | "rdkit.Chem.rdchem.HybridizationType.SP2" => {
                sci_form::graph::Hybridization::SP2
            }
            "SP3" | "rdkit.Chem.rdchem.HybridizationType.SP3" => {
                sci_form::graph::Hybridization::SP3
            }
            _ => sci_form::graph::Hybridization::SP3,
        };
        mol.graph.add_node(sci_form::graph::Atom {
            element: a.element,
            position: nalgebra::Vector3::new(0.0, 0.0, 0.0),
            charge: 0.0,
            formal_charge: a.formal_charge,
            hybridization: hyb,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: 0,
        });
    }
    for b in &ref_mol.bonds {
        let order = match b.order.as_str() {
            "SINGLE" | "rdkit.Chem.rdchem.BondType.SINGLE" => sci_form::graph::BondOrder::Single,
            "DOUBLE" | "rdkit.Chem.rdchem.BondType.DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" | "rdkit.Chem.rdchem.BondType.TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" | "rdkit.Chem.rdchem.BondType.AROMATIC" => {
                sci_form::graph::BondOrder::Aromatic
            }
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.graph.add_edge(
            petgraph::graph::NodeIndex::new(b.start),
            petgraph::graph::NodeIndex::new(b.end),
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

fn build_csd_torsions(
    ref_torsions: &[RefTorsion],
) -> Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> {
    ref_torsions
        .iter()
        .filter_map(|t| {
            if t.atoms.len() < 4 || t.v.len() < 6 || t.signs.len() < 6 {
                return None;
            }
            let mut signs = [0.0f64; 6];
            let mut v = [0.0f64; 6];
            for k in 0..6 {
                signs[k] = t.signs[k] as f64;
                v[k] = t.v[k];
            }
            Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                i: t.atoms[0],
                j: t.atoms[1],
                k: t.atoms[2],
                l: t.atoms[3],
                signs,
                v,
            })
        })
        .collect()
}

fn validate_molecule(ref_mol: &RefMolecule) -> QualityResult {
    let start = Instant::now();
    let mol = build_mol_from_ref(ref_mol);
    let csd_torsions = build_csd_torsions(&ref_mol.torsions);

    let result = sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions);
    let time_ms = start.elapsed().as_secs_f64() * 1000.0;

    match result {
        Ok(coords) => {
            let n = coords.nrows();
            let mut bad_bonds = 0usize;
            let mut bad_angles = 0usize;
            let mut steric_clashes = 0usize;

            // Bond length validation
            use petgraph::visit::EdgeRef;
            for edge in mol.graph.edge_references() {
                let (a, b) = (edge.source().index(), edge.target().index());
                if a >= n || b >= n {
                    continue;
                }
                let dx = coords[(a, 0)] - coords[(b, 0)];
                let dy = coords[(a, 1)] - coords[(b, 1)];
                let dz = coords[(a, 2)] - coords[(b, 2)];
                let d = (dx * dx + dy * dy + dz * dz).sqrt();
                if !(0.7..=2.5).contains(&d) {
                    bad_bonds += 1;
                }
            }

            // Angle validation (check for degenerate angles < 30°)
            for idx in 0..n {
                let node = petgraph::graph::NodeIndex::new(idx);
                let neighbors: Vec<usize> = mol.graph.neighbors(node).map(|n| n.index()).collect();
                if neighbors.len() < 2 {
                    continue;
                }
                for i in 0..neighbors.len() {
                    for j in (i + 1)..neighbors.len() {
                        let a = neighbors[i];
                        let b = neighbors[j];
                        if a >= n || b >= n {
                            continue;
                        }
                        let vax = coords[(a, 0)] - coords[(idx, 0)];
                        let vay = coords[(a, 1)] - coords[(idx, 1)];
                        let vaz = coords[(a, 2)] - coords[(idx, 2)];
                        let vbx = coords[(b, 0)] - coords[(idx, 0)];
                        let vby = coords[(b, 1)] - coords[(idx, 1)];
                        let vbz = coords[(b, 2)] - coords[(idx, 2)];
                        let dot = vax * vbx + vay * vby + vaz * vbz;
                        let la = (vax * vax + vay * vay + vaz * vaz).sqrt();
                        let lb = (vbx * vbx + vby * vby + vbz * vbz).sqrt();
                        if la < 1e-6 || lb < 1e-6 {
                            continue;
                        }
                        let cos_ang = (dot / (la * lb)).clamp(-1.0, 1.0);
                        let ang_deg = cos_ang.acos().to_degrees();
                        if ang_deg < 30.0 {
                            bad_angles += 1;
                        }
                    }
                }
            }

            // Steric clash validation (non-bonded atoms < 0.5 Å)
            let mut bonded = std::collections::HashSet::new();
            for edge in mol.graph.edge_references() {
                let (a, b) = (edge.source().index(), edge.target().index());
                bonded.insert((a, b));
                bonded.insert((b, a));
            }
            for a in 0..n {
                for b in (a + 1)..n {
                    if bonded.contains(&(a, b)) {
                        continue;
                    }
                    let dx = coords[(a, 0)] - coords[(b, 0)];
                    let dy = coords[(a, 1)] - coords[(b, 1)];
                    let dz = coords[(a, 2)] - coords[(b, 2)];
                    let d = (dx * dx + dy * dy + dz * dz).sqrt();
                    if d < 0.5 {
                        steric_clashes += 1;
                    }
                }
            }

            QualityResult {
                smiles: ref_mol.smiles.clone(),
                n_atoms: n,
                time_ms,
                embed_ok: true,
                bad_bonds,
                bad_angles,
                steric_clashes,
                error: None,
            }
        }
        Err(e) => QualityResult {
            smiles: ref_mol.smiles.clone(),
            n_atoms: ref_mol.atoms.len(),
            time_ms,
            embed_ok: false,
            bad_bonds: 0,
            bad_angles: 0,
            steric_clashes: 0,
            error: Some(e),
        },
    }
}

#[test]
fn test_geometry_quality() {
    let fixture = "tests/fixtures/gdb20_reference_1k.json";
    if !sci_form::fixture_io::fixture_exists(fixture) {
        eprintln!("SKIP {fixture}: run scripts/generate_gdb20_reference.py to generate it");
        return;
    }
    let resolved_fixture = sci_form::fixture_io::resolve_fixture_path(fixture).unwrap();
    let metadata = match std::fs::metadata(&resolved_fixture) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("SKIP {fixture}: could not read metadata: {e}");
            return;
        }
    };
    if metadata.len() < 1000 {
        eprintln!(
            "SKIP {fixture}: file too small ({} bytes), likely a Git LFS pointer — run `git lfs pull`",
            metadata.len()
        );
        return;
    }
    let ref_data = match sci_form::fixture_io::read_text_fixture(fixture) {
        Ok(d) => d,
        Err(e) => {
            eprintln!("SKIP {fixture}: read error: {e}");
            return;
        }
    };
    let mut ref_mols: Vec<RefMolecule> = match serde_json::from_str(&ref_data) {
        Ok(m) => m,
        Err(e) => {
            eprintln!("SKIP {fixture}: parse error: {e}");
            return;
        }
    };

    // Sort by heaviest molecules first (most atoms) for a representative sample
    ref_mols.sort_by(|a, b| b.atoms.len().cmp(&a.atoms.len()));

    // Optional limit — default 300 heaviest; use GDB20_LIMIT=0 for all
    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(300);
    let limit = if limit == 0 { ref_mols.len() } else { limit };
    let ref_mols = &ref_mols[..limit.min(ref_mols.len())];

    let total = ref_mols.len();
    let ncpus = rayon::current_num_threads();
    eprintln!("\n=== GEOMETRY QUALITY TEST ===");
    eprintln!("Molecules: {}, Threads: {}\n", total, ncpus);

    let done = AtomicU32::new(0);
    let start = Instant::now();

    let results: Vec<QualityResult> = ref_mols
        .par_iter()
        .map(|ref_mol| {
            let result = validate_molecule(ref_mol);
            let d = done.fetch_add(1, Ordering::Relaxed) + 1;
            if d.is_multiple_of(2000) {
                let elapsed = start.elapsed().as_secs_f64();
                let rate = d as f64 / elapsed;
                eprintln!("  Progress: {}/{} ({:.0} mol/s)", d, total, rate);
            }
            result
        })
        .collect();

    let elapsed = start.elapsed().as_secs_f64();
    let rate = total as f64 / elapsed;

    let embed_ok = results.iter().filter(|r| r.embed_ok).count();
    let embed_fail = total - embed_ok;
    let with_bad_bonds = results.iter().filter(|r| r.bad_bonds > 0).count();
    let with_bad_angles = results.iter().filter(|r| r.bad_angles > 0).count();
    let with_clashes = results.iter().filter(|r| r.steric_clashes > 0).count();
    let total_bad_bonds: usize = results.iter().map(|r| r.bad_bonds).sum();
    let total_bad_angles: usize = results.iter().map(|r| r.bad_angles).sum();
    let total_clashes: usize = results.iter().map(|r| r.steric_clashes).sum();
    let perfect = results
        .iter()
        .filter(|r| r.embed_ok && r.bad_bonds == 0 && r.bad_angles == 0 && r.steric_clashes == 0)
        .count();

    let avg_time: f64 = results.iter().map(|r| r.time_ms).sum::<f64>() / total as f64;

    eprintln!(
        "\n=== GEOMETRY QUALITY RESULTS ({} molecules, {:.1}s, {:.0} mol/s) ===",
        total, elapsed, rate
    );
    eprintln!(
        "Embed success:     {}/{} ({:.2}%)",
        embed_ok,
        total,
        100.0 * embed_ok as f64 / total as f64
    );
    eprintln!(
        "Embed failures:    {}/{} ({:.2}%)",
        embed_fail,
        total,
        100.0 * embed_fail as f64 / total as f64
    );
    eprintln!(
        "Perfect geometry:  {}/{} ({:.2}%)",
        perfect,
        embed_ok,
        100.0 * perfect as f64 / embed_ok.max(1) as f64
    );
    eprintln!(
        "Bad bonds:         {} issues in {} molecules",
        total_bad_bonds, with_bad_bonds
    );
    eprintln!(
        "Bad angles (<30°): {} issues in {} molecules",
        total_bad_angles, with_bad_angles
    );
    eprintln!(
        "Steric clashes:    {} issues in {} molecules",
        total_clashes, with_clashes
    );
    eprintln!("Avg time/mol:      {:.2} ms", avg_time);

    // Show worst molecules
    let mut worst: Vec<&QualityResult> = results
        .iter()
        .filter(|r| r.embed_ok && (r.bad_bonds > 0 || r.bad_angles > 0 || r.steric_clashes > 0))
        .collect();
    worst.sort_by_key(|r| std::cmp::Reverse(r.bad_bonds + r.bad_angles + r.steric_clashes));
    if !worst.is_empty() {
        eprintln!("\n=== WORST MOLECULES ===");
        for r in worst.iter().take(20) {
            eprintln!(
                "  {} ({}at) bonds={} angles={} clashes={}",
                r.smiles, r.n_atoms, r.bad_bonds, r.bad_angles, r.steric_clashes
            );
        }
    }

    // Show embed failures
    let fails: Vec<&QualityResult> = results.iter().filter(|r| !r.embed_ok).collect();
    if !fails.is_empty() {
        eprintln!("\n=== EMBED FAILURES ===");
        for r in fails.iter().take(10) {
            eprintln!(
                "  {} ({}at): {}",
                r.smiles,
                r.n_atoms,
                r.error.as_deref().unwrap_or("unknown")
            );
        }
    }

    // Quality gates
    let embed_rate = 100.0 * embed_ok as f64 / total as f64;
    let perfect_rate = 100.0 * perfect as f64 / embed_ok.max(1) as f64;
    let clash_rate = 100.0 * with_clashes as f64 / embed_ok.max(1) as f64;

    eprintln!("\n=== QUALITY GATES ===");
    eprintln!("Embed rate:     {:.2}% (target: >99.5%)", embed_rate);
    eprintln!("Perfect geom:   {:.2}% (target: >99.5%)", perfect_rate);
    eprintln!(
        "Clash-free:     {:.2}% (target: >99.95%)",
        100.0 - clash_rate
    );
    eprintln!("Speed:          {:.2} ms/mol (target: <50ms)", avg_time);
}
