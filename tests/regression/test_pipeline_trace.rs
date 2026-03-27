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
use serde::Deserialize;
/// Pipeline trace diagnostic: dump intermediate values for failing molecules
use std::fs;

#[derive(Deserialize)]
struct RefAtom {
    element: u8,
    x: f64,
    y: f64,
    z: f64,
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
    signs: Vec<f64>,
}

#[derive(Deserialize)]
struct RefMolecule {
    smiles: String,
    atoms: Vec<RefAtom>,
    bonds: Vec<RefBond>,
    #[serde(default)]
    torsions: Vec<RefTorsion>,
}

fn build_mol_from_ref(ref_mol: &RefMolecule) -> sci_form::graph::Molecule {
    let mut mol = sci_form::graph::Molecule::new(&ref_mol.smiles);
    let mut node_indices = Vec::with_capacity(ref_mol.atoms.len());
    for atom in &ref_mol.atoms {
        let hybridization = match atom.hybridization.as_str() {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        let new_atom = sci_form::graph::Atom {
            element: atom.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if atom.element == 1 || atom.element == 0 {
                1
            } else {
                0
            },
        };
        node_indices.push(mol.add_atom(new_atom));
    }
    for bond in &ref_mol.bonds {
        let order = match bond.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        mol.add_bond(
            node_indices[bond.start],
            node_indices[bond.end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    mol
}

fn pairwise_rmsd_f64(coords: &nalgebra::DMatrix<f64>, ref_atoms: &[RefAtom]) -> f64 {
    let n = ref_atoms.len();
    let mut sq_sum = 0.0f64;
    let mut npairs = 0u64;
    for a in 0..n {
        for b in (a + 1)..n {
            let dr = ((ref_atoms[a].x - ref_atoms[b].x).powi(2)
                + (ref_atoms[a].y - ref_atoms[b].y).powi(2)
                + (ref_atoms[a].z - ref_atoms[b].z).powi(2))
            .sqrt();
            let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
            .sqrt();
            sq_sum += (dr - du).powi(2);
            npairs += 1;
        }
    }
    if npairs > 0 {
        (sq_sum / npairs as f64).sqrt()
    } else {
        0.0
    }
}

/// Compare power method eigendecomposition with nalgebra
#[test]
fn test_trace_pipeline() {
    let ref_data = sci_form::fixture_io::read_text_fixture("tests/fixtures/gdb20_reference_1k.json")
        .unwrap();
    let ref_mols: Vec<RefMolecule> = serde_json::from_str(&ref_data).unwrap();

    let limit: usize = std::env::var("GDB20_LIMIT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let threshold = 0.5;

    let mut n_pass = 0;
    let mut n_fail = 0;
    let mut fail_eigenval_stats: Vec<(String, f64, f64, usize)> = Vec::new();

    for (mol_idx, ref_mol) in ref_mols.iter().take(limit).enumerate() {
        let mol = build_mol_from_ref(ref_mol);
        let n = mol.graph.node_count();

        // Compute bounds
        let mut bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
        if !sci_form::distgeom::triangle_smooth_tol(&mut bounds, 0.0) {
            bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, false);
            sci_form::distgeom::triangle_smooth_tol(&mut bounds, 0.0);
        }

        // Compare power vs nalgebra eigendecomposition for the FIRST attempt
        let mut rng = sci_form::distgeom::MinstdRand::new(42);
        let dists = sci_form::distgeom::pick_rdkit_distances(&mut rng, &bounds);
        let rng_state = rng.get_state();
        let coords_power = sci_form::distgeom::compute_initial_coords_rdkit(&mut rng, &dists, 3);
        let mut rng2 = sci_form::distgeom::MinstdRand::new(rng_state as u32);
        let coords_nalgebra =
            sci_form::distgeom::compute_initial_coords_nalgebra_f64(&mut rng2, &dists, 3);

        let max_dist_diff = match (&coords_power, &coords_nalgebra) {
            (Some(cp), Some(cn)) => {
                let mut max_dd = 0.0f64;
                for i in 0..n {
                    for j in (i + 1)..n {
                        let d_p = ((cp[(i, 0)] - cp[(j, 0)]).powi(2)
                            + (cp[(i, 1)] - cp[(j, 1)]).powi(2)
                            + (cp[(i, 2)] - cp[(j, 2)]).powi(2))
                        .sqrt();
                        let d_n = ((cn[(i, 0)] - cn[(j, 0)]).powi(2)
                            + (cn[(i, 1)] - cn[(j, 1)]).powi(2)
                            + (cn[(i, 2)] - cn[(j, 2)]).powi(2))
                        .sqrt();
                        let diff = (d_p - d_n).abs();
                        if diff > max_dd {
                            max_dd = diff;
                        }
                    }
                }
                max_dd
            }
            _ => f64::NAN,
        };

        // Always run the full pipeline to compute RMSD
        let csd_torsions: Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> = ref_mol
            .torsions
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
            .collect();

        let result =
            sci_form::conformer::generate_3d_conformer_with_torsions(&mol, 42, &csd_torsions);
        let rmsd = match &result {
            Ok(coords) => {
                let coords_f64 = coords.map(|v| v as f64);
                pairwise_rmsd_f64(&coords_f64, &ref_mol.atoms)
            }
            Err(_) => f64::NAN,
        };

        if rmsd > threshold || rmsd.is_nan() {
            n_fail += 1;
            // Compute eigenvalue spectrum for failing molecules
            let d_size = n * (n + 1) / 2;
            let mut sq_packed = vec![0.0f64; d_size];
            let mut sum_sq_all = 0.0f64;
            for i in 0..n {
                let id = i * (i + 1) / 2;
                for j in 0..=i {
                    let d = dists[(i, j)];
                    sq_packed[id + j] = d * d;
                    sum_sq_all += d * d;
                }
            }
            sum_sq_all /= (n * n) as f64;
            let mut d0 = vec![0.0f64; n];
            for i in 0..n {
                let mut row_sum = 0.0f64;
                for j in 0..n {
                    let idx = if i >= j {
                        i * (i + 1) / 2 + j
                    } else {
                        j * (j + 1) / 2 + i
                    };
                    row_sum += sq_packed[idx];
                }
                d0[i] = row_sum / n as f64 - sum_sq_all;
            }
            let mut t_full = nalgebra::DMatrix::from_element(n, n, 0.0f64);
            for i in 0..n {
                for j in 0..=i {
                    let sq_val = sq_packed[if i >= j {
                        i * (i + 1) / 2 + j
                    } else {
                        j * (j + 1) / 2 + i
                    }];
                    let val = 0.5 * (d0[i] + d0[j] - sq_val);
                    t_full[(i, j)] = val;
                    t_full[(j, i)] = val;
                }
            }
            let eigen = nalgebra::SymmetricEigen::new(t_full);
            let mut evals: Vec<f64> = eigen.eigenvalues.iter().copied().collect();
            evals.sort_by(|a, b| b.partial_cmp(a).unwrap());

            let neg_count = evals.iter().filter(|&&v| v < -1e-3).count();
            let min_eig = evals.last().copied().unwrap_or(0.0);

            fail_eigenval_stats.push((ref_mol.smiles.clone(), rmsd, min_eig, neg_count));

            if n_fail <= 15 {
                eprintln!("FAIL mol[{}] '{}': RMSD={:.4}, max_pw_dist_diff={:.6e}, neg_eigs={}, min_eig={:.4}",
                    mol_idx, ref_mol.smiles, rmsd, max_dist_diff, neg_count, min_eig);
                let top5: Vec<String> = evals.iter().take(5).map(|v| format!("{:.4}", v)).collect();
                eprintln!("  eigenvalues[0..5] = [{}]", top5.join(", "));
            }
        } else {
            n_pass += 1;
        }
    }

    eprintln!("\n=== SUMMARY ===");
    eprintln!(
        "Pass: {}, Fail: {} ({:.1}%)",
        n_pass,
        n_fail,
        100.0 * n_fail as f64 / (n_pass + n_fail) as f64
    );

    // Analyze correlation between negative eigenvalues and failures
    if !fail_eigenval_stats.is_empty() {
        let avg_min_eig: f64 =
            fail_eigenval_stats.iter().map(|s| s.2).sum::<f64>() / fail_eigenval_stats.len() as f64;
        let avg_neg_count: f64 = fail_eigenval_stats.iter().map(|s| s.3 as f64).sum::<f64>()
            / fail_eigenval_stats.len() as f64;
        eprintln!(
            "Failing molecules: avg_min_eigenvalue={:.4}, avg_neg_eigenvalue_count={:.1}",
            avg_min_eig, avg_neg_count
        );
    }
}
