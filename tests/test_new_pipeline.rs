//! Quick benchmark of the new ETKDG pipeline with retry-on-failure.
//! Tests `generate_3d_conformer` (single call, no brute force) against reference data.

use serde::Deserialize;
use std::collections::HashMap;
use std::fs;

fn calc_dihedral_f32(
    p1: &nalgebra::Vector3<f32>,
    p2: &nalgebra::Vector3<f32>,
    p3: &nalgebra::Vector3<f32>,
    p4: &nalgebra::Vector3<f32>,
) -> f32 {
    let b1 = p2 - p1;
    let b2 = p3 - p2;
    let b3 = p4 - p3;
    let n1 = b1.cross(&b2);
    let n2 = b2.cross(&b3);
    let n1l = n1.norm();
    let n2l = n2.norm();
    if n1l < 1e-6 || n2l < 1e-6 {
        return 0.0;
    }
    let n1u = n1 / n1l;
    let n2u = n2 / n2l;
    let cos_d = n1u.dot(&n2u).clamp(-1.0, 1.0);
    let sign = n1u.dot(&b3);
    let angle = cos_d.acos();
    if sign < 0.0 {
        -angle
    } else {
        angle
    }
}

#[derive(Deserialize, Debug)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize, Debug)]
struct OracleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize, Debug)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

#[derive(Deserialize, Debug)]
struct CsdTorsion {
    atoms: Vec<usize>,
    v: Vec<f64>,
    signs: Vec<i32>,
}

#[test]
fn test_new_pipeline() {
    let data = fs::read_to_string("tests/fixtures/reference_coords_no_mmff.json")
        .expect("Should be able to read reference JSON");
    let mut molecules: Vec<OracleMolecule> =
        serde_json::from_str(&data).expect("JSON was not well-formatted");

    // Load CSD torsion parameters
    let csd_torsions: HashMap<String, Vec<CsdTorsion>> = {
        let path = "tests/fixtures/torsion_params.json";
        match fs::read_to_string(path) {
            Ok(data) => serde_json::from_str(&data).unwrap_or_default(),
            Err(_) => HashMap::new(),
        }
    };
    println!(
        "Loaded CSD torsion params for {} molecules",
        csd_torsions.len()
    );

    // Shuffle and pick 100 random molecules
    use rand::seq::SliceRandom;
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(123);
    molecules.shuffle(&mut rng);
    molecules.truncate(100);

    let mut total_rmsd = 0.0f32;
    let mut max_rmsd = 0.0f32;
    let mut above_05 = 0u32;
    let mut failures = 0u32;
    let mut count = 0u32;

    let start_time = std::time::Instant::now();

    for mol in &molecules {
        let n = mol.atoms.len();
        // Build molecule from oracle data
        let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);
        let mut node_indices = Vec::new();

        for atom in &mol.atoms {
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
            node_indices.push(our_mol.add_atom(new_atom));
        }

        for bond in &mol.bonds {
            let order = match bond.order.as_str() {
                "DOUBLE" => sci_form::graph::BondOrder::Double,
                "TRIPLE" => sci_form::graph::BondOrder::Triple,
                "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
                _ => sci_form::graph::BondOrder::Single,
            };
            our_mol.add_bond(
                node_indices[bond.start],
                node_indices[bond.end],
                sci_form::graph::Bond {
                    order,
                    stereo: sci_form::graph::BondStereo::None,
                },
            );
        }

        // Build CSD torsion contribs for this molecule (if available)
        let csd_contribs: Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> =
            if let Some(torsions) = csd_torsions.get(&mol.smiles) {
                torsions
                    .iter()
                    .map(|t| {
                        let mut signs = [0.0f64; 6];
                        let mut v = [0.0f64; 6];
                        for k in 0..6 {
                            signs[k] = t.signs[k] as f64;
                            v[k] = t.v[k] as f64;
                        }
                        sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                            i: t.atoms[0],
                            j: t.atoms[1],
                            k: t.atoms[2],
                            l: t.atoms[3],
                            signs,
                            v,
                        }
                    })
                    .collect()
            } else {
                Vec::new()
            };

        // Call the new pipeline (single call — internally retries on failure)
        match sci_form::conformer::generate_3d_conformer_with_torsions(&our_mol, 42, &csd_contribs)
        {
            Ok(coords) => {
                // Compute pairwise distance RMSD vs reference
                let mut sq_sum = 0.0f32;
                let mut npairs = 0;
                for i in 0..n {
                    for j in (i + 1)..n {
                        let dr = ((mol.atoms[i].x - mol.atoms[j].x).powi(2)
                            + (mol.atoms[i].y - mol.atoms[j].y).powi(2)
                            + (mol.atoms[i].z - mol.atoms[j].z).powi(2))
                        .sqrt();
                        let du = ((coords[(i, 0)] - coords[(j, 0)]).powi(2)
                            + (coords[(i, 1)] - coords[(j, 1)]).powi(2)
                            + (coords[(i, 2)] - coords[(j, 2)]).powi(2))
                        .sqrt();
                        sq_sum += (dr - du).powi(2);
                        npairs += 1;
                    }
                }
                let rmsd = if npairs > 0 {
                    (sq_sum / npairs as f32).sqrt()
                } else {
                    0.0
                };
                total_rmsd += rmsd;
                if rmsd > max_rmsd {
                    max_rmsd = rmsd;
                }
                if rmsd > 0.5 {
                    above_05 += 1;
                }
                let has_csd = csd_torsions.contains_key(&mol.smiles);
                let tag = if rmsd > 0.5 { "***" } else { "   " };
                println!(
                    "{} Mol {:3} ({:40}) RMSD: {:.3} Å  n={:2} csd={}",
                    tag,
                    count,
                    &mol.smiles[..mol.smiles.len().min(40)],
                    rmsd,
                    n,
                    has_csd,
                );
                // For worst molecules, show detailed diagnostics
                if rmsd > 0.5 {
                    // Build reference DMatrix
                    let mut ref_coords = nalgebra::DMatrix::zeros(n, 3);
                    for ii in 0..n {
                        ref_coords[(ii, 0)] = mol.atoms[ii].x;
                        ref_coords[(ii, 1)] = mol.atoms[ii].y;
                        ref_coords[(ii, 2)] = mol.atoms[ii].z;
                    }

                    // Compare torsion angles for CSD contributions
                    if !csd_contribs.is_empty() {
                        for tc in &csd_contribs {
                            let (i, j, k, l) = (tc.i, tc.j, tc.k, tc.l);
                            // our dihedral
                            let op1 = nalgebra::Vector3::new(
                                coords[(i, 0)],
                                coords[(i, 1)],
                                coords[(i, 2)],
                            );
                            let op2 = nalgebra::Vector3::new(
                                coords[(j, 0)],
                                coords[(j, 1)],
                                coords[(j, 2)],
                            );
                            let op3 = nalgebra::Vector3::new(
                                coords[(k, 0)],
                                coords[(k, 1)],
                                coords[(k, 2)],
                            );
                            let op4 = nalgebra::Vector3::new(
                                coords[(l, 0)],
                                coords[(l, 1)],
                                coords[(l, 2)],
                            );
                            let our_a = calc_dihedral_f32(&op1, &op2, &op3, &op4);
                            // ref dihedral
                            let rp1 = nalgebra::Vector3::new(
                                ref_coords[(i, 0)],
                                ref_coords[(i, 1)],
                                ref_coords[(i, 2)],
                            );
                            let rp2 = nalgebra::Vector3::new(
                                ref_coords[(j, 0)],
                                ref_coords[(j, 1)],
                                ref_coords[(j, 2)],
                            );
                            let rp3 = nalgebra::Vector3::new(
                                ref_coords[(k, 0)],
                                ref_coords[(k, 1)],
                                ref_coords[(k, 2)],
                            );
                            let rp4 = nalgebra::Vector3::new(
                                ref_coords[(l, 0)],
                                ref_coords[(l, 1)],
                                ref_coords[(l, 2)],
                            );
                            let ref_a = calc_dihedral_f32(&rp1, &rp2, &rp3, &rp4);
                            let diff = (our_a - ref_a)
                                .abs()
                                .min((our_a - ref_a + 2.0 * std::f32::consts::PI).abs())
                                .min((our_a - ref_a - 2.0 * std::f32::consts::PI).abs());
                            if diff > 0.3 {
                                println!("      TORSION ({},{},{},{}) ours={:.1}° ref={:.1}° diff={:.1}° v={:?} signs={:?}",
                                    i, j, k, l,
                                    our_a.to_degrees(), ref_a.to_degrees(), diff.to_degrees(),
                                    tc.v.iter().map(|x| (*x * 10.0).round() / 10.0).collect::<Vec<_>>(),
                                    tc.signs.iter().map(|x| *x as i32).collect::<Vec<_>>());
                            }
                        }
                    }
                }
            }
            Err(e) => {
                failures += 1;
                println!(
                    "Mol {:3} ({:40}) FAILED: {}",
                    count,
                    &mol.smiles[..mol.smiles.len().min(40)],
                    e,
                );
            }
        }
        count += 1;
    }

    let elapsed = start_time.elapsed();
    let success_count = count - failures;
    let avg_rmsd = if success_count > 0 {
        total_rmsd / success_count as f32
    } else {
        0.0
    };

    println!("\n=== NEW PIPELINE TEST RESULTS ===");
    println!(
        "Molecules: {}, Successes: {}, Failures: {}",
        count, success_count, failures
    );
    println!(
        "Avg RMSD: {:.3} Å, Max RMSD: {:.3} Å, Above 0.5: {}",
        avg_rmsd, max_rmsd, above_05
    );
    println!(
        "Time: {:.2}s ({:.1} ms/mol)",
        elapsed.as_secs_f64(),
        elapsed.as_secs_f64() * 1000.0 / count as f64
    );
}
