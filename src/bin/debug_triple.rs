//! Debug: inspect torsion contributions during 3D conformer generation
//! for molecules with triple bonds.
//!
//! Category: debug

use sci_form::conformer::generate_3d_conformer_with_torsions;
use sci_form::forcefield::etkdg_3d::M6TorsionContrib;
use sci_form::graph::Molecule;

fn main() {
    let smi = "CCNCC#N";
    let mol = Molecule::from_smiles(smi).unwrap();

    // Print atom ordering
    for i in 0..mol.graph.node_count() {
        let a = &mol.graph[petgraph::graph::NodeIndex::new(i)];
        let nbs: Vec<_> = mol
            .graph
            .neighbors(petgraph::graph::NodeIndex::new(i))
            .map(|n| n.index())
            .collect();
        let elem = match a.element {
            1 => "H",
            6 => "C",
            7 => "N",
            _ => "?",
        };
        println!("Atom {}: {} (nbs: {:?})", i, elem, nbs);
    }

    // CSD torsion params from torsion_params.json
    let csd = vec![
        M6TorsionContrib {
            i: 0,
            j: 1,
            k: 2,
            l: 3,
            signs: [1.0, 1.0, 1.0, -1.0, 1.0, 1.0],
            v: [4.0, 3.1, 3.9, -0.8, 0.0, 0.7],
        },
        M6TorsionContrib {
            i: 4,
            j: 3,
            k: 2,
            l: 1,
            signs: [1.0, 1.0, 1.0, -1.0, 1.0, 1.0],
            v: [4.0, 3.1, 3.9, -0.8, 0.0, 0.7],
        },
    ];

    match generate_3d_conformer_with_torsions(&mol, 42, &csd) {
        Ok(coords) => {
            println!("\nGenerated coordinates:");
            for i in 0..mol.graph.node_count() {
                println!(
                    "Atom {}: ({:.4}, {:.4}, {:.4})",
                    i,
                    coords[(i, 0)],
                    coords[(i, 1)],
                    coords[(i, 2)]
                );
            }

            // Compute torsion angles for CSD bonds
            fn torsion_angle(
                coords: &nalgebra::DMatrix<f32>,
                i: usize,
                j: usize,
                k: usize,
                l: usize,
            ) -> f32 {
                let p = |idx: usize| {
                    nalgebra::Vector3::new(coords[(idx, 0)], coords[(idx, 1)], coords[(idx, 2)])
                };
                let b1 = p(j) - p(i);
                let b2 = p(k) - p(j);
                let b3 = p(l) - p(k);
                let n1 = b1.cross(&b2);
                let n2 = b2.cross(&b3);
                let m = n1.cross(&b2.normalize());
                let x = n1.dot(&n2);
                let y = m.dot(&n2);
                (-y).atan2(x).to_degrees()
            }

            println!(
                "\nTorsion 0-1-2-3: {:.1}°",
                torsion_angle(&coords, 0, 1, 2, 3)
            );
            println!(
                "Torsion 4-3-2-1: {:.1}°",
                torsion_angle(&coords, 4, 3, 2, 1)
            );

            // Compute pairwise distances and compare with RDKit reference
            let rdkit: [(f64, f64, f64); 14] = [
                (-2.0850, 0.4112, -0.1475),
                (-0.6112, 0.3211, 0.2951),
                (0.1450, -0.6726, -0.3864),
                (1.4583, -0.9298, 0.1811),
                (2.3875, 0.1756, 0.0559),
                (3.1279, 1.0676, -0.0355),
                (-2.7179, -0.1176, 0.5883),
                (-2.3963, 1.4886, -0.1075),
                (-2.2142, 0.0674, -1.1895),
                (-0.1698, 1.3385, -0.0045),
                (-0.5778, 0.2989, 1.3753),
                (0.3883, -0.3306, -1.3745),
                (1.8719, -1.7641, -0.4525),
                (1.3933, -1.3542, 1.2021),
            ];
            let n = rdkit.len();
            let mut sq_sum = 0.0f64;
            let mut npairs = 0u64;
            for a in 0..n {
                for b in (a + 1)..n {
                    let dr = ((rdkit[a].0 - rdkit[b].0).powi(2)
                        + (rdkit[a].1 - rdkit[b].1).powi(2)
                        + (rdkit[a].2 - rdkit[b].2).powi(2))
                    .sqrt();
                    let du = ((coords[(a, 0)] - coords[(b, 0)]).powi(2)
                        + (coords[(a, 1)] - coords[(b, 1)]).powi(2)
                        + (coords[(a, 2)] - coords[(b, 2)]).powi(2))
                    .sqrt() as f64;
                    sq_sum += (dr - du).powi(2);
                    npairs += 1;
                }
            }
            println!("\nPairwise RMSD: {:.4}", (sq_sum / npairs as f64).sqrt());
        }
        Err(e) => println!("Failed: {}", e),
    }
}
