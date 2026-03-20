use petgraph::visit::EdgeRef;

use sci_form::graph::BondOrder;
use sci_form::materials::{CoordinationGeometry, Sbu, Topology};

fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
    coords
        .chunks_exact(3)
        .map(|chunk| [chunk[0], chunk[1], chunk[2]])
        .collect()
}

fn ml_bonds(mol: &sci_form::graph::Molecule) -> Vec<(usize, usize, u8)> {
    mol.graph
        .edge_references()
        .map(|edge| {
            let order = match edge.weight().order {
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Aromatic => 4,
                BondOrder::Unknown => 0,
            };
            (edge.source().index(), edge.target().index(), order)
        })
        .collect()
}

#[test]
fn smoke_conformer_and_clustering() {
    let ethanol_a = sci_form::embed("CCO", 42);
    let ethanol_b = sci_form::embed("CCO", 43);
    let ethanol_c = sci_form::embed("CCO", 44);

    assert!(
        ethanol_a.error.is_none(),
        "embed failed: {:?}",
        ethanol_a.error
    );
    assert!(
        ethanol_b.error.is_none(),
        "embed failed: {:?}",
        ethanol_b.error
    );
    assert!(
        ethanol_c.error.is_none(),
        "embed failed: {:?}",
        ethanol_c.error
    );
    assert_eq!(ethanol_a.coords.len(), ethanol_a.num_atoms * 3);

    let conformers = vec![
        ethanol_a.coords.clone(),
        ethanol_b.coords.clone(),
        ethanol_c.coords.clone(),
    ];

    let rmsd_matrix = sci_form::compute_rmsd_matrix(&conformers);
    assert_eq!(rmsd_matrix.len(), conformers.len());
    assert!(rmsd_matrix.iter().all(|row| row.len() == conformers.len()));

    let cluster_result = sci_form::butina_cluster(&conformers, 0.75);
    assert_eq!(cluster_result.assignments.len(), conformers.len());
    assert!(!cluster_result.centroid_indices.is_empty());

    let self_rmsd = sci_form::compute_rmsd(&conformers[0], &conformers[0]);
    assert!(self_rmsd.abs() < 1e-8);

    #[cfg(feature = "parallel")]
    {
        let config = sci_form::ConformerConfig {
            seed: 42,
            num_threads: 2,
        };
        let batch = sci_form::embed_batch(&["CCO", "c1ccccc1", "CC(=O)O"], &config);
        assert_eq!(batch.len(), 3);
        assert!(batch.iter().all(|result| result.error.is_none()));
    }
}

#[test]
fn smoke_properties_and_solvation() {
    let ethanol = sci_form::embed("CCO", 42);
    assert!(ethanol.error.is_none(), "embed failed: {:?}", ethanol.error);

    let positions = positions_from_flat(&ethanol.coords);
    let charges = sci_form::compute_charges("CCO").expect("compute_charges");
    assert_eq!(charges.charges.len(), ethanol.num_atoms);
    assert!(charges.total_charge.abs() < 1e-6);

    let sasa = sci_form::compute_sasa(&ethanol.elements, &ethanol.coords, Some(1.4))
        .expect("compute_sasa");
    assert!(sasa.total_sasa > 0.0);
    assert_eq!(sasa.atom_sasa.len(), ethanol.num_atoms);

    let population =
        sci_form::compute_population(&ethanol.elements, &positions).expect("compute_population");
    assert_eq!(population.mulliken_charges.len(), ethanol.num_atoms);
    assert_eq!(population.lowdin_charges.len(), ethanol.num_atoms);

    let dipole = sci_form::compute_dipole(&ethanol.elements, &positions).expect("compute_dipole");
    assert!(dipole.magnitude.is_finite());

    let esp = sci_form::compute_esp(&ethanol.elements, &positions, 0.8, 2.0).expect("compute_esp");
    assert!(!esp.values.is_empty());

    let dos = sci_form::compute_dos(&ethanol.elements, &positions, 0.3, -20.0, 5.0, 64)
        .expect("compute_dos");
    assert_eq!(dos.total_dos.len(), 64);
    assert_eq!(dos.energies.len(), 64);

    let uff = sci_form::compute_uff_energy("CCO", &ethanol.coords).expect("compute_uff_energy");
    let mmff =
        sci_form::compute_mmff94_energy("CCO", &ethanol.coords).expect("compute_mmff94_energy");
    assert!(uff.is_finite());
    assert!(mmff.is_finite());

    let nonpolar = sci_form::compute_nonpolar_solvation(&ethanol.elements, &positions, Some(1.4));
    assert!(nonpolar.total_sasa > 0.0);
    assert_eq!(nonpolar.atom_contributions.len(), ethanol.num_atoms);

    let gb = sci_form::compute_gb_solvation(
        &ethanol.elements,
        &positions,
        &charges.charges,
        Some(78.5),
        Some(1.0),
        Some(1.4),
    );
    assert!(gb.total_energy_kcal_mol.is_finite());
}

#[test]
fn smoke_ml_quantum_and_framework() {
    let ethanol = sci_form::embed("CCO", 42);
    assert!(ethanol.error.is_none(), "embed failed: {:?}", ethanol.error);

    let parsed = sci_form::parse("CCO").expect("parse");
    let bonds = ml_bonds(&parsed);
    let descriptors = sci_form::compute_ml_descriptors(&ethanol.elements, &bonds, &[], &[]);
    let properties = sci_form::predict_ml_properties(&descriptors);
    assert!(properties.druglikeness.is_finite());

    let positions = positions_from_flat(&ethanol.coords);
    let pm3 = sci_form::compute_pm3(&ethanol.elements, &positions).expect("compute_pm3");
    let xtb = sci_form::compute_xtb(&ethanol.elements, &positions).expect("compute_xtb");
    assert!(pm3.converged);
    assert!(xtb.converged);
    assert_eq!(pm3.mulliken_charges.len(), ethanol.num_atoms);
    assert_eq!(xtb.mulliken_charges.len(), ethanol.num_atoms);

    let cell = sci_form::create_unit_cell(12.0, 12.0, 12.0, 90.0, 90.0, 90.0);
    let node = Sbu::metal_node(30, 2.0, CoordinationGeometry::Octahedral);
    let linker = Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
    let topology = Topology::pcu();
    let structure = sci_form::assemble_framework(&node, &linker, &topology, &cell);
    assert!(structure.num_atoms() > 0);
}

#[test]
fn smoke_stereo_rings_and_fingerprints() {
    let chiral = sci_form::embed("C[C@H](F)Cl", 42);
    assert!(chiral.error.is_none(), "embed failed: {:?}", chiral.error);

    let positions = positions_from_flat(&chiral.coords);
    let stereo = sci_form::analyze_stereo("C[C@H](F)Cl", &chiral.coords).expect("analyze_stereo");
    assert!(stereo.n_stereocenters >= 1);
    assert!(stereo
        .stereocenters
        .iter()
        .any(|center| center.configuration.is_some()));

    let rings = sci_form::compute_sssr("c1ccccc1").expect("compute_sssr");
    assert!(!rings.rings.is_empty());
    assert!(
        rings
            .ring_size_histogram
            .get(6)
            .copied()
            .unwrap_or_default()
            >= 1
    );

    let benzene_fp = sci_form::compute_ecfp("c1ccccc1", 2, 2048).expect("compute_ecfp");
    let toluene_fp = sci_form::compute_ecfp("Cc1ccccc1", 2, 2048).expect("compute_ecfp");
    let tanimoto = sci_form::compute_tanimoto(&benzene_fp, &toluene_fp);
    assert!((0.0..1.0).contains(&tanimoto));

    let _ = sci_form::compute_dipole(&chiral.elements, &positions).expect("compute_dipole");
}
