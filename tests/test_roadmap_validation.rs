/// Comprehensive validation tests for ALL roadmap features.
///
/// This test file validates each checklist item from the roadmap
/// (08-future-methods-and-metal-support.md) against reference values
/// and physical invariants.
///
/// Run: cargo test --test test_roadmap_validation

// ─── Phase 1: Graceful Metal Complex Support ─────────────────────────────────

#[test]
fn phase1_detect_transition_metals_automatically() {
    let elements_organic = [6u8, 1, 1, 1, 1]; // CH4
    let elements_metal = [26u8, 17, 17, 17, 17, 17, 17]; // FeCl6

    let support_organic = sci_form::get_eht_support(&elements_organic);
    assert!(!support_organic.has_transition_metals);
    assert_eq!(support_organic.level, sci_form::eht::SupportLevel::Supported);

    let support_metal = sci_form::get_eht_support(&elements_metal);
    assert!(support_metal.has_transition_metals);
    assert_eq!(
        support_metal.level,
        sci_form::eht::SupportLevel::Experimental
    );
    assert!(support_metal.provisional_elements.contains(&26));
}

#[test]
fn phase1_uff_fallback_for_unsupported_eht() {
    // Pt system with experimental EHT disabled → should route to UFF
    let elements = [78u8, 17, 17, 7, 7, 1, 1, 1, 1, 1, 1];
    let positions = [
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0],
        [-2.30, 0.0, 0.0],
        [0.0, 2.00, 0.0],
        [0.0, -2.00, 0.0],
        [0.90, 2.60, 0.0],
        [-0.90, 2.60, 0.0],
        [0.00, 2.00, 0.95],
        [0.90, -2.60, 0.0],
        [-0.90, -2.60, 0.0],
        [0.00, -2.00, -0.95],
    ];

    let result =
        sci_form::compute_eht_or_uff_fallback("[Pt](Cl)(Cl)([NH3])[NH3]", &elements, &positions, false);
    assert!(result.is_ok());
    match result.unwrap() {
        sci_form::ElectronicWorkflowResult::UffFallback { energy_kcal_mol, .. } => {
            assert!(energy_kcal_mol.is_finite());
        }
        _ => panic!("Expected UFF fallback when experimental EHT is disabled"),
    }
}

#[test]
fn phase1_capability_query_per_element() {
    let elements = [6u8, 1, 1, 1, 1];
    let caps = sci_form::get_system_capabilities(&elements);
    assert!(caps.embed.available);
    assert!(caps.uff.available);
    assert!(caps.eht.available);
    assert!(caps.population.available);
    assert!(caps.orbital_grid.available);
}

#[test]
fn phase1_low_confidence_mode_for_metals() {
    let elements = [78u8, 17, 17];
    let caps = sci_form::get_system_capabilities(&elements);
    assert_eq!(caps.eht.confidence, sci_form::eht::SupportLevel::Experimental);
    assert!(!caps.eht.warnings.is_empty());
}

// ─── Phase 2: EHT d-Orbital Extension ───────────────────────────────────────

#[test]
fn phase2_d_orbitals_in_basis_set() {
    // Fe should have 4s(1) + 4p(3) + 3d(5) = 9 basis functions
    let elements = [26u8];
    let positions = [[0.0, 0.0, 0.0]];
    let basis = sci_form::eht::basis::build_basis(&elements, &positions);
    assert_eq!(basis.len(), 9);

    // Verify we have s, p, and d orbitals
    let l_values: Vec<u8> = basis.iter().map(|ao| ao.l).collect();
    assert!(l_values.contains(&0), "Missing s orbital");
    assert!(l_values.contains(&1), "Missing p orbital");
    assert!(l_values.contains(&2), "Missing d orbital");
    assert_eq!(l_values.iter().filter(|&&l| l == 2).count(), 5, "Should have 5 d orbitals");
}

#[test]
fn phase2_first_row_tm_coverage() {
    // Sc(21) through Zn(30) should all be parameterized
    for z in 21..=30u8 {
        let params = sci_form::eht::params::get_params(z);
        assert!(
            params.is_some(),
            "Missing EHT parameters for Z={}",
            z
        );
    }
}

#[test]
fn phase2_second_row_tm_coverage() {
    // Y(39) through Cd(48) should all be parameterized
    for z in 39..=48u8 {
        let params = sci_form::eht::params::get_params(z);
        assert!(
            params.is_some(),
            "Missing EHT parameters for Z={}",
            z
        );
    }
}

#[test]
fn phase2_third_row_tm_subset_coverage() {
    // Hf(72) through Hg(80) should be parameterized
    for z in 72..=80u8 {
        let params = sci_form::eht::params::get_params(z);
        assert!(
            params.is_some(),
            "Missing EHT parameters for Z={}",
            z
        );
    }
}

#[test]
fn phase2_valence_electron_counting_5d_metals() {
    // Critical regression: 5d metals must have correct group-based electron counts
    // Hf(72) is group 4 with 4 electrons, not group 3 with 3 electrons
    let cisplatin_elements = [78u8, 17, 17, 7, 7, 1, 1, 1, 1, 1, 1];
    let positions: Vec<[f64; 3]> = vec![
        [0.0, 0.0, 0.0], [2.30, 0.0, 0.0], [-2.30, 0.0, 0.0],
        [0.0, 2.00, 0.0], [0.0, -2.00, 0.0],
        [0.90, 2.60, 0.0], [-0.90, 2.60, 0.0], [0.00, 2.00, 0.95],
        [0.90, -2.60, 0.0], [-0.90, -2.60, 0.0], [0.00, -2.00, -0.95],
    ];

    let result = sci_form::eht::solve_eht(&cisplatin_elements, &positions, None).unwrap();
    // Pt(10) + 2×Cl(7) + 2×N(5) + 6×H(1) = 40 (must be even)
    assert_eq!(result.n_electrons, 40, "Cisplatin must have 40 valence electrons (even)");
    assert_eq!(result.n_electrons % 2, 0, "Closed-shell system must have even electrons");
}

#[test]
fn phase2_ferrocene_has_nonzero_gap() {
    let mut elements = vec![26u8];
    let mut positions = vec![[0.0, 0.0, 0.0]];
    let ring_r_c = 1.65;
    let ring_r_h = 2.70;
    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians();
        elements.push(6);
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), 1.66]);
        elements.push(1);
        positions.push([ring_r_h * a.cos(), ring_r_h * a.sin(), 2.01]);
    }
    for i in 0..5 {
        let a = (i as f64) * 72.0_f64.to_radians() + 36.0_f64.to_radians();
        elements.push(6);
        positions.push([ring_r_c * a.cos(), ring_r_c * a.sin(), -1.66]);
        elements.push(1);
        positions.push([ring_r_h * a.cos(), ring_r_h * a.sin(), -2.01]);
    }

    let result = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
    assert_eq!(result.n_electrons, 58); // Fe(8) + 10C(40) + 10H(10)
    // With proper s/p/d basis, ferrocene should have a nonzero HOMO-LUMO gap
    assert!(result.gap > 0.001, "Ferrocene gap should be > 0.001 eV, got {}", result.gap);
}

#[test]
fn phase2_overlap_supports_spd_combinations() {
    let elements = [26u8, 17]; // Fe-Cl
    let positions = [[0.0, 0.0, 0.0], [2.30, 0.0, 0.0]];
    let basis = sci_form::eht::basis::build_basis(&elements, &positions);
    let s = sci_form::eht::build_overlap_matrix(&basis);

    // Overlap matrix should be symmetric
    for i in 0..s.nrows() {
        for j in 0..s.ncols() {
            assert!(
                (s[(i, j)] - s[(j, i)]).abs() < 1e-10,
                "Overlap matrix not symmetric at ({},{})",
                i, j
            );
        }
    }
    // Diagonal should be 1.0 (normalized)
    for i in 0..s.nrows() {
        assert!(
            (s[(i, i)] - 1.0).abs() < 1e-10,
            "Overlap diagonal not 1.0 at {}",
            i
        );
    }
}

#[test]
fn phase2_orbital_grid_supports_d_orbitals() {
    let elements = [26u8, 17, 17, 17, 17, 17, 17];
    let positions = [
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0], [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0], [0.0, -2.30, 0.0],
        [0.0, 0.0, 2.30], [0.0, 0.0, -2.30],
    ];

    let result = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
    let basis = sci_form::eht::basis::build_basis(&elements, &positions);
    // Should be able to generate an orbital grid for HOMO
    let grid = sci_form::eht::evaluate_orbital_on_grid(
        &basis, &result.coefficients, result.homo_index, &positions, 0.4, 3.0,
    );
    let mesh = sci_form::eht::marching_cubes(&grid, 0.02);
    assert!(mesh.vertices.len() > 0, "Orbital mesh should have vertices for d-orbital system");
}

#[test]
fn phase2_energy_ordering_nondecreasing() {
    let elements = [26u8, 17, 17, 17, 17, 17, 17];
    let positions = [
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0], [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0], [0.0, -2.30, 0.0],
        [0.0, 0.0, 2.30], [0.0, 0.0, -2.30],
    ];
    let result = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
    for w in result.energies.windows(2) {
        assert!(w[0] <= w[1] + 1e-10, "Energies not sorted: {} > {}", w[0], w[1]);
    }
}

// ─── Phase 3: Reactivity Analysis ───────────────────────────────────────────

#[test]
fn phase3_frontier_descriptors_exist() {
    let elements = [8u8, 1, 1]; // water
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_frontier_descriptors(&elements, &positions).unwrap();
    assert_eq!(result.num_atoms, 3);
    assert!(result.gap > 0.0);
    // HOMO contributions should sum to ~1
    let homo_sum: f64 = result.homo_atom_contributions.iter().sum();
    assert!((homo_sum - 1.0).abs() < 1e-6);
    // LUMO contributions should sum to ~1
    let lumo_sum: f64 = result.lumo_atom_contributions.iter().sum();
    assert!((lumo_sum - 1.0).abs() < 1e-6);
}

#[test]
fn phase3_homo_lumo_gap_automatically_computed() {
    let elements = [6u8, 1, 1, 1, 1]; // CH4
    let positions = [
        [0.0, 0.0, 0.0],
        [0.6, 0.6, 0.6], [-0.6, -0.6, 0.6],
        [0.6, -0.6, -0.6], [-0.6, 0.6, -0.6],
    ];
    let result = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
    // Gap should be non-negative
    assert!(result.gap >= 0.0);
    assert_eq!(result.gap, result.lumo_energy - result.homo_energy);
}

#[test]
fn phase3_fukui_descriptors_per_atom() {
    let elements = [8u8, 1, 1]; // water
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_fukui_descriptors(&elements, &positions).unwrap();

    assert_eq!(result.num_atoms, 3);
    assert_eq!(result.f_plus.len(), 3);
    assert_eq!(result.f_minus.len(), 3);
    assert_eq!(result.f_radical.len(), 3);
    assert_eq!(result.dual_descriptor.len(), 3);
    assert_eq!(result.condensed.len(), 3);

    // f+ sums to ~1 (normalized)
    let fp_sum: f64 = result.f_plus.iter().sum();
    assert!((fp_sum - 1.0).abs() < 1e-6);
    // f- sums to ~1 (normalized)
    let fm_sum: f64 = result.f_minus.iter().sum();
    assert!((fm_sum - 1.0).abs() < 1e-6);
    // f_radical = 0.5*(f+ + f-)
    for i in 0..3 {
        let expected = 0.5 * (result.f_plus[i] + result.f_minus[i]);
        assert!((result.f_radical[i] - expected).abs() < 1e-10);
    }
}

#[test]
fn phase3_condensed_fukui_descriptors_per_atom() {
    let elements = [6u8, 8, 8, 1]; // formic acid simplified
    let positions = [
        [0.0, 0.0, 0.0],
        [1.2, 0.0, 0.0],
        [-0.5, 1.0, 0.0],
        [-0.5, -1.0, 0.0],
    ];
    let result = sci_form::compute_fukui_descriptors(&elements, &positions).unwrap();

    for atom in &result.condensed {
        assert!(atom.f_plus >= 0.0, "f+ should be non-negative for atom {}", atom.atom_index);
        assert!(atom.f_minus >= 0.0, "f- should be non-negative for atom {}", atom.atom_index);
        assert!(atom.f_radical >= 0.0, "f_radical should be non-negative for atom {}", atom.atom_index);
        assert_eq!(atom.dual_descriptor, atom.f_plus - atom.f_minus);
    }
}

#[test]
fn phase3_reactivity_ranking_helpers() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_reactivity_ranking(&elements, &positions).unwrap();

    // Should have rankings for all 3 atoms for each attack type
    assert_eq!(result.nucleophilic_attack_sites.len(), 3);
    assert_eq!(result.electrophilic_attack_sites.len(), 3);
    assert_eq!(result.radical_attack_sites.len(), 3);

    // Rankings should be sorted descending by score
    for sites in [
        &result.nucleophilic_attack_sites,
        &result.electrophilic_attack_sites,
        &result.radical_attack_sites,
    ] {
        for w in sites.windows(2) {
            assert!(w[0].score >= w[1].score, "Rankings not sorted descending");
        }
    }
}

#[test]
fn phase3_empirical_pka_acetic_acid() {
    let result = sci_form::compute_empirical_pka("CC(=O)O").unwrap();
    // Acetic acid should have at least one acidic site
    assert!(
        !result.acidic_sites.is_empty(),
        "Acetic acid should detect acidic sites"
    );
    // The carboxylic acid pKa should be roughly in range [2, 7]
    let best_pka = result.acidic_sites[0].estimated_pka;
    assert!(
        best_pka >= 2.0 && best_pka <= 7.0,
        "Acetic acid pKa should be ~4.7, got {}",
        best_pka
    );
}

#[test]
fn phase3_empirical_pka_amine() {
    let result = sci_form::compute_empirical_pka("CCN").unwrap();
    // Ethylamine should have at least one basic site
    assert!(
        !result.basic_sites.is_empty(),
        "Ethylamine should detect basic sites"
    );
    let best_pka = result.basic_sites[0].estimated_pka;
    assert!(
        best_pka >= 8.0 && best_pka <= 14.0,
        "Ethylamine pKa should be ~10.6, got {}",
        best_pka
    );
}

#[test]
fn phase3_uv_vis_spectrum() {
    let conf = sci_form::embed("c1ccccc1", 42);
    assert!(conf.error.is_none());
    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let spectrum = sci_form::compute_uv_vis_spectrum(&conf.elements, &positions, 0.3, 0.0, 20.0, 100).unwrap();
    assert_eq!(spectrum.energies_ev.len(), 100);
    assert_eq!(spectrum.intensities.len(), 100);
    // Should have at least some peaks
    assert!(!spectrum.peaks.is_empty(), "Benzene should produce UV-Vis transitions");
    // All peaks should have finite values
    for peak in &spectrum.peaks {
        assert!(peak.energy_ev.is_finite());
        assert!(peak.wavelength_nm.is_finite());
        assert!(peak.intensity.is_finite());
        assert!(peak.intensity > 0.0);
    }
}

// ─── Phase 4: Bonding and Topology Analysis ─────────────────────────────────

#[test]
fn phase4_bond_order_analysis() {
    let elements = [1u8, 1]; // H2
    let positions = [[0.0, 0.0, 0.0], [0.74, 0.0, 0.0]];
    let result = sci_form::compute_bond_orders(&elements, &positions).unwrap();

    assert!(!result.bonds.is_empty());
    // H-H bond order should be positive and roughly ~1
    let h_h_bond = &result.bonds[0];
    assert_eq!(h_h_bond.atom_i, 0);
    assert_eq!(h_h_bond.atom_j, 1);
    assert!(h_h_bond.wiberg > 0.5, "H-H Wiberg BO should be > 0.5, got {}", h_h_bond.wiberg);
    assert!(h_h_bond.mayer > 0.5, "H-H Mayer BO should be > 0.5, got {}", h_h_bond.mayer);
}

#[test]
fn phase4_bond_order_ethylene_double_bond() {
    let conf = sci_form::embed("C=C", 42);
    assert!(conf.error.is_none());
    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let result = sci_form::compute_bond_orders(&conf.elements, &positions).unwrap();

    // Find C-C bond (elements 6-6)
    let cc_bond = result.bonds.iter().find(|b| {
        conf.elements[b.atom_i] == 6 && conf.elements[b.atom_j] == 6
    });
    assert!(cc_bond.is_some(), "Should find C=C bond");
    let cc = cc_bond.unwrap();
    // C=C double bond → Wiberg BO should be > 1.2
    assert!(cc.wiberg > 1.2, "C=C Wiberg BO should be > 1.2, got {}", cc.wiberg);
}

#[test]
fn phase4_aromaticity_detection_benzene() {
    let features = sci_form::analyze_graph_features("c1ccccc1").unwrap();
    // All 6 carbons should be aromatic
    let aromatic_count = features.aromaticity.aromatic_atoms.iter().filter(|&&a| a).count();
    assert!(aromatic_count >= 6, "Benzene should have ≥6 aromatic atoms, got {}", aromatic_count);
    // Should have aromatic bonds
    assert!(!features.aromaticity.aromatic_bonds.is_empty(), "Benzene should have aromatic bonds");
}

#[test]
fn phase4_aromaticity_detection_methane_not_aromatic() {
    let features = sci_form::analyze_graph_features("C").unwrap();
    let aromatic_count = features.aromaticity.aromatic_atoms.iter().filter(|&&a| a).count();
    assert_eq!(aromatic_count, 0, "Methane should have no aromatic atoms");
}

#[test]
fn phase4_stereocenter_detection() {
    // Alanine: C(N)(C)(=O)O has a stereocenter at the alpha carbon
    // Use a molecule with explicit chirality tag
    let features = sci_form::analyze_graph_features("C([C@@H](N)C)O").unwrap();
    assert!(
        !features.stereocenters.tagged_tetrahedral_centers.is_empty(),
        "Should detect tagged tetrahedral stereocenter"
    );
}

#[test]
fn phase4_coordination_geometry_octahedral() {
    let elements = [26u8, 17, 17, 17, 17, 17, 17]; // FeCl6
    let positions = [
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0], [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0], [0.0, -2.30, 0.0],
        [0.0, 0.0, 2.30], [0.0, 0.0, -2.30],
    ];
    let result = sci_form::compute_topology(&elements, &positions);
    assert_eq!(result.metal_centers.len(), 1);
    assert_eq!(result.metal_centers[0].coordination_number, 6);
    assert_eq!(
        result.metal_centers[0].geometry,
        sci_form::topology::CoordinationGeometryGuess::Octahedral
    );
    assert!(result.metal_centers[0].geometry_score > 0.8);
}

#[test]
fn phase4_coordination_geometry_square_planar() {
    let elements = [78u8, 17, 17, 7, 7, 1, 1, 1, 1, 1, 1]; // cisplatin
    let positions = [
        [0.0, 0.0, 0.0],
        [2.32, 0.0, 0.0], [-2.32, 0.0, 0.0],
        [0.0, 2.05, 0.0], [0.0, -2.05, 0.0],
        [0.90, 2.65, 0.0], [-0.90, 2.65, 0.0], [0.00, 2.05, 0.95],
        [0.90, -2.65, 0.0], [-0.90, -2.65, 0.0], [0.00, -2.05, -0.95],
    ];
    let result = sci_form::compute_topology(&elements, &positions);
    assert_eq!(result.metal_centers.len(), 1);
    assert_eq!(result.metal_centers[0].coordination_number, 4);
    assert_eq!(
        result.metal_centers[0].geometry,
        sci_form::topology::CoordinationGeometryGuess::SquarePlanar
    );
}

#[test]
fn phase4_coordination_geometry_tetrahedral() {
    let elements = [30u8, 17, 17, 17, 17]; // ZnCl4
    let dist = 2.25;
    let t = 0.5773502691896258; // 1/sqrt(3)
    let positions = [
        [0.0, 0.0, 0.0],
        [t * dist, t * dist, t * dist],
        [t * dist, -t * dist, -t * dist],
        [-t * dist, t * dist, -t * dist],
        [-t * dist, -t * dist, t * dist],
    ];
    let result = sci_form::compute_topology(&elements, &positions);
    assert_eq!(result.metal_centers.len(), 1);
    assert_eq!(result.metal_centers[0].coordination_number, 4);
    assert_eq!(
        result.metal_centers[0].geometry,
        sci_form::topology::CoordinationGeometryGuess::Tetrahedral
    );
}

#[test]
fn phase4_topology_output_is_machine_readable() {
    let elements = [26u8, 17, 17, 17, 17, 17, 17];
    let positions = [
        [0.0, 0.0, 0.0],
        [2.30, 0.0, 0.0], [-2.30, 0.0, 0.0],
        [0.0, 2.30, 0.0], [0.0, -2.30, 0.0],
        [0.0, 0.0, 2.30], [0.0, 0.0, -2.30],
    ];
    let result = sci_form::compute_topology(&elements, &positions);
    // Should be serializable to JSON
    let json = serde_json::to_string(&result).unwrap();
    assert!(json.contains("\"geometry\":\"octahedral\""));
    assert!(json.contains("\"coordination_number\":6"));
}

// ─── Phase 5: Dynamics and Conformational Sampling ──────────────────────────

#[test]
fn phase5_conformer_search_with_uff() {
    let result = sci_form::search_conformers_with_uff("CCCC", 5, 42, 0.5).unwrap();
    assert!(result.generated > 0);
    assert!(result.unique > 0);
    assert!(result.unique <= result.generated);
    // Conformers should be ranked by energy (lowest first)
    for w in result.conformers.windows(2) {
        assert!(w[0].energy_kcal_mol <= w[1].energy_kcal_mol + 1e-10);
    }
    // Should have clusters
    assert!(!result.clusters.is_empty());
}

#[test]
fn phase5_conformer_search_rmsd_clustering() {
    // Large threshold → should collapse duplicates
    let result = sci_form::search_conformers_with_uff("CCCC", 10, 42, 100.0).unwrap();
    assert!(result.unique <= result.generated);
    // Very large threshold → all should cluster into 1
    assert!(result.unique <= 3, "Very large RMSD threshold should aggressively deduplicate");
}

#[test]
fn phase5_md_velocity_verlet_nve() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());

    let traj = sci_form::compute_md_trajectory("CCO", &conf.coords, 10, 0.5, 42).unwrap();
    assert_eq!(traj.frames.len(), 11); // frame 0 + 10 steps
    assert_eq!(traj.dt_fs, 0.5);

    // Energy should be finite at all frames
    for frame in &traj.frames {
        assert!(frame.potential_energy_kcal_mol.is_finite());
        assert!(frame.kinetic_energy_kcal_mol.is_finite());
        assert!(frame.temperature_k.is_finite());
    }
}

#[test]
fn phase5_md_nvt_thermostat() {
    let conf = sci_form::embed("C", 42);
    assert!(conf.error.is_none());

    let traj = sci_form::compute_md_trajectory_nvt("C", &conf.coords, 20, 0.5, 42, 300.0, 10.0)
        .unwrap();
    assert_eq!(traj.frames.len(), 21);

    // With thermostat, temperature shouldn't diverge wildly
    for frame in traj.frames.iter().skip(5) {
        assert!(
            frame.temperature_k < 5000.0,
            "Temperature diverged to {} K at step {}",
            frame.temperature_k,
            frame.step
        );
    }
}

#[test]
fn phase5_simplified_neb_path() {
    let conf1 = sci_form::embed("CCO", 42);
    let conf2 = sci_form::embed("CCO", 123);
    assert!(conf1.error.is_none() && conf2.error.is_none());

    let result = sci_form::compute_simplified_neb_path(
        "CCO", &conf1.coords, &conf2.coords, 5, 3, 0.1, 0.001,
    ).unwrap();

    assert_eq!(result.images.len(), 5);
    // First and last images should be close to start/end
    // Energies should all be finite
    for image in &result.images {
        assert!(image.potential_energy_kcal_mol.is_finite());
    }
}

#[test]
fn phase5_conformer_ensemble_format() {
    let result = sci_form::search_conformers_with_uff("CCC", 3, 42, 0.5).unwrap();
    for conf in &result.conformers {
        assert!(conf.energy_kcal_mol.is_finite());
        assert!(!conf.coords.is_empty());
        assert!(conf.cluster_id.is_some());
    }
}

// ─── Phase 6: Methods Availability ──────────────────────────────────────────

#[test]
fn phase6_charges_available() {
    let result = sci_form::compute_charges("CCO").unwrap();
    assert!(!result.charges.is_empty());
    // Charge sum should be ~0 for neutral molecule
    let sum: f64 = result.charges.iter().sum();
    assert!(sum.abs() < 0.1, "Charge sum should be ~0 for neutral CCO, got {}", sum);
}

#[test]
fn phase6_sasa_available() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let result = sci_form::compute_sasa(&conf.elements, &conf.coords, None).unwrap();
    assert!(result.total_sasa > 0.0, "SASA should be positive");
    assert_eq!(result.atom_sasa.len(), conf.num_atoms);
}

#[test]
fn phase6_population_available() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_population(&elements, &positions).unwrap();
    assert_eq!(result.mulliken_charges.len(), 3);
    assert_eq!(result.lowdin_charges.len(), 3);
    // Oxygen should be more negative than hydrogen
    assert!(
        result.mulliken_charges[0] < result.mulliken_charges[1],
        "O should be more negative than H"
    );
}

#[test]
fn phase6_dipole_available() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_dipole(&elements, &positions).unwrap();
    // Water should have a nonzero dipole
    assert!(result.magnitude > 0.0, "Water dipole should be nonzero");
    assert!(result.magnitude < 10.0, "Water dipole should be < 10 D");
}

#[test]
fn phase6_esp_available() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_esp(&elements, &positions, 0.5, 2.0).unwrap();
    assert!(result.dims[0] > 0 && result.dims[1] > 0 && result.dims[2] > 0);
    assert_eq!(result.values.len(), result.dims[0] * result.dims[1] * result.dims[2]);
}

#[test]
fn phase6_dos_available() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let result = sci_form::compute_dos(&elements, &positions, 0.3, -30.0, 5.0, 100).unwrap();
    assert_eq!(result.energies.len(), 100);
    assert_eq!(result.total_dos.len(), 100);
    // DOS result contains smoothed spectrum; verify it has valid values
    assert!(result.total_dos.iter().all(|v| v.is_finite()));
}

#[test]
fn phase6_uff_energy_available() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let energy = sci_form::compute_uff_energy("CCO", &conf.coords).unwrap();
    assert!(energy.is_finite());
}

#[test]
fn phase6_rmsd_available() {
    let conf1 = sci_form::embed("CC", 42);
    let conf2 = sci_form::embed("CC", 123);
    assert!(conf1.error.is_none() && conf2.error.is_none());
    let rmsd = sci_form::compute_rmsd(&conf1.coords, &conf2.coords);
    assert!(rmsd >= 0.0);
    assert!(rmsd.is_finite());
}

// ─── Phase 7: API Method Selection ──────────────────────────────────────────

#[test]
fn phase7_method_selection_organic() {
    let elements = [6u8, 1, 1, 1, 1];
    let plan = sci_form::get_system_method_plan(&elements);
    // For organic: EHT should be recommended for orbitals with "supported" confidence
    assert_eq!(
        plan.orbitals.recommended,
        Some(sci_form::ScientificMethod::Eht)
    );
    let eht_meta = plan.orbitals.methods.iter().find(|m| m.method == sci_form::ScientificMethod::Eht).unwrap();
    assert_eq!(eht_meta.confidence, sci_form::eht::SupportLevel::Supported);
    assert!(eht_meta.confidence_score > 0.9);
}

#[test]
fn phase7_method_selection_metal() {
    let elements = [26u8, 17, 17];
    let plan = sci_form::get_system_method_plan(&elements);
    // For metals: EHT should still be available but experimental
    let eht_meta = plan.orbitals.methods.iter().find(|m| m.method == sci_form::ScientificMethod::Eht).unwrap();
    assert_eq!(eht_meta.confidence, sci_form::eht::SupportLevel::Experimental);
    assert!(eht_meta.confidence_score < 0.9);
}

#[test]
fn phase7_confidence_scoring() {
    let organic = [6u8, 1, 1, 1, 1];
    let metal = [26u8, 17, 17];

    let plan_organic = sci_form::get_system_method_plan(&organic);
    let plan_metal = sci_form::get_system_method_plan(&metal);

    // Organic confidence should be higher than metal
    let organic_eht = plan_organic.orbitals.methods.iter()
        .find(|m| m.method == sci_form::ScientificMethod::Eht).unwrap();
    let metal_eht = plan_metal.orbitals.methods.iter()
        .find(|m| m.method == sci_form::ScientificMethod::Eht).unwrap();

    assert!(organic_eht.confidence_score > metal_eht.confidence_score);
}

#[test]
fn phase7_compare_methods_workflow() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let comparison = sci_form::compare_methods("CCO", &conf.elements, &positions, true);
    assert!(!comparison.comparisons.is_empty());

    // At least UFF and EHT should be present
    let uff_entry = comparison.comparisons.iter().find(|c| c.method == sci_form::ScientificMethod::Uff);
    let eht_entry = comparison.comparisons.iter().find(|c| c.method == sci_form::ScientificMethod::Eht);
    assert!(uff_entry.is_some());
    assert!(eht_entry.is_some());
}

#[test]
fn phase7_method_limitations_exposed() {
    let elements = [26u8, 17, 17];
    let plan = sci_form::get_system_method_plan(&elements);
    let eht_meta = plan.orbitals.methods.iter()
        .find(|m| m.method == sci_form::ScientificMethod::Eht).unwrap();
    // Should have at least one limitation about experimental parameters
    assert!(
        eht_meta.limitations.iter().any(|l| l.contains("provisional") || l.contains("experimental")),
        "Metal EHT should expose limitation about provisional parameters"
    );
}

// ─── Physical Correctness Invariants ─────────────────────────────────────────

#[test]
fn invariant_charge_conservation_neutral_molecule() {
    let result = sci_form::compute_charges("c1ccccc1").unwrap();
    let sum: f64 = result.charges.iter().sum();
    assert!(
        sum.abs() < 0.05,
        "Charge sum for neutral benzene should be ~0, got {}",
        sum
    );
}

#[test]
fn invariant_population_charge_conservation() {
    let elements = [6u8, 8, 1, 1];
    let positions = [
        [0.0, 0.0, 0.0], [1.2, 0.0, 0.0],
        [-0.5, 0.9, 0.0], [-0.5, -0.9, 0.0],
    ];
    let pop = sci_form::compute_population(&elements, &positions).unwrap();
    let mulliken_sum: f64 = pop.mulliken_charges.iter().sum();
    assert!(
        mulliken_sum.abs() < 0.001,
        "Mulliken charge sum should be ~0 for neutral molecule, got {}",
        mulliken_sum
    );
}

#[test]
fn invariant_sasa_positive() {
    let conf = sci_form::embed("CCCCCC", 42);
    assert!(conf.error.is_none());
    let result = sci_form::compute_sasa(&conf.elements, &conf.coords, None).unwrap();
    assert!(result.total_sasa > 0.0);
    for &a in &result.atom_sasa {
        assert!(a >= 0.0, "Per-atom SASA should be non-negative");
    }
}

#[test]
fn invariant_dipole_vector_magnitude_consistency() {
    let elements = [8u8, 1, 1];
    let positions = [[0.0, 0.0, 0.0], [0.757, 0.586, 0.0], [-0.757, 0.586, 0.0]];
    let d = sci_form::compute_dipole(&elements, &positions).unwrap();
    let computed_mag = (d.vector[0].powi(2) + d.vector[1].powi(2) + d.vector[2].powi(2)).sqrt();
    assert!(
        (d.magnitude - computed_mag).abs() < 1e-10,
        "Dipole magnitude should match vector norm"
    );
}

#[test]
fn invariant_symmetric_molecule_dipole_near_zero() {
    let elements = [6u8, 1, 1, 1, 1];
    let positions = [
        [0.0, 0.0, 0.0],
        [0.6, 0.6, 0.6], [-0.6, -0.6, 0.6],
        [0.6, -0.6, -0.6], [-0.6, 0.6, -0.6],
    ];
    let d = sci_form::compute_dipole(&elements, &positions).unwrap();
    assert!(d.magnitude < 2.0, "Methane-like symmetry should give small dipole, got {}", d.magnitude);
}

#[test]
fn invariant_all_energies_finite() {
    let conf = sci_form::embed("CCO", 42);
    assert!(conf.error.is_none());
    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

    let eht = sci_form::eht::solve_eht(&conf.elements, &positions, None).unwrap();
    for e in &eht.energies {
        assert!(e.is_finite(), "Orbital energy must be finite");
    }
    assert!(eht.homo_energy.is_finite());
    assert!(eht.lumo_energy.is_finite());
    assert!(eht.gap.is_finite());
}
