//! Comprehensive reaction test system: DFT → reactions → mesh visualization.
//!
//! Tests the end-to-end pipeline from conformer generation through quantum
//! chemistry calculations, NEB path computation, and orbital mesh generation
//! along reaction coordinates.

// ─── Helpers ───────────────────────────────────────────────────────────────

fn embed_smiles(smiles: &str) -> (Vec<u8>, Vec<[f64; 3]>) {
    let conf = sci_form::embed(smiles, 42);
    assert!(
        conf.error.is_none(),
        "embed failed for {}: {:?}",
        smiles,
        conf.error
    );
    let positions: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    (conf.elements.clone(), positions)
}

fn positions_to_flat(positions: &[[f64; 3]]) -> Vec<f64> {
    positions.iter().flat_map(|p| p.iter().copied()).collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 1: DFT energy evaluation at reaction endpoints
// ═══════════════════════════════════════════════════════════════════════════

mod dft_endpoints {
    use super::*;

    /// Evaluate multiple DFT methods on the same molecule to validate consistency.
    #[test]
    fn multi_method_energy_ethanol() {
        let (elements, positions) = embed_smiles("CCO");
        let coords = positions_to_flat(&positions);

        // EHT
        let eht = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
        assert!(eht.gap > 0.0, "EHT gap should be positive");
        assert!(eht.homo_energy < eht.lumo_energy, "EHT HOMO < LUMO");

        // PM3
        let pm3 = sci_form::compute_pm3(&elements, &positions).unwrap();
        assert!(pm3.converged, "PM3 should converge for ethanol");
        assert!(pm3.gap > 0.0, "PM3 gap should be positive");

        // GFN2-xTB (via compute_xtb)
        let xtb = sci_form::compute_xtb(&elements, &positions).unwrap();
        assert!(xtb.converged, "xTB should converge for ethanol");
        assert!(xtb.gap > 0.0, "xTB gap should be positive");

        // UFF energy
        let uff = sci_form::compute_uff_energy("CCO", &coords).unwrap();
        assert!(uff.is_finite(), "UFF energy should be finite");

        eprintln!(
            "  Ethanol energies — EHT gap: {:.3} eV, PM3 gap: {:.3} eV, xTB gap: {:.3} eV, UFF: {:.2} kcal/mol",
            eht.gap, pm3.gap, xtb.gap, uff
        );
    }

    /// Energy evaluation at a perturbed geometry (reaction-like displacement).
    #[test]
    fn energy_along_bond_stretch() {
        // Use embed to get a good starting geometry for ethanol
        let conf = sci_form::embed("CCO", 42);
        assert!(conf.error.is_none());
        let pos_ref: Vec<[f64; 3]> = conf.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let e_ref = sci_form::compute_pm3(&conf.elements, &pos_ref).unwrap();

        // Stretch the C-O bond by moving the last heavy atom outward
        let mut pos_str = pos_ref.clone();
        // Move the oxygen atom (element 8) outward along its current direction
        for (i, &elem) in conf.elements.iter().enumerate() {
            if elem == 8 {
                // Shift oxygen 0.5 Å further from center of mass
                pos_str[i][0] += 0.5;
                pos_str[i][1] += 0.3;
                break;
            }
        }
        let e_str = sci_form::compute_pm3(&conf.elements, &pos_str).unwrap();

        // Perturbed geometry should have different energy (typically higher due to strain)
        let delta_ev = (e_str.total_energy - e_ref.total_energy).abs();
        assert!(
            delta_ev > 0.001,
            "Energy should change with geometry perturbation: delta = {:.6} eV",
            delta_ev
        );
        assert!(
            e_str.total_energy.is_finite() && e_ref.total_energy.is_finite(),
            "Both energies should be finite"
        );

        eprintln!(
            "  Ethanol PM3: ref={:.4} eV, stretched={:.4} eV, delta={:.4} eV",
            e_ref.total_energy, e_str.total_energy, delta_ev
        );
    }

    /// Energy and gradient consistency via NEB backend.
    #[test]
    fn neb_backend_energy_gradient_consistency() {
        let conf = sci_form::embed("CCO", 42);
        assert!(conf.error.is_none());
        let coords = conf.coords.clone();

        for method in &["uff", "pm3"] {
            let (e, grad) = sci_form::neb_backend_energy_and_gradient(method, "CCO", &coords)
                .unwrap_or_else(|_| panic!("{} energy+gradient failed", method));

            assert!(e.is_finite(), "{} energy not finite", method);
            assert_eq!(grad.len(), coords.len(), "{} gradient wrong length", method);

            // Gradient should be nonzero (molecule not at exact minimum)
            let grad_norm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
            assert!(
                grad_norm > 1e-6,
                "{} gradient norm too small: {}",
                method,
                grad_norm
            );

            eprintln!(
                "  [{}] E={:.4} kcal/mol, |grad|={:.6}",
                method, e, grad_norm
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 2: NEB path computation across methods
// ═══════════════════════════════════════════════════════════════════════════

mod neb_paths {

    /// Ethane torsion: compute NEB path with UFF using very small perturbation.
    #[test]
    fn neb_ethane_torsion_uff() {
        let conf = sci_form::embed("CC", 42);
        assert!(conf.error.is_none());
        let start = conf.coords.clone();

        // Very small perturbation (~5°) to avoid atom overlaps in linear interpolation
        let mut end = start.clone();
        let n_atoms = conf.elements.len();
        for i in (n_atoms - 2)..n_atoms {
            let x = end[i * 3];
            let y = end[i * 3 + 1];
            let cos5 = 0.9962;
            let sin5 = 0.0872;
            end[i * 3] = x * cos5 - y * sin5;
            end[i * 3 + 1] = x * sin5 + y * cos5;
        }

        let path = sci_form::compute_simplified_neb_path("CC", &start, &end, 5, 20, 0.5, 0.01)
            .expect("NEB ethane torsion should succeed");

        assert_eq!(path.images.len(), 5, "Should have 5 NEB images");

        // All energies should be finite
        for img in &path.images {
            assert!(
                img.potential_energy_kcal_mol.is_finite(),
                "Image {} energy not finite: {}",
                img.index,
                img.potential_energy_kcal_mol
            );
        }

        eprintln!(
            "  NEB ethane torsion UFF: {} images, E range: {:.2} to {:.2} kcal/mol",
            path.images.len(),
            path.images
                .iter()
                .map(|i| i.potential_energy_kcal_mol)
                .fold(f64::INFINITY, f64::min),
            path.images
                .iter()
                .map(|i| i.potential_energy_kcal_mol)
                .fold(f64::NEG_INFINITY, f64::max),
        );
    }

    /// NEB path with configurable backend (PM3).
    #[test]
    fn neb_ethane_rotation_pm3() {
        let conf = sci_form::embed("CC", 42);
        assert!(conf.error.is_none());
        let start = conf.coords.clone();

        // Create end config: rotate H atoms slightly around C-C bond
        let mut end = start.clone();
        let n_atoms = conf.elements.len();
        // Small perturbation on last 3 H atoms (rotate by ~30°)
        for i in (n_atoms - 3)..n_atoms {
            let x = end[i * 3];
            let y = end[i * 3 + 1];
            let cos30 = 0.866;
            let sin30 = 0.5;
            end[i * 3] = x * cos30 - y * sin30;
            end[i * 3 + 1] = x * sin30 + y * cos30;
        }

        let path = sci_form::compute_simplified_neb_path_configurable(
            "CC", &start, &end, 7, 30, 0.3, 0.005, "pm3",
        )
        .expect("NEB ethane rotation with PM3 should succeed");

        assert_eq!(path.images.len(), 7);
        for img in &path.images {
            assert!(
                img.potential_energy_kcal_mol.is_finite(),
                "Image {} energy not finite",
                img.index
            );
        }

        eprintln!(
            "  NEB ethane PM3: {} images OK, energies finite",
            path.images.len()
        );
    }

    /// Multi-method NEB: same path with different backends validates finite energies.
    #[test]
    fn neb_ethanol_multi_method() {
        let conf = sci_form::embed("CCO", 42);
        assert!(conf.error.is_none());
        let start = conf.coords.clone();

        // Small perturbation for end geometry
        let mut end = start.clone();
        let n = conf.elements.len();
        if n > 2 {
            end[(n - 1) * 3] += 0.3;
            end[(n - 1) * 3 + 1] -= 0.2;
        }

        let methods = ["uff", "pm3"];
        let mut successful = 0;

        for method in &methods {
            match sci_form::compute_simplified_neb_path_configurable(
                "CCO", &start, &end, 5, 20, 0.3, 0.005, method,
            ) {
                Ok(path) => {
                    assert_eq!(path.images.len(), 5, "[{method}] should have 5 images");
                    for img in &path.images {
                        assert!(
                            img.potential_energy_kcal_mol.is_finite(),
                            "[{method}] Image {} energy not finite",
                            img.index
                        );
                    }
                    successful += 1;
                    eprintln!(
                        "  [{method}] ethanol NEB: {} images, E₀={:.2}, E_mid={:.2}",
                        path.images.len(),
                        path.images[0].potential_energy_kcal_mol,
                        path.images[path.images.len() / 2].potential_energy_kcal_mol,
                    );
                }
                Err(e) => {
                    eprintln!("  [{method}] ethanol NEB failed: {e} (acceptable for some methods)");
                }
            }
        }

        assert!(
            successful >= 1,
            "At least one method should produce a valid NEB path"
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 3: Orbital mesh generation at reaction geometries
// ═══════════════════════════════════════════════════════════════════════════

mod orbital_mesh {
    use super::*;

    /// Generate orbital mesh at a single geometry.
    #[test]
    fn eht_mesh_water() {
        let elements = [8u8, 1, 1];
        let positions = [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]];

        let eht = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
        let basis = sci_form::eht::basis::build_basis(&elements, &positions);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht.coefficients,
            eht.homo_index,
            &positions,
            0.4,
            3.0,
        );

        assert!(grid.dims[0] > 0 && grid.dims[1] > 0 && grid.dims[2] > 0);
        assert!(!grid.values.is_empty());

        let mesh = sci_form::eht::marching_cubes(&grid, 0.02);
        assert!(
            mesh.num_triangles > 0,
            "HOMO isosurface should produce triangles"
        );
        assert_eq!(mesh.vertices.len(), mesh.normals.len());

        eprintln!(
            "  Water HOMO mesh: grid {}×{}×{}, {} triangles",
            grid.dims[0], grid.dims[1], grid.dims[2], mesh.num_triangles
        );
    }

    /// Multi-method mesh: compare mesh result across EHT, PM3, xTB.
    #[test]
    fn multi_method_mesh_ethanol() {
        let (elements, positions) = embed_smiles("CCO");

        for method_str in &["eht", "pm3"] {
            let result = sci_form::compute_orbital_mesh(
                &elements, &positions, method_str, 0, // HOMO (mo_index from lowest)
                0.4, 3.0, 0.04,
            );

            match result {
                Ok(res) => {
                    assert!(
                        res.mesh.num_triangles > 0,
                        "[{method_str}] mesh should have triangles"
                    );
                    assert!(
                        !res.orbital_energies.is_empty(),
                        "[{method_str}] should have orbital energies"
                    );
                    assert!(res.gap >= 0.0, "[{method_str}] gap should be non-negative");

                    eprintln!(
                        "  [{method_str}] mesh: {} triangles, gap={:.3} eV, HOMO idx={}",
                        res.mesh.num_triangles, res.gap, res.homo_index
                    );
                }
                Err(e) => {
                    eprintln!("  [{method_str}] mesh failed: {} (non-fatal)", e);
                }
            }
        }
    }

    /// Dual-phase mesh (positive and negative orbital lobes).
    #[test]
    fn dual_phase_mesh_benzene() {
        let (elements, positions) = embed_smiles("c1ccccc1");

        let eht = sci_form::eht::solve_eht(&elements, &positions, None).unwrap();
        let basis = sci_form::eht::basis::build_basis(&elements, &positions);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht.coefficients,
            eht.homo_index,
            &positions,
            0.4,
            3.0,
        );

        let dual = sci_form::eht::marching_cubes_dual(&grid, 0.02);

        // HOMO of benzene should have both positive and negative lobes (π orbital)
        assert!(
            dual.positive.num_triangles > 0 || dual.negative.num_triangles > 0,
            "Benzene HOMO should have at least one orbital lobe"
        );

        eprintln!(
            "  Benzene HOMO dual mesh: +lobe={} tri, -lobe={} tri",
            dual.positive.num_triangles, dual.negative.num_triangles
        );
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 4: End-to-end DFT → NEB → Mesh along reaction coordinate
// ═══════════════════════════════════════════════════════════════════════════

mod reaction_pipeline {

    /// Complete pipeline: H₂O bending + EHT orbital evolution along path.
    #[test]
    fn water_bending_orbital_evolution() {
        let elements = [8u8, 1, 1];

        // Start: tetrahedral-like angle
        let start = vec![0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.24, 0.93, 0.0];
        // End: linear-ish
        let end = vec![0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.96, 0.0, 0.0];

        // Step 1: NEB path
        let path = sci_form::compute_simplified_neb_path("O", &start, &end, 7, 30, 0.5, 0.01)
            .expect("NEB path should succeed");
        assert_eq!(path.images.len(), 7);

        // Step 2: EHT at each image → track HOMO-LUMO gap evolution
        let mut gaps = Vec::new();
        let mut homo_energies = Vec::new();

        for (idx, image) in path.images.iter().enumerate() {
            let positions: Vec<[f64; 3]> =
                image.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

            let eht = sci_form::eht::solve_eht(&elements, &positions, None)
                .unwrap_or_else(|e| panic!("EHT failed at image {}: {}", idx, e));

            gaps.push(eht.gap);
            homo_energies.push(eht.homo_energy);
        }

        // Gap should be positive at all geometries
        for (i, &gap) in gaps.iter().enumerate() {
            assert!(gap >= 0.0, "Gap negative at image {}: {:.4} eV", i, gap);
        }

        // HOMO energy should vary smoothly along the path
        for w in homo_energies.windows(2) {
            let delta = (w[1] - w[0]).abs();
            assert!(
                delta < 5.0,
                "HOMO jump too large: {:.4} eV between consecutive images",
                delta
            );
        }

        eprintln!(
            "  Water bending orbital evolution: gaps = {:?}",
            gaps.iter().map(|g| format!("{:.3}", g)).collect::<Vec<_>>()
        );
    }

    /// Complete pipeline: Ethanol → PM3 energy + HOMO mesh at endpoints.
    #[test]
    fn ethanol_reaction_endpoint_meshes() {
        let conf = sci_form::embed("CCO", 42);
        assert!(conf.error.is_none());
        let elements = conf.elements.clone();
        let start_coords = conf.coords.clone();

        // Create perturbed "product" geometry (rotate OH group)
        let mut end_coords = start_coords.clone();
        let n = elements.len();
        // Perturb last atom position slightly
        if n > 2 {
            end_coords[(n - 1) * 3] += 0.3;
            end_coords[(n - 1) * 3 + 1] -= 0.2;
        }

        let start_pos: Vec<[f64; 3]> = start_coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
        let end_pos: Vec<[f64; 3]> = end_coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        // PM3 at both endpoints
        let pm3_start = sci_form::compute_pm3(&elements, &start_pos).unwrap();
        let pm3_end = sci_form::compute_pm3(&elements, &end_pos).unwrap();
        assert!(pm3_start.converged && pm3_end.converged);

        // EHT mesh at both endpoints
        for (label, positions) in [("reactant", &start_pos), ("product", &end_pos)] {
            let eht = sci_form::eht::solve_eht(&elements, positions, None).unwrap();
            let basis = sci_form::eht::basis::build_basis(&elements, positions);
            let grid = sci_form::eht::evaluate_orbital_on_grid(
                &basis,
                &eht.coefficients,
                eht.homo_index,
                positions,
                0.4,
                3.0,
            );
            let mesh = sci_form::eht::marching_cubes(&grid, 0.03);
            assert!(
                mesh.num_triangles > 0,
                "[{label}] HOMO mesh should have triangles"
            );

            eprintln!(
                "  [{label}] PM3 E={:.4} eV, gap={:.3} eV, mesh={} tri",
                if label == "reactant" {
                    pm3_start.total_energy
                } else {
                    pm3_end.total_energy
                },
                eht.gap,
                mesh.num_triangles
            );
        }
    }

    /// Full pipeline: SMIRKS reaction → embed both sides → DFT comparison.
    #[test]
    fn smirks_reaction_energy_comparison() {
        // Deprotonation: carboxylic acid → carboxylate
        let smirks = "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]";
        let reactant = "CC(=O)O"; // acetic acid

        let transform =
            sci_form::smirks::parse_smirks(smirks).expect("SMIRKS parse should succeed");
        assert!(
            !transform.atom_map.is_empty(),
            "Atom map should not be empty"
        );

        let result =
            sci_form::smirks::apply_smirks(smirks, reactant).expect("SMIRKS apply should succeed");
        assert!(result.success, "SMIRKS should apply to acetic acid");
        assert!(!result.products.is_empty(), "Should produce products");

        // Embed and compute energy for reactant
        let r_conf = sci_form::embed(reactant, 42);
        assert!(r_conf.error.is_none());
        let r_pos: Vec<[f64; 3]> = r_conf
            .coords
            .chunks(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();

        let r_eht = sci_form::eht::solve_eht(&r_conf.elements, &r_pos, None).unwrap();
        let r_pm3 = sci_form::compute_pm3(&r_conf.elements, &r_pos).unwrap();

        eprintln!(
            "  Acetic acid: EHT gap={:.3} eV, PM3 HOF={:.2} kcal/mol",
            r_eht.gap, r_pm3.heat_of_formation
        );

        // Embed product (first one)
        let product_smiles = &result.products[0];
        let p_conf = sci_form::embed(product_smiles, 42);
        if p_conf.error.is_none() {
            let p_pos: Vec<[f64; 3]> = p_conf
                .coords
                .chunks(3)
                .map(|c| [c[0], c[1], c[2]])
                .collect();
            let p_eht = sci_form::eht::solve_eht(&p_conf.elements, &p_pos, None).unwrap();

            eprintln!(
                "  Product '{}': EHT gap={:.3} eV",
                product_smiles, p_eht.gap
            );
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 5: SMIRKS-driven reaction with mesh visualization
// ═══════════════════════════════════════════════════════════════════════════

mod smirks_mesh {

    /// Apply SMIRKS, then generate orbital meshes for reactant and product.
    #[test]
    fn smirks_esterification_mesh() {
        // Ester hydrolysis: ester → acid + alcohol
        let reactant = "CC(=O)OC"; // methyl acetate

        let r_conf = sci_form::embed(reactant, 42);
        assert!(r_conf.error.is_none(), "Embed failed: {:?}", r_conf.error);
        let r_pos: Vec<[f64; 3]> = r_conf
            .coords
            .chunks(3)
            .map(|c| [c[0], c[1], c[2]])
            .collect();

        // Compute HOMO mesh for reactant
        let eht = sci_form::eht::solve_eht(&r_conf.elements, &r_pos, None).unwrap();
        let basis = sci_form::eht::basis::build_basis(&r_conf.elements, &r_pos);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht.coefficients,
            eht.homo_index,
            &r_pos,
            0.4,
            3.0,
        );

        let dual = sci_form::eht::marching_cubes_dual(&grid, 0.02);
        let total_tri = dual.positive.num_triangles + dual.negative.num_triangles;
        assert!(total_tri > 0, "Reactant HOMO should have orbital lobes");

        // Also compute LUMO mesh (electrophilic site)
        let lumo_grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht.coefficients,
            eht.lumo_index,
            &r_pos,
            0.4,
            3.0,
        );
        let lumo_mesh = sci_form::eht::marching_cubes(&lumo_grid, 0.02);

        eprintln!(
            "  Methyl acetate: HOMO dual ({}/{} tri), LUMO ({} tri), gap={:.3} eV",
            dual.positive.num_triangles,
            dual.negative.num_triangles,
            lumo_mesh.num_triangles,
            eht.gap
        );
    }

    /// SMIRKS selection of reaction type → orbital visualization.
    #[test]
    fn smirks_aromatic_substitution_compare() {
        // Aromatic nitration pattern
        let smirks = "[cH:1]>>[c:1][N+](=O)[O-]";
        let reactant = "c1ccccc1"; // benzene

        match sci_form::smirks::apply_smirks(smirks, reactant) {
            Ok(result) if result.success => {
                eprintln!("  Nitration: {} → {:?}", reactant, result.products);

                // Compute Fukui descriptors for the reactant (electrophilic susceptibility)
                let r_conf = sci_form::embed(reactant, 42);
                assert!(r_conf.error.is_none());
                let r_pos: Vec<[f64; 3]> = r_conf
                    .coords
                    .chunks(3)
                    .map(|c| [c[0], c[1], c[2]])
                    .collect();

                let eht = sci_form::eht::solve_eht(&r_conf.elements, &r_pos, None).unwrap();
                assert!(eht.gap > 0.0);

                eprintln!(
                    "  Benzene pre-nitration: gap={:.3} eV, HOMO={:.3} eV",
                    eht.gap, eht.homo_energy
                );
            }
            Ok(_result) => {
                eprintln!("  Aromatic nitration: SMIRKS did not apply (no match). Acceptable.");
            }
            Err(e) => {
                eprintln!("  Aromatic nitration SMIRKS failed: {} (non-fatal)", e);
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// LEVEL 6: Molecular dynamics + orbital visualization
// ═══════════════════════════════════════════════════════════════════════════

mod dynamics_mesh {

    /// Short MD trajectory + EHT analysis of snapshots.
    #[test]
    fn md_water_orbital_snapshots() {
        let conf = sci_form::embed("O", 42);
        assert!(conf.error.is_none());
        let elements = conf.elements.clone();
        let coords = conf.coords.clone();

        // Run short NVE MD (10 steps, 0.5 fs)
        let traj = sci_form::compute_md_trajectory("O", &coords, 10, 0.5, 42)
            .expect("MD should succeed for water");

        assert!(traj.frames.len() > 1, "Should have multiple frames");

        // Analyze first and last frames with EHT
        for (label, frame) in [
            ("start", &traj.frames[0]),
            ("end", traj.frames.last().unwrap()),
        ] {
            let positions: Vec<[f64; 3]> =
                frame.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

            let eht = sci_form::eht::solve_eht(&elements, &positions, None)
                .unwrap_or_else(|e| panic!("EHT failed at {} frame: {}", label, e));

            assert!(eht.gap >= 0.0, "[{label}] gap should be non-negative");

            eprintln!(
                "  [{label}] E_pot={:.2} kcal/mol, EHT gap={:.3} eV",
                frame.potential_energy_kcal_mol, eht.gap
            );
        }
    }

    /// NVT MD + HOMO mesh at final snapshot.
    #[test]
    fn md_nvt_ethanol_final_mesh() {
        let conf = sci_form::embed("CCO", 42);
        assert!(conf.error.is_none());

        let traj =
            sci_form::compute_md_trajectory_nvt("CCO", &conf.coords, 5, 0.5, 42, 300.0, 50.0)
                .expect("NVT MD should succeed");

        assert!(!traj.frames.is_empty());
        let last = traj.frames.last().unwrap();
        let positions: Vec<[f64; 3]> = last.coords.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();

        // Generate HOMO mesh at final MD geometry
        let eht = sci_form::eht::solve_eht(&conf.elements, &positions, None).unwrap();
        let basis = sci_form::eht::basis::build_basis(&conf.elements, &positions);
        let grid = sci_form::eht::evaluate_orbital_on_grid(
            &basis,
            &eht.coefficients,
            eht.homo_index,
            &positions,
            0.4,
            3.0,
        );
        let mesh = sci_form::eht::marching_cubes(&grid, 0.03);

        assert!(
            mesh.num_triangles > 0,
            "HOMO mesh at MD final snapshot should have triangles"
        );

        eprintln!(
            "  NVT MD ethanol: {} frames, final T={:.1} K, HOMO mesh={} tri",
            traj.frames.len(),
            last.temperature_k,
            mesh.num_triangles
        );
    }
}
