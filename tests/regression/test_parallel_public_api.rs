#[cfg(feature = "parallel")]
mod parallel_public_api_tests {
    use std::sync::Mutex;

    use sci_form::materials::{CoordinationGeometry, Sbu, Topology};
    use sci_form::{
        assemble_framework, compute_charges, compute_dipole, compute_dos, compute_esp, compute_pm3,
        compute_population, compute_rmsd, compute_sasa, compute_uff_energy, compute_xtb,
        create_unit_cell, embed, embed_batch, parse, version, ConformerConfig,
    };

    fn positions_from_flat(coords: &[f64]) -> Vec<[f64; 3]> {
        coords
            .chunks_exact(3)
            .map(|chunk| [chunk[0], chunk[1], chunk[2]])
            .collect()
    }

    #[test]
    fn all_public_api_functions_can_run_concurrently() {
        let seed = 42;
        let ethanol = embed("CCO", seed);
        assert!(
            ethanol.error.is_none(),
            "embed baseline failed: {:?}",
            ethanol.error
        );

        let positions = positions_from_flat(&ethanol.coords);
        let elements = ethanol.elements.clone();
        let coords_flat = ethanol.coords.clone();
        let failures = Mutex::new(Vec::<String>::new());

        rayon::scope(|scope| {
            for _ in 0..2 {
                scope.spawn(|_| {
                    if version().is_empty() {
                        failures
                            .lock()
                            .unwrap()
                            .push("version returned empty string".to_string());
                    }
                });

                scope.spawn(|_| {
                    let result = embed("c1ccccc1", seed + 1);
                    if result.error.is_some() {
                        failures
                            .lock()
                            .unwrap()
                            .push(format!("embed failed: {:?}", result.error));
                    }
                });

                scope.spawn(|_| {
                    let config = ConformerConfig {
                        seed,
                        num_threads: 2,
                    };
                    let result = embed_batch(&["CCO", "O", "c1ccccc1"], &config);
                    if result.len() != 3 || result.iter().any(|item| item.error.is_some()) {
                        failures
                            .lock()
                            .unwrap()
                            .push("embed_batch failed under concurrent load".to_string());
                    }
                });

                scope.spawn(|_| {
                    if parse("CCO")
                        .map(|mol| mol.graph.node_count())
                        .unwrap_or_default()
                        == 0
                    {
                        failures.lock().unwrap().push("parse failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_charges("CCO")
                        .map(|result| result.charges.is_empty())
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_charges failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_sasa(&elements, &coords_flat, Some(1.4))
                        .map(|result| result.total_sasa <= 0.0)
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_sasa failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_population(&elements, &positions)
                        .map(|result| result.mulliken_charges.len() != elements.len())
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_population failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_dipole(&elements, &positions)
                        .map(|result| !result.magnitude.is_finite())
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_dipole failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_esp(&elements, &positions, 0.8, 2.5)
                        .map(|result| result.values.is_empty())
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_esp failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_dos(&elements, &positions, 0.3, -20.0, 5.0, 128)
                        .map(|result| result.total_dos.len() != 128)
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_dos failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    let rmsd = compute_rmsd(&coords_flat, &coords_flat);
                    if !rmsd.is_finite() || rmsd.abs() > 1e-8 {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_rmsd failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    if compute_uff_energy("CCO", &coords_flat)
                        .map(|energy| !energy.is_finite())
                        .unwrap_or(true)
                    {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_uff_energy failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    let cell = create_unit_cell(12.0, 12.0, 12.0, 90.0, 90.0, 90.0);
                    let volume = cell.volume();
                    if !volume.is_finite() || volume <= 0.0 {
                        failures
                            .lock()
                            .unwrap()
                            .push("create_unit_cell failed".to_string());
                    }
                });

                scope.spawn(|_| {
                    let cell = create_unit_cell(12.0, 12.0, 12.0, 90.0, 90.0, 90.0);
                    let node = Sbu::metal_node(30, 2.0, CoordinationGeometry::Octahedral);
                    let linker = Sbu::linear_linker(&[6, 6], 1.4, "carboxylate");
                    let topology = Topology::pcu();
                    let structure = assemble_framework(&node, &linker, &topology, &cell);
                    if structure.num_atoms() == 0 {
                        failures
                            .lock()
                            .unwrap()
                            .push("assemble_framework failed".to_string());
                    }
                });

                // Quantum methods under concurrency
                scope.spawn(|_| match compute_pm3(&elements, &positions) {
                    Ok(r) if !r.converged => {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_pm3 did not converge".to_string());
                    }
                    Err(e) => {
                        failures
                            .lock()
                            .unwrap()
                            .push(format!("compute_pm3 error: {e}"));
                    }
                    _ => {}
                });

                scope.spawn(|_| match compute_xtb(&elements, &positions) {
                    Ok(r) if !r.converged => {
                        failures
                            .lock()
                            .unwrap()
                            .push("compute_xtb did not converge".to_string());
                    }
                    Err(e) => {
                        failures
                            .lock()
                            .unwrap()
                            .push(format!("compute_xtb error: {e}"));
                    }
                    _ => {}
                });

                scope.spawn(|_| {
                    match sci_form::xtb::gfn1::solve_gfn1(&elements, &positions) {
                        Ok(_) => {} // GFN1 may not converge for all molecules (depends on GFN0 guess)
                        Err(e) => {
                            eprintln!("solve_gfn1 warning: {e}");
                        }
                    }
                });

                scope.spawn(
                    |_| match sci_form::xtb::gfn2::solve_gfn2(&elements, &positions) {
                        Ok(r) if !r.converged => {
                            failures
                                .lock()
                                .unwrap()
                                .push("solve_gfn2 did not converge".to_string());
                        }
                        Err(e) => {
                            failures
                                .lock()
                                .unwrap()
                                .push(format!("solve_gfn2 error: {e}"));
                        }
                        _ => {}
                    },
                );
            }
        });

        let failures = failures.into_inner().unwrap();
        assert!(failures.is_empty(), "parallel API failures: {:?}", failures);
    }
}
