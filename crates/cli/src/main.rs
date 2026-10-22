mod args;
mod calc_cmds;
mod embed_cmds;
mod experimental_cmds;
mod format;

use args::{Cli, Commands};
use clap::Parser;

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Embed {
            smiles,
            seed,
            format,
        } => embed_cmds::cmd_embed(&smiles, seed, &format),
        Commands::Batch {
            input,
            output,
            seed,
            threads,
            format,
        } => embed_cmds::cmd_batch(input, output, seed, threads, &format),
        Commands::Parse { smiles } => embed_cmds::cmd_parse(&smiles),
        Commands::Info => embed_cmds::cmd_info(),
        Commands::Charges { smiles } => embed_cmds::cmd_charges(&smiles),
        Commands::Uff { smiles, coords } => embed_cmds::cmd_uff(&smiles, &coords),
        Commands::SimplifiedNebPath {
            smiles,
            start_coords,
            end_coords,
            n_images,
            n_iter,
            spring_k,
            step_size,
            method,
        } => experimental_cmds::cmd_simplified_neb_path(
            &smiles,
            &start_coords,
            &end_coords,
            n_images,
            n_iter,
            spring_k,
            step_size,
            &method,
        ),
        Commands::NebEnergy {
            smiles,
            coords,
            method,
        } => experimental_cmds::cmd_neb_energy(&smiles, &coords, &method),
        Commands::NebGradient {
            smiles,
            coords,
            method,
        } => experimental_cmds::cmd_neb_gradient(&smiles, &coords, &method),
        Commands::Rmsd { coords, reference } => embed_cmds::cmd_rmsd(&coords, &reference),
        Commands::Eht {
            elements,
            coords,
            k,
        } => calc_cmds::cmd_eht(&elements, &coords, k),
        Commands::Sasa {
            elements,
            coords,
            probe_radius,
        } => calc_cmds::cmd_sasa(&elements, &coords, probe_radius),
        Commands::Population { elements, coords } => calc_cmds::cmd_population(&elements, &coords),
        Commands::Dipole { elements, coords } => calc_cmds::cmd_dipole(&elements, &coords),
        Commands::Esp {
            elements,
            coords,
            spacing,
            padding,
        } => calc_cmds::cmd_esp(&elements, &coords, spacing, padding),
        Commands::Dos {
            elements,
            coords,
            sigma,
            e_min,
            e_max,
            n_points,
        } => calc_cmds::cmd_dos(&elements, &coords, sigma, e_min, e_max, n_points),
        Commands::Cell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => calc_cmds::cmd_cell(a, b, c, alpha, beta, gamma),
        Commands::Assemble {
            topology,
            a,
            metal,
            geometry,
            supercell,
        } => calc_cmds::cmd_assemble(&topology, a, metal, &geometry, supercell),
        Commands::Ani { elements, coords } => calc_cmds::cmd_ani(&elements, &coords),
        Commands::Pm3 { elements, coords } => calc_cmds::cmd_pm3(&elements, &coords),
        Commands::Xtb { elements, coords } => calc_cmds::cmd_xtb(&elements, &coords),
        Commands::Gfn1 { elements, coords } => calc_cmds::cmd_gfn1(&elements, &coords),
        Commands::Gfn2 { elements, coords } => calc_cmds::cmd_gfn2(&elements, &coords),
        Commands::Hf3c { elements, coords } => calc_cmds::cmd_hf3c(&elements, &coords),
        Commands::Stereo { smiles, coords } => calc_cmds::cmd_stereo(&smiles, coords.as_deref()),
        Commands::Solvation {
            elements,
            coords,
            charges,
            probe_radius,
        } => calc_cmds::cmd_solvation(&elements, &coords, charges.as_deref(), probe_radius),
        Commands::Sssr { smiles } => calc_cmds::cmd_sssr(&smiles),
        Commands::Ecfp {
            smiles,
            radius,
            n_bits,
        } => calc_cmds::cmd_ecfp(&smiles, radius, n_bits),
        Commands::Tanimoto {
            smiles1,
            smiles2,
            radius,
            n_bits,
        } => calc_cmds::cmd_tanimoto(&smiles1, &smiles2, radius, n_bits),

        // Cross-surface parity commands
        Commands::EmbedDiverse {
            smiles,
            n_conformers,
            rmsd_cutoff,
            base_seed,
        } => embed_cmds::cmd_embed_diverse(&smiles, n_conformers, rmsd_cutoff, base_seed),
        Commands::CanonicalSmiles { smiles } => embed_cmds::cmd_canonical_smiles(&smiles),
        Commands::EcfpBatch {
            smiles_list,
            radius,
            n_bits,
        } => calc_cmds::cmd_ecfp_batch(&smiles_list, radius, n_bits),
        Commands::ChargesConfigured {
            smiles,
            max_iter,
            initial_damping,
            convergence_threshold,
        } => calc_cmds::cmd_charges_configured(
            &smiles,
            max_iter,
            initial_damping,
            convergence_threshold,
        ),
        Commands::EnsembleJCouplings {
            smiles,
            conformer_coords,
            energies_kcal,
            temperature_k,
        } => calc_cmds::cmd_ensemble_j_couplings(
            &smiles,
            &conformer_coords,
            &energies_kcal,
            temperature_k,
        ),
        Commands::IrSpectrumBroadened {
            elements,
            coords,
            method,
            gamma,
            wn_min,
            wn_max,
            n_points,
            broadening,
            step_size,
            smiles,
        } => calc_cmds::cmd_ir_spectrum_broadened(
            &elements,
            &coords,
            &method,
            gamma,
            wn_min,
            wn_max,
            n_points,
            &broadening,
            step_size,
            smiles.as_deref(),
        ),
        Commands::Cisd {
            elements,
            coords,
            n_states,
        } => calc_cmds::cmd_cisd(&elements, &coords, n_states),
        Commands::DosMultimethod {
            elements,
            coords,
            method,
            sigma,
            e_min,
            e_max,
            n_points,
        } => calc_cmds::cmd_dos_multimethod(
            &elements, &coords, &method, sigma, e_min, e_max, n_points,
        ),
        Commands::EhtGradient { elements, coords } => {
            calc_cmds::cmd_eht_gradient(&elements, &coords)
        }
        Commands::OptimizeGeometryEht {
            elements,
            coords,
            max_steps,
            grad_threshold,
            energy_threshold,
            step_size,
        } => calc_cmds::cmd_optimize_geometry_eht(
            &elements,
            &coords,
            max_steps,
            grad_threshold,
            energy_threshold,
            step_size,
        ),
        Commands::AevsTm { elements, coords } => calc_cmds::cmd_aevs_tm(&elements, &coords),
        Commands::Whim {
            elements,
            coords,
            weighting,
        } => calc_cmds::cmd_whim(&elements, &coords, &weighting),
        Commands::Rdf {
            elements,
            coords,
            charges,
            beta,
            r_max,
            n_points,
        } => calc_cmds::cmd_rdf(
            &elements,
            &coords,
            charges.as_deref(),
            beta,
            r_max,
            n_points,
        ),
        Commands::Getaway {
            elements,
            coords,
            bonds,
            max_lag,
        } => calc_cmds::cmd_getaway(&elements, &coords, bonds.as_deref(), max_lag),
        Commands::TrainRandomForest {
            features,
            targets,
            n_trees,
            max_depth,
            min_samples_split,
            min_samples_leaf,
            sample_fraction,
            seed,
        } => calc_cmds::cmd_train_random_forest(
            &features,
            &targets,
            n_trees,
            max_depth,
            min_samples_split,
            min_samples_leaf,
            sample_fraction,
            seed,
        ),
        Commands::TrainGradientBoosting {
            features,
            targets,
            n_estimators,
            learning_rate,
            max_depth,
            min_samples_split,
            min_samples_leaf,
            subsample,
            seed,
        } => calc_cmds::cmd_train_gradient_boosting(
            &features,
            &targets,
            n_estimators,
            learning_rate,
            max_depth,
            min_samples_split,
            min_samples_leaf,
            subsample,
            seed,
        ),
        Commands::AssignPeaks {
            frequencies,
            intensities,
            elements,
            mode_displacements,
        } => calc_cmds::cmd_assign_peaks(
            &frequencies,
            &intensities,
            &elements,
            mode_displacements.as_deref(),
        ),
        Commands::SpaceGroup { number, symbol } => {
            calc_cmds::cmd_space_group(number, symbol.as_deref())
        }
        Commands::AllSpaceGroups => calc_cmds::cmd_all_space_groups(),
        Commands::PeriodicMolecule {
            elements,
            frac_coords,
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            bond_tolerance,
        } => calc_cmds::cmd_periodic_molecule(
            &elements,
            &frac_coords,
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
            bond_tolerance,
        ),
        Commands::Hapticity {
            elements,
            coords,
            bonds,
        } => calc_cmds::cmd_hapticity(&elements, &coords, bonds.as_deref()),
        Commands::OptimizeFramework {
            smiles,
            coords,
            method,
            max_iter,
            force_tol,
            energy_tol,
            max_step,
            fixed_atoms,
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => {
            let cell = match (a, b, c) {
                (Some(a), Some(b), Some(c)) => Some([a, b, c, alpha, beta, gamma]),
                _ => None,
            };
            calc_cmds::cmd_optimize_framework(
                &smiles,
                &coords,
                &method,
                max_iter,
                force_tol,
                energy_tol,
                max_step,
                fixed_atoms.as_deref(),
                cell,
            )
        }
        Commands::FilterDiverse {
            conformers,
            rmsd_cutoff,
        } => embed_cmds::cmd_filter_diverse(&conformers, rmsd_cutoff),

        // Experimental commands
        Commands::Eeq {
            elements,
            coords,
            total_charge,
        } => experimental_cmds::cmd_eeq(&elements, &coords, total_charge),
        #[cfg(feature = "experimental-alpb")]
        Commands::Alpb {
            elements,
            coords,
            charges,
            dielectric,
        } => experimental_cmds::cmd_alpb(&elements, &coords, charges.as_deref(), dielectric),
        #[cfg(feature = "experimental-d4")]
        Commands::D4 {
            elements,
            coords,
            three_body,
        } => experimental_cmds::cmd_d4(&elements, &coords, three_body),
        #[cfg(feature = "experimental-cpm")]
        Commands::Cpm {
            elements,
            coords,
            mu_ev,
            dielectric,
        } => experimental_cmds::cmd_cpm(&elements, &coords, mu_ev, dielectric),

        #[cfg(feature = "alpha-edl")]
        Commands::EdlProfile {
            potential,
            model,
            ionic_strength,
            temperature,
            n_points,
            extent,
        } => experimental_cmds::cmd_edl_profile(
            potential,
            &model,
            ionic_strength,
            temperature,
            n_points,
            extent,
        ),
        #[cfg(feature = "alpha-edl")]
        Commands::EdlCapacitanceScan {
            v_min,
            v_max,
            n_points,
            ionic_strength,
            model,
        } => experimental_cmds::cmd_edl_capacitance_scan(
            v_min,
            v_max,
            n_points,
            ionic_strength,
            &model,
        ),
        #[cfg(feature = "alpha-edl")]
        Commands::EdlHelmholtz {
            potential,
            ionic_strength,
            temperature,
        } => experimental_cmds::cmd_edl_helmholtz(potential, ionic_strength, temperature),
        #[cfg(feature = "alpha-edl")]
        Commands::EdlGouyChapman {
            potential,
            ionic_strength,
            temperature,
        } => experimental_cmds::cmd_edl_gouy_chapman(potential, ionic_strength, temperature),
        #[cfg(feature = "alpha-edl")]
        Commands::EdlGcs {
            potential,
            ionic_strength,
            temperature,
            stern_thickness,
        } => {
            experimental_cmds::cmd_edl_gcs(potential, ionic_strength, temperature, stern_thickness)
        }
        #[cfg(feature = "alpha-kinetics")]
        Commands::HtstRate {
            activation_free_energy,
            reaction_free_energy,
            temperature,
        } => experimental_cmds::cmd_htst_rate(
            activation_free_energy,
            reaction_free_energy,
            temperature,
        ),
        #[cfg(feature = "alpha-periodic-linear")]
        Commands::Kmesh {
            nx,
            ny,
            nz,
            centering,
        } => experimental_cmds::cmd_kmesh([nx, ny, nz], &centering),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
        Commands::EdlChart {
            potential,
            ionic_strength,
        } => experimental_cmds::cmd_edl_chart(potential, ionic_strength),
        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-edl"))]
        Commands::CapacitanceChart {
            v_min,
            v_max,
            n_points,
            ionic_strength,
            model,
        } => {
            experimental_cmds::cmd_capacitance_chart(v_min, v_max, n_points, ionic_strength, &model)
        }

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
        Commands::ArrheniusChart {
            temperatures,
            rates,
        } => experimental_cmds::cmd_arrhenius_chart(&temperatures, &rates),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
        Commands::ThermalPropChart {
            activation_free_energy,
            reaction_free_energy,
            temperatures,
        } => experimental_cmds::cmd_thermal_prop_chart(
            activation_free_energy,
            reaction_free_energy,
            &temperatures,
        ),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-kinetics"))]
        Commands::PhasePortraitChart { network, n_steps } => {
            experimental_cmds::cmd_phase_portrait_chart(network.as_deref(), n_steps)
        }

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
        Commands::BandStructureChart {
            elements,
            coords,
            lattice,
            n_kpoints_per_segment,
        } => experimental_cmds::cmd_band_structure_chart(
            &elements,
            &coords,
            &lattice,
            n_kpoints_per_segment,
        ),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
        Commands::DosChart {
            energies,
            dos_values,
        } => experimental_cmds::cmd_dos_chart(&energies, &dos_values),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
        Commands::KpointPathChart {
            n_points_per_segment,
        } => experimental_cmds::cmd_kpoint_path_chart(n_points_per_segment),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-periodic-linear"))]
        Commands::FermiSurfaceChart {
            blocks,
            fermi_energy_ev,
            n_kx,
            n_ky,
        } => experimental_cmds::cmd_fermi_surface_chart(&blocks, fermi_energy_ev, n_kx, n_ky),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
        Commands::TrajectoryChart {
            smiles,
            reactant_coords,
            product_coords,
            method,
        } => experimental_cmds::cmd_trajectory_chart(
            &smiles,
            &reactant_coords,
            &product_coords,
            &method,
        ),

        #[cfg(all(feature = "alpha-render-bridge", feature = "alpha-gsm"))]
        Commands::ReactionCoordinateChart {
            smiles,
            reactant_coords,
            product_coords,
            method,
        } => experimental_cmds::cmd_reaction_coordinate_chart(
            &smiles,
            &reactant_coords,
            &product_coords,
            &method,
        ),

        // ─── PERIODIC LINEAR ───
        #[cfg(feature = "alpha-periodic-linear")]
        Commands::ComputePeriodicDos {
            hamiltonian_blocks,
            nx,
            ny,
            nz,
            centering,
            order,
            e_min,
            e_max,
            n_points,
            fermi_energy_ev,
        } => experimental_cmds::cmd_compute_periodic_dos(
            &hamiltonian_blocks,
            [nx, ny, nz],
            &centering,
            order,
            e_min,
            e_max,
            n_points,
            fermi_energy_ev,
        ),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::SolvePeriodicRandnla {
            hamiltonian_blocks,
            overlap_blocks,
            n_electrons,
            nx,
            ny,
            nz,
            centering,
            sketch_size,
        } => experimental_cmds::cmd_solve_periodic_randnla(
            &hamiltonian_blocks,
            &overlap_blocks,
            n_electrons,
            [nx, ny, nz],
            &centering,
            sketch_size,
        ),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::BlochPhase {
            k_x,
            k_y,
            k_z,
            t_x,
            t_y,
            t_z,
        } => experimental_cmds::cmd_bloch_phase(k_x, k_y, k_z, t_x, t_y, t_z),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::BuildBlochHamiltonian {
            blocks,
            k_x,
            k_y,
            k_z,
        } => experimental_cmds::cmd_build_bloch_hamiltonian(&blocks, k_x, k_y, k_z),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::AssemblePeriodicOperators {
            hamiltonian_blocks,
            overlap_blocks,
            nx,
            ny,
            nz,
            centering,
        } => experimental_cmds::cmd_assemble_periodic_operators(
            &hamiltonian_blocks,
            overlap_blocks.as_deref(),
            [nx, ny, nz],
            &centering,
        ),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::ValidateElectronCount {
            n_electrons,
            n_bands,
            n_kpoints,
        } => experimental_cmds::cmd_validate_electron_count(n_electrons, n_bands, n_kpoints),

        #[cfg(feature = "alpha-gsm")]
        Commands::GsmBackendPlan { smiles } => experimental_cmds::cmd_gsm_backend_plan(&smiles),

        #[cfg(feature = "alpha-gsm")]
        Commands::GsmCompareBackends {
            smiles,
            coords,
            methods,
        } => experimental_cmds::cmd_gsm_compare_backends(&smiles, &coords, methods.as_deref()),

        // ─── KINETICS ──────────
        #[cfg(feature = "alpha-kinetics")]
        Commands::ExtractKineticsDiagnostics { network } => {
            experimental_cmds::cmd_extract_kinetics_diagnostics(network.as_deref())
        }

        #[cfg(feature = "alpha-kinetics")]
        Commands::SolveMicrokineticNetwork { network, n_steps } => {
            experimental_cmds::cmd_solve_microkinetic_network(network.as_deref(), n_steps)
        }

        #[cfg(feature = "alpha-kinetics")]
        Commands::SolveMicrokineticSteadyState {
            network,
            n_steps,
            tol,
        } => {
            experimental_cmds::cmd_solve_microkinetic_steady_state(network.as_deref(), n_steps, tol)
        }

        #[cfg(feature = "alpha-kinetics")]
        Commands::AnalyzeGsmMbhHtstStep {
            smiles,
            reactant_coords,
            product_coords,
            method,
            step_id,
            temperature_k,
            pressure_bar,
            n_nodes,
            mbh_fd_step,
        } => experimental_cmds::cmd_analyze_gsm_mbh_htst_step(
            &smiles,
            &reactant_coords,
            &product_coords,
            &method,
            &step_id,
            temperature_k,
            pressure_bar,
            n_nodes,
            mbh_fd_step,
        ),

        // Render Bridge commands
        Commands::ChartToJson { json_input } => {
            experimental_cmds::cmd_chart_to_json(json_input.unwrap_or_default())
        }
        Commands::ChartFromJson { json_input } => {
            experimental_cmds::cmd_chart_from_json(json_input.unwrap_or_default())
        }
        Commands::ValidateChartSchema { json_input } => {
            experimental_cmds::cmd_validate_chart_schema(json_input.unwrap_or_default())
        }

        #[cfg(feature = "alpha-reaction-dynamics")]
        Commands::ReactionDynamics3d {
            reactant_smiles,
            product_smiles,
            method,
            n_images,
            n_approach,
            n_departure,
            smirks,
            seed,
        } => {
            let reactants: Vec<String> =
                serde_json::from_str(&reactant_smiles).unwrap_or_else(|e| {
                    eprintln!("bad reactant_smiles: {e}");
                    std::process::exit(1);
                });
            let products: Vec<String> = serde_json::from_str(&product_smiles).unwrap_or_else(|e| {
                eprintln!("bad product_smiles: {e}");
                std::process::exit(1);
            });
            let config = sci_form::alpha::reaction_dynamics::ReactionDynamics3DConfig {
                method,
                n_images,
                n_approach_frames: n_approach,
                n_departure_frames: n_departure,
                smirks,
                seed,
                ..Default::default()
            };
            let r_refs: Vec<&str> = reactants.iter().map(String::as_str).collect();
            let p_refs: Vec<&str> = products.iter().map(String::as_str).collect();
            match sci_form::alpha::reaction_dynamics::compute_reaction_dynamics_3d(
                &r_refs, &p_refs, &config,
            ) {
                Ok(result) => println!("{}", serde_json::to_string_pretty(&result).unwrap()),
                Err(e) => {
                    eprintln!("Error: {e}");
                    std::process::exit(1);
                }
            }
        }
    }
}
