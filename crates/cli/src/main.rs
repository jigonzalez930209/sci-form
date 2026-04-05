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
        Commands::Stereo { smiles, coords } => calc_cmds::cmd_stereo(&smiles, &coords),
        Commands::Solvation {
            elements,
            coords,
            charges,
            probe_radius,
        } => calc_cmds::cmd_solvation(&elements, &coords, &charges, probe_radius),
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

        // Experimental commands
        #[cfg(feature = "experimental-eeq")]
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
        } => experimental_cmds::cmd_alpb(&elements, &coords, &charges, dielectric),
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

        #[cfg(feature = "alpha-render-bridge")]
        Commands::EdlChart {
            potential,
            ionic_strength,
        } => experimental_cmds::cmd_edl_chart(potential, ionic_strength),
        #[cfg(feature = "alpha-render-bridge")]
        Commands::CapacitanceChart {
            v_min,
            v_max,
            n_points,
            ionic_strength,
            model,
        } => {
            experimental_cmds::cmd_capacitance_chart(v_min, v_max, n_points, ionic_strength, &model)
        }

        #[cfg(feature = "alpha-render-bridge")]
        Commands::ArrheniusChart {
            temperatures,
            rates,
        } => experimental_cmds::cmd_arrhenius_chart(&temperatures, &rates),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::BandStructureChart {
            k_points,
            n_bands,
            e_min,
            e_max,
        } => experimental_cmds::cmd_band_structure_chart(k_points, n_bands, e_min, e_max),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::DosChart {
            energies,
            dos_values,
        } => experimental_cmds::cmd_dos_chart(&energies, &dos_values),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::TrajectoryChart => experimental_cmds::cmd_trajectory_chart(),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::KpointPathChart => experimental_cmds::cmd_kpoint_path_chart(),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::FermiSurfaceChart => experimental_cmds::cmd_fermi_surface_chart(),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::PhasePortraitChart => experimental_cmds::cmd_phase_portrait_chart(),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::ReactionCoordinateChart => experimental_cmds::cmd_reaction_coordinate_chart(),

        #[cfg(feature = "alpha-render-bridge")]
        Commands::ThermalPropChart => experimental_cmds::cmd_thermal_prop_chart(),

        // ─── PERIODIC LINEAR ───
        #[cfg(feature = "alpha-periodic-linear")]
        Commands::ComputePeriodicDos { n_kpoints, order } => {
            experimental_cmds::cmd_compute_periodic_dos(n_kpoints, order)
        }

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::SolvePeriodicRandnla { n_kpoints } => {
            experimental_cmds::cmd_solve_periodic_randnla(n_kpoints)
        }

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::BlochPhase { k_x, k_y, k_z } => experimental_cmds::cmd_bloch_phase(k_x, k_y, k_z),

        #[cfg(feature = "alpha-periodic-linear")]
        Commands::BuildBlochHamiltonian { n_basis } => {
            experimental_cmds::cmd_build_bloch_hamiltonian(n_basis)
        }

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
        } => experimental_cmds::cmd_gsm_compare_backends(&smiles, &coords, &methods),

        // ─── KINETICS ──────────
        #[cfg(feature = "alpha-kinetics")]
        Commands::ExtractKineticsDiagnostics => {
            experimental_cmds::cmd_extract_kinetics_diagnostics()
        }

        #[cfg(feature = "alpha-kinetics")]
        Commands::SolveMicrokineticNetwork { n_steps } => {
            experimental_cmds::cmd_solve_microkinetic_network(n_steps)
        }

        #[cfg(feature = "alpha-kinetics")]
        Commands::SolveMicrokineticSteadyState { n_steps } => {
            experimental_cmds::cmd_solve_microkinetic_steady_state(n_steps)
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
        Commands::ChartToJson => experimental_cmds::cmd_chart_to_json(),
        Commands::ChartFromJson { json_input } => {
            experimental_cmds::cmd_chart_from_json(json_input.unwrap_or_default())
        }
        Commands::ValidateChartSchema { json_input } => {
            experimental_cmds::cmd_validate_chart_schema(json_input.unwrap_or_default())
        }

        // Auxiliary commands
        Commands::ValidateExperimentalResult => {
            experimental_cmds::cmd_validate_experimental_result()
        }
        Commands::MergeExperimentResults => experimental_cmds::cmd_merge_experiment_results(),
        Commands::ExperimentalResultToSdf => experimental_cmds::cmd_experimental_result_to_sdf(),
        Commands::ExperimentalResultToJson => experimental_cmds::cmd_experimental_result_to_json(),
        Commands::BenchmarkFunction => experimental_cmds::cmd_benchmark_function(),
        Commands::TraceFunctionCalls => experimental_cmds::cmd_trace_function_calls(),
        Commands::ReportGenerator => experimental_cmds::cmd_report_generator(),
        Commands::PipelineValidator => experimental_cmds::cmd_pipeline_validator(),
        Commands::CacheResults => experimental_cmds::cmd_cache_results(),

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
