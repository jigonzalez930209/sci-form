#[path = "regression/test_algorithm_validation.rs"]
mod test_algorithm_validation;
#[path = "regression/test_benchmark_100.rs"]
mod test_benchmark_100;
#[path = "regression/test_bounds_compare.rs"]
mod test_bounds_compare;
#[path = "regression/test_charges.rs"]
mod test_charges;
#[path = "regression/test_diverse_molecules.rs"]
mod test_diverse_molecules;
#[path = "regression/test_eht.rs"]
mod test_eht;
#[path = "regression/test_eht_metal_references.rs"]
mod test_eht_metal_references;
#[path = "regression/test_embedding_trace.rs"]
mod test_embedding_trace;
#[path = "regression/test_energy_comparison.rs"]
mod test_energy_comparison;
#[path = "regression/test_ensemble_rmsd.rs"]
mod test_ensemble_rmsd;
#[path = "regression/test_gdb20.rs"]
mod test_gdb20;
#[path = "regression/test_gdb20_rmsd.rs"]
mod test_gdb20_rmsd;
#[path = "regression/test_geometry_quality.rs"]
mod test_geometry_quality;
#[path = "regression/test_grad_isolate.rs"]
mod test_grad_isolate;
#[path = "regression/test_grad_torsion_detail.rs"]
mod test_grad_torsion_detail;
#[path = "regression/test_gradient_check.rs"]
mod test_gradient_check;
#[path = "regression/test_gradient_perterm.rs"]
mod test_gradient_perterm;
#[path = "regression/test_metal_experimental_geometry.rs"]
mod test_metal_experimental_geometry;
#[path = "regression/test_new_features.rs"]
mod test_new_features;
#[path = "regression/test_new_pipeline.rs"]
mod test_new_pipeline;
#[path = "regression/test_parallel_public_api.rs"]
mod test_parallel_public_api;
#[path = "regression/test_phase_c.rs"]
mod test_phase_c;
#[path = "regression/test_phase_c_validation.rs"]
mod test_phase_c_validation;
#[path = "regression/test_pipeline_trace.rs"]
mod test_pipeline_trace;
#[path = "regression/test_reaction_dft_mesh.rs"]
mod test_reaction_dft_mesh;
#[path = "regression/test_roadmap_validation.rs"]
mod test_roadmap_validation;
#[path = "regression/test_sasa.rs"]
mod test_sasa;
#[path = "regression/test_session_algorithms.rs"]
mod test_session_algorithms;
#[path = "regression/test_spectroscopy.rs"]
mod test_spectroscopy;
#[path = "regression/test_step_by_step.rs"]
mod test_step_by_step;
#[path = "regression/test_tet_centers.rs"]
mod test_tet_centers;

#[path = "experimental/test_d4.rs"]
mod test_d4;

// Benchmarks — always compiled (EEQ/D4/ALPB are core, rest behind alpha/beta flags)
#[path = "regression/test_experimental_comparison.rs"]
mod test_experimental_comparison;
#[path = "regression/test_extended_molecules.rs"]
mod test_extended_molecules;
#[path = "regression/test_gpu_candidates.rs"]
mod test_gpu_candidates;

// The `tests/experimental/` suite (EEQ, ALPB, D4, CPM and the alpha/beta
// gated modules) is aggregated by `tests/experimental_gates.rs` so it can run
// as its own CI target.
