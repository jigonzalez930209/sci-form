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

#[cfg(feature = "alpha-cga")]
#[path = "experimental/test_cga.rs"]
mod test_cga;

#[cfg(feature = "beta-randnla")]
#[path = "experimental/test_randnla.rs"]
mod test_randnla;

#[cfg(feature = "beta-riemannian")]
#[path = "experimental/test_riemannian.rs"]
mod test_riemannian;

#[cfg(feature = "beta-kpm")]
#[path = "experimental/test_kpm.rs"]
mod test_kpm;

// EEQ, ALPB, D4 promoted to core — always compiled
#[path = "experimental/test_eeq.rs"]
mod test_eeq;

#[path = "experimental/test_alpb.rs"]
mod test_alpb;

#[path = "experimental/test_d4.rs"]
mod test_d4;

#[cfg(feature = "alpha-sdr")]
#[path = "experimental/test_sdr.rs"]
mod test_sdr;

#[cfg(feature = "beta-mbh")]
#[path = "experimental/test_mbh.rs"]
mod test_mbh;

#[cfg(feature = "beta-cpm")]
#[path = "experimental/test_cpm.rs"]
mod test_cpm;

#[cfg(feature = "alpha-gsm")]
#[path = "experimental/test_gsm.rs"]
mod test_gsm;

// Benchmarks — always compiled (EEQ/D4/ALPB are core, rest behind alpha/beta flags)
#[path = "experimental/test_benchmarks.rs"]
mod test_benchmarks;
#[path = "regression/test_experimental_comparison.rs"]
mod test_experimental_comparison;
#[path = "regression/test_extended_molecules.rs"]
mod test_extended_molecules;
#[path = "regression/test_gpu_candidates.rs"]
mod test_gpu_candidates;
