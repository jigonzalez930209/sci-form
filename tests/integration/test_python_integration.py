#!/usr/bin/env python3
"""Python binding integration tests for sci-form.

Run with pytest against an installed `sci_form` wheel:

    maturin build --release --out dist
    pip install dist/sciforma-*.whl
    pytest tests/integration/test_python_integration.py -v

Every public entry point of the `sci_form` surface is exercised at least
once, so a binding that silently changes shape or returns a fabricated value
is caught here rather than by a user.
"""

from __future__ import annotations

import math

import pytest

sci_form = pytest.importorskip("sci_form")

# Small water-like geometry: elements, flat coordinates, and a 3-vector form.
E = [6, 6, 8]
C = [0.0, 0.0, 0.0, 1.2, 0.6, 0.0, 2.4, 0.6, 0.0]
# Ethanol (9 atoms) with 27 coordinates.
CCO = [
    0.0, 0.0, 0.0,
    1.2, 0.6, 0.0,
    2.4, 0.6, 0.0,
    0.0, 0.1, 1.5,
    0.0, 0.2, 2.4,
    0.0, -0.1, 1.5,
    0.0, -0.2, 2.4,
    -0.1, 1.5, 0.6,
    -0.2, 2.4, 0.6,
]


def finite(*values) -> bool:
    return all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)


# ─── Core embedding and parsing ──────────────────────────────────────────────


def test_version():
    assert sci_form.version().startswith("sci-form")


def test_embed():
    result = sci_form.embed("CCO", seed=42)
    assert result.error is None
    assert result.num_atoms == 9
    assert len(result.coords) == 27
    # `Vec<u8>` is exposed as `bytes` by PyO3; the element list is
    # `list(result.elements)`.
    assert list(result.elements) == [6, 6, 8, 1, 1, 1, 1, 1, 1]
    positions = result.get_positions()
    assert len(positions) == 9
    assert all(len(p) == 3 for p in positions)


def test_embed_batch():
    results = sci_form.embed_batch(["CCO", "CCC"], seed=42)
    assert len(results) == 2
    assert all(r.num_atoms > 0 for r in results)


def test_embed_diverse():
    conformers = sci_form.embed_diverse("CCO", 3, 1.0, 42)
    assert 1 <= len(conformers) <= 3


def test_parse():
    # The parser adds implicit hydrogens, so benzene has 6 C + 6 H.
    parsed = sci_form.parse("c1ccccc1")
    assert parsed["num_atoms"] == 12
    assert parsed["num_bonds"] == 12
    assert len(parsed["atoms"]) == 12


def test_canonical_smiles_is_order_independent():
    assert sci_form.to_canonical_smiles("OCC") == sci_form.to_canonical_smiles("CCO")


# ─── Charges, surface, electrostatics ────────────────────────────────────────


def test_charges():
    result = sci_form.charges("CCO")
    assert len(result.charges) == 9
    assert finite(*result.charges)
    assert abs(result.total_charge) < 1e-6


def test_charges_configured():
    result = sci_form.compute_charges_configured("CCO", max_iter=12)
    assert len(result.charges) == 9
    assert finite(*result.charges)


def test_sasa():
    result = sci_form.sasa(E, C)
    assert result.total_sasa > 0.0


def test_population():
    result = sci_form.population(E, C)
    assert len(result.mulliken_charges) == 3
    assert finite(*result.mulliken_charges)


def test_dipole_returns_finite_3_vector():
    result = sci_form.dipole(E, C)
    assert len(result.vector) == 3
    assert finite(*result.vector)
    assert finite(result.magnitude)


def test_esp_grid_has_origin_dims_and_values():
    grid = sci_form.esp(E, C, spacing=0.5, padding=1.0)
    assert len(grid.origin) == 3
    assert finite(*grid.origin)
    assert len(grid.dims) == 3
    assert len(grid.values) == grid.dims[0] * grid.dims[1] * grid.dims[2]
    assert finite(*grid.values[:8])


def test_dos():
    result = sci_form.dos(E, C)
    assert len(result.energies) == len(result.total_dos)


def test_dos_multimethod():
    result = sci_form.compute_dos_multimethod(E, C)
    assert len(result.energies) == len(result.total_dos)
    assert result.gap > 0.0


def test_rmsd():
    result = sci_form.rmsd([0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0])
    assert result.rmsd == pytest.approx(0.0, abs=1e-12)


def test_rmsd_matrix():
    matrix = sci_form.rmsd_matrix([[0, 0, 0, 1, 0, 0], [0, 0, 0, 1.01, 0, 0]])
    assert len(matrix) == 2
    assert matrix[0][0] == pytest.approx(0.0)


# ─── Force fields and electronic structure ───────────────────────────────────


def test_uff_energy():
    assert finite(sci_form.uff_energy("CCO", CCO))


def test_mmff94_energy():
    assert finite(sci_form.mmff94_energy("CCO", CCO))


def test_pm3():
    assert finite(sci_form.pm3_calculate(E, C).total_energy)


def test_gfn1():
    assert finite(sci_form.gfn1_calculate(E, C).total_energy)


def test_gfn2():
    assert finite(sci_form.gfn2_calculate(E, C).total_energy)


def test_hf3c():
    assert finite(sci_form.hf3c_calculate(E, C).energy)


def test_cisd():
    result = sci_form.compute_cisd(E, C, n_states=3)
    assert len(result.excitations) <= 3
    assert finite(result.correlation_energy)


def test_eht():
    result = sci_form.eht_calculate(E, C)
    assert result.homo_energy < result.lumo_energy
    assert finite(result.homo_energy, result.lumo_energy, result.gap)


def test_eht_support():
    support = sci_form.eht_support(E)
    assert len(support.supported_elements) > 0
    assert not support.has_transition_metals


def test_eht_orbital_mesh():
    mesh = sci_form.eht_orbital_mesh(E, C, 0)
    assert set(mesh) >= {"vertices", "normals", "indices"}
    assert len(mesh["vertices"]) == len(mesh["normals"])


def test_eht_gradient():
    assert finite(sci_form.compute_eht_gradient(E, C).rms_gradient)


def test_optimize_geometry_eht():
    result = sci_form.optimize_geometry_eht(E, C)
    assert finite(result.energy)


def test_eht_or_uff_fallback():
    result = sci_form.eht_or_uff_fallback("CCO", E, C)
    assert result["mode"] in ("eht", "uff")
    assert finite(result["homo_energy"], result["lumo_energy"], result["gap"])


# ─── ML, descriptors and ANI ─────────────────────────────────────────────────


def test_ml_descriptors_and_prediction():
    assert finite(sci_form.ml_descriptors("CCO").molecular_weight)
    assert finite(sci_form.ml_predict("CCO").logp)


def test_whim():
    assert finite(sci_form.compute_whim(E, C).total_spread)


def test_rdf():
    assert len(sci_form.compute_rdf(E, C).radii) > 0


def test_getaway():
    getaway = sci_form.compute_getaway(E, C)
    assert finite(getaway.hat_total, getaway.r_total, getaway.h_mean)


def test_ani():
    result = sci_form.ani_calculate(E, C)
    assert len(result.forces) == 3


def test_aevs_tm():
    result = sci_form.compute_aevs_tm(E, C)
    assert result.n_atoms == 3
    assert result.aev_length > 0
    assert len(result.aevs) == 3


def test_train_random_forest():
    model = sci_form.train_random_forest([[0.0], [1.0], [2.0], [3.0]], [0.0, 1.0, 2.0, 3.0])
    assert finite(model.train_rmse)
    assert len(model.predictions) == 4


def test_train_gradient_boosting():
    model = sci_form.train_gradient_boosting([[0.0], [1.0], [2.0], [3.0]], [0.0, 1.0, 2.0, 3.0])
    assert finite(model.train_rmse)


# ─── Spectroscopy, IR and NMR ────────────────────────────────────────────────


def test_vibrational_analysis():
    result = sci_form.vibrational_analysis(E, C)
    assert result.n_atoms == 3
    assert len(result.modes) == 3 * 3


def test_ir_spectrum():
    result = sci_form.ir_spectrum(E, C)
    assert len(result.wavenumbers) == len(result.intensities)


def test_ir_spectrum_broadened():
    result = sci_form.compute_ir_spectrum_broadened(E, C, gamma=10.0, n_points=200)
    assert len(result.wavenumbers) == 200


def test_assign_peaks():
    result = sci_form.assign_peaks([3600.0], [80.0], [1, 1])
    assert len(result.assignments) == 1


def test_nmr_shifts():
    result = sci_form.nmr_shifts("CCO")
    assert len(result.h_shifts) > 0
    assert len(result.c_shifts) > 0


def test_nmr_couplings():
    assert len(sci_form.nmr_couplings("CCO", CCO)) > 0


def test_ensemble_j_couplings():
    assert len(sci_form.compute_ensemble_j_couplings("CCO", [C], [0.0])) > 0


def test_nmr_spectrum():
    assert len(sci_form.nmr_spectrum("CCO").ppm_axis) > 0


def test_hose_codes():
    codes = sci_form.hose_codes("CCO")
    assert len(codes) == 9
    assert all(len(c) == 3 for c in codes)


def test_stda_uvvis():
    assert len(sci_form.stda_uvvis(E, C).wavelengths_nm) > 0


def test_fluorescence():
    assert len(sci_form.fluorescence_spectrum(E, C).wavelengths_nm) > 0


# ─── Solvation ───────────────────────────────────────────────────────────────


def test_nonpolar_solvation():
    result = sci_form.nonpolar_solvation(E, C)
    assert finite(result.energy_kcal_mol)
    assert result.total_sasa > 0.0


def test_gb_solvation():
    result = sci_form.gb_solvation(E, C, [0.1, -0.1, 0.0])
    assert finite(result.total_energy_kcal_mol)
    assert finite(result.electrostatic_energy_kcal_mol)


def test_gb_solvation_rejects_a_zero_dielectric():
    """C-H4: `solvent_dielectric=0` used to divide by zero and return NaN."""
    with pytest.raises(ValueError):
        sci_form.compute_gb_solvation(E, C, [0.1, -0.1, 0.0], solvent_dielectric=0.0)


def test_gb_solvation_rejects_a_charge_vector_of_the_wrong_length():
    """C-M1: a short charge vector used to index past its end and abort."""
    with pytest.raises(ValueError):
        sci_form.compute_gb_solvation(E, C, [0.1])


def test_gb_solvation_rejects_a_coordinate_vector_of_the_wrong_length():
    """C-M1: a short coordinate array used to panic in the flat->triple split."""
    with pytest.raises(ValueError):
        sci_form.compute_gb_solvation([8, 1, 1], [0.0], [0.1, -0.1, 0.0])


def test_gb_solvation_rejects_a_non_finite_charge():
    with pytest.raises(ValueError):
        sci_form.compute_gb_solvation(E, C, [float("nan"), -0.1, 0.0])


def test_gb_solvation_rejects_a_negative_probe_radius():
    with pytest.raises(ValueError):
        sci_form.compute_gb_solvation(E, C, [0.1, -0.1, 0.0], probe_radius=-1.0)


def test_gb_solvation_keeps_the_vacuum_limit_working():
    result = sci_form.compute_gb_solvation(E, C, [0.1, -0.1, 0.0], solvent_dielectric=1.0)
    assert finite(result.total_energy_kcal_mol)


# ─── Open-shell SCF, gradients, ML and CIF (C-H2) ───────────────────────────


def test_uhf_open_shell_scf():
    """A methyl radical (doublet) must converge and report its orbitals."""
    CH3 = [6, 1, 1, 1]
    CH3_COORDS = [
        0.0, 0.0, 0.0,
        1.09, 0.0, 0.0,
        -0.36, 1.03, 0.0,
        -0.36, -0.52, 0.89,
    ]
    result = sci_form.compute_uhf(CH3, CH3_COORDS, 0, 2)
    assert result.converged
    assert result.total_energy < 0.0
    assert len(result.alpha_orbital_energies) == result.n_basis
    assert len(result.beta_orbital_energies) == result.n_basis
    assert len(result.mulliken_charges) == len(CH3)
    assert finite(result.total_energy, result.s2_expectation)
    # The dense coefficients are exposed row by row.
    assert len(result.alpha_coefficients) == result.n_basis
    assert len(result.alpha_coefficients[0]) == result.n_basis


def test_rohf_open_shell_scf():
    CH3 = [6, 1, 1, 1]
    CH3_COORDS = [
        0.0, 0.0, 0.0,
        1.09, 0.0, 0.0,
        -0.36, 1.03, 0.0,
        -0.36, -0.52, 0.89,
    ]
    result = sci_form.compute_rohf(CH3, CH3_COORDS, 0, 2)
    assert result.converged
    assert finite(result.total_energy)


def test_uhf_rejects_an_impossible_multiplicity():
    # A methyl radical has 9 electrons: a singlet (multiplicity 1) is
    # impossible, and a multiplicity of 0 does not exist.
    CH3 = [6, 1, 1, 1]
    CH3_COORDS = [
        0.0, 0.0, 0.0,
        1.09, 0.0, 0.0,
        -0.36, 1.03, 0.0,
        -0.36, -0.52, 0.89,
    ]
    with pytest.raises((ValueError, RuntimeError)):
        sci_form.compute_uhf(CH3, CH3_COORDS, 0, 1)
    with pytest.raises((ValueError, RuntimeError)):
        sci_form.compute_uhf(CH3, CH3_COORDS, 0, 0)


def test_pm3_gradient():
    result = sci_form.compute_pm3_gradient(E, C)
    assert len(result.gradients) == len(E)
    assert finite(result.energy, result.heat_of_formation)
    assert all(finite(*g) for g in result.gradients)


def test_xtb_gradient():
    result = sci_form.compute_xtb_gradient(E, C)
    assert len(result.gradients) == len(E)
    assert finite(result.energy)
    assert all(finite(*g) for g in result.gradients)


def test_uff_energy_breakdown():
    terms = sci_form.compute_uff_energy_breakdown("CCO", CCO)
    assert terms
    assert all(isinstance(kind, str) for kind, _ in terms)
    total = sum(value for _, value in terms)
    assert abs(total - sci_form.compute_uff_energy("CCO", CCO)) < 1e-6


def test_mmff94_energy_breakdown():
    terms = sci_form.compute_mmff94_energy_breakdown("CCO", CCO)
    assert terms
    total = sum(value for _, value in terms)
    assert abs(total - sci_form.compute_mmff94_energy("CCO", CCO)) < 1e-6


def test_tpsa():
    # Ethanol: one hydroxyl (20.23 Å² for the OH fragment set).
    assert sci_form.compute_tpsa("CCO") == pytest.approx(20.23, abs=0.01)
    assert sci_form.compute_tpsa("CCCC") == pytest.approx(0.0, abs=0.01)


def test_predict_ensemble():
    result = sci_form.predict_ensemble("CCO")
    assert finite(result.logp, result.logp_std, result.tpsa, result.confidence)
    assert 0.0 <= result.confidence <= 1.0
    assert result.veber.passes in (True, False)
    assert result.pka_acidic is None or finite(result.pka_acidic)


CIF_TEXT = """data_test
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Na1 Na 0.0 0.0 0.0 1.0
"""


def test_parse_cif():
    structure = sci_form.parse_cif(CIF_TEXT)
    assert structure.data_block == "test"
    assert structure.a == pytest.approx(5.0)
    assert len(structure.atom_sites) == 1
    assert structure.atom_sites[0].label == "Na1"
    assert structure.atom_sites[0].atomic_number == 11
    assert structure.atom_sites[0].occupancy == pytest.approx(1.0)


def test_write_cif_roundtrip():
    structure = sci_form.parse_cif(CIF_TEXT)
    text = sci_form.write_cif(structure)
    assert "data_test" in text
    again = sci_form.parse_cif(text)
    assert again.data_block == structure.data_block
    assert len(again.atom_sites) == len(structure.atom_sites)
    assert again.atom_sites[0].label == structure.atom_sites[0].label


def test_parse_cif_rejects_garbage():
    with pytest.raises(ValueError):
        sci_form.parse_cif("not a cif file")


# ─── Method names (C-H3) ────────────────────────────────────────────────────


def test_ir_spectrum_rejects_a_misspelled_line_shape():
    """The broadening used to fall back to Lorentzian for anything that was
    not "gaussian", so the typo silently changed the spectrum."""
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_ir_spectrum_broadened(
            E, C, method="pm3", gamma=10.0, wn_min=400.0, wn_max=4000.0,
            n_points=100, broadening="gauassian",
        )
    message = str(excinfo.value)
    assert "gaussian" in message
    assert "lorentzian" in message


def test_dos_multimethod_rejects_an_unknown_method():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_dos_multimethod(E, C, "ccsd")
    message = str(excinfo.value)
    for method in ("eht", "pm3", "xtb", "gfn1", "gfn2", "hf3c"):
        assert method in message


def test_orbital_mesh_rejects_an_unknown_method():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_orbital_mesh(E, C, "dft", 0, 0.4, 3.0, 0.02)
    message = str(excinfo.value)
    for method in ("eht", "pm3", "xtb", "hf3c"):
        assert method in message


def test_whim_rejects_an_unknown_weighting():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_whim(E, C, "gravitational")
    message = str(excinfo.value)
    assert "mass" in message
    assert "volume" in message


# ─── Rings, fingerprints, clustering ─────────────────────────────────────────


def test_sssr():
    assert len(sci_form.sssr("c1ccccc1").rings) == 1


def test_ecfp_and_tanimoto():
    a = sci_form.ecfp("CCO")
    assert a.n_bits > 0
    assert sci_form.tanimoto(a, a) == pytest.approx(1.0)


def test_ecfp_batch():
    assert len(sci_form.compute_ecfp_batch(["CCO", "CCC"])) == 2


def test_butina_cluster_and_filter_diverse():
    conformers = [[0, 0, 0, 1, 0, 0], [0, 0, 0, 1.01, 0, 0]]
    assert len(sci_form.butina_cluster(conformers, 1.0).centroid_indices) >= 1
    assert len(sci_form.filter_diverse(conformers, 0.5)) >= 1


# ─── Stereochemistry, graph features, reactivity ─────────────────────────────


def test_stereo_analysis():
    assert sci_form.stereo_analysis("C[C@H](F)Cl").n_stereocenters == 1


def test_graph_features():
    features = sci_form.graph_features("c1ccccc1")
    # 6 aromatic carbons + 6 implicit hydrogens.
    assert len(features.aromatic_atoms) == 12
    assert all(features.aromatic_atoms[:6])
    assert len(features.aromatic_bonds) == 6


def test_frontier_and_fukui():
    frontier = sci_form.frontier_descriptors(E, C)
    assert finite(frontier.homo_energy, frontier.lumo_energy)
    assert len(frontier.dual_descriptor) == 3
    assert len(frontier.homo_atom_contributions) == 3
    assert len(sci_form.fukui_descriptors(E, C).f_plus) == 3


def test_reactivity_ranking():
    ranking = sci_form.reactivity_ranking(E, C)
    assert len(ranking.nucleophilic_attack_sites) == 3
    assert len(ranking.electrophilic_attack_sites) == 3


def test_uvvis_and_empirical_pka():
    assert len(sci_form.uvvis_spectrum(E, C).peaks) > 0
    pka = sci_form.empirical_pka("CC(=O)O")
    assert len(pka.acidic_sites) + len(pka.basic_sites) > 0


def test_bond_orders_wiberg():
    assert len(sci_form.bond_orders(E, C).wiberg) == 3


def test_bond_orders():
    orders = sci_form.bond_orders(E, C)
    assert len(orders.wiberg) == 3
    assert len(orders.mayer) == 3


# ─── Materials and periodic systems ──────────────────────────────────────────


def test_unit_cell():
    cell = sci_form.unit_cell(10.0, 10.0, 10.0)
    assert len(cell.lattice) == 3


def test_assemble_framework():
    structure = sci_form.assemble_framework("pcu", 30, "octahedral", 10.0, 1)
    assert structure.num_atoms > 0


def test_space_groups():
    assert sci_form.space_group_by_number(225).symbol
    assert sci_form.space_group_by_symbol("P 21/c").number == 14
    assert len(sci_form.all_space_groups()) == 230


def test_periodic_molecule_and_hapticity():
    assert sci_form.build_periodic_molecule([6], [[0.5, 0.5, 0.5]], 10, 10, 10).n_atoms == 1
    analysis = sci_form.detect_hapticity(
        [26] + [6] * 10, [0.0, 0.0, 0.0] + [i * 1.2 for i in range(1, 11)]
    )
    assert analysis.n_interactions >= 0


def test_optimize_framework():
    assert finite(sci_form.optimize_framework("CCO", CCO).energy)


# ─── Reactions and transport ─────────────────────────────────────────────────


def test_smirks_roundtrip():
    transform = sci_form.parse_smirks("[C:1][OH:2]>>[C:1][Cl:2]")
    assert transform is not None
    result = sci_form.apply_smirks("[C:1][OH:2]>>[C:1][Cl:2]", "CCO")
    assert result.success
    assert len(result.products) == 1


def test_transport():
    batch = sci_form.pack_conformers([sci_form.embed("CCO", seed=42)])
    assert "coords" in batch.column_names
    tasks = sci_form.split_worker_tasks(["CCO", "CCC"], 4, 42)
    assert len(tasks) >= 1
    assert sci_form.estimate_workers(100) >= 1


# ─── Dynamics and planning ───────────────────────────────────────────────────


def test_md_trajectory():
    trajectory = sci_form.md_trajectory("CCO", CCO, 10, 1.0, 42)
    # The trajectory stores the initial frame plus one frame per step.
    assert len(trajectory.frames) == 11
    # Flat coordinate buffer: 3 floats per atom (9 atoms).
    assert len(trajectory.frames[0].coords) == 27


def test_search_conformers():
    result = sci_form.search_conformers("CCO", 3, 42, 0.5)
    assert result.generated >= 1


def test_simplified_neb_path():
    path = sci_form.simplified_neb_path("CCO", CCO, CCO, 5)
    assert len(path.images) == 5


def test_system_capabilities_and_plan():
    capabilities = sci_form.system_capabilities(E)
    assert capabilities.embed.available
    plan = sci_form.system_method_plan(E)
    assert plan.geometry.recommended is not None


def test_compare_methods():
    result = sci_form.compare_methods("CCO", E, C)
    assert len(result.comparisons) >= 2


# ─── Submodules ──────────────────────────────────────────────────────────────


def test_potentials_submodule():
    assert finite(sci_form.potentials.muller_brown_2d(0.0, 0.0))
    assert finite(sci_form.potentials.double_well_1d(0.0))


def test_alpha_and_beta_submodules_are_attributes():
    # Regression: the submodules used to be registered under the invalid
    # identifier "sci_form.alpha" / "sci_form.beta", which made
    # `sci_form.alpha` unreachable as an attribute.
    assert hasattr(sci_form, "alpha")
    assert hasattr(sci_form, "beta")
    assert hasattr(sci_form, "potentials")


def test_alpha_submodule_import():
    # The `alpha` tier is feature-gated: a wheel built without the
    # `alpha-*` features ships the (empty) submodule but not its functions.
    alpha = pytest.importorskip("sci_form.alpha")
    if not hasattr(alpha, "alpha_boys_function"):
        pytest.skip("the alpha tier is not built into this wheel")

    assert finite(alpha.alpha_boys_function(0, 2.0))


def test_beta_submodule_import():
    # The `beta` tier is feature-gated in the same way.
    beta = pytest.importorskip("sci_form.beta")
    if not hasattr(beta, "kpm_dos"):
        pytest.skip("the beta tier is not built into this wheel")

    assert len(beta.kpm_dos(E, C).energies) > 0
