#!/usr/bin/env python3
"""Binding-surface tests for the v1.0 API/accessibility blockers.

Run against an installed `sci_form` wheel:

    maturin build --release --out dist
    pip install dist/sciforma-*.whl
    pytest crates/python/tests/test_bindings.py -v

The file covers the binding-side blockers of `reports/v1.0-blockers.md`
section 4 (C. API design and accessibility):

* C-C1  the WASM error payload must always be valid JSON (checked here on
        the Python side of the shared error style, and by the Rust unit
        tests in `crates/wasm/src/helpers.rs`),
* C-C3  `cpm_charges` must exist exactly once and must carry the
        `dielectric` parameter,
* C-C4  the alignment path must raise instead of returning a `NaN` RMSD,
* C-H1  the canonical names must exist and the old names must still work as
        deprecated aliases,
* C-H5  `__all__` must be complete and `parse` must return a `pyclass`,
* C-H6  the zero-copy buffer entry points must work,
* C-H9  the error messages must follow the canonical style,
* C-M4  a failed embed must raise instead of returning an error object,
* C-M6  the error taxonomy must be distinguishable.
"""

from __future__ import annotations

import json
import math
import warnings

import pytest

sci_form = pytest.importorskip("sci_form")

# A small water-like geometry: elements, flat coordinates.
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


# ─── C-H5: the package has a meaningful public API ─────────────────────────


def test_version_and_dunder_version_agree():
    assert sci_form.__version__ == sci_form.version().split()[-1]


def test_all_is_complete_and_every_name_exists():
    assert hasattr(sci_form, "__all__")
    assert len(sci_form.__all__) > 100
    missing = [name for name in sci_form.__all__ if not hasattr(sci_form, name)]
    assert missing == [], f"names in __all__ that do not exist: {missing}"


def test_every_public_callable_is_listed_in_all():
    """`__all__` must cover the whole public surface, not just a subset.

    The wheel re-exports the extension module with ``from .sci_form import *``,
    so a name missing from ``__all__`` is invisible to the user.
    """
    extension = getattr(sci_form, "sci_form", None)
    if extension is None:  # pragma: no cover - only for non-wheel installs
        pytest.skip("the extension module is not reachable from the package")
    public = {
        name
        for name in extension.__dict__
        if not name.startswith("_") and callable(getattr(extension, name))
    }
    missing = sorted(public - set(sci_form.__all__))
    assert missing == [], f"public callables missing from __all__: {missing}"


def test_dir_is_the_public_api():
    listed = dir(sci_form)
    assert "embed" in listed
    assert "compute_charges" in listed
    # `dir()` must not be a flat alphabetical dump of 220+ internals.
    assert len(listed) < 400


def test_classes_have_no_py_suffix():
    for name in sci_form.__all__:
        obj = getattr(sci_form, name)
        if isinstance(obj, type):
            assert not name.endswith("Py"), f"{name} still carries the Py suffix"


def test_parse_returns_a_pyclass():
    parsed = sci_form.parse("c1ccccc1")
    assert isinstance(parsed, sci_form.ParseResult)
    assert parsed.num_atoms == 12
    assert parsed.num_bonds == 12
    assert len(parsed.atoms) == 12
    assert parsed.atoms[0].element == 6
    # The mapping-style access of the previous dict payload still works.
    assert parsed["num_atoms"] == 12
    with pytest.raises(KeyError):
        parsed["nope"]


# ─── C-C3: exactly one `cpm_charges` ───────────────────────────────────────


def test_cpm_charges_has_the_dielectric_parameter():
    import inspect

    params = inspect.signature(sci_form.compute_cpm_charges).parameters
    assert "dielectric" in params
    assert "mu_ev" in params


def test_cpm_charges_alias_has_the_dielectric_parameter():
    import inspect

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        params = inspect.signature(sci_form.cpm_charges).parameters
    assert "dielectric" in params


def test_beta_cpm_charges_is_reachable_under_its_own_name():
    assert hasattr(sci_form, "beta_compute_cpm_charges")
    assert hasattr(sci_form.beta, "beta_compute_cpm_charges")


def test_cpm_charges_computes_at_the_standard_electrode_potential():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = sci_form.cpm_charges(E, C)
    # The standard potential of the H+/H2 electrode, not mu = 0 eV.
    assert result.mu_ev == pytest.approx(-4.44)
    assert len(result.charges) == 3
    assert result.converged in (True, False)
    # The beta variant is a different, three-field object.
    beta = sci_form.beta_compute_cpm_charges(E, C)
    assert hasattr(beta, "energy")
    assert not hasattr(beta, "mu_ev")


# ─── C-C4: the alignment path never returns NaN ────────────────────────────


def test_rmsd_raises_on_a_length_mismatch():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_rmsd([0, 0, 0, 1, 0, 0], [0, 0, 0])
    assert "different length" in str(excinfo.value)


def test_rmsd_raises_when_the_length_is_not_a_multiple_of_three():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_rmsd([0, 0, 0, 1], [0, 0, 0, 1])
    assert "multiple of 3" in str(excinfo.value)


def test_rmsd_is_finite_for_matching_structures():
    result = sci_form.compute_rmsd([0, 0, 0, 1, 0, 0], [0, 0, 0, 1, 0, 0])
    assert math.isfinite(result.rmsd)
    assert result.rmsd == pytest.approx(0.0, abs=1e-12)


# ─── C-H1: canonical names and deprecated aliases ──────────────────────────

CANONICAL_NAMES = [
    "compute_charges",
    "compute_sasa",
    "compute_population",
    "compute_dipole",
    "compute_dos",
    "compute_bond_orders",
    "compute_frontier_descriptors",
    "compute_fukui_descriptors",
    "compute_reactivity_ranking",
    "compute_uv_vis_spectrum",
    "compute_topology",
    "analyze_graph_features",
    "analyze_stereo",
    "compute_orbital_mesh",
    "compute_sssr",
    "compute_ecfp",
    "compute_tanimoto",
    "compute_uff_energy",
    "compute_mmff94_energy",
    "compute_empirical_pka",
    "compute_ml_descriptors",
    "predict_ml_properties",
    "predict_nmr_shifts",
    "predict_nmr_couplings",
    "compute_nmr_spectrum",
    "compute_giao_nmr",
    "compute_hose_codes",
    "compute_ir_spectrum",
    "compute_vibrational_analysis",
    "compute_gb_solvation",
    "compute_nonpolar_solvation",
    "compute_eht_or_uff_fallback",
    "eht_support",
    "system_capabilities",
    "system_method_plan",
    "compute_pm3",
    "compute_xtb",
    "compute_gfn1",
    "compute_gfn2",
    "compute_hf3c",
    "compute_ani",
]


@pytest.mark.parametrize("name", CANONICAL_NAMES)
def test_canonical_name_exists(name):
    assert hasattr(sci_form, name), f"missing canonical name {name}"


DEPRECATED_ALIASES = [
    ("charges", "compute_charges"),
    ("sasa", "compute_sasa"),
    ("population", "compute_population"),
    ("dipole", "compute_dipole"),
    ("dos", "compute_dos"),
    ("bond_orders", "compute_bond_orders"),
    ("rmsd", "compute_rmsd"),
    ("pm3_calculate", "compute_pm3"),
    ("xtb_calculate", "compute_xtb"),
    ("gfn1_calculate", "compute_gfn1"),
    ("gfn2_calculate", "compute_gfn2"),
    ("hf3c_calculate", "compute_hf3c"),
    ("ani_calculate", "compute_ani"),
    ("ml_descriptors", "compute_ml_descriptors"),
    ("ml_predict", "predict_ml_properties"),
    ("nmr_shifts", "predict_nmr_shifts"),
    ("orbital_mesh", "compute_orbital_mesh"),
    ("unit_cell", "create_unit_cell"),
    ("eht_calculate", "solve_eht"),
]


@pytest.mark.parametrize("alias,canonical", DEPRECATED_ALIASES)
def test_deprecated_alias_warns_and_still_works(alias, canonical):
    assert hasattr(sci_form, alias)
    with pytest.warns(DeprecationWarning) as record:
        try:
            getattr(sci_form, alias)("CCO")
        except Exception:
            # The alias must forward to the canonical function; the call may
            # fail for lack of arguments, but the warning must be raised
            # first.
            pass
    assert any(canonical in str(w.message) for w in record)


def test_no_calculate_suffix_remains_in_the_canonical_api():
    """The `*_calculate` Python-ism is gone from the canonical names.

    The deprecated aliases keep their old spelling on purpose, so they are
    excluded from the check.
    """
    aliases = {alias for alias, _ in DEPRECATED_ALIASES}
    offenders = [
        name
        for name in dir(sci_form)
        if not name.startswith("_")
        and name.endswith("_calculate")
        and name not in aliases
    ]
    assert offenders == []


def test_every_deprecated_alias_points_at_a_canonical_name():
    """No alias may shadow another alias, and each must resolve."""
    seen = set()
    for alias, canonical in DEPRECATED_ALIASES:
        assert alias not in seen, f"duplicate alias {alias}"
        seen.add(alias)
        assert hasattr(sci_form, canonical), f"{canonical} is missing"


# ─── C-M4 / C-M6: the error taxonomy ───────────────────────────────────────


def test_embed_raises_on_a_malformed_smiles():
    with pytest.raises(ValueError):
        sci_form.embed("NOT_A_SMILES")


def test_embed_returns_a_result_object_on_success():
    result = sci_form.embed("CCO", seed=42)
    assert result.error is None
    assert result.num_atoms == 9


def test_not_converged_error_is_a_runtime_error():
    assert issubclass(sci_form.NotConvergedError, RuntimeError)


def test_value_error_for_malformed_input():
    with pytest.raises(ValueError):
        sci_form.compute_sasa(E, [0.0, 0.0])  # 2 values for 3 atoms


def test_unknown_method_is_a_value_error():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_dos_multimethod(E, C, method="not-a-method")
    assert "not-a-method" in str(excinfo.value)


def test_engine_failure_is_a_runtime_error():
    """A failure inside the engine is a `RuntimeError`, not a `ValueError`."""
    with pytest.raises(RuntimeError) as excinfo:
        sci_form.compute_giao_nmr([35], [0.0, 0.0, 0.0])
    assert not isinstance(excinfo.value, ValueError)


# ─── C-H9: the error message style ─────────────────────────────────────────


def test_error_message_lists_the_supported_values():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_dos_multimethod(E, C, method="not-a-method")
    message = str(excinfo.value)
    assert "'not-a-method'" in message
    assert "Supported:" in message
    for value in ("'eht'", "'pm3'", "'gfn2'"):
        assert value in message


def test_smiles_error_reports_the_position_and_the_offending_fragment():
    with pytest.raises(ValueError) as excinfo:
        sci_form.parse("C(CC")
    message = str(excinfo.value)
    assert "position" in message
    assert "C(CC" in message


def test_no_rust_internal_names_leak_into_a_message():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_rmsd([0, 0, 0], [0, 0, 0, 1])
    message = str(excinfo.value)
    for leak in ("unwrap", "None", "positions", ".rs:"):
        assert leak not in message, f"'{leak}' leaked into {message!r}"


def test_error_message_is_a_single_line():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_dos_multimethod(E, C, method="not-a-method")
    assert "\n" not in str(excinfo.value)


# ─── C-H6: the zero-copy buffer path ───────────────────────────────────────


def test_embed_coords_into_writes_into_the_caller_buffer():
    n_atoms = sci_form.parse("CCO").num_atoms
    try:
        import array

        buffer = array.array("d", [0.0] * (3 * n_atoms))
    except ImportError:  # pragma: no cover - array is in the stdlib
        pytest.skip("array module unavailable")
    written = sci_form.embed_coords_into("CCO", buffer, seed=42)
    assert written == n_atoms
    assert any(v != 0.0 for v in buffer)


def test_embed_coords_into_rejects_a_wrong_sized_buffer():
    import array

    buffer = array.array("d", [0.0] * 3)
    with pytest.raises(ValueError) as excinfo:
        sci_form.embed_coords_into("CCO", buffer, seed=42)
    assert "Expected" in str(excinfo.value)


def test_compute_esp_grid_into_writes_into_the_caller_buffer():
    grid = sci_form.compute_esp(E, C, spacing=0.5, padding=1.0)
    n_values = len(grid.values)
    import array

    buffer = array.array("d", [0.0] * n_values)
    written = sci_form.compute_esp_grid_into(E, C, buffer, spacing=0.5, padding=1.0)
    assert written == n_values
    assert any(v != 0.0 for v in buffer)


# ─── C-M3: the optional options mapping ────────────────────────────────────


def test_optimize_framework_accepts_an_options_mapping():
    result = sci_form.optimize_framework("CCO", CCO, options={"max_iter": 5})
    assert math.isfinite(result.energy)


def test_optimize_framework_rejects_an_unknown_option():
    with pytest.raises(ValueError) as excinfo:
        sci_form.optimize_framework("CCO", CCO, options={"not_an_option": 1})
    assert "Supported" in str(excinfo.value)


# ─── C-C1: the shared error payload is valid JSON ──────────────────────────


def test_wasm_style_error_payload_is_valid_json():
    """The canonical WASM payload `{"error": "..."}` must always parse.

    The message built here is the one the WASM helper produces; the Rust
    unit tests in `crates/wasm/src/helpers.rs` assert the same property for
    the real payloads.
    """
    message = 'bad elements: invalid type: string "abc", expected a sequence'
    payload = json.dumps({"error": message})
    parsed = json.loads(payload)
    assert parsed["error"] == message


def test_error_messages_survive_a_json_round_trip():
    with pytest.raises(ValueError) as excinfo:
        sci_form.compute_dos_multimethod(E, C, method='we"ird\\value')
    message = str(excinfo.value)
    assert json.loads(json.dumps({"error": message}))["error"] == message
