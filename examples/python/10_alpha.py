"""Alpha-tier experimental modules. """
from sci_form import embed
from sci_form.alpha import (
    dft_calculate,
    reaxff_energy,
    reaxff_gradient,
    eem_charges,
    alpha_compute_mlff,
    alpha_compute_aevs,
    alpha_boys_function,
    alpha_eri_ssss,
    alpha_schwarz_bound,
    rotate_dihedral_cga,
    refine_torsion_cga,
    gsm_interpolate,
    gsm_backend_plan,
    gsm_compare_backends,
    gsm_find_ts,
    sdr_embed_py,
    reaction_dynamics_3d,
)
from common import require_ok


def run():
    print("=== 10 alpha ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    dft = dft_calculate(conf.elements, conf.coords, method="pbe")
    print("DFT", dft.energy, dft.gap, dft.converged, dft.scf_iterations, dft.n_basis, dft.xc_energy, dft.nuclear_repulsion)

    rgrad = reaxff_gradient(conf.elements, conf.coords)
    print("ReaxFF energy", rgrad.energy_kcal_mol, "grad len", len(rgrad.gradient))

    re = reaxff_energy(conf.elements, conf.coords)
    print("ReaxFF energy breakdown", re.bonded, re.coulomb, re.van_der_waals, re.total)

    eem = eem_charges(conf.elements, conf.coords)
    print("EEM charges", len(eem.charges), "total", eem.total_charge)

    aevs = alpha_compute_aevs(conf.elements, conf.coords, radial_cutoff=6.0, angular_cutoff=3.5)
    print("AEV", aevs.n_atoms, aevs.aev_length, len(aevs.aevs))

    # MLFF requires trained element nets; AEVs computed separately above

    print("Boys:", alpha_boys_function(2, 1.5))
    print("ERI:", alpha_eri_ssss(1.0, [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0], 1.0, [0.0, 1.0, 0.0], 1.0, [0.0, 0.0, 1.0]))
    print("Schwarz:", alpha_schwarz_bound(1.0, [0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]))

    # CGA: rotate a dihedral
    coords = conf.coords
    axis_a = [coords[0], coords[1], coords[2]]
    axis_b = [coords[3], coords[4], coords[5]]
    newcoords = rotate_dihedral_cga(coords, [0, 3, 4, 5], axis_a, axis_b, 0.5)
    print("rotated coords len", len(newcoords))
    # refine not run (needs bonds); call via None?


    # GSM interpolation simple
    interp = gsm_interpolate(conf.coords, conf.coords, 0.5)
    print("gsm interpolate len", len(interp))

    bp = gsm_backend_plan("CCO")
    print("gsm backend plan", bp)
    cmpb = gsm_compare_backends("CCO", conf.coords, None)
    print("gsm compare", len(cmpb))
    print("gsm compare", cmpb)
    ts = gsm_find_ts("CCO", conf.coords, conf.coords, n_nodes=5, method="uff")
    print("gsm find ts", ts)
    # gsm_find_ts_with_method not exposed; skip

    sdr = sdr_embed_py([(0,1,1.5)], 9)
    print("SDR", sdr.coords[:3], sdr.max_distance_error)
    print("SDR", sdr)

    rd = reaction_dynamics_3d(["CC"], ["CC"], "uff", 5, 2, 2, None, 42)
    print("reaction dynamics", rd.ts_frame_index, rd.activation_energy_kcal_mol, len(rd.frames))
    print("reaction dynamics", rd)


if __name__ == "__main__":
    run()
