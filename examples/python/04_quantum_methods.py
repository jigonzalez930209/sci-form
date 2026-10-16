"""Semi-empirical quantum methods."""
from sci_form import (
    embed,
    pm3_calculate,
    xtb_calculate,
    gfn1_calculate,
    gfn2_calculate,
    hf3c_calculate,
)
from common import require_ok


def run():
    print("=== 04 quantum ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    pm3 = pm3_calculate(conf.elements, conf.coords)
    print("PM3 energy", pm3.total_energy, "gap", pm3.gap, "converged", pm3.converged, "iterations", pm3.scf_iterations)
    print("HOF", pm3.heat_of_formation, "n_basis", pm3.n_basis, "n_electrons", pm3.n_electrons)
    print("mulliken", pm3.mulliken_charges[:3])

    xtb = xtb_calculate(conf.elements, conf.coords)
    print("xTB energy", xtb.total_energy, "gap", xtb.gap, "converged", xtb.converged)

    g1 = gfn1_calculate(conf.elements, conf.coords)
    print("GFN1 energy", g1.total_energy, "dispersion", g1.dispersion_energy, "converged", g1.converged)

    g2 = gfn2_calculate(conf.elements, conf.coords)
    print("GFN2 energy", g2.total_energy, "dispersion", g2.dispersion_energy, "halogen", g2.halogen_bond_energy, "converged", g2.converged)

    hf = hf3c_calculate(conf.elements, conf.coords, max_scf_iter=50, n_cis_states=3, corrections=True)
    print("HF-3c", hf)


if __name__ == "__main__":
    run()
