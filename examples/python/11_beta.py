"""Beta-tier modules. """
from sci_form import embed
from sci_form.beta import (
    kpm_dos,
    kpm_mulliken,
    mbh_frequencies,
    eht_randnla,
    cpm_charges,
    beta_psd_distance,
    beta_psd_projection,
)
from common import require_ok


def run():
    print("=== 11 beta ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    kpm = kpm_dos(conf.elements, conf.coords, order=128)
    print("KPM points", len(kpm.energies), "order", kpm.order, "energies[0]", kpm.energies[0])

    km = kpm_mulliken(conf.elements, conf.coords)
    print("KPM mulliken keys", [a for a in dir(km) if not a.startswith("_")])


    mbh = mbh_frequencies(conf.elements, conf.coords, "CCO", fd_step=0.005)
    print("MBH", mbh)

    rn = eht_randnla(conf.elements, conf.coords, sketch_size=16)
    print("RandNLA orbitals", len(rn.orbital_energies), "residual", rn.residual_error)

    cpm = cpm_charges(conf.elements, conf.coords, potential=0.0)
    print("CPM", cpm.charges[:3], cpm.total_charge, cpm.energy)

    dist = beta_psd_distance([1.0, 0.0, 0.0, 1.0], [2.0, 0.0, 0.0, 2.0], 2)
    print("PSD dist", dist)
    proj = beta_psd_projection([-1.0, 0.0, 0.0, 1.0], 2)
    print("PSD proj", proj)


if __name__ == "__main__":
    run()
