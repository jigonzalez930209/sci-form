"""Charges, SASA, dipole, population, DOS, RMSD."""
from sci_form import embed, charges, sasa, dipole, dos, rmsd, population
from common import require_ok


def run():
    print("=== 02 charges/sasa/dipole ===")
    cr = charges("CCO")
    print(cr)
    print("charges:", cr.charges)
    print("total_charge:", cr.total_charge, "iterations:", cr.iterations)

    cr2 = charges("[NH4+]")
    print("charges [NH4+] total:", cr2.total_charge, cr2.iterations, cr2.charges)
    conf = require_ok(embed("CCO", seed=42), "embed CCO")

    s = sasa(conf.elements, conf.coords, probe_radius=1.4)
    print(s, "atom_sasa len", len(s.atom_sasa))

    d = dipole(conf.elements, conf.coords)
    print(d, "vector", d.vector, "magnitude", d.magnitude, d.unit)

    pop = population(conf.elements, conf.coords)
    print(pop)
    print("mulliken:", pop.mulliken_charges)
    print("lowdin:", pop.lowdin_charges)
    print("total_mulliken", pop.total_charge_mulliken)

    dosr = dos(conf.elements, conf.coords, sigma=0.2, e_min=-25.0, e_max=5.0, n_points=300)
    print("dos points", len(dosr.energies), "sigma", dosr.sigma)
    print("max dos", max(dosr.total_dos))

    c2 = require_ok(embed("CCO", seed=123), "embed CCO 2")
    al = rmsd(conf.coords, c2.coords)
    print(al, "rmsd", al.rmsd, "aligned len", len(al.aligned_coords))


if __name__ == "__main__":
    run()
