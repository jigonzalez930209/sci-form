"""Solvation and materials: nonpolar, GB, unit cell, frameworks. """
from sci_form import embed, nonpolar_solvation, gb_solvation, charges as gasteiger_charges, unit_cell, assemble_framework

from common import require_ok


def run():
    print("=== 08 solvation/materials ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    q = gasteiger_charges("CCO").charges
    np_solv = nonpolar_solvation(conf.elements, conf.coords)
    print("nonpolar", np_solv.energy_kcal_mol, "total_sasa", np_solv.total_sasa, "atom contrib", len(np_solv.atom_contributions))

    gb = gb_solvation(conf.elements, conf.coords, q)
    print("gb", gb.total_energy_kcal_mol, gb.electrostatic_energy_kcal_mol, gb.nonpolar_energy_kcal_mol, gb.born_radii[:3])

    cell = unit_cell(10.0, 10.0, 10.0)
    print("cubic cell volume", cell.volume, "lattice", cell.lattice)
    cell2 = unit_cell(5.0, 5.0, 12.0, 90, 90, 120)
    print("hex cell volume", cell2.volume)

    mof = assemble_framework(topology="pcu", metal=30, geometry="octahedral", lattice_a=26.3)
    print("MOF atoms", mof.num_atoms, "elements len", len(mof.elements), "cart coords len", len(mof.cart_coords), "frac", len(mof.frac_coords))


if __name__ == "__main__":
    run()
