"""EHT, bond orders, frontier, UFF/MMFF, ANI, ESP, orbitals."""
from sci_form import (
    embed,
    eht_calculate,
    eht_orbital_mesh,
    bond_orders,
    frontier_descriptors,
    uff_energy,
    mmff94_energy,
    ani_calculate,
    esp,
    orbital_mesh,
)
from common import require_ok


def run():
    print("=== 03 eht/forcefields ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    eht = eht_calculate(conf.elements, conf.coords, k=1.75)
    print(eht)
    print("homo", eht.homo_energy, "lumo", eht.lumo_energy, "gap", eht.gap)
    print("n_electrons", eht.n_electrons, "homo_index", eht.homo_index, "lumo_index", eht.lumo_index)

    mesh = eht_orbital_mesh(conf.elements, conf.coords, eht.homo_index, isovalue=0.02, spacing=0.3)
    print("mesh keys", list(mesh.keys()))
    print("triangles", mesh["num_triangles"], "vertices len", len(mesh["vertices"]))

    bo = bond_orders(conf.elements, conf.coords)
    print(bo)
    print("pairs", len(bo.atom_pairs), "wiberg", bo.wiberg[:3])

    fd = frontier_descriptors(conf.elements, conf.coords)
    print(fd)
    print("gap", fd.gap, "homo_atom_contrib", len(fd.homo_atom_contributions))

    print("UFF:", uff_energy("CCO", conf.coords))
    print("MMFF94:", mmff94_energy("CCO", conf.coords))

    ani = ani_calculate(conf.elements, conf.coords)
    print(ani)

    esp_res = esp(conf.elements, conf.coords)
    print(esp_res)

    om = orbital_mesh(conf.elements, conf.coords, method="eht", mo_index=eht.homo_index, isovalue=0.02)
    print("orbital mesh", om.method, om.mo_index, om.num_triangles, "vertices", len(om.vertices), "gap", om.gap)


if __name__ == "__main__":
    run()
