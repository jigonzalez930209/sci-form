"""Reactivity, ML descriptors, stereochemistry, rings, fingerprints."""
from sci_form import (
    embed,
    fukui_descriptors,
    reactivity_ranking,
    empirical_pka,
    ml_descriptors,
    ml_predict,
    stereo_analysis,
    sssr,
    ecfp,
    tanimoto,
    bond_orders,
    frontier_descriptors,
)
from common import require_ok


def run():
    print("=== 05 reactivity/ml/stereo ===")
    conf = require_ok(embed("CCO", seed=42), "embed")

    fd = frontier_descriptors(conf.elements, conf.coords)
    fk = fukui_descriptors(conf.elements, conf.coords)
    print("fukui condensed f_plus", fk.condensed_f_plus)
    print("fukui validity notes", fk.validity_notes)

    rr = reactivity_ranking(conf.elements, conf.coords)
    print(rr)

    pka = empirical_pka("CC(=O)O")
    print(pka)

    md = ml_descriptors("CCO")
    print("ml_descriptors", md)
    print("MW", md.molecular_weight, "n_heavy", md.n_heavy_atoms, "fsp3", md.fsp3)

    pred = ml_predict("CCO")
    print("ml_predict", pred)

    sa2 = stereo_analysis("C/C=C/C")
    print("stereo C/C=C/C double bonds", sa2.double_bonds)
    sa3 = stereo_analysis("C(F)(Cl)(Br)I")
    print("stereocenters (topology)", sa3.n_stereocenters, sa3.stereocenters)
    conf_chiral = require_ok(embed("C(F)(Cl)(Br)I", seed=42), "embed chiral")
    sa4 = stereo_analysis("C(F)(Cl)(Br)I", coords=conf_chiral.coords)
    print("chiral R/S", sa4.n_stereocenters, [ (s.atom_index, s.configuration) for s in sa4.stereocenters])

    r = sssr("c1ccc2ccccc2c1")
    print("naphthalene rings", len(r.rings), "histogram", r.ring_size_histogram)

    fp1 = ecfp("c1ccccc1", radius=2, n_bits=2048)
    fp2 = ecfp("Cc1ccccc1", radius=2, n_bits=2048)
    print("tanimoto", tanimoto(fp1, fp2), "self", tanimoto(fp1, fp1))


if __name__ == "__main__":
    run()
