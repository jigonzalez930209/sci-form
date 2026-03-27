#!/usr/bin/env python3

from sci_form import (
    charges,
    dipole,
    ecfp,
    embed,
    embed_batch,
    nonpolar_solvation,
    parse,
    pm3_calculate,
    population,
    sssr,
    stereo_analysis,
    tanimoto,
    xtb_calculate,
)
from sci_form.alpha import alpha_compute_aevs, dft_calculate, reaxff_gradient
from sci_form.beta import cpm_charges, eht_randnla, kpm_dos


def require_ok(result):
    if not result.is_ok():
        raise RuntimeError(result.error)
    return result


def main():
    ethanol = require_ok(embed("CCO", seed=42))
    print(f"ethanol atoms: {ethanol.num_atoms}")

    batch = embed_batch(["O", "CCO", "c1ccccc1"], seed=42, num_threads=0)
    successes = sum(1 for item in batch if item.is_ok())
    print(f"batch success: {successes}/{len(batch)}")

    parsed = parse("CC(=O)O")
    print(f"acetic acid topology: {parsed['num_atoms']} atoms, {parsed['num_bonds']} bonds")

    charge_result = charges("CC(=O)O")
    print(f"Gasteiger total charge: {charge_result.total_charge:.4f}")

    population_result = population(ethanol.elements, ethanol.coords)
    print(f"Mulliken total charge: {population_result.total_charge_mulliken:.4f}")

    dipole_result = dipole(ethanol.elements, ethanol.coords)
    print(f"dipole magnitude: {dipole_result.magnitude:.3f} {dipole_result.unit}")

    pm3_result = pm3_calculate(ethanol.elements, ethanol.coords)
    print(f"PM3 gap: {pm3_result.gap:.3f} eV, converged: {pm3_result.converged}")

    xtb_result = xtb_calculate(ethanol.elements, ethanol.coords)
    print(f"xTB gap: {xtb_result.gap:.3f} eV, converged: {xtb_result.converged}")

    chiral = require_ok(embed("C(F)(Cl)(Br)I", seed=42))
    stereo = stereo_analysis("C(F)(Cl)(Br)I", coords=chiral.coords)
    print(f"stereocenters: {stereo.n_stereocenters}")

    solvation = nonpolar_solvation(ethanol.elements, ethanol.coords)
    print(f"non-polar solvation: {solvation.energy_kcal_mol:.3f} kcal/mol")

    rings = sssr("c1ccccc1")
    print(f"benzene rings: {len(rings.rings)}")

    benzene = ecfp("c1ccccc1", radius=2, n_bits=2048)
    toluene = ecfp("Cc1ccccc1", radius=2, n_bits=2048)
    print(f"benzene/toluene tanimoto: {tanimoto(benzene, toluene):.3f}")

    dft_result = dft_calculate(ethanol.elements, ethanol.coords, method="pbe")
    print(f"DFT gap: {dft_result.gap:.3f} eV, converged: {dft_result.converged}")

    reaxff = reaxff_gradient(ethanol.elements, ethanol.coords)
    print(f"ReaxFF energy: {reaxff.energy_kcal_mol:.3f} kcal/mol, gradient length: {len(reaxff.gradient)}")

    aevs = alpha_compute_aevs(ethanol.elements, ethanol.coords)
    print(f"AEV atoms: {aevs.n_atoms}, descriptor length: {aevs.aev_length}")

    kpm = kpm_dos(ethanol.elements, ethanol.coords, order=128)
    print(f"KPM points: {len(kpm.energies)}, order: {kpm.order}")

    randnla = eht_randnla(ethanol.elements, ethanol.coords, sketch_size=16)
    print(f"RandNLA orbitals: {len(randnla.orbital_energies)}, residual: {randnla.residual_error:.3e}")

    cpm = cpm_charges(ethanol.elements, ethanol.coords, potential=0.0)
    print(f"CPM total charge: {cpm.total_charge:.4f}")


if __name__ == "__main__":
    main()
