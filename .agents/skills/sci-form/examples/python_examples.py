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


def require_ok(result):
    if not result.is_ok():
        raise RuntimeError(result.error)
    return result


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

# ─── 16. Stereochemistry ──────────────────────────────────────────────────────

# Tetrahedral R/S — with 3D coords
chiral = sci_form.embed("C(F)(Cl)(Br)I", seed=42)
stereo = sci_form.stereo_analysis("C(F)(Cl)(Br)I", coords=chiral.coords)
print(f"Stereocenters: {stereo.n_stereocenters}  double bonds: {stereo.n_double_bonds}")
for sc in stereo.stereocenters:
    print(f"  atom {sc.atom_index} ({chr(sc.element+64 if sc.element else 63)}): "
          f"config={sc.configuration}  priorities={sc.priorities}")

# E/Z double bond — with 3D coords
alk = sci_form.embed("CC=CC", seed=42)
ez = sci_form.stereo_analysis("CC=CC", coords=alk.coords)
for db in ez.double_bonds:
    print(f"  double bond {db.atom1}-{db.atom2}: config={db.configuration}")

# Topology-only mode (no 3D, limited assignment)
stereo_topo = sci_form.stereo_analysis("C(F)(Cl)(Br)I")
print(f"Topology-only stereocenters identified: {stereo_topo.n_stereocenters}")

# ─── 17. Solvation ────────────────────────────────────────────────────────────

water = sci_form.embed("O", seed=42)

# Non-polar (SASA-based ASP model)
np_solv = sci_form.nonpolar_solvation(water.elements, water.coords)
print(f"Non-polar solvation: {np_solv.energy_kcal_mol:.4f} kcal/mol  SASA={np_solv.total_sasa:.2f} Å²")
for i, (contrib, sasa_i) in enumerate(zip(np_solv.atom_contributions, np_solv.atom_sasa)):
    print(f"  atom {i}: ΔG={contrib:.4f} kcal/mol  sasa={sasa_i:.2f} Å²")

# Generalized Born (HCT model)
ch = sci_form.charges("O")
gb = sci_form.gb_solvation(water.elements, water.coords, ch.charges,
                            solvent_dielectric=78.5, solute_dielectric=1.0, probe_radius=1.4)
print(f"GB electrostatic: {gb.electrostatic_energy_kcal_mol:.4f} kcal/mol")
print(f"GB nonpolar:      {gb.nonpolar_energy_kcal_mol:.4f} kcal/mol")
print(f"GB total:         {gb.total_energy_kcal_mol:.4f} kcal/mol")
print(f"Born radii: {[f'{r:.3f}' for r in gb.born_radii]}")

# ─── 18. Rings & Fingerprints ─────────────────────────────────────────────────

# SSSR — Smallest Set of Smallest Rings
benz_sssr = sci_form.sssr(BENZENE_SMILES)
print(f"Benzene rings: {len(benz_sssr.rings)}, histogram: {benz_sssr.ring_size_histogram}")
for ring in benz_sssr.rings:
    print(f"  ring size {ring.size}, atoms={ring.atoms}, aromatic={ring.is_aromatic}")

naph_sssr = sci_form.sssr("c1ccc2ccccc2c1")  # naphthalene
print(f"Naphthalene: {len(naph_sssr.rings)} rings, ring_count/atom: {naph_sssr.atom_ring_count}")

# ECFP fingerprints
fp_benz = sci_form.ecfp(BENZENE_SMILES, radius=2, n_bits=2048)
fp_tol  = sci_form.ecfp("Cc1ccccc1", radius=2, n_bits=2048)
fp_acid = sci_form.ecfp(ACETIC_ACID, radius=2, n_bits=2048)
print(f"Benzene fp: {fp_benz.n_bits} bits, {len(fp_benz.on_bits)} on, density={fp_benz.density:.4f}")

sim_bt = sci_form.tanimoto(fp_benz, fp_tol)
sim_ba = sci_form.tanimoto(fp_benz, fp_acid)
print(f"Tanimoto benzene/toluene: {sim_bt:.4f}")
print(f"Tanimoto benzene/acetic:  {sim_ba:.4f}")

# ECFP4 vs ECFP2
fp4 = sci_form.ecfp(BENZENE_SMILES, radius=2, n_bits=2048)
fp2 = sci_form.ecfp(BENZENE_SMILES, radius=1, n_bits=2048)
print(f"ECFP4 on-bits: {len(fp4.on_bits)}  ECFP2 on-bits: {len(fp2.on_bits)}")

# ─── 19. Clustering ───────────────────────────────────────────────────────────

smiles_set = ["c1ccccc1", "Cc1ccccc1", "Clc1ccccc1", "CC(=O)O", "CCO"]
conformers = [sci_form.embed(s, seed=42).coords for s in smiles_set]

# Butina RMSD clustering
clusters = sci_form.butina_cluster(conformers, rmsd_cutoff=2.0)
print(f"Clusters (cutoff=2.0 Å): {clusters.n_clusters}")
print(f"  assignments: {clusters.assignments}")
print(f"  centroid indices: {clusters.centroid_indices}")
print(f"  cluster sizes: {clusters.cluster_sizes}")

# Full RMSD matrix
mat = sci_form.rmsd_matrix(conformers)
print(f"RMSD matrix shape: {len(mat)}×{len(mat[0])}")
for i, row in enumerate(mat):
    print(f"  row {i}: {[f'{v:.3f}' for v in row]}")

# Diverse subset selection
diverse_idx = sci_form.filter_diverse(conformers, rmsd_cutoff=1.0)
print(f"Diverse indices (cutoff=1.0 Å): {diverse_idx}")

# ─── Extra: Ensemble NMR J-couplings ──────────────────────────────────────────

conf_list = [sci_form.embed("CCC", seed=s).coords for s in range(5)]
energies  = [0.0, 0.3, 0.6, 1.2, 2.5]  # kcal/mol relative
ens_j = sci_form.nmr_ensemble_couplings("CCC", conf_list, energies, temperature_k=298.15)
print(f"Ensemble J-couplings: {len(ens_j)} pairs")
for c in ens_j:
    print(f"  J({c.h1_index},{c.h2_index}) = {c.j_hz:.2f} Hz  {c.n_bonds}J  type={c.coupling_type}")

# ─── Extra: IR Gaussian broadened ─────────────────────────────────────────────

vib3 = sci_form.vibrational_analysis(WATER_ELEMENTS, WATER_COORDS, method="eht")
ir_g = sci_form.ir_spectrum_broadened(vib3, gamma=20.0, wn_min=400.0, wn_max=4000.0,
                                       n_points=1000, broadening="gaussian")
print(f"IR Gaussian: {len(ir_g.wavenumbers)} points, {len(ir_g.peaks)} peaks")

# ─── Extra: Charges with custom config ────────────────────────────────────────

ch_cfg = sci_form.charges_configured(ACETIC_ACID, max_iter=10, initial_damping=0.4)
print(f"Configured charges: {[f'{c:.3f}' for c in ch_cfg.charges]}")
