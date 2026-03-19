#!/usr/bin/env python3
"""
sci-form Python Examples — Complete API Reference
Covers all 50+ functions across all API groups.
Build with: maturin develop --features parallel
Install   : pip install sci-form  (if published)
"""

import sci_form

# ─── Helpers ──────────────────────────────────────────────────────────────────

WATER_ELEMENTS = [8, 1, 1]
WATER_COORDS   = [0.0, 0.0, 0.0, 0.96, 0.0, 0.0, -0.24, 0.93, 0.0]
BENZENE_SMILES = "c1ccccc1"
ACETIC_ACID    = "CC(=O)O"

# ─── 1. Geometry / Embedding ──────────────────────────────────────────────────

# Single conformer
result = sci_form.embed("C1=CC=CC=C1", seed=42)
assert result.is_ok(), f"Failed: {result.error}"
print(f"Benzene: {result.num_atoms} atoms, time={result.time_ms:.1f} ms")
positions = result.get_positions()  # -> list of (x, y, z) tuples
print(f"  First atom: {positions[0]}")

# Parallel batch — num_threads=0 uses all available cores
smiles_list = ["O", "CCO", "c1ccccc1", "CC(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"]
batch = sci_form.embed_batch(smiles_list, seed=42, num_threads=0)
ok = sum(1 for r in batch if r.is_ok())
print(f"Batch: {ok}/{len(batch)} succeeded")

# Topology parse (no 3D)
info = sci_form.parse("CC(=O)N")
print(f"Acetamide: {info['num_atoms']} atoms, {info['num_bonds']} bonds")
for atom in info['atoms']:
    print(f"  Z={atom['element']} hyb={atom['hybridization']} q={atom['formal_charge']}")

# RMSD / Kabsch alignment
r1 = sci_form.embed("C", seed=1)
r2 = sci_form.embed("C", seed=2)
aligned = sci_form.rmsd(r1.coords, r2.coords)
print(f"Methane RMSD: {aligned.rmsd:.4f} Å")

# ─── 2. Force Fields ──────────────────────────────────────────────────────────

benz = sci_form.embed(BENZENE_SMILES, seed=42)
e_uff = sci_form.uff_energy(BENZENE_SMILES, benz.coords)
print(f"UFF benzene: {e_uff:.3f} kcal/mol")

e_mmff = sci_form.mmff94_energy(BENZENE_SMILES, benz.coords)
print(f"MMFF94 benzene: {e_mmff:.3f} kcal/mol")

e_heur = sci_form.uff_energy_with_aromatic_heuristics(BENZENE_SMILES, benz.coords)
print(f"UFF corrected: {e_heur.corrected_energy_kcal_mol:.3f} kcal/mol "
      f"(aromatic bonds: {e_heur.aromatic_bond_count})")
print(f"  Notes: {e_heur.notes}")

# ─── 3. EHT / Quantum ─────────────────────────────────────────────────────────

eht = sci_form.eht_calculate(WATER_ELEMENTS, WATER_COORDS, k=1.75)
print(f"H2O EHT: {eht}")
print(f"  HOMO={eht.homo_energy:.3f} eV  LUMO={eht.lumo_energy:.3f} eV  gap={eht.gap:.3f} eV")
print(f"  Support: {eht.support_level}  Warnings: {eht.warnings}")

# EHT or UFF fallback (smart routing)
result_eff = sci_form.eht_or_uff_fallback("O", WATER_ELEMENTS, WATER_COORDS, allow_experimental_eht=False)
print(f"EHT/UFF fallback mode: {result_eff['mode']}")

# Orbital isosurface mesh (HOMO by default)
mesh = sci_form.eht_orbital_mesh(WATER_ELEMENTS, WATER_COORDS,
                                  mo_index=eht.homo_index, spacing=0.2, isovalue=0.02)
print(f"HOMO mesh: {mesh['num_triangles']} triangles")

# Density of States
dos = sci_form.dos(WATER_ELEMENTS, WATER_COORDS, sigma=0.3, e_min=-30.0, e_max=5.0, n_points=500)
print(f"DOS: {len(dos.energies)} energy points, sigma={dos.sigma} eV")

# EHT support query (pre-flight check)
support = sci_form.eht_support(WATER_ELEMENTS)
print(f"EHT support: {support.level}, TM={support.has_transition_metals}")

# ─── 4. Charges & Surface ─────────────────────────────────────────────────────

charges = sci_form.charges(ACETIC_ACID)
print(f"Gasteiger charges: {[f'{c:.3f}' for c in charges.charges]}")
print(f"  Total: {charges.total_charge:.4f}  Iterations: {charges.iterations}")

sasa = sci_form.sasa(WATER_ELEMENTS, WATER_COORDS, probe_radius=1.4)
print(f"SASA: {sasa.total_sasa:.2f} Ų  (per-atom: {[f'{a:.2f}' for a in sasa.atom_sasa]})")

dipole = sci_form.dipole(WATER_ELEMENTS, WATER_COORDS)
print(f"Dipole: {dipole.magnitude:.3f} {dipole.unit}  vector={dipole.vector}")

# ─── 5. Population & Bond Orders ──────────────────────────────────────────────

pop = sci_form.population(WATER_ELEMENTS, WATER_COORDS)
print(f"Mulliken: {[f'{c:.3f}' for c in pop.mulliken_charges]}")
print(f"Löwdin:   {[f'{c:.3f}' for c in pop.lowdin_charges]}")

bo = sci_form.bond_orders(WATER_ELEMENTS, WATER_COORDS)
print(f"Bond orders: {len(bo.atom_pairs)} pairs")
for (a, b), w, m in zip(bo.atom_pairs, bo.wiberg, bo.mayer):
    print(f"  {a}-{b}: Wiberg={w:.3f}  Mayer={m:.3f}")

# ─── 6. Reactivity ────────────────────────────────────────────────────────────

phenol = sci_form.embed("Oc1ccccc1", seed=42)
ph_el = phenol.elements
ph_co = phenol.coords

frontier = sci_form.frontier_descriptors(ph_el, ph_co)
print(f"Frontier: gap={frontier.gap:.3f} eV")
print(f"  HOMO contributions: {[f'{c:.3f}' for c in frontier.homo_atom_contributions[:4]]}")

fukui = sci_form.fukui_descriptors(ph_el, ph_co)
print(f"Fukui: gap={fukui.gap:.3f} eV  validity={fukui.validity_notes}")
for i, (fp, fm) in enumerate(zip(fukui.condensed_f_plus, fukui.condensed_f_minus)):
    print(f"  Atom {fukui.condensed_atom_indices[i]}: f+={fp:.3f}  f-={fm:.3f}")

ranking = sci_form.reactivity_ranking(ph_el, ph_co)
top_nuc = ranking.nucleophilic_attack_sites[0] if ranking.nucleophilic_attack_sites else None
top_elec = ranking.electrophilic_attack_sites[0] if ranking.electrophilic_attack_sites else None
if top_nuc:  print(f"Top nucleophilic site: atom {top_nuc.atom_index} (score={top_nuc.score:.3f})")
if top_elec: print(f"Top electrophilic site: atom {top_elec.atom_index} (score={top_elec.score:.3f})")

pka = sci_form.empirical_pka(ACETIC_ACID)
for site in pka.acidic_sites:
    print(f"Acidic pKa: {site.estimated_pka:.1f} at atom {site.atom_index} ({site.environment})")
for site in pka.basic_sites:
    print(f"Basic pKa: {site.estimated_pka:.1f} at atom {site.atom_index}")

# ─── 7. UV-Vis Spectroscopy ───────────────────────────────────────────────────

uvvis = sci_form.uvvis_spectrum(WATER_ELEMENTS, WATER_COORDS,
                                sigma=0.25, e_min=0.5, e_max=8.0, n_points=600)
print(f"UV-Vis (EHT): {len(uvvis.energies_ev)} points, {len(uvvis.peaks)} peaks")
for peak in uvvis.peaks[:3]:
    print(f"  {peak.wavelength_nm:.1f} nm ({peak.energy_ev:.2f} eV)  I={peak.intensity:.4f}  MO:{peak.from_mo}->{peak.to_mo}")

stda = sci_form.stda_uvvis(WATER_ELEMENTS, WATER_COORDS,
                            sigma=0.3, e_min=1.0, e_max=8.0, n_points=500, broadening="gaussian")
print(f"sTDA: {len(stda.excitations)} excitations")
for exc in stda.excitations[:3]:
    print(f"  {exc.wavelength_nm:.1f} nm  f={exc.oscillator_strength:.4f}  μ={exc.transition_dipole:.4f}")

# ─── 8. IR Spectroscopy ───────────────────────────────────────────────────────

vib = sci_form.vibrational_analysis(WATER_ELEMENTS, WATER_COORDS,
                                     method="eht", step_size=0.005)
print(f"Vibrational: {vib.n_atoms} atoms, {len(vib.modes)} modes, ZPVE={vib.zpve_ev:.4f} eV")
for mode in vib.modes[:5]:
    print(f"  {mode.frequency_cm1:.1f} cm⁻¹  IR={mode.ir_intensity:.4f}  real={mode.is_real}")

ir = sci_form.ir_spectrum(WATER_ELEMENTS, WATER_COORDS,
                           method="eht", gamma=15.0, wn_min=400.0, wn_max=4000.0, n_points=1000)
print(f"IR: {len(ir.wavenumbers)} points, {len(ir.peaks)} peaks")
for peak in ir.peaks[:3]:
    print(f"  Peak at {peak.frequency_cm1:.1f} cm⁻¹  I={peak.ir_intensity:.4f}")

# ─── 9. NMR Spectroscopy ──────────────────────────────────────────────────────

shifts = sci_form.nmr_shifts(ACETIC_ACID)
print(f"¹H shifts for {ACETIC_ACID}:")
for s in shifts.h_shifts:
    print(f"  atom {s.atom_index}: {s.shift_ppm:.2f} ppm ({s.environment})  conf={s.confidence:.2f}")
print(f"¹³C shifts:")
for s in shifts.c_shifts:
    print(f"  atom {s.atom_index}: {s.shift_ppm:.1f} ppm ({s.environment})")

couplings = sci_form.nmr_couplings("CCC", coords=[])
for c in couplings:
    print(f"J({c.h1_index},{c.h2_index}): {c.j_hz:.1f} Hz  {c.n_bonds}J  type={c.coupling_type}")

nmr_1h = sci_form.nmr_spectrum(ACETIC_ACID, nucleus="1H", gamma=0.02, ppm_min=0.0, ppm_max=12.0, n_points=1000)
print(f"¹H spectrum: {len(nmr_1h.ppm_axis)} points, {len(nmr_1h.peaks)} peaks")
for p in nmr_1h.peaks:
    print(f"  {p.shift_ppm:.2f} ppm  mult={p.multiplicity}  env={p.environment}")

nmr_13c = sci_form.nmr_spectrum(ACETIC_ACID, nucleus="13C", gamma=0.5, ppm_min=0.0, ppm_max=220.0, n_points=1000)
print(f"¹³C spectrum: {len(nmr_13c.peaks)} peaks")

hose = sci_form.hose_codes("c1ccccc1", max_radius=2)
for atom_idx, element, code in hose[:3]:
    print(f"HOSE atom {atom_idx} (Z={element}): {code}")

# ─── 10. PM3 and xTB ──────────────────────────────────────────────────────────

pm3 = sci_form.pm3_calculate(WATER_ELEMENTS, WATER_COORDS)
print(f"PM3: Hf={pm3.heat_of_formation:.3f} kcal/mol  gap={pm3.gap:.3f} eV  "
      f"SCF={pm3.scf_iterations} iters  converged={pm3.converged}")

xtb = sci_form.xtb_calculate(WATER_ELEMENTS, WATER_COORDS)
print(f"xTB: E={xtb.total_energy:.6f} eV  gap={xtb.gap:.3f} eV  "
      f"SCC={xtb.scf_iterations} iters  converged={xtb.converged}")

# ─── 11. ML Properties ────────────────────────────────────────────────────────

desc = sci_form.ml_descriptors("Cc1ccccc1")  # Toluene
print(f"Toluene descriptors:")
print(f"  MW={desc.molecular_weight:.2f}  HBA={desc.n_hba}  HBD={desc.n_hbd}  "
      f"Fsp3={desc.fsp3:.2f}  nRings={desc.n_rings}")

ml = sci_form.ml_predict("Cc1ccccc1")
print(f"  LogP={ml.logp:.2f}  MR={ml.molar_refractivity:.2f}  "
      f"logS={ml.log_solubility:.2f}  Lipinski={ml.lipinski_passes}  "
      f"druglikeness={ml.druglikeness:.2f}")

# ─── 12. Topology & Graph ─────────────────────────────────────────────────────

gf = sci_form.graph_features(BENZENE_SMILES)
n_arom = sum(gf.aromatic_atoms)
print(f"Benzene: {n_arom} aromatic atoms, {len(gf.tagged_tetrahedral_centers)} tagged stereocenters")

# Topology analysis (metal coordination)
porphyrin_metal = [26, 7, 7, 7, 7]  # Fe + 4 N
porphyrin_coords = [0.0,0.0,0.0, 2.0,0.0,0.0, -2.0,0.0,0.0, 0.0,2.0,0.0, 0.0,-2.0,0.0]
topo = sci_form.topology_analysis(porphyrin_metal, porphyrin_coords)
for center in topo.metal_centers:
    print(f"Metal Z={center.element}: CN={center.coordination_number}  geom={center.geometry}  score={center.geometry_score:.3f}")

# ─── 13. Materials / Crystallography ──────────────────────────────────────────

cell = sci_form.unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0)
print(f"Cubic cell: a={cell.a}Å  vol={cell.volume:.2f} ų")

mof = sci_form.assemble_framework(topology="pcu", metal=30, geometry="octahedral",
                                   lattice_a=10.0, supercell=1)
print(f"MOF: {mof.num_atoms} atoms, labels={mof.labels[:3]}")

# 2x2x2 Supercell
mof_super = sci_form.assemble_framework(topology="pcu", metal=30, geometry="octahedral",
                                        lattice_a=10.0, supercell=2)
print(f"Supercell: {mof_super.num_atoms} atoms")

# ─── 14. Transport / Streaming ────────────────────────────────────────────────

all_smiles = ["O", "CCO", "c1ccccc1"] * 100
n_workers = sci_form.estimate_workers(len(all_smiles), max_workers=8)
tasks = sci_form.split_worker_tasks(all_smiles, n_workers=n_workers, seed=42)
print(f"Transport: {n_workers} workers, {len(tasks)} tasks")

batch_results = sci_form.embed_batch(all_smiles, seed=42, num_threads=0)
record = sci_form.pack_conformers(batch_results)
print(f"RecordBatch: {record.num_rows} rows, {record.num_columns} cols, {record.byte_size} bytes")

# ─── 15. System Planning ──────────────────────────────────────────────────────

caps = sci_form.system_capabilities(WATER_ELEMENTS)
print(f"H2O capabilities:")
print(f"  embed={caps.embed.available}  uff={caps.uff.available}  eht={caps.eht.available}")

plan = sci_form.system_method_plan(WATER_ELEMENTS)
print(f"Recommended for orbitals: {plan.orbitals.recommended}")
print(f"Fallback for orbitals:    {plan.orbitals.fallback}")
print(f"Rationale: {plan.orbitals.rationale}")

cmp = sci_form.compare_methods("O", WATER_ELEMENTS, WATER_COORDS, allow_experimental_eht=False)
for entry in cmp.comparisons:
    print(f"  {entry.method}: status={entry.status}  available={entry.available}")
