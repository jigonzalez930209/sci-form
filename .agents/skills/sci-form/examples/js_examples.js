/**
 * sci-form WASM/JavaScript Examples — Complete API Reference
 * 
 * Works in: Browser (Webpack/Vite), Node.js (with WASM support)
 * Package:  npm install sci-form-wasm  (or build with wasm-pack)
 * 
 * ⚠️  CRITICAL for parallel operations:
 *     await initThreadPool(navigator.hardwareConcurrency)
 *     must be called ONCE before any parallel function.
 */

import init, {
  initThreadPool,
  version,
  // Geometry
  embed, embed_coords, embed_coords_typed, embed_batch, parse_smiles,
  // EHT
  eht_calculate, eht_support, eht_or_uff_fallback,
  eht_orbital_mesh, eht_orbital_grid_typed, eht_orbital_grid_from_coefficients_typed,
  // Charges & Surface
  compute_charges, compute_sasa, compute_dipole,
  compute_esp_grid_typed, compute_esp_grid_info,
  // Population & Bond Orders
  compute_population, compute_bond_orders,
  // Reactivity
  compute_frontier_descriptors, compute_fukui_descriptors,
  compute_reactivity_ranking, compute_empirical_pka,
  // UV-Vis
  compute_uv_vis_spectrum, compute_stda_uvvis,
  // IR
  compute_vibrational_analysis, compute_ir_spectrum,
  // NMR
  predict_nmr_shifts, predict_nmr_couplings, compute_nmr_spectrum, compute_hose_codes,
  // PM3 / xTB
  compute_pm3, compute_xtb,
  // ML
  compute_molecular_descriptors, compute_ml_properties,
  // Force fields
  compute_uff_energy, compute_uff_energy_with_aromatic_heuristics, compute_mmff94_energy,
  // Graph / Materials
  analyze_graph_features, compute_topology,
  create_unit_cell, assemble_framework,
  // Transport
  pack_batch_arrow, split_worker_tasks, estimate_workers,
  // System Planning
  system_capabilities, system_method_plan, compare_methods,
  // DOS / RMSD
  compute_dos, compute_rmsd,
  // Stereochemistry
  analyze_stereo,
  // Solvation
  compute_nonpolar_solvation, compute_gb_solvation,
  // Rings & Fingerprints
  compute_sssr, compute_ecfp, compute_tanimoto,
  // Clustering
  butina_cluster, compute_rmsd_matrix,
} from "sci-form-wasm";

// ─── JSON helpers ──────────────────────────────────────────────────────────────
const fromJSON = (s) => JSON.parse(s);
const toJSON   = (v) => JSON.stringify(v);

// ─── Test data ─────────────────────────────────────────────────────────────────

const ELEMENTS    = "[8,1,1]";          // Water O, H, H
const COORDS      = "[0,0,0, 0.96,0,0, -0.24,0.93,0]";
const BENZENE     = "c1ccccc1";
const ACETIC_ACID = "CC(=O)O";

// ─── Main async entrypoint ────────────────────────────────────────────────────

async function main() {

  // 1️⃣ Initialize WASM module (required first step)
  await init();

  // 2️⃣ Initialize thread pool for parallel processing (STRONGLY RECOMMENDED)
  //    This requires SharedArrayBuffer + COOP/COEP headers in the browser.
  const cores = (typeof navigator !== "undefined")
    ? navigator.hardwareConcurrency ?? 4
    : 4;
  await initThreadPool(cores);
  console.log(`sci-form ${version()} — ${cores} threads ready`);


  // ─── 1. Geometry / Embedding ───────────────────────────────────────────────

  // Full conformer result as JSON
  const embedResult = fromJSON(embed(BENZENE, 42));
  console.log("Benzene:", embedResult.num_atoms, "atoms");
  if (embedResult.error) throw new Error(embedResult.error);

  // Just coords as JSON
  const embedCoords = fromJSON(embed_coords(BENZENE, 42));
  // { coords: [...], num_atoms: N }

  // 🚀 Zero-copy typed array (no JSON overhead, ideal for GPU/WebGL)
  const coordsTA = embed_coords_typed(BENZENE, 42); // Float64Array
  console.log("Benzene coords typed array length:", coordsTA.length); // 3 * N

  // Batch embed (uses Rayon thread pool initialized above)
  const smilesList = "C1=CC=CC=C1\nO\nCCO\nCC(=O)O\nCN1C=NC2=C1C(=O)N(C(=O)N2C)C";
  const batchJSON = fromJSON(embed_batch(smilesList, 42));
  const ok = batchJSON.filter(r => !r.error).length;
  console.log(`Batch: ${ok}/${batchJSON.length} succeeded`);

  // Parse (topology only)
  const parsed = fromJSON(parse_smiles("CC(=O)N"));
  console.log("Acetamide:", parsed.num_atoms, "atoms,", parsed.num_bonds, "bonds");


  // ─── 2. Force Fields ───────────────────────────────────────────────────────

  // Get benzene coords first
  const benzCoords = toJSON(fromJSON(embed_coords(BENZENE, 42)).coords);

  const uffResult = fromJSON(compute_uff_energy(BENZENE, benzCoords));
  console.log("UFF:", uffResult.energy.toFixed(3), uffResult.unit);

  const uffHeur = fromJSON(compute_uff_energy_with_aromatic_heuristics(BENZENE, benzCoords));
  console.log("UFF corrected:", uffHeur.corrected_energy_kcal_mol.toFixed(3), "kcal/mol",
    "| aromatic bonds:", uffHeur.aromatic_bond_count);

  const mmff = fromJSON(compute_mmff94_energy(BENZENE, benzCoords));
  console.log("MMFF94:", mmff.energy.toFixed(3), mmff.unit);


  // ─── 3. EHT / Quantum ──────────────────────────────────────────────────────

  // Pre-flight check
  const support = fromJSON(eht_support(ELEMENTS));
  console.log("EHT support:", support.level, "| TM:", support.has_transition_metals);

  // Run EHT
  const eht = fromJSON(eht_calculate(ELEMENTS, COORDS, 1.75));
  if (eht.error) throw new Error(eht.error);
  console.log(`H2O EHT: HOMO=${eht.homo_energy.toFixed(3)} eV, gap=${eht.gap.toFixed(3)} eV`);

  // Smart routing: EHT if supported, else UFF
  const workflow = fromJSON(eht_or_uff_fallback("O", ELEMENTS, COORDS, false));
  console.log("Workflow mode:", workflow.mode);

  // Orbital mesh (vertices/normals for WebGL)
  const mesh = fromJSON(eht_orbital_mesh(ELEMENTS, COORDS, eht.homo_index, 0.2, 0.02));
  console.log("HOMO mesh:", mesh.num_triangles, "triangles");

  // 🚀 Zero-copy orbital grid (Float32Array for GPU volume rendering)
  const orbGrid = eht_orbital_grid_typed(ELEMENTS, COORDS, eht.homo_index, 0.2);
  console.log("Orbital grid Float32Array length:", orbGrid.length);

  // Orbital grid from precomputed coefficients
  const coeffJSON = JSON.stringify(eht.coefficients);
  const orbGrid2 = eht_orbital_grid_from_coefficients_typed(ELEMENTS, COORDS, coeffJSON, eht.homo_index, 0.2);

  // Density of States
  const dos = fromJSON(compute_dos(ELEMENTS, COORDS, 0.3, -30.0, 5.0, 500));
  console.log("DOS: sigma=", dos.sigma, "| points:", dos.energies.length);


  // ─── 4. Charges & Surface ─────────────────────────────────────────────────

  const charges = fromJSON(compute_charges(ACETIC_ACID));
  console.log("Gasteiger charges:", charges.charges.map(c => c.toFixed(3)));
  console.log("  Total:", charges.total_charge.toFixed(4), "Iterations:", charges.iterations);

  const sasa = fromJSON(compute_sasa(ELEMENTS, COORDS, 1.4));
  console.log("SASA:", sasa.total_sasa.toFixed(2), "Ų");

  const dipole = fromJSON(compute_dipole(ELEMENTS, COORDS));
  console.log("Dipole:", dipole.magnitude.toFixed(3), dipole.unit, "| vector:", dipole.vector);

  // 🚀 ESP grid zero-copy (huge datasets)
  const espValues = compute_esp_grid_typed(ELEMENTS, COORDS, 0.5, 3.0); // Float64Array
  const espInfo   = fromJSON(compute_esp_grid_info(ELEMENTS, COORDS, 0.5, 3.0));
  console.log("ESP grid:", espInfo.dims, "=", espValues.length, "values");


  // ─── 5. Population & Bond Orders ──────────────────────────────────────────

  const pop = fromJSON(compute_population(ELEMENTS, COORDS));
  console.log("Mulliken:", pop.mulliken_charges.map(c => c.toFixed(3)));
  console.log("Löwdin:  ", pop.lowdin_charges.map(c => c.toFixed(3)));

  const bo = fromJSON(compute_bond_orders(ELEMENTS, COORDS));
  bo.bonds.slice(0, 3).forEach(b =>
    console.log(`Bond ${b.atom_i}-${b.atom_j}: Wiberg=${b.wiberg.toFixed(3)} Mayer=${b.mayer.toFixed(3)}`)
  );


  // ─── 6. Reactivity ────────────────────────────────────────────────────────

  // Get phenol coordinates first
  const phenol = fromJSON(embed("Oc1ccccc1", 42));
  const phEl   = JSON.stringify(phenol.elements);
  const phCo   = JSON.stringify(phenol.coords);

  const frontier = fromJSON(compute_frontier_descriptors(phEl, phCo));
  console.log("Frontier gap:", frontier.gap.toFixed(3), "eV");
  console.log("HOMO contributions:", frontier.homo_atom_contributions.slice(0, 4).map(c => c.toFixed(3)));

  const fukui = fromJSON(compute_fukui_descriptors(phEl, phCo));
  console.log("Fukui validity:", fukui.validity_notes);
  fukui.condensed.slice(0, 3).forEach(row =>
    console.log(`  Atom ${row.atom_index}: f+=${row.f_plus.toFixed(3)} f-=${row.f_minus.toFixed(3)}`)
  );

  const ranking = fromJSON(compute_reactivity_ranking(phEl, phCo));
  if (ranking.nucleophilic_attack_sites[0])
    console.log("Top nucleophilic:", ranking.nucleophilic_attack_sites[0]);

  const pka = fromJSON(compute_empirical_pka(ACETIC_ACID));
  pka.acidic_sites.forEach(s =>
    console.log(`pKa acidic: ${s.estimated_pka.toFixed(1)} at atom ${s.atom_index} (${s.environment})`)
  );


  // ─── 7. UV-Vis Spectroscopy ───────────────────────────────────────────────

  const uvvis = fromJSON(compute_uv_vis_spectrum(ELEMENTS, COORDS, 0.25, 0.5, 8.0, 600));
  console.log("UV-Vis (EHT):", uvvis.peaks.length, "peaks");
  uvvis.peaks.slice(0, 3).forEach(p =>
    console.log(`  ${p.wavelength_nm.toFixed(1)} nm  I=${p.intensity.toFixed(4)}  MO:${p.from_mo}->${p.to_mo}`)
  );

  const stda = fromJSON(compute_stda_uvvis(ELEMENTS, COORDS, 0.3, 1.0, 8.0, 500, "gaussian"));
  console.log("sTDA excitations:", stda.excitations.length);
  stda.excitations.slice(0, 3).forEach(e =>
    console.log(`  ${e.wavelength_nm.toFixed(1)} nm  f=${e.oscillator_strength.toFixed(4)}`)
  );


  // ─── 8. IR Spectroscopy ───────────────────────────────────────────────────

  // Step 1: Vibrational analysis (Numerical Hessian)
  const vibJSON = compute_vibrational_analysis(ELEMENTS, COORDS, "eht", 0.005);
  const vib = fromJSON(vibJSON);
  console.log("Vibrational modes:", vib.modes.length, "| ZPVE:", vib.zpve_ev.toFixed(4), "eV");
  vib.modes.slice(0, 5).forEach(m =>
    console.log(`  ${m.frequency_cm1.toFixed(1)} cm⁻¹  IR=${m.ir_intensity.toFixed(4)}  real=${m.is_real}`)
  );

  // Step 2: IR spectrum from the analysis JSON (chained call pattern)
  const ir = fromJSON(compute_ir_spectrum(vibJSON, 15.0, 400.0, 4000.0, 1000));
  console.log("IR spectrum:", ir.wavenumbers.length, "points,", ir.peaks.length, "peaks");


  // ─── 9. NMR Spectroscopy ──────────────────────────────────────────────────

  const shifts = fromJSON(predict_nmr_shifts(ACETIC_ACID));
  console.log("¹H NMR shifts:", shifts.h_shifts.map(s => `${s.shift_ppm.toFixed(2)} ppm`));
  console.log("¹³C NMR shifts:", shifts.c_shifts.map(s => `${s.shift_ppm.toFixed(1)} ppm`));

  const couplings = fromJSON(predict_nmr_couplings("CCC", "[]"));
  couplings.forEach(c =>
    console.log(`J(${c.h1_index},${c.h2_index}): ${c.j_hz.toFixed(1)} Hz  ${c.n_bonds}J  ${c.coupling_type}`)
  );

  const nmr1h = fromJSON(compute_nmr_spectrum(ACETIC_ACID, "1H", 0.02, 0.0, 12.0, 1000));
  console.log("¹H spectrum:", nmr1h.ppm_axis.length, "points,", nmr1h.peaks.length, "peaks");

  const nmr13c = fromJSON(compute_nmr_spectrum(ACETIC_ACID, "13C", 0.5, 0.0, 220.0, 1000));
  console.log("¹³C spectrum:", nmr13c.peaks.length, "peaks");

  const hose = fromJSON(compute_hose_codes(BENZENE, 2));
  hose.slice(0, 3).forEach(h => console.log(`HOSE atom ${h.atom_index}: ${h.full_code}`));


  // ─── 10. PM3 / xTB ────────────────────────────────────────────────────────

  const pm3 = fromJSON(compute_pm3(ELEMENTS, COORDS));
  console.log(`PM3: Hf=${pm3.heat_of_formation.toFixed(3)} kcal/mol  gap=${pm3.gap.toFixed(3)} eV  converged=${pm3.converged}`);

  const xtb = fromJSON(compute_xtb(ELEMENTS, COORDS));
  console.log(`xTB: E=${xtb.total_energy.toFixed(6)} eV  gap=${xtb.gap.toFixed(3)} eV  converged=${xtb.converged}`);


  // ─── 11. ML Properties ────────────────────────────────────────────────────

  const desc = fromJSON(compute_molecular_descriptors("Cc1ccccc1"));
  console.log(`Toluene: MW=${desc.molecular_weight.toFixed(2)}  HBA=${desc.n_hba}  HBD=${desc.n_hbd}  Fsp3=${desc.fsp3.toFixed(2)}`);

  const ml = fromJSON(compute_ml_properties("Cc1ccccc1"));
  console.log(`  LogP=${ml.logp.toFixed(2)}  logS=${ml.log_solubility.toFixed(2)}  Lipinski=${ml.lipinski_passes}  DL=${ml.druglikeness.toFixed(2)}`);


  // ─── 12. Topology & Graph ─────────────────────────────────────────────────

  const gf = fromJSON(analyze_graph_features(BENZENE));
  const nArom = gf.aromaticity.aromatic_atoms.filter(Boolean).length;
  console.log("Benzene: aromatic atoms:", nArom);

  const topo = fromJSON(compute_topology(ELEMENTS, COORDS));
  topo.metal_centers.forEach(c =>
    console.log(`Metal Z=${c.element}: CN=${c.coordination_number}  geom=${c.geometry}`)
  );


  // ─── 13. Materials / Crystallography ──────────────────────────────────────

  const cell = fromJSON(create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0));
  console.log("Unit cell: a=", cell.a, "Å  vol=", cell.volume.toFixed(2), "ų");

  // Assemble MOF (PCU topology, Zn octahedral, 10 Å cell)
  const framework = fromJSON(assemble_framework("pcu", 30, "octahedral", 10.0, 1));
  console.log("MOF atoms:", framework.num_atoms);

  // Supercell 2×2×2
  const superFW = fromJSON(assemble_framework("pcu", 30, "octahedral", 10.0, 2));
  console.log("Supercell atoms:", superFW.num_atoms);


  // ─── 14. Transport & Streaming ────────────────────────────────────────────

  // Arrow columnar format for streaming large batches
  const batchResults = fromJSON(embed_batch(smilesList, 42));
  const arrow = fromJSON(pack_batch_arrow(JSON.stringify(batchResults)));
  console.log("Arrow batch:", arrow.num_rows, "rows,", arrow.num_columns, "cols");

  // Split into worker tasks for Web Worker dispatch
  const bigSmiles = JSON.stringify(["O","CCO","c1ccccc1","CC(=O)O","CCC","CCCC","CCCCC","CCCCCC"]);
  const nWorkers = estimate_workers(8, 4);
  const tasks = fromJSON(split_worker_tasks(bigSmiles, nWorkers, 42));
  console.log(`Transport: ${nWorkers} workers, ${tasks.length} tasks`);


  // ─── 15. RMSD ─────────────────────────────────────────────────────────────

  const r1 = JSON.stringify(fromJSON(embed_coords("C", 1)).coords);
  const r2 = JSON.stringify(fromJSON(embed_coords("C", 2)).coords);
  const rmsd = fromJSON(compute_rmsd(r1, r2));
  console.log("Methane RMSD:", rmsd.rmsd.toFixed(4), "Å");


  // ─── 16. System Planning ──────────────────────────────────────────────────

  const caps = fromJSON(system_capabilities(ELEMENTS));
  console.log("H2O capabilities:", {
    embed: caps.embed.available,
    uff:   caps.uff.available,
    eht:   caps.eht.available,
  });

  const plan = fromJSON(system_method_plan(ELEMENTS));
  console.log("Recommended for orbitals:", plan.orbitals.recommended);
  console.log("Rationale:", plan.orbitals.rationale);

  const cmp = fromJSON(compare_methods("O", ELEMENTS, COORDS, false));
  cmp.comparisons.forEach(e =>
    console.log(`  ${e.method}: status=${e.status}  available=${e.available}`)
  );


  // ─── 17. Stereochemistry ───────────────────────────────────────────

  // R/S chirality with 3D coords
  const chiralResult = fromJSON(embed("C(F)(Cl)(Br)I", 42));
  const chiralElem   = JSON.stringify(chiralResult.elements);
  const chiralCoords = JSON.stringify(chiralResult.coords);
  const stereo = fromJSON(analyze_stereo("C(F)(Cl)(Br)I", chiralCoords));
  console.log(`Stereocenters: ${stereo.n_stereocenters}  double bonds: ${stereo.n_double_bonds}`);
  stereo.stereocenters.forEach(sc =>
    console.log(`  atom ${sc.atom_index}: config=${sc.configuration}  priorities=${JSON.stringify(sc.priorities)}`)
  );

  // E/Z with 3D coords
  const alkResult = fromJSON(embed("CC=CC", 42));
  const ezStereo = fromJSON(analyze_stereo("CC=CC", JSON.stringify(alkResult.coords)));
  ezStereo.double_bonds.forEach(db =>
    console.log(`  double bond ${db.atom1}-${db.atom2}: config=${db.configuration}`)
  );

  // Topology-only (empty coords)
  const stereoTopo = fromJSON(analyze_stereo("C(F)(Cl)(Br)I", "[]"));
  console.log("Topology stereocenters:", stereoTopo.n_stereocenters);


  // ─── 18. Solvation ────────────────────────────────────────────────

  // Non-polar solvation (SASA ASP model)
  const npSolv = fromJSON(compute_nonpolar_solvation(ELEMENTS, COORDS, 1.4));
  console.log(`Non-polar solvation: ${npSolv.energy_kcal_mol.toFixed(4)} kcal/mol  SASA=${npSolv.total_sasa.toFixed(2)} Å²`);
  npSolv.atom_contributions.forEach((c, i) =>
    console.log(`  atom ${i}: ΔG=${c.toFixed(4)} kcal/mol  sasa=${npSolv.atom_sasa[i].toFixed(2)} Å²`)
  );

  // Generalized Born (HCT model)
  const chargesW = fromJSON(compute_charges("O"));
  const chargesJSON = JSON.stringify(chargesW.charges);
  const gb = fromJSON(compute_gb_solvation(ELEMENTS, COORDS, chargesJSON, 78.5, 1.0, 1.4));
  console.log(`GB electrostatic: ${gb.electrostatic_energy_kcal_mol.toFixed(4)} kcal/mol`);
  console.log(`GB nonpolar:      ${gb.nonpolar_energy_kcal_mol.toFixed(4)} kcal/mol`);
  console.log(`GB total:         ${gb.total_energy_kcal_mol.toFixed(4)} kcal/mol`);
  console.log(`Born radii: ${gb.born_radii.map(r => r.toFixed(3))}`);


  // ─── 19. Rings & Fingerprints ─────────────────────────────────────────

  // SSSR — Smallest Set of Smallest Rings
  const benzRings = fromJSON(compute_sssr(BENZENE));
  console.log(`Benzene rings: ${benzRings.rings.length}  histogram: ${JSON.stringify(benzRings.ring_size_histogram)}`);
  benzRings.rings.forEach(r =>
    console.log(`  ring size ${r.size}  atoms=${JSON.stringify(r.atoms)}  aromatic=${r.is_aromatic}`)
  );

  const naphRings = fromJSON(compute_sssr("c1ccc2ccccc2c1")); // naphthalene
  console.log(`Naphthalene: ${naphRings.rings.length} rings`);

  // ECFP fingerprints
  const fpBenz = fromJSON(compute_ecfp(BENZENE, 2, 2048));
  const fpTol  = fromJSON(compute_ecfp("Cc1ccccc1", 2, 2048));
  const fpAcid = fromJSON(compute_ecfp(ACETIC_ACID, 2, 2048));
  console.log(`Benzene fp: ${fpBenz.n_bits} bits, ${fpBenz.on_bits.length} on`);

  // Tanimoto similarity
  const simBT = fromJSON(compute_tanimoto(
    JSON.stringify(fpBenz), JSON.stringify(fpTol)
  ));
  const simBA = fromJSON(compute_tanimoto(
    JSON.stringify(fpBenz), JSON.stringify(fpAcid)
  ));
  console.log(`Tanimoto benzene/toluene: ${simBT.tanimoto.toFixed(4)}`);
  console.log(`Tanimoto benzene/acetic:  ${simBA.tanimoto.toFixed(4)}`);


  // ─── 20. Clustering ─────────────────────────────────────────────────

  const smilesList4Cluster = ["c1ccccc1", "Cc1ccccc1", "Clc1ccccc1", "CC(=O)O", "CCO"];
  const conformerCoords = smilesList4Cluster.map(s => fromJSON(embed(s, 42)).coords);
  const conformersJSON = JSON.stringify(conformerCoords);

  // Butina RMSD clustering
  const clustResult = fromJSON(butina_cluster(conformersJSON, 2.0));
  console.log(`Clusters (cutoff=2.0Å): ${clustResult.n_clusters}`);
  console.log(`  assignments: ${JSON.stringify(clustResult.assignments)}`);
  console.log(`  centroids:   ${JSON.stringify(clustResult.centroid_indices)}`);
  console.log(`  sizes:       ${JSON.stringify(clustResult.cluster_sizes)}`);

  // Pair-wise RMSD matrix
  const rmsdMat = fromJSON(compute_rmsd_matrix(conformersJSON));
  console.log(`RMSD matrix: ${rmsdMat.length}×${rmsdMat[0].length}`);
  rmsdMat.forEach((row, i) =>
    console.log(`  row ${i}: ${row.map(v => v.toFixed(3))}`))
  );
}

main().catch(console.error);


// ─── Web Worker Pattern (parallel batch in browser) ───────────────────────────
//
// main.js:
//   import init, { initThreadPool, embed_batch } from "sci-form-wasm";
//   await init();
//   await initThreadPool(navigator.hardwareConcurrency);
//   const results = JSON.parse(embed_batch(smiles.join("\n"), 42));
//
// NOTE: Requires server headers:
//   Cross-Origin-Opener-Policy: same-origin
//   Cross-Origin-Embedder-Policy: require-corp
//
// Vite config:
//   server: { headers: { "Cross-Origin-Opener-Policy": "same-origin",
//                        "Cross-Origin-Embedder-Policy": "require-corp" } }
