import init, {
  analyze_stereo,
  compute_charges,
  compute_dipole,
  compute_ecfp,
  compute_nonpolar_solvation,
  compute_pm3,
  compute_population,
  compute_sssr,
  compute_tanimoto,
  compute_xtb,
  create_unit_cell,
  embed,
  embed_batch,
  embed_coords_typed,
  initThreadPool,
  parse_smiles,
  version,
} from "sci-form-wasm";

const fromJSON = (value) => JSON.parse(value);
const toJSON = (value) => JSON.stringify(value);

async function main() {
  await init();

  if (typeof SharedArrayBuffer !== "undefined") {
    const cores = typeof navigator === "undefined"
      ? 4
      : navigator.hardwareConcurrency ?? 4;
    await initThreadPool(cores);
  }

  console.log(`sci-form-wasm ${version()}`);

  const ethanol = fromJSON(embed("CCO", 42));
  if (ethanol.error) {
    throw new Error(ethanol.error);
  }
  console.log("ethanol atoms:", ethanol.num_atoms);

  const coordsTyped = embed_coords_typed("CCO", 42);
  console.log("typed coordinate length:", coordsTyped.length);

  const batch = fromJSON(embed_batch("O\nCCO\nc1ccccc1", 42));
  const successes = batch.filter((entry) => !entry.error).length;
  console.log(`batch success: ${successes}/${batch.length}`);

  const parsed = fromJSON(parse_smiles("CC(=O)O"));
  console.log("acetic acid topology:", parsed.num_atoms, "atoms", parsed.num_bonds, "bonds");

  const chargeResult = fromJSON(compute_charges("CC(=O)O"));
  console.log("Gasteiger total charge:", chargeResult.total_charge.toFixed(4));

  const elements = toJSON(ethanol.elements);
  const coords = toJSON(ethanol.coords);

  const population = fromJSON(compute_population(elements, coords));
  console.log("Mulliken total charge:", population.total_charge_mulliken.toFixed(4));

  const dipole = fromJSON(compute_dipole(elements, coords));
  console.log("dipole magnitude:", dipole.magnitude.toFixed(3), dipole.unit);

  const pm3 = fromJSON(compute_pm3(elements, coords));
  console.log("PM3 gap:", pm3.gap.toFixed(3), "eV", "converged:", pm3.converged);

  const xtb = fromJSON(compute_xtb(elements, coords));
  console.log("xTB gap:", xtb.gap.toFixed(3), "eV", "converged:", xtb.converged);

  const chiral = fromJSON(embed("C(F)(Cl)(Br)I", 42));
  if (chiral.error) {
    throw new Error(chiral.error);
  }
  const stereo = fromJSON(analyze_stereo("C(F)(Cl)(Br)I", toJSON(chiral.coords)));
  console.log("stereocenters:", stereo.n_stereocenters);

  const solvation = fromJSON(compute_nonpolar_solvation(elements, coords, 1.4));
  console.log("non-polar solvation:", solvation.energy_kcal_mol.toFixed(3), "kcal/mol");

  const rings = fromJSON(compute_sssr("c1ccccc1"));
  console.log("benzene rings:", rings.rings.length);

  const benzene = compute_ecfp("c1ccccc1", 2, 2048);
  const toluene = compute_ecfp("Cc1ccccc1", 2, 2048);
  const similarity = fromJSON(compute_tanimoto(benzene, toluene));
  console.log("benzene/toluene tanimoto:", similarity.tanimoto.toFixed(3));

  const cell = fromJSON(create_unit_cell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0));
  console.log("cubic cell volume:", cell.volume.toFixed(2));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});

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
