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
import {
  alpha_compute_aevs,
  alpha_compute_dft,
  alpha_compute_reaxff_gradient,
} from "sci-form-wasm/alpha";
import {
  beta_compute_cpm_charges,
  beta_compute_kpm_dos,
  beta_solve_eht_randnla,
} from "sci-form-wasm/beta";

const fromJSON = (value) => JSON.parse(value);
const toJSON = (value) => JSON.stringify(value);

async function main() {
  await init();

  if (typeof SharedArrayBuffer !== "undefined") {
    const cores = typeof navigator === "undefined" ? 4 : navigator.hardwareConcurrency ?? 4;
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

  const dft = fromJSON(alpha_compute_dft(elements, coords, "pbe"));
  console.log("DFT gap:", dft.gap.toFixed(3), "eV", "converged:", dft.converged);

  const reaxff = fromJSON(alpha_compute_reaxff_gradient(elements, coords));
  console.log("ReaxFF energy:", reaxff.energy_kcal_mol.toFixed(3), "gradient length:", reaxff.gradient.length);

  const aevs = fromJSON(alpha_compute_aevs(elements, coords, "{}"));
  console.log("AEV atoms:", aevs.n_atoms, "descriptor length:", aevs.aev_length);

  const kpm = fromJSON(beta_compute_kpm_dos(elements, coords, toJSON({ order: 128, n_points: 256 })));
  console.log("KPM points:", kpm.energies.length, "order:", kpm.order);

  const randnla = fromJSON(beta_solve_eht_randnla(elements, coords, toJSON({ sketch_size: 16 })));
  console.log("RandNLA orbitals:", randnla.orbital_energies.length, "residual:", randnla.residual_error);

  const cpm = fromJSON(beta_compute_cpm_charges(elements, coords, 0.0));
  console.log("CPM total charge:", cpm.total_charge.toFixed(4));
}

main().catch((error) => {
  console.error(error);
  process.exitCode = 1;
});
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
