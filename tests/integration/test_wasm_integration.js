#!/usr/bin/env node
/**
 * sci-form WASM/TypeScript Integration Tests
 *
 * Tests that WASM bindings produce correct results across all API functions.
 *
 * Requires a built WASM package:
 *   npm run wasm:build:node     (emits crates/wasm/pkg-node/sci_form_wasm.js)
 *   npm run wasm:build          (emits crates/wasm/pkg/sci_form_wasm.js, ESM)
 */

const fs = require('node:fs');
const path = require('node:path');

const rootDir = path.resolve(__dirname, '..', '..');
const CANDIDATES = [
  path.join(rootDir, 'crates', 'wasm', 'pkg-node', 'sci_form_wasm.js'),
  path.join(rootDir, 'crates', 'wasm', 'pkg', 'sci_form_wasm.js'),
];

const wasmPath = CANDIDATES.find((candidate) => fs.existsSync(candidate));
if (!wasmPath) {
  console.error(
    'FATAL: no built WASM package found.\n' +
      `Looked in:\n${CANDIDATES.map((c) => `  - ${c}`).join('\n')}\n` +
      'Build it first with: npm run wasm:build:node'
  );
  process.exit(1);
}
console.log(`Using WASM package: ${path.relative(rootDir, wasmPath)}`);
const wasm = require(wasmPath);

let passed = 0;
let failed = 0;

function assert(condition, msg) {
  if (!condition) throw new Error(msg);
}

function test(name, fn) {
  try {
    fn();
    passed++;
    console.log(`  [PASS] ${name}`);
  } catch (e) {
    failed++;
    console.log(`  [FAIL] ${name}: ${e.message}`);
  }
}

console.log('='.repeat(60));
console.log('sci-form WASM/TypeScript Integration Tests');
console.log('='.repeat(60));

test('Version', () => {
  const v = wasm.version();
  assert(v.startsWith('sci-form'), `Bad version: ${v}`);
  console.log(`    version: ${v}`);
});

test('Single embed - ethanol', () => {
  const r = JSON.parse(wasm.embed('CCO', 42));
  assert(!r.error, `Embed failed: ${r.error}`);
  assert(r.num_atoms === 9, `Expected 9 atoms, got ${r.num_atoms}`);
  assert(r.coords.length === 27, `Expected 27 coords, got ${r.coords.length}`);
  assert(r.elements.length === 9, 'Expected 9 elements');
  assert(r.elements[0] === 6, 'First atom should be Carbon');
  assert(r.elements[2] === 8, 'Third atom should be Oxygen');
  assert(r.bonds.length > 0, 'Should have bonds');
  console.log(`    CCO: ${r.num_atoms} atoms, ${r.bonds.length} bonds`);
});

test('Single embed - benzene', () => {
  const r = JSON.parse(wasm.embed('c1ccccc1', 42));
  assert(!r.error, `Embed failed: ${r.error}`);
  assert(r.num_atoms === 12, `Expected 12 atoms, got ${r.num_atoms}`);
  console.log(`    benzene: ${r.num_atoms} atoms`);
});

test('Single embed - aspirin', () => {
  const r = JSON.parse(wasm.embed('CC(=O)Oc1ccccc1C(=O)O', 42));
  assert(!r.error, `Embed failed: ${r.error}`);
  assert(r.num_atoms === 21, `Expected 21 atoms, got ${r.num_atoms}`);
  console.log(`    aspirin: ${r.num_atoms} atoms`);
});

test('Reproducibility', () => {
  const r1 = JSON.parse(wasm.embed('CCO', 42));
  const r2 = JSON.parse(wasm.embed('CCO', 42));
  assert(JSON.stringify(r1.coords) === JSON.stringify(r2.coords), 'Same seed should give same result');
  console.log('    same seed = same coordinates: OK');
});

test('Multiple molecules', () => {
  const mols = {
    'C': 5, 'N': 4, 'O': 3, 'C#N': 3, 'C=C': 6,
    'c1ccncc1': 11, 'CC(=O)O': 8,
  };
  for (const [smi, expected] of Object.entries(mols)) {
    const r = JSON.parse(wasm.embed(smi, 42));
    assert(!r.error, `Failed for ${smi}: ${r.error}`);
    assert(r.num_atoms === expected, `${smi}: expected ${expected} atoms, got ${r.num_atoms}`);
  }
  console.log(`    ${Object.keys(mols).length} molecules: all OK`);
});

test('Embed coords (compact)', () => {
  const r = JSON.parse(wasm.embed_coords('CCO', 42));
  assert(!r.error, `embed_coords failed: ${r.error}`);
  assert(r.num_atoms === 9, `Expected 9 atoms`);
  assert(r.coords.length === 27, `Expected 27 coords`);
  console.log(`    embed_coords: ${r.num_atoms} atoms, ${r.coords.length} coord values`);
});

test('Batch embed', () => {
  const results = JSON.parse(wasm.embed_batch('CCO\nc1ccccc1\nCC(=O)O', 42));
  assert(results.length === 3, `Expected 3 results, got ${results.length}`);
  results.forEach(r => {
    assert(!r.error, `Batch failed for ${r.smiles}: ${r.error}`);
  });

  // Verify batch matches single results
  const single = JSON.parse(wasm.embed('CCO', 42));
  const batchCCO = results.find(r => r.smiles === 'CCO');
  assert(JSON.stringify(single.coords) === JSON.stringify(batchCCO.coords),
    'Batch/single mismatch for CCO');
  console.log(`    batch 3 molecules: all match single results`);
});

test('Parse SMILES', () => {
  const r = JSON.parse(wasm.parse_smiles('CCO'));
  assert(r.num_atoms === 9, `Expected 9 atoms, got ${r.num_atoms}`);
  assert(r.num_bonds === 8, `Expected 8 bonds, got ${r.num_bonds}`);
  console.log(`    parse CCO: ${r.num_atoms} atoms, ${r.num_bonds} bonds`);
});

test('Error handling', () => {
  const r = JSON.parse(wasm.embed('INVALID_XYZ', 42));
  assert(r.error !== null && r.error !== undefined, 'Should have error for invalid SMILES');
  console.log(`    invalid SMILES error: "${r.error}"`);
});

test('Coordinate sanity', () => {
  const r = JSON.parse(wasm.embed('CCO', 42));
  for (const [a, b, order] of r.bonds) {
    const dx = r.coords[a * 3] - r.coords[b * 3];
    const dy = r.coords[a * 3 + 1] - r.coords[b * 3 + 1];
    const dz = r.coords[a * 3 + 2] - r.coords[b * 3 + 2];
    const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
    assert(dist > 0.5 && dist < 3.0, `Bond ${a}-${b} length ${dist.toFixed(2)} out of range`);
  }
  console.log(`    all bond lengths in [0.5, 3.0] Å`);
});

test('Data roundtrip JSON', () => {
  const r = JSON.parse(wasm.embed('CCO', 42));
  const json = JSON.stringify(r);
  const back = JSON.parse(json);
  assert(back.num_atoms === r.num_atoms, 'num_atoms mismatch after roundtrip');
  assert(JSON.stringify(back.coords) === JSON.stringify(r.coords), 'coords mismatch');
  console.log('    JSON roundtrip: OK');
});

// C-H4: a zero dielectric constant used to divide by zero and return a
// `NaN` energy; the binding must report it as a JSON error instead.
test('GB solvation rejects a zero dielectric constant', () => {
  const elements = JSON.stringify([8, 1, 1]);
  const coords = JSON.stringify([0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692]);
  const charges = JSON.stringify([0.5, -0.25, -0.25]);
  const r = JSON.parse(
    wasm.compute_gb_solvation(elements, coords, charges, 0.0, 1.0, 1.4),
  );
  assert(r.error !== null && r.error !== undefined, 'Should have an error for epsilon = 0');
  assert(!String(r.error).includes('NaN'), `The error must not be a NaN: ${r.error}`);
  console.log(`    zero dielectric error: "${r.error}"`);
});

// C-M1: a charge vector shorter than the atom list used to index past its
// end inside the GB kernel.
test('GB solvation rejects a charge vector of the wrong length', () => {
  const elements = JSON.stringify([8, 1, 1]);
  const coords = JSON.stringify([0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692]);
  const charges = JSON.stringify([0.5]);
  const r = JSON.parse(
    wasm.compute_gb_solvation(elements, coords, charges, 78.5, 1.0, 1.4),
  );
  assert(r.error !== null && r.error !== undefined, 'Should have an error for 1 charge / 3 atoms');
  console.log(`    charge length error: "${r.error}"`);
});

// C-H2: the capabilities that used to be Rust-only.
test('TPSA and the ensemble prediction', () => {
  const tpsa = JSON.parse(wasm.compute_tpsa('CCO'));
  assert(Math.abs(tpsa.tpsa - 20.23) < 0.01, `TPSA = ${tpsa.tpsa}`);
  const ensemble = JSON.parse(wasm.predict_ensemble('CCO'));
  assert(typeof ensemble.logp === 'number', 'missing logp');
  assert(ensemble.confidence >= 0 && ensemble.confidence <= 1, 'confidence out of range');
  console.log(`    TPSA = ${tpsa.tpsa} Å², logP = ${ensemble.logp.toFixed(3)}`);
});

test('Open-shell SCF and the analytical gradients', () => {
  const ch3 = JSON.stringify([6, 1, 1, 1]);
  const ch3Coords = JSON.stringify([
    0, 0, 0,
    1.09, 0, 0,
    -0.36, 1.03, 0,
    -0.36, -0.52, 0.89,
  ]);
  const uhf = JSON.parse(wasm.compute_uhf(ch3, ch3Coords, 0, 2));
  assert(!uhf.error, `UHF failed: ${uhf.error}`);
  assert(uhf.converged === true, 'the UHF SCF did not converge');
  assert(uhf.total_energy < 0, `UHF energy = ${uhf.total_energy}`);
  const rohf = JSON.parse(wasm.compute_rohf(ch3, ch3Coords, 0, 2));
  assert(!rohf.error, `ROHF failed: ${rohf.error}`);
  const water = JSON.stringify([8, 1, 1]);
  const waterCoords = JSON.stringify([0, 0, 0.1173, 0, 0.7572, -0.4692, 0, -0.7572, -0.4692]);
  const grad = JSON.parse(wasm.compute_pm3_gradient(water, waterCoords));
  assert(!grad.error, `PM3 gradient failed: ${grad.error}`);
  assert(grad.gradients.length === 3, 'one gradient per atom');
  console.log(`    UHF E = ${uhf.total_energy.toFixed(4)} Ha, ROHF E = ${rohf.total_energy.toFixed(4)} Ha`);
});

test('CIF import and export', () => {
  const cif = [
    'data_test',
    '_cell_length_a 5.0',
    '_cell_length_b 5.0',
    '_cell_length_c 5.0',
    '_cell_angle_alpha 90.0',
    '_cell_angle_beta 90.0',
    '_cell_angle_gamma 90.0',
    'loop_',
    '_atom_site_label',
    '_atom_site_type_symbol',
    '_atom_site_fract_x',
    '_atom_site_fract_y',
    '_atom_site_fract_z',
    '_atom_site_occupancy',
    'Na1 Na 0.0 0.0 0.0 1.0',
  ].join('\n');
  const structure = JSON.parse(wasm.parse_cif(cif));
  assert(!structure.error, `parse_cif failed: ${structure.error}`);
  assert(structure.data_block === 'test', `data block = ${structure.data_block}`);
  assert(structure.atom_sites.length === 1, 'one atom site');
  const back = wasm.write_cif(JSON.stringify(structure));
  assert(back.includes('data_test'), 'the written CIF must carry the data block');
  console.log(`    CIF roundtrip: ${structure.atom_sites[0].label}`);
});

console.log('='.repeat(60));
console.log(`Results: ${passed} passed, ${failed} failed`);
console.log('='.repeat(60));
process.exit(failed > 0 ? 1 : 0);
