#!/usr/bin/env node
/**
 * sci-form WASM/TypeScript Integration Tests
 *
 * Tests that WASM bindings produce correct results across all API functions.
 */

const wasm = require('../crates/wasm/pkg/sci_form_wasm.js');

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

console.log('='.repeat(60));
console.log(`Results: ${passed} passed, ${failed} failed`);
console.log('='.repeat(60));
process.exit(failed > 0 ? 1 : 0);
