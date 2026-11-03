#!/usr/bin/env node
// Runs every live example of the documentation against the COMPILED WASM
// package (crates/wasm/pkg-docs, built by `crates/wasm/build.sh --docs-only`).
// An example passes when it runs without throwing and does not return
// {"error": ...}. It also reports the public functions that have no example.
//
//   node scripts/test-live-examples.mjs [--list-uncovered]

import { readFileSync } from 'node:fs';
import { fileURLToPath, pathToFileURL } from 'node:url';
import path from 'node:path';
import { makeScope, runExample } from '../docs/examples/runtime.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const pkg = path.join(root, 'crates/wasm/pkg-docs');
const catalog = JSON.parse(readFileSync(path.join(root, 'docs/examples/catalog.json'), 'utf8'));

const mod = await import(pathToFileURL(path.join(pkg, 'sci_form_wasm.js')).href);
await mod.default({ module_or_path: readFileSync(path.join(pkg, 'sci_form_wasm_bg.wasm')) });
const scope = makeScope(mod);

let failed = 0;
for (const ex of catalog) {
  try {
    const { value, ms } = await runExample(scope, ex.code);
    if (value && typeof value === 'object' && !Array.isArray(value) && 'error' in value && value.error) {
      throw new Error(String(value.error));
    }
    if (value === undefined) throw new Error('the example returned nothing');
    console.log(`ok   ${ex.id.padEnd(46)} ${ms.toFixed(0).padStart(6)} ms`);
  } catch (e) {
    failed++;
    console.log(`FAIL ${ex.id}: ${e.message}`);
  }
}

// Coverage: every public function should have an example, except the
// variants that only change the seed width, the output container, or that are
// internal batching helpers.
const skip = /(_seed64$|_typed$|_wasm$|^pack_|^split_worker|^eht_orbital_grid|^embed_coords|^compute_esp_grid$|^init|^start|^alpha_|^beta_|^__|^LiveSimulation|^Wasm|^wbg|^compute_md_trajectory_nvt_)/;
const have = new Set(catalog.map((e) => e.id));
const exports = Object.keys(mod).filter((k) => typeof mod[k] === 'function' && !skip.test(k) && k !== 'default');
const uncovered = exports.filter((k) => !have.has(k));
console.log(`\n${catalog.length} examples, ${failed} failed; ${uncovered.length} public functions without an example`);
if (process.argv.includes('--list-uncovered') || uncovered.length) console.log('uncovered:', uncovered.join(', '));
process.exit(failed ? 1 : 0);
