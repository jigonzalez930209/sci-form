#!/usr/bin/env node
/**
 * post-wasm-build.mjs
 *
 * Run after `wasm-pack build` to install alpha/beta subpath exports
 * into both `pkg/` (ESM/browser) and `pkg-node/` (CJS/Node) output dirs.
 *
 * Usage:
 *   node scripts/post-wasm-build.mjs
 *
 * This script:
 *  1. Copies crates/wasm/js/alpha/ and crates/wasm/js/beta/ into pkg/ and pkg-node/
 *  2. Patches pkg/package.json and pkg-node/package.json with:
 *     - An `exports` map that includes `./alpha` and `./beta` subpaths
 *     - A `typesVersions` map for TypeScript resolution
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const root = path.resolve(__dirname, '..');

const JS_SRC = path.join(root, 'crates', 'wasm', 'js');
const PKG = path.join(root, 'pkg');
const PKG_NODE = path.join(root, 'pkg-node');

// ─── 1. Copy JS/TS re-export files into pkg directories ─────────────────────

function copySubpath(subdir, destRoot, mainEntry) {
  const src = path.join(JS_SRC, subdir);
  const dest = path.join(destRoot, subdir);
  if (!fs.existsSync(src)) return;
  fs.mkdirSync(dest, { recursive: true });

  for (const file of fs.readdirSync(src)) {
    let content = fs.readFileSync(path.join(src, file), 'utf8');
    // Adjust relative import path: in pkg/ the main module is at ../sci_form_wasm.js
    // in pkg-node/ it is at ../sci_form.js
    if (mainEntry !== '../sci_form_wasm.js') {
      content = content.replaceAll('../sci_form_wasm.js', mainEntry);
    }
    fs.writeFileSync(path.join(dest, file), content);
  }
  console.log(`  copied ${subdir}/ → ${path.relative(root, dest)}/`);
}

for (const [pkgDir, mainEntry] of [
  [PKG, '../sci_form_wasm.js'],
  [PKG_NODE, '../sci_form.js'],
]) {
  if (!fs.existsSync(pkgDir)) {
    console.warn(`  ⚠ ${path.relative(root, pkgDir)} not found — run wasm-pack build first`);
    continue;
  }
  console.log(`Patching ${path.relative(root, pkgDir)}/`);
  copySubpath('alpha', pkgDir, mainEntry);
  copySubpath('beta', pkgDir, mainEntry);
}

// ─── 2. Patch package.json exports fields ───────────────────────────────────

function patchPackageJson(pkgDir, mainJs, mainDts, isCjs) {
  const pkgJson = path.join(pkgDir, 'package.json');
  if (!fs.existsSync(pkgJson)) return;

  const pkg = JSON.parse(fs.readFileSync(pkgJson, 'utf8'));

  // Root export
  const rootExport = isCjs
    ? { require: `./${mainJs}`, types: `./${mainDts}` }
    : { import: `./${mainJs}`, types: `./${mainDts}` };

  const makeSubpathExport = (subdir, file) =>
    isCjs
      ? { require: `./${subdir}/${file}.js`, types: `./${subdir}/${file}.d.ts` }
      : { import: `./${subdir}/${file}.js`, types: `./${subdir}/${file}.d.ts` };

  pkg.exports = {
    '.': rootExport,
    './alpha': makeSubpathExport('alpha', 'index'),
    './beta': makeSubpathExport('beta', 'index'),
  };

  pkg.typesVersions = {
    '*': {
      alpha: [`./alpha/index.d.ts`],
      beta: [`./beta/index.d.ts`],
    },
  };

  fs.writeFileSync(pkgJson, JSON.stringify(pkg, null, 2) + '\n');
  console.log(`  patched ${path.relative(root, pkgJson)}`);
}

if (fs.existsSync(PKG)) {
  patchPackageJson(PKG, 'sci_form_wasm.js', 'sci_form_wasm.d.ts', false);
}
if (fs.existsSync(PKG_NODE)) {
  patchPackageJson(PKG_NODE, 'sci_form.js', 'sci_form.d.ts', true);
}

console.log('\n✓ alpha/beta subpath exports installed');
console.log('  Usage: import { alpha_compute_dft } from \'sci-form-wasm/alpha\'');
console.log('  Usage: import { beta_compute_kpm_dos } from \'sci-form-wasm/beta\'');
