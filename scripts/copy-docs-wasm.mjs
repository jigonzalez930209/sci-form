#!/usr/bin/env node
// Copies the single-threaded WASM build used by the live examples
// (crates/wasm/pkg-docs, from `crates/wasm/build.sh --docs-only`) into the
// documentation site (docs/public/wasm). Without it the site still builds and
// the examples show a note instead of a result.
import { cpSync, existsSync, mkdirSync, rmSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const src = path.join(root, 'crates/wasm/pkg-docs');
const dst = path.join(root, 'docs/public/wasm');
if (!existsSync(path.join(src, 'sci_form_wasm_bg.wasm'))) {
  console.error('crates/wasm/pkg-docs not found: run `bash crates/wasm/build.sh --docs-only` first');
  process.exit(1);
}
rmSync(dst, { recursive: true, force: true });
mkdirSync(dst, { recursive: true });
cpSync(src, dst, { recursive: true });
console.log('copied', src, '->', dst);
