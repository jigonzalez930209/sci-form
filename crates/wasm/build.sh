#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# crates/wasm/build.sh — Canonical build script for the sci-form-wasm npm pkg.
#
# ALWAYS use this script instead of calling wasm-pack directly, because:
#   1. wasm-pack regenerates pkg/package.json from scratch on every run,
#      discarding the "snippets" entry needed by wasm-bindgen-rayon.
#   2. This script re-adds "snippets" and "README.md" to the files allowlist
#      and copies the root README so npm consumers get documentation.
#
# Usage:
#   ./build.sh                   # web (parallel) + nodejs targets
#   ./build.sh --web-only        # skip nodejs build
#   ./build.sh --profile dev     # dev profile (faster, unoptimised)
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PROFILE="release"
WEB_ONLY=false

for arg in "$@"; do
  case $arg in
    --web-only)  WEB_ONLY=true ;;
    --profile)   shift; PROFILE="$1" ;;
    --profile=*) PROFILE="${arg#*=}" ;;
  esac
done

PROFILE_FLAG="--$PROFILE"

# ── 1. Web target (browser / Vite) with Rayon parallelisation ────────────────
echo "→ Building WASM — web target (parallel enabled, profile: $PROFILE)..."
wasm-pack build --target web "$PROFILE_FLAG" --features parallel

# ── 2. Patch pkg/package.json ────────────────────────────────────────────────
# wasm-pack overwrites pkg/package.json each run.  We must re-add:
#   • "snippets"  — required at runtime by the wasm-bindgen-rayon import in
#                   sci_form_wasm.js; without it the file is omitted from the
#                   npm tarball and Vite throws "Failed to resolve import".
#   • "README.md" — so npm/pnpm consumers see documentation.
echo "→ Patching pkg/package.json..."
node -e "
const fs = require('fs');
const pkgPath = 'pkg/package.json';
const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));

['snippets', 'README.md'].forEach(f => {
  if (!pkg.files.includes(f)) pkg.files.push(f);
});
pkg.readme = 'README.md';

fs.writeFileSync(pkgPath, JSON.stringify(pkg, null, 2) + '\n');
console.log('  files:', JSON.stringify(pkg.files));
"

# ── 3. Copy root README into pkg/ ────────────────────────────────────────────
echo "→ Copying README.md → pkg/README.md..."
cp ../../README.md pkg/README.md

# ── 4. Node.js target (Jest / server-side) without parallelisation ───────────
if [ "$WEB_ONLY" = false ]; then
  echo "→ Building WASM — nodejs target (no parallelization, profile: $PROFILE)..."
  wasm-pack build --target nodejs "$PROFILE_FLAG" --out-dir pkg-node --no-default-features
fi

echo "✓ sci-form-wasm build complete."
echo "  Web pkg : $SCRIPT_DIR/pkg/"
if [ "$WEB_ONLY" = false ]; then
  echo "  Node pkg: $SCRIPT_DIR/pkg-node/"
fi
