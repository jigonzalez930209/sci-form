#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# crates/wasm/build.sh — Canonical build script for the sci-form-wasm npm pkg.
#
# ALWAYS use this script instead of calling wasm-pack directly, because:
#   1. wasm-pack regenerates pkg/package.json from scratch on every run,
#      discarding the "snippets" entry needed by wasm-bindgen-rayon.
#   2. This script re-adds "snippets" and "README.md" to the files allowlist
#      and copies the root README so npm consumers get documentation.
#   3. It installs the ./alpha and ./beta subpath packages in BOTH output
#      directories (ESM for pkg/, CommonJS for pkg-node/).
#
# The npm package is built twice:
#   pkg/      — ESM, "type": "module", for browsers (Vite, webpack, ...)
#   pkg-node/ — CommonJS, for Node.js (require(...))
#
# The *published* package is `pkg/`, and it is dual: its `exports` map sends
# `import` to the ESM half at the root and `require` to the CommonJS half in
# `pkg/cjs/`.  Both halves expose the same API (the `parallel` and
# `experimental-gpu` features change behaviour, not the export list); the
# CommonJS half simply has no rayon threads and no WebGPU, because Node.js
# has neither.
#
# Usage:
#   ./build.sh                      # web (parallel) + nodejs targets
#   ./build.sh --web-only           # skip the nodejs build
#   ./build.sh --node-only          # skip the web build
#   ./build.sh --profile dev        # dev profile (faster, unoptimised)
#
# Options:
#   --web-features "<space separated feature list>"
#   --node-features "<space separated feature list>"
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

PROFILE="release"
WEB_ONLY=false
NODE_ONLY=false
WEB_FEATURES="parallel"
NODE_FEATURES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --web-only)
      WEB_ONLY=true
      shift
      ;;
    --node-only)
      NODE_ONLY=true
      shift
      ;;
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --profile=*)
      PROFILE="${1#*=}"
      shift
      ;;
    --web-features)
      WEB_FEATURES="$2"
      shift 2
      ;;
    --web-features=*)
      WEB_FEATURES="${1#*=}"
      shift
      ;;
    --node-features)
      NODE_FEATURES="$2"
      shift 2
      ;;
    --node-features=*)
      NODE_FEATURES="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ "$WEB_ONLY" == true && "$NODE_ONLY" == true ]]; then
  echo "error: --web-only and --node-only are mutually exclusive" >&2
  exit 1
fi

PROFILE_FLAG="--$PROFILE"

# ─────────────────────────────────────────────────────────────────────────────
# Patch helpers.
#
# wasm-pack rewrites <out-dir>/package.json on every run, so the "snippets"
# entry (required at runtime by the wasm-bindgen-rayon import) and the
# alpha/beta subpath exports have to be re-installed afterwards.
#
#   pkg/      — ESM:  "type": "module", `import ... from 'sci-form-wasm/alpha'`
#   pkg-node/ — CJS:  `require('sci-form-wasm/alpha')`
# ─────────────────────────────────────────────────────────────────────────────
patch_package_json() {
  local dir="$1" module_kind="$2"
  node -e "
const fs = require('fs');
const path = require('path');
const pkgPath = process.argv[1];
const kind = process.argv[2];
const pkg = JSON.parse(fs.readFileSync(pkgPath, 'utf8'));

// 'snippets' holds the wasm-bindgen-rayon worker glue: without it the
// tarball is missing the file that sci_form_wasm.js imports at runtime and
// the bundler fails with 'Failed to resolve import'.
for (const f of ['snippets', 'README.md', 'alpha', 'beta']) {
  if (!pkg.files.includes(f)) pkg.files.push(f);
}
pkg.readme = 'README.md';

if (kind === 'esm') {
  pkg.type = 'module';
  // The condition order matters: Node matches the first satisfied
  // condition, and "default" always matches.  "default" therefore has to
  // come last, otherwise it shadows "require" and every CommonJS consumer
  // would load the ESM half.
  //
  // Dual package: when the CommonJS half was assembled into pkg/cjs/ the
  // "require" condition points at it, so require('sci-form-wasm') and
  // import ... from 'sci-form-wasm' both work from the same package.
  const dual = fs.existsSync(path.join(path.dirname(pkgPath), 'cjs'));
  if (dual && !pkg.files.includes('cjs')) pkg.files.push('cjs');
  const entry = (esm, cjs) => ({
    types: esm + '.d.ts',
    import: esm + '.js',
    ...(dual ? { require: cjs + '.js' } : {}),
    default: esm + '.js'
  });
  pkg.exports = {
    '.': entry('./sci_form_wasm', './cjs/sci_form_wasm'),
    './alpha': entry('./alpha/index', './cjs/alpha/index'),
    './beta': entry('./beta/index', './cjs/beta/index')
  };
} else {
  // CommonJS: the nodejs target of wasm-pack emits a CommonJS glue file.
  // Writing an ESM export statement here would be a syntax error for
  // require(...).
  pkg.type = 'commonjs';
  pkg.main = './sci_form_wasm.js';
  pkg.exports = {
    '.': {
      types: './sci_form_wasm.d.ts',
      require: './sci_form_wasm.js',
      default: './sci_form_wasm.js'
    },
    './alpha': {
      types: './alpha/index.d.ts',
      require: './alpha/index.js',
      default: './alpha/index.js'
    },
    './beta': {
      types: './beta/index.d.ts',
      require: './beta/index.js',
      default: './beta/index.js'
    }
  };
}

// TypeScript resolution for the subpaths (TS < 4.7 does not read 'exports').
pkg.typesVersions = {
  '*': {
    alpha: ['./alpha/index.d.ts'],
    beta: ['./beta/index.d.ts']
  }
};

fs.writeFileSync(pkgPath, JSON.stringify(pkg, null, 2) + '\n');
console.log('  files:', JSON.stringify(pkg.files));
console.log('  exports:', JSON.stringify(pkg.exports));
" "$dir/package.json" "$module_kind"
}

# ── 1. Web target (browser / Vite) with Rayon parallelisation ────────────────
if [[ "$NODE_ONLY" == false ]]; then
  echo "→ Building WASM — web target (features: $WEB_FEATURES, profile: $PROFILE)..."
  WEB_RUSTFLAGS="${RUSTFLAGS:-}"
  if [[ " $WEB_FEATURES " == *" parallel "* ]]; then
    WEB_RUSTFLAGS="${WEB_RUSTFLAGS:+$WEB_RUSTFLAGS }-C target-feature=+atomics,+bulk-memory,+mutable-globals"
  fi
  CARGO_TARGET_WASM32_UNKNOWN_UNKNOWN_RUSTFLAGS="$WEB_RUSTFLAGS" RUSTFLAGS="$WEB_RUSTFLAGS" wasm-pack build --target web "$PROFILE_FLAG" --features "$WEB_FEATURES"

  # ── 2. Copy the README and the alpha/beta ESM subpath shims ────────────
  echo "→ Copying README.md → pkg/README.md..."
  cp ../../README.md pkg/README.md

  echo "→ Copying alpha/beta ESM subpath shims into pkg/..."
  mkdir -p pkg/alpha pkg/beta
  cp js/alpha/index.js pkg/alpha/index.js
  cp js/alpha/index.d.ts pkg/alpha/index.d.ts
  cp js/beta/index.js pkg/beta/index.js
  cp js/beta/index.d.ts pkg/beta/index.d.ts
fi

# ── 4. Node.js target (server-side / CommonJS) ───────────────────────────────
if [[ "$WEB_ONLY" == false ]]; then
  echo "→ Building WASM — nodejs target (features: ${NODE_FEATURES:-<none>}, profile: $PROFILE)..."
  if [[ -n "${NODE_FEATURES// }" ]]; then
    wasm-pack build --target nodejs "$PROFILE_FLAG" --out-dir pkg-node --no-default-features --features "$NODE_FEATURES"
  else
    wasm-pack build --target nodejs "$PROFILE_FLAG" --out-dir pkg-node --no-default-features
  fi

  echo "→ Copying README.md → pkg-node/README.md..."
  cp ../../README.md pkg-node/README.md
  # The nodejs target of wasm-pack emits a self-contained CommonJS glue file
  # and therefore no `snippets/` directory of its own.  The npm `files`
  # allowlist of pkg-node/package.json lists "snippets" (it is the entry
  # wasm-bindgen-rayon needs at runtime in the ESM package), so the directory
  # is copied over from the web build to keep the two packages in sync and
  # the allowlist valid.
  if [[ -d pkg/snippets ]]; then
    echo "→ Copying snippets/ → pkg-node/snippets/..."
    rm -rf pkg-node/snippets
    cp -r pkg/snippets pkg-node/snippets
  fi

  # The ESM shims under js/ use `export { ... } from '../sci_form_wasm.js'`;
  # pkg-node/ is CommonJS, so the subpath entry points are plain
  # `module.exports = require(...)` re-exports.
  echo "→ Creating alpha/beta CommonJS subpath shims in pkg-node/..."
  mkdir -p pkg-node/alpha pkg-node/beta
  cat > pkg-node/alpha/index.js <<'EOF'
'use strict';
module.exports = require('../sci_form_wasm.js');
EOF
  cat > pkg-node/beta/index.js <<'EOF'
'use strict';
module.exports = require('../sci_form_wasm.js');
EOF
  cp js/alpha/index.d.ts pkg-node/alpha/index.d.ts
  cp js/beta/index.d.ts pkg-node/beta/index.d.ts

  # ── 5. Assemble the dual package: pkg/cjs/ (CommonJS) ───────────────────
  # The published package is `pkg/`.  Its `exports` map sends `import` to the
  # ESM half at the root and `require` to the CommonJS half in `pkg/cjs/`, so
  # `require('sci-form-wasm')` works in Node.js — the path the TypeScript
  # guide documents — and `import ... from 'sci-form-wasm'` keeps working in
  # browsers and bundlers.
  #
  # The nested `pkg/cjs/package.json` marks that subtree as CommonJS; the
  # root package is `"type": "module"`, so without it every `.js` file below
  # `cjs/` would be parsed as ESM and `require()` would fail.
  if [[ "$NODE_ONLY" == false ]]; then
    echo "→ Assembling the dual package: pkg/cjs/ (CommonJS) ..."
    rm -rf pkg/cjs
    mkdir -p pkg/cjs
    cp pkg-node/sci_form_wasm.js pkg/cjs/sci_form_wasm.js
    cp pkg-node/sci_form_wasm_bg.wasm pkg/cjs/sci_form_wasm_bg.wasm
    cp pkg-node/sci_form_wasm.d.ts pkg/cjs/sci_form_wasm.d.ts
    cp pkg-node/sci_form_wasm_bg.wasm.d.ts pkg/cjs/sci_form_wasm_bg.wasm.d.ts
    cp -r pkg-node/alpha pkg/cjs/alpha
    cp -r pkg-node/beta pkg/cjs/beta
    printf '{\n  "type": "commonjs"\n}\n' > pkg/cjs/package.json
  fi
fi

# ── 6. Patch the package manifests ───────────────────────────────────────────
# The patch runs last: the ESM manifest only grows the `require` conditions
# when the CommonJS half is actually present.
if [[ "$NODE_ONLY" == false ]]; then
  echo "→ Patching pkg/package.json..."
  patch_package_json pkg esm
fi
if [[ "$WEB_ONLY" == false ]]; then
  echo "→ Patching pkg-node/package.json..."
  patch_package_json pkg-node cjs
fi

# ── 7. Self-check: the CommonJS half must load ──────────────────────────────
if [[ "$NODE_ONLY" == false && "$WEB_ONLY" == false ]]; then
  echo "→ Checking that the CommonJS half of the dual package loads..."
  node -e "
const sci = require('./pkg/cjs/sci_form_wasm.js');
const alpha = require('./pkg/cjs/alpha/index.js');
const beta = require('./pkg/cjs/beta/index.js');
const r = JSON.parse(sci.embed('CCO', 42));
if (r.num_atoms !== 9) throw new Error('expected 9 atoms, got ' + r.num_atoms);
if (typeof alpha.alpha_modules_info !== 'function') throw new Error('missing alpha subpath export');
if (typeof beta.beta_modules_info !== 'function') throw new Error('missing beta subpath export');
console.log('  ✓ pkg/cjs loads and embeds');
"
fi

echo "✓ sci-form-wasm build complete."
if [[ "$NODE_ONLY" == false ]]; then
  echo "  Web pkg : $SCRIPT_DIR/pkg/ (ESM, + pkg/cjs/ CommonJS)"
fi
if [[ "$WEB_ONLY" == false ]]; then
  echo "  Node pkg: $SCRIPT_DIR/pkg-node/"
fi
