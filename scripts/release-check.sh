#!/usr/bin/env bash
# The release gate. Every test of the project runs HERE, on the maintainer's
# machine; GitHub only builds and publishes (see .github/workflows/).
#
#   scripts/release-check.sh              # the whole gate
#   scripts/release-check.sh --heavy      # + the #[ignore]d heavy tests (long)
#   scripts/release-check.sh --skip-wasm  # when the wasm toolchain is not installed
#
# On success it writes .release-check.ok with the commit it verified;
# scripts/release.sh refuses to tag anything else.
set -uo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"

HEAVY=false; SKIP_WASM=false
for a in "$@"; do case "$a" in --heavy) HEAVY=true;; --skip-wasm) SKIP_WASM=true;; esac; done

RUST_FEATS="$(python3 scripts/_features.py RUST_EXPERIMENTAL_FEATURES)"
PY_FEATS="$(python3 scripts/_features.py PYTHON_BINDING_FEATURES --comma)"
WASM_WEB="$(python3 scripts/_features.py WASM_WEB_FEATURES)"
WASM_NODE="$(python3 scripts/_features.py WASM_NODE_FEATURES)"
CLI_FEATS="$(sed -n '/^\[features\]/,/^\[/p' crates/cli/Cargo.toml | grep -E '^(experimental|alpha|beta)-' | cut -d' ' -f1 | paste -sd, -)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

FAILED=()
step() { # step "name" command...
  local name="$1"; shift
  printf '\n\033[1m== %s\033[0m\n' "$name"
  if "$@"; then printf '\033[32m   ok: %s\033[0m\n' "$name"; else printf '\033[31m   FAILED: %s\033[0m\n' "$name"; FAILED+=("$name"); fi
}

# ── 1. Versions ────────────────────────────────────────────────────────────
check_versions() {
  python3 - <<'PY'
import re, json, sys
from pathlib import Path
def toml_ver(p): return re.search(r'^version = "(.*)"$', Path(p).read_text(), re.M).group(1)
v = {p: toml_ver(p) for p in ["Cargo.toml","crates/cli/Cargo.toml","crates/python/Cargo.toml","crates/wasm/Cargo.toml","crates/python/pyproject.toml"]}
v["package.json"] = json.loads(Path("package.json").read_text())["version"]
for k, x in v.items(): print(f"   {x:<10} {k}")
sys.exit(0 if len(set(v.values())) == 1 else 1)
PY
}
step "versions agree in every manifest" check_versions
step "CI / release feature lists agree" python3 scripts/check_feature_parity.py

# ── 2. Rust ────────────────────────────────────────────────────────────────
step "cargo fmt" cargo fmt --all -- --check
step "clippy -D warnings" cargo clippy --workspace --all-targets --features "$RUST_FEATS" -- -D warnings
step "rust tests (workspace, release)" cargo test --workspace --exclude sci-form-python --all-targets --release --no-fail-fast --features "$RUST_FEATS"
step "doctests" cargo test --doc --release --features "$RUST_FEATS"
step "CLI guide with every command compiled in" cargo test --release -p sci-form-cli --features "$CLI_FEATS" --test guide --target-dir "$WORK/clifull"
if $HEAVY; then
  step "heavy tests (#[ignore])" cargo test --release --features "$RUST_FEATS" -- --ignored --skip bench
fi

# ── 3. The CLI as users get it ─────────────────────────────────────────────
cli_package() {
  cargo build --release -p sci-form-cli || return 1
  local host; host="$(rustc -vV | sed -n 's/^host: //p')"
  rm -rf "$WORK/dist" && bash scripts/package-cli.sh "$host" "$WORK/dist" >/dev/null || return 1
  local tgz; tgz="$(ls "$WORK"/dist/*.tar.gz 2>/dev/null | head -1)"
  [ -n "$tgz" ] || { echo "no tar.gz produced"; return 1; }
  mkdir -p "$WORK/unpack" && tar -xzf "$tgz" -C "$WORK/unpack"
  local exe; exe="$(find "$WORK/unpack" -name sci-form -type f | head -1)"
  local guide runex
  "$exe" --version || return 1
  guide="$("$exe")" && grep -q "Try, in this order" <<<"$guide" || { echo "bare guide missing"; return 1; }
  runex="$("$exe" examples tpsa --run)" && grep -q '20.2' <<<"$runex" || { echo "examples --run failed"; return 1; }
  ls "$WORK"/dist
}
step "CLI packages (tar.gz/deb/rpm) install and run" cli_package

# ── 4. Python ──────────────────────────────────────────────────────────────
python_tests() {
  (cd crates/python && maturin build --release --features "$PY_FEATS" --out "$WORK/wheels") >"$WORK/wheel.log" 2>&1 || { tail -20 "$WORK/wheel.log"; return 1; }
  python3 -m venv "$WORK/venv" && "$WORK/venv/bin/pip" install -q "$WORK"/wheels/*.whl pytest numpy || return 1
  (cd /tmp && "$WORK/venv/bin/python" -I -m pytest -q --no-header -p no:cacheprovider \
      "$ROOT/tests/integration/test_python_integration.py" "$ROOT/crates/python/tests/test_bindings.py") || return 1
  (cd /tmp && "$WORK/venv/bin/python" -I "$ROOT/tests/integration/test_smirks_reactions.py") || return 1
}
step "python wheel + pytest + smirks" python_tests

# ── 5. WASM, docs examples ────────────────────────────────────────────────
wasm_tests() {
  (cd crates/wasm && bash build.sh --web-features "$WASM_WEB" --node-features "$WASM_NODE" && bash build.sh --docs-only --node-features "$WASM_NODE") >"$WORK/wasm.log" 2>&1 || { tail -20 "$WORK/wasm.log"; return 1; }
  node tests/integration/test_wasm_integration.js || return 1
  node scripts/test-live-examples.mjs | tail -4 || return 1
  # the published package, installed by name, as a user would
  (cd crates/wasm/pkg && npm pack --pack-destination "$WORK" >/dev/null) && mkdir "$WORK/npmt" && (cd "$WORK/npmt" && npm init -y >/dev/null && npm install --no-audit --no-fund "$WORK"/sci-form-wasm-*.tgz >/dev/null && node -e "const s=require('sci-form-wasm'); const r=JSON.parse(s.embed('CCO',42)); if(r.num_atoms!==9) throw new Error('npm package broken'); console.log('npm package ok')")
}
if $SKIP_WASM; then echo "(skipping wasm)"; else step "wasm build + node tests + live examples + npm install" wasm_tests; fi

cross_language() {
  # needs the wheel venv, the release CLI and the node wasm package
  [ -x "$WORK/venv/bin/python" ] && [ -d crates/wasm/pkg-node ] || { echo "skipped: needs the python and wasm steps"; return 0; }
  "$WORK/venv/bin/python" -I tests/integration/test_cross_language.py | tail -3
}
$SKIP_WASM || step "python / cli / wasm agree (cross-language)" cross_language

# ── 6. Documentation ───────────────────────────────────────────────────────
docs_build() {
  $SKIP_WASM || node scripts/copy-docs-wasm.mjs
  npm run docs:build 2>&1 | tail -3
}
step "documentation site builds" docs_build

# ── 7. Workflows ───────────────────────────────────────────────────────────
if command -v actionlint >/dev/null 2>&1; then step "actionlint" actionlint .github/workflows/*.yml; fi

# ── Result ─────────────────────────────────────────────────────────────────
echo
if [ ${#FAILED[@]} -eq 0 ]; then
  if [ -n "$(git status --porcelain)" ]; then
    echo "All checks passed, but the working tree has uncommitted changes:"
    echo "commit them and run the gate again (the stamp is per commit)."
    exit 1
  fi
  git rev-parse HEAD > .release-check.ok
  printf '\033[32mALL CHECKS PASSED for %s\033[0m\nNext: scripts/release.sh <version>\n' "$(git rev-parse --short HEAD)"
else
  printf '\033[31mFAILED (%d):\033[0m\n' "${#FAILED[@]}"; printf '  - %s\n' "${FAILED[@]}"
  rm -f .release-check.ok
  exit 1
fi
