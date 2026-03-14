#!/usr/bin/env bash
# Tune WASM binary size with wasm-opt.
# Usage: ./scripts/wasm_opt_tune.sh [pkg-dir]
#
# Prerequisite: install wasm-opt via `npm i -g wasm-opt` or `brew install binaryen`
set -e

PKG_DIR="${1:-pkg}"
WASM_FILE="${PKG_DIR}/sci_form_bg.wasm"

if [ ! -f "$WASM_FILE" ]; then
  echo "Error: $WASM_FILE not found. Run 'wasm-pack build --target web' first."
  exit 1
fi

if ! command -v wasm-opt &>/dev/null; then
  echo "Error: wasm-opt not found. Install via: npm i -g wasm-opt  or  brew install binaryen"
  exit 1
fi

ORIG_SIZE=$(wc -c < "$WASM_FILE")
echo "Original: ${ORIG_SIZE} bytes ($(( ORIG_SIZE / 1024 )) KiB)"

# Optimization levels: -O1, -O2, -O3, -O4, -Os (size), -Oz (aggressive size)
for OPT in "-Os" "-Oz" "-O3"; do
  OUT="${WASM_FILE%.wasm}.opt${OPT#-}.wasm"
  wasm-opt "$WASM_FILE" "$OPT" -o "$OUT" 2>/dev/null || continue
  OPT_SIZE=$(wc -c < "$OUT")
  SAVINGS=$(( (ORIG_SIZE - OPT_SIZE) * 100 / ORIG_SIZE ))
  echo "${OPT}: ${OPT_SIZE} bytes ($(( OPT_SIZE / 1024 )) KiB) — ${SAVINGS}% smaller"
done

# Apply best optimization (-Oz) in-place
BEST="${WASM_FILE%.wasm}.optOz.wasm"
if [ -f "$BEST" ]; then
  cp "$BEST" "$WASM_FILE"
  FINAL_SIZE=$(wc -c < "$WASM_FILE")
  echo ""
  echo "Applied -Oz: ${FINAL_SIZE} bytes ($(( FINAL_SIZE / 1024 )) KiB)"
  # Clean up
  rm -f "${WASM_FILE%.wasm}".opt*.wasm
fi
