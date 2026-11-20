#!/usr/bin/env bash
# Run cargo inside a memory-capped systemd scope so a heavy build or test run
# cannot push the machine into an OOM kill of the whole desktop session.
#
#   scripts/cargo-guard.sh test --profile numeric --features alpha-v2 --test v2_validation
#
# The scope throttles at MemoryHigh and, past MemoryMax, the kernel kills only
# processes inside the scope (cargo/rustc/test binaries), never the session.
#
# Overrides (environment):
#   GUARD_HIGH_PCT  throttle threshold, % of RAM   (default 60)
#   GUARD_MAX_PCT   hard limit, % of RAM           (default 75)
#   GUARD_SWAP      swap the scope may use         (default 1G)
#   CARGO_BUILD_JOBS parallel rustc jobs           (default 8)
set -euo pipefail

total_kib=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
high_kib=$(( total_kib * ${GUARD_HIGH_PCT:-60} / 100 ))
max_kib=$(( total_kib * ${GUARD_MAX_PCT:-75} / 100 ))
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-8}"

if ! command -v systemd-run >/dev/null 2>&1; then
    echo "cargo-guard: systemd-run not found; running with CARGO_BUILD_JOBS=$CARGO_BUILD_JOBS only" >&2
    exec nice -n 10 cargo "$@"
fi

echo "cargo-guard: MemoryHigh=$((high_kib / 1024))M MemoryMax=$((max_kib / 1024))M" \
     "swap=${GUARD_SWAP:-1G} jobs=$CARGO_BUILD_JOBS" >&2
exec systemd-run --user --scope --quiet --collect \
    -p MemoryHigh="${high_kib}K" \
    -p MemoryMax="${max_kib}K" \
    -p MemorySwapMax="${GUARD_SWAP:-1G}" \
    nice -n 10 cargo "$@"
