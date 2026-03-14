#!/usr/bin/env bash
# Usage:
#   ./scripts/bump_version.sh          # auto-increment patch (0.1.4 → 0.1.5)
#   ./scripts/bump_version.sh 0.2.0    # set explicit version

set -euo pipefail

# ── Resolve new version ───────────────────────────────────────────────────────
CURRENT=$(grep '^version' Cargo.toml | head -1 | sed 's/version = "\(.*\)"/\1/')

if [[ $# -ge 1 ]]; then
  NEW_VERSION="$1"
  # Basic semver format check
  if ! [[ "$NEW_VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "Error: version must be in X.Y.Z format, got: $NEW_VERSION" >&2
    exit 1
  fi
else
  # Auto-increment patch
  MAJOR=$(echo "$CURRENT" | cut -d. -f1)
  MINOR=$(echo "$CURRENT" | cut -d. -f2)
  PATCH=$(echo "$CURRENT" | cut -d. -f3)
  NEW_VERSION="$MAJOR.$MINOR.$((PATCH + 1))"
fi

echo "Bumping: $CURRENT → $NEW_VERSION"

# ── Update files ──────────────────────────────────────────────────────────────

# Cargo.toml files — replace exactly the `version = "X.Y.Z"` line in [package]
for toml in \
    Cargo.toml \
    crates/cli/Cargo.toml \
    crates/wasm/Cargo.toml \
    crates/python/Cargo.toml; do
  if [[ -f "$toml" ]]; then
    # Only update the first occurrence (the [package] version, not dependency versions)
    sed -i "0,/^version = \"[0-9]*\.[0-9]*\.[0-9]*\"/s//version = \"$NEW_VERSION\"/" "$toml"
    echo "  updated $toml"
  fi
done

# pyproject.toml — update version = "X.Y.Z" under [project]
PYPROJECT="crates/python/pyproject.toml"
if [[ -f "$PYPROJECT" ]]; then
  sed -i "s/^version = \"[0-9]*\.[0-9]*\.[0-9]*\"/version = \"$NEW_VERSION\"/" "$PYPROJECT"
  echo "  updated $PYPROJECT"
fi

# package.json files — update "version": "X.Y.Z"
for pkg in \
    crates/wasm/pkg/package.json \
    pkg/package.json \
    pkg-node/package.json; do
  if [[ -f "$pkg" ]]; then
    sed -i "s/\"version\": \"[0-9]*\.[0-9]*\.[0-9]*\"/\"version\": \"$NEW_VERSION\"/" "$pkg"
    echo "  updated $pkg"
  fi
done

# ── Git commit + tag + push ───────────────────────────────────────────────────
git add \
  Cargo.toml \
  crates/cli/Cargo.toml \
  crates/wasm/Cargo.toml \
  crates/python/Cargo.toml \
  crates/python/pyproject.toml \
  crates/wasm/pkg/package.json \
  pkg/package.json \
  pkg-node/package.json \
  Cargo.lock 2>/dev/null || true

git commit -m "chore: bump version to $NEW_VERSION"

TAG="v$NEW_VERSION"
git tag "$TAG"
echo "Tagged: $TAG"

git push origin main
git push origin "$TAG"
echo "Pushed main + $TAG"
echo "Done! Version is now $NEW_VERSION"
