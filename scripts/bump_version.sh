#!/usr/bin/env bash
# Usage:
#   ./scripts/bump_version.sh          # auto-increment patch (current → next patch)
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

# Cargo.toml files — update [package] version (first occurrence) + inline dep versions
for toml in \
    Cargo.toml \
    crates/cli/Cargo.toml \
    crates/wasm/Cargo.toml \
    crates/python/Cargo.toml; do
  if [[ -f "$toml" ]]; then
    # Update [package] version (first standalone `version = "X.Y.Z"` line)
    sed -i "0,/^version = \"[0-9]*\.[0-9]*\.[0-9]*\"/s//version = \"$NEW_VERSION\"/" "$toml"
    # Update inline dep version references: version = "X.Y.Z" inside { } blocks
    sed -i "s/version = \"$CURRENT\"/version = \"$NEW_VERSION\"/g" "$toml"
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
if [[ -f "package.json" ]]; then
  sed -i "s/\"version\": \"[0-9]*\.[0-9]*\.[0-9]*\"/\"version\": \"$NEW_VERSION\"/" package.json
  echo "  updated package.json"
fi

for pkg in \
    crates/wasm/pkg/package.json \
    pkg/package.json \
    pkg-node/package.json; do
  if [[ -f "$pkg" ]]; then
    sed -i "s/\"version\": \"[0-9]*\.[0-9]*\.[0-9]*\"/\"version\": \"$NEW_VERSION\"/" "$pkg"
    echo "  updated $pkg"
  fi
done

if [[ -f "package-lock.json" ]]; then
  npm install --package-lock-only --ignore-scripts --no-audit --no-fund
  echo "  updated package-lock.json"
fi

# ── Git commit + tag + push ───────────────────────────────────────────────────
git add \
  package.json \
  package-lock.json \
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
