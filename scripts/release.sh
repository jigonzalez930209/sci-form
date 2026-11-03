#!/usr/bin/env bash
# Tag and push a release — only after the local gate passed for this commit.
#
#   scripts/release.sh 1.0.0
#
# GitHub then builds the CLI packages, wheels and npm package, publishes them
# and deploys the documentation (.github/workflows/release.yml).
set -euo pipefail
cd "$(dirname "$0")/.."

VERSION="${1:?usage: scripts/release.sh <version>   (e.g. 1.0.0)}"
TAG="v${VERSION#v}"; VERSION="${TAG#v}"

[ -z "$(git status --porcelain)" ] || { echo "The working tree is not clean."; exit 1; }
[ -f .release-check.ok ] || { echo "Run scripts/release-check.sh first (no stamp found)."; exit 1; }
[ "$(cat .release-check.ok)" = "$(git rev-parse HEAD)" ] || { echo "The gate passed for a different commit. Run scripts/release-check.sh again."; exit 1; }
for m in Cargo.toml crates/cli/Cargo.toml crates/python/Cargo.toml crates/wasm/Cargo.toml crates/python/pyproject.toml; do
  grep -q "^version = \"$VERSION\"" "$m" || { echo "$m is not at $VERSION"; exit 1; }
done
git rev-parse "$TAG" >/dev/null 2>&1 && { echo "$TAG already exists."; exit 1; }

echo "About to tag $(git rev-parse --short HEAD) as $TAG and push it to origin."
echo "This triggers the public release: crates.io, PyPI, npm, GitHub Release, docs."
read -r -p "Type the version ($VERSION) to continue: " ans
[ "$ans" = "$VERSION" ] || { echo "Aborted."; exit 1; }

git tag -a "$TAG" -m "sci-form $VERSION"
git push origin "$(git rev-parse --abbrev-ref HEAD)" "$TAG"
echo "Pushed $TAG. Watch: https://github.com/jigonzalez930209/sci-form/actions"
