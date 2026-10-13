#!/usr/bin/env bash
# Install MOPAC (the F42 PM3 arbiter) from the official openmopac releases.
#
# Usage: scripts/validation/install_mopac.sh [version]
#
# Downloads the prebuilt Linux tarball (no Fortran toolchain needed), unpacks
# it under $HOME/.local/opt/mopac-<version>/ and prints the wrapper path.
# The wrapper sets LD_LIBRARY_PATH so `mopac <file.mop>` just works.
set -euo pipefail

VERSION="${1:-23.2.5}"
DEST="$HOME/.local/opt/mopac-${VERSION}"
URL="https://github.com/openmopac/mopac/releases/download/v${VERSION}/mopac-${VERSION}-linux.tar.gz"

if [ ! -x "$DEST/bin/mopac" ]; then
  mkdir -p "$DEST"
  echo "Downloading $URL"
  curl -sL --fail -o /tmp/mopac.tar.gz "$URL"
  tar xzf /tmp/mopac.tar.gz -C "$DEST" --strip-components=1
  chmod +x "$DEST/bin/"*
fi

mkdir -p "$HOME/.local/bin"
cat > "$HOME/.local/bin/mopac" <<EOF
#!/usr/bin/env bash
export LD_LIBRARY_PATH="$DEST/lib:$DEST/bin:\${LD_LIBRARY_PATH:-}"
exec "$DEST/bin/mopac" "\$@"
EOF
chmod +x "$HOME/.local/bin/mopac"

echo "MOPAC installed: $DEST"
echo "Wrapper: $HOME/.local/bin/mopac"
echo
echo "Note: a .mop file needs a blank line after the title line, otherwise the"
echo "first geometry line is silently dropped (MOPAC format quirk)."
echo "Reference outputs for validation live in tests/fixtures/validation/mopac_pm3.json."
