#!/usr/bin/env python3
"""Fail when the CI and release workflows disagree about the feature sets.

The published artefacts must expose exactly the capabilities that CI tests
and that the documentation promises.  The two workflows used to carry their
own copy of the feature lists, and the release copy silently fell behind:
the published wheel and npm package were missing the `alpha-edl`,
`alpha-periodic-linear`, `alpha-kinetics`, `alpha-render-bridge` and
`alpha-reaction-dynamics` modules while CI built and tested all of them.

One list, one place: this script compares the lists and fails the build when
they drift apart.

Usage:
    python3 scripts/check_feature_parity.py [--ci PATH] [--release PATH]

Exit codes:
    0 — the lists agree.
    1 — they differ (the differences are printed).
    2 — a list is missing from one of the workflows.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# The variables that must be identical in both workflows.
SHARED_VARS = (
    "PYTHON_BINDING_FEATURES",
    "WASM_WEB_FEATURES",
    "WASM_NODE_FEATURES",
)


def extract_env_vars(text: str) -> dict[str, str]:
    """Return the workflow-level ``env:`` variables and their values.

    A tiny parser instead of a YAML import: the workflows use folded
    scalars (``>-``) whose value is the concatenation of the indented
    continuation lines, and the script must run on a bare python3.
    """
    env: dict[str, str] = {}
    lines = text.splitlines()
    in_env = False
    block_indent = -1
    current: str | None = None
    current_indent = -1

    for line in lines:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indent = len(line) - len(line.lstrip())

        # A folded scalar continues while the lines stay indented deeper
        # than the key that introduced it.
        if current is not None:
            if indent > current_indent:
                env[current] = f"{env[current]} {line.strip()}".strip()
                continue
            current = None

        if in_env and indent <= block_indent and re.match(r"^\S", line):
            in_env = False

        if re.match(r"^env:\s*$", line):
            in_env = True
            block_indent = indent
            continue

        if in_env:
            match = re.match(r"^\s+([A-Z0-9_]+):\s*(.*)$", line)
            if match:
                name, value = match.group(1), match.group(2).strip()
                if value in (">-", ">", "|", "|-"):
                    env[name] = ""
                    current = name
                    current_indent = indent
                else:
                    env[name] = value.strip("'\"")
    return {k: " ".join(v.split()) for k, v in env.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ci", default=".github/workflows/ci.yml")
    parser.add_argument("--release", default=".github/workflows/release.yml")
    args = parser.parse_args()

    ci = extract_env_vars(Path(args.ci).read_text(encoding="utf-8"))
    release = extract_env_vars(Path(args.release).read_text(encoding="utf-8"))

    failures: list[str] = []
    for var in SHARED_VARS:
        if var not in ci:
            failures.append(f"{args.ci} does not define {var}")
            continue
        if var not in release:
            failures.append(f"{args.release} does not define {var}")
            continue
        if ci[var] != release[var]:
            failures.append(
                f"{var} differs\n"
                f"    ci:     {ci[var]}\n"
                f"    release: {release[var]}"
            )

    if failures:
        print("The CI and release feature lists must be identical:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        return 1

    for var in SHARED_VARS:
        print(f"OK  {var}: {ci[var]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
