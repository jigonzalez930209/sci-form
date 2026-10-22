"""Locate the upstream reference data used by the validation generators.

The upstream sources (dftd4, grimme-lab/multicharge, grimme-lab/mctc-lib,
tblite, xtb, MOPAC, LAMMPS, YAeHMOP, pormake) are **not** vendored in this
repository: they are large third-party trees with their own licences.

Point `SCIFORM_REF_DATA` at the directory that holds them::

    export SCIFORM_REF_DATA=/path/to/refdata

The default is the repo-local, git-ignored cache directory::

    <repo>/.cache/refdata

Individual scripts may also be pointed at a specific checkout with their own
command-line flag (e.g. `--source-dir`, `--multicharge-dir`).
"""

from __future__ import annotations

import os
from pathlib import Path

#: Repository root (`<repo>/scripts/validation/refdata.py` -> `<repo>`).
REPO = Path(__file__).resolve().parents[2]

#: Default location of the upstream reference data (git-ignored).
DEFAULT_REF_DATA = REPO / ".cache" / "refdata"

#: Exit code used by the `--check` mode when the upstream input is not
#: available on this machine.  `0` means "the committed table is up to date",
#: `1` means "the committed table is OUT OF DATE" (a hard failure) and `2`
#: means "the check could not run" (reported, but not a failure).
EXIT_SKIPPED = 2


def ref_data_dir() -> Path:
    """Return the directory that holds the upstream reference data.

    Honours `$SCIFORM_REF_DATA`; falls back to `<repo>/.cache/refdata`.
    """
    return Path(os.environ.get("SCIFORM_REF_DATA", DEFAULT_REF_DATA)).expanduser()


def skip(reason: str) -> int:
    """Report a skipped `--check` run and return the `EXIT_SKIPPED` code."""
    print(f"SKIPPED: {reason}")
    return EXIT_SKIPPED

