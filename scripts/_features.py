#!/usr/bin/env python3
"""Print a feature list from the env block of .github/workflows/ci.yml.

    scripts/_features.py RUST_EXPERIMENTAL_FEATURES [--comma]
"""
import re
import sys
from pathlib import Path

name = sys.argv[1]
text = (Path(__file__).resolve().parent.parent / ".github/workflows/ci.yml").read_text()
m = re.search(rf"{name}: >-\n((?:    .*\n)+)", text)
if not m:
    sys.exit(f"{name} not found in ci.yml")
feats = m.group(1).split()
print((",".join if "--comma" in sys.argv else " ".join)(feats))
