#!/usr/bin/env python3
"""Check tblite orbital occupations and estimate entropy.

Environment: SCIFORM_CLI, SCIFORM_VENV, GFN2_LD_PRELOAD (see
`scripts/gfn2_component_compare.py`).
"""
import os, sys, json, subprocess
from pathlib import Path
import numpy as np

# The reference python packages (tblite) come from the active interpreter or
# from the virtualenv pointed at by SCIFORM_VENV — never from a hard-coded
# path that only exists on one machine.
_VENV = os.environ.get('SCIFORM_VENV')
if _VENV:
    sys.path.insert(0, str(Path(_VENV) / 'lib' /
                            f'python{sys.version_info.major}.{sys.version_info.minor}' /
                            'site-packages'))
# The tblite wheel links against libgfortran; only preload it when the
# distribution actually ships it there.
_PRELOAD = os.environ.get('GFN2_LD_PRELOAD', '/usr/lib/x86_64-linux-gnu/libgfortran.so.5')
if Path(_PRELOAD).exists():
    os.environ['LD_PRELOAD'] = _PRELOAD

from tblite.interface import Calculator, Result

AATOAU = 1.8897259886
CLI = os.environ.get('SCIFORM_CLI',
                    str(Path(__file__).resolve().parents[1] / 'target' / 'release' / 'sci-form'))
EV2HA = 1.0 / 27.21138505

kt = 9.500e-4  # Ha

def embed(smiles, seed=42):
    r = subprocess.run([CLI, 'embed', smiles, '-s', str(seed)], capture_output=True, text=True)
    data = json.loads(r.stdout)
    return data['elements'], data['coords']

# Test benzene (larger gap but more orbitals)
for smiles in ["C", "c1ccccc1", "CCCCCC"]:
    elems, coords = embed(smiles)
    numbers = np.array(elems, dtype=np.int32)
    positions = np.array(coords).reshape(-1, 3) * AATOAU
    calc = Calculator("GFN2-xTB", numbers, positions)
    res = calc.singlepoint(Result())
    
    energy = res.get("energy")
    occs = res.get("orbital-occupations")
    orb_e = res.get("orbital-energies")
    
    print(f"\n=== {smiles} ===")
    print(f"  Total energy: {energy:.10f} Ha")
    print(f"  Occupations: {occs}")
    print(f"  Orbital energies (Ha): {orb_e}")
    
    # Compute what entropy would be
    ts = 0.0
    for f in occs:
        if f > 1e-15 and (1.0 - f) > 1e-15:
            ts += f * np.log(f) + (1.0 - f) * np.log(1.0 - f)
    ts *= kt
    print(f"  Estimated ts (per spin): {ts:.12f} Ha")
    print(f"  Estimated 2*ts (total): {2*ts:.12f} Ha")
    
    # Check: are any occupations non-integer?
    non_int = [(i, f) for i, f in enumerate(occs) if abs(f - round(f)) > 1e-10]
    if non_int:
        print(f"  Non-integer occupations: {non_int}")
    else:
        print(f"  All occupations are integer (no smearing effect)")
