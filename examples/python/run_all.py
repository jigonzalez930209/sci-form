#!/usr/bin/env python3
import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import traceback

modules = [
    "01_embed",
    "02_charges_sasa_dipole",
    "03_eht_forcefields",
    "04_quantum_methods",
    "05_reactivity_ml_stereo",
    "06_spectroscopy",
    "07_nmr",
    "08_solvation_materials",
    "09_clustering_transport",
    "10_alpha",
    "11_beta",
]

ok = 0
for name in modules:
    print("\n" + "#" * 72)
    print(f"# RUNNING {name}")
    print("#" * 72)
    try:
        m = importlib.import_module(name)
        m.run()
        ok += 1
    except Exception as e:
        print(f"FAIL in {name}: {type(e).__name__}: {e}")
        traceback.print_exc()
print("\n" + "=" * 72)
print(f"Passed example modules: {ok}/{len(modules)}")
