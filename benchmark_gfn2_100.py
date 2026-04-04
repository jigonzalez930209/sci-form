#!/usr/bin/env python3
"""100-SMILES GFN2-xTB benchmark: sci-form vs tblite."""
import subprocess, json, sys, os, time
import numpy as np

os.environ.setdefault("LD_PRELOAD", "/usr/lib/x86_64-linux-gnu/libgfortran.so.5")
from tblite.interface import Calculator

EV_PER_HARTREE = 27.21138505
KCAL_PER_HARTREE = 627.509474
SCIFORM_BIN = "./target/release/sci-form"


def run_tblite_gfn2(elements, positions):
    numbers = np.array(elements)
    pos_bohr = np.array(positions) / 0.529177
    calc = Calculator("GFN2-xTB", numbers, pos_bohr)
    old_fd = os.dup(1)
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1)
    try:
        res = calc.singlepoint()
    finally:
        os.dup2(old_fd, 1)
        os.close(old_fd)
        os.close(devnull)
    return res.get("energy")


def run_sciform_gfn2(smiles, positions):
    flat = []
    for p in positions:
        flat.extend(p)
    coords_json = json.dumps(flat)
    cmd = [SCIFORM_BIN, "neb-energy", smiles, coords_json, "--method", "gfn2"]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd="/home/lestad/github/sci-form", timeout=60)
    try:
        data = json.loads(result.stdout)
        if isinstance(data, list):
            data = data[0]
        e_kcal = data.get("energy_kcal_mol")
        if e_kcal is not None:
            return e_kcal / KCAL_PER_HARTREE
        return None
    except:
        return None


def embed_sciform(smiles):
    cmd = [SCIFORM_BIN, "embed", smiles]
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd="/home/lestad/github/sci-form")
    try:
        data = json.loads(result.stdout)
        elements = data["elements"]
        coords = data["coords"]
        positions = []
        for i in range(0, len(coords), 3):
            positions.append([coords[i], coords[i+1], coords[i+2]])
        return elements, positions
    except:
        return None, None


# Read first 100 SMILES from GDB20
smi_file = "/home/lestad/github/sci-form/GDB20.1000.smi"
with open(smi_file) as f:
    all_smiles = [line.strip().split()[0] for line in f if line.strip()][:100]

print(f"GFN2-xTB 100-SMILES Benchmark: sci-form vs tblite", flush=True)
print(f"{'='*90}", flush=True)

results = []
t0 = time.time()

for idx, smi in enumerate(all_smiles):
    elements, positions = embed_sciform(smi)
    if elements is None:
        print(f"[{idx+1:3d}/100] {smi:25s}  EMBED FAILED", flush=True)
        results.append({"smi": smi, "status": "embed_fail"})
        continue

    e_tblite = run_tblite_gfn2(elements, positions)
    e_sciform = run_sciform_gfn2(smi, positions)

    if e_sciform is None:
        print(f"[{idx+1:3d}/100] {smi:25s}  SCIFORM FAILED", flush=True)
        results.append({"smi": smi, "status": "sciform_fail"})
        continue

    diff = abs(e_sciform - e_tblite)
    pct = abs(diff / e_tblite) * 100
    n = len(elements)
    elem_set = sorted(set(elements))
    status = "OK" if pct < 0.1 else ("CLOSE" if pct < 0.5 else "FAIL")
    print(f"[{idx+1:3d}/100] {smi:25s}  ({n:2d} at)  "
          f"tblite={e_tblite:12.6f}  sci={e_sciform:12.6f}  "
          f"diff={diff:.6f}  ({pct:.4f}%)  [{status}]", flush=True)
    results.append({
        "smi": smi, "status": status, "pct": pct, "diff_ha": diff,
        "e_tblite": e_tblite, "e_sciform": e_sciform,
        "n_atoms": n, "elements": elem_set,
    })

elapsed = time.time() - t0

# Summary
print(f"\n{'='*90}")
ok   = [r for r in results if r.get("status") == "OK"]
close = [r for r in results if r.get("status") == "CLOSE"]
fail  = [r for r in results if r.get("status") == "FAIL"]
errs  = [r for r in results if r.get("status") in ("embed_fail", "sciform_fail")]
pcts  = [r["pct"] for r in results if "pct" in r]

print(f"OK (<0.1%): {len(ok)}  CLOSE (0.1-0.5%): {len(close)}  FAIL (>0.5%): {len(fail)}  ERR: {len(errs)}")
if pcts:
    print(f"Mean: {np.mean(pcts):.4f}%  Median: {np.median(pcts):.4f}%  "
          f"Max: {np.max(pcts):.4f}%  Min: {np.min(pcts):.4f}%")
    under01 = sum(1 for p in pcts if p < 0.1)
    under05 = sum(1 for p in pcts if p < 0.5)
    print(f"Under 0.1%: {under01}/{len(pcts)}  Under 0.5%: {under05}/{len(pcts)}")

# Worst 10
worst = sorted([r for r in results if "pct" in r], key=lambda r: -r["pct"])[:10]
if worst:
    print(f"\nTop 10 worst:")
    for r in worst:
        print(f"  {r['smi']:25s}  {r['pct']:.4f}%  elements={r['elements']}")

# Element analysis
from collections import defaultdict
elem_errors = defaultdict(list)
for r in results:
    if "pct" in r:
        for e in r["elements"]:
            elem_errors[e].append(r["pct"])

print(f"\nPer-element average error:")
SYMS = {1:"H",5:"B",6:"C",7:"N",8:"O",9:"F",14:"Si",15:"P",16:"S",17:"Cl",35:"Br"}
for e in sorted(elem_errors):
    vals = elem_errors[e]
    sym = SYMS.get(e, str(e))
    print(f"  {sym:3s} (Z={e:2d}): mean={np.mean(vals):.4f}%  max={np.max(vals):.4f}%  n={len(vals)}")

print(f"\nElapsed: {elapsed:.1f}s")
