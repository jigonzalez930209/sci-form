import subprocess, json, numpy as np, os, sys
from tblite.interface import Calculator

EV_PER_HARTREE = 27.21138505

molecules = {
    "H2": {"elements": [1,1], "coords": [0,0,0, 0,0,0.74]},
    "HF": {"elements": [1,9], "coords": [0,0,0, 0,0,0.917]},
    "HCl": {"elements": [1,17], "coords": [0,0,0, 0,0,1.275]},
    "HBr": {"elements": [1,35], "coords": [0,0,0, 0,0,1.414]},
    "HI": {"elements": [1,53], "coords": [0,0,0, 0,0,1.609]},
    "H2O": {"elements": [8,1,1], "coords": [0,0,0.117, 0,0.757,-0.470, 0,-0.757,-0.470]},
    "CH4": {"elements": [6,1,1,1,1], "coords": [0,0,0, 0.629,0.629,0.629, -0.629,-0.629,0.629, -0.629,0.629,-0.629, 0.629,-0.629,-0.629]},
    "NH3": {"elements": [7,1,1,1], "coords": [0,0,0.117, 0,0.939,-0.381, 0.813,-0.470,-0.381, -0.813,-0.470,-0.381]},
    "CO": {"elements": [6,8], "coords": [0,0,0, 0,0,1.128]},
    "N2": {"elements": [7,7], "coords": [0,0,0, 0,0,1.098]},
    "F2": {"elements": [9,9], "coords": [0,0,0, 0,0,1.412]},
    "Cl2": {"elements": [17,17], "coords": [0,0,0, 0,0,1.988]},
    "Br2": {"elements": [35,35], "coords": [0,0,0, 0,0,2.281]},
    "I2": {"elements": [53,53], "coords": [0,0,0, 0,0,2.666]},
    "CH3F": {"elements": [6,9,1,1,1], "coords": [0,0,0.63, 0,0,-0.75, 1.03,0,1.03, -0.515,0.89,1.03, -0.515,-0.89,1.03]},
    "CH3Cl": {"elements": [6,17,1,1,1], "coords": [0,0,0.67, 0,0,-1.11, 1.03,0,1.07, -0.515,0.89,1.07, -0.515,-0.89,1.07]},
    "CH3Br": {"elements": [6,35,1,1,1], "coords": [0,0,0.67, 0,0,-1.27, 1.03,0,1.07, -0.515,0.89,1.07, -0.515,-0.89,1.07]},
    "CH3I": {"elements": [6,53,1,1,1], "coords": [0,0,0.67, 0,0,-1.49, 1.03,0,1.07, -0.515,0.89,1.07, -0.515,-0.89,1.07]},
    "PH3": {"elements": [15,1,1,1], "coords": [0,0,0.127, 0,1.186,-0.593, 1.027,-0.593,-0.593, -1.027,-0.593,-0.593]},
    "H2S": {"elements": [16,1,1], "coords": [0,0,0.103, 0,0.958,-0.822, 0,-0.958,-0.822]},
    "SiH4": {"elements": [14,1,1,1,1], "coords": [0,0,0, 0.854,0.854,0.854, -0.854,-0.854,0.854, -0.854,0.854,-0.854, 0.854,-0.854,-0.854]},
    "CCl4": {"elements": [6,17,17,17,17], "coords": [0,0,0, 1.02,1.02,1.02, -1.02,-1.02,1.02, -1.02,1.02,-1.02, 1.02,-1.02,-1.02]},
    "CBr4": {"elements": [6,35,35,35,35], "coords": [0,0,0, 1.10,1.10,1.10, -1.10,-1.10,1.10, -1.10,1.10,-1.10, 1.10,-1.10,-1.10]},
}

# Suppress tblite Fortran output by redirecting fd 6 (Fortran stdout)
devnull = os.open(os.devnull, os.O_WRONLY)
old_stdout_fd = os.dup(1)

results = []
for name, mol in molecules.items():
    els = np.array(mol["elements"])
    cds = np.array(mol["coords"]).reshape(-1, 3) * 1.8897259886
    
    # Suppress stdout for tblite
    os.dup2(devnull, 1)
    try:
        calc = Calculator("GFN2-xTB", els, cds)
        res = calc.singlepoint()
        e_ref = float(res.get("energy"))
    except Exception as e:
        os.dup2(old_stdout_fd, 1)
        results.append((name, None, None, str(e)))
        continue
    os.dup2(old_stdout_fd, 1)
    
    el_json = json.dumps(mol["elements"])
    co_json = json.dumps(mol["coords"])
    try:
        out = subprocess.run(["./target/release/sci-form", "gfn2", el_json, co_json],
                           capture_output=True, text=True, timeout=30)
        r = json.loads(out.stdout)
        e_our = r["total_energy"] / EV_PER_HARTREE
    except Exception as e:
        results.append((name, None, e_ref, str(e)))
        continue
    
    results.append((name, e_our, e_ref, None))

os.close(devnull)
os.close(old_stdout_fd)

print(f"{'Molecule':<10} {'Our(Ha)':>12} {'tblite(Ha)':>12} {'Err(Ha)':>10} {'Err%':>8}")
print("-" * 58)

errors = []
for name, e_our, e_ref, err_msg in results:
    if err_msg:
        print(f"{name:<10} ERROR: {err_msg}")
        continue
    err = abs(e_our - e_ref)
    pct = abs(err / e_ref) * 100
    errors.append((name, pct))
    marker = " <<<" if pct > 1.0 else (" *" if pct > 0.3 else "")
    print(f"{name:<10} {e_our:>12.6f} {e_ref:>12.6f} {err:>10.6f} {pct:>7.3f}%{marker}")

print("-" * 58)
avg = np.mean([e[1] for e in errors])
maxe = max(errors, key=lambda x: x[1])
mine = min(errors, key=lambda x: x[1])
print(f"Average: {avg:.4f}%, Max: {maxe[0]} ({maxe[1]:.4f}%), Min: {mine[0]} ({mine[1]:.4f}%)")
gt1 = [e for e in errors if e[1] > 1.0]
print(f"Molecules >1% error: {len(gt1)}/{len(errors)}" + (f" [{', '.join(e[0] for e in gt1)}]" if gt1 else ""))
