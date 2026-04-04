#!/usr/bin/env python3
"""Benchmark GFN2-xTB: our Rust implementation vs tblite for 100+ organic molecules."""
import subprocess, json, sys, os, numpy as np

# 100+ diverse organic molecules with SMILES, elements, and coordinates
# We generate 3D coords using sci-form embed, then compare GFN2 energies
molecules_smiles = [
    "C", "CC", "CCC", "CCCC", "CCCCC",          # alkanes
    "C=C", "CC=C", "CC=CC", "C=CC=C",            # alkenes
    "C#C", "CC#C", "CC#CC",                        # alkynes
    "C1CC1", "C1CCC1", "C1CCCC1", "C1CCCCC1",    # cycloalkanes
    "c1ccccc1", "c1ccncc1", "c1ccoc1", "c1ccsc1", # aromatics
    "c1ccc(cc1)O", "c1ccc(cc1)N", "c1ccc(cc1)F",  # substituted benzene
    "c1ccc(cc1)Cl", "c1ccc(cc1)C",                 # more substituted
    "CO", "CCO", "CCCO", "CC(C)O", "OCC=O",       # alcohols/aldehydes
    "C=O", "CC=O", "CCC=O", "CC(=O)C",            # ketones
    "OC=O", "CC(=O)O", "CCC(=O)O",                # carboxylic acids
    "N", "CN", "CCN", "CCCN", "NC(=O)C",          # amines/amides
    "C#N", "CC#N", "CCC#N",                        # nitriles
    "O=C(N)N", "NC(=O)NC",                         # ureas
    "CS", "CCS", "CCSC",                           # thioethers
    "CF", "CCF", "FC(F)F", "FC(F)(F)F",           # fluorinated
    "CCl", "CCCl", "ClC(Cl)Cl",                    # chlorinated
    "CBr", "CCBr",                                  # brominated
    "C(=O)OC", "CC(=O)OC", "CC(=O)OCC",           # esters
    "COC", "COCC", "C1CCOC1",                      # ethers
    "c1ccc2ccccc2c1",                               # naphthalene
    "C1=CC=CC=C1O",                                 # phenol
    "C1=CC=CC=C1N",                                 # aniline
    "CC(=O)Nc1ccc(O)cc1",                           # acetaminophen
    "OC(=O)c1ccccc1",                               # benzoic acid
    "c1ccc(cc1)c1ccccc1",                           # biphenyl
    "O=C1CCCCC1",                                   # cyclohexanone
    "C1CCNCC1",                                     # piperidine
    "C1CCOCC1",                                     # tetrahydropyran
    "c1cn[nH]c1",                                   # pyrazole
    "c1c[nH]cn1",                                   # imidazole
    "c1ccnc(c1)N",                                  # 2-aminopyridine
    "c1ccc(cc1)C=O",                                # benzaldehyde
    "c1ccc(cc1)C(=O)O",                             # benzoic acid
    "CC(C)(C)C",                                    # neopentane
    "C(C)(C)=O",                                    # acetone
    "C(F)(Cl)Br",                                   # halomethane
    "CC(=O)OC(=O)C",                                # acetic anhydride
    "c1cc(ccc1O)O",                                 # catechol
    "NC(=O)c1ccccc1",                               # benzamide
    "OC(=O)CC(=O)O",                                # malonic acid
    "c1ccc(c(c1)O)O",                               # catechol isomer
    "FN=O",                                         # nitrosyl fluoride
    "O=S=O",                                        # SO2
    "CS(C)=O",                                      # DMSO
    "FC=O",                                         # formyl fluoride
    "O=CO",                                         # formic acid
    "NCC(=O)O",                                     # glycine
    "CC(N)C(=O)O",                                  # alanine
    "O=C(O)C(O)=O",                                 # oxalic acid
    "c1ccc(cc1)OC",                                 # anisole
    "c1cc(cc(c1)O)O",                               # resorcinol
    "OC(CO)CO",                                     # glycerol
    "C1CC1C1CC1",                                   # bicyclopropyl
    "C12CC1CC2",                                    # spiropentane
    "C1=CC2=CC=CC=C2C=C1",                          # azulene
    "c1ccc2[nH]ccc2c1",                             # indole
    "CC(=O)NCC(=O)O",                               # aceturic acid
    "C(#N)C#N",                                     # malononitrile
    "N=C=O",                                        # isocyanic acid
    "O=C=O",                                        # CO2
    "N#N",                                          # N2
    "O=O",                                          # O2
    "[H][H]",                                       # H2
    "FF",                                           # F2
    "ClCl",                                         # Cl2
    "BrBr",                                         # Br2
    "O",                                            # water
    "N",                                            # ammonia
    "P",                                            # PH3
    "S",                                            # H2S
    "C(Cl)(Cl)(Cl)Cl",                              # CCl4
    "C(F)(F)(F)Cl",                                 # CFC-13
    "c1cnc2ccccc2n1",                               # quinoxaline
]

env = os.environ.copy()
env.pop('GFN2_DEBUG', None)

results = []
tblite_results = []

from tblite.interface import Calculator

for i, smi in enumerate(molecules_smiles):
    # Step 1: Generate 3D coords with sci-form embed
    try:
        r_embed = subprocess.run(
            ["./target/release/sci-form", "embed", smi, "-s", "42"],
            capture_output=True, text=True, timeout=30, env=env, cwd="/home/lestad/github/sci-form"
        )
        embed_data = json.loads(r_embed.stdout)
        if embed_data.get("error"):
            results.append((smi, None, False, -1, "embed_error"))
            tblite_results.append((smi, None))
            continue
        elements = embed_data["elements"]
        coords = embed_data["coords"]
    except Exception as e:
        results.append((smi, None, False, -1, f"embed_fail: {e}"))
        tblite_results.append((smi, None))
        continue

    # Step 2: Run our GFN2
    elems_str = json.dumps(elements)
    coords_str = json.dumps(coords)
    try:
        r_gfn2 = subprocess.run(
            ["./target/release/sci-form", "gfn2", elems_str, coords_str],
            capture_output=True, text=True, timeout=120, env=env, cwd="/home/lestad/github/sci-form"
        )
        data = json.loads(r_gfn2.stdout)
        our_eV = data["total_energy"]
        conv = data.get("converged", False)
        iters = data.get("scc_iterations", -1)
        results.append((smi, our_eV, conv, iters, "ok"))
    except Exception as e:
        results.append((smi, None, False, -1, f"gfn2_fail: {e}"))
        tblite_results.append((smi, None))
        continue

    # Step 3: Run tblite
    n = len(elements)
    positions = np.array(coords).reshape(n, 3) * 1.8897259886  # Å to bohr
    numbers = np.array(elements, dtype=int)
    try:
        calc = Calculator("GFN2-xTB", numbers, positions)
        res = calc.singlepoint()
        tbl_eV = res.get("energy") * 27.21138505
        tblite_results.append((smi, tbl_eV))
    except Exception as e:
        tblite_results.append((smi, None))

# Print results
print(f"\n{'#':>3} {'SMILES':<30} {'Our (eV)':>14} {'tblite (eV)':>14} {'Diff (eV)':>12} {'Error%':>10} {'Conv':>5} {'Iter':>5}")
print("-" * 100)
errors_pct = []
for idx, ((smi, our, conv, iters, status), (_, tbl)) in enumerate(zip(results, tblite_results)):
    if our is not None and tbl is not None:
        diff = our - tbl
        pct = abs(diff / tbl) * 100
        errors_pct.append(pct)
        flag = " ✗" if pct > 0.01 else ""
        print(f"{idx+1:>3} {smi:<30} {our:>14.6f} {tbl:>14.6f} {diff:>12.6f} {pct:>9.4f}%{flag} {'Y' if conv else 'N':>5} {iters:>5}")
    else:
        print(f"{idx+1:>3} {smi:<30} {'N/A':>14} {'N/A':>14} {'N/A':>12} {'N/A':>10} {status}")

print(f"\n--- Summary ---")
if errors_pct:
    print(f"Total molecules tested: {len(errors_pct)}")
    print(f"Below 0.01%: {sum(1 for e in errors_pct if e < 0.01)}")
    print(f"Below 0.02%: {sum(1 for e in errors_pct if e < 0.02)}")
    print(f"Below 0.05%: {sum(1 for e in errors_pct if e < 0.05)}")
    print(f"Below 0.10%: {sum(1 for e in errors_pct if e < 0.10)}")
    print(f"Average error: {np.mean(errors_pct):.4f}%")
    print(f"Max error: {max(errors_pct):.4f}%")
    print(f"Median error: {np.median(errors_pct):.4f}%")
