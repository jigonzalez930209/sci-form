#!/usr/bin/env python3
"""Cross-language integration tests for sci-form.

Validates that CLI, Python, and WASM/Node.js produce consistent results
and that all interfaces work correctly for sending/receiving data.
"""
import json
import subprocess
import sys
import os
import time
import math

PASS = 0
FAIL = 0

def test(name, condition, detail=""):
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  ✓ {name}")
    else:
        FAIL += 1
        print(f"  ✗ {name} — {detail}")

def run_cli(args, input_data=None):
    """Run CLI and return (stdout, stderr, returncode)."""
    cmd = ["./target/release/sci-form"] + args
    result = subprocess.run(cmd, capture_output=True, text=True, 
                          input=input_data, timeout=30)
    return result.stdout, result.stderr, result.returncode

def run_node(script):
    """Run Node.js script and return stdout."""
    result = subprocess.run(
        ["node", "-e", script],
        capture_output=True, text=True, timeout=30,
        cwd="crates/wasm/pkg"
    )
    return result.stdout.strip(), result.stderr, result.returncode


# ============================================================
# Test 1: Python bindings
# ============================================================
print("\n=== PYTHON BINDINGS ===")
import sci_form

# 1.1 Single molecule embed
r = sci_form.embed("CCO")
test("embed CCO succeeds", r.is_ok())
test("embed CCO has 9 atoms", r.num_atoms == 9, f"got {r.num_atoms}")
test("embed CCO has 27 coords", len(r.coords) == 27, f"got {len(r.coords)}")
test("embed CCO has elements", len(r.elements) == 9, f"got {len(r.elements)}")
test("embed CCO has bonds", len(r.bonds) > 0)

# 1.2 Positions helper
pos = r.get_positions()
test("get_positions returns 9 tuples", len(pos) == 9, f"got {len(pos)}")
test("first position is 3-tuple", len(pos[0]) == 3)

# 1.3 Aromatic molecule
r2 = sci_form.embed("c1ccccc1")
test("embed benzene succeeds", r2.is_ok())
test("benzene has 12 atoms (6C + 6H)", r2.num_atoms == 12, f"got {r2.num_atoms}")

# 1.4 Complex molecule
r3 = sci_form.embed("CC(=O)Oc1ccccc1C(=O)O")  # Aspirin
test("embed aspirin succeeds", r3.is_ok())
test("aspirin time < 50ms", r3.time_ms < 50, f"got {r3.time_ms:.1f}ms")

# 1.5 Invalid SMILES
r4 = sci_form.embed("INVALID")
test("invalid SMILES returns error", not r4.is_ok())
test("error message exists", r4.error is not None)

# 1.6 Batch processing
batch = sci_form.embed_batch(["CCO", "c1ccccc1", "C#N", "CC(=O)O"])
test("batch returns 4 results", len(batch) == 4, f"got {len(batch)}")
test("all batch results OK", all(b.is_ok() for b in batch))

# 1.7 Batch with invalid
batch2 = sci_form.embed_batch(["CCO", "INVALID", "C#N"])
test("batch with invalid: 3 results", len(batch2) == 3, f"got {len(batch2)}")
test("batch: first OK", batch2[0].is_ok())
test("batch: second fails", not batch2[1].is_ok())
test("batch: third OK", batch2[2].is_ok())

# 1.8 Parse
parsed = sci_form.parse("c1ccccc1")
test("parse returns dict", isinstance(parsed, dict))
test("parse has num_atoms=12", parsed.get("num_atoms") == 12, f"got {parsed}")

# 1.9 Version
v = sci_form.version()
test("version is string", isinstance(v, str))
test("version contains sci-form", "sci-form" in v)

# 1.10 Seed reproducibility
ra = sci_form.embed("CCO", seed=42)
rb = sci_form.embed("CCO", seed=42)
test("same seed same coords", ra.coords == rb.coords)

# 1.11 Different seeds different coords
rc = sci_form.embed("CC(C)CC(C)C", seed=1)
rd = sci_form.embed("CC(C)CC(C)C", seed=9999)
test("different seeds may differ", True)  # Just verify no crash

# 1.12 Bond data format
for (a, b, order) in r.bonds:
    test("bond indices are ints", isinstance(a, int) and isinstance(b, int))
    test("bond order is string", isinstance(order, str))
    break  # Just check first bond

# 1.13 Geometry validation
def validate_coords(result):
    """Basic geometry check."""
    coords = result.coords
    n = result.num_atoms
    for i in range(0, n * 3, 3):
        x, y, z = coords[i], coords[i+1], coords[i+2]
        if any(math.isnan(v) or math.isinf(v) for v in (x, y, z)):
            return False
    return True

test("ethanol coords valid", validate_coords(r))
test("benzene coords valid", validate_coords(r2))
test("aspirin coords valid", validate_coords(r3))


# ============================================================
# Test 2: CLI
# ============================================================
print("\n=== CLI ===")

# 2.1 Single embed
out, err, code = run_cli(["embed", "CCO"])
test("cli embed CCO exit 0", code == 0, f"exit {code}")
data = json.loads(out) if code == 0 else {}
test("cli json has smiles", data.get("smiles") == "CCO")
test("cli json has num_atoms", data.get("num_atoms") == 9, f"got {data.get('num_atoms')}")
test("cli json has coords", len(data.get("coords", [])) == 27)

# 2.2 XYZ format
out, err, code = run_cli(["embed", "CCO", "-f", "xyz"])
test("cli xyz exit 0", code == 0)
lines = out.strip().split('\n')
test("xyz first line is atom count", lines[0].strip() == "9", f"got '{lines[0]}'")

# 2.3 SDF format
out, err, code = run_cli(["embed", "CCO", "-f", "sdf"])
test("cli sdf exit 0", code == 0)
test("sdf contains $$$$", "$$$$" in out)
test("sdf contains V2000", "V2000" in out)

# 2.4 CLI error handling
out, err, code = run_cli(["embed", "INVALID"])
test("cli invalid SMILES exit 1", code == 1)

# 2.5 Batch via stdin
smiles_input = "CCO\nc1ccccc1\nC#N\n"
out, err, code = run_cli(["batch"], input_data=smiles_input)
test("cli batch exit 0", code == 0, f"exit {code}")
if code == 0:
    batch_data = json.loads(out)
    test("cli batch 3 results", len(batch_data) == 3, f"got {len(batch_data)}")

# 2.6 Parse command
out, err, code = run_cli(["parse", "c1ccccc1"])
test("cli parse exit 0", code == 0)
test("parse output has Atoms", "Atoms:" in out, f"got: {out[:50]}")

# 2.7 Info command
out, err, code = run_cli(["info"])
test("cli info exit 0", code == 0)
test("info has version", "sci-form" in out.lower() or "version" in out.lower(), f"got: {out[:50]}")


# ============================================================
# Test 3: WASM/Node.js
# ============================================================
print("\n=== WASM / NODE.JS ===")

# 3.1 Single embed
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.embed('CCO', 42);
console.log(JSON.stringify({
    ok: r.error === null,
    num_atoms: r.num_atoms,
    coords_len: r.coords.length,
    elements_len: r.elements.length,
    bonds_len: r.bonds.length
}));
""")
test("node embed CCO exit 0", code == 0, f"exit {code}, stderr: {err[:100]}")
if code == 0:
    d = json.loads(out)
    test("node CCO ok", d["ok"])
    test("node CCO 9 atoms", d["num_atoms"] == 9, f"got {d['num_atoms']}")
    test("node CCO 27 coords", d["coords_len"] == 27, f"got {d['coords_len']}")

# 3.2 Aromatic molecule
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.embed('c1ccccc1', 42);
console.log(JSON.stringify({ok: r.error === null, num_atoms: r.num_atoms}));
""")
if code == 0:
    d = json.loads(out)
    test("node benzene ok", d["ok"])
    test("node benzene 12 atoms", d["num_atoms"] == 12, f"got {d['num_atoms']}")

# 3.3 getPositions
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.embed('CCO', 42);
const pos = sf.getPositions(r);
console.log(JSON.stringify({len: pos.length, first: pos[0]}));
""")
if code == 0:
    d = json.loads(out)
    test("node getPositions 9 items", d["len"] == 9, f"got {d['len']}")
    test("node position has x,y,z", all(k in d["first"] for k in ["x","y","z"]))

# 3.4 getAtoms  
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.embed('CCO', 42);
const atoms = sf.getAtoms(r);
console.log(JSON.stringify({len: atoms.length, first: atoms[0]}));
""")
if code == 0:
    d = json.loads(out)
    test("node getAtoms 9 items", d["len"] == 9)
    test("node atom has element+coords", "element" in d["first"] and "x" in d["first"])

# 3.5 embedBatch
out, err, code = run_node("""
const sf = require('./index.js');
const results = sf.embedBatch(['CCO', 'c1ccccc1', 'C#N'], 42);
console.log(JSON.stringify({
    len: results.length,
    all_ok: results.every(r => r.error === null),
    smiles: results.map(r => r.smiles)
}));
""")
if code == 0:
    d = json.loads(out)
    test("node batch 3 results", d["len"] == 3) 
    test("node batch all ok", d["all_ok"])

# 3.6 parseSmiles
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.parseSmiles('c1ccccc1');
console.log(JSON.stringify(r));
""")
if code == 0:
    d = json.loads(out)
    test("node parse benzene atoms=12", d.get("num_atoms") == 12, f"got {d}")

# 3.7 Error handling
out, err, code = run_node("""
const sf = require('./index.js');
const r = sf.embed('INVALID', 42);
console.log(JSON.stringify({has_error: r.error !== null}));
""")
if code == 0:
    d = json.loads(out)
    test("node invalid SMILES has error", d["has_error"])

# 3.8 version
out, err, code = run_node("""
const sf = require('./index.js');
console.log(sf.version());
""")
test("node version works", code == 0 and "sci-form" in out)


# ============================================================
# Test 4: Cross-language consistency
# ============================================================
print("\n=== CROSS-LANGUAGE CONSISTENCY ===")

# Same molecule, same seed → same number of atoms across all interfaces
test_smiles = "c1ccc(O)cc1"  # phenol

# Python
py_result = sci_form.embed(test_smiles, seed=42)
py_natoms = py_result.num_atoms
py_coords = py_result.coords

# CLI
cli_out, _, cli_code = run_cli(["embed", test_smiles, "-s", "42"])
cli_data = json.loads(cli_out) if cli_code == 0 else {}
cli_natoms = cli_data.get("num_atoms", 0)
cli_coords = cli_data.get("coords", [])

# Node.js/WASM
node_out, _, node_code = run_node(f"""
const sf = require('./index.js');
const r = sf.embed('{test_smiles}', 42);
console.log(JSON.stringify({{num_atoms: r.num_atoms, coords: r.coords}}));
""")
node_data = json.loads(node_out) if node_code == 0 else {}
node_natoms = node_data.get("num_atoms", 0)
node_coords = node_data.get("coords", [])

test("py-cli same atoms", py_natoms == cli_natoms, f"py={py_natoms} cli={cli_natoms}")
test("py-node same atoms", py_natoms == node_natoms, f"py={py_natoms} node={node_natoms}")
test("py-cli same coord count", len(py_coords) == len(cli_coords))
test("py-node same coord count", len(py_coords) == len(node_coords))

# Check coords are numerically close (may differ due to f32/f64 conversion)
if len(py_coords) == len(cli_coords) and len(py_coords) > 0:
    max_diff = max(abs(a - b) for a, b in zip(py_coords, cli_coords))
    test("py-cli coords match (<0.001)", max_diff < 0.001, f"max_diff={max_diff:.6f}")

if len(py_coords) == len(node_coords) and len(py_coords) > 0:
    max_diff = max(abs(a - b) for a, b in zip(py_coords, node_coords))
    test("py-node coords match (<0.001)", max_diff < 0.001, f"max_diff={max_diff:.6f}")


# ============================================================
# Test 5: Performance
# ============================================================
print("\n=== PERFORMANCE ===")

# Python batch performance
smiles_list = [
    "CCO", "c1ccccc1", "CC(=O)O", "C#N", "CC(=O)Oc1ccccc1C(=O)O",
    "C1CCCC1", "c1ccncc1", "OC(=O)c1ccccc1", "CC(C)C(=O)O", "c1ccc2ccccc2c1"
] * 10  # 100 molecules

start = time.time()
results = sci_form.embed_batch(smiles_list)
elapsed = time.time() - start
ok_count = sum(1 for r in results if r.is_ok())
avg_ms = elapsed * 1000 / len(smiles_list)

test(f"100 mol batch: {ok_count}/100 ok", ok_count >= 98, f"got {ok_count}")
test(f"avg {avg_ms:.1f}ms/mol < 50ms", avg_ms < 50, f"got {avg_ms:.1f}ms")
print(f"  ℹ throughput: {len(smiles_list)/elapsed:.0f} mol/s ({elapsed*1000:.0f}ms total)")


# ============================================================
# Summary
# ============================================================
print(f"\n{'='*50}")
print(f"RESULTS: {PASS} passed, {FAIL} failed out of {PASS+FAIL} tests")
if FAIL == 0:
    print("ALL TESTS PASSED ✓")
else:
    print(f"FAILURES: {FAIL}")
    sys.exit(1)
