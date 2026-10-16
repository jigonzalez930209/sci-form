"""Tier 0: conformer generation and SMILES parsing."""
from sci_form import embed, embed_batch, parse
from common import require_ok


def run():
    print("=== 01 embed ===")
    result = require_ok(embed("CCO", seed=42), "embed")
    print(f"ConformerResult: {result}")
    print("atoms:", result.num_atoms, "elements:", result.elements)
    print("coords len:", len(result.coords))
    print("bonds:", result.bonds)
    print("time_ms:", result.time_ms)
    print("positions[0]:", result.get_positions()[0])

    batch = embed_batch(["O", "CCO", "c1ccccc1"], seed=7, num_threads=0)
    print("batch successes:", sum(1 for r in batch if r.is_ok()), "/", len(batch))
    for r in batch:
        if r.is_ok():
            print("  ", r.smiles, r.num_atoms, "atoms", f"{r.time_ms:.1f}ms")

    parsed = parse("CC(=O)O")
    print("parsed acetic acid:", parsed)
    print("num_atoms", parsed["num_atoms"], "num_bonds", parsed["num_bonds"])
    print("atoms[0]", parsed["atoms"][0])


if __name__ == "__main__":
    run()
