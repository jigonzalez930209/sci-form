"""Shared helpers for sci-form Python examples."""
from sci_form import embed


def require_ok(result, label="operation"):
    if not result.is_ok():
        raise RuntimeError(f"{label} failed: {result.error}")
    return result


def ethanol(seed=42):
    return require_ok(embed("CCO", seed=seed), "embed ethanol")


def benzene(seed=42):
    return require_ok(embed("c1ccccc1", seed=seed), "embed benzene")


def square_sections():
    """Print a separator."""
    print("-" * 72)
