#!/usr/bin/env python3
"""RDKit oracle for sci-form's geometry / force-field parity gates.

Covers the four distance-geometry / force-field cross-library gates that were
missing (audit `reports/auditoria-total.md`):

- F71 — distance-bounds matrix parity
        `sci_form::distgeom::calculate_bounds_matrix_opts` vs
        `rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix`.
- F72 — ETKDG embedding parity
        `sci_form::conformer::generate_3d_conformer` vs
        `rdkit.Chem.AllChem.EmbedMolecule` (ETKDGv3).
- F73 — ETKDG-3D force-field parity
        `sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions` vs
        RDKit's `construct3DForceField` term set.  RDKit 2026.09.1 does **not**
        export `construct3DForceField` to Python, so the oracle rebuilds the
        term set from the upstream C++ (RDKit `Release_2026_09_1`,
        `Code/DistGeom/DistGeomUtils.cpp`) and evaluates the distance / angle
        terms with RDKit's own `ForceField` engine.
- F74 — BFGS minimisation parity
        `sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit` vs RDKit's
        `ForceField::minimize` (same BFGS code path that RDKit uses inside
        `EmbedMolecule`) and vs `rdForceFieldHelpers.UFFOptimizeMolecule`.

Every subcommand:
  1. computes reference values with RDKit on a frozen corpus,
  2. writes a JSON fixture under tests/fixtures/validation/,
  3. emits a report row into reports/validation/summary.json.

The Rust side (`tests/rdkit_parity.rs`) recomputes the sci-form values from
the same fixtures and asserts the tolerances, so CI does not need RDKit.

Usage:
    python scripts/validation/oracle_rdkit_geometry.py bounds|embed|etkdg3d|bfgs|all
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

from common import (
    REPO,
    Report,
    emit_fixture_report,
    kabsch_rmsd,
    read_smiles_corpus,
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDistGeom, rdForceFieldHelpers
except ImportError as exc:  # pragma: no cover
    print(f"RDKit is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

RDKIT_VERSION = Chem.rdBase.rdkitVersion
REFERENCE = f"RDKit {RDKIT_VERSION}"

# Fixed seeds — every number in the fixtures must be reproducible.
ETKDG_SEED = 20261021
PERTURB_SEED = 20261022

# ---------------------------------------------------------------------------
# Frozen corpus
# ---------------------------------------------------------------------------

#: Molecules used for the geometry gates (bounds matrix, embedding, BFGS).
#: Chosen to cover: alkanes, alkenes, alkynes, carbonyls, carboxylic acids,
#: nitriles, alcohols, ethers, amines, aromatics (benzene / pyridine / pyrrole /
#: thiophene), saturated rings, cyclic ketones, amides, esters, halides,
#: sulfoxides, a branched aromatic and a pyranose sugar.
GEOMETRY_CORPUS = [
    "C",
    "CC",
    "CCC",
    "CCCC",
    "CC(C)C",
    "C=C",
    "C#C",
    "CC=O",
    "CC(=O)O",
    "CC#N",
    "CCO",
    "COC",
    "CN",
    "CCN",
    "c1ccccc1",
    "Cc1ccccc1",
    "c1ccncc1",
    "c1cc[nH]c1",
    "c1ccsc1",
    "C1CCCCC1",
    "O=C1CCCCC1",
    "C1CCNCC1",
    "CC(=O)NC",
    "COC(C)=O",
    "CCCl",
    "CCBr",
    "CS(C)=O",
    "CC(C)(C)c1ccccc1",
    "OCC1OC(O)C(O)C(O)C1O",
]

#: Molecules used for the ETKDG-3D force-field gate (F73).
ETKDG3D_CORPUS = [
    "CC",
    "CCC",
    "CCCC",
    "CC(C)C",
    "C=C",
    "C#C",
    "CC=O",
    "CC#N",
    "CCO",
    "COC",
    "CN",
    "c1ccccc1",
    "Cc1ccccc1",
    "c1ccncc1",
    "C1CCCCC1",
    "CC(=O)NC",
    "COC(C)=O",
    "CCCl",
    "CS(C)=O",
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def with_hs(smiles: str) -> Chem.Mol | None:
    """RDKit molecule with explicit hydrogens (matches sci-form's parser)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    if len(Chem.GetMolFrags(mol)) > 1:
        return None  # sci-form keeps only the largest fragment
    return Chem.AddHs(mol)


def element_list(mol: Chem.Mol) -> list[int]:
    return [a.GetAtomicNum() for a in mol.GetAtoms()]


def pair_list(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def bounds_from_rdkit(mol: Chem.Mol, smooth: bool) -> tuple[list[float], list[float]]:
    """RDKit bounds matrix as (lower[], upper[]) in row-major (i<j) order.

    RDKit stores the *lower* bound in the lower triangle and the *upper* bound
    in the upper triangle; sci-form uses the same convention
    (`bounds[(j, i)]` = lower, `bounds[(i, j)]` = upper for `i < j`).
    """
    bm = rdDistGeom.GetMoleculeBoundsMatrix(
        mol,
        True,  # set15bounds
        False,  # scaleVDW (ignored by RDKit)
        smooth,  # doTriangleSmoothing
        False,  # useMacrocycle14config
        True,  # forceTransAmides
        True,  # set14bounds
        True,  # set13bounds
    )
    n = mol.GetNumAtoms()
    lower, upper = [], []
    for i, j in pair_list(n):
        lower.append(float(bm[j, i]))
        upper.append(float(bm[i, j]))
    return lower, upper


def topo_distances(mol: Chem.Mol) -> list[int]:
    """Topological (graph) distance for every (i<j) pair."""
    dm = Chem.GetDistanceMatrix(mol)
    n = mol.GetNumAtoms()
    return [int(round(dm[i, j])) for i, j in pair_list(n)]


def conformer_coords(mol: Chem.Mol, conf_id: int = 0) -> list[float]:
    conf = mol.GetConformer(conf_id)
    out: list[float] = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        out.extend([p.x, p.y, p.z])
    return out


def perturb(coords: list[float], seed: int, sigma: float = 0.25) -> list[float]:
    """Deterministic random displacement (same on both sides via the fixture)."""
    rng = np.random.default_rng(seed)
    arr = np.asarray(coords, dtype=float).reshape(-1, 3)
    return (arr + rng.normal(scale=sigma, size=arr.shape)).reshape(-1).tolist()


def embed_rdkit(smiles: str, seed: int = ETKDG_SEED):
    """ETKDGv3 embedding with a fixed seed (returns mol or None)."""
    mol = with_hs(smiles)
    if mol is None:
        return None
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = False
    cid = AllChem.EmbedMolecule(mol, params)
    if cid < 0:
        return None
    return mol


# ---------------------------------------------------------------------------
# F71 — distance-bounds matrix parity
# ---------------------------------------------------------------------------

#: Per-topological-distance-class tolerances (Å).  Measured against
#: RDKit 2026.09.1 on the frozen corpus (see `reports/validation/summary.json`
#: for the live numbers produced by `tests/rdkit_parity.rs`).
#:
#: * class 1 (1-2 bonds): sci-form reproduces RDKit's UFF bond-length formula
#:   exactly (0.0 Å); the 0.03 Å slack covers the C-Br/Halogen bond length
#:   difference (0.026 Å measured on CCBr).
#: * class 2 (1-3): RDKit derives the distance window from an *angle* window
#:   (`compute13Dist(lb12, lb12, angle-tol)`) while sci-form uses the ideal
#:   angle with a fixed ±0.04 Å distance tolerance, so the windows differ by up
#:   to 0.34 Å for aromatic heteroatom rings (thiophene).
#: * class 3 (1-4): same reason (torsion cis/trans distances) — up to 0.37 Å.
#: * class 4 (1-5): **known structural divergence** — sci-form enumerates 1-4
#:   paths *through ring bonds* (so it can derive a 1-5 bound for e.g. two
#:   hydrogens of a 5-membered aromatic ring) whereas RDKit only enumerates
#:   ring paths along the ring itself and therefore leaves the upper bound at
#:   1000 Å.  The tolerance therefore only constrains the pairs where *both*
#:   sides set a finite upper bound; the number of divergent pairs is pinned
#:   separately (see `n_divergent_1_5`).
#: * class 5/6 (1-6 and beyond): VDW-only lower bounds, up to 0.54 / 0.95 Å
#:   measured (glucose, a pyranose with five hydroxyls).
BOUNDS_TOLERANCES = {
    "1": 0.03,
    "2": 0.35,
    "3": 0.40,
    "4": 1.00,
    "5": 0.60,
    "6": 1.00,
}


def run_bounds() -> Report:
    cases = []
    for smi in GEOMETRY_CORPUS:
        mol = with_hs(smi)
        if mol is None:
            continue
        lower, upper = bounds_from_rdkit(mol, smooth=False)
        cases.append(
            {
                "smiles": smi,
                "n_atoms": mol.GetNumAtoms(),
                "elements": element_list(mol),
                "lower": lower,
                "upper": upper,
                "topo": topo_distances(mol),
            }
        )
    payload = {
        "reference": REFERENCE,
        "rdkit_api": (
            "rdkit.Chem.rdDistGeom.GetMoleculeBoundsMatrix(mol, set15bounds=True,"
            " scaleVDW=False, doTriangleSmoothing=False, useMacrocycle14config=False,"
            " forceTransAmides=True, set14bounds=True, set13bounds=True)"
        ),
        "sci_form_api": "sci_form::distgeom::calculate_bounds_matrix_opts(mol, use_set15=true)",
        "tolerances": BOUNDS_TOLERANCES,
        "tolerance_notes": {
            "1": "1-2 bond lengths: UFF formula, exact match except halogen bonds (0.026 A)",
            "2": "1-3: RDKit uses an angle window, sci-form a distance window (max 0.34 A)",
            "3": "1-4: cis/trans distance windows differ slightly (max 0.37 A)",
            "4": "1-5: sci-form derives extra 1-5 upper bounds for paths through ring bonds; "
            "only pairs bounded on both sides are compared (see n_divergent_1_5)",
            "5": "1-6: VDW lower bound only, max 0.54 A (glucose)",
            "6": "1-7 and beyond: VDW lower bound only, max 0.95 A (glucose)",
        },
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F71",
        observable="bounds_matrix_max_abs",
        reference_tool=f"{REFERENCE} GetMoleculeBoundsMatrix",
        units="angstrom",
        tolerance=BOUNDS_TOLERANCES["2"],
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "per-pair lower/upper distance bounds; tolerances are per topological "
            "distance class and are read from the fixture"
        ),
    )
    emit_fixture_report("bounds_matrix_rdkit.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F72 — ETKDG embedding parity
# ---------------------------------------------------------------------------

#: Maximum deviation allowed between the heavy-atom bond lengths of the two
#: independent embeddings, expressed as a fraction of the distance-bounds
#: window that RDKit allows for that bond.  RDKit's ETKDGv3 and sci-form's
#: ETKDG port use the same UFF bond lengths but different random streams, so
#: the individual bond lengths differ slightly while both stay inside the
#: same allowed window.
EMBED_BOND_LENGTH_TOL = 1.0
#: Upper bound on the Kabsch RMSD between the two independent embeddings (Å).
#: Two independent embeddings of the same molecule are *not* expected to be
#: identical (different random streams); the RMSD is recorded as a measured
#: statistic and only bounded from above.
EMBED_RMSD_TOL = 3.0


def run_embed() -> Report:
    cases = []
    for smi in GEOMETRY_CORPUS:
        mol = embed_rdkit(smi)
        if mol is None:
            continue
        coords = conformer_coords(mol)
        elements = element_list(mol)
        heavy = [i for i, z in enumerate(elements) if z != 1]
        heavy_coords = [coords[3 * i + d] for i in heavy for d in range(3)]
        lower, upper = bounds_from_rdkit(mol, smooth=False)
        bounds_lookup = {
            pair: (lo, up)
            for pair, lo, up in zip(pair_list(mol.GetNumAtoms()), lower, upper)
        }
        bonds = []
        for b in mol.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if elements[i] == 1 or elements[j] == 1:
                continue
            order = {
                Chem.BondType.SINGLE: 1,
                Chem.BondType.DOUBLE: 2,
                Chem.BondType.TRIPLE: 3,
                Chem.BondType.AROMATIC: 4,
            }.get(b.GetBondType(), 1)
            bonds.append([min(i, j), max(i, j), order])
        bond_lengths = []
        bond_bounds = []
        for i, j, _ in bonds:
            d = np.linalg.norm(
                np.asarray(coords[3 * i : 3 * i + 3]) - np.asarray(coords[3 * j : 3 * j + 3])
            )
            bond_lengths.append(float(d))
            lo, up = bounds_lookup[(min(i, j), max(i, j))]
            bond_bounds.append([lo, up])
        cases.append(
            {
                "smiles": smi,
                "n_atoms": mol.GetNumAtoms(),
                "elements": elements,
                "heavy_atoms": heavy,
                "heavy_bonds": bonds,
                "rdkit_coords": coords,
                "rdkit_heavy_coords": heavy_coords,
                "rdkit_heavy_bond_lengths": bond_lengths,
                "rdkit_heavy_bond_bounds": bond_bounds,
            }
        )
    payload = {
        "reference": f"{REFERENCE} (AllChem.ETKDGv3, randomSeed={ETKDG_SEED})",
        "seed": ETKDG_SEED,
        "tolerance_bond_length": EMBED_BOND_LENGTH_TOL,
        "tolerance_rmsd": EMBED_RMSD_TOL,
        "tolerance_units": (
            "bond length: fraction of the RDKit distance-bounds window; "
            "RMSD: angstrom"
        ),
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F72",
        observable="embedding_bond_length_max_abs",
        reference_tool=f"{REFERENCE} AllChem.EmbedMolecule (ETKDGv3)",
        units="fraction of the bounds window",
        tolerance=EMBED_BOND_LENGTH_TOL,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "the two embedders use independent random streams, so coordinates are "
            "NOT expected to match; the gate checks connectivity, atom count, "
            "heavy-atom bond-length agreement and a bounded Kabsch RMSD"
        ),
    )
    emit_fixture_report("embedding_rdkit.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F73 — ETKDG-3D force-field parity
# ---------------------------------------------------------------------------

#: Tolerance on the ETKDG-3D force-field energy (kcal/mol) evaluated on the
#: shared perturbed geometry.
ETKDG3D_ENERGY_TOL = 5.0

# RDKit's `ETKDGForceConsts` (Code/GraphMol/ForceFieldHelpers/CrystalFF/
# TorsionPreferences.h): kTermAngle = 1.0, kTermImproper = 10.0.
K_TERM_ANGLE = 1.0
K_TERM_IMPROPER = 10.0
KNOWN_DIST_TOL = 0.01
KNOWN_DIST_FORCE_CONSTANT = 100.0
BOUNDS_MAT_FORCE_SCALING = 1.0


def rdkit_angles(mol: Chem.Mol) -> list[tuple[int, int, int, bool]]:
    """RDKit's `collectBondsAndAngles` (BoundsMatrixBuilder.cpp).

    Returns (i, j, k, is_sp) tuples where j is the central atom.
    """
    angles: list[tuple[int, int, int, bool]] = []
    nb = mol.GetNumBonds()
    for bi in range(nb):
        b1 = mol.GetBondWithIdx(bi)
        aid11, aid12 = b1.GetBeginAtomIdx(), b1.GetEndAtomIdx()
        for bj in range(bi + 1, nb):
            b2 = mol.GetBondWithIdx(bj)
            aid21, aid22 = b2.GetBeginAtomIdx(), b2.GetEndAtomIdx()
            if aid11 != aid21 and aid11 != aid22 and aid12 != aid21 and aid12 != aid22:
                continue
            if aid12 == aid21:
                a1, a2, a3 = aid11, aid12, aid22
            elif aid12 == aid22:
                a1, a2, a3 = aid11, aid12, aid21
            elif aid11 == aid21:
                a1, a2, a3 = aid12, aid11, aid22
            else:
                a1, a2, a3 = aid12, aid11, aid21
            is_sp = False
            if b1.GetBondType() == Chem.BondType.TRIPLE or b2.GetBondType() == Chem.BondType.TRIPLE:
                is_sp = True
            elif (
                b1.GetBondType() == Chem.BondType.DOUBLE
                and b2.GetBondType() == Chem.BondType.DOUBLE
                and mol.GetAtomWithIdx(a2).GetDegree() == 2
            ):
                is_sp = True
            angles.append((a1, a2, a3, is_sp))
    return angles


def rdkit_improper_centers(mol: Chem.Mol) -> list[tuple[int, bool]]:
    """RDKit's improper (out-of-plane) centres: C/N/O, SP2, degree 3.

    Returns (center_idx, is_c_bound_to_sp2_o) pairs.
    """
    out = []
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        if z not in (6, 7, 8):
            continue
        if atom.GetHybridization() != Chem.HybridizationType.SP2:
            continue
        if atom.GetDegree() != 3:
            continue
        bound_to_o = any(
            nb.GetAtomicNum() == 8
            and nb.GetHybridization() == Chem.HybridizationType.SP2
            for nb in atom.GetNeighbors()
        )
        out.append((atom.GetIdx(), bool(z == 6 and bound_to_o)))
    return out


def uff_inversion_energy(
    coords: np.ndarray, center: int, neighbors: list[int], is_c_bound_to_o: bool
) -> float:
    """RDKit's `ForceFields::UFF::InversionContribs` energy (UFF/Utils.cpp)."""
    k_base = 50.0 if is_c_bound_to_o else 6.0
    force_constant = k_base / 3.0 * K_TERM_IMPROPER
    c0, c1, c2 = 1.0, -1.0, 0.0
    p = coords
    total = 0.0
    a0, a1, a2 = neighbors
    # RDKit's three permutations: (0,1,2,3), (0,1,3,2), (2,1,3,0) of
    # (a0, center, a1, a2).
    for (i1, i3, i4) in ((a0, a1, a2), (a0, a2, a1), (a1, a2, a0)):
        p1, p2, p3, p4 = p[i1], p[center], p[i3], p[i4]
        rji = p1 - p2
        rjk = p3 - p2
        rjl = p4 - p2
        l2ji = float(rji.dot(rji))
        l2jk = float(rjk.dot(rjk))
        l2jl = float(rjl.dot(rjl))
        if l2ji < 1e-16 or l2jk < 1e-16 or l2jl < 1e-16:
            continue
        n = np.cross(rji, rjk)
        n = n / (np.sqrt(l2ji) * np.sqrt(l2jk))
        l2n = float(n.dot(n))
        if l2n < 1e-16:
            continue
        cos_y = float(n.dot(rjl)) / (np.sqrt(l2jl) * np.sqrt(l2n))
        cos_y = max(-1.0, min(1.0, cos_y))
        sin_y = np.sqrt(max(0.0, 1.0 - cos_y * cos_y))
        cos_2w = 2.0 * sin_y * sin_y - 1.0
        total += force_constant * (c0 + c1 * sin_y + c2 * cos_2w)
    return total


def build_rdkit_etkdg3d_forcefield(
    mol: Chem.Mol, conf_id: int = 0
):
    """Rebuild RDKit's `construct3DForceField` distance/angle term set.

    Follows `Code/DistGeom/DistGeomUtils.cpp` (RDKit Release_2026_09_1):
      * add12Terms  — bond distances d ± 0.01 with k = 100
      * add13Terms  — angle distances (d ± 0.01, k = 100) unless the central
                      atom is an improper centre, in which case the bounds
                      matrix is used; SP centres get a 179-180 deg angle
                      constraint with k = kTermAngle
      * addLongRangeDistanceConstraints — every remaining pair gets a
                      bounds-matrix distance constraint with k = 10
    """
    n = mol.GetNumAtoms()
    lower, upper = bounds_from_rdkit(mol, smooth=True)
    lb = {}
    ub = {}
    for (i, j), lo, up in zip(pair_list(n), lower, upper):
        lb[(i, j)] = lo
        ub[(i, j)] = up

    marked: set[tuple[int, int]] = set()
    dist12: list[tuple[int, int, float, float, float]] = []
    dist13: list[tuple[int, int, float, float, float]] = []
    long_terms: list[tuple[int, int, float, float, float]] = []
    angle_terms: list[tuple[int, int, int, float, float, float]] = []

    def mark(a: int, b: int) -> None:
        marked.add((a, b) if a < b else (b, a))

    # --- 1-2 terms -------------------------------------------------------
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        mark(i, j)
        d = float(
            np.linalg.norm(
                np.asarray(mol.GetConformer(conf_id).GetAtomPosition(i))
                - np.asarray(mol.GetConformer(conf_id).GetAtomPosition(j))
            )
        )
        dist12.append(
            (i, j, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL, KNOWN_DIST_FORCE_CONSTANT)
        )

    # --- 1-3 terms -------------------------------------------------------
    improper = rdkit_improper_centers(mol)
    improper_centers = {c for c, _ in improper}
    for (i, j, k, is_sp) in rdkit_angles(mol):
        mark(i, k)
        if is_sp:
            angle_terms.append((i, j, k, 179.0, 180.0, K_TERM_ANGLE))
            continue
        if j in improper_centers:
            key = (i, k) if i < k else (k, i)
            dist13.append((i, k, lb[key], ub[key], KNOWN_DIST_FORCE_CONSTANT))
            continue
        d = float(
            np.linalg.norm(
                np.asarray(mol.GetConformer(conf_id).GetAtomPosition(i))
                - np.asarray(mol.GetConformer(conf_id).GetAtomPosition(k))
            )
        )
        dist13.append(
            (i, k, d - KNOWN_DIST_TOL, d + KNOWN_DIST_TOL, KNOWN_DIST_FORCE_CONSTANT)
        )

    # --- long-range terms ------------------------------------------------
    for (i, j) in pair_list(n):
        if (i, j) in marked:
            continue
        long_terms.append(
            (i, j, lb[(i, j)], ub[(i, j)], BOUNDS_MAT_FORCE_SCALING * 10.0)
        )

    dist_terms = dist12 + dist13 + long_terms
    ff = rdForceFieldHelpers.CreateEmptyForceFieldForMol(mol, conf_id)
    for (i, j, lo, hi, k) in dist_terms:
        ff.AddDistanceConstraint(i, j, lo, hi, k)
    for (i, j, k, lo, hi, fc) in angle_terms:
        ff.UFFAddAngleConstraint(i, j, k, False, lo, hi, fc)
    # A second force field with the 1-2/1-3/angle terms only: the long-range
    # term set depends on the experimental-torsion bookkeeping, which differs
    # between RDKit and sci-form, so the two classes are reported separately.
    ff_1213 = rdForceFieldHelpers.CreateEmptyForceFieldForMol(mol, conf_id)
    for (i, j, lo, hi, k) in dist12 + dist13:
        ff_1213.AddDistanceConstraint(i, j, lo, hi, k)
    for (i, j, k, lo, hi, fc) in angle_terms:
        ff_1213.UFFAddAngleConstraint(i, j, k, False, lo, hi, fc)
    ff_long = rdForceFieldHelpers.CreateEmptyForceFieldForMol(mol, conf_id)
    for (i, j, lo, hi, k) in long_terms:
        ff_long.AddDistanceConstraint(i, j, lo, hi, k)
    return (
        ff,
        ff_1213,
        ff_long,
        dist12,
        dist13,
        long_terms,
        angle_terms,
        improper,
        lower,
        upper,
    )


def run_etkdg3d() -> Report:
    cases = []
    for smi in ETKDG3D_CORPUS:
        mol = embed_rdkit(smi)
        if mol is None:
            continue
        # The conformer itself satisfies every constraint (energy ~ 0), so the
        # force field is evaluated on a deterministic perturbation of it: both
        # sides then see identical, non-trivial energies.
        base = conformer_coords(mol)
        eval_coords = perturb(base, PERTURB_SEED + len(cases))
        mol2 = Chem.Mol(mol)
        conf = mol2.GetConformer(0)
        for i in range(mol2.GetNumAtoms()):
            conf.SetAtomPosition(i, tuple(eval_coords[3 * i : 3 * i + 3]))
        (
            ff,
            ff_1213,
            ff_long,
            dist12,
            dist13,
            long_terms,
            angle_terms,
            improper,
            lower,
            upper,
        ) = build_rdkit_etkdg3d_forcefield(mol2, 0)
        coords = np.asarray(eval_coords, dtype=float).reshape(-1, 3)
        improper_energy = 0.0
        improper_terms = []
        for center, bound_o in improper:
            nbrs = [a.GetIdx() for a in mol2.GetAtomWithIdx(center).GetNeighbors()]
            improper_energy += uff_inversion_energy(coords, center, nbrs, bound_o)
            k_base = 50.0 if bound_o else 6.0
            k = k_base / 3.0 * K_TERM_IMPROPER
            a0, a1, a2 = nbrs
            improper_terms.append([a0, center, a1, a2, k])
            improper_terms.append([a0, center, a2, a1, k])
            improper_terms.append([a1, center, a2, a0, k])
        n_bonds = len(dist12)
        n_angles = len(rdkit_angles(mol2))
        cases.append(
            {
                "smiles": smi,
                "n_atoms": mol2.GetNumAtoms(),
                "elements": element_list(mol2),
                "eval_coords": eval_coords,
                "bounds_lower": lower,
                "bounds_upper": upper,
                "rdkit_dist12": [list(t) for t in dist12],
                "rdkit_dist13": [list(t) for t in dist13],
                "rdkit_long": [list(t) for t in long_terms],
                "rdkit_angle": [list(t) for t in angle_terms],
                "rdkit_inversion": improper_terms,
                "rdkit_energy": float(ff.CalcEnergy()),
                "rdkit_energy_12_13": float(ff_1213.CalcEnergy()),
                "rdkit_energy_long": float(ff_long.CalcEnergy()),
                "rdkit_n_dist12": n_bonds,
                "rdkit_n_dist13": len(dist13),
                "rdkit_n_long": len(long_terms),
                "rdkit_n_angle": len(angle_terms),
                "rdkit_n_improper": len(improper),
                "rdkit_improper_energy": improper_energy,
            }
        )
    payload = {
        "reference": f"{REFERENCE} ForceField engine",
        "rdkit_api": (
            "rdkit.Chem.rdForceFieldHelpers.CreateEmptyForceFieldForMol + "
            "AddDistanceConstraint / UFFAddAngleConstraint; term set rebuilt from "
            "RDKit Release_2026_09_1 Code/DistGeom/DistGeomUtils.cpp "
            "(construct3DForceField is not exported to Python)"
        ),
        "sci_form_api": (
            "sci_form::forcefield::etkdg_3d::build_etkdg_3d_ff_with_torsions "
            "+ etkdg_3d_energy_f64"
        ),
        "perturbation_seed": PERTURB_SEED,
        "perturbation_sigma": 0.25,
        "tolerance": ETKDG3D_ENERGY_TOL,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F73",
        observable="etkdg3d_ff_energy_abs_error",
        reference_tool=f"{REFERENCE} ForceField (ETKDG-3D term set)",
        units="kcal/mol",
        tolerance=ETKDG3D_ENERGY_TOL,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "1-2/1-3 distance constraints, SP angle constraints and the UFF "
            "inversion (improper) terms are compared on a shared perturbed "
            "geometry; the long-range constraints are reported but not gated "
            "because sci-form's experimental-torsion bookkeeping differs from "
            "RDKit's"
        ),
    )
    emit_fixture_report("etkdg3d_ff_rdkit.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# F74 — BFGS minimisation parity
# ---------------------------------------------------------------------------

#: Tolerance on the Kabsch RMSD between sci-form's bounds-BFGS minimum and
#: RDKit's ForceField::minimize minimum (Å).  Both minimisers enforce the
#: *same* distance-bounds constraint set, but they use different local
#: optimisers (RDKit's flat-bottom harmonic vs sci-form's bounds-violation
#: energy), so the two minima are close but not identical.
#: Measured max: 0.136 Å (thiophene).
BFGS_RMSD_TOL = 0.25
#: Maximum residual bounds-violation *energy* accepted after minimisation.
#: sci-form's BFGS stops on a gradient criterion (1e-3) instead of driving the
#: violation to zero; the measured worst case is 5.7e-3.
BFGS_ENERGY_TOL = 1e-2
#: Maximum residual distance-bound violation accepted after minimisation (Å).
#: Measured max: 0.083 Å (a 1-5 pair whose sci-form bound is tighter than
#: RDKit's, see F71).
BFGS_VIOLATION_TOL = 0.1


def run_bfgs() -> Report:
    cases = []
    for smi in GEOMETRY_CORPUS:
        mol = embed_rdkit(smi)
        if mol is None:
            continue
        start = conformer_coords(mol)
        # The *unsmoothed* bounds matrix is used on both sides: RDKit's
        # `ForceField::minimize` cannot always converge on the triangle-smoothed
        # matrix (the smoothed lower bounds make the constraint set infeasible
        # for some molecules and RDKit raises "bad direction in linearSearch").
        lower, upper = bounds_from_rdkit(mol, smooth=False)

        # (a) RDKit's own BFGS (`ForceField::minimize`) on flat-bottom distance
        #     constraints built from RDKit's bounds matrix.
        ff = rdForceFieldHelpers.CreateEmptyForceFieldForMol(mol, 0)
        n_bounds = 0
        for (i, j), lo, up in zip(pair_list(mol.GetNumAtoms()), lower, upper):
            if lo <= 0.0 or up >= 1000.0:
                continue
            ff.AddDistanceConstraint(i, j, lo, up, KNOWN_DIST_FORCE_CONSTANT)
            n_bounds += 1
        conv = -1
        rdkit_bfgs_coords = list(start)
        try:
            conv = int(ff.Minimize(500, 1e-4, 1e-6))
            rdkit_bfgs_coords = conformer_coords(mol, 0)
        except Exception:  # pragma: no cover - RDKit BFGS invariant violation
            conv = -1
            rdkit_bfgs_coords = list(start)

        # (b) RDKit's UFF minimiser from the same start (the API named in the
        #     task description).
        mol3 = Chem.Mol(mol)
        uff_before = None
        uff_after = None
        uff_conv = -1
        uff_coords = list(start)
        try:
            uff_ff = rdForceFieldHelpers.UFFGetMoleculeForceField(mol3, 0)
            if uff_ff is not None:
                uff_before = float(uff_ff.CalcEnergy())
                uff_conv = int(rdForceFieldHelpers.UFFOptimizeMolecule(mol3, 1000, 10.0, 0))
                uff_ff2 = rdForceFieldHelpers.UFFGetMoleculeForceField(mol3, 0)
                uff_after = float(uff_ff2.CalcEnergy()) if uff_ff2 is not None else None
                uff_coords = conformer_coords(mol3, 0)
        except Exception:
            pass

        cases.append(
            {
                "smiles": smi,
                "n_atoms": mol.GetNumAtoms(),
                "elements": element_list(mol),
                "start_coords": start,
                "bounds_lower": lower,
                "bounds_upper": upper,
                "rdkit_bfgs_coords": rdkit_bfgs_coords,
                "rdkit_bfgs_converged": conv,
                "rdkit_n_bounds": n_bounds,
                "rdkit_uff_coords": uff_coords,
                "rdkit_uff_converged": uff_conv,
                "rdkit_uff_energy_before": uff_before,
                "rdkit_uff_energy_after": uff_after,
            }
        )
    payload = {
        "reference": f"{REFERENCE} ForceField::minimize (BFGS) + UFFOptimizeMolecule",
        "rdkit_api": (
            "rdForceFieldHelpers.CreateEmptyForceFieldForMol + AddDistanceConstraint "
            "+ Minimize; rdForceFieldHelpers.UFFOptimizeMolecule"
        ),
        "sci_form_api": "sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit",
        "tolerance_rmsd": BFGS_RMSD_TOL,
        "tolerance_energy": BFGS_ENERGY_TOL,
        "tolerance_violation": BFGS_VIOLATION_TOL,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F74",
        observable="bfgs_minimised_rmsd",
        reference_tool=f"{REFERENCE} ForceField::minimize",
        units="angstrom",
        tolerance=BFGS_RMSD_TOL,
        level="L2",
        n_cases=len(cases),
        metric="max_abs",
        value=float("nan"),
        notes=(
            "both sides start from the same coordinates and enforce the same "
            "distance-bounds constraint set; the RMSD measures the separation of "
            "the two local minima reached by the two optimisers"
        ),
    )
    emit_fixture_report("bfgs_opt_rdkit.json", payload, report)
    return report


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs="?",
        default="all",
        choices=["bounds", "embed", "etkdg3d", "bfgs", "all"],
    )
    args = parser.parse_args()
    commands = {
        "bounds": run_bounds,
        "embed": run_embed,
        "etkdg3d": run_etkdg3d,
        "bfgs": run_bfgs,
    }
    if args.command == "all":
        for name, fn in commands.items():
            print(f"--- {name} ---")
            fn()
    else:
        commands[args.command]()
    print(f"\nfixtures -> {REPO / 'tests' / 'fixtures' / 'validation'}")
    print(f"summary  -> {REPO / 'reports' / 'validation' / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
