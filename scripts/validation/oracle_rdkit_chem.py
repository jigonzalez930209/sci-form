#!/usr/bin/env python3
"""RDKit oracle for sci-form cross-validation (fichas F81, F82, F83).

F81 — CIP stereochemistry (R/S at stereocentres, E/Z at double bonds) against
      RDKit's CIP labeller.  Both sides derive the label from the *same* 3-D
      geometry, so the comparison isolates the R/S sign convention and the
      CIP priority implementation.  The fixture also carries the label RDKit
      assigns from the SMILES `@`/`@@` tag alone, which documents whether
      sci-form consumes the tag.
F82 — ML descriptors (`src/ml/descriptors.rs::compute_descriptors`) and the
      Ertl TPSA (`src/ml/ensemble.rs::compute_tpsa`) against RDKit.
F83 — Canonical SMILES (`src/canonical_smiles.rs::to_canonical_smiles`)
      against `Chem.MolToSmiles(canonical=True)`.

Each subcommand writes a JSON fixture under tests/fixtures/validation/ and
emits a report row into reports/validation/summary.json.  The Rust side
(tests/rdkit_chem_parity.rs) recomputes the sci-form values from the same
fixtures and asserts the tolerances, so CI never needs RDKit.

Usage:
    python scripts/validation/oracle_rdkit_chem.py stereo|descriptors|canonical|all
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from common import (
    FIXTURE_DIR,
    REPO,
    Report,
    emit_fixture_report,
    read_smiles_corpus,
)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
except ImportError as exc:  # pragma: no cover
    print(f"RDKit is required for this oracle: {exc}", file=sys.stderr)
    sys.exit(1)

RDKIT_VERSION = Chem.rdBase.rdkitVersion
REFERENCE = f"RDKit {RDKIT_VERSION}"

BOHR_PER_ANGSTROM = 1.0 / 0.529177210903


# ═══════════════════════════════════════════════════════════════════════════
# helpers
# ═══════════════════════════════════════════════════════════════════════════


def embed(smiles: str, seed: int = 0x5EED) -> tuple[list[int], list[float]]:
    """Deterministic ETKDGv3 embedding (heavy atoms + explicit H)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = False
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError(f"embedding failed for {smiles}")
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
    conf = mol.GetConformer()
    elements = [a.GetAtomicNum() for a in mol.GetAtoms()]
    coords: list[float] = []
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        coords += [p.x, p.y, p.z]
    return elements, coords


def _label_from_stereo_info(mol: Chem.Mol) -> list[list]:
    """CIP labels (R/S per centre, E/Z per double bond) of `mol`.

    `mol` must already carry stereo information (either from the SMILES tag
    or assigned from a 3-D conformer).
    """
    centers: list[list] = []
    bonds: list[list] = []
    for info in Chem.FindPotentialStereo(mol):
        if info.type != Chem.StereoType.Atom_Tetrahedral:
            continue
        if info.specified != Chem.StereoSpecified.Specified:
            continue
        idx = info.centeredOn
        label = None
        for center_idx, tag in Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False
        ):
            if center_idx == idx:
                label = tag
                break
        if label in ("R", "S"):
            centers.append([idx, label])
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo = bond.GetStereo()
        if stereo == Chem.BondStereo.STEREOE:
            bonds.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), "E"])
        elif stereo == Chem.BondStereo.STEREOZ:
            bonds.append([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), "Z"])
    centers.sort()
    bonds.sort()
    return centers, bonds


def _assign_from_structure(mol: Chem.Mol) -> Chem.Mol:
    """Copy of `mol` with stereo assigned from its 3-D conformer."""
    work = Chem.Mol(mol)
    Chem.AssignAtomChiralTagsFromStructure(work, confId=0, replaceExistingTags=True)
    Chem.SetBondStereoFromDirections(work)
    Chem.AssignStereochemistry(work, cleanIt=True, force=True)
    return work


# ═══════════════════════════════════════════════════════════════════════════
# F81 — CIP stereochemistry
# ═══════════════════════════════════════════════════════════════════════════

# name, 3-D source SMILES, expected structural feature
STEREO_CORPUS = [
    ("alanine", "C[C@H](N)C(=O)O", "center"),
    ("alanine-mirror", "C[C@H](N)C(=O)O", "center", "mirror"),
    ("2-butanol", "CC[C@H](C)O", "center"),
    ("lactic-acid", "C[C@@H](O)C(=O)O", "center"),
    ("bromo-chlorofluoromethane", "FC(Cl)Br", "center"),
    ("two-centres", "C[C@H](O)[C@@H](C)O", "center"),
    ("two-centres-mirror", "C[C@H](O)[C@@H](C)O", "center", "mirror"),
    ("cysteine", "N[C@@H](CS)C(=O)O", "center"),
    ("glutamic-acid", "N[C@@H](CCC(=O)O)C(=O)O", "center"),
    ("trans-2-butene", "C/C=C/C", "bond"),
    ("cis-2-butene", "C/C=C\\C", "bond"),
    ("trans-2-pentene", "CC/C=C/CC", "bond"),
    ("cis-2-pentene", "CC/C=C\\CC", "bond"),
    ("1,2-difluoroethene-trans", "F/C=C/F", "bond"),
    ("1,2-difluoroethene-cis", "F/C=C\\F", "bond"),
    ("trans-1,3-pentadiene", "C/C=C/C=C/C", "bond"),
    ("allene-penta-2,3-diene", "CC=C=CC", "allene"),
    ("allene-hexa-2,3-diene", "CCC=C=CC", "allene"),
    ("atropisomer-dimethylbiphenyl", "Cc1ccccc1-c1ccccc1C", "atropisomer"),
    (
        "atropisomer-dinitrobiphenyl",
        "O=[N+]([O-])c1ccccc1-c1ccccc1[N+](=O)[O-]",
        "atropisomer",
    ),
]


def run_stereo() -> Report:
    cases = []
    for entry in STEREO_CORPUS:
        name, smiles = entry[0], entry[1]
        feature = entry[2]
        mirror = len(entry) > 3 and entry[3] == "mirror"

        elements, coords = embed(smiles)
        if mirror:
            # Mirror the geometry through the yz plane: the 3-D configuration
            # is inverted while the SMILES tag stays the same.  A 3-D driven
            # analyser must report the *inverted* label, a tag driven one the
            # original.
            mirrored: list[float] = []
            for i in range(len(coords) // 3):
                mirrored += [-coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]]
            coords = mirrored

        heavy = Chem.MolFromSmiles(smiles)
        heavy_h = Chem.AddHs(heavy)
        with_h = Chem.Mol(heavy_h)
        with_h.RemoveAllConformers()
        conf = Chem.Conformer(with_h.GetNumAtoms())
        for i in range(with_h.GetNumAtoms()):
            p = coords[3 * i : 3 * i + 3]
            conf.SetAtomPosition(i, (p[0], p[1], p[2]))
        with_h.AddConformer(conf)
        assert with_h.GetNumAtoms() == len(elements)

        from_3d = _assign_from_structure(with_h)
        centers_3d, bonds_3d = _label_from_stereo_info(from_3d)

        # Labels RDKit assigns from the SMILES tag alone (no 3-D involved).
        tagged = Chem.MolFromSmiles(smiles)
        Chem.AssignStereochemistry(tagged, cleanIt=True, force=True)
        centers_tag, bonds_tag = _label_from_stereo_info(tagged)

        cases.append(
            {
                "name": name,
                "smiles": smiles,
                "feature": feature,
                "mirrored": mirror,
                "elements": elements,
                "coords": coords,
                "rdkit_centers_from_3d": centers_3d,
                "rdkit_bonds_from_3d": bonds_3d,
                "rdkit_centers_from_tag": centers_tag,
                "rdkit_bonds_from_tag": bonds_tag,
            }
        )

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F81",
        observable="cip_stereo_labels_vs_rdkit",
        reference_tool=f"{REFERENCE} CIP labeller",
        units="count",
        tolerance=0,
        level="L2",
        n_cases=len(cases),
        metric="mismatched_labels",
        value=float("nan"),
        notes=(
            "R/S and E/Z derived from the same 3-D geometry on both sides; "
            "the fixture also stores RDKit's SMILES-tag-only labels."
        ),
    )
    emit_fixture_report("rdkit_stereo.json", payload, report)
    return report


# ═══════════════════════════════════════════════════════════════════════════
# F82 — ML descriptors
# ═══════════════════════════════════════════════════════════════════════════


def _wiener_index(mol: Chem.Mol) -> float:
    """Wiener index from RDKit's topological distance matrix (heavy atoms)."""
    dist = Chem.GetDistanceMatrix(mol, useBO=False, useAtomWts=False, force=True)
    n = len(dist)
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            total += float(dist[i][j])
    return total


def _descriptor_corpus(limit: int) -> list[str]:
    """Curated SMILES bank used for the descriptor comparison.

    The curated `diverse_molecules.smi` bank is used (every entry is parsable
    by sci-form — see the Track G guard in tests/cross_validation.rs) and is
    topped up with the ChEMBL sample.
    """
    return read_smiles_corpus(limit=limit)


def run_descriptors(limit: int = 260) -> Report:
    cases = []
    skipped = 0
    for smi in _descriptor_corpus(limit):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            skipped += 1
            continue
        if len(Chem.GetMolFrags(mol)) > 1:
            skipped += 1  # sci-form keeps the largest fragment by design
            continue
        with_h = Chem.AddHs(mol)
        n_h = with_h.GetNumAtoms() - mol.GetNumAtoms()
        try:
            tpsa = rdMolDescriptors.CalcTPSA(mol)
            tpsa_sp = rdMolDescriptors.CalcTPSA(mol, False, True)
            balaban = Descriptors.BalabanJ(mol)
        except Exception:
            skipped += 1
            continue
        cases.append(
            {
                "name": smi,
                "smiles": smi,
                "n_heavy": mol.GetNumAtoms(),
                "n_hydrogens": n_h,
                "elements_with_h": [a.GetAtomicNum() for a in with_h.GetAtoms()],
                "n_bonds_heavy": mol.GetNumBonds(),
                "n_bonds_with_h": with_h.GetNumBonds(),
                "mol_wt": Descriptors.MolWt(mol),
                "heavy_atom_count": Descriptors.HeavyAtomCount(mol),
                "num_hbd": Descriptors.NumHDonors(mol),
                "num_hba": Descriptors.NumHAcceptors(mol),
                "num_rotatable_strict": int(
                    rdMolDescriptors.CalcNumRotatableBonds(
                        mol, rdMolDescriptors.NumRotatableBondsOptions.Strict
                    )
                ),
                "num_rotatable_nonstrict": int(
                    rdMolDescriptors.CalcNumRotatableBonds(
                        mol, rdMolDescriptors.NumRotatableBondsOptions.NonStrict
                    )
                ),
                "fraction_csp3": Descriptors.FractionCSP3(mol),
                "num_rings": int(rdMolDescriptors.CalcNumRings(mol)),
                "n_sssr": len(Chem.GetSSSR(mol)),
                "num_aromatic_rings": int(Descriptors.NumAromaticRings(mol)),
                "num_aromatic_atoms": int(
                    sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
                ),
                "wiener_index": _wiener_index(mol),
                "balaban_j": balaban,
                "tpsa": tpsa,
                "tpsa_with_s_p": tpsa_sp,
            }
        )

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "skipped": skipped,
        "cases": cases,
    }
    report = Report(
        ficha="F82",
        observable="ml_descriptors_vs_rdkit",
        reference_tool=f"{REFERENCE} Descriptors / rdMolDescriptors",
        units="mixed",
        tolerance=0.0,
        level="L2",
        n_cases=len(cases),
        metric="mismatched_cases",
        value=float("nan"),
        notes="descriptor-by-descriptor comparison; tolerances live with the measured rows",
    )
    emit_fixture_report("rdkit_descriptors.json", payload, report)
    return report


# ═══════════════════════════════════════════════════════════════════════════
# F83 — canonical SMILES
# ═══════════════════════════════════════════════════════════════════════════

CANONICAL_CORPUS = [
    ("ethanol", "CCO"),
    ("acetic-acid", "CC(=O)O"),
    ("benzene", "c1ccccc1"),
    ("toluene", "Cc1ccccc1"),
    ("phenol", "Oc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("branching", "CC(C)(C)C"),
    ("isobutane", "CC(C)C"),
    ("neopentane", "CC(C)(C)C"),
    ("long-chain", "CCCCCCCCCC"),
    ("symmetric-diol", "OCC(O)CO"),
    ("symmetric-ester", "CC(=O)OC(C)=O"),
    ("two-rings", "C1CCC2(CC1)CCCCC2"),
    ("ring-numbering", "C1CC2CCC1CC2"),
    ("fused-rings", "c1ccc2c(c1)ccc1ccccc12"),
    ("spiro", "C1CCC2(CC1)CCCC2"),
    ("nitro", "C[N+](=O)[O-]"),
    ("amide", "CC(=O)NC"),
    ("sulfonamide", "CS(=O)(=O)N"),
    ("halides", "FC(Cl)Br"),
    ("ether", "CCOCC"),
    ("ester", "CCOC(=O)C"),
    ("amine", "CCN(C)C"),
    ("phenol-ether", "COc1ccccc1O"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("caffeine", "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("ibuprofen", "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
    ("glucose", "OCC1OC(O)C(O)C(O)C1O"),
    ("biphenyl", "c1ccc(-c2ccccc2)cc1"),
    ("triple-bond", "C#CC#C"),
    ("allene", "CC=C=CC"),
    ("cyclohexane", "C1CCCCC1"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
]


def _reordered_smiles(mol: Chem.Mol, seed: int = 0xC0FFEE) -> str | None:
    """A SMILES with the same molecule but a different atom ordering."""
    n = mol.GetNumAtoms()
    if n < 2:
        return None
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    work = Chem.RenumberAtoms(mol, perm.tolist())
    return Chem.MolToSmiles(work, canonical=True)


def run_canonical() -> Report:
    cases = []
    for name, smi in CANONICAL_CORPUS:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        canonical = Chem.MolToSmiles(mol, canonical=True)
        reordered = [
            _reordered_smiles(mol, seed) for seed in (0xC0FFEE, 0xBEEF, 0x1234)
        ]
        reordered = [r for r in reordered if r]
        # Round-trip check: RDKit's canonical form re-parses to the same
        # canonical form (used as the ground truth for "same molecule").
        cases.append(
            {
                "name": name,
                "smiles": smi,
                "rdkit_canonical": canonical,
                "rdkit_canonical_reordered": reordered[0] if reordered else None,
                "rdkit_canonical_reorderings": reordered,
                "n_heavy_atoms": mol.GetNumAtoms(),
            }
        )

    payload = {
        "reference": REFERENCE,
        "generated": time.strftime("%Y-%m-%d"),
        "cases": cases,
    }
    report = Report(
        ficha="F83",
        observable="canonical_smiles_vs_rdkit",
        reference_tool=f"{REFERENCE} MolToSmiles(canonical=True)",
        units="count",
        tolerance=0,
        level="L2",
        n_cases=len(cases),
        metric="mismatched_cases",
        value=float("nan"),
        notes="canonical-string equality plus round-trip / reordering invariance",
    )
    emit_fixture_report("rdkit_canonical_smiles.json", payload, report)
    return report


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs="?",
        default="all",
        choices=["stereo", "descriptors", "canonical", "all"],
    )
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    if args.command in ("stereo", "all"):
        print("--- stereo (F81) ---")
        run_stereo()
    if args.command in ("descriptors", "all"):
        print("--- descriptors (F82) ---")
        run_descriptors(args.limit or 260)
    if args.command in ("canonical", "all"):
        print("--- canonical (F83) ---")
        run_canonical()

    print(f"\nfixtures -> {FIXTURE_DIR}")
    print(f"summary  -> {REPO / 'reports' / 'validation' / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
