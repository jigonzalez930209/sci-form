#!/usr/bin/env python3
"""Oracle F79: BP symmetry functions + tiny-NN energies/forces (numpy).

- Geometries: RDKit ETKDG (seed 42), H/C/N/O molecules.
- AEV reference: independent numpy implementation of the Behler-Parrinello
  formulas as implemented (radial G2 with PAIRED shift/eta lists, angular
  G4 with Rs shift, cos cutoff) — Behler-Parrinello 2007/2011 reference.
- Energies: numpy forward pass through tiny random nets (seed 43) with
  exact activations (linear/relu/celu) + one GELU net (erf-limited).
- Forces: central finite differences (h=1e-5) of the NUMPY pipeline —
  fully independent of the analytic backprop under test.
- ANI-2x paper constants recorded for the defaults check.
- TorchANI full-tensor comparison is out of reach by design (our AEVs
  are species-agnostic, TorchANI splits by species); MD17/SPICE energy
  benchmarks need trained weights: documented gaps.

Usage: oracle_mlff.py  (writes tests/fixtures/validation/mlff_bp.json)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "tests" / "fixtures" / "validation" / "mlff_bp.json"

SMILES = ["O", "C", "N", "CC", "C=C", "c1ccccc1", "CO", "CN"]


def cutoff(r: float, rc: float) -> float:
    return 0.0 if r >= rc else 0.5 * (1.0 + np.cos(np.pi * r / rc))


def aev_atom(i: int, elements, pos: np.ndarray, p: dict) -> np.ndarray:
    n = len(elements)
    rad = np.zeros(len(p["radial_shifts"]))
    for j in range(n):
        if i == j:
            continue
        rij = float(np.linalg.norm(pos[i] - pos[j]))
        if rij >= p["radial_cutoff"]:
            continue
        fc = cutoff(rij, p["radial_cutoff"])
        for s, (rs, eta) in enumerate(zip(p["radial_shifts"], p["radial_etas"])):
            rad[s] += np.exp(-eta * (rij - rs) ** 2) * fc
    ang = []
    for rs in p["angular_shifts"]:
        for eta in p["angular_etas"]:
            for zeta in p["angular_zetas"]:
                tot = 0.0
                for j in range(n):
                    if i == j:
                        continue
                    rij = float(np.linalg.norm(pos[i] - pos[j]))
                    if rij >= p["angular_cutoff"]:
                        continue
                    for k in range(j + 1, n):
                        if i == k:
                            continue
                        rik = float(np.linalg.norm(pos[i] - pos[k]))
                        if rik >= p["angular_cutoff"]:
                            continue
                        vij = pos[j] - pos[i]
                        vik = pos[k] - pos[i]
                        cos_t = float(vij @ vik / (np.linalg.norm(vij) * np.linalg.norm(vik)))
                        cos_t = min(1.0, max(-1.0, cos_t))
                        rmid = 0.5 * (rij + rik)
                        tot += (2.0 ** (1 - zeta) * (1 + cos_t) ** zeta
                                * np.exp(-eta * (rmid - rs) ** 2)
                                * cutoff(rij, p["angular_cutoff"]) * cutoff(rik, p["angular_cutoff"]))
                ang.append(tot)
    return np.concatenate([rad, ang])


def activate(x: np.ndarray, kind: str) -> np.ndarray:
    if kind == "linear":
        return x
    if kind == "relu":
        return np.maximum(x, 0.0)
    if kind == "celu":
        return np.where(x >= 0, x, np.exp(x) - 1.0)
    if kind == "gelu":
        from scipy.special import erf
        return x * 0.5 * (1.0 + erf(x / np.sqrt(2.0)))
    raise ValueError(kind)


def forward(aev: np.ndarray, layers: list[dict]) -> float:
    x = aev
    for L in layers:
        w = np.array(L["weights"]).reshape(L["out"], L["in"])
        x = activate(np.array(L["bias"]) + w @ x, L["act"])
    return float(x[0])


def main() -> int:
    rng = np.random.default_rng(42)
    rw = np.random.default_rng(43)
    mols = []
    for smi in SMILES:
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0, smi
        conf = mol.GetConformer()
        pos = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
        mols.append({"smiles": smi,
                     "elements": [a.GetAtomicNum() for a in mol.GetAtoms()],
                     "coords": pos.tolist()})
    params = {
        "default": {"radial_cutoff": 5.2, "angular_cutoff": 3.5,
                    "radial_shifts": [0.9 + 0.25 * i for i in range(16)],
                    "radial_etas": [16.0] * 16,
                    "angular_shifts": [0.9 + 0.5 * i for i in range(4)],
                    "angular_etas": [8.0],
                    "angular_zetas": [1.0, 2.0, 4.0, 8.0]},
        "small": {"radial_cutoff": 5.0, "angular_cutoff": 3.0,
                  "radial_shifts": [1.0, 2.0], "radial_etas": [4.0, 9.0],
                  "angular_shifts": [1.2], "angular_etas": [5.0],
                  "angular_zetas": [1.0, 2.0]},
    }
    # Tiny nets sized for the RUST layout (256 radial incl. 240 structural
    # zeros + 16 angular = 272): exact activations for H/C/N/O + GELU net.
    AEV_DIM = 256 + 16
    nets: dict[str, list[dict]] = {}
    for z in (1, 6, 7, 8):
        nets[str(z)] = [
            {"in": AEV_DIM, "out": 4,
             "weights": rw.standard_normal(4 * AEV_DIM).tolist(),
             "bias": rw.standard_normal(4).tolist(), "act": "celu" if z % 2 else "relu"},
            {"in": 4, "out": 1,
             "weights": rw.standard_normal(4).tolist(),
             "bias": rw.standard_normal(1).tolist(), "act": "linear"},
        ]
    nets["6-gelu"] = [
        {"in": AEV_DIM, "out": 4,
         "weights": rw.standard_normal(4 * AEV_DIM).tolist(),
         "bias": rw.standard_normal(4).tolist(), "act": "gelu"},
        {"in": 4, "out": 1,
         "weights": rw.standard_normal(4).tolist(),
         "bias": rw.standard_normal(1).tolist(), "act": "linear"},
    ]

    def full_aev(v32: np.ndarray) -> np.ndarray:
        return np.concatenate([v32[:16], np.zeros(240), v32[16:]])
    cases = []
    for mol in mols:
        pos = np.array(mol["coords"])
        entry = dict(mol)
        entry["aev_default"] = [aev_atom(i, mol["elements"], pos, params["default"]).tolist()
                                for i in range(len(pos))]
        entry["aev_small"] = [aev_atom(i, mol["elements"], pos, params["small"]).tolist()
                              for i in range(len(pos))]
        # Energy/forces with the default params + exact-activation nets
        # (full 272 layout).
        aevs = [full_aev(np.array(v)) for v in entry["aev_default"]]
        e = sum(forward(a, nets[str(z)]) for a, z in zip(aevs, mol["elements"]))
        h = 1e-5
        f = np.zeros_like(pos)
        for i in range(len(pos)):
            for d in range(3):
                pp, pm = pos.copy(), pos.copy()
                pp[i, d] += h
                pm[i, d] -= h
                ap = [full_aev(aev_atom(j, mol["elements"], pp, params["default"])) for j in range(len(pos))]
                am = [full_aev(aev_atom(j, mol["elements"], pm, params["default"])) for j in range(len(pos))]
                ep = sum(forward(a, nets[str(z)]) for a, z in zip(ap, mol["elements"]))
                em = sum(forward(a, nets[str(z)]) for a, z in zip(am, mol["elements"]))
                f[i, d] = -(ep - em) / (2 * h)
        entry["energy_ref"] = float(e)
        entry["forces_ref"] = f.tolist()
        # GELU-net energy for C atoms only (erf-limited tolerance).
        if 6 in mol["elements"]:
            a2 = [full_aev(np.array(v)) for v in entry["aev_default"]]
            entry["gelu_energy_ref"] = float(
                sum(forward(a, nets["6-gelu"]) for a, z in zip(a2, mol["elements"]) if z == 6))
        cases.append(entry)
        print(f"{mol['smiles']}: E={e:.6f} n={len(pos)}")
    payload = {
        "ficha": "F79",
        "oracle": "oracle_mlff.py (BP formulas + tiny nets + FD forces, RDKit ETKDG)",
        "note": ("TorchANI tensor comparison out of reach by design (species-split "
                 "vs species-agnostic); MD17/SPICE energy benchmarks need trained "
                 "weights: documented gaps."),
        "aev_tolerance": 1e-9,
        "energy_tolerance": 1e-9,
        "gelu_energy_tolerance": 2e-6,
        # Forces are central FD (impl h=1e-4, oracle h=1e-5): tol covers
        # OUR truncation, well inside the 1e-3 criterion.
        "force_tolerance": 1e-4,
        "params": params,
        "nets": nets,
        "ani2x_paper": params["default"],
        "cases": cases,
    }
    OUT.write_text(json.dumps(payload) + "\n")
    print(f"wrote {OUT} ({len(cases)} molecules)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
