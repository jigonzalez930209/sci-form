#!/usr/bin/env python3
"""
Layered reaction comparison harness.

This script keeps the reaction benchmark honest by splitting it into four layers:
1. Reaction semantics: sci-form vs RDKit on SMIRKS transforms.
2. Path workflow smoke tests: sci-form simplified NEB (any method) vs optional
     geomeTRIC NEB and ASE NEB with tblite on conformer-derived endpoint pairs.
3. Common-geometry energetics: configurable methods on the same endpoint geometries.
4. Cross-validation: sci-form GFN1/GFN2 vs tblite GFN1/GFN2 and sci-form HF-3c
     vs PySCF HF/STO-3G single-point energies on identical geometries.

Optional external dependencies:
- RDKit for SMIRKS semantics.
- geomeTRIC for external NEB optimization, wired through a sci-form CLI engine.
- tblite + ASE for truly external GFN1/GFN2 NEB and single-point validation.
- PySCF for ab-initio HF reference single-point energies.

Examples:
    python scripts/compare_reaction_layers.py
    python scripts/compare_reaction_layers.py --methods uff,pm3,xtb
    python scripts/compare_reaction_layers.py --methods uff,pm3,xtb,gfn2 --external
    python scripts/compare_reaction_layers.py --json
    python scripts/compare_reaction_layers.py --cli './target/release/sci-form'
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import logging
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from compare_smirks_reactions import compare_reactions, generate_report


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CLI = [
    "cargo",
    "run",
    "-p",
    "sci-form-cli",
    "--features",
    "alpha-gsm,alpha-kinetics",
    "--",
]


@dataclass(frozen=True)
class PathBenchmarkCase:
    name: str
    smiles: str
    seed_a: int
    seed_b: int
    atom_symbols: List[str]
    methods: List[str]
    n_images: int = 5


PATH_TS_CASES: List[PathBenchmarkCase] = [
    PathBenchmarkCase(
        name="n_butane_rotamers",
        smiles="C(C(C(C([H])([H])[H])([H])[H])([H])[H])([H])([H])[H]",
        seed_a=2,
        seed_b=5,
        atom_symbols=["C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
        methods=["uff", "pm3"],
    ),
    PathBenchmarkCase(
        name="n_propanol_rotamers",
        smiles="C(C(C(O[H])([H])[H])([H])[H])([H])([H])[H]",
        seed_a=2,
        seed_b=3,
        atom_symbols=["C", "C", "C", "O", "H", "H", "H", "H", "H", "H", "H", "H"],
        methods=["uff", "pm3"],
    ),
]


HARTREE_TO_KCAL_MOL = 627.5094740631
KCAL_MOL_TO_HARTREE = 1.0 / HARTREE_TO_KCAL_MOL
BOHR_TO_ANGSTROM = 0.529177210903
KCAL_MOL_ANGSTROM_TO_HARTREE_BOHR = BOHR_TO_ANGSTROM / HARTREE_TO_KCAL_MOL


def detect_cli(cli_arg: Optional[str]) -> List[str]:
    if cli_arg:
        return cli_arg.split()
    for candidate in (ROOT / "target" / "debug" / "sci-form", ROOT / "target" / "release" / "sci-form"):
        if candidate.exists():
            return [str(candidate)]
    if shutil.which("cargo"):
        return DEFAULT_CLI
    raise RuntimeError("No CLI launcher available. Pass --cli or install cargo.")


def cli_json(cli_prefix: List[str], args: List[str]) -> Dict[str, Any]:
    command = cli_prefix + args
    completed = subprocess.run(
        command,
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        return {
            "success": False,
            "command": command,
            "stdout": completed.stdout,
            "stderr": completed.stderr.strip(),
        }
    try:
        return {"success": True, "payload": json.loads(completed.stdout)}
    except json.JSONDecodeError as exc:
        return {
            "success": False,
            "command": command,
            "stdout": completed.stdout,
            "stderr": f"invalid JSON output: {exc}",
        }


def interpolate_images(start_coords: Sequence[float], end_coords: Sequence[float], n_images: int) -> List[List[float]]:
    images: List[List[float]] = []
    for image_index in range(n_images):
        fraction = image_index / (n_images - 1)
        images.append(
            [
                (1.0 - fraction) * start_coords[idx] + fraction * end_coords[idx]
                for idx in range(len(start_coords))
            ]
        )
    return images


def summarize_neb_images(images: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    energies = [float(image["potential_energy_kcal_mol"]) for image in images]
    peak_index = max(range(len(energies)), key=energies.__getitem__)
    return {
        "n_images": len(images),
        "barrier_kcal_mol": energies[peak_index] - energies[0],
        "reaction_energy_kcal_mol": energies[-1] - energies[0],
        "peak_image_index": peak_index,
        "start_energy_kcal_mol": energies[0],
        "peak_energy_kcal_mol": energies[peak_index],
        "end_energy_kcal_mol": energies[-1],
    }


def align_coords_kabsch(start_coords: Sequence[float], end_coords: Sequence[float]) -> List[float]:
    try:
        import numpy as np
    except Exception as exc:
        raise RuntimeError(f"NumPy unavailable for endpoint alignment: {exc}") from exc

    start_xyz = np.asarray(start_coords, dtype=float).reshape(-1, 3)
    end_xyz = np.asarray(end_coords, dtype=float).reshape(-1, 3)
    start_centroid = start_xyz.mean(axis=0)
    end_centroid = end_xyz.mean(axis=0)
    start_centered = start_xyz - start_centroid
    end_centered = end_xyz - end_centroid
    covariance = end_centered.T @ start_centered
    u_matrix, _, vt_matrix = np.linalg.svd(covariance)
    rotation = u_matrix @ vt_matrix
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1
        rotation = u_matrix @ vt_matrix
    aligned = end_centered @ rotation + start_centroid
    return aligned.reshape(-1).tolist()


def rmsd_between_coords(start_coords: Sequence[float], end_coords: Sequence[float]) -> float:
    try:
        import numpy as np
    except Exception as exc:
        raise RuntimeError(f"NumPy unavailable for RMSD: {exc}") from exc

    start_xyz = np.asarray(start_coords, dtype=float).reshape(-1, 3)
    end_xyz = np.asarray(end_coords, dtype=float).reshape(-1, 3)
    deltas = start_xyz - end_xyz
    return float(np.sqrt((deltas * deltas).sum() / start_xyz.shape[0]))


def prepare_path_case(case: PathBenchmarkCase) -> Dict[str, Any]:
    try:
        import sci_form
    except Exception as exc:
        return {
            "success": False,
            "case": case.name,
            "stderr": f"sci_form unavailable for endpoint preparation: {exc}",
        }

    start = sci_form.embed(case.smiles, case.seed_a)
    end = sci_form.embed(case.smiles, case.seed_b)
    if getattr(start, "error", None) or getattr(end, "error", None):
        return {
            "success": False,
            "case": case.name,
            "stderr": f"embed failed: {getattr(start, 'error', None) or getattr(end, 'error', None)}",
        }

    aligned_end = align_coords_kabsch(start.coords, end.coords)
    return {
        "success": True,
        "name": case.name,
        "smiles": case.smiles,
        "start_coords": list(start.coords),
        "end_coords": aligned_end,
        "methods": case.methods,
        "n_images": case.n_images,
        "endpoint_rmsd_angstrom": rmsd_between_coords(start.coords, aligned_end),
        "selection_note": (
            f"curated sci-form embed seeds {case.seed_a} and {case.seed_b}, then rigidly aligned before path interpolation"
        ),
        "atom_symbols": case.atom_symbols,
    }


def direct_cli_binary(cli_prefix: List[str]) -> Optional[List[str]]:
    if len(cli_prefix) == 1:
        return cli_prefix
    for candidate in (ROOT / "target" / "debug" / "sci-form", ROOT / "target" / "release" / "sci-form"):
        if candidate.exists():
            return [str(candidate)]
    return None


def cli_uff_energy(cli_prefix: List[str], smiles: str, coords: Sequence[float]) -> float:
    result = cli_json(cli_prefix, ["uff", smiles, json.dumps(list(coords))])
    if not result.get("success"):
        raise RuntimeError(result.get("stderr") or "UFF CLI call failed")
    return float(result["payload"]["energy"])


def cli_neb_energy(cli_prefix: List[str], smiles: str, coords: Sequence[float], method: str = "uff") -> float:
    """Single-point energy via the neb-energy CLI command with any backend."""
    result = cli_json(cli_prefix, ["neb-energy", smiles, json.dumps(list(coords)), "--method", method])
    if not result.get("success"):
        raise RuntimeError(result.get("stderr") or f"neb-energy ({method}) failed")
    return float(result["payload"]["energy_kcal_mol"])


def cli_neb_gradient(cli_prefix: List[str], smiles: str, coords: Sequence[float], method: str = "uff") -> tuple:
    """Energy + gradient via the neb-gradient CLI command with any backend."""
    result = cli_json(cli_prefix, ["neb-gradient", smiles, json.dumps(list(coords)), "--method", method])
    if not result.get("success"):
        raise RuntimeError(result.get("stderr") or f"neb-gradient ({method}) failed")
    payload = result["payload"]
    return float(payload["energy_kcal_mol"]), payload["gradient_kcal_mol_ang"]


def run_sci_form_neb(cli_prefix: List[str], prepared: Dict[str, Any], method: str = "uff") -> Dict[str, Any]:
    result = cli_json(
        cli_prefix,
        [
            "simplified-neb-path",
            prepared["smiles"],
            json.dumps(prepared["start_coords"]),
            json.dumps(prepared["end_coords"]),
            "--n-images",
            str(prepared["n_images"]),
            "--n-iter",
            "20",
            "--spring-k",
            "0.01",
            "--step-size",
            "1e-6",
            "--method",
            method,
        ],
    )
    if result.get("success"):
        result["summary"] = summarize_neb_images(result["payload"]["images"])
    return result


def run_geometric_neb(cli_prefix: List[str], prepared: Dict[str, Any], method: str = "uff") -> Dict[str, Any]:
    direct_cli = direct_cli_binary(cli_prefix)
    if direct_cli is None:
        return {
            "success": False,
            "stderr": "geomeTRIC adapter requires a built sci-form CLI binary; run cargo build first or pass --cli ./target/debug/sci-form",
        }

    try:
        import numpy as np
        import geometric.engine
        import geometric.molecule
        import geometric.neb as geometric_neb
    except Exception as exc:
        return {"success": False, "stderr": f"geomeTRIC/NumPy unavailable: {exc}"}

    class SciFormCliEngine(geometric.engine.Engine):
        def __init__(self, molecule: Any, smiles: str, cli_cmd: List[str], engine_method: str = "uff", fd_step: float = 1e-3):
            super().__init__(molecule)
            self.smiles = smiles
            self.cli_cmd = cli_cmd
            self.engine_method = engine_method
            self.fd_step = fd_step

        def calc_new(self, coords: Any, dirname: str) -> Dict[str, Any]:
            angstrom = np.asarray(coords, dtype=float) * BOHR_TO_ANGSTROM
            try:
                energy_kcal, grad_flat = cli_neb_gradient(self.cli_cmd, self.smiles, angstrom.tolist(), self.engine_method)
                gradient = np.asarray(grad_flat, dtype=float)
            except RuntimeError:
                energy_kcal = cli_neb_energy(self.cli_cmd, self.smiles, angstrom.tolist(), self.engine_method)
                gradient = np.zeros_like(angstrom)
                for index in range(len(angstrom)):
                    forward = angstrom.copy()
                    backward = angstrom.copy()
                    forward[index] += self.fd_step
                    backward[index] -= self.fd_step
                    e_plus = cli_neb_energy(self.cli_cmd, self.smiles, forward.tolist(), self.engine_method)
                    e_minus = cli_neb_energy(self.cli_cmd, self.smiles, backward.tolist(), self.engine_method)
                    gradient[index] = (e_plus - e_minus) / (2.0 * self.fd_step)
            return {
                "energy": energy_kcal * KCAL_MOL_TO_HARTREE,
                "gradient": gradient * KCAL_MOL_ANGSTROM_TO_HARTREE_BOHR,
            }

    def xyz_text() -> str:
        lines: List[str] = []
        for image_index, image in enumerate(
            interpolate_images(prepared["start_coords"], prepared["end_coords"], prepared["n_images"])
        ):
            lines.append(str(len(prepared["atom_symbols"])))
            lines.append(f"{prepared['name']} image {image_index}")
            for atom_index, symbol in enumerate(prepared["atom_symbols"]):
                base = atom_index * 3
                lines.append(
                    f"{symbol} {image[base]:.10f} {image[base + 1]:.10f} {image[base + 2]:.10f}"
                )
        return "\n".join(lines) + "\n"

    with tempfile.TemporaryDirectory(prefix=f"geom_neb_{prepared['name']}_") as tmpdir:
        xyz_path = Path(tmpdir) / "chain.xyz"
        xyz_path.write_text(xyz_text(), encoding="utf-8")

        try:
            geometric_neb.logger.setLevel(logging.CRITICAL)
            geometric_neb.logger.disabled = True
            molecule = geometric.molecule.Molecule(str(xyz_path))
            engine = SciFormCliEngine(molecule, prepared["smiles"], direct_cli, engine_method=method)
            params = geometric_neb.NEBParams(
                images=prepared["n_images"],
                neb_maxcyc=2,
                maxg=0.50,
                avgg=0.25,
                nebk=0.05,
                trust=0.05,
                tmax=0.10,
                climb=0.0,
                align=True,
                verbose=False,
            )
            chain = geometric_neb.ElasticBand(molecule, engine, tmpdir, params)
            optimized_chain, opt_cycles = geometric_neb.OptimizeChain(chain, engine, params)
        except Exception as exc:
            return {"success": False, "stderr": f"geomeTRIC NEB error: {exc}"}

        images = []
        for image_index, structure in enumerate(optimized_chain.Structures):
            coords = np.asarray(structure.M.xyzs[0], dtype=float).reshape(-1).tolist()
            energy_kcal_mol = float(structure.energy * HARTREE_TO_KCAL_MOL)
            images.append(
                {
                    "index": image_index,
                    "coords": coords,
                    "potential_energy_kcal_mol": energy_kcal_mol,
                }
            )

        return {
            "success": True,
            "payload": {
                "images": images,
                "opt_cycles": int(opt_cycles),
                "converged": bool(
                    optimized_chain.maxg <= params.maxg and optimized_chain.avgg <= params.avgg
                ),
                "maxg_ev_ang": float(optimized_chain.maxg),
                "avgg_ev_ang": float(optimized_chain.avgg),
            },
            "summary": summarize_neb_images(images),
        }


def tool_availability() -> Dict[str, bool]:
    availability = {"rdkit": False, "geometric": False, "sci_form": False, "tblite": False, "ase": False, "pyscf": False}
    for name, mod in [("rdkit", "rdkit"), ("geometric", "geometric"), ("sci_form", "sci_form"), ("pyscf", "pyscf")]:
        try:
            __import__(mod)
            availability[name] = True
        except Exception:
            pass
    try:
        from ase.mep.neb import NEB  # noqa: F401
        availability["ase"] = True
    except Exception:
        pass
    try:
        from tblite.ase import TBLite  # noqa: F401
        availability["tblite"] = True
    except Exception:
        pass
    return availability


# ── External engine adapters ─────────────────────────────────────────────────

EV_TO_KCAL_MOL = 23.0605
ATOMIC_SYMBOLS = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br",
]


def _elements_from_case(case: "PathBenchmarkCase", label: str) -> Optional[List[int]]:
    """Extract atomic numbers from a prepared case via sci_form.embed."""
    try:
        import sci_form
        seed = case.seed_a if label == "reactant" else case.seed_b
        conf = sci_form.embed(case.smiles, seed)
        return list(conf.elements) if conf.is_ok() else None
    except Exception:
        return None


def _symbols_to_numbers(symbols: List[str]) -> List[int]:
    """Convert atomic symbols to numbers using local table."""
    lut = {s: i for i, s in enumerate(ATOMIC_SYMBOLS)}
    return [lut.get(s, 0) for s in symbols]


def tblite_single_point(elements: List[int], coords: List[float], method: str = "GFN2-xTB") -> Dict[str, Any]:
    """Compute tblite xTB single-point energy for cross-validation."""
    try:
        import numpy as np
        from ase import Atoms
        from tblite.ase import TBLite
    except Exception as exc:
        return {"success": False, "stderr": f"tblite/ASE unavailable: {exc}"}

    positions = np.array(coords).reshape(-1, 3)
    atoms = Atoms(numbers=elements, positions=positions)
    atoms.calc = TBLite(method=method)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            energy_ev = atoms.get_potential_energy()
            forces = atoms.get_forces()
    except Exception as exc:
        return {"success": False, "stderr": f"tblite calculation error: {exc}"}
    return {
        "success": True,
        "energy_ev": float(energy_ev),
        "energy_kcal_mol": float(energy_ev * EV_TO_KCAL_MOL),
        "max_force_ev_ang": float(np.max(np.abs(forces))),
        "method": method,
    }


def pyscf_single_point(elements: List[int], coords: List[float], basis: str = "sto-3g") -> Dict[str, Any]:
    """Compute PySCF RHF single-point energy for ab-initio reference."""
    try:
        import numpy as np
        from pyscf import gto, scf
    except Exception as exc:
        return {"success": False, "stderr": f"PySCF unavailable: {exc}"}

    positions = np.array(coords).reshape(-1, 3)
    atom_str = "; ".join(
        f"{ATOMIC_SYMBOLS[z] if z < len(ATOMIC_SYMBOLS) else 'X'} {positions[i, 0]:.10f} {positions[i, 1]:.10f} {positions[i, 2]:.10f}"
        for i, z in enumerate(elements)
    )
    try:
        mol = gto.Mole()
        mol.atom = atom_str
        mol.basis = basis
        mol.unit = "Angstrom"
        mol.verbose = 0
        mol.build()
        mf = scf.RHF(mol)
        mf.verbose = 0
        with contextlib.redirect_stdout(io.StringIO()):
            energy_hartree = mf.kernel()
    except Exception as exc:
        return {"success": False, "stderr": f"PySCF SCF error: {exc}"}
    return {
        "success": True,
        "energy_hartree": float(energy_hartree),
        "energy_kcal_mol": float(energy_hartree * HARTREE_TO_KCAL_MOL),
        "converged": bool(mf.converged),
        "basis": basis,
    }


def run_ase_neb_tblite(prepared: Dict[str, Any], gfn_method: str = "GFN2-xTB", max_steps: int = 10) -> Dict[str, Any]:
    """Run ASE NEB with tblite calculator — truly external energy surface."""
    try:
        import numpy as np
        from ase import Atoms
        from ase.mep.neb import NEB
        from ase.optimize import BFGS
        from tblite.ase import TBLite
    except Exception as exc:
        return {"success": False, "stderr": f"ASE/tblite unavailable: {exc}"}

    atom_symbols = prepared["atom_symbols"]
    n_images = prepared["n_images"]
    start_pos = np.array(prepared["start_coords"]).reshape(-1, 3)
    end_pos = np.array(prepared["end_coords"]).reshape(-1, 3)

    start_atoms = Atoms(symbols=atom_symbols, positions=start_pos)
    start_atoms.calc = TBLite(method=gfn_method)
    end_atoms = Atoms(symbols=atom_symbols, positions=end_pos)
    end_atoms.calc = TBLite(method=gfn_method)

    images = [start_atoms]
    for _ in range(n_images - 2):
        image = start_atoms.copy()
        image.calc = TBLite(method=gfn_method)
        images.append(image)
    images.append(end_atoms)

    neb = NEB(images, method="improvedtangent")
    neb.interpolate()

    try:
        optimizer = BFGS(neb, logfile=None)
        with contextlib.redirect_stdout(io.StringIO()):
            optimizer.run(fmax=0.5, steps=max_steps)
    except Exception as exc:
        return {"success": False, "stderr": f"ASE NEB optimization error: {exc}"}

    result_images = []
    for i, image in enumerate(images):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                energy_ev = image.get_potential_energy()
        except Exception:
            energy_ev = 0.0
        result_images.append({
            "index": i,
            "coords": image.get_positions().reshape(-1).tolist(),
            "potential_energy_kcal_mol": float(energy_ev * EV_TO_KCAL_MOL),
        })

    return {
        "success": True,
        "payload": {"images": result_images, "method": gfn_method},
        "summary": summarize_neb_images(result_images),
    }


def run_cross_validation_layer(cli_prefix: List[str], methods: List[str]) -> Dict[str, Any]:
    """Cross-validate sci-form internal energies against external libraries.

    Always tests GFN1/GFN2/xTB single-point (fast) plus whatever methods are passed.
    The key comparison is sci-form GFN1 vs tblite GFN1 and sci-form GFN2 vs tblite GFN2.
    """
    # Always include key xTB methods in cross-validation (single-point is fast)
    internal_methods = list(dict.fromkeys(methods + ["xtb", "gfn1", "gfn2"]))
    rows: List[Dict[str, Any]] = []
    for case in PATH_TS_CASES:
        prepared = prepare_path_case(case)
        if not prepared.get("success"):
            continue
        elements = _symbols_to_numbers(case.atom_symbols)
        for label, coords in [("reactant", prepared["start_coords"]), ("product", prepared["end_coords"])]:
            row: Dict[str, Any] = {"case": case.name, "geometry": label}
            # Internal methods via neb-energy (single-point, fast for all methods)
            for method in internal_methods:
                try:
                    e = cli_neb_energy(cli_prefix, prepared["smiles"], coords, method)
                    row[f"sci_form_{method}_kcal"] = e
                except Exception as exc:
                    row[f"sci_form_{method}_error"] = str(exc)
            # tblite GFN1
            tb1 = tblite_single_point(elements, coords, "GFN1-xTB")
            if tb1.get("success"):
                row["tblite_gfn1_kcal"] = tb1["energy_kcal_mol"]
            else:
                row["tblite_gfn1_error"] = tb1.get("stderr", "unknown")
            # tblite GFN2
            tb2 = tblite_single_point(elements, coords, "GFN2-xTB")
            if tb2.get("success"):
                row["tblite_gfn2_kcal"] = tb2["energy_kcal_mol"]
            else:
                row["tblite_gfn2_error"] = tb2.get("stderr", "unknown")
            # PySCF HF/STO-3G
            pyscf_res = pyscf_single_point(elements, coords)
            if pyscf_res.get("success"):
                row["pyscf_hf_kcal"] = pyscf_res["energy_kcal_mol"]
                row["pyscf_hf_converged"] = pyscf_res["converged"]
            else:
                row["pyscf_hf_error"] = pyscf_res.get("stderr", "unknown")
            rows.append(row)
    return {
        "rows": rows,
        "note": (
            "Cross-validation compares sci-form internal energies against tblite GFN1/GFN2-xTB "
            "and PySCF HF/STO-3G on identical geometries. "
            "sci-form GFN1 vs tblite GFN1 and sci-form GFN2 vs tblite GFN2 should agree closely."
        ),
    }


def run_semantics_layer() -> Dict[str, Any]:
    with contextlib.redirect_stdout(io.StringIO()):
        results = compare_reactions()
    comparable = [
        row
        for row in results
        if row["sci_form"].get("success") is not None and row["rdkit"].get("success") is not None
    ]
    agreements = sum(
        1
        for row in comparable
        if row["sci_form"].get("success") == row["rdkit"].get("success")
    )
    return {
        "total_cases": len(results),
        "comparable_cases": len(comparable),
        "agreements": agreements,
        "disagreements": len(comparable) - agreements,
        "results": results,
    }


def run_path_layer(cli_prefix: List[str], methods: List[str], run_external: bool = False) -> Dict[str, Any]:
    cases: List[Dict[str, Any]] = []
    for case in PATH_TS_CASES:
        prepared = prepare_path_case(case)
        if not prepared.get("success"):
            cases.append({"name": case.name, "input_smiles": case.smiles, "preparation": prepared})
            continue

        plan = cli_json(cli_prefix, ["gsm-backend-plan", prepared["smiles"]])

        # Internal NEB: run each configured method
        neb_by_method: Dict[str, Any] = {}
        for method in methods:
            neb_by_method[method] = run_sci_form_neb(cli_prefix, prepared, method=method)

        # geomeTRIC NEB with configurable engine (use first method)
        geometric_neb = run_geometric_neb(cli_prefix, prepared, method=methods[0])

        # External ASE NEB with tblite (truly independent surface)
        ase_neb_result = None
        if run_external:
            ase_neb_result = run_ase_neb_tblite(prepared, gfn_method="GFN2-xTB", max_steps=10)

        entry: Dict[str, Any] = {
            "name": case.name,
            "input_smiles": case.smiles,
            "smiles": prepared["smiles"],
            "preparation": {
                "success": True,
                "endpoint_rmsd_angstrom": prepared["endpoint_rmsd_angstrom"],
                "selection_note": prepared["selection_note"],
            },
            "plan": plan,
            "sci_form_neb": neb_by_method,
            "geometric_neb": geometric_neb,
        }
        if ase_neb_result is not None:
            entry["ase_neb_tblite"] = ase_neb_result
        cases.append(entry)

    return {
        "cases": cases,
        "methods": methods,
        "note": (
            "Path smoke cases use curated sci-form embed seed pairs with rigid endpoint alignment. "
            "The internal path uses simplified NEB with configurable backends (methods: "
            + ", ".join(methods)
            + "). The geomeTRIC reference uses ElasticBand wired to the sci-form CLI engine. "
            + ("ASE NEB with tblite GFN2-xTB provides a truly external energy-surface reference. " if run_external else "")
            + "GSM+MBH+HTST remains available through the CLI."
        ),
    }


def run_common_geometry_energy_layer(cli_prefix: List[str], methods: List[str]) -> Dict[str, Any]:
    rows = []
    for case in PATH_TS_CASES:
        prepared = prepare_path_case(case)
        if not prepared.get("success"):
            rows.append(
                {
                    "case": case.name,
                    "geometry": "endpoints",
                    "comparison": prepared,
                }
            )
            continue
        for label, coords in (("reactant", prepared["start_coords"]), ("product", prepared["end_coords"])):
            # Use gsm-compare-backends with the chosen methods
            rows.append(
                {
                    "case": case.name,
                    "geometry": label,
                    "smiles": prepared["smiles"],
                    "comparison": cli_json(
                        cli_prefix,
                        [
                            "gsm-compare-backends",
                            prepared["smiles"],
                            json.dumps(coords),
                            "--methods",
                            json.dumps(methods),
                        ],
                    ),
                }
            )
    return {
        "rows": rows,
        "methods": methods,
        "note": f"Energy comparisons reuse the same curated, aligned endpoint geometries. Methods: {', '.join(methods)}. Heavier methods available through gsm-compare-backends on demand.",
    }


def render_text(summary: Dict[str, Any]) -> None:
    print("=" * 88)
    print("Layered Reaction Comparison")
    print("=" * 88)
    print("Tool availability:")
    for name, available in summary["tool_availability"].items():
        print(f"  {name:12s}: {'yes' if available else 'no'}")

    semantics = summary["semantics"]
    print("\nSemantics layer:")
    print(
        f"  comparable={semantics['comparable_cases']} agreements={semantics['agreements']} disagreements={semantics['disagreements']}"
    )

    path_layer = summary["path_ts"]
    methods_used = path_layer.get("methods", ["uff"])
    print(f"\nPath / TS layer (methods: {', '.join(methods_used)}):")
    for case in path_layer["cases"]:
        if not case["preparation"].get("success"):
            print(f"  {case['name']}: endpoint preparation error")
            continue
        rmsd = case["preparation"]["endpoint_rmsd_angstrom"]
        print(f"  {case['name']}: endpoint RMSD={rmsd:.3f} A")

        # sci-form NEB per method
        neb_data = case.get("sci_form_neb", {})
        for method, neb_result in (neb_data.items() if isinstance(neb_data, dict) else []):
            state = "ok" if neb_result.get("success") else "error"
            if neb_result.get("summary"):
                barrier = neb_result["summary"]["barrier_kcal_mol"]
                print(f"    sci-form NEB ({method}): {state}, barrier={barrier:.3f} kcal/mol")
            else:
                print(f"    sci-form NEB ({method}): {state}")

        # geomeTRIC
        geometric_state = "ok" if case["geometric_neb"].get("success") else "error"
        if case["geometric_neb"].get("summary"):
            barrier = case["geometric_neb"]["summary"]["barrier_kcal_mol"]
            print(f"    geomeTRIC NEB: {geometric_state}, barrier={barrier:.3f} kcal/mol")
        else:
            print(f"    geomeTRIC NEB: {geometric_state}")

        # ASE NEB with tblite (if run)
        ase_neb = case.get("ase_neb_tblite")
        if ase_neb is not None:
            ase_state = "ok" if ase_neb.get("success") else "error"
            if ase_neb.get("summary"):
                barrier = ase_neb["summary"]["barrier_kcal_mol"]
                gfn = ase_neb.get("payload", {}).get("method", "GFN2-xTB")
                print(f"    ASE NEB ({gfn}): {ase_state}, barrier={barrier:.3f} kcal/mol")
            else:
                stderr = ase_neb.get("stderr", "")
                print(f"    ASE NEB (tblite): {ase_state} {stderr}")
    print(f"  note: {path_layer['note']}")

    energy_layer = summary["common_geometry_energies"]
    energy_methods = energy_layer.get("methods", [])
    print(f"\nCommon-geometry energy layer (methods: {', '.join(energy_methods)}):")
    for row in energy_layer["rows"]:
        state = "ok" if row["comparison"].get("success") else "error"
        print(f"  {row['case']} / {row['geometry']}: {state}")
    print(f"  note: {energy_layer['note']}")

    # Cross-validation layer
    cross_val = summary.get("cross_validation")
    if cross_val:
        print("\nCross-validation layer (sci-form vs tblite/PySCF):")
        for row in cross_val.get("rows", []):
            parts = [f"{row['case']}/{row['geometry']}"]
            for key, val in sorted(row.items()):
                if key in ("case", "geometry"):
                    continue
                if key.endswith("_kcal") and isinstance(val, (int, float)):
                    parts.append(f"{key}={val:.2f}")
                elif key.endswith("_error"):
                    parts.append(f"{key}={val}")
            print(f"  {' | '.join(parts)}")
        print(f"  note: {cross_val['note']}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Layered reaction comparison harness")
    parser.add_argument("--cli", help="CLI launcher, e.g. './target/release/sci-form'")
    parser.add_argument("--json", action="store_true", help="Emit only JSON summary")
    parser.add_argument(
        "--methods",
        default="uff,pm3",
        help="Comma-separated list of internal backends: uff,mmff94,pm3,xtb,gfn1,gfn2,hf3c (default: uff,pm3)",
    )
    parser.add_argument(
        "--external",
        action="store_true",
        help="Enable external validation (tblite ASE NEB + PySCF single-point cross-validation)",
    )
    args = parser.parse_args()

    cli_prefix = detect_cli(args.cli)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    summary: Dict[str, Any] = {
        "tool_availability": tool_availability(),
        "methods": methods,
        "external": args.external,
        "semantics": run_semantics_layer(),
        "path_ts": run_path_layer(cli_prefix, methods, run_external=args.external),
        "common_geometry_energies": run_common_geometry_energy_layer(cli_prefix, methods),
    }
    if args.external:
        summary["cross_validation"] = run_cross_validation_layer(cli_prefix, methods)

    if args.json:
        print(json.dumps(summary, indent=2))
    else:
        render_text(summary)
        print()
        generate_report(summary["semantics"]["results"])
    return 0


if __name__ == "__main__":
    sys.exit(main())