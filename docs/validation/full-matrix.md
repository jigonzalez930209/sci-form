# Full validation matrix

This appendix is generated from `reports/validation/summary.json`. It is the
raw record behind the curated [scientific validation report](./scientific-validation.md).

- Rows: **129**
- Passing gates: **122**
- Recorded method/provenance gaps: **7**
- Failing gates: **0**
- Reference tools: RDKit 2026.09.1, PySCF 2.14.0, tblite 0.7.0,
  dftd4 4.3.0, xtb 6.7.1, LAMMPS, ASE 3.29.0, pymatgen 2026.9.24,
  pormake 0.2.3, FreeSASA, NumPy/SciPy, and the independent oracles.

`GAP` means the row documents a method or provenance gap; it is not a
silent failure. The per-row notes and the owning test are kept in
`reports/validation/summary.json`; the curated report explains the
limitations behind the gaps.

| ID | Observable | Reference | Metric | Measured | Tolerance | Level | Status |
|---|---|---|---|---|---|---|---|
| F01 | <code>parse_<wbr>invariants_<wbr>mismatches</code> | RDKit 2026.09.1 | <code>mismatched_<wbr>cases</code> | 0 | 0 | L2 | PASS |
| F05 | <code>sssr_<wbr>multiset_<wbr>mismatches</code> | RDKit 2026.09.1 | <code>mismatched_<wbr>cases</code> | 0 | 0 | L2 | PASS |
| F13 | <code>kabsch_<wbr>rmsd_<wbr>max_<wbr>abs</code> | NumPy SVD Kabsch | <code>max_<wbr>abs</code> | 2.220e-15 | 1.000e-06 | L2 | PASS |
| F14 | <code>butina_<wbr>partition_<wbr>mismatches</code> | RDKit 2026.09.1 Butina | <code>mismatched_<wbr>cases</code> | 0 | 0 | L2 | PASS |
| F16 | <code>uff_<wbr>energy_<wbr>abs_<wbr>error</code> | RDKit 2026.09.1 UFFGetMoleculeForceField | <code>max_<wbr>abs</code> | 2.274e-13 | 0.01 | L2 | PASS |
| F17 | <code>uff_<wbr>atom_<wbr>type_<wbr>mismatches</code> | RDKit 2026.09.1 UFF AtomTyper | <code>mismatched_<wbr>labels</code> | 0 | 0 | L2 | PASS |
| F21 | <code>gasteiger_<wbr>charges_<wbr>max_<wbr>abs</code> | RDKit 2026.09.1 | <code>max_<wbr>abs</code> | 8.327e-17 | 0.02 | L2 | PASS |
| F22 | <code>eeqbc_<wbr>charges_<wbr>vs_<wbr>multicharge</code> | — | <code>max_<wbr>abs</code> | 6.902e-12 | 1.000e-10 | L2 | PASS |
| F22 | <code>eeqbc_<wbr>cn_<wbr>vs_<wbr>multicharge</code> | — | <code>max_<wbr>abs</code> | 1.908e-11 | 1.000e-10 | L2 | PASS |
| F23 | <code>mulliken_<wbr>charges_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 8.564e-08 | 0.005 | L2 | PASS |
| F26 | <code>scf_<wbr>mulliken_<wbr>dipole_<wbr>parity_<wbr>debye</code> | — | <code>max_<wbr>abs</code> | 1.829e-07 | 1.000e-06 | L2 | PASS |
| F28 | <code>sasa_<wbr>max_<wbr>rel_<wbr>error</code> | FreeSASA python-bindings calcCoord (same radii/probe/npoints) | <code>max_<wbr>rel</code> | 0.0153304 | 0.02 | L2 | PASS |
| F31 | <code>eht_<wbr>metal_<wbr>params_<wbr>vs_<wbr>yaehmop</code> | — | <code>max_<wbr>abs</code> | 0 | 0.02 | L2 | PASS |
| F36 | <code>scf_<wbr>matrices_<wbr>energy_<wbr>mulliken</code> | PySCF 2.14.0 | <code>max_<wbr>abs_<wbr>e_<wbr>tot</code> | 1.224e-08 | 1.000e-07 | L2 | PASS |
| F37 | <code>two_<wbr>electron_<wbr>integrals_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 2.139e-11 | 1.000e-10 | L2 | PASS |
| F38 | <code>uhf_<wbr>rohf_<wbr>energy_<wbr>s2</code> | PySCF 2.14.0 | <code>max_<wbr>abs_<wbr>e_<wbr>tot</code> | 6.502e-08 | 1.000e-06 | L2 | PASS |
| F39 | <code>ao_<wbr>to_<wbr>mo_<wbr>transform_<wbr>max_<wbr>abs</code> | PySCF 2.14.0 ao2mo.full | <code>max_<wbr>abs</code> | 1.834e-14 | 1.000e-10 | L2 | PASS |
| F41 | <code>ao_<wbr>dipole_<wbr>integrals_<wbr>max_<wbr>abs_<wbr>au</code> | — | <code>max_<wbr>abs</code> | 1.100e-09 | 1.000e-08 | L2 | PASS |
| F41 | <code>cis_<wbr>oscillator_<wbr>strength_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.004068 | 0.005 | L2 | PASS |
| F41 | <code>cis_<wbr>singlet_<wbr>excitations_<wbr>ev</code> | PySCF 2.14.0 tdscf.TDA (singlets) | <code>max_<wbr>abs</code> | 1.724e-06 | 0.01 | L2 | PASS |
| F41 | <code>cis_<wbr>singlet_<wbr>excitations_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 2.622e-06 | 0.01 | L2 | PASS |
| F41 | <code>cisd_<wbr>correlation_<wbr>vs_<wbr>pyscf</code> | — | <code>max_<wbr>abs</code> | 2.073e-08 | 1.000e-06 | L2 | PASS |
| F42 | <code>mopac_<wbr>pm3_<wbr>max_<wbr>abs</code> | MOPAC 23.2.5 PM3 1SCF | <code>max_<wbr>abs</code> | 0.00834273 | 0.02 | L2 | PASS |
| F42 | <code>pm3_<wbr>heat_<wbr>of_<wbr>formation_<wbr>guard</code> | NIST WebBook gas-phase ΔHf°(298 K) | <code>max_<wbr>abs_<wbr>vs_<wbr>experiment</code> | 8.60829 | 5000 | L4 | PASS |
| F43 | <code>gfn0_<wbr>single_<wbr>point_<wbr>ev</code> | xtb GFN0-xTB | <code>max_<wbr>abs</code> | 1.650e-04 | 0.01 | L4 | PASS |
| F44 | <code>gfn1_<wbr>gfn2_<wbr>fixture_<wbr>provenance</code> | tblite 0.7.0 | <code>max_<wbr>abs</code> | — | 0.0217 | L4 | GAP |
| F44 | <code>gfn1_<wbr>single_<wbr>point_<wbr>ev</code> | tblite 0.7.0 | <code>max_<wbr>abs</code> | 8.063e-05 | 0.01 | L2 | PASS |
| F44-UKS | <code>uks_<wbr>gfn1_<wbr>nospin_<wbr>vs_<wbr>tblite_<wbr>uhf_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 0.0335832 | 0.05 | L2 | PASS |
| F45 | <code>gfn2_<wbr>single_<wbr>point_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 2.253e-05 | 0.0217 | L2 | PASS |
| F45-UKS | <code>uks_<wbr>nospin_<wbr>vs_<wbr>tblite_<wbr>uhf_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 0.0345862 | 0.05 | L2 | PASS |
| F46 | <code>d3_<wbr>hf3c_<wbr>max_<wbr>abs_<wbr>hartree</code> | — | <code>max_<wbr>abs</code> | 6.772e-09 | 1.594e-06 | L2 | PASS |
| F46 | <code>d4_<wbr>pbe_<wbr>standalone_<wbr>max_<wbr>abs_<wbr>hartree</code> | — | <code>max_<wbr>abs</code> | 6.198e-10 | 1.594e-06 | L2 | PASS |
| F46 | <code>d4d3_<wbr>dispersion_<wbr>fixture_<wbr>provenance</code> | dftd4 4.3.0 | <code>max_<wbr>abs</code> | — | 1.594e-06 | L4 | GAP |
| F47 | <code>eas_<wbr>site_<wbr>ranking</code> | EAS textbook orientation rules (DFT-confirmed) | <code>fraction_<wbr>correct</code> | 0.875 | 0.8 | L4 | PASS |
| F48 | <code>hessian_<wbr>diagonalization_<wbr>vs_<wbr>xtb</code> | xtb GFN2-xTB vibspectrum | <code>max_<wbr>abs</code> | 0.00891019 | 2 | L2 | PASS |
| F49 | <code>ir_<wbr>group_<wbr>recall_<wbr>vs_<wbr>nist</code> | — | <code>recall</code> | 1 | 0.9 | L2 | PASS |
| F54 | <code>hla_<wbr>baseline_<wbr>vs_<wbr>reference</code> | — | <code>max_<wbr>abs</code> | 0.004 | 0.01 | L2 | PASS |
| F54 | <code>hla_<wbr>substituent_<wbr>vs_<wbr>reference</code> | — | <code>max_<wbr>abs</code> | 0.420899 | 0.5 | L2 | PASS |
| F59 | <code>ecfp_<wbr>tanimoto_<wbr>rank_<wbr>correlation</code> | RDKit 2026.09.1 Morgan/ECFP4 | <code>spearman</code> | 0.755541 | 0.7 | L2 | PASS |
| F64 | <code>space_<wbr>groups_<wbr>vs_<wbr>spglib</code> | — | <code>max_<wbr>abs_<wbr>translation</code> | 0 | 1.000e-12 | L2 | PASS |
| F65 | <code>cif_<wbr>roundtrip_<wbr>vs_<wbr>pymatgen</code> | pymatgen 2026.9.24 | <code>max_<wbr>abs_<wbr>lattice_<wbr>param</code> | 0 | 1.000e-06 | L2 | PASS |
| F66 | <code>minimum_<wbr>image_<wbr>distance_<wbr>vs_<wbr>ase</code> | ASE 3.29.0 | <code>max_<wbr>abs</code> | 8.882e-16 | 1.000e-10 | L0/L2 | PASS |
| F67 | <code>mof5_<wbr>assembly_<wbr>vs_<wbr>pormake</code> | pormake Builder (pcu/N16/E14) | <code>max_<wbr>abs</code> | 0 | 0 | L2 | PASS |
| F67 | <code>topology_<wbr>degree_<wbr>mismatches</code> | — | <code>mismatches</code> | 0 | 0 | L2 | PASS |
| F68 | <code>mof5_<wbr>porosity_<wbr>lcd_<wbr>vs_<wbr>pyzeo</code> | — | <code>max_<wbr>abs</code> | 0.00372613 | 0.3 | L2 | PASS |
| F68 | <code>mof5_<wbr>porosity_<wbr>pld_<wbr>vs_<wbr>pyzeo</code> | — | <code>max_<wbr>abs</code> | 0.290373 | 0.5 | L2 | PASS |
| F68 | <code>mof5_<wbr>porosity_<wbr>vf_<wbr>vs_<wbr>pyzeo</code> | — | <code>max_<wbr>abs</code> | 7.769e-04 | 0.01 | L2 | PASS |
| F69 | <code>powder_<wbr>xrd_<wbr>vs_<wbr>pymatgen</code> | pymatgen 2026.9.24 | <code>max_<wbr>abs_<wbr>two_<wbr>theta</code> | 4.263e-14 | 0.01 | L2 | PASS |
| F70 | <code>uff_<wbr>opt_<wbr>energy_<wbr>vs_<wbr>rdkit</code> | RDKit 2025.09 UFFOptimizeMolecule | <code>max_<wbr>abs</code> | 7.844e-10 | 1.000e-06 | L2 | PASS |
| F71 | <code>bounds_<wbr>matrix_<wbr>max_<wbr>abs</code> | RDKit 2026.09.1 GetMoleculeBoundsMatrix | <code>max_<wbr>abs</code> | — | 0.35 | L2 | GAP |
| F71 | <code>bounds_<wbr>topo_<wbr>1_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.0258985 | 0.03 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>2_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.33652 | 0.35 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>3_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.373165 | 0.4 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>4_<wbr>divergent_<wbr>pairs</code> | — | <code>count</code> | 8 | 388 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>4_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 1.776e-15 | 1 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>5_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.5375 | 0.6 | L2 | PASS |
| F71 | <code>bounds_<wbr>topo_<wbr>6_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.95 | 1 | L2 | PASS |
| F72 | <code>embedding_<wbr>bond_<wbr>length_<wbr>max_<wbr>abs</code> | RDKit 2026.09.1 AllChem.EmbedMolecule (ETKDGv3) | <code>max_<wbr>abs</code> | 0.505243 | 1 | L2 | PASS |
| F72 | <code>embedding_<wbr>heavy_<wbr>rmsd_<wbr>mean</code> | — | <code>mean</code> | 0.0285526 | 3 | L2 | PASS |
| F73 | <code>etkdg3d_<wbr>ff_<wbr>energy_<wbr>abs_<wbr>error</code> | RDKit 2026.09.1 ForceField (ETKDG-3D term set) | <code>max_<wbr>abs</code> | 9.095e-13 | 5 | L2 | PASS |
| F73 | <code>etkdg3d_<wbr>ff_<wbr>inversion_<wbr>energy_<wbr>abs_<wbr>error</code> | — | <code>max_<wbr>abs</code> | 4.263e-14 | 5 | L2 | PASS |
| F74 | <code>bfgs_<wbr>bounds_<wbr>energy_<wbr>after_<wbr>min</code> | — | <code>max_<wbr>abs</code> | 0.00569787 | 0.01 | L2 | PASS |
| F74 | <code>bfgs_<wbr>minimised_<wbr>rmsd</code> | RDKit 2026.09.1 ForceField::minimize | <code>max_<wbr>abs</code> | 0.135508 | 0.25 | L2 | PASS |
| F75 | <code>smarts_<wbr>match_<wbr>agreement</code> | RDKit 2026.09.1 GetSubstructMatches | <code>fraction_<wbr>geq</code> | 1 | 0.9 | L2 | PASS |
| F75 | <code>smarts_<wbr>raw_<wbr>count_<wbr>agreement</code> | — | <code>fraction_<wbr>geq</code> | 0.999403 | 0.9 | L2 | PASS |
| F76 | <code>csd_<wbr>torsion_<wbr>table_<wbr>drift</code> | RDKit Release_2026_09_1 torsionPreferences*.in | <code>differing_<wbr>entries</code> | 0 | 0 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>impl_<wbr>parity_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 1.814e-06 | 0.001 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>impl_<wbr>parity_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 5.634e-08 | 1.000e-05 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>term_<wbr>effect_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 18.2198 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>fullspace_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 9.73173 | 10 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>fullspace_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 1.0161 | 5 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 9.73173 | 10 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 1.0161 | 5 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tdhf_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 9.25628 | 10 | L2 | PASS |
| F77 | <code>stda_<wbr>exchange_<wbr>vs_<wbr>pyscf_<wbr>tdhf_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 0.806881 | 5 | L2 | PASS |
| F77 | <code>stda_<wbr>impl_<wbr>parity_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 3.264e-06 | 0.001 | L2 | PASS |
| F77 | <code>stda_<wbr>impl_<wbr>parity_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 4.876e-08 | 1.000e-05 | L2 | PASS |
| F77 | <code>stda_<wbr>missing_<wbr>exchange_<wbr>gap_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 18.2198 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>monopole_<wbr>coulomb_<wbr>gap_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 9.73173 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>scf_<wbr>orbital_<wbr>energy_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 1.253e-07 | 1.000e-06 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>fullspace_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 16.0437 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>fullspace_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 0.927451 | 5 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>max_<wbr>abs_<wbr>ev</code> | PySCF 2.14.0 tdscf.TDA / tdscf.TDHF | <code>max_<wbr>abs</code> | 16.0437 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tda_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 0.927451 | 5 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tdhf_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 16.071 | 20 | L2 | PASS |
| F77 | <code>stda_<wbr>vs_<wbr>pyscf_<wbr>tdhf_<wbr>max_<wbr>abs_<wbr>oscillator</code> | — | <code>max_<wbr>abs</code> | 0.789967 | 5 | L2 | PASS |
| F78 | <code>esp_<wbr>far_<wbr>field_<wbr>max_<wbr>abs_<wbr>dev</code> | PySCF 2.14.0 RHF/STO-3G density ESP (cubegen scheme) | <code>max_<wbr>abs</code> | 0.00439039 | 0.006 | L2 | PASS |
| F78 | <code>esp_<wbr>far_<wbr>field_<wbr>max_<wbr>rel_<wbr>dev</code> | — | <code>max_<wbr>abs</code> | 0.0320234 | 0.05 | L2 | PASS |
| F78 | <code>esp_<wbr>grid_<wbr>formula_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 1.110e-16 | 1.000e-12 | L2 | PASS |
| F78 | <code>esp_<wbr>monopole_<wbr>residual_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 0.00923947 | 0.02 | L2 | PASS |
| F78 | <code>esp_<wbr>near_<wbr>field_<wbr>max_<wbr>abs_<wbr>dev</code> | — | <code>max_<wbr>abs</code> | 0.0497805 | 1 | L2 | PASS |
| F79b | <code>ir_<wbr>frequencies_<wbr>max_<wbr>abs_<wbr>cm1</code> | — | <code>max_<wbr>abs</code> | 0.44626 | 5 | L2 | PASS |
| F79b | <code>ir_<wbr>intensity_<wbr>vs_<wbr>exact_<wbr>dipole_<wbr>reference_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 37.1477 | 60 | L2 | PASS |
| F79b | <code>ir_<wbr>intensity_<wbr>vs_<wbr>mulliken_<wbr>reference_<wbr>max_<wbr>abs</code> | PySCF 2.14.0 analytic Hessian + SCF dipole derivatives | <code>max_<wbr>abs</code> | 0.576716 | 5 | L2 | PASS |
| F80 | <code>ani_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 8.318e-07 | 1.000e-06 | L2 | PASS |
| F80 | <code>eht_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 3.055e-08 | 1.493e-04 | L2 | PASS |
| F80 | <code>mmff94_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 8.385e-08 | 4.845e-05 | L2 | PASS |
| F80 | <code>mmff94_<wbr>term_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 8.328e-08 | 1.000e-04 | L2 | PASS |
| F80 | <code>scf_<wbr>nuclear_<wbr>repulsion_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 3.885e-10 | 9.485e-07 | L2 | PASS |
| F80 | <code>uff_<wbr>gradient_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 4.008e-08 | 5.414e-05 | L2 | PASS |
| F81 | <code>cip_<wbr>stereo_<wbr>labels_<wbr>vs_<wbr>rdkit</code> | RDKit 2026.09.1 CIP labeller | <code>mismatched_<wbr>labels</code> | 0 | 0 | L2 | PASS |
| F81 | <code>smiles_<wbr>chiral_<wbr>tag_<wbr>handling</code> | — | <code>mismatched_<wbr>cases</code> | 0 | 0 | L2 | PASS |
| F82 | <code>ertl_<wbr>tpsa_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 5.684e-14 | 0.01 | L2 | PASS |
| F82 | <code>hba_<wbr>definition_<wbr>difference</code> | — | <code>max_<wbr>abs</code> | 27 | 27 | L2 | PASS |
| F82 | <code>ml_<wbr>descriptors_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 5.230e-12 | 0.01 | L2 | PASS |
| F82 | <code>ml_<wbr>descriptors_<wbr>vs_<wbr>rdkit</code> | RDKit 2026.09.1 Descriptors / rdMolDescriptors | <code>mismatched_<wbr>cases</code> | — | 0 | L2 | GAP |
| F83 | <code>canonical_<wbr>smiles_<wbr>identical_<wbr>to_<wbr>rdkit</code> | — | <code>fraction_<wbr>identical</code> | 0.25 | 0.25 | L2 | PASS |
| F83 | <code>canonical_<wbr>smiles_<wbr>roundtrip_<wbr>reordering</code> | — | <code>mismatched_<wbr>cases</code> | 0 | 0 | L2 | PASS |
| F83 | <code>canonical_<wbr>smiles_<wbr>vs_<wbr>rdkit</code> | RDKit 2026.09.1 MolToSmiles(canonical=True) | <code>mismatched_<wbr>cases</code> | — | 0 | L2 | GAP |
| F84 | <code>alpb_<wbr>electrostatic_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 3.14229 | 5 | L2 | PASS |
| F84 | <code>alpb_<wbr>solvation_<wbr>vs_<wbr>tblite</code> | tblite 0.7.0 GB/ALPB solvation | <code>max_<wbr>abs</code> | — | 0 | L2 | GAP |
| F84 | <code>alpb_<wbr>total_<wbr>vs_<wbr>tblite_<wbr>alpb_<wbr>reported</code> | — | <code>max_<wbr>abs</code> | 11.6572 | 11.6572 | L2 | PASS |
| F85 | <code>variable_<wbr>cell_<wbr>relaxed_<wbr>lattice_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 3.779e-04 | 0.001 | L2 | PASS |
| F85 | <code>variable_<wbr>cell_<wbr>stress_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 4.291e-07 | 1.000e-05 | L2 | PASS |
| F85 | <code>variable_<wbr>cell_<wbr>stress_<wbr>vs_<wbr>ase</code> | ASE 3.29.0 LennardJones stress + UnitCellFilter | <code>max_<wbr>abs</code> | — | 0 | L2 | GAP |
| F86-EHT | <code>eht_<wbr>derived_<wbr>vsip_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 4.441e-16 | 1.000e-12 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>derived_<wbr>zeta_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 5.551e-17 | 1.000e-12 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>gap_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 3.213e-12 | 1.000e-06 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>hamiltonian_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 1.420e-11 | 1.000e-10 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>homo_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 8.562e-13 | 1.000e-06 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>lumo_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 3.299e-12 | 1.000e-06 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>orbital_<wbr>energy_<wbr>max_<wbr>abs_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 1.252e-11 | 1.000e-06 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>overlap_<wbr>diagonal_<wbr>convention_<wbr>effect_<wbr>ev</code> | — | <code>max_<wbr>abs</code> | 4.097 | 1 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>overlap_<wbr>diagonal_<wbr>normalisation_<wbr>error</code> | — | <code>max_<wbr>abs</code> | 0.333333 | 0.4 | L2 | PASS |
| F86-EHT | <code>eht_<wbr>overlap_<wbr>max_<wbr>abs</code> | oracle_eht.py (independent numpy EHT) | <code>max_<wbr>abs</code> | 7.535e-13 | 1.000e-10 | L2 | PASS |
| F87-EEQ2019 | <code>eeq2019_<wbr>charge_<wbr>max_<wbr>abs</code> | oracle_eeq2019.py (independent numpy EEQ2019) | <code>max_<wbr>abs</code> | 6.082e-07 | 1.000e-05 | L2 | PASS |
| F87-EEQ2019 | <code>eeq2019_<wbr>total_<wbr>charge_<wbr>max_<wbr>abs</code> | — | <code>max_<wbr>abs</code> | 8.327e-17 | 1.000e-09 | L2 | PASS |
| F100 | <code>mmff94_<wbr>full_<wbr>energy_<wbr>max_<wbr>abs</code> | RDKit 2026.09.1 MMFF94 (MMFFGetMoleculeForceField) | <code>max_<wbr>abs</code> | 2.806e-11 | 1.000e-06 | L2 | PASS |
| F100 | <code>mmff94_<wbr>partial_<wbr>charge_<wbr>max_<wbr>abs</code> | RDKit 2026.09.1 MMFF94 (MMFFGetMoleculeForceField) | <code>max_<wbr>abs</code> | 0 | 1.000e-09 | L2 | PASS |

## Reproduce

```bash
# Rebuild the machine-readable summary from the external campaign
python3 scripts/validation/collect.py

# Re-render this appendix
python3 scripts/docs/generate_validation_matrix.py
```
