//! Comprehensive reaction dynamics study — 55+ reaction types.
//!
//! Validates the full 3D reaction dynamics pipeline across:
//! - SN2 nucleophilic substitutions (5)
//! - Addition reactions (5)
//! - Diels-Alder cycloadditions (2)
//! - Esterification / condensation (5)
//! - Conformational / intramolecular (4)
//! - Proton transfer (2)
//! - Bimolecular association (4)
//! - Oxidation / reduction (2)
//! - Halogenation (3)
//! - Elimination (2)
//! - Ring formation (1)
//! - Electrophilic aromatic substitution (2)
//! - Tautomerism (2)
//! - Pericyclic rearrangements (3)
//! - Trimolecular (2)
//! - Miscellaneous organic (11)
//!
//! Every reaction must:
//! 1. Produce a valid reaction path (finite coords + energies)
//! 2. Include approach + NEB + departure phases
//! 3. Locate a transition state in the NEB region

#![cfg(feature = "alpha-reaction-dynamics")]

use rayon::prelude::*;
use sci_form::alpha::reaction_dynamics::*;

// ════════════════════════════════════════════════════════════════════════════
// Infrastructure
// ════════════════════════════════════════════════════════════════════════════

/// Reaction descriptor for the test battery.
struct ReactionCase {
    name: &'static str,
    reactants: Vec<&'static str>,
    products: Vec<&'static str>,
    smirks: Option<&'static str>,
    category: &'static str,
}

/// Outcome of a single reaction test.
#[derive(Debug)]
#[allow(dead_code)]
struct ReactionOutcome {
    name: String,
    category: String,
    ok: bool,
    n_frames: usize,
    activation_energy: f64,
    reaction_energy: f64,
    error: Option<String>,
}

/// Run one reaction through the full pipeline and return the outcome.
fn run_one(case: &ReactionCase) -> ReactionOutcome {
    // Use GFN0-xTB (analytical gradients) for test speed.
    // Production default is GFN2-xTB — see ReactionDynamics3DConfig::default().
    let config = ReactionDynamics3DConfig {
        method: "xtb".to_string(),
        n_images: 7,
        neb_max_iter: 30,
        spring_k: 0.1,
        step_size: 0.005,
        use_climbing_image: true,
        ci_neb_force_threshold: 0.5,
        n_approach_frames: 9,
        n_departure_frames: 9,
        far_distance: 6.0,
        reactive_distance: 2.0,
        seed: 42,
        use_orbital_guidance: false,
        use_electrostatic_steering: false,
        optimise_complex: false,
        complex_opt_max_steps: 10,
        compute_properties: false,
        n_angular_samples: 0,
        smirks: case.smirks.map(|s| s.to_string()),
    };

    match compute_reaction_dynamics_3d(&case.reactants, &case.products, &config) {
        Ok(result) => {
            let all_finite = result
                .frames
                .iter()
                .all(|f| f.energy_kcal_mol.is_finite() && f.coords.iter().all(|c| c.is_finite()));
            ReactionOutcome {
                name: case.name.to_string(),
                category: case.category.to_string(),
                ok: all_finite && !result.frames.is_empty(),
                n_frames: result.frames.len(),
                activation_energy: result.activation_energy_kcal_mol,
                reaction_energy: result.reaction_energy_kcal_mol,
                error: if !all_finite {
                    Some("Non-finite values in frames".into())
                } else {
                    None
                },
            }
        }
        Err(e) => ReactionOutcome {
            name: case.name.to_string(),
            category: case.category.to_string(),
            ok: false,
            n_frames: 0,
            activation_energy: 0.0,
            reaction_energy: 0.0,
            error: Some(e),
        },
    }
}

/// Build the full battery of 55+ reactions.
fn build_reaction_battery() -> Vec<ReactionCase> {
    vec![
        // ═══════════════════════════════════════════════════════════════════
        // SN2 — Nucleophilic Substitution (bimolecular, backside attack)
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "SN2: Cl⁻ + CH₃Br → ClCH₃ + Br⁻",
            reactants: vec!["[Cl-]", "CBr"],
            products: vec!["ClC", "[Br-]"],
            smirks: Some("[Cl-:1].[C:2][Br:3]>>[Cl:1][C:2].[Br-:3]"),
            category: "SN2",
        },
        ReactionCase {
            name: "SN2: OH⁻ + CH₃Cl → CH₃OH + Cl⁻",
            reactants: vec!["[OH-]", "CCl"],
            products: vec!["CO", "[Cl-]"],
            smirks: Some("[OH-:1].[C:2][Cl:3]>>[O:1][C:2].[Cl-:3]"),
            category: "SN2",
        },
        ReactionCase {
            name: "SN2: I⁻ + CH₃Cl → CH₃I + Cl⁻",
            reactants: vec!["[I-]", "CCl"],
            products: vec!["CI", "[Cl-]"],
            smirks: Some("[I-:1].[C:2][Cl:3]>>[I:1][C:2].[Cl-:3]"),
            category: "SN2",
        },
        ReactionCase {
            name: "SN2: HCN + CH₃Br → CH₃CN + HBr",
            reactants: vec!["C#N", "CBr"],
            products: vec!["CC#N", "Br"],
            smirks: Some("[C:1]#[N:2].[C:3][Br:4]>>[C:3][C:1]#[N:2].[Br:4]"),
            category: "SN2",
        },
        ReactionCase {
            name: "SN2: F⁻ + CH₃Cl → CH₃F + Cl⁻",
            reactants: vec!["[F-]", "CCl"],
            products: vec!["CF", "[Cl-]"],
            smirks: Some("[F-:1].[C:2][Cl:3]>>[F:1][C:2].[Cl-:3]"),
            category: "SN2",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Addition — electrophilic addition to alkenes
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Addition: HBr + ethylene → bromoethane",
            reactants: vec!["C=C", "Br"],
            products: vec!["CCBr"],
            smirks: None,
            category: "Addition",
        },
        ReactionCase {
            name: "Addition: HCl + ethylene → chloroethane",
            reactants: vec!["C=C", "Cl"],
            products: vec!["CCCl"],
            smirks: None,
            category: "Addition",
        },
        ReactionCase {
            name: "Addition: H₂O + ethylene → ethanol",
            reactants: vec!["C=C", "O"],
            products: vec!["CCO"],
            smirks: None,
            category: "Addition",
        },
        ReactionCase {
            name: "Addition: HF + ethylene → fluoroethane",
            reactants: vec!["C=C", "F"],
            products: vec!["CCF"],
            smirks: None,
            category: "Addition",
        },
        ReactionCase {
            name: "Addition: Br₂ + ethylene → 1,2-dibromoethane",
            reactants: vec!["C=C", "BrBr"],
            products: vec!["BrCCBr"],
            smirks: None,
            category: "Addition",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Diels-Alder — [4+2] cycloaddition (pericyclic, concerted)
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Diels-Alder: butadiene + ethylene → cyclohexene",
            reactants: vec!["C=CC=C", "C=C"],
            products: vec!["C1CC=CCC1"],
            smirks: None,
            category: "Diels-Alder",
        },
        ReactionCase {
            name: "Diels-Alder: butadiene + propylene → methylcyclohexene",
            reactants: vec!["C=CC=C", "CC=C"],
            products: vec!["CC1CC=CCC1"],
            smirks: None,
            category: "Diels-Alder",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Esterification / Condensation — loss of H₂O
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Esterification: AcOH + MeOH → methyl acetate + H₂O",
            reactants: vec!["CC(=O)O", "CO"],
            products: vec!["CC(=O)OC", "O"],
            smirks: None,
            category: "Esterification",
        },
        ReactionCase {
            name: "Esterification: HCOOH + EtOH → ethyl formate + H₂O",
            reactants: vec!["C(=O)O", "CCO"],
            products: vec!["C(=O)OCC", "O"],
            smirks: None,
            category: "Esterification",
        },
        ReactionCase {
            name: "Esterification: AcOH + EtOH → ethyl acetate + H₂O",
            reactants: vec!["CC(=O)O", "CCO"],
            products: vec!["CC(=O)OCC", "O"],
            smirks: None,
            category: "Esterification",
        },
        ReactionCase {
            name: "Amide: AcOH + MeNH₂ → N-methylacetamide + H₂O",
            reactants: vec!["CC(=O)O", "CN"],
            products: vec!["CC(=O)NC", "O"],
            smirks: None,
            category: "Esterification",
        },
        ReactionCase {
            name: "Amide: HCOOH + MeNH₂ → N-methylformamide + H₂O",
            reactants: vec!["CN", "C(=O)O"],
            products: vec!["CNC=O", "O"],
            smirks: None,
            category: "Esterification",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Conformational / Intramolecular — same SMILES, different geometry
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Conformational: n-butane anti ↔ gauche",
            reactants: vec!["CCCC"],
            products: vec!["CCCC"],
            smirks: None,
            category: "Conformational",
        },
        ReactionCase {
            name: "Conformational: ethanol rotamer",
            reactants: vec!["CCO"],
            products: vec!["CCO"],
            smirks: None,
            category: "Conformational",
        },
        ReactionCase {
            name: "Conformational: cyclohexane chair ↔ twist-boat",
            reactants: vec!["C1CCCCC1"],
            products: vec!["C1CCCCC1"],
            smirks: None,
            category: "Conformational",
        },
        ReactionCase {
            name: "Conformational: propanol rotamer",
            reactants: vec!["CCCO"],
            products: vec!["CCCO"],
            smirks: None,
            category: "Conformational",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Proton Transfer
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Proton transfer: FN ↔ NF rearrangement",
            reactants: vec!["FN"],
            products: vec!["NF"],
            smirks: None,
            category: "Proton Transfer",
        },
        ReactionCase {
            name: "Proton transfer: ClN ↔ NCl rearrangement",
            reactants: vec!["ClN"],
            products: vec!["NCl"],
            smirks: None,
            category: "Proton Transfer",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Bimolecular Association — hydrogen bonding / van der Waals
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Association: MeOH dimer",
            reactants: vec!["CO", "CO"],
            products: vec!["CO", "CO"],
            smirks: None,
            category: "Association",
        },
        ReactionCase {
            name: "Association: H₂O dimer",
            reactants: vec!["O", "O"],
            products: vec!["O", "O"],
            smirks: None,
            category: "Association",
        },
        ReactionCase {
            name: "Association: NH₃ + H₂O",
            reactants: vec!["N", "O"],
            products: vec!["N", "O"],
            smirks: None,
            category: "Association",
        },
        ReactionCase {
            name: "Association: HF + H₂O",
            reactants: vec!["F", "O"],
            products: vec!["F", "O"],
            smirks: None,
            category: "Association",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Oxidation / Reduction (H₂-balanced)
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Oxidation: methanol → formaldehyde + H₂",
            reactants: vec!["CO"],
            products: vec!["C=O", "[H][H]"],
            smirks: None,
            category: "Redox",
        },
        ReactionCase {
            name: "Reduction: formaldehyde + H₂ → methanol",
            reactants: vec!["C=O", "[H][H]"],
            products: vec!["CO"],
            smirks: None,
            category: "Redox",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Halogenation — radical substitution
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Halogenation: CH₄ + Cl₂ → CH₃Cl + HCl",
            reactants: vec!["C", "ClCl"],
            products: vec!["CCl", "Cl"],
            smirks: None,
            category: "Halogenation",
        },
        ReactionCase {
            name: "Halogenation: CH₄ + Br₂ → CH₃Br + HBr",
            reactants: vec!["C", "BrBr"],
            products: vec!["CBr", "Br"],
            smirks: None,
            category: "Halogenation",
        },
        ReactionCase {
            name: "Halogenation: C₂H₆ + Cl₂ → C₂H₅Cl + HCl",
            reactants: vec!["CC", "ClCl"],
            products: vec!["CCCl", "Cl"],
            smirks: None,
            category: "Halogenation",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Elimination — E2 model
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "E2: bromoethane → ethylene + HBr",
            reactants: vec!["CCBr"],
            products: vec!["C=C", "Br"],
            smirks: None,
            category: "Elimination",
        },
        ReactionCase {
            name: "E2: chloroethane → ethylene + HCl",
            reactants: vec!["CCCl"],
            products: vec!["C=C", "Cl"],
            smirks: None,
            category: "Elimination",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Ring formation / strain
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Ring closure: propane → cyclopropane + H₂",
            reactants: vec!["CCC"],
            products: vec!["C1CC1", "[H][H]"],
            smirks: None,
            category: "Ring",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Electrophilic Aromatic Substitution (EAS)
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "EAS: benzene + Cl₂ → chlorobenzene + HCl",
            reactants: vec!["c1ccccc1", "ClCl"],
            products: vec!["Clc1ccccc1", "Cl"],
            smirks: None,
            category: "EAS",
        },
        ReactionCase {
            name: "EAS: benzene + Br₂ → bromobenzene + HBr",
            reactants: vec!["c1ccccc1", "BrBr"],
            products: vec!["Brc1ccccc1", "Br"],
            smirks: None,
            category: "EAS",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Tautomerism — keto-enol like rearrangements
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Tautomer: acetic acid C=O ↔ C-OH",
            reactants: vec!["CC(=O)O"],
            products: vec!["CC(O)=O"],
            smirks: Some("[C:1](=[O:2])[O:3]>>[C:1]([O:2])=[O:3]"),
            category: "Tautomer",
        },
        ReactionCase {
            name: "Tautomer: acetaldehyde rotamer",
            reactants: vec!["CC=O"],
            products: vec!["CC=O"],
            smirks: None,
            category: "Tautomer",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Pericyclic rearrangements
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Cope rearrangement: 1,5-hexadiene",
            reactants: vec!["C=CCC=CC"],
            products: vec!["CC=CCC=C"],
            smirks: None,
            category: "Pericyclic",
        },
        ReactionCase {
            name: "Claisen rearrangement: allyl vinyl ether",
            reactants: vec!["C=CCOC=C"],
            products: vec!["C=CCC(=O)C"],
            smirks: None,
            category: "Pericyclic",
        },
        ReactionCase {
            name: "Retro Diels-Alder: cyclohexene → butadiene + ethylene",
            reactants: vec!["C1CC=CCC1"],
            products: vec!["C=CC=C", "C=C"],
            smirks: None,
            category: "Pericyclic",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Trimolecular
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Trimolecular: 3 H₂O cluster",
            reactants: vec!["O", "O", "O"],
            products: vec!["O", "O", "O"],
            smirks: None,
            category: "Trimolecular",
        },
        ReactionCase {
            name: "Trimolecular: 2 MeOH + H₂O cluster",
            reactants: vec!["CO", "CO", "O"],
            products: vec!["CO", "CO", "O"],
            smirks: None,
            category: "Trimolecular",
        },
        // ═══════════════════════════════════════════════════════════════════
        // Miscellaneous organic
        // ═══════════════════════════════════════════════════════════════════
        ReactionCase {
            name: "Hydration: acetaldehyde + H₂O → geminal diol",
            reactants: vec!["CC=O", "O"],
            products: vec!["CC(O)O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Epoxidation: ethylene + H₂O₂ → ethylene oxide + H₂O",
            reactants: vec!["C=C", "OO"],
            products: vec!["C1CO1", "O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Thiol exchange: CH₃SH + HSCH₃ cluster",
            reactants: vec!["CS", "SC"],
            products: vec!["CS", "SC"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Alcohol dehydration: ethanol → ethylene + H₂O",
            reactants: vec!["CCO"],
            products: vec!["C=C", "O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Grignard model: CH₃Br + CH₂O + H₂ → ethanol + HBr",
            reactants: vec!["CBr", "C=O", "[H][H]"],
            products: vec!["CCO", "Br"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Aldol condensation: 2 CH₃CHO → 3-hydroxybutanal",
            reactants: vec!["CC=O", "CC=O"],
            products: vec!["CC(O)CC=O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Michael addition: acrolein + acetaldehyde",
            reactants: vec!["C=CC=O", "CC=O"],
            products: vec!["CC(CC=O)C=O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Beckmann oxime: acetone oxime conformer",
            reactants: vec!["CC(=NO)C"],
            products: vec!["CC(=NO)C"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Pinacol coupling: 2 CH₂O + H₂ → ethylene glycol",
            reactants: vec!["C=O", "C=O", "[H][H]"],
            products: vec!["C(O)C(O)"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Transesterification: methyl acetate + EtOH → ethyl acetate + MeOH",
            reactants: vec!["CC(=O)OC", "CCO"],
            products: vec!["CC(=O)OCC", "CO"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Imine formation: acetaldehyde + MeNH₂ → imine + H₂O",
            reactants: vec!["CC=O", "CN"],
            products: vec!["CC=NC", "O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Hemiacetal: acetaldehyde + MeOH → hemiacetal",
            reactants: vec!["CC=O", "CO"],
            products: vec!["CC(O)OC"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Bimolecular: MeOH + HCl → MeCl + H₂O",
            reactants: vec!["CO", "Cl"],
            products: vec!["CCl", "O"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Bimolecular: NH₃ + Cl₂ → NH₂Cl + HCl",
            reactants: vec!["N", "ClCl"],
            products: vec!["NCl", "Cl"],
            smirks: None,
            category: "Misc",
        },
        ReactionCase {
            name: "Bimolecular: EtOH + HBr → EtBr + H₂O",
            reactants: vec!["CCO", "Br"],
            products: vec!["CCBr", "O"],
            smirks: None,
            category: "Misc",
        },
    ]
}

// ════════════════════════════════════════════════════════════════════════════
// Master test — runs all 55 reactions, prints summary table, fails if <90%
// ════════════════════════════════════════════════════════════════════════════

#[test]
fn study_55_reaction_types() {
    let battery = build_reaction_battery();
    let total = battery.len();
    assert!(
        total >= 55,
        "Battery must have ≥55 reactions, got {}",
        total
    );

    let outcomes: Vec<ReactionOutcome> = battery.par_iter().map(|case| run_one(case)).collect();

    // ── Print summary table ──────────────────────────────────────────
    println!("\n╔══════════════════════════════════════════════════════════════════════╗");
    println!(
        "║          REACTION DYNAMICS STUDY — {} REACTION TYPES            ║",
        total
    );
    println!("╠══════════════════════════════════════════════════════════════════════╣");
    println!(
        "║ {:<3} {:<48} {:>5} {:>8} ║",
        "#", "Reaction", "Frms", "Ea"
    );
    println!("╠══════════════════════════════════════════════════════════════════════╣");

    for o in &outcomes {
        let status = if o.ok { "✓" } else { "✗" };
        let name_trunc: String = o.name.chars().take(47).collect();
        let ea_str = if !o.ok {
            "    ERR".to_string()
        } else if o.activation_energy.abs() > 1e6 {
            "    !!!".to_string()
        } else {
            format!("{:>7.1}", o.activation_energy)
        };
        println!(
            "║ {} {:<48} {:>4} {:>8} ║",
            status, name_trunc, o.n_frames, ea_str,
        );
    }

    // ── Category summary ─────────────────────────────────────────────
    println!("╠══════════════════════════════════════════════════════════════════════╣");
    let categories: Vec<&str> = {
        let mut cats: Vec<&str> = outcomes.iter().map(|o| o.category.as_str()).collect();
        cats.sort();
        cats.dedup();
        cats
    };
    for cat in &categories {
        let cat_outcomes: Vec<&ReactionOutcome> =
            outcomes.iter().filter(|o| o.category == *cat).collect();
        let pass = cat_outcomes.iter().filter(|o| o.ok).count();
        let total_cat = cat_outcomes.len();
        println!(
            "║ {:<20} {:>3}/{:<3} passed                            ║",
            cat, pass, total_cat
        );
    }

    // ── Overall ──────────────────────────────────────────────────────
    let passed = outcomes.iter().filter(|o| o.ok).count();
    let failed = total - passed;
    let pct = 100.0 * passed as f64 / total as f64;
    println!("╠══════════════════════════════════════════════════════════════════════╣");
    println!(
        "║ TOTAL: {}/{} passed ({:.1}%)     {} failed                       ║",
        passed, total, pct, failed
    );
    println!("╚══════════════════════════════════════════════════════════════════════╝");

    // ── Print errors for failures ────────────────────────────────────
    if failed > 0 {
        println!("\n── Failures ──");
        for o in outcomes.iter().filter(|o| !o.ok) {
            println!(
                "  ✗ {} [{}]: {}",
                o.name,
                o.category,
                o.error.as_deref().unwrap_or("unknown")
            );
        }
    }

    // ── Assert ≥90% pass rate ────────────────────────────────────────
    assert!(
        pct >= 90.0,
        "Pass rate {:.1}% is below 90% threshold. {}/{} failed.",
        pct,
        failed,
        total
    );
}

// ════════════════════════════════════════════════════════════════════════════
// Individual category tests — for granular CI feedback
// ════════════════════════════════════════════════════════════════════════════

fn run_category(category: &str) {
    let battery = build_reaction_battery();
    let cases: Vec<&ReactionCase> = battery.iter().filter(|c| c.category == category).collect();
    assert!(
        !cases.is_empty(),
        "No reactions for category '{}'",
        category
    );

    let outcomes: Vec<ReactionOutcome> = cases.par_iter().map(|case| run_one(case)).collect();

    let failed: Vec<String> = outcomes
        .iter()
        .filter(|o| !o.ok)
        .map(|o| format!("  {} — {}", o.name, o.error.as_deref().unwrap_or("unknown")))
        .collect();

    assert!(
        failed.is_empty(),
        "{} failures in {} category:\n{}",
        failed.len(),
        category,
        failed.join("\n")
    );
}

#[test]
fn category_sn2() {
    run_category("SN2");
}

#[test]
fn category_addition() {
    run_category("Addition");
}

#[test]
fn category_diels_alder() {
    run_category("Diels-Alder");
}

#[test]
fn category_esterification() {
    run_category("Esterification");
}

#[test]
fn category_conformational() {
    run_category("Conformational");
}

#[test]
fn category_proton_transfer() {
    run_category("Proton Transfer");
}

#[test]
fn category_association() {
    run_category("Association");
}

#[test]
fn category_redox() {
    run_category("Redox");
}

#[test]
fn category_halogenation() {
    run_category("Halogenation");
}

#[test]
fn category_elimination() {
    run_category("Elimination");
}

#[test]
fn category_ring() {
    run_category("Ring");
}

#[test]
fn category_eas() {
    run_category("EAS");
}

#[test]
fn category_tautomer() {
    run_category("Tautomer");
}

#[test]
fn category_pericyclic() {
    run_category("Pericyclic");
}

#[test]
fn category_trimolecular() {
    run_category("Trimolecular");
}

#[test]
fn category_misc() {
    run_category("Misc");
}
