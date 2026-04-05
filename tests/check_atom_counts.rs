//! Temporary test to validate atom counts for reaction test suite.

#[test]
fn check_all_reaction_atom_counts() {
    let reactions: Vec<(&str, Vec<&str>, Vec<&str>)> = vec![
        // ═══ SN2 Reactions ═══
        ("SN2_Cl_CH3Br", vec!["[Cl-]", "CBr"], vec!["ClC", "[Br-]"]),
        ("SN2_OH_CH3Cl", vec!["[OH-]", "CCl"], vec!["CO", "[Cl-]"]),
        ("SN2_I_CH3Cl", vec!["[I-]", "CCl"], vec!["CI", "[Cl-]"]),
        ("SN2_HCN_CH3Br", vec!["C#N", "CBr"], vec!["CC#N", "Br"]),
        ("SN2_F_CH3Cl", vec!["[F-]", "CCl"], vec!["CF", "[Cl-]"]),
        // ═══ Addition Reactions ═══
        ("Addition_HBr_ethylene", vec!["C=C", "Br"], vec!["CCBr"]),
        ("Addition_HCl_ethylene", vec!["C=C", "Cl"], vec!["CCCl"]),
        ("Addition_H2O_ethylene", vec!["C=C", "O"], vec!["CCO"]),
        ("Addition_HF_ethylene", vec!["C=C", "F"], vec!["CCF"]),
        ("Addition_Br2_ethylene", vec!["C=C", "BrBr"], vec!["BrCCBr"]),
        // ═══ Diels-Alder ═══
        (
            "DielsAlder_butadiene_ethylene",
            vec!["C=CC=C", "C=C"],
            vec!["C1CC=CCC1"],
        ),
        (
            "DielsAlder_butadiene_propylene",
            vec!["C=CC=C", "CC=C"],
            vec!["CC1CC=CCC1"],
        ),
        // ═══ Esterification ═══
        (
            "Esterification_AcOH_MeOH",
            vec!["CC(=O)O", "CO"],
            vec!["CC(=O)OC", "O"],
        ),
        (
            "Esterification_HCOOH_EtOH",
            vec!["C(=O)O", "CCO"],
            vec!["C(=O)OCC", "O"],
        ),
        (
            "Esterification_AcOH_EtOH",
            vec!["CC(=O)O", "CCO"],
            vec!["CC(=O)OCC", "O"],
        ),
        // ═══ Conformational / Intramolecular ═══
        ("Conformational_butane", vec!["CCCC"], vec!["CCCC"]),
        ("Conformational_ethanol", vec!["CCO"], vec!["CCO"]),
        (
            "Conformational_cyclohexane",
            vec!["C1CCCCC1"],
            vec!["C1CCCCC1"],
        ),
        ("Conformational_propanol", vec!["CCCO"], vec!["CCCO"]),
        // ═══ Proton Transfer ═══
        ("ProtonTransfer_FN", vec!["FN"], vec!["NF"]),
        ("ProtonTransfer_ClN", vec!["ClN"], vec!["NCl"]),
        // ═══ Bimolecular Association ═══
        ("Association_MeOH_MeOH", vec!["CO", "CO"], vec!["CO", "CO"]),
        ("Association_H2O_H2O", vec!["O", "O"], vec!["O", "O"]),
        ("Association_NH3_H2O", vec!["N", "O"], vec!["N", "O"]),
        ("Association_HF_H2O", vec!["F", "O"], vec!["F", "O"]),
        // ═══ Oxidation / Reduction (model — must include H₂ to balance) ═══
        ("Oxidation_methanol", vec!["CO"], vec!["C=O", "[H][H]"]),
        ("Reduction_formaldehyde", vec!["C=O", "[H][H]"], vec!["CO"]),
        // ═══ Halogenation ═══
        (
            "Halogenation_methane_Cl",
            vec!["C", "ClCl"],
            vec!["CCl", "Cl"],
        ),
        (
            "Halogenation_methane_Br",
            vec!["C", "BrBr"],
            vec!["CBr", "Br"],
        ),
        (
            "Halogenation_ethane_Cl",
            vec!["CC", "ClCl"],
            vec!["CCCl", "Cl"],
        ),
        // ═══ Elimination ═══
        ("E2_model_bromoethane", vec!["CCBr"], vec!["C=C", "Br"]),
        ("E2_model_chloroethane", vec!["CCCl"], vec!["C=C", "Cl"]),
        // ═══ Ring formation (balance H₂) ═══
        (
            "Ring_propane_cyclopropane",
            vec!["CCC"],
            vec!["C1CC1", "[H][H]"],
        ),
        // ═══ Electrophilic aromatic substitution (balance HCl/HBr) ═══
        (
            "EAS_benzene_Cl",
            vec!["c1ccccc1", "ClCl"],
            vec!["Clc1ccccc1", "Cl"],
        ),
        (
            "EAS_benzene_Br",
            vec!["c1ccccc1", "BrBr"],
            vec!["Brc1ccccc1", "Br"],
        ),
        // ═══ Tautomerism ═══
        ("Tautomer_AcOH", vec!["CC(=O)O"], vec!["CC(O)=O"]),
        ("Tautomer_acetaldehyde", vec!["CC=O"], vec!["CC=O"]),
        // ═══ Condensation ═══
        (
            "Condensation_amine_acid",
            vec!["CC(=O)O", "CN"],
            vec!["CC(=O)NC", "O"],
        ),
        // ═══ Various bimolecular ═══
        ("Bimolecular_CH3OH_HCl", vec!["CO", "Cl"], vec!["CCl", "O"]),
        ("Bimolecular_NH3_HCl", vec!["N", "ClCl"], vec!["NCl", "Cl"]),
        ("Bimolecular_EtOH_HBr", vec!["CCO", "Br"], vec!["CCBr", "O"]),
        // ═══ Trimolecular ═══
        (
            "Trimolecular_3_H2O",
            vec!["O", "O", "O"],
            vec!["O", "O", "O"],
        ),
        (
            "Trimolecular_2MeOH_H2O",
            vec!["CO", "CO", "O"],
            vec!["CO", "CO", "O"],
        ),
        // ═══ More reactions to reach 50+ ═══
        ("Hydration_acetaldehyde", vec!["CC=O", "O"], vec!["CC(O)O"]),
        (
            "Amide_MeNH2_HCOOH",
            vec!["CN", "C(=O)O"],
            vec!["CNC=O", "O"],
        ),
        (
            "Epoxide_ethylene_peracid",
            vec!["C=C", "OO"],
            vec!["C1CO1", "O"],
        ),
        ("Thiol_exchange", vec!["CS", "SC"], vec!["CS", "SC"]),
        ("Carbene_insertion", vec!["C=C"], vec!["C=C"]), // conformational model
        ("Aziridine_formation", vec!["C1CN1"], vec!["C1CN1"]), // conformational ring-strain model
        ("Cope_hexadiene", vec!["C=CCC=CC"], vec!["CC=CCC=C"]),
        ("Claisen_allyl_vinyl", vec!["C=CCOC=C"], vec!["C=CCC(=O)C"]),
        // ═══ Additional to reach 50+ ═══
        ("Alcohol_dehydration", vec!["CCO"], vec!["C=C", "O"]),
        (
            "Grignard_model",
            vec!["CBr", "C=O", "[H][H]"],
            vec!["CCO", "Br"],
        ),
        (
            "Aldol_condensation",
            vec!["CC=O", "CC=O"],
            vec!["CC(O)CC=O"],
        ),
        (
            "Michael_addition",
            vec!["C=CC=O", "CC=O"],
            vec!["CC(CC=O)C=O"],
        ),
        ("Beckmann_oxime", vec!["CC(=NO)C"], vec!["CC(=NO)C"]),
        (
            "Pinacol_coupling",
            vec!["C=O", "C=O", "[H][H]"],
            vec!["C(O)C(O)"],
        ),
        ("Retro_DielsAlder", vec!["C1CC=CCC1"], vec!["C=CC=C", "C=C"]),
        (
            "Transesterification",
            vec!["CC(=O)OC", "CCO"],
            vec!["CC(=O)OCC", "CO"],
        ),
        ("Imine_formation", vec!["CC=O", "CN"], vec!["CC=NC", "O"]),
        ("Hemiacetal_formation", vec!["CC=O", "CO"], vec!["CC(O)OC"]),
    ];

    let mut failures = Vec::new();
    for (name, reactants, products) in &reactions {
        let r_conf: Vec<_> = reactants.iter().map(|s| sci_form::embed(s, 42)).collect();
        let p_conf: Vec<_> = products.iter().map(|s| sci_form::embed(s, 42)).collect();

        let r_errs: Vec<_> = r_conf
            .iter()
            .enumerate()
            .filter(|(_, c)| c.error.is_some())
            .map(|(i, c)| {
                format!(
                    "R[{}] '{}': {:?}",
                    i,
                    reactants[i],
                    c.error.as_ref().unwrap()
                )
            })
            .collect();
        let p_errs: Vec<_> = p_conf
            .iter()
            .enumerate()
            .filter(|(_, c)| c.error.is_some())
            .map(|(i, c)| {
                format!(
                    "P[{}] '{}': {:?}",
                    i,
                    products[i],
                    c.error.as_ref().unwrap()
                )
            })
            .collect();

        if !r_errs.is_empty() || !p_errs.is_empty() {
            failures.push(format!(
                "EMBED_FAIL {}: {} {}",
                name,
                r_errs.join(", "),
                p_errs.join(", ")
            ));
            continue;
        }

        let r_total: usize = r_conf.iter().map(|c| c.num_atoms).sum();
        let p_total: usize = p_conf.iter().map(|c| c.num_atoms).sum();

        if r_total != p_total {
            let r_detail: Vec<String> = r_conf
                .iter()
                .enumerate()
                .map(|(i, c)| format!("{}={}", reactants[i], c.num_atoms))
                .collect();
            let p_detail: Vec<String> = p_conf
                .iter()
                .enumerate()
                .map(|(i, c)| format!("{}={}", products[i], c.num_atoms))
                .collect();
            failures.push(format!(
                "MISMATCH {}: R({})={} vs P({})={}",
                name,
                r_detail.join("+"),
                r_total,
                p_detail.join("+"),
                p_total
            ));
        } else {
            println!("OK {:40} : {} atoms", name, r_total);
        }
    }

    if !failures.is_empty() {
        for f in &failures {
            println!("FAIL: {}", f);
        }
        panic!(
            "{} reactions failed atom count validation:\n{}",
            failures.len(),
            failures.join("\n")
        );
    }
}
