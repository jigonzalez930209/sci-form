pub struct UffAtomParams {
    pub label: &'static str,
    pub r1: f64,     // Covalent radius (Angstroms)
    pub theta0: f64, // Equilibrium angle (Radians)
    pub x1: f64,     // VDW distance (Angstroms)
    pub d1: f64,     // VDW depth (kcal/mol)
    pub zeta: f64,   // Valence effective charge
    pub z_star: f64, // Effective charge for bonds
    pub chi: f64,    // Electronegativity
}

pub fn get_uff_params(label: &str) -> Option<UffAtomParams> {
    match label {
        "C_3" => Some(UffAtomParams {
            label: "C_3",
            r1: 0.706,
            theta0: 109.47,
            x1: 3.851,
            d1: 0.105,
            zeta: 1.912,
            z_star: 1.912,
            chi: 5.343,
        }),
        "H_" => Some(UffAtomParams {
            label: "H_",
            r1: 0.354,
            theta0: 180.0,
            x1: 2.886,
            d1: 0.044,
            zeta: 0.712,
            z_star: 0.712,
            chi: 4.528,
        }),
        _ => None,
    }
}
