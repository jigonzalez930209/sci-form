//! Output formatting and shared parse helpers.

pub fn format_xyz(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "0\nERROR: {}\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let mut out = format!("{}\n{}\n", n, result.smiles);
    for i in 0..n {
        let elem = element_symbol(result.elements[i]);
        let (x, y, z) = (
            result.coords[i * 3],
            result.coords[i * 3 + 1],
            result.coords[i * 3 + 2],
        );
        out.push_str(&format!("{:2} {:12.6} {:12.6} {:12.6}\n", elem, x, y, z));
    }
    out
}

pub fn format_sdf(result: &sci_form::ConformerResult) -> String {
    if result.error.is_some() || result.coords.is_empty() {
        return format!(
            "\n  sci-form\n\nERROR: {}\n$$$$\n",
            result.error.as_deref().unwrap_or("unknown")
        );
    }
    let n = result.num_atoms;
    let nb = result.bonds.len();
    let mut out = format!("{}\n  sci-form\n\n", result.smiles);
    out.push_str(&format!(
        "{:3}{:3}  0  0  0  0  0  0  0  0999 V2000\n",
        n, nb
    ));
    for i in 0..n {
        let (x, y, z) = (
            result.coords[i * 3],
            result.coords[i * 3 + 1],
            result.coords[i * 3 + 2],
        );
        let elem = element_symbol(result.elements[i]);
        out.push_str(&format!(
            "{:10.4}{:10.4}{:10.4} {:3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            x, y, z, elem
        ));
    }
    for (a, b, order) in &result.bonds {
        let bond_type = match order.as_str() {
            "SINGLE" => 1,
            "DOUBLE" => 2,
            "TRIPLE" => 3,
            "AROMATIC" => 4,
            _ => 1,
        };
        out.push_str(&format!("{:3}{:3}{:3}  0\n", a + 1, b + 1, bond_type));
    }
    out.push_str("M  END\n$$$$\n");
    out
}

pub fn element_symbol(atomic_num: u8) -> &'static str {
    match atomic_num {
        1 => "H",
        6 => "C",
        7 => "N",
        8 => "O",
        9 => "F",
        15 => "P",
        16 => "S",
        17 => "Cl",
        35 => "Br",
        53 => "I",
        5 => "B",
        14 => "Si",
        34 => "Se",
        _ => "X",
    }
}

/// Parse JSON elements + coords, exit on error.
pub fn parse_elems_coords(elements: &str, coords: &str) -> (Vec<u8>, Vec<f64>, Vec<[f64; 3]>) {
    let elems: Vec<u8> = serde_json::from_str(elements).unwrap_or_else(|e| {
        eprintln!("Bad elements JSON: {}", e);
        std::process::exit(1);
    });
    let flat: Vec<f64> = serde_json::from_str(coords).unwrap_or_else(|e| {
        eprintln!("Bad coords JSON: {}", e);
        std::process::exit(1);
    });
    let positions: Vec<[f64; 3]> = flat.chunks_exact(3).map(|c| [c[0], c[1], c[2]]).collect();
    (elems, flat, positions)
}
