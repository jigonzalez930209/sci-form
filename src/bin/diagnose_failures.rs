//! Quick diagnostic: find which SMILES fail embedding in the first 1000 molecules.

use std::io::BufRead;

fn main() {
    let path = std::path::Path::new("data/chembl_1k_practical_asc.smi");
    let file = std::fs::File::open(path).unwrap();
    let reader = std::io::BufReader::new(file);

    let mut total = 0usize;
    let mut failures = Vec::new();

    for line in reader.lines().take(1000) {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        total += 1;
        let smi = parts[0];
        let id = parts[1];

        let result = sci_form::embed(smi, 42);
        if let Some(err) = result.error {
            failures.push((smi.to_string(), id.to_string(), err));
        }
    }

    println!("Total: {total}, Failures: {}", failures.len());
    for (smi, id, err) in &failures {
        println!("  FAIL: {id} | {smi}");
        println!("    Error: {err}");
    }
}
