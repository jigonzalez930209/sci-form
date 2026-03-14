use sci_form::conformer::generate_3d_conformer;
use sci_form::graph::Molecule;

fn main() {
    let smiles = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();
    println!("N={}", n);

    match generate_3d_conformer(&mol, 42) {
        Ok(coords) => {
            for i in 0..n {
                println!(
                    "{:2}: {:15.10} {:15.10} {:15.10}",
                    i,
                    coords[(i, 0)],
                    coords[(i, 1)],
                    coords[(i, 2)]
                );
            }
        }
        Err(e) => println!("Error: {}", e),
    }
}
