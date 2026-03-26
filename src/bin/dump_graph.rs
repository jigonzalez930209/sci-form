//! Debug: parse a molecule and print the graph topology (atoms, bonds, properties).
//!
//! Category: debug

use petgraph::visit::EdgeRef;
use sci_form::graph::Molecule;

fn main() {
    let smiles = "C#CCOC(C)CC1CC2C3CCC(C)C(O)(C3)C2O1";
    let mol = Molecule::from_smiles(smiles).unwrap();
    let n = mol.graph.node_count();
    println!("Atoms: {}", n);
    println!("Atom details:");
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let atom = &mol.graph[ni];
        let nbs: Vec<_> = mol.graph.neighbors(ni).map(|x| x.index()).collect();
        println!(
            "  {:2}: {:3} hyb={:?} nbs={:?}",
            i, atom.element, atom.hybridization, nbs
        );
    }
    println!("Bonds: {}", mol.graph.edge_count());
    for edge in mol.graph.edge_references() {
        let bond = &mol.graph[edge.id()];
        println!(
            "  {}-{}: {:?}",
            edge.source().index(),
            edge.target().index(),
            bond.order
        );
    }
}
