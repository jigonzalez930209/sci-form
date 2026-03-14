use sci_form::smiles::SmilesParser;
use sci_form::graph::Molecule;
use sci_form::distgeom::bounds::calculate_bounds_matrix;

fn main() {
    let smiles_list = vec![
        "CCNCC#N",
        "CCC(=O)CCC=O",
        "CNc1ccccc1",
        "CC#CC1OC1C",
        "CC1(C)OC2CC1C2",
    ];

    for smi in smiles_list {
        let mut mol = Molecule::new("");
        let mut parser = SmilesParser::new(smi, &mut mol);
        parser.parse().unwrap();

        let bounds = calculate_bounds_matrix(&mol);
        let n = mol.graph.node_count();

        // Get heavy atom indices
        let heavy: Vec<usize> = (0..n)
            .filter(|&i| {
                let ni = petgraph::graph::NodeIndex::new(i);
                mol.graph[ni].element != 1
            })
            .collect();

        println!("OUR_BOUNDS {}", smi);
        for ii in 0..heavy.len() {
            for jj in (ii + 1)..heavy.len() {
                let i = heavy[ii];
                let j = heavy[jj];
                let (mi, ma) = if i < j { (i, j) } else { (j, i) };
                let lb = bounds[(ma, mi)];
                let ub = bounds[(mi, ma)];
                println!("  {} {} {:.6} {:.6}", ii, jj, lb, ub);
            }
        }
        println!();
    }
}
