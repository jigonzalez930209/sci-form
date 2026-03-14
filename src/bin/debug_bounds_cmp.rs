use sci_form::graph::Molecule;

fn main() {
    let smi = "C#CCOCC#C";
    let mol = Molecule::from_smiles(smi).unwrap();
    let n = mol.graph.node_count();
    println!("Our bounds for {} (n={}):", smi, n);

    let bounds = sci_form::distgeom::calculate_bounds_matrix_opts(&mol, true);
    let mut b = bounds.clone();
    sci_form::distgeom::triangle_smooth_tol(&mut b, 0.05);

    // Print heavy atom bounds (skip H=1)
    let heavy: Vec<usize> = (0..n)
        .filter(|&i| {
            let ni = petgraph::graph::NodeIndex::new(i);
            mol.graph[ni].element != 1
        })
        .collect();

    println!("Heavy atoms: {:?}", heavy);
    for ii in 0..heavy.len() {
        for jj in (ii + 1)..heavy.len() {
            let i = heavy[ii];
            let j = heavy[jj];
            let ni = petgraph::graph::NodeIndex::new(i);
            let nj = petgraph::graph::NodeIndex::new(j);
            let lo = b[(j, i)];
            let hi = b[(i, j)];
            let elem_i = match mol.graph[ni].element {
                6 => "C",
                7 => "N",
                8 => "O",
                _ => "?",
            };
            let elem_j = match mol.graph[nj].element {
                6 => "C",
                7 => "N",
                8 => "O",
                _ => "?",
            };
            println!(
                "  [{},{}] ({}-{}) [{:.3}, {:.3}]",
                i, j, elem_i, elem_j, lo, hi
            );
        }
    }
}
