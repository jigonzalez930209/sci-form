use petgraph::visit::EdgeRef;
use sci_form::distgeom::{
    calculate_bounds_matrix_opts, compute_initial_coords_rdkit, identify_chiral_sets,
    pick_rdkit_distances, triangle_smooth_tol, MinstdRand,
};
use sci_form::forcefield::bounds_ff::minimize_bfgs_rdkit;
use sci_form::forcefield::etkdg_3d::{build_etkdg_3d_ff_with_torsions, minimize_etkdg_3d_bfgs};
use sci_form::graph::Molecule;

fn dist(coords: &nalgebra::DMatrix<f64>, i: usize, j: usize) -> f64 {
    let dx = coords[(i, 0)] - coords[(j, 0)];
    let dy = coords[(i, 1)] - coords[(j, 1)];
    let dz = coords[(i, 2)] - coords[(j, 2)];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn dist_f32(coords: &nalgebra::DMatrix<f32>, i: usize, j: usize) -> f32 {
    let dx = coords[(i, 0)] - coords[(j, 0)];
    let dy = coords[(i, 1)] - coords[(j, 1)];
    let dz = coords[(i, 2)] - coords[(j, 2)];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn main() {
    let smi = "C#CCOCC#C";
    let mol = Molecule::from_smiles(smi).unwrap();
    let n = mol.graph.node_count();

    // Bonds to track
    let bond_pairs: Vec<(usize, usize)> = mol
        .graph
        .edge_references()
        .filter(|e| mol.graph[e.source()].element != 1 && mol.graph[e.target()].element != 1)
        .map(|e| (e.source().index(), e.target().index()))
        .collect();

    println!("Molecule: {} (n={})", smi, n);

    // Step 1: Bounds
    let bounds = {
        let raw = calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        triangle_smooth_tol(&mut b, 0.05);
        b
    };

    println!("\nBounds:");
    for &(i, j) in &bond_pairs {
        let lo = bounds[(j, i)];
        let hi = bounds[(i, j)];
        println!("  {}-{}: [{:.4}, {:.4}]", i, j, lo, hi);
    }

    // Print all heavy-atom pairwise bounds
    let heavy: Vec<usize> = (0..n)
        .filter(|&i| mol.graph[petgraph::graph::NodeIndex::new(i)].element != 1)
        .collect();
    println!("\nAll heavy-atom bounds:");
    for a in 0..heavy.len() {
        for b in (a + 1)..heavy.len() {
            let i = heavy[a];
            let j = heavy[b];
            let lo = bounds[(j.max(i), j.min(i))];
            let hi = bounds[(j.min(i), j.max(i))];
            println!("  d({},{}) [{:.6}, {:.6}]", i, j, lo, hi);
        }
    }

    // Step 2: Sample distances and embed
    let chiral_sets = identify_chiral_sets(&mol);
    let mut rng = MinstdRand::new(42);
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };

    let dists = pick_rdkit_distances(&mut rng, &bounds);
    println!("\nSampled distances:");
    for &(i, j) in &bond_pairs {
        let d = if i < j { dists[(i, j)] } else { dists[(j, i)] };
        println!("  {}-{}: {:.4}", i, j, d);
    }

    let mut coords = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim).unwrap();
    println!("\nAfter embedding:");
    for &(i, j) in &bond_pairs {
        println!("  {}-{}: {:.4}", i, j, dist(&coords, i, j));
    }

    // Step 3: Bounds FF minimization
    {
        let mut need_more = 1;
        let mut restarts = 0;
        while need_more != 0 && restarts < 5 {
            need_more =
                minimize_bfgs_rdkit(&mut coords, &bounds, &chiral_sets, 400, 1e-3, 5.0, 0.1, 1.0);
            restarts += 1;
        }
    }
    println!("\nAfter bounds FF:");
    for &(i, j) in &bond_pairs {
        println!("  {}-{}: {:.4}", i, j, dist(&coords, i, j));
    }

    // Step 4: Drop to 3D
    let coords3d = coords.columns(0, 3).into_owned();
    let coords3d_f32 = coords3d.map(|v| v as f32);

    // Step 5: ETKDG 3D FF — load CSD torsions
    let csd_torsions: Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> = {
        let data = std::fs::read_to_string("tests/fixtures/torsion_params.json").unwrap();
        let all: std::collections::HashMap<String, Vec<serde_json::Value>> =
            serde_json::from_str(&data).unwrap();
        if let Some(torsions) = all.get(smi) {
            torsions
                .iter()
                .map(|t| {
                    let atoms: Vec<usize> = t["atoms"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|a| a.as_u64().unwrap() as usize)
                        .collect();
                    let v: Vec<f64> = t["v"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_f64().unwrap())
                        .collect();
                    let signs: Vec<f64> = t["signs"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|x| x.as_f64().unwrap())
                        .collect();
                    sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                        i: atoms[0],
                        j: atoms[1],
                        k: atoms[2],
                        l: atoms[3],
                        signs: [signs[0], signs[1], signs[2], signs[3], signs[4], signs[5]],
                        v: [v[0], v[1], v[2], v[3], v[4], v[5]],
                    }
                })
                .collect()
        } else {
            vec![]
        }
    };
    println!("\nCSD torsions: {}", csd_torsions.len());
    for t in &csd_torsions {
        println!(
            "  [{},{},{},{}] V={:?} s={:?}",
            t.i, t.j, t.k, t.l, t.v, t.signs
        );
    }
    let ff = build_etkdg_3d_ff_with_torsions(&mol, &coords3d, &bounds, &csd_torsions);
    println!("\n3D FF 1-2 constraints:");
    for c in ff
        .dist_12
        .iter()
        .chain(ff.dist_13.iter())
        .chain(ff.dist_long.iter())
    {
        if bond_pairs.contains(&(c.i, c.j)) || bond_pairs.contains(&(c.j, c.i)) {
            println!(
                "  {}-{}: [{:.4}, {:.4}] k={}",
                c.i, c.j, c.min_len, c.max_len, c.k
            );
        }
    }

    // Energy before 3D FF
    let e_before = sci_form::forcefield::etkdg_3d::etkdg_3d_energy(&coords3d_f32, &mol, &ff);
    println!("\nEnergy before 3D FF: {:.4}", e_before);

    let refined = minimize_etkdg_3d_bfgs(&mol, &coords3d, &ff, 300, 1e-3);
    let refined_f32 = refined.map(|v| v as f32);

    let e_after = sci_form::forcefield::etkdg_3d::etkdg_3d_energy(&refined_f32, &mol, &ff);
    println!("Energy after 3D FF: {:.4}", e_after);

    println!("\nFinal coordinates:");
    for i in 0..n {
        println!(
            "  atom {:2} ({}): {:.4} {:.4} {:.4}",
            i,
            mol.graph[petgraph::graph::NodeIndex::new(i)].element,
            refined_f32[(i, 0)],
            refined_f32[(i, 1)],
            refined_f32[(i, 2)]
        );
    }

    println!("\nHeavy atom pairwise distances:");
    let heavy: Vec<usize> = (0..n)
        .filter(|&i| mol.graph[petgraph::graph::NodeIndex::new(i)].element != 1)
        .collect();
    for a in 0..heavy.len() {
        for b in (a + 1)..heavy.len() {
            let i = heavy[a];
            let j = heavy[b];
            let d = dist_f32(&refined_f32, i, j);
            print!("d({}-{})={:.4} ", i, j, d);
        }
        println!();
    }

    println!("\nFF details:");
    println!(
        "  dist constraints: {}",
        ff.dist_12.len() + ff.dist_13.len() + ff.dist_long.len()
    );
    println!("  angle constraints: {}", ff.angle_constraints.len());
    println!("  torsion contribs: {}", ff.torsion_contribs.len());
    println!("  inversion contribs: {}", ff.inversion_contribs.len());
    for tc in &ff.torsion_contribs {
        if tc.v.iter().any(|&x| x > 0.0) {
            println!(
                "  torsion [{},{},{},{}] V={:?} s={:?}",
                tc.i, tc.j, tc.k, tc.l, tc.v, tc.signs
            );
        }
    }

    println!("\nAfter 3D FF:");
    for &(i, j) in &bond_pairs {
        let d = dist_f32(&refined_f32, i, j);
        let lo = bounds[(j.max(i), j.min(i))];
        let hi = bounds[(j.min(i), j.max(i))];
        println!(
            "  {}-{}: {:.4} (bounds: [{:.4}, {:.4}], diff: {:.4})",
            i,
            j,
            d,
            lo,
            hi,
            d as f64 - (lo + hi) / 2.0
        );
    }

    // Also check if any 1-2 constraint is being violated
    println!("\n1-2 constraint violations:");
    for c in ff
        .dist_12
        .iter()
        .chain(ff.dist_13.iter())
        .chain(ff.dist_long.iter())
    {
        if c.k >= 99.0 {
            let d = dist_f32(&refined_f32, c.i, c.j);
            if (d as f64) < c.min_len || (d as f64) > c.max_len {
                println!(
                    "  {}-{}: d={:.4}, constraint=[{:.4},{:.4}], violation={:.4}",
                    c.i,
                    c.j,
                    d,
                    c.min_len,
                    c.max_len,
                    if (d as f64) < c.min_len {
                        c.min_len - d as f64
                    } else {
                        d as f64 - c.max_len
                    }
                );
            }
        }
    }
}
