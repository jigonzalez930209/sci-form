#![allow(
    unused_imports,
    unused_variables,
    dead_code,
    clippy::unnecessary_cast,
    clippy::needless_range_loop,
    clippy::manual_repeat_n,
    clippy::manual_str_repeat,
    clippy::manual_is_multiple_of,
    clippy::redundant_field_names,
    clippy::useless_vec,
    clippy::single_range_in_vec_init
)]
use nalgebra::{DMatrix, Vector3};
use petgraph::visit::EdgeRef;
use rand::SeedableRng;
use serde::Deserialize;

#[derive(Deserialize)]
struct OracleAtom {
    element: u8,
    x: f32,
    y: f32,
    z: f32,
    formal_charge: i8,
    hybridization: String,
}

#[derive(Deserialize)]
struct OracleBond {
    start: usize,
    end: usize,
    order: String,
}

#[derive(Deserialize)]
struct OracleMolecule {
    smiles: String,
    atoms: Vec<OracleAtom>,
    bonds: Vec<OracleBond>,
}

fn build_mol(mol: &OracleMolecule) -> sci_form::graph::Molecule {
    let mut our_mol = sci_form::graph::Molecule::new(&mol.smiles);
    let mut node_indices = Vec::new();
    for atom in &mol.atoms {
        let hybridization = match atom.hybridization.as_str() {
            "SP" => sci_form::graph::Hybridization::SP,
            "SP2" => sci_form::graph::Hybridization::SP2,
            "SP3" => sci_form::graph::Hybridization::SP3,
            "SP3D" => sci_form::graph::Hybridization::SP3D,
            "SP3D2" => sci_form::graph::Hybridization::SP3D2,
            _ => sci_form::graph::Hybridization::Unknown,
        };
        let new_atom = sci_form::graph::Atom {
            element: atom.element,
            position: nalgebra::Vector3::zeros(),
            charge: 0.0,
            formal_charge: atom.formal_charge,
            hybridization,
            chiral_tag: sci_form::graph::ChiralType::Unspecified,
            explicit_h: if atom.element == 1 || atom.element == 0 {
                1
            } else {
                0
            },
        };
        node_indices.push(our_mol.add_atom(new_atom));
    }
    for bond in &mol.bonds {
        let order = match bond.order.as_str() {
            "DOUBLE" => sci_form::graph::BondOrder::Double,
            "TRIPLE" => sci_form::graph::BondOrder::Triple,
            "AROMATIC" => sci_form::graph::BondOrder::Aromatic,
            _ => sci_form::graph::BondOrder::Single,
        };
        our_mol.add_bond(
            node_indices[bond.start],
            node_indices[bond.end],
            sci_form::graph::Bond {
                order,
                stereo: sci_form::graph::BondStereo::None,
            },
        );
    }
    our_mol
}

/// Calculate energy breakdown by component
fn energy_breakdown(
    coords: &DMatrix<f32>,
    mol: &sci_form::graph::Molecule,
    params: &sci_form::forcefield::FFParams,
    bounds_matrix: &DMatrix<f64>,
) -> (f32, f32, f32, f32, f32, f32) {
    let n = mol.graph.node_count();
    let mut e_bonds = 0.0f32;
    let mut e_angles = 0.0f32;
    let mut e_bounds = 0.0f32;
    let mut e_oop = 0.0f32;
    let mut e_torsion_uff = 0.0f32;
    let mut e_torsion_etkdg = 0.0f32;

    // Bonds
    for edge in mol.graph.edge_references() {
        let idx1 = edge.source().index();
        let idx2 = edge.target().index();
        let p1 = Vector3::new(coords[(idx1, 0)], coords[(idx1, 1)], coords[(idx1, 2)]);
        let p2 = Vector3::new(coords[(idx2, 0)], coords[(idx2, 1)], coords[(idx2, 2)]);
        let r_eq = sci_form::distgeom::get_bond_length(mol, edge.source(), edge.target()) as f32;
        e_bonds += sci_form::forcefield::energy::bond_stretch_energy(&p1, &p2, params.kb, r_eq);
    }

    // Angles
    for i in 0..n {
        let ni = petgraph::graph::NodeIndex::new(i);
        let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
        for j in 0..nbs.len() {
            for k in (j + 1)..nbs.len() {
                let p1 = Vector3::new(
                    coords[(nbs[j].index(), 0)],
                    coords[(nbs[j].index(), 1)],
                    coords[(nbs[j].index(), 2)],
                );
                let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
                let p2 = Vector3::new(
                    coords[(nbs[k].index(), 0)],
                    coords[(nbs[k].index(), 1)],
                    coords[(nbs[k].index(), 2)],
                );
                let ideal =
                    sci_form::graph::get_corrected_ideal_angle(mol, ni, nbs[j], nbs[k]) as f32;
                e_angles += sci_form::forcefield::energy::angle_bend_energy(
                    &p1,
                    &pc,
                    &p2,
                    params.k_theta,
                    ideal,
                );
            }
        }
    }

    // Bounds
    for i in 0..n {
        for j in (i + 1)..n {
            let p1 = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p2 = Vector3::new(coords[(j, 0)], coords[(j, 1)], coords[(j, 2)]);
            let upper = bounds_matrix[(i, j)] as f32;
            let lower = bounds_matrix[(j, i)] as f32;
            e_bounds += sci_form::forcefield::energy::bounds_energy(
                &p1,
                &p2,
                lower,
                upper,
                params.k_bounds,
            );
        }
    }

    // OOP
    if params.k_oop.abs() > 1e-8 {
        for i in 0..n {
            let ni = petgraph::graph::NodeIndex::new(i);
            if mol.graph[ni].hybridization != sci_form::graph::Hybridization::SP2 {
                continue;
            }
            let nbs: Vec<_> = mol.graph.neighbors(ni).collect();
            if nbs.len() != 3 {
                continue;
            }
            let pc = Vector3::new(coords[(i, 0)], coords[(i, 1)], coords[(i, 2)]);
            let p1 = Vector3::new(
                coords[(nbs[0].index(), 0)],
                coords[(nbs[0].index(), 1)],
                coords[(nbs[0].index(), 2)],
            );
            let p2 = Vector3::new(
                coords[(nbs[1].index(), 0)],
                coords[(nbs[1].index(), 1)],
                coords[(nbs[1].index(), 2)],
            );
            let p3 = Vector3::new(
                coords[(nbs[2].index(), 0)],
                coords[(nbs[2].index(), 1)],
                coords[(nbs[2].index(), 2)],
            );
            let v1 = p1 - pc;
            let v2 = p2 - pc;
            let v3 = p3 - pc;
            let vol = v1.dot(&v2.cross(&v3));
            e_oop += params.k_oop * vol * vol;
        }
    }

    // UFF Torsions
    if params.k_omega.abs() > 1e-8 && n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            let hyb_u = mol.graph[u].hybridization;
            let hyb_v = mol.graph[v].hybridization;
            if hyb_u == sci_form::graph::Hybridization::SP
                || hyb_v == sci_form::graph::Hybridization::SP
            {
                continue;
            }
            let (n_fold, gamma, weight) =
                sci_form::forcefield::energy::torsion_params(hyb_u, hyb_v);
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            for &nu in &neighbors_u {
                for &nv in &neighbors_v {
                    let p1 = Vector3::new(
                        coords[(nu.index(), 0)],
                        coords[(nu.index(), 1)],
                        coords[(nu.index(), 2)],
                    );
                    let p2 = Vector3::new(
                        coords[(u.index(), 0)],
                        coords[(u.index(), 1)],
                        coords[(u.index(), 2)],
                    );
                    let p3 = Vector3::new(
                        coords[(v.index(), 0)],
                        coords[(v.index(), 1)],
                        coords[(v.index(), 2)],
                    );
                    let p4 = Vector3::new(
                        coords[(nv.index(), 0)],
                        coords[(nv.index(), 1)],
                        coords[(nv.index(), 2)],
                    );
                    e_torsion_uff += sci_form::forcefield::energy::torsional_energy(
                        &p1,
                        &p2,
                        &p3,
                        &p4,
                        params.k_omega * weight,
                        n_fold,
                        gamma,
                    );
                }
            }
        }
    }

    // ETKDG torsions
    if n >= 4 {
        for edge in mol.graph.edge_references() {
            let u = edge.source();
            let v = edge.target();
            if sci_form::graph::min_path_excluding2(mol, u, v, u, v, 7).is_some() {
                continue;
            }
            let m6 =
                sci_form::forcefield::etkdg_lite::infer_etkdg_parameters(mol, u.index(), v.index());
            if m6.v.iter().all(|&x| x.abs() < 1e-6) {
                continue;
            }
            let neighbors_u: Vec<_> = mol.graph.neighbors(u).filter(|&x| x != v).collect();
            let neighbors_v: Vec<_> = mol.graph.neighbors(v).filter(|&x| x != u).collect();
            if neighbors_u.is_empty() || neighbors_v.is_empty() {
                continue;
            }
            let nu = neighbors_u[0];
            let nv = neighbors_v[0];
            let p1 = Vector3::new(
                coords[(nu.index(), 0)],
                coords[(nu.index(), 1)],
                coords[(nu.index(), 2)],
            );
            let p2 = Vector3::new(
                coords[(u.index(), 0)],
                coords[(u.index(), 1)],
                coords[(u.index(), 2)],
            );
            let p3 = Vector3::new(
                coords[(v.index(), 0)],
                coords[(v.index(), 1)],
                coords[(v.index(), 2)],
            );
            let p4 = Vector3::new(
                coords[(nv.index(), 0)],
                coords[(nv.index(), 1)],
                coords[(nv.index(), 2)],
            );
            e_torsion_etkdg +=
                sci_form::forcefield::etkdg_lite::calc_torsion_energy_m6(&p1, &p2, &p3, &p4, &m6);
        }
    }

    (
        e_bonds,
        e_angles,
        e_bounds,
        e_oop,
        e_torsion_uff,
        e_torsion_etkdg,
    )
}

/// Calculate gradient norm breakdown by component
fn gradient_norms(
    coords: &DMatrix<f32>,
    mol: &sci_form::graph::Molecule,
    params: &sci_form::forcefield::FFParams,
    bounds_matrix: &DMatrix<f64>,
) -> (f32, f32, f32, f32, f32, f32) {
    let n = mol.graph.node_count();

    // Bond gradient only
    let bond_params = sci_form::forcefield::FFParams {
        kb: params.kb,
        k_theta: 0.0,
        k_omega: 0.0,
        k_oop: 0.0,
        k_bounds: 0.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };
    let g_bond = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        &bond_params,
        bounds_matrix,
    );

    // Angle gradient only
    let angle_params = sci_form::forcefield::FFParams {
        kb: 0.0,
        k_theta: params.k_theta,
        k_omega: 0.0,
        k_oop: 0.0,
        k_bounds: 0.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };
    let g_angle = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        &angle_params,
        bounds_matrix,
    );

    // Bounds gradient only
    let bounds_params = sci_form::forcefield::FFParams {
        kb: 0.0,
        k_theta: 0.0,
        k_omega: 0.0,
        k_oop: 0.0,
        k_bounds: params.k_bounds,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };
    let g_bounds = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        &bounds_params,
        bounds_matrix,
    );

    // OOP gradient only
    let oop_params = sci_form::forcefield::FFParams {
        kb: 0.0,
        k_theta: 0.0,
        k_omega: 0.0,
        k_oop: params.k_oop,
        k_bounds: 0.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };
    let g_oop = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        &oop_params,
        bounds_matrix,
    );

    // UFF torsion gradient only
    let torsion_params = sci_form::forcefield::FFParams {
        kb: 0.0,
        k_theta: 0.0,
        k_omega: params.k_omega,
        k_oop: 0.0,
        k_bounds: 0.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };
    let g_torsion = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        &torsion_params,
        bounds_matrix,
    );

    // ETKDG torsion gradient — need to isolate. Since compute_analytical_gradient always includes ETKDG M6,
    // calculate it by: full - (bond+angle+bounds+oop+torsion)
    let g_full = sci_form::forcefield::gradients::compute_analytical_gradient(
        coords,
        mol,
        params,
        bounds_matrix,
    );
    let g_etkdg = &g_full - &g_bond - &g_angle - &g_bounds - &g_oop - &g_torsion;

    (
        g_bond.norm(),
        g_angle.norm(),
        g_bounds.norm(),
        g_oop.norm(),
        g_torsion.norm(),
        g_etkdg.norm(),
    )
}

#[test]
fn debug_ff_energy_breakdown() {
    let data = std::fs::read_to_string("tests/fixtures/reference_coords.json").unwrap();
    let mut molecules: Vec<OracleMolecule> = serde_json::from_str(&data).unwrap();

    use rand::seq::SliceRandom;
    let mut rng = rand::rngs::StdRng::seed_from_u64(123);
    molecules.shuffle(&mut rng);
    molecules.truncate(100);

    let params = sci_form::forcefield::FFParams {
        kb: 300.0,
        k_theta: 200.0,
        k_omega: 10.0,
        k_oop: 20.0,
        k_bounds: 10.0,
        k_chiral: 0.0,
        k_vdw: 0.0,
    };

    // Test on worst molecules: 62, 48, 28, 92 (from previous analysis)
    let test_indices = vec![62, 48, 28, 92, 0, 10];

    for &mol_idx in &test_indices {
        if mol_idx >= molecules.len() {
            continue;
        }
        let mol_data = &molecules[mol_idx];
        let our_mol = build_mol(mol_data);
        let n = our_mol.graph.node_count();

        let bounds = sci_form::distgeom::calculate_bounds_matrix(&our_mol);
        let smoothed = sci_form::distgeom::smooth_bounds_matrix(bounds);

        // Generate embedding with seed 0
        let mut rng = rand::rngs::StdRng::seed_from_u64(13);
        let dists = sci_form::distgeom::pick_random_distances(&mut rng, &smoothed);
        let metric = sci_form::distgeom::compute_metric_matrix(&dists);
        let mut coords4d = sci_form::distgeom::generate_nd_coordinates(&mut rng, &metric, 4);

        sci_form::forcefield::minimize_bounds_lbfgs_ex(
            &mut coords4d,
            &smoothed,
            &[],
            400,
            1e-4,
            5.0,
            0.1,
        );
        sci_form::forcefield::minimize_bounds_lbfgs_ex(
            &mut coords4d,
            &smoothed,
            &[],
            200,
            1e-4,
            5.0,
            1.0,
        );

        let coords_before = coords4d.columns(0, 3).into_owned();

        println!("\n=== Mol {} ({}) n={} ===", mol_idx, &mol_data.smiles, n);

        // Energy breakdown BEFORE FF
        let (e_b, e_a, e_bnd, e_oop, e_t_uff, e_t_etkdg) =
            energy_breakdown(&coords_before, &our_mol, &params, &smoothed);
        println!("BEFORE FF:");
        println!("  Bonds:    {:.4}", e_b);
        println!("  Angles:   {:.4}", e_a);
        println!("  Bounds:   {:.4}", e_bnd);
        println!("  OOP:      {:.4}", e_oop);
        println!("  Tors UFF: {:.4}", e_t_uff);
        println!("  Tors ETKDG: {:.4}", e_t_etkdg);
        println!(
            "  TOTAL:    {:.4}",
            e_b + e_a + e_bnd + e_oop + e_t_uff + e_t_etkdg
        );

        // Gradient norms
        let (gn_b, gn_a, gn_bnd, gn_oop, gn_t_uff, gn_t_etkdg) =
            gradient_norms(&coords_before, &our_mol, &params, &smoothed);
        println!("Gradient norms:");
        println!("  Bonds:    {:.4}", gn_b);
        println!("  Angles:   {:.4}", gn_a);
        println!("  Bounds:   {:.4}", gn_bnd);
        println!("  OOP:      {:.4}", gn_oop);
        println!("  Tors UFF: {:.4}", gn_t_uff);
        println!("  Tors ETKDG: {:.4}", gn_t_etkdg);

        // Full gradient max component
        let g_full = sci_form::forcefield::gradients::compute_analytical_gradient(
            &coords_before,
            &our_mol,
            &params,
            &smoothed,
        );
        println!("  Full norm:    {:.4}", g_full.norm());
        println!("  Full max:     {:.4}", g_full.abs().max());

        // Run NEW 3D ETKDG minimizer
        let coords_after =
            sci_form::forcefield::minimize_etkdg_3d(&our_mol, &coords_before, &smoothed, 200, 1e-4);

        // Gradient consistency check for worst molecule
        if mol_idx == 62 {
            let ff = sci_form::forcefield::build_etkdg_3d_ff(
                &our_mol,
                &coords_before.map(|v| v as f64),
                &smoothed,
            );
            let anal_g = sci_form::forcefield::etkdg_3d_gradient(&coords_before, &our_mol, &ff);
            let e0 = sci_form::forcefield::etkdg_3d_energy(&coords_before, &our_mol, &ff);
            let eps = 1e-4f32;

            // Check individual energy components
            println!("  COMPONENT GRADIENT CHECK (Mol 62):");

            // Test distance constraints alone
            let mut g_dist = nalgebra::DMatrix::zeros(n, 3);
            for c in ff
                .dist_12
                .iter()
                .chain(ff.dist_13.iter())
                .chain(ff.dist_long.iter())
            {
                let p1 = nalgebra::Vector3::new(
                    coords_before[(c.i, 0)],
                    coords_before[(c.i, 1)],
                    coords_before[(c.i, 2)],
                );
                let p2 = nalgebra::Vector3::new(
                    coords_before[(c.j, 0)],
                    coords_before[(c.j, 1)],
                    coords_before[(c.j, 2)],
                );
                sci_form::forcefield::gradients::analytical_grad_distance_constraint(
                    &p1,
                    &p2,
                    c.min_len as f32,
                    c.max_len as f32,
                    c.k as f32,
                    &mut g_dist,
                    c.i,
                    c.j,
                );
            }

            // Numerical gradient for dist constraints only
            for atom in 0..2 {
                for d in 0..3 {
                    let mut cp = coords_before.clone();
                    let mut cm = coords_before.clone();
                    cp[(atom, d)] += eps;
                    cm[(atom, d)] -= eps;
                    let mut ep = 0.0f32;
                    let mut em = 0.0f32;
                    for c in ff
                        .dist_12
                        .iter()
                        .chain(ff.dist_13.iter())
                        .chain(ff.dist_long.iter())
                    {
                        let p1p = nalgebra::Vector3::new(cp[(c.i, 0)], cp[(c.i, 1)], cp[(c.i, 2)]);
                        let p2p = nalgebra::Vector3::new(cp[(c.j, 0)], cp[(c.j, 1)], cp[(c.j, 2)]);
                        ep += sci_form::forcefield::energy::distance_constraint_energy(
                            &p1p,
                            &p2p,
                            c.min_len as f32,
                            c.max_len as f32,
                            c.k as f32,
                        );
                        let p1m = nalgebra::Vector3::new(cm[(c.i, 0)], cm[(c.i, 1)], cm[(c.i, 2)]);
                        let p2m = nalgebra::Vector3::new(cm[(c.j, 0)], cm[(c.j, 1)], cm[(c.j, 2)]);
                        em += sci_form::forcefield::energy::distance_constraint_energy(
                            &p1m,
                            &p2m,
                            c.min_len as f32,
                            c.max_len as f32,
                            c.k as f32,
                        );
                    }
                    let num_g = (ep - em) / (2.0 * eps);
                    let anal = g_dist[(atom, d)];
                    let axis = ["x", "y", "z"][d];
                    println!(
                        "    DIST atom={} {}: anal={:.4} num={:.4} diff={:.4}",
                        atom,
                        axis,
                        anal,
                        num_g,
                        (anal - num_g).abs()
                    );
                }
            }

            // Test single torsion term numerically
            println!("    --- SINGLE TORSION TEST ---");
            // Find first SP3-SP3 bond's first torsion
            'outer: for edge in our_mol.graph.edge_references() {
                let u = edge.source();
                let v = edge.target();
                let hyb_u = our_mol.graph[u].hybridization;
                let hyb_v = our_mol.graph[v].hybridization;
                if hyb_u == sci_form::graph::Hybridization::SP
                    || hyb_v == sci_form::graph::Hybridization::SP
                {
                    continue;
                }
                let (n_fold, gamma, weight) =
                    sci_form::forcefield::energy::torsion_params(hyb_u, hyb_v);
                let neighbors_u: Vec<_> = our_mol.graph.neighbors(u).filter(|&x| x != v).collect();
                let neighbors_v: Vec<_> = our_mol.graph.neighbors(v).filter(|&x| x != u).collect();
                if neighbors_u.is_empty() || neighbors_v.is_empty() {
                    continue;
                }
                let nu = neighbors_u[0];
                let nv = neighbors_v[0];
                let vv = 10.0 * weight;
                println!(
                    "    Torsion: {}-{}-{}-{} hyb={:?}-{:?} n={} gamma={:.2} V={:.2}",
                    nu.index(),
                    u.index(),
                    v.index(),
                    nv.index(),
                    hyb_u,
                    hyb_v,
                    n_fold,
                    gamma,
                    vv
                );

                let mk = |coords: &nalgebra::DMatrix<f32>, idx: usize| {
                    nalgebra::Vector3::new(coords[(idx, 0)], coords[(idx, 1)], coords[(idx, 2)])
                };

                // Analytical gradient for this single torsion
                let mut single_grad = nalgebra::DMatrix::zeros(n, 3);
                sci_form::forcefield::gradients::analytical_grad_torsion(
                    &mk(&coords_before, nu.index()),
                    &mk(&coords_before, u.index()),
                    &mk(&coords_before, v.index()),
                    &mk(&coords_before, nv.index()),
                    vv,
                    n_fold,
                    gamma,
                    &mut single_grad,
                    nu.index(),
                    u.index(),
                    v.index(),
                    nv.index(),
                );

                // Numerical gradient for this single torsion
                let indices = [nu.index(), u.index(), v.index(), nv.index()];
                for &idx in &indices {
                    for d in 0..3 {
                        let mut cp = coords_before.clone();
                        let mut cm = coords_before.clone();
                        cp[(idx, d)] += eps;
                        cm[(idx, d)] -= eps;
                        let ep = sci_form::forcefield::energy::torsional_energy(
                            &mk(&cp, nu.index()),
                            &mk(&cp, u.index()),
                            &mk(&cp, v.index()),
                            &mk(&cp, nv.index()),
                            vv,
                            n_fold,
                            gamma,
                        );
                        let em = sci_form::forcefield::energy::torsional_energy(
                            &mk(&cm, nu.index()),
                            &mk(&cm, u.index()),
                            &mk(&cm, v.index()),
                            &mk(&cm, nv.index()),
                            vv,
                            n_fold,
                            gamma,
                        );
                        let num_g = (ep - em) / (2.0 * eps);
                        let anal = single_grad[(idx, d)];
                        let axis = ["x", "y", "z"][d];
                        let ok = if num_g.abs() > 0.01 {
                            (anal / num_g - 1.0).abs() < 0.1
                        } else {
                            (anal - num_g).abs() < 0.1
                        };
                        if !ok {
                            println!(
                                "    FAIL atom={} {} anal={:.6} num={:.6} ratio={:.4}",
                                idx,
                                axis,
                                anal,
                                num_g,
                                if num_g.abs() > 1e-6 {
                                    anal / num_g
                                } else {
                                    0.0
                                }
                            );
                        }
                    }
                }
                println!("    Single torsion check done");
                break 'outer;
            }

            // Test ETKDG M6 torsion gradient
            println!("    --- ETKDG M6 ---");
            // Get ETKDG-only gradient by using no-other-terms params
            let no_params = sci_form::forcefield::FFParams {
                kb: 0.0,
                k_theta: 0.0,
                k_omega: 0.0,
                k_oop: 0.0,
                k_bounds: 0.0,
                k_chiral: 0.0,
                k_vdw: 0.0,
            };
            let g_etkdg = sci_form::forcefield::gradients::compute_analytical_gradient(
                &coords_before,
                &our_mol,
                &no_params,
                &smoothed,
            );
            for atom in 0..2 {
                for d in 0..3 {
                    let mut cp = coords_before.clone();
                    let mut cm = coords_before.clone();
                    cp[(atom, d)] += eps;
                    cm[(atom, d)] -= eps;
                    let mut ep = 0.0f32;
                    let mut em = 0.0f32;
                    for edge in our_mol.graph.edge_references() {
                        let u = edge.source();
                        let v = edge.target();
                        if sci_form::graph::min_path_excluding2(&our_mol, u, v, u, v, 7).is_some() {
                            continue;
                        }
                        let m6 = sci_form::forcefield::etkdg_lite::infer_etkdg_parameters(
                            &our_mol,
                            u.index(),
                            v.index(),
                        );
                        if m6.v.iter().all(|&x| x.abs() < 1e-6) {
                            continue;
                        }
                        let neighbors_u: Vec<_> =
                            our_mol.graph.neighbors(u).filter(|&x| x != v).collect();
                        let neighbors_v: Vec<_> =
                            our_mol.graph.neighbors(v).filter(|&x| x != u).collect();
                        if neighbors_u.is_empty() || neighbors_v.is_empty() {
                            continue;
                        }
                        let nu = neighbors_u[0];
                        let nv = neighbors_v[0];
                        let mk = |coords: &nalgebra::DMatrix<f32>, idx: usize| {
                            nalgebra::Vector3::new(
                                coords[(idx, 0)],
                                coords[(idx, 1)],
                                coords[(idx, 2)],
                            )
                        };
                        ep += sci_form::forcefield::etkdg_lite::calc_torsion_energy_m6(
                            &mk(&cp, nu.index()),
                            &mk(&cp, u.index()),
                            &mk(&cp, v.index()),
                            &mk(&cp, nv.index()),
                            &m6,
                        );
                        em += sci_form::forcefield::etkdg_lite::calc_torsion_energy_m6(
                            &mk(&cm, nu.index()),
                            &mk(&cm, u.index()),
                            &mk(&cm, v.index()),
                            &mk(&cm, nv.index()),
                            &m6,
                        );
                    }
                    let num_g = (ep - em) / (2.0 * eps);
                    let anal = g_etkdg[(atom, d)];
                    let axis = ["x", "y", "z"][d];
                    println!(
                        "    M6 atom={} {}: anal={:.4} num={:.4} diff={:.4}",
                        atom,
                        axis,
                        anal,
                        num_g,
                        (anal - num_g).abs()
                    );
                }
            }

            // Test OOP gradient
            println!("    --- OOP ---");
            let oop_params = sci_form::forcefield::FFParams {
                kb: 0.0,
                k_theta: 0.0,
                k_omega: 0.0,
                k_oop: 10.0,
                k_bounds: 0.0,
                k_chiral: 0.0,
                k_vdw: 0.0,
            };
            let g_oop = sci_form::forcefield::gradients::compute_analytical_gradient(
                &coords_before,
                &our_mol,
                &oop_params,
                &smoothed,
            );
            for atom in 0..2 {
                for d in 0..3 {
                    let mut cp = coords_before.clone();
                    let mut cm = coords_before.clone();
                    cp[(atom, d)] += eps;
                    cm[(atom, d)] -= eps;
                    // OOP energy
                    let mut ep = 0.0f32;
                    let mut em = 0.0f32;
                    for i in 0..n {
                        let ni = petgraph::graph::NodeIndex::new(i);
                        if our_mol.graph[ni].hybridization != sci_form::graph::Hybridization::SP2 {
                            continue;
                        }
                        let nbs: Vec<_> = our_mol.graph.neighbors(ni).collect();
                        if nbs.len() != 3 {
                            continue;
                        }
                        let mk = |coords: &nalgebra::DMatrix<f32>, idx: usize| {
                            nalgebra::Vector3::new(
                                coords[(idx, 0)],
                                coords[(idx, 1)],
                                coords[(idx, 2)],
                            )
                        };
                        let (pcp, p1p, p2p, p3p) = (
                            mk(&cp, i),
                            mk(&cp, nbs[0].index()),
                            mk(&cp, nbs[1].index()),
                            mk(&cp, nbs[2].index()),
                        );
                        let v1p = p1p - pcp;
                        let v2p = p2p - pcp;
                        let v3p = p3p - pcp;
                        let volp = v1p.dot(&v2p.cross(&v3p));
                        ep += 10.0 * volp * volp;
                        let (pcm, p1m, p2m, p3m) = (
                            mk(&cm, i),
                            mk(&cm, nbs[0].index()),
                            mk(&cm, nbs[1].index()),
                            mk(&cm, nbs[2].index()),
                        );
                        let v1m = p1m - pcm;
                        let v2m = p2m - pcm;
                        let v3m = p3m - pcm;
                        let volm = v1m.dot(&v2m.cross(&v3m));
                        em += 10.0 * volm * volm;
                    }
                    let num_g = (ep - em) / (2.0 * eps);
                    let anal = g_oop[(atom, d)];
                    let axis = ["x", "y", "z"][d];
                    println!(
                        "    OOP atom={} {}: anal={:.4} num={:.4} diff={:.4}",
                        atom,
                        axis,
                        anal,
                        num_g,
                        (anal - num_g).abs()
                    );
                }
            }
        }

        // Energy breakdown AFTER FF
        let (e_b2, e_a2, e_bnd2, e_oop2, e_t_uff2, e_t_etkdg2) =
            energy_breakdown(&coords_after, &our_mol, &params, &smoothed);
        println!("AFTER FF:");
        println!("  Bonds:    {:.4} (delta {:.4})", e_b2, e_b2 - e_b);
        println!("  Angles:   {:.4} (delta {:.4})", e_a2, e_a2 - e_a);
        println!("  Bounds:   {:.4} (delta {:.4})", e_bnd2, e_bnd2 - e_bnd);
        println!("  OOP:      {:.4} (delta {:.4})", e_oop2, e_oop2 - e_oop);
        println!(
            "  Tors UFF: {:.4} (delta {:.4})",
            e_t_uff2,
            e_t_uff2 - e_t_uff
        );
        println!(
            "  Tors ETKDG: {:.4} (delta {:.4})",
            e_t_etkdg2,
            e_t_etkdg2 - e_t_etkdg
        );
        println!(
            "  TOTAL:    {:.4} (delta {:.4})",
            e_b2 + e_a2 + e_bnd2 + e_oop2 + e_t_uff2 + e_t_etkdg2,
            (e_b2 + e_a2 + e_bnd2 + e_oop2 + e_t_uff2 + e_t_etkdg2)
                - (e_b + e_a + e_bnd + e_oop + e_t_uff + e_t_etkdg)
        );

        // Coordinate changes
        let diff = &coords_after - &coords_before;
        let mut max_atom_move = 0.0f32;
        let mut total_move = 0.0f32;
        for i in 0..n {
            let d = (diff[(i, 0)] * diff[(i, 0)]
                + diff[(i, 1)] * diff[(i, 1)]
                + diff[(i, 2)] * diff[(i, 2)])
                .sqrt();
            max_atom_move = max_atom_move.max(d);
            total_move += d;
        }
        println!("Coord changes:");
        println!("  Max atom move: {:.6} Å", max_atom_move);
        println!("  Avg atom move: {:.6} Å", total_move / n as f32);
    }
}
