/// Trace the pipeline for a single molecule, dumping intermediate coordinates
/// to identify where divergence from RDKit occurs.
///
/// Usage:  cargo run --release --bin trace_pipeline -- "SMILES"
///
/// Outputs coordinates at each stage:
///   1. After embedding (initial coords)
///   2. After bounds FF
///   3. After 3D FF (final coords)
///
/// These can be compared with RDKit reference to identify the divergence point.

use nalgebra::DMatrix;
use sci_form::distgeom::{
    calculate_bounds_matrix_opts, triangle_smooth_tol,
    identify_chiral_sets, identify_tetrahedral_centers,
    pick_rdkit_distances, compute_initial_coords_rdkit,
    MinstdRand, check_tetrahedral_centers, check_chiral_centers,
    check_double_bond_geometry,
};
use sci_form::forcefield::bounds_ff::{minimize_bfgs_rdkit, ChiralSet};
use sci_form::forcefield::etkdg_3d::{
    build_etkdg_3d_ff_with_torsions, minimize_etkdg_3d_bfgs,
};
use sci_form::graph::Molecule;

const BASIN_THRESH: f32 = 5.0;
const FORCE_TOL: f32 = 1e-3;
const ERROR_TOL: f64 = 1e-5;
const MAX_MINIMIZED_E_PER_ATOM: f64 = 0.05;

fn compute_bounds_energy(
    coords: &DMatrix<f64>, bounds: &DMatrix<f64>,
    chiral_sets: &[ChiralSet],
    basin_thresh: f64, w4d: f64, wchiral: f64,
) -> f64 {
    let n = coords.nrows();
    let dim = coords.ncols();
    let mut energy = 0.0;
    for i in 1..n {
        for j in 0..i {
            let ub = bounds[(j, i)];
            let lb = bounds[(i, j)];
            if ub - lb > basin_thresh { continue; }
            let mut d2 = 0.0;
            for d in 0..dim {
                let diff = coords[(i, d)] - coords[(j, d)];
                d2 += diff * diff;
            }
            let ub2 = ub * ub;
            let lb2 = lb * lb;
            let val = if d2 > ub2 {
                d2 / ub2 - 1.0
            } else if d2 < lb2 {
                2.0 * lb2 / (lb2 + d2) - 1.0
            } else {
                0.0
            };
            if val > 0.0 { energy += val * val; }
        }
    }
    energy
}

fn pairwise_rmsd(a: &DMatrix<f64>, b: &DMatrix<f64>) -> f64 {
    let n = a.nrows();
    let mut sum = 0.0;
    let mut count = 0;
    for i in 0..n {
        for j in (i+1)..n {
            let da = ((a[(i,0)]-a[(j,0)]).powi(2) + (a[(i,1)]-a[(j,1)]).powi(2) + (a[(i,2)]-a[(j,2)]).powi(2)).sqrt();
            let db = ((b[(i,0)]-b[(j,0)]).powi(2) + (b[(i,1)]-b[(j,1)]).powi(2) + (b[(i,2)]-b[(j,2)]).powi(2)).sqrt();
            sum += (da - db).powi(2);
            count += 1;
        }
    }
    (sum / count as f64).sqrt()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let smiles = if args.len() > 1 { &args[1] } else { "N#Cc1ccc(Br)c(CC2CCCCC2)c1C=O" };
    
    eprintln!("=== PIPELINE TRACE: {} ===\n", smiles);
    
    let mol = Molecule::from_smiles(smiles).expect("Failed to parse SMILES");
    let n = mol.graph.node_count();
    eprintln!("Atoms: {}", n);
    
    // Load reference data
    let ref_data: Vec<serde_json::Value> = serde_json::from_str(
        &std::fs::read_to_string("tests/fixtures/gdb20_reference.json").expect("open ref")
    ).expect("parse ref");
    
    let ref_mol = ref_data.iter().find(|r| r["smiles"].as_str() == Some(smiles));
    let ref_coords = ref_mol.map(|r| {
        let atoms = r["atoms"].as_array().unwrap();
        let mut m = DMatrix::zeros(atoms.len(), 3);
        for (i, a) in atoms.iter().enumerate() {
            m[(i, 0)] = a["x"].as_f64().unwrap();
            m[(i, 1)] = a["y"].as_f64().unwrap();
            m[(i, 2)] = a["z"].as_f64().unwrap();
        }
        m
    });
    
    let csd_torsions: Vec<sci_form::forcefield::etkdg_3d::M6TorsionContrib> = ref_mol.map(|r| {
        r["torsions"].as_array().unwrap().iter().filter_map(|t| {
            let atoms = t["atoms"].as_array()?;
            let v_arr = t["v"].as_array()?;
            let s_arr = t["signs"].as_array()?;
            let mut v = [0.0f64; 6];
            let mut signs = [0.0f64; 6];
            for k in 0..6.min(v_arr.len()) { v[k] = v_arr[k].as_f64().unwrap_or(0.0); }
            for k in 0..6.min(s_arr.len()) { signs[k] = s_arr[k].as_f64().unwrap_or(0.0); }
            Some(sci_form::forcefield::etkdg_3d::M6TorsionContrib {
                i: atoms[0].as_u64()? as usize,
                j: atoms[1].as_u64()? as usize,
                k: atoms[2].as_u64()? as usize,
                l: atoms[3].as_u64()? as usize,
                signs, v,
            })
        }).collect()
    }).unwrap_or_default();
    
    // Build bounds matrix
    let bounds = {
        let raw = calculate_bounds_matrix_opts(&mol, true);
        let mut b = raw;
        if triangle_smooth_tol(&mut b, 0.0) { b }
        else {
            let raw2 = calculate_bounds_matrix_opts(&mol, false);
            let mut b2 = raw2.clone();
            if triangle_smooth_tol(&mut b2, 0.0) { b2 }
            else { let mut b3 = raw2; triangle_smooth_tol(&mut b3, 0.05); b3 }
        }
    };
    
    let chiral_sets = identify_chiral_sets(&mol);
    let tet_centers = identify_tetrahedral_centers(&mol);
    let use_4d = !chiral_sets.is_empty();
    let embed_dim = if use_4d { 4 } else { 3 };
    
    eprintln!("Chiral sets: {}, Use 4D: {}", chiral_sets.len(), use_4d);
    
    let mut rng = MinstdRand::new(42u32);
    let max_iters = 10 * n;
    
    for attempt in 0..max_iters {
        let dists = pick_rdkit_distances(&mut rng, &bounds);
        let coords_opt = compute_initial_coords_rdkit(&mut rng, &dists, embed_dim);
        let mut coords = match coords_opt {
            Some(c) => c,
            None => { continue; },
        };
        
        // === Stage 1: Initial coords ===
        let initial_coords = coords.columns(0, 3).into_owned();
        if let Some(ref rc) = ref_coords {
            let rmsd = pairwise_rmsd(&initial_coords, rc);
            eprintln!("  attempt {} | STAGE 1 (initial): RMSD vs ref = {:.6}", attempt, rmsd);
        }
        
        // === Stage 2: Bounds FF ===
        let init_energy = compute_bounds_energy(&coords, &bounds, &chiral_sets, 5.0, 0.1, 1.0);
        if init_energy > ERROR_TOL {
            let mut need_more = 1;
            let mut iters = 0;
            while need_more != 0 {
                need_more = minimize_bfgs_rdkit(
                    &mut coords, &bounds, &chiral_sets,
                    400, FORCE_TOL as f64, BASIN_THRESH, 0.1, 1.0,
                );
                iters += 1;
            }
            let post_energy = compute_bounds_energy(&coords, &bounds, &chiral_sets, 5.0, 0.1, 1.0);
            eprintln!("  attempt {} | STAGE 2 (boundsFF): energy {:.8} → {:.8} ({} rounds)", 
                attempt, init_energy, post_energy, iters);
        }
        
        let post_bounds_coords = coords.columns(0, 3).into_owned();
        if let Some(ref rc) = ref_coords {
            let rmsd = pairwise_rmsd(&post_bounds_coords, rc);
            eprintln!("  attempt {} | STAGE 2 (boundsFF): RMSD vs ref = {:.6}", attempt, rmsd);
        }
        
        // Energy check
        let energy = compute_bounds_energy(&coords, &bounds, &chiral_sets, 5.0, 0.1, 1.0);
        if energy / n as f64 >= MAX_MINIMIZED_E_PER_ATOM {
            eprintln!("  attempt {} → energy check FAILED ({:.6}/atom)", attempt, energy / n as f64);
            continue;
        }
        
        // Tet check
        if !check_tetrahedral_centers(&coords, &tet_centers) {
            eprintln!("  attempt {} → tet check FAILED", attempt);
            continue;
        }
        
        // Chiral check
        if !chiral_sets.is_empty() && !check_chiral_centers(&coords, &chiral_sets) {
            eprintln!("  attempt {} → chiral check FAILED", attempt);
            continue;
        }
        
        // Second bounds FF (4D only)
        if use_4d {
            let e2 = compute_bounds_energy(&coords, &bounds, &chiral_sets, 5.0, 1.0, 0.2);
            if e2 > ERROR_TOL {
                let mut need_more = 1;
                while need_more != 0 {
                    need_more = minimize_bfgs_rdkit(
                        &mut coords, &bounds, &chiral_sets,
                        200, FORCE_TOL as f64, BASIN_THRESH, 1.0, 0.2,
                    );
                }
            }
        }
        
        let coords3d = coords.columns(0, 3).into_owned();
        
        // === Stage 3: 3D FF ===
        let ff = build_etkdg_3d_ff_with_torsions(&mol, &coords3d, &bounds, &csd_torsions);
        
        let flat: Vec<f64> = (0..n).flat_map(|a| vec![coords3d[(a,0)], coords3d[(a,1)], coords3d[(a,2)]]).collect();
        let e3d_before = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&flat, n, &mol, &ff);
        eprintln!("  attempt {} | STAGE 3 (3dFF): energy before = {:.8}", attempt, e3d_before);
        
        let refined = if e3d_before > ERROR_TOL {
            minimize_etkdg_3d_bfgs(&mol, &coords3d, &ff, 300, FORCE_TOL)
        } else {
            coords3d.clone()
        };
        
        let flat_after: Vec<f64> = (0..n).flat_map(|a| vec![refined[(a,0)], refined[(a,1)], refined[(a,2)]]).collect();
        let e3d_after = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&flat_after, n, &mol, &ff);
        
        if let Some(ref rc) = ref_coords {
            let rmsd = pairwise_rmsd(&refined, rc);
            eprintln!("  attempt {} | STAGE 3 (3dFF): energy {:.8} → {:.8}, RMSD vs ref = {:.6}", 
                attempt, e3d_before, e3d_after, rmsd);
        }
        
        // Planarity check
        let n_improper = ff.inversion_contribs.len() / 3;
        let plan_energy = sci_form::forcefield::etkdg_3d::planarity_check_energy_f64(&flat_after, n, &ff);
        let plan_thresh = n_improper as f64 * 0.7;
        if plan_energy > plan_thresh {
            eprintln!("  attempt {} → planarity FAILED ({:.4} > {:.4})", attempt, plan_energy, plan_thresh);
            continue;
        }
        
        // Double bond check
        if !check_double_bond_geometry(&mol, &refined) {
            eprintln!("  attempt {} → double bond FAILED", attempt);
            continue;
        }
        
        eprintln!("  attempt {} → SUCCESS", attempt);
        
        // Print final coordinate comparison
        if let Some(ref rc) = ref_coords {
            let final_rmsd = pairwise_rmsd(&refined, rc);
            eprintln!("\n=== FINAL RESULT ===");
            eprintln!("Final RMSD vs RDKit: {:.6} Å", final_rmsd);
            eprintln!("Bounds FF RMSD vs RDKit: {:.6} Å", pairwise_rmsd(&post_bounds_coords, rc));
            
            // Also compute RMSD at initial coords 
            eprintln!("Initial coords RMSD vs RDKit: {:.6} Å", pairwise_rmsd(&initial_coords, rc));
            
            // Evaluate 3D FF energy at RDKit coords
            let ref_flat: Vec<f64> = (0..n).flat_map(|a| vec![rc[(a,0)], rc[(a,1)], rc[(a,2)]]).collect();
            let e_at_rdkit = sci_form::forcefield::etkdg_3d::etkdg_3d_energy_f64(&ref_flat, n, &mol, &ff);
            eprintln!("\n3D FF energy at OUR coords: {:.8}", e3d_after);
            eprintln!("3D FF energy at RDKIT coords: {:.8}", e_at_rdkit);
            eprintln!("Difference: {:.8} (negative = ours lower = better?)", e3d_after - e_at_rdkit);
        }
        
        return;
    }
    
    eprintln!("All attempts failed!");
}
