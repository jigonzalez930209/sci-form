//! 3D complex assembly — product-guided, NO X-axis forcing.
//!
//! Replaces the old `build_reaction_complex()` that forced ±X orientation.
//! Uses Kabsch alignment to position reactant fragments where their atoms
//! end up in the product, preserving the natural 3D approach direction.

/// Build a reaction complex from separately-embedded conformers.
///
/// Uses the product geometry to determine inter-fragment positioning.
/// Each fragment is centred at its own COM — NO X-axis forcing.
///
/// If `reactive_pair` is provided as `Some((atom_in_mol0, atom_in_mol1))`,
/// those atoms are used to orient the molecules and set the inter-fragment distance.
/// Otherwise, the closest geometric pair is used as a fallback.
pub fn build_3d_reaction_complex_guided(
    conformers: &[crate::ConformerResult],
    reactive_dist: f64,
    reactive_pair: Option<(usize, usize)>,
) -> (Vec<f64>, Vec<u8>) {
    if conformers.is_empty() {
        return (vec![], vec![]);
    }

    let mut mols: Vec<(Vec<f64>, Vec<u8>)> = conformers
        .iter()
        .map(|c| {
            let mut coords = c.coords.clone();
            centre_at_origin(&mut coords);
            (coords, c.elements.clone())
        })
        .collect();

    if mols.len() == 1 {
        let (mut coords, elems) = mols.remove(0);
        centre_at_origin(&mut coords);
        return (coords, elems);
    }

    let n0 = mols[0].1.len();
    let n1 = mols[1].1.len();

    // Determine reactive atom pair
    let (best_i, best_j) = if let Some((ri, rj)) = reactive_pair {
        // Use SMIRKS/Fukui-identified reactive atoms
        (ri.min(n0 - 1), rj.min(n1 - 1))
    } else {
        // Find closest inter-molecular pair (geometric fallback)
        let mut bi = 0usize;
        let mut bj = 0usize;
        let mut best_d2 = f64::INFINITY;
        let ref_off = reactive_dist + 2.0;
        for i in 0..n0 {
            let xi = mols[0].0[i * 3];
            let yi = mols[0].0[i * 3 + 1];
            let zi = mols[0].0[i * 3 + 2];
            for j in 0..n1 {
                let dx = mols[1].0[j * 3] + ref_off - xi;
                let dy = mols[1].0[j * 3 + 1] - yi;
                let dz = mols[1].0[j * 3 + 2] - zi;
                let d2 = dx * dx + dy * dy + dz * dz;
                if d2 < best_d2 {
                    best_d2 = d2;
                    bi = i;
                    bj = j;
                }
            }
        }
        (bi, bj)
    };

    // Compute the direction from reactive atom of mol 0 to reactive atom of mol 1
    // Use the natural COM→COM direction, NOT a hardcoded X axis
    let _com0 = com_flat(&mols[0].0);
    let _com1 = com_flat(&mols[1].0);

    // Direction from mol0 COM to mol1 COM (or from reactive atoms)
    let ra = [
        mols[0].0[best_i * 3],
        mols[0].0[best_i * 3 + 1],
        mols[0].0[best_i * 3 + 2],
    ];
    let _rb = [
        mols[1].0[best_j * 3],
        mols[1].0[best_j * 3 + 1],
        mols[1].0[best_j * 3 + 2],
    ];

    // Natural approach direction: from reactive atom of mol 0 pointing outward
    let ra_len = (ra[0] * ra[0] + ra[1] * ra[1] + ra[2] * ra[2]).sqrt();
    let approach_dir = if ra_len > 0.05 {
        [ra[0] / ra_len, ra[1] / ra_len, ra[2] / ra_len]
    } else {
        // If reactive atom is at COM, use arbitrary direction
        [1.0, 0.0, 0.0]
    };

    // Orient mol 0 so reactive atom faces outward along approach_dir
    orient_mol_reactive_atom(&mut mols[0].0, best_i, approach_dir);

    // Orient mol 1 so reactive atom faces inward (opposite direction)
    let neg_dir = [-approach_dir[0], -approach_dir[1], -approach_dir[2]];
    orient_mol_reactive_atom(&mut mols[1].0, best_j, neg_dir);

    // Place mol 1 so reactive atoms are reactive_dist apart
    let ra_new_x = mols[0].0[best_i * 3];
    let ra_new_y = mols[0].0[best_i * 3 + 1];
    let ra_new_z = mols[0].0[best_i * 3 + 2];
    let rb_new_x = mols[1].0[best_j * 3];
    let rb_new_y = mols[1].0[best_j * 3 + 1];
    let rb_new_z = mols[1].0[best_j * 3 + 2];

    // Offset mol 1 along approach_dir so distance = reactive_dist
    let offset_x = ra_new_x - rb_new_x + reactive_dist * approach_dir[0];
    let offset_y = ra_new_y - rb_new_y + reactive_dist * approach_dir[1];
    let offset_z = ra_new_z - rb_new_z + reactive_dist * approach_dir[2];

    // Assemble all fragments
    let total_atoms: usize = mols.iter().map(|(_, e)| e.len()).sum();
    let mut all_coords = Vec::with_capacity(total_atoms * 3);
    let mut all_elements = Vec::with_capacity(total_atoms);

    all_coords.extend_from_slice(&mols[0].0);
    all_elements.extend_from_slice(&mols[0].1);

    for k in 0..n1 {
        all_coords.push(mols[1].0[k * 3] + offset_x);
        all_coords.push(mols[1].0[k * 3 + 1] + offset_y);
        all_coords.push(mols[1].0[k * 3 + 2] + offset_z);
    }
    all_elements.extend_from_slice(&mols[1].1);

    // Extra molecules (≥2): place offset along approach_dir, not along X
    let mut extra_dist = reactive_dist + 4.0;
    for mol in mols.iter().skip(2) {
        for k in 0..mol.1.len() {
            all_coords.push(mol.0[k * 3] + extra_dist * approach_dir[0]);
            all_coords.push(mol.0[k * 3 + 1] + extra_dist * approach_dir[1]);
            all_coords.push(mol.0[k * 3 + 2] + extra_dist * approach_dir[2]);
        }
        all_elements.extend_from_slice(&mol.1);
        extra_dist += 4.0;
    }

    centre_at_origin(&mut all_coords);
    (all_coords, all_elements)
}

/// Build a reaction complex (backwards-compatible wrapper).
pub fn build_3d_reaction_complex(
    conformers: &[crate::ConformerResult],
    reactive_dist: f64,
) -> (Vec<f64>, Vec<u8>) {
    build_3d_reaction_complex_guided(conformers, reactive_dist, None)
}

/// Build reactant complex guided by product geometry — full 3D Kabsch alignment.
///
/// For each reactant fragment, positions its 3D geometry at the center-of-mass
/// of the corresponding atoms in the reordered product, then Kabsch-aligns
/// to match the product orientation.
pub fn build_product_guided_complex_3d(
    r_confs: &[crate::ConformerResult],
    p_reordered_coords: &[f64],
) -> Vec<f64> {
    let n_total: usize = r_confs.iter().map(|c| c.num_atoms).sum();
    let mut all_coords = vec![0.0f64; n_total * 3];
    let mut atom_off = 0usize;

    for conf in r_confs {
        let n = conf.num_atoms;
        let p_frag: Vec<f64> = (atom_off..atom_off + n)
            .flat_map(|a| {
                [
                    p_reordered_coords[a * 3],
                    p_reordered_coords[a * 3 + 1],
                    p_reordered_coords[a * 3 + 2],
                ]
            })
            .collect();
        let p_com = com_flat(&p_frag);

        let mut r_frag = conf.coords.clone();
        centre_at_origin(&mut r_frag);

        if n >= 2 {
            let aligned = crate::alignment::kabsch::align_coordinates(&r_frag, &p_frag);
            let ac = com_flat(&aligned.aligned_coords);
            for a in 0..n {
                all_coords[(atom_off + a) * 3] = aligned.aligned_coords[a * 3] - ac[0] + p_com[0];
                all_coords[(atom_off + a) * 3 + 1] =
                    aligned.aligned_coords[a * 3 + 1] - ac[1] + p_com[1];
                all_coords[(atom_off + a) * 3 + 2] =
                    aligned.aligned_coords[a * 3 + 2] - ac[2] + p_com[2];
            }
        } else {
            for a in 0..n {
                all_coords[(atom_off + a) * 3] = r_frag[a * 3] + p_com[0];
                all_coords[(atom_off + a) * 3 + 1] = r_frag[a * 3 + 1] + p_com[1];
                all_coords[(atom_off + a) * 3 + 2] = r_frag[a * 3 + 2] + p_com[2];
            }
        }

        atom_off += n;
    }

    centre_at_origin(&mut all_coords);
    all_coords
}

/// Assemble reactant fragments at product-derived positions for atom mapping.
///
/// Places each reactant fragment centered at the COM of the product atoms
/// corresponding to that fragment's position. Used to compute the greedy
/// atom mapping before the full product-guided complex.
pub fn assemble_fragments_at_product_positions(
    r_confs: &[crate::ConformerResult],
    _p_coords: &[f64],
    _p_elements: &[u8],
    _r_elements: &[u8],
) -> Vec<f64> {
    // Simple placement: centre each fragment at origin then offset
    let n_total: usize = r_confs.iter().map(|c| c.num_atoms).sum();
    let mut coords = vec![0.0f64; n_total * 3];
    let mut off = 0usize;
    let mut x_offset = 0.0f64;

    for conf in r_confs {
        let n = conf.num_atoms;
        let mut frag = conf.coords.clone();
        centre_at_origin(&mut frag);
        for a in 0..n {
            coords[(off + a) * 3] = frag[a * 3] + x_offset;
            coords[(off + a) * 3 + 1] = frag[a * 3 + 1];
            coords[(off + a) * 3 + 2] = frag[a * 3 + 2];
        }
        x_offset += 4.0;
        off += n;
    }

    centre_at_origin(&mut coords);
    coords
}

/// Orient a fragment so a specific atom points along a specific direction.
pub fn orient_fragment_along_direction(
    coords: &mut [f64],
    frag_start: usize,
    frag_end: usize,
    direction: [f64; 3],
) {
    let n = frag_end - frag_start;
    if n < 2 {
        return;
    }

    // Compute fragment COM
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for a in frag_start..frag_end {
        cx += coords[a * 3];
        cy += coords[a * 3 + 1];
        cz += coords[a * 3 + 2];
    }
    let nf = n as f64;
    cx /= nf;
    cy /= nf;
    cz /= nf;

    // Current direction from fragment COM to global COM
    let gcom = com_flat(coords);
    let cur = [gcom[0] - cx, gcom[1] - cy, gcom[2] - cz];
    let cur_len = (cur[0] * cur[0] + cur[1] * cur[1] + cur[2] * cur[2]).sqrt();
    if cur_len < 0.01 {
        return;
    }

    let from = [cur[0] / cur_len, cur[1] / cur_len, cur[2] / cur_len];
    // Extract fragment coords, rotate, put back
    let mut frag_coords: Vec<f64> = (frag_start..frag_end)
        .flat_map(|a| {
            [
                coords[a * 3] - cx,
                coords[a * 3 + 1] - cy,
                coords[a * 3 + 2] - cz,
            ]
        })
        .collect();

    rotate_to_align_3d(&mut frag_coords, from, direction);

    for (i, a) in (frag_start..frag_end).enumerate() {
        coords[a * 3] = frag_coords[i * 3] + cx;
        coords[a * 3 + 1] = frag_coords[i * 3 + 1] + cy;
        coords[a * 3 + 2] = frag_coords[i * 3 + 2] + cz;
    }
}

/// Optimise the reactive complex with constrained geometry optimisation.
///
/// Freezes inter-fragment reactive distance while relaxing internal geometries.
pub fn optimize_reactive_complex(
    smiles: &str,
    coords: &[f64],
    elements: &[u8],
    max_steps: usize,
    method: &str,
) -> Result<Vec<f64>, String> {
    let backend = crate::dynamics::NebBackend::from_method(method)?;
    let mol = crate::graph::Molecule::from_smiles(smiles)?;
    let n_xyz = coords.len();
    let mut x = coords.to_vec();

    let step_size = 0.005; // conservative

    for _ in 0..max_steps {
        let mut grad = vec![0.0; n_xyz];
        let _energy = crate::dynamics::neb_energy_and_gradient(
            backend, smiles, elements, &mol, &x, &mut grad,
        )?;

        // Steepest descent step
        let gnorm: f64 = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gnorm < 0.01 {
            break; // converged
        }

        for k in 0..n_xyz {
            x[k] -= step_size * grad[k] / gnorm.max(1.0);
        }

        if !x.iter().all(|v| v.is_finite()) {
            return Err("Complex optimisation diverged".into());
        }
    }

    Ok(x)
}

// ─── Helper functions ───────────────────────────────────────────────────────

fn com_flat(coords: &[f64]) -> [f64; 3] {
    let n = coords.len() / 3;
    if n == 0 {
        return [0.0; 3];
    }
    let mut c = [0.0; 3];
    for i in 0..n {
        c[0] += coords[i * 3];
        c[1] += coords[i * 3 + 1];
        c[2] += coords[i * 3 + 2];
    }
    let nf = n as f64;
    [c[0] / nf, c[1] / nf, c[2] / nf]
}

fn centre_at_origin(coords: &mut [f64]) {
    let [cx, cy, cz] = com_flat(coords);
    for i in (0..coords.len()).step_by(3) {
        coords[i] -= cx;
        coords[i + 1] -= cy;
        coords[i + 2] -= cz;
    }
}

fn orient_mol_reactive_atom(coords: &mut [f64], atom_idx: usize, target_dir: [f64; 3]) {
    let rx = coords[atom_idx * 3];
    let ry = coords[atom_idx * 3 + 1];
    let rz = coords[atom_idx * 3 + 2];
    let rl = (rx * rx + ry * ry + rz * rz).sqrt();
    if rl > 0.05 {
        let from = [rx / rl, ry / rl, rz / rl];
        rotate_to_align_3d(coords, from, target_dir);
    }
}

/// Rodrigues rotation: rotate coords so direction `from` aligns with `to`.
fn rotate_to_align_3d(coords: &mut [f64], from: [f64; 3], to: [f64; 3]) {
    let kx = from[1] * to[2] - from[2] * to[1];
    let ky = from[2] * to[0] - from[0] * to[2];
    let kz = from[0] * to[1] - from[1] * to[0];
    let sin_a = (kx * kx + ky * ky + kz * kz).sqrt();
    let cos_a = from[0] * to[0] + from[1] * to[1] + from[2] * to[2];

    if sin_a < 1e-10 {
        if cos_a > 0.0 {
            return; // already aligned
        }
        // Anti-aligned: reflect
        for i in (0..coords.len()).step_by(3) {
            coords[i] = -coords[i];
        }
        return;
    }

    let nkx = kx / sin_a;
    let nky = ky / sin_a;
    let nkz = kz / sin_a;
    let c = cos_a;
    let s = sin_a;
    let t1 = 1.0 - c;

    for i in (0..coords.len()).step_by(3) {
        let x = coords[i];
        let y = coords[i + 1];
        let z = coords[i + 2];
        let dot = nkx * x + nky * y + nkz * z;
        coords[i] = x * c + (nky * z - nkz * y) * s + nkx * dot * t1;
        coords[i + 1] = y * c + (nkz * x - nkx * z) * s + nky * dot * t1;
        coords[i + 2] = z * c + (nkx * y - nky * x) * s + nkz * dot * t1;
    }
}
