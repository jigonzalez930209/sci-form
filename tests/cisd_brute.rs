//! Cross-check of the CISD implementation against an explicit
//! second-quantized reference.
//!
//! The reference builds the Hamiltonian by applying the fermionic operators
//! `a†_p a_q` and `a†_p a†_q a_s a_r` to each determinant (bitmask
//! representation, canonical signs) — no Slater–Condon formulas — and is
//! itself validated against PySCF `ci.ucisd` on H₂O/STO-3G to 1.5e-12 Ha.
//! The library matrix must reproduce it element by element (the synthetic
//! 8-orbital/8-electron case below is sensitive to the creation-operator
//! ordering that the Slater–Condon formulas must follow).
use nalgebra::{DMatrix, SymmetricEigen};

/// Spin-orbital index convention: `idx = spin * n + spatial` (spin 0 = alpha).
fn brute_force_correlation(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    space: &[(u64, u64)],
) -> f64 {
    let h = brute_force_matrix(
        orbital_energies,
        coefficients,
        eris,
        n_basis,
        n_occupied,
        space,
    );
    let eigen = SymmetricEigen::<f64, nalgebra::Dyn>::new(h);
    let mut lowest = f64::INFINITY;
    for e in eigen.eigenvalues.iter() {
        lowest = lowest.min(*e);
    }
    let so = |p: usize, q: usize| -> f64 {
        spin_orbital_h(
            orbital_energies,
            coefficients,
            eris,
            n_basis,
            n_occupied,
            p,
            q,
        )
    };
    let se = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        spin_orbital_eri(coefficients, eris, n_basis, p, q, r, s)
    };
    let ref_det = space[0];
    let mut e_ref = 0.0;
    let occ: Vec<usize> = (0..2 * n_basis)
        .filter(|i| {
            let b = 1u64 << (i % n_basis);
            if i / n_basis == 0 {
                ref_det.0 & b != 0
            } else {
                ref_det.1 & b != 0
            }
        })
        .collect();
    for &p in &occ {
        e_ref += so(p, p);
        for &q in &occ {
            e_ref += 0.5 * (se(p, p, q, q) - se(p, q, q, p));
        }
    }
    lowest - e_ref
}

fn spin_orbital_h(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    p: usize,
    q: usize,
) -> f64 {
    use sci_form::hf::mo_transform::ao_to_mo_transform;
    if p / n_basis != q / n_basis {
        return 0.0;
    }
    let mo = ao_to_mo_transform(coefficients, eris, n_basis, n_basis);
    let sp = p % n_basis;
    let sq = q % n_basis;
    let mut fock = 0.0;
    for j in 0..n_occupied {
        fock += 2.0 * mo.get(sp, sq, j, j) - mo.get(sp, j, j, sq);
    }
    let h = if p == q { orbital_energies[sp] } else { 0.0 };
    h - fock
}

fn spin_orbital_eri(
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    p: usize,
    q: usize,
    r: usize,
    s: usize,
) -> f64 {
    use sci_form::hf::mo_transform::ao_to_mo_transform;
    if p / n_basis != q / n_basis || r / n_basis != s / n_basis {
        return 0.0;
    }
    let mo = ao_to_mo_transform(coefficients, eris, n_basis, n_basis);
    mo.get(p % n_basis, q % n_basis, r % n_basis, s % n_basis)
}

fn brute_force_matrix(
    orbital_energies: &[f64],
    coefficients: &DMatrix<f64>,
    eris: &[f64],
    n_basis: usize,
    n_occupied: usize,
    space: &[(u64, u64)],
) -> DMatrix<f64> {
    use sci_form::hf::mo_transform::ao_to_mo_transform;
    let mo = ao_to_mo_transform(coefficients, eris, n_basis, n_basis);
    let n_spin = 2 * n_basis;
    let spin = |i: usize| i / n_basis;
    let spatial = |i: usize| i % n_basis;
    let so = |p: usize, q: usize| -> f64 {
        if spin(p) == spin(q) {
            let h = if p == q {
                orbital_energies[spatial(p)]
            } else {
                0.0
            };
            // h_mo via the canonical Fock condition
            let mut fock = 0.0;
            for j in 0..n_occupied {
                fock += 2.0 * mo.get(spatial(p), spatial(q), j, j)
                    - mo.get(spatial(p), j, j, spatial(q));
            }
            h - fock
        } else {
            0.0
        }
    };
    let se = |p: usize, q: usize, r: usize, s: usize| -> f64 {
        if spin(p) == spin(q) && spin(r) == spin(s) {
            mo.get(spatial(p), spatial(q), spatial(r), spatial(s))
        } else {
            0.0
        }
    };

    let nd = space.len();
    let mut h = DMatrix::zeros(nd, nd);
    let mask = |det: &(u64, u64)| -> (u64, u64) { *det };

    for (col, det) in space.iter().enumerate() {
        // one-electron: sum_pq h_pq a+_p a_q
        for p in 0..n_spin {
            for q in 0..n_spin {
                let val = so(p, q);
                if val.abs() < 1e-14 {
                    continue;
                }
                if let Some(st1) = annihilate(det, q, n_basis) {
                    if let Some(st2) = create(&st1, p, n_basis) {
                        if let Some(row) = space.iter().position(|d| mask(d) == st2) {
                            let sign =
                                phase_annihilate(det, q, n_basis) * phase_create(&st1, p, n_basis);
                            h[(row, col)] += sign * val;
                        }
                    }
                }
            }
        }
        // two-electron: 1/2 sum (pr|qs) a+_p a+_q a_s a_r
        for p in 0..n_spin {
            for q in 0..n_spin {
                for r in 0..n_spin {
                    for s in 0..n_spin {
                        let val = se(p, r, q, s);
                        if val.abs() < 1e-14 {
                            continue;
                        }
                        let Some(st1) = annihilate(det, r, n_basis) else {
                            continue;
                        };
                        let s1 = phase_annihilate(det, r, n_basis);
                        let Some(st2) = annihilate(&st1, s, n_basis) else {
                            continue;
                        };
                        let s2 = phase_annihilate(&st1, s, n_basis);
                        let Some(st3) = create(&st2, q, n_basis) else {
                            continue;
                        };
                        let s3 = phase_create(&st2, q, n_basis);
                        let Some(st4) = create(&st3, p, n_basis) else {
                            continue;
                        };
                        let s4 = phase_create(&st3, p, n_basis);
                        if let Some(row) = space.iter().position(|d| mask(d) == st4) {
                            h[(row, col)] += 0.5 * s1 * s2 * s3 * s4 * val;
                        }
                    }
                }
            }
        }
    }
    h
}

fn annihilate(det: &(u64, u64), p: usize, n: usize) -> Option<(u64, u64)> {
    let b = 1u64 << (p % n);
    let (m, other) = if p / n == 0 {
        (det.0, det.1)
    } else {
        (det.1, det.0)
    };
    if m & b == 0 {
        return None;
    }
    let new = m & !b;
    Some(if p / n == 0 {
        (new, other)
    } else {
        (other, new)
    })
}

fn create(det: &(u64, u64), p: usize, n: usize) -> Option<(u64, u64)> {
    let b = 1u64 << (p % n);
    let (m, other) = if p / n == 0 {
        (det.0, det.1)
    } else {
        (det.1, det.0)
    };
    if m & b != 0 {
        return None;
    }
    let new = m | b;
    Some(if p / n == 0 {
        (new, other)
    } else {
        (other, new)
    })
}

fn phase_annihilate(det: &(u64, u64), p: usize, n: usize) -> f64 {
    let m = if p / n == 0 { det.0 } else { det.1 };
    let below = (m & ((1u64 << (p % n)) - 1)).count_ones();
    if below.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    }
}

fn phase_create(det: &(u64, u64), p: usize, n: usize) -> f64 {
    // Canonical ordering: a creation operator anticommutes past the occupied
    // orbitals *below* its index.
    let m = if p / n == 0 { det.0 } else { det.1 };
    let below = (m & ((1u64 << (p % n)) - 1)).count_ones();
    if below.is_multiple_of(2) {
        1.0
    } else {
        -1.0
    }
}

/// Rebuild the CISD space exactly like the library does.
fn build_space(n_basis: usize, n_occupied: usize) -> Vec<(u64, u64)> {
    let n_virtual = n_basis - n_occupied;
    let ref_mask = (1u64 << n_occupied) - 1;
    let reference = (ref_mask, ref_mask);
    let mut space = vec![reference];
    let mut occ: Vec<usize> = (0..n_occupied).collect();
    occ.extend((0..n_occupied).map(|k| n_basis + k));
    let mut virt: Vec<usize> = (0..2 * n_virtual)
        .map(|k| {
            let spin = k / n_virtual;
            let spatial = k % n_virtual + n_occupied;
            spin * n_basis + spatial
        })
        .collect();
    virt.sort_unstable();
    let spin_of = |i: usize| i / n_basis;
    let set = |det: &mut (u64, u64), idx: usize, value: bool| {
        let b = 1u64 << (idx % n_basis);
        if idx / n_basis == 0 {
            if value {
                det.0 |= b
            } else {
                det.0 &= !b
            }
        } else if value {
            det.1 |= b
        } else {
            det.1 &= !b
        }
    };
    for &i in &occ {
        for &a in &virt {
            if spin_of(i) != spin_of(a) {
                continue;
            }
            let mut d = reference;
            set(&mut d, i, false);
            set(&mut d, a, true);
            space.push(d);
        }
    }
    for (ii, &i) in occ.iter().enumerate() {
        for &j in occ.iter().skip(ii + 1) {
            let na = (spin_of(i) == 0) as usize + (spin_of(j) == 0) as usize;
            for (aa, &a) in virt.iter().enumerate() {
                for &b in virt.iter().skip(aa + 1) {
                    let nb = (spin_of(a) == 0) as usize + (spin_of(b) == 0) as usize;
                    if na != nb {
                        continue;
                    }
                    let mut d = reference;
                    set(&mut d, i, false);
                    set(&mut d, j, false);
                    set(&mut d, a, true);
                    set(&mut d, b, true);
                    space.push(d);
                }
            }
        }
    }
    space
}

#[test]
fn brute_force_vs_library() {
    // Synthetic 8-orbital / 8-electron system with deterministic "integrals".
    let n = 8usize;
    let nocc = 4usize;
    let n_pair = n * (n + 1) / 2;
    let mut eris = vec![0.0; n_pair * (n_pair + 1) / 2];
    let mut seed = 12345u64;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 33) as f64 / (1u64 << 31) as f64) - 0.5
    };
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    let v = next() * 0.3;
                    eris[sci_form::hf::integrals::eri_index(i, j, k, l, n)] = v;
                }
            }
        }
    }
    let orbital_energies: Vec<f64> = (0..n).map(|i| -0.5 + 0.2 * i as f64).collect();
    let coefficients = DMatrix::identity(n, n);
    let space = build_space(n, nocc);
    let e_brute = brute_force_correlation(&orbital_energies, &coefficients, &eris, n, nocc, &space);
    let res = sci_form::hf::cisd::compute_cisd(&orbital_energies, &coefficients, &eris, n, nocc, 3);
    println!(
        "n_dets={} library={:.10} brute={:.10} Δ={:.2e}",
        space.len(),
        res.correlation_energy,
        e_brute,
        res.correlation_energy - e_brute
    );
    assert!(
        (res.correlation_energy - e_brute).abs() < 1e-8,
        "library {:.10} vs brute {:.10}",
        res.correlation_energy,
        e_brute
    );
}

#[test]
fn diff_matrices_synthetic() {
    let n = 8usize;
    let nocc = 4usize;
    let n_pair = n * (n + 1) / 2;
    let mut eris = vec![0.0; n_pair * (n_pair + 1) / 2];
    let mut seed = 12345u64;
    let mut next = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 33) as f64 / (1u64 << 31) as f64) - 0.5
    };
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                for l in 0..n {
                    eris[sci_form::hf::integrals::eri_index(i, j, k, l, n)] = next() * 0.3;
                }
            }
        }
    }
    let orbital_energies: Vec<f64> = (0..n).map(|i| -0.5 + 0.2 * i as f64).collect();
    let coefficients = DMatrix::identity(n, n);
    let h_lib = sci_form::hf::cisd::debug_matrix(&orbital_energies, &coefficients, &eris, n, nocc);
    let space = sci_form::hf::cisd::debug_space(n, nocc);
    let h_brute = brute_force_matrix(&orbital_energies, &coefficients, &eris, n, nocc, &space);
    let mut worst = 0.0f64;
    let mut worst_pair = None;
    let mut n_bad = 0;
    for row in 0..space.len() {
        for col in 0..=row {
            let d = (h_lib[(row, col)] - h_brute[(row, col)]).abs();
            if d > 1e-9 {
                n_bad += 1;
            }
            if d > worst {
                worst = d;
                worst_pair = Some((
                    row,
                    col,
                    space[row],
                    space[col],
                    h_lib[(row, col)],
                    h_brute[(row, col)],
                ));
            }
        }
    }
    println!("pares con |Δ|>1e-9: {n_bad}  peor = {worst:.3e}");
    if let Some((r, c, d1, d2, a, b)) = worst_pair {
        println!(
            "peor par: row={r} col={c} d1=({:#b},{:#b}) d2=({:#b},{:#b}) lib={a:+.8} brute={b:+.8}",
            d1.0, d1.1, d2.0, d2.1
        );
    }
    assert_eq!(n_bad, 0, "matrices differ");
}
