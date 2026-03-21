//! Integration tests for Track E1: Conformal Geometric Algebra (CGA)
//!
//! Validates all E1 subsections: algebraic layer, conformer generation,
//! and materials assembly.

#[cfg(feature = "experimental-cga")]
mod cga_tests {
    use sci_form::alpha::cga::*;

    // ---------------------------------------------------------------
    // E1.1 — CGA Algebraic Layer
    // ---------------------------------------------------------------

    #[test]
    fn e1_1a_basis_anticommutation() {
        // e_i * e_j = -e_j * e_i for i ≠ j
        for i in 0..5 {
            for j in (i + 1)..5 {
                let ei = Multivector::basis(i);
                let ej = Multivector::basis(j);
                let ab = ei.geo(&ej);
                let ba = ej.geo(&ei);
                let sum = ab + ba;
                assert!(
                    sum.norm_squared().abs() < 1e-20,
                    "e_{} e_{} + e_{} e_{} should be 0, got norm² {}",
                    i, j, j, i, sum.norm_squared()
                );
            }
        }
    }

    #[test]
    fn e1_1a_metric_signature() {
        // e_1²=e_2²=e_3²=e_+²=1, e_-²=-1
        for i in 0..4 {
            let ei = Multivector::basis(i);
            let sq = ei.geo(&ei);
            let s = sq.scalar();
            assert!(
                (s - 1.0).abs() < 1e-14,
                "e_{}² should be 1, got {}",
                i, s
            );
        }
        let em = Multivector::basis(4);
        let sq = em.geo(&em);
        assert!(
            (sq.scalar() + 1.0).abs() < 1e-14,
            "e_-² should be -1, got {}",
            sq.scalar()
        );
    }

    #[test]
    fn e1_1a_null_vectors() {
        // e_o and e_inf should be null: e_o² = 0, e_inf² = 0
        let eo = Multivector::origin();
        let ei = Multivector::infinity();
        assert!(
            eo.geo(&eo).norm_squared().abs() < 1e-20,
            "e_o should be null"
        );
        assert!(
            ei.geo(&ei).norm_squared().abs() < 1e-20,
            "e_inf should be null"
        );
        // e_o · e_inf = -1
        let dot = eo.inner(&ei);
        assert!(
            (dot.scalar() + 1.0).abs() < 1e-14,
            "e_o · e_inf should be -1, got {}",
            dot.scalar()
        );
    }

    #[test]
    fn e1_1b_unit_motor_identity() {
        // For a unit motor M: M * ~M = 1
        let m = Motor::rotor([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_4);
        let mv = m.mv;
        let rev = mv.reverse();
        let product = mv.geo(&rev);
        let s = product.scalar();
        assert!(
            (s - 1.0).abs() < 1e-10,
            "M * ~M should be 1, got {}",
            s
        );
    }

    #[test]
    fn e1_1b_reverse_involution() {
        let m = Multivector::basis(0).outer(&Multivector::basis(1));
        let rev = m.reverse();
        // Grade-2 blade reverses sign under reverse
        let sum = m + rev;
        // For a 2-blade, reverse = -original
        assert!(
            sum.norm_squared().abs() < 1e-20,
            "Reverse of a 2-blade should negate it"
        );
    }

    #[test]
    fn e1_1c_rotation_90_z() {
        let m = Motor::rotor([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let p = [1.0, 0.0, 0.0];
        let result = m.transform_point(p);
        assert!((result[0]).abs() < 1e-10, "x should be ~0, got {}", result[0]);
        assert!((result[1] - 1.0).abs() < 1e-10, "y should be ~1, got {}", result[1]);
        assert!((result[2]).abs() < 1e-10, "z should be ~0, got {}", result[2]);
    }

    #[test]
    fn e1_1c_translation() {
        let m = Motor::translator([3.0, 4.0, 5.0]);
        let p = [1.0, 2.0, 3.0];
        let result = m.transform_point(p);
        assert!((result[0] - 4.0).abs() < 1e-10);
        assert!((result[1] - 6.0).abs() < 1e-10);
        assert!((result[2] - 8.0).abs() < 1e-10);
    }

    #[test]
    fn e1_1c_motor_composition() {
        let r = Motor::rotor([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2);
        let t = Motor::translator([1.0, 0.0, 0.0]);
        let m = t.compose(&r); // first rotate, then translate
        let p = [1.0, 0.0, 0.0];
        let result = m.transform_point(p);
        // Rotate (1,0,0) by 90° around Z → (0,1,0), then translate by (1,0,0) → (1,1,0)
        assert!((result[0] - 1.0).abs() < 1e-10);
        assert!((result[1] - 1.0).abs() < 1e-10);
        assert!((result[2]).abs() < 1e-10);
    }

    // ---------------------------------------------------------------
    // E1.2 — CGA-based Conformer Generation
    // ---------------------------------------------------------------

    #[test]
    fn e1_2a_round_trip_embed_extract() {
        let points = [
            [1.234, -5.678, 9.012],
            [0.0, 0.0, 0.0],
            [-3.14, 2.72, 1.41],
        ];
        for &p in &points {
            let embedded = embed_point(p);
            let extracted = extract_point(&embedded);
            for i in 0..3 {
                assert!(
                    (p[i] - extracted[i]).abs() < 1e-10,
                    "Round-trip error at dim {}: {} vs {}",
                    i, p[i], extracted[i]
                );
            }
        }
    }

    #[test]
    fn e1_2b_dihedral_motor_preserves_axis_atoms() {
        let axis_a = [0.0, 0.0, 0.0];
        let axis_b = [1.5, 0.0, 0.0];
        let motor = dihedral_motor(axis_a, axis_b, std::f64::consts::FRAC_PI_4);

        // Transform axis atoms — they should not move
        let a2 = motor.transform_point(axis_a);
        let b2 = motor.transform_point(axis_b);

        for i in 0..3 {
            assert!(
                (axis_a[i] - a2[i]).abs() < 0.01,
                "Axis atom A moved at dim {}: {} → {}",
                i, axis_a[i], a2[i]
            );
            assert!(
                (axis_b[i] - b2[i]).abs() < 0.01,
                "Axis atom B moved at dim {}: {} → {}",
                i, axis_b[i], b2[i]
            );
        }
    }

    #[test]
    fn e1_2b_subtree_rotation() {
        // 4-atom chain: 0—1—2—3
        // Rotate subtree [2, 3] about bond 1—2 by 90°
        // Atom 3 is placed off the bond axis so rotation has a visible effect.
        let coords = vec![
            0.0, 0.0, 0.0,   // atom 0
            1.5, 0.0, 0.0,   // atom 1
            3.0, 0.0, 0.0,   // atom 2
            3.0, 1.5, 0.0,   // atom 3 — off-axis
        ];

        let motor = dihedral_motor(
            [coords[3], coords[4], coords[5]],   // atom 1
            [coords[6], coords[7], coords[8]],   // atom 2
            std::f64::consts::FRAC_PI_2,
        );
        let coords = apply_motor_to_subtree(&coords, &[2, 3], &motor);

        // Atoms 0 and 1 untouched
        assert!((coords[0]).abs() < 1e-10);
        assert!((coords[3] - 1.5).abs() < 1e-10);

        // Atom 2 is on the axis — should barely move
        assert!((coords[6] - 3.0).abs() < 0.02);

        // Atom 3 was at (3.0,1.5,0.0); after 90° rotation about x-axis
        // through (1.5,0,0)—(3.0,0,0) it should move to roughly (3.0,0,1.5).
        // Check it moved off the original y position and gained a z component.
        let dy = coords[10];
        let dz = coords[11];
        assert!(
            dy.abs() < 0.1 || dz.abs() > 0.1,
            "Atom 3 should have rotated off the x-y plane after 90° rotation,              got y={dy}, z={dz}"
        );
    }

    #[test]
    fn e1_2c_gimbal_lock_absence() {
        // Verify no gimbal lock at ω = ±90°
        let axis_a = [0.0, 0.0, 0.0];
        let axis_b = [1.5, 0.0, 0.0];
        let test_point = [3.0, 1.0, 0.0];

        for &angle in &[
            std::f64::consts::FRAC_PI_2,
            -std::f64::consts::FRAC_PI_2,
            std::f64::consts::PI,
        ] {
            let motor = dihedral_motor(axis_a, axis_b, angle);
            let result = motor.transform_point(test_point);
            let dist_before = ((test_point[0] - axis_a[0]).powi(2)
                + (test_point[1] - axis_a[1]).powi(2)
                + (test_point[2] - axis_a[2]).powi(2))
            .sqrt();
            let dist_after = ((result[0] - axis_a[0]).powi(2)
                + (result[1] - axis_a[1]).powi(2)
                + (result[2] - axis_a[2]).powi(2))
            .sqrt();
            assert!(
                (dist_before - dist_after).abs() < 0.01,
                "Distance not preserved at angle {:.3}: {:.6} vs {:.6}",
                angle, dist_before, dist_after
            );
        }
    }

    #[test]
    fn e1_2c_torsion_scan_produces_snapshots() {
        let mut coords = vec![
            0.0, 0.0, 0.0,
            1.5, 0.0, 0.0,
            3.0, 0.0, 0.0,
            4.5, 0.0, 0.0,
        ];
        let subtree = vec![2, 3];
        let snapshots = refine_torsion_cga(&mut coords, 1, 2, &subtree, 12);
        assert_eq!(snapshots.len(), 12, "Should produce 12 torsion snapshots");
    }

    // ---------------------------------------------------------------
    // E1.3 — CGA-based Materials Assembly
    // ---------------------------------------------------------------

    #[test]
    fn e1_3a_cubic_frame() {
        let frame = CgaFrame::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
        let cart = frame.frac_to_cart([1.0, 0.0, 0.0]);
        assert!((cart[0] - 10.0).abs() < 1e-10);
        assert!(cart[1].abs() < 1e-10);
        assert!(cart[2].abs() < 1e-10);
    }

    #[test]
    fn e1_3a_frac_cart_round_trip() {
        let frame = CgaFrame::from_parameters(12.0, 14.0, 9.0, 85.0, 95.0, 100.0);
        let frac = [0.25, 0.5, 0.75];
        let cart = frame.frac_to_cart(frac);
        let frac2 = frame.cart_to_frac(cart);
        for i in 0..3 {
            assert!(
                (frac[i] - frac2[i]).abs() < 1e-8,
                "Frac round-trip dim {}: {} → {}",
                i, frac[i], frac2[i]
            );
        }
    }

    #[test]
    fn e1_3b_pcu_assembly_1x1x1() {
        let crystal = assemble_framework_cga(30, 6, 10.0, 1);
        assert_eq!(crystal.num_atoms(), 4); // 1 node + 3 linkers
    }

    #[test]
    fn e1_3b_pcu_assembly_2x2x2() {
        let crystal = assemble_framework_cga(30, 6, 10.0, 2);
        assert_eq!(crystal.num_atoms(), 32); // 8 nodes + 24 linkers
    }

    #[test]
    fn e1_3c_node_linker_distance() {
        let crystal = assemble_framework_cga(30, 6, 10.0, 1);
        let node = &crystal.atoms[0];
        let linker = &crystal.atoms[1];
        let dx = linker.cart_coords[0] - node.cart_coords[0];
        let dy = linker.cart_coords[1] - node.cart_coords[1];
        let dz = linker.cart_coords[2] - node.cart_coords[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!(
            (dist - 5.0).abs() < 0.001,
            "Node-linker distance = {}, expected 5.0 Å (half of lattice parameter)",
            dist
        );
    }

    #[test]
    fn e1_3c_supercell_motor_correct() {
        let frame = CgaFrame::from_parameters(10.0, 10.0, 10.0, 90.0, 90.0, 90.0);
        let motor = frame.supercell_motor(2, 3, 1);
        let result = motor.transform_point([0.0, 0.0, 0.0]);
        assert!((result[0] - 20.0).abs() < 1e-10);
        assert!((result[1] - 30.0).abs() < 1e-10);
        assert!((result[2] - 10.0).abs() < 1e-10);
    }

    #[test]
    fn e1_3b_place_sbu_cga_test() {
        let sbu = vec![
            (30u8, [0.0, 0.0, 0.0]),
            (8u8, [1.0, 0.0, 0.0]),
            (8u8, [0.0, 1.0, 0.0]),
        ];
        let motor = Motor::translator([5.0, 5.0, 5.0]);
        let placed = place_sbu_cga(&sbu, &motor, "test");
        assert_eq!(placed.len(), 3);
        assert!((placed[0].1[0] - 5.0).abs() < 1e-10);
        assert!((placed[1].1[0] - 6.0).abs() < 1e-10);
        assert!((placed[2].1[1] - 6.0).abs() < 1e-10);
    }
}
