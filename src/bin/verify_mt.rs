use sci_form::distgeom::Mt19937;

fn main() {
    let mut mt = Mt19937::new(42);
    let expected_u32: [u32; 10] = [
        1608637542, 3421126067, 4083286876, 787846414, 3143890026, 3348747335, 2571218620,
        2563451924, 670094950, 1914837113,
    ];

    println!("Verifying MT19937 u32 output (seed=42):");
    let mut all_ok = true;
    for i in 0..10 {
        let v = mt.next_u32();
        let ok = v == expected_u32[i];
        println!(
            "  [{}] got={}, expected={}, {}",
            i,
            v,
            expected_u32[i],
            if ok { "OK" } else { "FAIL" }
        );
        if !ok {
            all_ok = false;
        }
    }

    let mut mt2 = Mt19937::new(42);
    let expected_double: [f64; 5] = [
        3.74540114309638739e-01,
        7.96542984200641513e-01,
        9.50714311562478542e-01,
        1.83434787672013044e-01,
        7.31993938330560923e-01,
    ];

    println!("\nVerifying MT19937 double output (seed=42):");
    for i in 0..5 {
        let d = mt2.next_double();
        let diff = (d - expected_double[i]).abs();
        let ok = diff < 1e-17;
        println!(
            "  [{}] got={:.17e}, expected={:.17e}, diff={:.2e} {}",
            i,
            d,
            expected_double[i],
            diff,
            if ok { "OK" } else { "FAIL" }
        );
        if !ok {
            all_ok = false;
        }
    }

    if all_ok {
        println!("\nAll MT19937 values match! ✓");
    } else {
        println!("\nSome values don't match! ✗");
        std::process::exit(1);
    }
}
