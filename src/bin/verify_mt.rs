use sci_form::distgeom::Mt19937;

fn main() {
    let mut mt = Mt19937::new(42);
    let expected_u32: [u32; 10] = [
        1608637542, 3421126067, 4083286876, 787846414, 3143890026, 3348747335, 2571218620,
        2563451924, 670094950, 1914837113,
    ];

    println!("Verifying MT19937 u32 output (seed=42):");
    let mut all_ok = true;
    for (i, &expected) in expected_u32.iter().enumerate() {
        let v = mt.next_u32();
        let ok = v == expected;
        println!(
            "  [{}] got={}, expected={}, {}",
            i,
            v,
            expected,
            if ok { "OK" } else { "FAIL" }
        );
        if !ok {
            all_ok = false;
        }
    }

    let mut mt2 = Mt19937::new(42);
    let expected_double: [f64; 5] = [
        3.745_401_143e-1,
        7.965_429_842e-1,
        9.507_143_116e-1,
        1.834_347_877e-1,
        7.319_939_383e-1,
    ];

    println!("\nVerifying MT19937 double output (seed=42):");
    for (i, &expected) in expected_double.iter().enumerate() {
        let d = mt2.next_double();
        let diff = (d - expected).abs();
        let ok = diff < 1e-17;
        println!(
            "  [{}] got={:.17e}, expected={:.17e}, diff={:.2e} {}",
            i,
            d,
            expected,
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
