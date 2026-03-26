//! Debug: exercise a single molecule embedding with LOG_ATTEMPTS enabled
//! to diagnose failures in the conformer pipeline.
//!
//! Category: debug

use std::time::Instant;

fn main() {
    std::env::set_var("LOG_ATTEMPTS", "1");
    // Only test one molecule for diagnosis (single seed, no retry)
    let smi = "O=c1c2cccc3c(Nc4c(F)cccc4F)ccc(c32)c2nc3ccccc3n12";
    let mol = sci_form::parse(smi).unwrap();
    let start = Instant::now();
    let result = sci_form::conformer::generate_3d_conformer(&mol, 42);
    let elapsed = start.elapsed().as_secs_f64();
    match result {
        Ok(_) => println!("OK in {:.1}s", elapsed),
        Err(e) => {
            println!("FAIL in {:.1}s: {}", elapsed, e);
            // Try seed 997
            let start2 = Instant::now();
            let result2 = sci_form::conformer::generate_3d_conformer(&mol, 997);
            match result2 {
                Ok(_) => println!("OK with seed 997 in {:.1}s", start2.elapsed().as_secs_f64()),
                Err(e2) => println!(
                    "FAIL with seed 997 in {:.1}s: {}",
                    start2.elapsed().as_secs_f64(),
                    e2
                ),
            }
        }
    }
}
