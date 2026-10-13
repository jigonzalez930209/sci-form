//! Latency probe for the ETKDG pipeline (profiling target, validation §14).
//!
//! With no argument: reports per-molecule single-call latency (median of N).
//! With `--profile`: runs a fixed workload for `perf record` attribution.
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cases = [
        ("CCO", "ethanol"),
        ("c1ccccc1", "benzene"),
        ("CC(=O)Oc1ccccc1C(=O)O", "aspirin"),
        ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "caffeine"),
        (
            "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C",
            "cholesterol",
        ),
    ];
    if args.iter().any(|a| a == "--cold") {
        // First call in a fresh process: includes one-time initialisation
        // (SMARTS torsion library, tables) plus parse and embed.
        for (smi, name) in cases {
            let t = Instant::now();
            let _ = sci_form::embed(smi, 42);
            println!(
                "{name:12} first-call {:8.2} ms",
                t.elapsed().as_secs_f64() * 1e3
            );
        }
        return;
    }
    if args.iter().any(|a| a == "--profile") {
        let mut total = 0usize;
        for _ in 0..40 {
            for (smi, _) in cases {
                total += sci_form::embed(smi, 42).coords.len();
            }
        }
        println!("probe done, coords touched: {total}");
        return;
    }
    let reps: usize = args.iter().find_map(|a| a.parse().ok()).unwrap_or(20);
    println!("{:12} {:>10} {:>10}", "molecule", "median ms", "min ms");
    for (smi, name) in cases {
        // warm-up
        let _ = sci_form::embed(smi, 42);
        let mut times = Vec::with_capacity(reps);
        for _ in 0..reps {
            let t = Instant::now();
            let _ = sci_form::embed(smi, 42);
            times.push(t.elapsed().as_secs_f64() * 1e3);
        }
        times.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "{:12} {:>10.2} {:>10.2}",
            name,
            times[times.len() / 2],
            times[0]
        );
    }
}
