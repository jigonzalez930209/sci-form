/// Benchmark: diverse molecule conformer generation
///
/// Tests sci-form against a curated set of molecules spanning all chemical
/// functional groups, halogens, heteroatoms, stereo, and ring systems.
/// Reports parse rate, embed rate, and geometry quality per category.
use std::collections::HashMap;
use std::fs;
use std::time::Instant;

struct MolEntry {
    smiles: String,
    name: String,
    category: String,
}

fn load_diverse_molecules() -> Vec<MolEntry> {
    let content = fs::read_to_string("tests/fixtures/diverse_molecules.smi")
        .expect("Cannot read tests/fixtures/diverse_molecules.smi");
    content
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let parts: Vec<&str> = l.split('\t').collect();
            if parts.len() >= 3 {
                Some(MolEntry {
                    smiles: parts[0].to_string(),
                    name: parts[1].to_string(),
                    category: parts[2].to_string(),
                })
            } else {
                None
            }
        })
        .collect()
}

fn bond_length(coords: &[f64], a: usize, b: usize) -> f64 {
    let dx = coords[a * 3] - coords[b * 3];
    let dy = coords[a * 3 + 1] - coords[b * 3 + 1];
    let dz = coords[a * 3 + 2] - coords[b * 3 + 2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn check_geometry(result: &sci_form::ConformerResult) -> (bool, Vec<String>) {
    let mut issues = Vec::new();
    if result.error.is_some() || result.coords.is_empty() {
        return (false, vec!["embed failed".to_string()]);
    }

    let n = result.num_atoms;

    // Check for NaN/Inf
    for v in &result.coords {
        if !v.is_finite() {
            issues.push("NaN/Inf in coordinates".to_string());
            return (false, issues);
        }
    }

    // Check bond lengths
    for (a, b, _order) in &result.bonds {
        let bl = bond_length(&result.coords, *a, *b);
        if !(0.7..=2.5).contains(&bl) {
            issues.push(format!("bond {}-{}: {:.3} Å", a, b, bl));
        }
    }

    // Check for atom clashes (non-bonded < 0.5 Å)
    let bonded: std::collections::HashSet<(usize, usize)> = result
        .bonds
        .iter()
        .flat_map(|(a, b, _)| vec![(*a, *b), (*b, *a)])
        .collect();

    for i in 0..n {
        for j in (i + 1)..n {
            if bonded.contains(&(i, j)) {
                continue;
            }
            let dx = result.coords[i * 3] - result.coords[j * 3];
            let dy = result.coords[i * 3 + 1] - result.coords[j * 3 + 1];
            let dz = result.coords[i * 3 + 2] - result.coords[j * 3 + 2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist < 0.5 {
                issues.push(format!("clash {}-{}: {:.3} Å", i, j, dist));
            }
        }
    }

    (issues.is_empty(), issues)
}

#[test]
fn test_diverse_molecules() {
    let entries = load_diverse_molecules();
    let total = entries.len();
    println!(
        "\n=== DIVERSE MOLECULE BENCHMARK ({} molecules) ===\n",
        total
    );

    let mut parse_ok = 0;
    let mut embed_ok = 0;
    let mut geom_ok = 0;
    let mut category_stats: HashMap<String, (usize, usize, usize, usize)> = HashMap::new();
    let mut failures: Vec<(String, String, String)> = Vec::new();
    let mut geom_issues: Vec<(String, String, Vec<String>)> = Vec::new();

    let start = Instant::now();

    for entry in &entries {
        let cat = &entry.category;
        let stat = category_stats.entry(cat.clone()).or_insert((0, 0, 0, 0));
        stat.0 += 1; // total

        // Parse
        match sci_form::parse(&entry.smiles) {
            Ok(_) => {
                parse_ok += 1;
                stat.1 += 1;
            }
            Err(e) => {
                failures.push((
                    entry.name.clone(),
                    entry.smiles.clone(),
                    format!("parse: {}", e),
                ));
                continue;
            }
        }

        // Embed
        let result = sci_form::embed(&entry.smiles, 42);
        if result.error.is_some() {
            embed_ok += 0;
            failures.push((
                entry.name.clone(),
                entry.smiles.clone(),
                format!("embed: {}", result.error.as_deref().unwrap_or("?")),
            ));
            continue;
        }
        embed_ok += 1;
        stat.2 += 1;

        // Geometry
        let (good, issues) = check_geometry(&result);
        if good {
            geom_ok += 1;
            stat.3 += 1;
        } else {
            geom_issues.push((entry.name.clone(), entry.smiles.clone(), issues));
        }
    }

    let elapsed = start.elapsed();

    // Summary
    println!("Total: {}", total);
    println!(
        "Parse OK: {}/{} ({:.1}%)",
        parse_ok,
        total,
        100.0 * parse_ok as f64 / total as f64
    );
    println!(
        "Embed OK: {}/{} ({:.1}%)",
        embed_ok,
        total,
        100.0 * embed_ok as f64 / total as f64
    );
    println!(
        "Geometry OK: {}/{} ({:.1}%)",
        geom_ok,
        total,
        100.0 * geom_ok as f64 / total as f64
    );
    println!(
        "Time: {:.2}s ({:.1} mol/s)\n",
        elapsed.as_secs_f64(),
        total as f64 / elapsed.as_secs_f64()
    );

    // Per-category breakdown
    println!(
        "{:<16} {:>5} {:>5} {:>5} {:>5}",
        "Category", "Total", "Parse", "Embed", "Geom"
    );
    println!("{}", "-".repeat(52));
    let mut cats: Vec<_> = category_stats.iter().collect();
    cats.sort_by_key(|(k, _)| (*k).clone());
    for (cat, (tot, par, emb, geo)) in &cats {
        println!("{:<16} {:>5} {:>5} {:>5} {:>5}", cat, tot, par, emb, geo);
    }

    // Failures
    if !failures.is_empty() {
        println!("\n--- FAILURES ({}) ---", failures.len());
        for (name, smi, reason) in &failures {
            println!("  {} [{}]: {}", name, smi, reason);
        }
    }

    // Geometry issues
    if !geom_issues.is_empty() {
        println!("\n--- GEOMETRY ISSUES ({}) ---", geom_issues.len());
        for (name, smi, issues) in &geom_issues {
            println!("  {} [{}]:", name, smi);
            for issue in issues {
                println!("    {}", issue);
            }
        }
    }

    // Assertions
    let parse_rate = parse_ok as f64 / total as f64;
    let embed_rate = embed_ok as f64 / total as f64;
    let geom_rate = geom_ok as f64 / total as f64;

    assert!(
        parse_rate >= 0.95,
        "Parse rate {:.1}% < 95%",
        parse_rate * 100.0
    );
    assert!(
        embed_rate >= 0.90,
        "Embed rate {:.1}% < 90%",
        embed_rate * 100.0
    );
    assert!(
        geom_rate >= 0.85,
        "Geometry rate {:.1}% < 85%",
        geom_rate * 100.0
    );
}
