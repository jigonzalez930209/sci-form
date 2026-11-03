//! The built-in guide must describe the binary truthfully: every command has a
//! catalogue entry, and every example in the catalogue runs and exits 0.

use std::process::Command;

const CATALOG: &str = include_str!("../src/catalog.tsv");

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_sci-form")
}

fn commands_in_binary() -> Vec<String> {
    let out = Command::new(bin()).arg("--help").output().unwrap();
    let text = String::from_utf8_lossy(&out.stdout).to_string();
    let mut names = Vec::new();
    for line in text.lines() {
        // grouped list: two spaces, a name, padding, a description
        if let Some(rest) = line.strip_prefix("  ") {
            if let Some(name) = rest.split_whitespace().next() {
                if rest.len() > name.len() + 2
                    && name
                        .chars()
                        .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == '-')
                    && !rest.starts_with("sci-form")
                    && !rest.starts_with('-')
                {
                    names.push(name.to_string());
                }
            }
        }
    }
    names
}

#[test]
fn every_command_has_a_catalogue_entry() {
    let in_catalogue: Vec<&str> = CATALOG
        .lines()
        .filter_map(|l| l.split('\t').nth(1))
        .collect();
    let listed = commands_in_binary();
    assert!(listed.len() > 50, "the grouped list is missing: {listed:?}");
    for name in listed {
        assert!(
            in_catalogue.contains(&name.as_str()),
            "command '{name}' has no catalogue entry"
        );
    }
}

#[test]
fn every_example_runs() {
    let dir = std::path::Path::new(bin()).parent().unwrap().to_path_buf();
    let path = format!(
        "{}:{}",
        dir.display(),
        std::env::var("PATH").unwrap_or_default()
    );
    let present = commands_in_binary();
    let mut checked = 0;
    for line in CATALOG.lines().filter(|l| !l.trim().is_empty()) {
        let mut f = line.splitn(3, '\t');
        let (_topic, name, example) = (f.next().unwrap(), f.next().unwrap(), f.next().unwrap());
        if !present.iter().any(|p| p == name) {
            continue; // feature-gated command not built into this binary
        }
        let out = Command::new("sh")
            .arg("-c")
            .arg(example)
            .env("PATH", &path)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "example for '{name}' failed: {example}\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        checked += 1;
    }
    assert!(checked > 50, "only {checked} examples were checked");
}

#[test]
fn bare_invocation_prints_the_guide() {
    let out = Command::new(bin()).output().unwrap();
    assert!(out.status.success());
    let text = String::from_utf8_lossy(&out.stdout);
    assert!(text.contains("Try, in this order"));
    assert!(text.contains("sci-form embed"));
}

#[test]
fn completions_are_generated() {
    let out = Command::new(bin())
        .args(["completions", "bash"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stdout).contains("sci-form"));
}
