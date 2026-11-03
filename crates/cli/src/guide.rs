//! The built-in guide: grouped command list, per-command examples, overview.
//!
//! The catalogue (`catalog.tsv`: topic, command, example) is the single source
//! of truth. `tests/guide.rs` checks that every command of the binary has an
//! entry and that every example exits successfully, so the help can never
//! describe a command that does not exist or show an example that fails.

use clap::{Command, CommandFactory};
use std::io::Write;

use crate::args::Cli;

const CATALOG: &str = include_str!("catalog.tsv");

/// One catalogue row.
pub struct Entry {
    pub topic: &'static str,
    pub name: &'static str,
    pub example: &'static str,
}

/// All catalogue rows, in display order.
pub fn entries() -> Vec<Entry> {
    CATALOG
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let mut it = l.splitn(3, '\t');
            Entry {
                topic: it.next().unwrap_or(""),
                name: it.next().unwrap_or(""),
                example: it.next().unwrap_or(""),
            }
        })
        .collect()
}

fn summary(cmd: &Command, name: &str) -> String {
    cmd.find_subcommand(name)
        .and_then(|c| c.get_about())
        .map(|a| a.to_string())
        .unwrap_or_default()
}

/// The commands compiled into this binary, grouped by topic.
fn grouped(cmd: &Command) -> Vec<(&'static str, Vec<(&'static str, String)>)> {
    let mut groups: Vec<(&'static str, Vec<(&'static str, String)>)> = Vec::new();
    for e in entries() {
        if cmd.find_subcommand(e.name).is_none() {
            continue; // feature-gated command not built into this binary
        }
        let about = summary(cmd, e.name);
        match groups.iter_mut().find(|(t, _)| *t == e.topic) {
            Some((_, v)) => v.push((e.name, about)),
            None => groups.push((e.topic, vec![(e.name, about)])),
        }
    }
    groups
}

fn command_list(cmd: &Command, topic: Option<&str>) -> String {
    let mut out = String::new();
    for (t, items) in grouped(cmd) {
        if let Some(q) = topic {
            if !t.to_lowercase().contains(&q.to_lowercase()) {
                continue;
            }
        }
        out.push_str(&format!("{t}\n"));
        for (n, about) in items {
            out.push_str(&format!("  {n:<26}{about}\n"));
        }
        out.push('\n');
    }
    out
}

/// Text appended to `sci-form --help`.
pub fn root_after_help(cmd: &Command) -> String {
    let n: usize = grouped(cmd).iter().map(|(_, v)| v.len()).sum();
    format!(
        "Start here:\n  sci-form embed \"CCO\"            3D structure of ethanol\n  \
         sci-form commands               all {n} commands, grouped by topic\n  \
         sci-form examples pm3 --run     a ready-to-run example (and its output)\n  \
         sci-form help <command>         every option of one command\n\n\
         Most calculations take a geometry as two JSON arguments, <ELEMENTS> and\n\
         <COORDS> (flat x,y,z in Å). `sci-form embed` prints both.\n\n\
         Documentation: https://sci-form.dev"
    )
}

/// The banner printed by a bare `sci-form`.
pub fn overview() {
    let cmd = Cli::command();
    println!(
        "sci-form {} — computational chemistry from the command line\n",
        env!("CARGO_PKG_VERSION")
    );
    println!("Try, in this order:\n");
    println!("  sci-form embed \"CCO\"                       build a 3D structure from SMILES");
    println!("  sci-form charges \"CCO\"                     partial charges");
    println!(
        "  sci-form examples pm3 --run                a quantum calculation, with its output\n"
    );
    println!("{}", command_list(&cmd, None));
    println!("More:");
    println!("  sci-form commands [topic]    this list, optionally one topic (e.g. quantum)");
    println!("  sci-form examples [command]  a ready-to-run example for each command");
    println!("  sci-form help <command>      every option of one command");
    println!("  sci-form completions bash    shell completion (bash, zsh, fish, powershell)");
    println!("  sci-form --version           version\n");
    println!("Documentation and live examples: https://sci-form.dev");
}

pub fn cmd_commands(topic: Option<String>) {
    let cmd = Cli::command();
    let text = command_list(&cmd, topic.as_deref());
    if text.is_empty() {
        eprintln!(
            "No topic matches '{}'. Topics: {}.",
            topic.unwrap_or_default(),
            grouped(&cmd)
                .iter()
                .map(|(t, _)| *t)
                .collect::<Vec<_>>()
                .join("; ")
        );
        std::process::exit(2);
    }
    print!("{text}");
}

pub fn cmd_examples(command: Option<String>, run: bool) {
    let cmd = Cli::command();
    let all = entries();
    let chosen: Vec<&Entry> = match &command {
        Some(c) => {
            let v: Vec<&Entry> = all.iter().filter(|e| e.name == c).collect();
            if v.is_empty() {
                eprintln!("No example for '{c}'. Run `sci-form commands` for the list.");
                std::process::exit(2);
            }
            v
        }
        None => all
            .iter()
            .filter(|e| cmd.find_subcommand(e.name).is_some())
            .collect(),
    };
    for e in chosen {
        println!("# {}\n{}", summary(&cmd, e.name), e.example);
        if run {
            println!();
            let _ = std::io::stdout().flush();
            run_example(e.example);
        }
        println!();
    }
}

fn run_example(example: &str) {
    // Run through the shell with this very binary first on the PATH, so the
    // example behaves exactly as typed.
    let exe = std::env::current_exe().unwrap_or_default();
    let dir = exe.parent().map(|p| p.to_path_buf()).unwrap_or_default();
    let mut path = std::ffi::OsString::from(&dir);
    if let Some(old) = std::env::var_os("PATH") {
        path.push(if cfg!(windows) { ";" } else { ":" });
        path.push(old);
    }
    let status = if cfg!(windows) {
        eprintln!("(--run needs a POSIX shell; copy the example into your shell instead)");
        return;
    } else {
        std::process::Command::new("sh")
            .arg("-c")
            .arg(example)
            .env("PATH", path)
            .status()
    };
    if let Err(e) = status {
        eprintln!("could not run the example: {e}");
    }
}

pub fn cmd_completions(shell: clap_complete::Shell) {
    let mut cmd = Cli::command();
    clap_complete::generate(shell, &mut cmd, "sci-form", &mut std::io::stdout());
}

/// The command with the guide wired in: the grouped list as `--help` footer
/// and one example under each command's own help.
pub fn build_command() -> Command {
    let mut cmd = Cli::command();
    let after = root_after_help(&cmd);
    let list = command_list(&cmd, None)
        .replace('{', "{{")
        .replace('}', "}}");
    cmd = cmd.after_help(after).help_template(format!(
        "{{about-with-newline}}\n{{usage-heading}} {{usage}}\n\nCommands:\n{list}\
         {{options}}\n\n{{after-help}}"
    ));
    for e in entries() {
        if cmd.find_subcommand(e.name).is_some() && !e.example.is_empty() {
            let text = format!("Example:\n  {}", e.example.replace('\n', "\n  "));
            cmd = cmd.mut_subcommand(e.name, |c| c.after_help(text));
        }
    }
    cmd
}
