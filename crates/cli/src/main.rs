mod args;
mod calc_cmds;
mod embed_cmds;
mod experimental_cmds;
mod format;

use args::{Cli, Commands};
use clap::Parser;

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Embed {
            smiles,
            seed,
            format,
        } => embed_cmds::cmd_embed(&smiles, seed, &format),
        Commands::Batch {
            input,
            output,
            seed,
            threads,
            format,
        } => embed_cmds::cmd_batch(input, output, seed, threads, &format),
        Commands::Parse { smiles } => embed_cmds::cmd_parse(&smiles),
        Commands::Info => embed_cmds::cmd_info(),
        Commands::Charges { smiles } => embed_cmds::cmd_charges(&smiles),
        Commands::Uff { smiles, coords } => embed_cmds::cmd_uff(&smiles, &coords),
        Commands::Rmsd { coords, reference } => embed_cmds::cmd_rmsd(&coords, &reference),
        Commands::Eht {
            elements,
            coords,
            k,
        } => calc_cmds::cmd_eht(&elements, &coords, k),
        Commands::Sasa {
            elements,
            coords,
            probe_radius,
        } => calc_cmds::cmd_sasa(&elements, &coords, probe_radius),
        Commands::Population { elements, coords } => calc_cmds::cmd_population(&elements, &coords),
        Commands::Dipole { elements, coords } => calc_cmds::cmd_dipole(&elements, &coords),
        Commands::Esp {
            elements,
            coords,
            spacing,
            padding,
        } => calc_cmds::cmd_esp(&elements, &coords, spacing, padding),
        Commands::Dos {
            elements,
            coords,
            sigma,
            e_min,
            e_max,
            n_points,
        } => calc_cmds::cmd_dos(&elements, &coords, sigma, e_min, e_max, n_points),
        Commands::Cell {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        } => calc_cmds::cmd_cell(a, b, c, alpha, beta, gamma),
        Commands::Assemble {
            topology,
            a,
            metal,
            geometry,
            supercell,
        } => calc_cmds::cmd_assemble(&topology, a, metal, &geometry, supercell),
        Commands::Ani { elements, coords } => calc_cmds::cmd_ani(&elements, &coords),
        Commands::Pm3 { elements, coords } => calc_cmds::cmd_pm3(&elements, &coords),
        Commands::Xtb { elements, coords } => calc_cmds::cmd_xtb(&elements, &coords),
        Commands::Gfn1 { elements, coords } => calc_cmds::cmd_gfn1(&elements, &coords),
        Commands::Gfn2 { elements, coords } => calc_cmds::cmd_gfn2(&elements, &coords),
        Commands::Hf3c { elements, coords } => calc_cmds::cmd_hf3c(&elements, &coords),
        Commands::Stereo { smiles, coords } => calc_cmds::cmd_stereo(&smiles, &coords),
        Commands::Solvation {
            elements,
            coords,
            charges,
            probe_radius,
        } => calc_cmds::cmd_solvation(&elements, &coords, &charges, probe_radius),
        Commands::Sssr { smiles } => calc_cmds::cmd_sssr(&smiles),
        Commands::Ecfp {
            smiles,
            radius,
            n_bits,
        } => calc_cmds::cmd_ecfp(&smiles, radius, n_bits),
        Commands::Tanimoto {
            smiles1,
            smiles2,
            radius,
            n_bits,
        } => calc_cmds::cmd_tanimoto(&smiles1, &smiles2, radius, n_bits),

        // Experimental commands
        #[cfg(feature = "experimental-eeq")]
        Commands::Eeq {
            elements,
            coords,
            total_charge,
        } => experimental_cmds::cmd_eeq(&elements, &coords, total_charge),
        #[cfg(feature = "experimental-alpb")]
        Commands::Alpb {
            elements,
            coords,
            charges,
            dielectric,
        } => experimental_cmds::cmd_alpb(&elements, &coords, &charges, dielectric),
        #[cfg(feature = "experimental-d4")]
        Commands::D4 {
            elements,
            coords,
            three_body,
        } => experimental_cmds::cmd_d4(&elements, &coords, three_body),
        #[cfg(feature = "experimental-cpm")]
        Commands::Cpm {
            elements,
            coords,
            mu_ev,
            dielectric,
        } => experimental_cmds::cmd_cpm(&elements, &coords, mu_ev, dielectric),
    }
}
