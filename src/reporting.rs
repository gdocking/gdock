//! Run-output formatting: terminal display and report files (`metrics.tsv`,
//! `epoch.tsv.gz`, `sampling.tsv`).

use colored::*;
use flate2::write::GzEncoder;
use flate2::Compression;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::chromosome::Chromosome;
use crate::constants::{EnergyWeights, CONVERGENCE_THRESHOLD, CONVERGENCE_WINDOW};
use crate::evaluator::Metrics;
use crate::hall_of_fame::HallOfFameEntry;

/// Truecolor keeps this a fixed crimson regardless of the terminal's ANSI theme.
pub(crate) fn brand(s: &str) -> ColoredString {
    s.truecolor(222, 45, 38).bold()
}

/// CAPRI quality bands: high ≥0.80, medium ≥0.49, acceptable ≥0.23, incorrect below.
pub(crate) fn colorize_dockq(dockq: f64) -> ColoredString {
    let s = format!("{:.3}", dockq);
    if dockq >= 0.80 {
        s.green()
    } else if dockq >= 0.49 {
        s.yellow()
    } else if dockq >= 0.23 {
        s.bright_yellow()
    } else {
        s.red()
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn header(
    version: &str,
    debug_mode: bool,
    receptor_file: &str,
    ligand_file: &str,
    reference_file: &Option<String>,
    n_pairs: usize,
    num_restraints: usize,
    weights: &EnergyWeights,
) {
    let rule = "━".repeat(80);
    println!("{}", rule.truecolor(160, 160, 160));
    println!("{}", rule.truecolor(160, 160, 160));
    println!(
        "\n{} {} {}",
        brand("🧬 gdock"),
        format!("v{version}").truecolor(160, 160, 160),
        "· information-driven protein-protein docking".truecolor(160, 160, 160)
    );
    if debug_mode {
        println!(
            "\n{}",
            "   ⚠️  DEBUG MODE: using DockQ as the fitness function\n"
                .yellow()
                .bold()
        );
    }
    println!(
        "  {}  {}",
        "Receptor  ".truecolor(160, 160, 160),
        receptor_file
    );
    println!(
        "  {}  {}",
        "Ligand    ".truecolor(160, 160, 160),
        ligand_file
    );
    match reference_file {
        Some(ref_file) => println!(
            "  {}  {} {}",
            "Reference ".truecolor(160, 160, 160),
            ref_file,
            "(DockQ mode)".truecolor(160, 160, 160)
        ),
        None => println!(
            "  {}  {}",
            "Reference ".truecolor(160, 160, 160),
            "none (score-only mode)".yellow()
        ),
    }
    println!(
        "  {}  {} pairs {} {} ambiguous restraints",
        "Restraints".truecolor(160, 160, 160),
        n_pairs.to_string().bold(),
        "→".truecolor(160, 160, 160),
        num_restraints.to_string().bold()
    );
    println!(
        "  {}  {} {:.2}   {} {:.2}   {} {:.2}   {} {:.2}   {} {:.2}",
        "Weights   ".truecolor(160, 160, 160),
        "vdw".truecolor(160, 160, 160),
        weights.vdw,
        "elec".truecolor(160, 160, 160),
        weights.elec,
        "desolv".truecolor(160, 160, 160),
        weights.desolv,
        "air".truecolor(160, 160, 160),
        weights.air,
        "clash".truecolor(160, 160, 160),
        weights.clash
    );
    println!("\n{}", rule.truecolor(160, 160, 160));
    println!("{}\n", rule.truecolor(160, 160, 160));
}

pub(crate) fn new_progress_bar(total: u64) -> ProgressBar {
    let progress = ProgressBar::new(total);
    progress.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.160} [{bar:40.160/244}] {pos}/{len} gens | {msg}")
            .unwrap()
            .progress_chars("█▓░"),
    );
    progress
}

pub(crate) fn evolution_start(max_generations: u64, population_size: u64) {
    println!(
        "{} Evolving {} generations {} {} individuals",
        brand("🧬"),
        max_generations.to_string().bold(),
        "·".truecolor(160, 160, 160),
        population_size.to_string().bold()
    );
}

/// `mean_fitness` is only used when `metrics` is `None` (no reference given).
pub(crate) fn progress_message(
    debug_mode: bool,
    best_fitness: f64,
    metrics: Option<&[Metrics]>,
    best_idx: usize,
    mean_fitness: f64,
) -> String {
    match metrics {
        Some(metrics) => {
            let best_dockq = metrics[best_idx].dockq;
            if debug_mode {
                let mean_dockq =
                    metrics.iter().map(|m| m.dockq).sum::<f64>() / metrics.len() as f64;
                format!(
                    "best DockQ {} · mean {:.3}",
                    colorize_dockq(best_dockq),
                    mean_dockq
                )
            } else {
                format!(
                    "best score {:.1} · DockQ {}",
                    best_fitness,
                    colorize_dockq(best_dockq)
                )
            }
        }
        None => format!("best score {:.1} · mean {:.1}", best_fitness, mean_fitness),
    }
}

pub(crate) fn finish_progress(bar: &ProgressBar, generations_run: u64, converged_early: bool) {
    if converged_early {
        // abandon_with_message doesn't snap position to len like finish_with_message does.
        bar.set_position(generations_run);
        bar.abandon_with_message(format!(
            "{} converged at gen {} — no improvement > {}% over {} gens",
            "✓".green(),
            generations_run,
            CONVERGENCE_THRESHOLD * 100.0,
            CONVERGENCE_WINDOW
        ));
    } else {
        bar.finish_with_message(format!(
            "{} completed {} generations",
            "✓".green(),
            generations_run
        ));
    }
    println!();
}

pub(crate) fn clustering(hof_len: usize) {
    println!(
        "\n{} Clustering {} Hall-of-Fame structures",
        brand("🔬"),
        hof_len.to_string().bold()
    );
}

pub(crate) fn saving_results() {
    println!("\n{}", brand("💾 Saving Results"));
}

pub(crate) fn final_metrics(best_score: &Metrics, best_dockq: &Metrics) {
    println!("\n{}", brand("📊 Final Metrics"));
    println!(
        "  {} DockQ={:.3} RMSD={:.2}Å iRMSD={:.2}Å FNAT={:.3}",
        "Best by score:".green(),
        best_score.dockq,
        best_score.rmsd,
        best_score.irmsd,
        best_score.fnat
    );
    println!(
        "  {} DockQ={:.3} RMSD={:.2}Å iRMSD={:.2}Å FNAT={:.3}",
        "Best by DockQ:".green(),
        best_dockq.dockq,
        best_dockq.rmsd,
        best_dockq.irmsd,
        best_dockq.fnat
    );
}

pub(crate) fn saved_file(path: &Path) {
    println!("  {} {}", "✓".green(), path.display());
}

/// Columns must match what `model_row` prints for the same `has_*` flags.
pub(crate) fn model_section_header(
    title: &str,
    qualifier: &str,
    has_score: bool,
    has_dockq: bool,
    has_cluster: bool,
) {
    println!("\n{} {}", brand(title), qualifier.truecolor(160, 160, 160));

    let mut cols = vec![format!("{:<8}", "model")];
    if has_score {
        cols.push(format!("{:>8}", "score"));
    }
    if has_dockq {
        cols.push(format!("{:<5}", "DockQ"));
    }
    if has_cluster {
        cols.push(format!("{:>7}", "cluster"));
    }
    println!(
        "{}",
        format!("  {}", cols.join("  ")).truecolor(160, 160, 160)
    );
}

pub(crate) fn model_row(
    name: &str,
    score: Option<f64>,
    dockq: Option<f64>,
    cluster: Option<usize>,
) {
    let mut parts = vec![format!("{:<8}", name).bold().to_string()];
    if let Some(s) = score {
        parts.push(format!("{:>8.1}", s));
    }
    if let Some(d) = dockq {
        parts.push(colorize_dockq(d).to_string());
    }
    if let Some(c) = cluster {
        parts.push(format!("{:>7}", c).bold().to_string());
    }
    println!("  {}", parts.join("  "));
}

pub(crate) fn summary(n_models: usize, epoch: bool, out_dir: &Path) {
    let epoch_note = if epoch { " + epoch.tsv.gz" } else { "" };
    println!(
        "\n{} {} models + metrics.tsv{} written to {}",
        brand("💾"),
        n_models.to_string().bold(),
        epoch_note,
        out_dir.display().to_string().truecolor(160, 160, 160)
    );
}

pub(crate) fn done() {
    println!("\n{}", "✨ Done!".bold().green());
}

pub(crate) fn sampling_header() {
    println!("\n{}", brand("📦 Sampling output"));
}

pub(crate) fn sampling_done(count: usize, dir: &Path, tsv_path: &Path) {
    println!(
        "  {} {} structures written to {}",
        "✓".green(),
        count,
        dir.display()
    );
    println!("  {} {}", "✓".green(), tsv_path.display());
}

// Report files (metrics.tsv, epoch.tsv.gz, sampling.tsv)

fn write_epoch_header<W: Write>(w: &mut W, has_eval: bool) {
    let cols = "generation\tindividual\talpha\tbeta\tgamma\tdx\tdy\tdz\tscore\tvdw\telec\tdesolv\tair\tclash\trestraint_penalty";
    if has_eval {
        writeln!(w, "{cols}\tdockq\tirmsd\trmsd\tfnat").unwrap();
    } else {
        writeln!(w, "{cols}").unwrap();
    }
}

fn write_epoch_row<W: Write>(
    w: &mut W,
    gen: u64,
    individual: usize,
    chr: &Chromosome,
    metrics: Option<&Metrics>,
) {
    write!(
        w,
        "{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        gen,
        individual,
        chr.genes[0],
        chr.genes[1],
        chr.genes[2],
        chr.genes[3],
        chr.genes[4],
        chr.genes[5],
        chr.fitness,
        chr.vdw,
        chr.elec,
        chr.desolv,
        chr.air,
        chr.clash,
        chr.restraint_penalty
    )
    .unwrap();
    match metrics {
        Some(m) => writeln!(
            w,
            "\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
            m.dockq, m.irmsd, m.rmsd, m.fnat
        )
        .unwrap(),
        None => writeln!(w).unwrap(),
    }
}

pub(crate) struct EpochWriter {
    writer: BufWriter<GzEncoder<File>>,
}

impl EpochWriter {
    pub(crate) fn create(out_dir: &Path, has_eval: bool) -> Self {
        let path = out_dir.join("epoch.tsv.gz");
        let file = File::create(&path).expect("Failed to create epoch.tsv.gz");
        let mut writer = BufWriter::new(GzEncoder::new(file, Compression::default()));
        write_epoch_header(&mut writer, has_eval);
        EpochWriter { writer }
    }

    pub(crate) fn write_row(
        &mut self,
        gen: u64,
        individual: usize,
        chr: &Chromosome,
        metrics: Option<&Metrics>,
    ) {
        write_epoch_row(&mut self.writer, gen, individual, chr, metrics);
    }

    /// Must be called to flush the buffer and write the gzip trailer.
    pub(crate) fn finish(self) {
        let encoder = self
            .writer
            .into_inner()
            .expect("Failed to flush epoch buffer");
        encoder.finish().expect("Failed to finish epoch.tsv.gz");
    }
}

pub(crate) fn write_noclust_metrics_header<W: Write>(w: &mut W) {
    writeln!(
        w,
        "model\tdockq\trmsd\tirmsd\tfnat\tscore\tgenerations_run\tconverged_early"
    )
    .unwrap();
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn write_noclust_metrics_row<W: Write>(
    w: &mut W,
    label: &str,
    metrics: &Metrics,
    fitness: f64,
    generations_run: u64,
    converged_early: bool,
) {
    writeln!(
        w,
        "{label}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{generations_run}\t{converged_early}",
        metrics.dockq, metrics.rmsd, metrics.irmsd, metrics.fnat, fitness
    )
    .unwrap();
}

pub(crate) fn write_clustered_metrics_header<W: Write>(w: &mut W, has_eval: bool) {
    if has_eval {
        writeln!(
            w,
            "model\tcluster_size\tscore\tdockq\trmsd\tirmsd\tfnat\tgenerations_run\tconverged_early"
        )
        .unwrap();
    } else {
        writeln!(
            w,
            "model\tcluster_size\tscore\tgenerations_run\tconverged_early"
        )
        .unwrap();
    }
}

/// `cluster_size: None` (ranked models) is written as `-`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_clustered_metrics_row<W: Write>(
    w: &mut W,
    model_name: &str,
    cluster_size: Option<usize>,
    fitness: f64,
    metrics: Option<&Metrics>,
    generations_run: u64,
    converged_early: bool,
) {
    let cluster_str = cluster_size
        .map(|c| c.to_string())
        .unwrap_or_else(|| "-".to_string());
    match metrics {
        Some(m) => writeln!(
            w,
            "{model_name}\t{cluster_str}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{generations_run}\t{converged_early}",
            fitness, m.dockq, m.rmsd, m.irmsd, m.fnat
        )
        .unwrap(),
        None => writeln!(
            w,
            "{model_name}\t{cluster_str}\t{:.4}\t{generations_run}\t{converged_early}",
            fitness
        )
        .unwrap(),
    }
}

pub(crate) fn write_sampling_header<W: Write>(w: &mut W) {
    writeln!(w, "model\tscore\tvdw\telec\tdesolv\tair\tclash").unwrap();
}

pub(crate) fn write_sampling_row<W: Write>(w: &mut W, model_name: &str, entry: &HallOfFameEntry) {
    writeln!(
        w,
        "{model_name}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        entry.fitness, entry.vdw, entry.elec, entry.desolv, entry.air, entry.clash
    )
    .unwrap();
}
