---
title: 'gdock: Information-driven protein-protein docking using genetic algorithms'
tags:
  - Rust
  - bioinformatics
  - protein docking
  - genetic algorithm
  - structural biology
authors:
  - name: Rodrigo V. Honorato
    orcid: 0000-0001-5267-3002
    corresponding: true
    affiliation: 1
affiliations:
  - name: Independent Researcher, Netherlands
    index: 1
date: 28 January 2026
bibliography: paper.bib
---

# Summary

_placeholder_

# Statement of need

Information-driven docking incorporates experimental data—from mutagenesis,
cross-linking mass spectrometry, or NMR—to guide protein complex structure
prediction [@vannoort2021information]. HADDOCK [@dominguez2003haddock] pioneered
this approach; other tools including ClusPro [@kozakov2017cluspro], ZDOCK
[@pierce2014zdock], and LightDock [@lightdock2018] also support restraint
integration to varying degrees.

`gdock` contributes to this ecosystem as a fast, minimal implementation:

- **Single binary**: No runtime dependencies or environment setup
- **Speed**: Median 16 seconds per complex on standard hardware
- **Rust implementation**: Native energy functions, memory-safe, easily auditable
- **CLI-first**: Designed for scripting and pipeline integration

The tool provides an accessible option for researchers who need restraint-driven
docking without the overhead of larger software packages.

# Software design

`gdock` is a rewrite in Rust of an earlier Python prototype. The ~7,000 line
codebase compiles to a single statically-linked binary with no external runtime
dependencies.

**Search algorithm.** A genetic algorithm explores rigid-body transformations
of the ligand relative to the receptor. Each chromosome encodes six genes:
three Euler angles (α, β, γ) for rotation and three displacement values
(x, y, z) for translation. A generation consists of a population of
chromosomes that evolves through tournament selection, uniform crossover,
creep mutation (Gaussian perturbations for local refinement), and elitism.
Fitness evaluation is parallelized across the population. The search
terminates early upon convergence.

**Scoring function.** The energy function combines four terms:

$$E_{total} = w_{vdw} E_{vdw} + w_{elec} E_{elec} + w_{desolv} E_{desolv} + w_{air} E_{air}$$

- $E_{vdw}$: Soft-core Lennard-Jones potential that remains finite at short
  distances, allowing the search to explore conformations with minor clashes
- $E_{elec}$: Coulombic interactions with distance-dependent dielectric
  ($\varepsilon = r$) to dampen long-range effects
- $E_{desolv}$: Empirical atomic solvation parameters penalizing burial of
  polar atoms and rewarding burial of hydrophobic atoms
- $E_{air}$: Flat-bottom harmonic potential on Cα–Cα distances between
  user-specified residue pairs (no penalty within 0–7 Å), inspired by
  HADDOCK's ambiguous interaction restraints [@dominguez2003haddock]

**Weight calibration.** The weights $w_{vdw}$, $w_{elec}$, and $w_{desolv}$ were
calibrated using the Dockground decoy set [@dockground2008], which provides 100
decoy structures per complex with exactly one near-native conformation. A grid
search tested weight combinations by re-scoring all decoys and measuring how
often the near-native structure ranked in the top 50. The final weights
($w_{vdw}=0.4$, $w_{elec}=0.05$, $w_{desolv}=3.4$) maximize this ranking
performance. The restraint weight $w_{air}$ is fixed at 1.0 since calibration
was performed without restraints.

**Output.** Final models are clustered using Fraction of Common Contacts (FCC)
[@rodrigues2012fcc], re-implemented natively in Rust, and ranked by score,
providing both diverse and top-scoring solutions.

**Code quality.** The codebase includes 174 unit tests covering parsing, energy
calculations, and algorithm behavior. Continuous integration enforces code
formatting (`rustfmt`), linting (`clippy` with warnings as errors), and test
passage on every commit. Rust's ownership model provides compile-time guarantees
against data races and use-after-free errors.

# Benchmarking

`gdock` was validated on 271 complexes from the Protein-Protein Docking
Benchmark v5 [@vreven2015updates], a standard dataset for assessing docking
methods. Restraints were derived from native contacts using a 5 Å distance
cutoff, simulating ideal information-driven scenarios.

Using the DockQ metric [@basu2016dockq] to assess model quality, `gdock`
achieved a 95.9% success rate (260/271 complexes with at least one acceptable
model, DockQ ≥ 0.23). Medium-quality models (DockQ ≥ 0.49) were obtained for
55.7% of complexes, and high-quality models (DockQ ≥ 0.80) for 3.7%
(\autoref{fig:dockq}).

![Distribution of docking quality across 271 benchmark complexes. Each complex
is categorized by its best DockQ score among 10 output
models.\label{fig:dockq}](plot_dockq_categories.pdf)

Performance benchmarks on a 94-core machine show a median docking time of 15.7
seconds per complex (mean: 52.7 seconds), with 56% of cases completing within
20 seconds (\autoref{fig:timing}).

![Distribution of docking times across benchmark complexes. The right-skewed
distribution reflects that larger complexes require more computational
time.\label{fig:timing}](plot_timing_histogram.pdf)

All scripts required to reproduce the calibration and benchmarking experiments
are available in a separate repository:
<https://github.com/rvhonorato/gdock-benchmark>.

# Acknowledgements

# References
