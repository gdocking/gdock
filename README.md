# gdock

**Information-driven protein-protein docking using a genetic algorithm**

[![ci](https://github.com/rvhonorato/gdock/actions/workflows/ci.yml/badge.svg)](https://github.com/rvhonorato/gdock/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/23671025da8a4334a754d8d5af76a34a)](https://app.codacy.com/gh/rvhonorato/gdock/dashboard)

![License](https://img.shields.io/badge/license-0BSD-blue)

<img src="imgs/gdock_logo.png" width="350">

gdock is a fast protein-protein docking tool written in Rust that uses
restraints and energy components to guide the docking process. It combines a
genetic algorithm with physics-based scoring to find optimal protein-protein
complexes.

> **Note**: This project is currently under review for publication in the
> [Journal of Open Source Software (JOSS)](https://joss.theoj.org/). A stable
> v2.0.0 release will follow upon acceptance.

## Features

- **Fast**: Genetic algorithm with early stopping and elitism
- **Information-driven**: Uses residue restraints to guide docking
- **Flexible scoring**: Configurable energy weights (VDW, electrostatics,
  desolvation, restraints)
- **Quality metrics**: Optional DockQ calculation when reference structure is
  provided
- **Clustering**: FCC-based clustering to group similar solutions

## Quick Start

```bash
# Clone and build
git clone https://github.com/rvhonorato/gdock
cd gdock
cargo build --release

# Run docking with example data
./target/release/gdock run \
  --receptor data/2oob_A.pdb \
  --ligand data/2oob_B.pdb \
  --restraints 933:6,936:8,940:42,941:44,946:45,950:46
```

Most docking runs complete in ~15 seconds on standard hardware.

## Requirements

- [Rust](https://www.rust-lang.org/tools/install) (1.70 or later)

## Installation

```bash
git clone https://github.com/rvhonorato/gdock
cd gdock
cargo build --release
```

The binary will be available at `./target/release/gdock`.

Upon stable release, pre-built binaries and `cargo install gdock` will also be
available.

## Usage

gdock has three subcommands: `run`, `score`, and `restraints`.

### Docking (`run`)

Run the full genetic algorithm docking:

```bash
gdock run \
  --receptor receptor.pdb \
  --ligand ligand.pdb \
  --restraints 933:6,936:8,940:42
```

With a reference structure for DockQ calculation:

```bash
gdock run \
  --receptor receptor.pdb \
  --ligand ligand.pdb \
  --restraints 933:6,936:8,940:42 \
  --reference native.pdb
```

Additional options:

- `-o, --output-dir <DIR>`: Output directory (default: current directory)
- `-n, --nproc <NUM>`: Number of processors (default: total - 2)
- `--no-clust`: Disable clustering
- `--w_vdw`, `--w_elec`, `--w_desolv`, `--w_air`: Custom energy weights

### Scoring (`score`)

Calculate energy components without running the GA:

```bash
gdock score \
  --receptor receptor.pdb \
  --ligand ligand.pdb \
  --restraints 933:6,936:8,940:42
```

### Generate restraints (`restraints`)

Generate restraints from interface contacts in a native structure:

```bash
gdock restraints \
  --receptor receptor_ref.pdb \
  --ligand ligand_ref.pdb \
  --cutoff 5.0
```

## Input Format

### PDB Files

- **Receptor**: PDB file containing the receptor protein (single chain)
- **Ligand**: PDB file containing the ligand protein (single chain)
- **Reference** (optional): PDB file containing the native complex

### Restraints

Comma-separated list of residue pairs in `receptor:ligand` format:

```text
933:6,936:8,940:42
```

These indicate which residues should be in contact, based on experimental data
or other information sources.

## Output

- `model_X.pdb`: Cluster representatives (unless `--no-clust`)
- `ranked_X.pdb`: Top 5 models ranked by score
- `metrics.tsv`: Tab-separated file with scores and metrics

Output structures can be visualized with molecular viewers such as
[PyMOL](https://pymol.org/) or [ChimeraX](https://www.cgl.ucsf.edu/chimerax/).

## Algorithm

gdock uses:

- **Genetic Algorithm**: Population of 150, elitism (top 5), tournament
  selection
- **Energy Function**: VDW + Electrostatics + Desolvation + AIR restraints
- **Restraints**: Flat-bottom potential (0-7 Angstrom) for specified residue pairs
- **Early Stopping**: Converges when no improvement for 10 generations
- **Clustering**: FCC-based clustering of final population

## Testing

Run the test suite:

```bash
cargo test
```

The test suite includes 174 tests covering parsing, energy calculations, and
algorithm behavior.

## Example

Using the test data included in the repository:

```bash
gdock run \
  --receptor data/2oob_A.pdb \
  --ligand data/2oob_B.pdb \
  --restraints 933:6,936:8,940:42,941:44,946:45,950:46 \
  --reference data/2oob.pdb \
  --output-dir results/
```

This will produce:

- `results/model_*.pdb` — Cluster representatives ranked by cluster size
- `results/ranked_*.pdb` — Top 5 models ranked by score
- `results/metrics.tsv` — Scores and DockQ values for all models

## Relevant repositories

- [`gdock-benchmark`](https://github.com/rvhonorato/gdock-benchmark): repository
containing all scripts and raw data relevant to benchmarking the performance
of `gdock`
- [`gdock-website`](https://github.com/rvhonorato/gdock-website): source
code for [gdock.org](https://gdock.org)

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests
on [GitHub](https://github.com/rvhonorato/gdock).

Before submitting a pull request, please ensure:

- All tests pass (`cargo test`)
- Code is formatted (`cargo fmt`)
- Linting passes (`cargo clippy`)

## Citation

Coming soon.

## License

BSD Zero Clause License. See [LICENSE](LICENSE) file.
