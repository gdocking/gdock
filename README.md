# gdock

**Information-driven protein-protein docking using a genetic algorithm**

[![ci](https://github.com/rvhonorato/gdock/actions/workflows/ci.yml/badge.svg)](https://github.com/rvhonorato/gdock/actions/workflows/ci.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/23671025da8a4334a754d8d5af76a34a)](https://app.codacy.com/gh/rvhonorato/gdock/dashboard)

![Crates.io License](https://img.shields.io/crates/l/gdock)
![Crates.io Version](https://img.shields.io/crates/v/gdock)
![Crates.io Total Downloads](https://img.shields.io/crates/d/gdock)

<img src="imgs/gdock_logo.png" width="350">

gdock is a fast protein-protein docking tool written in Rust that uses
restraints and energy components to guide the docking process. It combines a
genetic algorithm with physics-based scoring to find optimal protein-protein
complexes.

## Features

- **Fast**: Genetic algorithm with early stopping and elitism
- **Information-driven**: Uses residue restraints to guide docking
- **Flexible scoring**: Configurable energy weights (VDW, electrostatics,
  desolvation, restraints)
- **Quality metrics**: Optional DockQ calculation when reference structure is
  provided
- **Clustering**: FCC-based clustering to group similar solutions

## Installation

### From crates.io

```bash
cargo install gdock
```

### From GitHub releases

Download pre-built binaries from the
[releases page](https://github.com/rvhonorato/gdock/releases).

### Build from source

```bash
git clone https://github.com/rvhonorato/gdock
cd gdock
cargo build --release
```

The binary will be available at `./target/release/gdock`.

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

## Algorithm

gdock uses:

- **Genetic Algorithm**: Population of 150, elitism (top 5), tournament
  selection
- **Energy Function**: VDW + Electrostatics + Desolvation + AIR restraints
- **Restraints**: Flat-bottom potential (0-7 Angstrom) for specified residue pairs
- **Early Stopping**: Converges when no improvement for 10 generations
- **Clustering**: FCC-based clustering of final population

## Example

Using the test data included in the repository:

```bash
gdock run \
  --receptor data/2oob_A.pdb \
  --ligand data/2oob_B.pdb \
  --restraints 933:6,936:8,940:42,941:44,946:45,950:46 \
  --reference data/2oob.pdb \
  --output-dir example/
```

## Relevant repositories

- [`gdock-benchmark`](https://github.com/rvhonorato/gdock-benchmark): repository
containing all scripts and raw data relevant to benchmarking the performance
of `gdock`
- [`gdock-website`](https://github.com/rvhonorato/gdock-website): source
code for [gdock.org](https://gdock.org)

## Citation

Coming soon.

## License

BSD Zero Clause License. See [LICENSE](LICENSE) file.
