# Contributing to gdock

Thank you for your interest in contributing to gdock! This document provides
guidelines to help you get started.

## Code of Conduct

This project follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).
By participating, you are expected to uphold this code. Please report
unacceptable behavior to <email@gdock.org>.

## How to Contribute

### Reporting Bugs

Open an [issue](https://github.com/rvhonorato/gdock/issues) with:

- A clear description of the problem
- Steps to reproduce the behavior
- Input files (PDB and restraints) if applicable
- Expected vs actual output
- Rust version (`rustc --version`) and OS

### Suggesting Features

Open an [issue](https://github.com/rvhonorato/gdock/issues) describing the
feature, its use case, and how it fits into the existing docking workflow.

### Submitting Pull Requests

1. Fork the repository and create a branch from `main`
2. Make your changes
3. Ensure all checks pass (see below)
4. Open a pull request against `main`

## Development Setup

```bash
git clone https://github.com/<your-username>/gdock.git
cd gdock
cargo build
```

Requires [Rust](https://www.rust-lang.org/tools/install) 1.70 or later.

## Quality Checks

All pull requests must pass the CI checks. Run these locally before submitting:

```bash
# Run the test suite
cargo test

# Check formatting
cargo fmt --all -- --check

# Run linting
cargo clippy -- -D warnings
```

To auto-format your code:

```bash
cargo fmt
```

## Project Structure

| Directory | Contents |
|-----------|----------|
| `src/` | Rust source code |
| `data/` | Test PDB files |
| `paper/` | JOSS paper draft |
| `scripts/` | Utility scripts |

See the [README](README.md) for a description of the subcommands and usage.

## License

By contributing, you agree that your contributions will be licensed under the
[BSD Zero Clause License](LICENSE).
