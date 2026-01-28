---
title: 'gdock: Information-driven protein-protein docking using genetic
  algorithms'
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

`gdock` is a protein-protein docking tool that uses a genetic algorithm
combined with physics-based scoring to predict the structure of protein
complexes. It is information-driven, meaning it uses experimental restraints
(Ambiguous Interaction Restraints, AIRs) to guide the docking process toward
biologically relevant solutions.

# Statement of need

Protein-protein interactions are fundamental to virtually all biological
processes. Determining the three-dimensional structure of protein complexes is
crucial for understanding molecular mechanisms and designing therapeutics.
While experimental methods like X-ray crystallography and cryo-EM can resolve
complex structures, they are time-consuming and not always feasible.
Computational docking provides a complementary approach.

`gdock` addresses the need for a fast, accessible docking tool that can
incorporate experimental data to guide predictions.

# State of the field

Several protein-protein docking tools exist, including HADDOCK
[@dominguez2003haddock], ClusPro [@kozakov2017cluspro], and ZDOCK
[@pierce2014zdock]. `gdock` was developed as a lightweight, fast alternative
that...

# Software design

`gdock` is implemented in Rust for performance and memory safety. The software
consists of several key components:

- **Genetic Algorithm**: Explores the conformational space through selection,
  crossover, and mutation
- **Energy Function**: Combines van der Waals, electrostatic, desolvation, and
  restraint terms
- **Clustering**: Uses Fraction of Common Contacts (FCC) to group similar
  solutions

The scoring function combines multiple energy terms:

$$E_{total} = w_{vdw} E_{vdw} + w_{elec} E_{elec} + w_{desolv} E_{desolv} + w_{air} E_{air}$$

where the weights were optimized through benchmarking on standard datasets.

# Research impact statement

# AI usage disclosure

# Acknowledgements

# References
