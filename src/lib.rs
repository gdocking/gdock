pub mod chromosome;
pub mod clustering;
pub mod constants;
pub mod evaluator;
pub mod fitness;
pub mod ga;
pub mod hall_of_fame;
pub mod population;
pub mod restraints;
pub mod runner;
pub mod structure;
pub mod toppar;
pub mod utils;

// Modules that use std::fs / CLI deps — not available on wasm32
#[cfg(not(target_arch = "wasm32"))]
pub mod commands;
#[cfg(not(target_arch = "wasm32"))]
pub mod reporting;
#[cfg(not(target_arch = "wasm32"))]
pub mod scoring;
