use crate::core::domain::Cluster;
use anyhow::Result;

/// The result of a physical evaluation.
#[derive(Debug, Clone)]
pub struct EvaluationResult {
    /// The potential energy (eV).
    pub energy: f64,
    /// The gradient norm (eV/Å) - useful for checking convergence quality.
    pub gradient_norm: Option<f64>,
    /// The updated cluster with relaxed coordinates.
    /// Returns None if the geometry exploded or failed.
    pub relaxed_cluster: Option<Cluster>,
}

/// A generic interface for physics engines.
/// Implementations must be Thread-Safe (Sync).
pub trait Evaluator: Send + Sync {
    /// Takes a raw cluster, relaxes it, and returns the result.
    fn evaluate(&self, cluster: &Cluster) -> Result<EvaluationResult>;

    /// Returns the name of the engine (e.g., "GULP 6.1").
    fn name(&self) -> &str;
}
