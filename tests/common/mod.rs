use anyhow::Result;
use klmc_ultimate::core::domain::Cluster;
use klmc_ultimate::engine::evaluator::{EvaluationResult, Evaluator};

pub struct MockEvaluator;

impl Evaluator for MockEvaluator {
    fn evaluate(&self, cluster: &Cluster) -> Result<EvaluationResult> {
        // Simple mock energy function: Sum of distances from origin
        // This favors compact clusters near origin
        let mut energy = 0.0;
        for atom in &cluster.atoms {
            energy += atom.position.coords.norm();
        }

        Ok(EvaluationResult {
            energy,
            gradient_norm: Some(0.1),
            relaxed_cluster: Some(cluster.clone()),
        })
    }

    fn name(&self) -> &str {
        "Mock Evaluator"
    }
}
