use klmc_ultimate::core::domain::{Params, AlgorithmType, Species, Cluster};
use klmc_ultimate::core::chemistry::InteractionGrid;
use klmc_ultimate::solvers::ga::GeneticAlgorithm;
use klmc_ultimate::solvers::SolverEvent;
use klmc_ultimate::engine::evaluator::{Evaluator, EvaluationResult};
use crossbeam_channel::unbounded;
use std::sync::Arc;
use std::time::Instant;
use anyhow::Result;

pub struct MockEvaluator;

impl Evaluator for MockEvaluator {
    fn evaluate(&self, cluster: &Cluster) -> Result<EvaluationResult> {
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

    fn name(&self) -> &str { "Mock Evaluator" }
}

fn main() {
    let params = Params {
        algorithm: AlgorithmType::GeneticAlgorithm,
        atom_count: 50,
        atom_counts: vec![25, 25],
        population_size: 100,
        max_steps: 100,
        elitism_count: 10,
        ..Default::default()
    };

    let species = vec![
        Species { symbol: "A".into(), radius_covalent: 0.5, ..Default::default() },
        Species { symbol: "B".into(), radius_covalent: 0.5, ..Default::default() },
    ];
    let grid = Arc::new(InteractionGrid::new(&species, 0.5));
    let evaluator = Arc::new(MockEvaluator);

    let ga = GeneticAlgorithm::new(evaluator, grid, params);

    let (tx, rx) = unbounded();

    let start = Instant::now();
    ga.solve(tx);

    for msg in rx {
        if let SolverEvent::Finished = msg {
            break;
        }
    }
    let duration = start.elapsed();

    println!("GA took: {:?}", duration);
}
