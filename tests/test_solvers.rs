use klmc_ultimate::core::domain::{Params, AlgorithmType, Cluster, Species};
use klmc_ultimate::core::chemistry::InteractionGrid;
use klmc_ultimate::solvers::ga::GeneticAlgorithm;
use klmc_ultimate::solvers::bh::BasinHopping;
use klmc_ultimate::solvers::SolverEvent;
use crossbeam_channel::unbounded;
use std::sync::Arc;
use crate::common::MockEvaluator;

mod common;

#[test]
fn test_ga_flow() {
    let params = Params {
        algorithm: AlgorithmType::GeneticAlgorithm,
        atom_count: 4,
        atom_counts: vec![2, 2],
        population_size: 10,
        max_steps: 5,
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
    ga.solve(tx);

    let mut finished = false;
    let mut received_stats = false;

    for msg in rx {
        match msg {
            SolverEvent::Finished => finished = true,
            SolverEvent::GenerationUpdate(_) => received_stats = true,
            _ => {}
        }
    }

    assert!(finished, "GA did not finish");
    assert!(received_stats, "GA did not send stats");
}

#[test]
fn test_bh_flow() {
    let params = Params {
        algorithm: AlgorithmType::BasinHopping,
        atom_count: 4,
        atom_counts: vec![2, 2],
        max_steps: 5,
        ..Default::default()
    };

    let species = vec![
        Species { symbol: "A".into(), radius_covalent: 0.5, ..Default::default() },
        Species { symbol: "B".into(), radius_covalent: 0.5, ..Default::default() },
    ];
    let grid = Arc::new(InteractionGrid::new(&species, 0.5));
    let evaluator = Arc::new(MockEvaluator);

    let bh = BasinHopping::new(evaluator, grid.clone(), params.clone());

    let (tx, rx) = unbounded();
    let mut rng = rand::thread_rng();
    let start_cluster = Cluster::new_random(&params.atom_counts, params.box_size, &grid, &mut rng)
        .expect("Failed to create start cluster");

    bh.solve(start_cluster, tx);

    let mut finished = false;
    for msg in rx {
        match msg {
            SolverEvent::Finished => finished = true,
            _ => {}
        }
    }

    assert!(finished, "BH did not finish");
}
