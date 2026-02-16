use std::sync::Arc;
use std::time::Instant;
use crossbeam_channel::Sender;
use rand::Rng;

use crate::core::domain::{Cluster, Params, ClusterStatus};
use crate::engine::evaluator::Evaluator;
use crate::engine::operators::Mutator;
use crate::core::spatial;
use crate::core::chemistry::InteractionGrid;
use crate::solvers::{SolverEvent, GenStats};

pub struct BasinHopping {
    evaluator: Arc<dyn Evaluator>,
    grid: Arc<InteractionGrid>,
    params: Params,
}

impl BasinHopping {
    pub fn new(
        evaluator: Arc<dyn Evaluator>,
        grid: Arc<InteractionGrid>,
        params: Params
    ) -> Self {
        Self {
            evaluator,
            grid,
            params,
        }
    }

    /// Runs the Basin Hopping loop (Monte Carlo Minimization).
    /// Tracks a single "Walker" cluster across the energy landscape.
    pub fn solve(&self, mut current: Cluster, tx: Sender<SolverEvent>) {
        let mut rng = rand::thread_rng();
        let kb_ev = 8.617333262e-5; // Boltzmann constant
        
        // Defensive: Validate inputs
        if self.params.bh_steps == 0 {
            let _ = tx.send(SolverEvent::Log("BH steps set to 0. Exiting.".to_string()));
            let _ = tx.send(SolverEvent::Finished);
            return;
        }

        // 1. Initial Relaxation
        if current.energy.is_none() {
            let _ = tx.send(SolverEvent::Log("Relaxing initial structure...".to_string()));
            match self.evaluator.evaluate(&current) {
                Ok(res) => {
                    current.energy = Some(res.energy);
                    
                    // Safe Geometry Update
                    if let Some(geom) = res.relaxed_cluster {
                        if geom.atoms.len() == current.atoms.len() {
                            for (orig, new) in current.atoms.iter_mut().zip(geom.atoms.iter()) {
                                orig.position = new.position;
                            }
                            if geom.lattice.is_some() { current.lattice = geom.lattice; }
                            spatial::wrap_or_center(&mut current);
                        }
                    }
                    current.status = ClusterStatus::Evaluated;
                },
                Err(e) => {
                    let _ = tx.send(SolverEvent::Log(format!("Initial relaxation failed: {}", e)));
                    let _ = tx.send(SolverEvent::Finished);
                    return;
                }
            }
        }

        let mut best = current.clone();
        if let Some(_) = best.energy {
            let _ = tx.send(SolverEvent::NewBest(best.clone()));
        }
        
        let start_time = Instant::now();
        let mut accepted_count = 0;

        // 2. Main Loop
        for i in 1..=self.params.bh_steps {
            // A. Perturb
            // Standard BH move: Random translation + slight rotation to escape shallow wells
            let mut trial = Mutator::new()
                .translate(self.params.step_size)
                .rotate(0.2) 
                .apply(&current, &mut rng);
            
            trial.origin = format!("BH_{}", i);

            // B. Pre-check Geometry
            // If the move creates an overlap (collision), reject immediately (infinite energy)
            if !spatial::check_overlap(&trial, &self.grid) {
                // Report current state (Rejection)
                self.report_step(&tx, i, &current);
                continue;
            }

            // C. Local Minimization
            match self.evaluator.evaluate(&trial) {
                Ok(res) => {
                    // Update trial with relaxed energy/geometry
                    trial.energy = Some(res.energy);
                    
                    if let Some(geom) = res.relaxed_cluster {
                        if geom.atoms.len() == trial.atoms.len() {
                            for (orig, new) in trial.atoms.iter_mut().zip(geom.atoms.iter()) {
                                orig.position = new.position;
                            }
                            if geom.lattice.is_some() { trial.lattice = geom.lattice; }
                            spatial::wrap_or_center(&mut trial);
                        } else {
                            // Atom count mismatch from engine -> Reject
                            self.report_step(&tx, i, &current);
                            continue;
                        }
                    }
                    trial.status = ClusterStatus::Evaluated;

                    // D. Metropolis Acceptance
                    let e_new = res.energy;
                    let e_old = current.energy.unwrap_or(f64::MAX);

                    let accepted = if e_new < e_old {
                        true
                    } else {
                        let delta = e_new - e_old;
                        let temp = self.params.temperature;
                        if temp <= 1e-9 {
                            false // Quench only
                        } else {
                            let prob = (-delta / (kb_ev * temp)).exp();
                            rng.gen::<f64>() < prob
                        }
                    };

                    if accepted {
                        accepted_count += 1;
                        current = trial; // Move walker

                        // Check Global Best
                        if let Some(best_e) = best.energy {
                            if e_new < best_e {
                                best = current.clone();
                                let _ = tx.send(SolverEvent::NewBest(best.clone()));
                            }
                        }
                    }

                    // Report stats (Current position of walker)
                    self.report_step(&tx, i, &current);
                },
                Err(_) => {
                    // Physics engine failed (e.g. SCF did not converge) -> Reject move
                    self.report_step(&tx, i, &current);
                }
            }
        }

        let duration = start_time.elapsed().as_secs_f64();
        let rate = if duration > 0.0 { self.params.bh_steps as f64 / duration } else { 0.0 };
        
        let _ = tx.send(SolverEvent::Log(format!("BH Finished. Acceptance: {}/{}", accepted_count, self.params.bh_steps)));
        let _ = tx.send(SolverEvent::WorkerHeartbeat(rate));
        let _ = tx.send(SolverEvent::Finished);
    }

    fn report_step(&self, tx: &Sender<SolverEvent>, iter: usize, cluster: &Cluster) {
        let e = cluster.energy.unwrap_or(0.0);
        
        // Map single walker to population stats
        let stats = GenStats {
            generation: iter,
            best_energy: e,
            avg_energy: e,
            worst_energy: e,
            diversity: 1.0, // A population of 1 is always 100% diverse relative to itself
            valid_count: 1,
            pop_size: 1,
            mutation_rate: 0.0, // BH uses fixed step_size, not mutation rate
        };

        let _ = tx.send(SolverEvent::GenerationUpdate(stats));
    }
}