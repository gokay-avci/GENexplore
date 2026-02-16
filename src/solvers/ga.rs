use std::collections::HashSet;
use std::sync::{Arc, Mutex};
use std::time::Instant;
use crossbeam_channel::Sender;
use rand::prelude::*;
use rayon::prelude::*;

use crate::core::domain::{Cluster, Params, ClusterStatus};
use crate::core::spatial;
use crate::core::chemistry::InteractionGrid;
use crate::engine::evaluator::Evaluator;
use crate::engine::operators::{Mutator, crossover_cut_splice};
use crate::analysis::topology;
use crate::solvers::{SolverEvent, GenStats};

pub struct GeneticAlgorithm {
    evaluator: Arc<dyn Evaluator>,
    grid: Arc<InteractionGrid>,
    params: Params,
}

impl GeneticAlgorithm {
    pub fn new(
        evaluator: Arc<dyn Evaluator>,
        grid: Arc<InteractionGrid>,
        params: Params
    ) -> Self {
        Self { evaluator, grid, params }
    }

    pub fn solve(&self, tx: Sender<SolverEvent>) {
        let _start_time = Instant::now();
        let mut rng = rand::thread_rng();

        // 1. Initialization Phase
        let _ = tx.send(SolverEvent::Log("Initializing Population...".to_string()));
        
        let mut population = self.generate_initial_population();
        
        if population.is_empty() {
            let _ = tx.send(SolverEvent::Log("CRITICAL: Failed to generate valid initial population.".to_string()));
            let _ = tx.send(SolverEvent::Finished);
            return;
        }

        self.evaluate_batch(&mut population);
        self.rank_population(&mut population);
        
        if let Some(best) = population.first() {
            if best.energy.is_some() {
                let _ = tx.send(SolverEvent::NewBest(best.clone()));
            }
        }

        // State Tracking
        let mut stagnation_counter = 0;
        let mut extinction_cooldown = 0;
        let mut last_global_best_e = population.first().and_then(|c| c.energy).unwrap_or(f64::MAX);
        let mut total_evals = 0;
        
        // Dynamic Parameters
        let mut current_mutation_rate = self.params.mutation_rate;

        // 2. Evolution Loop
        for gen in 1..=self.params.max_steps {
            let gen_start = Instant::now();

            // A. Elitism
            let mut next_gen = Vec::with_capacity(self.params.population_size);
            let elites = population.iter()
                .take(self.params.elitism_count)
                .cloned()
                .collect::<Vec<_>>();
            next_gen.extend(elites);

            // B. Breeding
            if !population.is_empty() {
                while next_gen.len() < self.params.population_size {
                    let p1 = self.tournament_select(&population, &mut rng);
                    let p2 = self.tournament_select(&population, &mut rng);

                    // Crossover
                    let mut child = if rng.gen::<f64>() < self.params.crossover_rate {
                        crossover_cut_splice(p1, p2, &mut rng).unwrap_or_else(|| p1.clone())
                    } else {
                        p1.clone()
                    };

                    // Mutation
                    if rng.gen::<f64>() < current_mutation_rate {
                        let rattle_mag = if stagnation_counter > 20 { 0.3 } else { 0.1 };
                        
                        let mut mutator = Mutator::new()
                            .rotate(0.5)
                            .rattle(rattle_mag)
                            .swap(1);
                        
                        if rng.gen_bool(0.2) {
                            mutator = mutator.breathing(0.05);
                        }
                        
                        child = mutator.apply(&child, &mut rng);
                    }

                    if spatial::check_overlap(&child, &self.grid) {
                        child.status = ClusterStatus::Born;
                        child.generation = gen as u64;
                        next_gen.push(child);
                    }
                }
            }

            // C. Evaluation
            let evals_this_gen = self.evaluate_batch(&mut next_gen);
            total_evals += evals_this_gen;

            // D. Topology & Diversity
            for c in &mut next_gen {
                if c.energy.is_some() {
                    // Use updated topology with PMOI + Graph Spectrum
                    c.hash_key = Some(topology::generate_hash_key(c, 1.5)); 
                }
            }

            // Deduplicate (Remove Isomers)
            let (mut unique_pop, diversity) = self.deduplicate_population(next_gen);

            // --- SMART REFILL STRATEGY ---
            // If deduplication removed individuals, fill the gap with Mutated Survivors
            // instead of random trash. This maintains energy quality while forcing geometric diversity.
            let target_size = self.params.population_size;
            
            if unique_pop.is_empty() {
                // Catastrophic collapse (should not happen with elitism, but safe fallback)
                unique_pop = self.generate_initial_population();
                self.evaluate_batch(&mut unique_pop);
            } else if unique_pop.len() < target_size {
                let needed = target_size - unique_pop.len();
                let mut refill = Vec::with_capacity(needed);
                
                // Cycle through survivors (Best -> Worst -> Best...)
                let source_pool = unique_pop.clone();
                let mut source_iter = source_pool.iter().cycle();

                while refill.len() < needed {
                    if let Some(parent) = source_iter.next() {
                        // Apply HEAVY mutation to force it into a new topological basin
                        // Twist + Rotate + Rattle
                        let mut child = Mutator::new()
                            .rotate(3.14)   // Full rotation potential
                            .twist(0.5)     // Significant twist
                            .rattle(0.2)    // Shake atoms
                            .apply(parent, &mut rng);
                        
                        child.origin = "Refill".to_string();
                        child.status = ClusterStatus::Born;
                        child.energy = None; // Force re-eval
                        refill.push(child);
                    }
                }

                // Evaluate the refill batch
                self.evaluate_batch(&mut refill);
                
                // Calculate hashes for refill to ensure they are tracked correctly next gen
                for c in &mut refill {
                    if c.energy.is_some() {
                        c.hash_key = Some(topology::generate_hash_key(c, 1.5));
                    }
                }
                
                unique_pop.extend(refill);
            }

            population = unique_pop;

            // E. Stagnation Logic
            self.rank_population(&mut population);
            
            let current_best_e = population.first().and_then(|c| c.energy).unwrap_or(f64::MAX);

            if current_best_e < last_global_best_e - 1e-5 {
                stagnation_counter = 0;
                extinction_cooldown = 0;
                current_mutation_rate = self.params.mutation_rate;
                last_global_best_e = current_best_e;
                if let Some(best) = population.first() {
                    let _ = tx.send(SolverEvent::NewBest(best.clone()));
                }
            } else {
                stagnation_counter += 1;
            }

            // Adaptive State Machine
            if extinction_cooldown > 0 {
                extinction_cooldown -= 1;
            } else {
                let catastrophic_stagnation = stagnation_counter > 50;
                let premature_convergence = gen > 20 && stagnation_counter > 20 && diversity < 0.1;

                if catastrophic_stagnation || premature_convergence {
                    let reason = if catastrophic_stagnation { "Stagnation 50+" } else { "Low Diversity" };
                    let _ = tx.send(SolverEvent::Log(format!("Mass Extinction ({}) -> Reseeding", reason)));
                    
                    let keep = self.params.elitism_count;
                    if population.len() > keep {
                        population.truncate(keep);
                    }
                    
                    // Fill with TRUE randoms to reset the gene pool completely
                    let mut attempts = 0;
                    while population.len() < self.params.population_size && attempts < self.params.population_size * 100 {
                        attempts += 1;
                        if let Some(mut r) = Cluster::new_random(
                            &self.params.atom_counts, 
                            self.params.box_size, 
                            &self.grid, 
                            &mut rng
                        ) {
                            if let Ok(res) = self.evaluator.evaluate(&r) {
                                r.energy = Some(res.energy);
                                r.status = ClusterStatus::Evaluated;
                                population.push(r);
                            }
                        }
                    }
                    self.rank_population(&mut population);
                    
                    stagnation_counter = 0;
                    extinction_cooldown = 50;
                    current_mutation_rate = self.params.mutation_rate;

                } else if stagnation_counter > 20 {
                    if current_mutation_rate < 0.5 {
                        let _ = tx.send(SolverEvent::Log("Stagnation (20+) -> Hyper-Mutation".to_string()));
                        current_mutation_rate = 0.5;
                    }
                }
            }

            // F. Telemetry
            let valid_clusters: Vec<&Cluster> = population.iter().filter(|c| c.energy.is_some()).collect();
            let valid_count = valid_clusters.len();
            
            let (best_e, worst_e, avg_e) = if valid_count > 0 {
                let best = valid_clusters.first().unwrap().energy.unwrap();
                let worst = valid_clusters.last().unwrap().energy.unwrap();
                let sum: f64 = valid_clusters.iter().map(|c| c.energy.unwrap()).sum();
                (best, worst, sum / valid_count as f64)
            } else {
                (0.0, 0.0, 0.0)
            };

            let _ = tx.send(SolverEvent::GenerationUpdate(GenStats {
                generation: gen,
                best_energy: best_e,
                avg_energy: avg_e,
                worst_energy: worst_e,
                diversity,
                valid_count,
                pop_size: population.len(),
                mutation_rate: current_mutation_rate,
            }));

            let duration = gen_start.elapsed().as_secs_f64();
            if duration > 0.0 {
                let ops = evals_this_gen as f64 / duration;
                let _ = tx.send(SolverEvent::WorkerHeartbeat(ops));
            }
        }

        let _ = tx.send(SolverEvent::Log(format!("GA Finished. Total Evals: {}", total_evals)));
        let _ = tx.send(SolverEvent::Finished);
    }

    // --- Helpers ---

    fn generate_initial_population(&self) -> Vec<Cluster> {
        let mut pop = Vec::new();
        let mut rng = rand::thread_rng();
        let attempts = self.params.population_size * 50;
        
        for _ in 0..attempts {
            if pop.len() >= self.params.population_size { break; }
            if let Some(c) = Cluster::new_random(
                &self.params.atom_counts, 
                self.params.box_size, 
                &self.grid, 
                &mut rng
            ) {
                pop.push(c);
            }
        }
        pop
    }

    fn evaluate_batch(&self, pop: &mut [Cluster]) -> usize {
        let eval_ref = &self.evaluator;
        let count = Arc::new(Mutex::new(0));

        pop.par_iter_mut()
            .filter(|c| c.status == ClusterStatus::Born)
            .for_each(|cluster| {
                match eval_ref.evaluate(cluster) {
                    Ok(res) => {
                        cluster.energy = Some(res.energy);
                        
                        if let Some(geom) = res.relaxed_cluster {
                            if geom.atoms.len() == cluster.atoms.len() {
                                for (orig, new) in cluster.atoms.iter_mut().zip(geom.atoms.iter()) {
                                    orig.position = new.position;
                                }
                                if geom.lattice.is_some() {
                                    cluster.lattice = geom.lattice;
                                }
                                spatial::wrap_or_center(cluster);
                            } else {
                                cluster.status = ClusterStatus::Discarded;
                                cluster.energy = None;
                                return;
                            }
                        }
                        
                        cluster.status = ClusterStatus::Evaluated;
                        if let Ok(mut c) = count.lock() { *c += 1; }
                    },
                    Err(_) => {
                        cluster.status = ClusterStatus::Discarded;
                        cluster.energy = None;
                    }
                }
            });
        
        let final_count = *count.lock().unwrap();
        final_count
    }

    fn tournament_select<'a>(&self, pop: &'a [Cluster], rng: &mut impl Rng) -> &'a Cluster {
        if pop.is_empty() {
            panic!("Tournament selection called on empty population");
        }

        let mut best = &pop[rng.gen_range(0..pop.len())];
        let mut best_e = best.energy.unwrap_or(f64::MAX);

        for _ in 0..1 {
            let candidate = &pop[rng.gen_range(0..pop.len())];
            let cand_e = candidate.energy.unwrap_or(f64::MAX);
            if cand_e < best_e {
                best = candidate;
                best_e = cand_e;
            }
        }
        best
    }

    fn rank_population(&self, pop: &mut Vec<Cluster>) {
        pop.sort_by(|a, b| {
            match (a.energy, b.energy) {
                (Some(ea), Some(eb)) => ea.partial_cmp(&eb).unwrap_or(std::cmp::Ordering::Equal),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => std::cmp::Ordering::Equal,
            }
        });
    }

    fn deduplicate_population(&self, pop: Vec<Cluster>) -> (Vec<Cluster>, f64) {
        let initial_count = pop.len();
        if initial_count == 0 { return (pop, 0.0); }

        let mut unique = Vec::new();
        let mut seen_hashes = HashSet::new();

        for c in pop {
            if c.status == ClusterStatus::Discarded || c.energy.is_none() { continue; }
            
            if let Some(hash) = &c.hash_key {
                if hash == "INVALID" || hash.contains("NAN") { 
                    unique.push(c);
                    continue; 
                }
                
                if !seen_hashes.contains(hash) {
                    seen_hashes.insert(hash.clone());
                    unique.push(c);
                }
            } else {
                unique.push(c);
            }
        }

        let diversity = if initial_count > 0 { 
            unique.len() as f64 / initial_count as f64 
        } else { 0.0 };
        
        (unique, diversity)
    }
}