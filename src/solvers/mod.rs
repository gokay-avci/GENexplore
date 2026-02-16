use crate::core::domain::Cluster;

/// Detailed statistics for a single generation/step.
/// Used for telemetry and UI visualization.
#[derive(Debug, Clone)]
pub struct GenStats {
    pub generation: usize,
    pub best_energy: f64,
    pub avg_energy: f64,
    pub worst_energy: f64,
    pub diversity: f64,     // 0.0 to 1.0 (Unique Isomers / Population Size)
    pub valid_count: usize, // Number of structures that passed geometry + GULP checks
    pub pop_size: usize,    // Total active population size

    /// The current mutation rate being applied.
    /// Allows the UI to visualize adaptive strategies (e.g. Hyper-Mutation spikes).
    pub mutation_rate: f64,
}

impl Default for GenStats {
    fn default() -> Self {
        Self {
            generation: 0,
            best_energy: 0.0,
            avg_energy: 0.0,
            worst_energy: 0.0,
            diversity: 0.0,
            valid_count: 0,
            pop_size: 0,
            mutation_rate: 0.0,
        }
    }
}

/// Events emitted by solvers to the main thread.
#[derive(Debug, Clone)]
pub enum SolverEvent {
    /// Diagnostic log message.
    Log(String),

    /// High-level heartbeat for async diagnostics (e.g., "Worker is alive").
    /// Contains (ops_per_sec estimate from worker perspective).
    WorkerHeartbeat(f64),

    /// A completed generation or basin-hopping step with full statistics.
    GenerationUpdate(GenStats),

    /// A structure that beats the current global best (Energy Record).
    NewBest(Cluster),

    /// Solver has finished its run.
    Finished,
}

pub mod bh;
pub mod ga;
