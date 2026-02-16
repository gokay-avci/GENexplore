use std::collections::VecDeque;
use std::time::{Duration, Instant};
use crossbeam_channel::{Receiver, TryRecvError};
use ratatui::widgets::TableState;

use crate::core::domain::{Cluster, Params};
use crate::solvers::{SolverEvent, GenStats};

// --- Constants ---
const HISTORY_CAPACITY: usize = 1000;
const LOG_CAPACITY: usize = 200;
const HOF_CAPACITY: usize = 50;

// --- Enums ---

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AppMode {
    Dashboard,
    Analysis,
    HallOfFame,
    StructureViewer,
    Help,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkerStatus {
    Idle,
    Starting,
    Running,
    Paused,
    Finished,
    Error,
}

// --- Telemetry & Analytics ---

#[derive(Debug, Clone)]
pub struct Telemetry {
    // History Queues for Sparklines
    pub best_energy_history: VecDeque<(f64, f64)>, // (Iter, Energy)
    pub avg_energy_history: VecDeque<(f64, f64)>,
    pub diversity_history: VecDeque<(f64, f64)>,   // (Iter, Diversity %)
    pub mutation_history: VecDeque<(f64, f64)>,    // (Iter, Rate)
    
    // Global Bounds for Chart Scaling
    pub global_min_energy: f64,
    pub global_max_energy: f64,
}

impl Telemetry {
    pub fn new() -> Self {
        Self {
            best_energy_history: VecDeque::with_capacity(HISTORY_CAPACITY),
            avg_energy_history: VecDeque::with_capacity(HISTORY_CAPACITY),
            diversity_history: VecDeque::with_capacity(HISTORY_CAPACITY),
            mutation_history: VecDeque::with_capacity(HISTORY_CAPACITY),
            global_min_energy: f64::MAX,
            global_max_energy: f64::MIN,
        }
    }

    pub fn ingest(&mut self, stats: &GenStats) {
        // Enforce Capacity
        if self.best_energy_history.len() >= HISTORY_CAPACITY {
            self.best_energy_history.pop_front();
            self.avg_energy_history.pop_front();
            self.diversity_history.pop_front();
            self.mutation_history.pop_front();
        }

        // Update Bounds
        if stats.best_energy < self.global_min_energy { self.global_min_energy = stats.best_energy; }
        if stats.worst_energy > self.global_max_energy { self.global_max_energy = stats.worst_energy; }

        // Sanity check for first point to prevent chart crash on zero range
        if (self.global_max_energy - self.global_min_energy).abs() < 1e-6 {
            self.global_max_energy += 1.0;
            self.global_min_energy -= 1.0;
        }

        let x = stats.generation as f64;
        
        self.best_energy_history.push_back((x, stats.best_energy));
        self.avg_energy_history.push_back((x, stats.avg_energy));
        self.diversity_history.push_back((x, stats.diversity * 100.0));
        self.mutation_history.push_back((x, stats.mutation_rate));
    }
}

/// Visualization State for 3D Viewer
#[derive(Debug, Clone)]
pub struct Viewport {
    pub azimuth: f64,
    pub elevation: f64,
    pub zoom: f64,
    pub auto_rotate: bool,
    pub last_tick: Instant,
}

impl Viewport {
    pub fn new() -> Self {
        Self {
            azimuth: 0.0,
            elevation: 0.0,
            zoom: 1.0,
            auto_rotate: true,
            last_tick: Instant::now(),
        }
    }

    pub fn tick(&mut self) {
        if self.auto_rotate {
            let now = Instant::now();
            let dt = now.duration_since(self.last_tick).as_secs_f64();
            self.azimuth += 0.5 * dt; // Rotate 0.5 rad/s
            self.last_tick = now;
        } else {
            self.last_tick = Instant::now();
        }
    }
}

// --- The Master State ---

pub struct AppState {
    // System
    pub should_quit: bool,
    pub mode: AppMode,
    pub params: Params,
    
    // Worker
    pub rx: Option<Receiver<SolverEvent>>, 
    pub worker_status: WorkerStatus,
    
    // Simulation Data
    pub total_iterations: usize,
    pub start_time: Instant,
    pub current_best: Option<Cluster>,
    pub hall_of_fame: Vec<Cluster>, 
    pub active_cluster: Option<Cluster>, 
    
    // Analytics
    pub telemetry: Telemetry,
    pub logs: VecDeque<String>,
    
    // UI Elements
    pub hof_state: TableState,
    pub viewport: Viewport,
    
    // Performance Metrics
    pub ops_counter: usize,
    pub ops_per_second: f64,
    last_ops_check: Instant,
}

impl AppState {
    pub fn new(default_params: Params) -> Self {
        let mut hof_state = TableState::default();
        hof_state.select(Some(0));

        Self {
            should_quit: false,
            mode: AppMode::Dashboard,
            params: default_params,
            rx: None,
            worker_status: WorkerStatus::Idle,
            total_iterations: 0,
            start_time: Instant::now(),
            current_best: None,
            hall_of_fame: Vec::with_capacity(HOF_CAPACITY),
            active_cluster: None,
            telemetry: Telemetry::new(),
            logs: VecDeque::with_capacity(LOG_CAPACITY),
            hof_state,
            viewport: Viewport::new(),
            ops_counter: 0,
            ops_per_second: 0.0,
            last_ops_check: Instant::now(),
        }
    }

    pub fn set_channel(&mut self, rx: Receiver<SolverEvent>) {
        self.rx = Some(rx);
        self.worker_status = WorkerStatus::Starting;
        self.start_time = Instant::now();
    }

    pub fn tick(&mut self) {
        self.viewport.tick();
        self.calc_metrics();

        // Process Events
        if let Some(rx) = self.rx.clone() {
            for _ in 0..100 {
                match rx.try_recv() {
                    Ok(evt) => self.handle_event(evt),
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => {
                        self.worker_status = WorkerStatus::Finished;
                        self.log("Worker disconnected.");
                        self.rx = None;
                        break;
                    }
                }
            }
        }
    }

    fn handle_event(&mut self, event: SolverEvent) {
        match event {
            SolverEvent::Log(msg) => self.log(msg),
            
            SolverEvent::WorkerHeartbeat(ops) => {
                self.worker_status = WorkerStatus::Running;
                if ops > 0.0 { self.ops_per_second = ops; }
            },

            SolverEvent::GenerationUpdate(stats) => {
                self.worker_status = WorkerStatus::Running;
                self.total_iterations = stats.generation;
                self.ops_counter += stats.valid_count;
                self.telemetry.ingest(&stats);
            },

            SolverEvent::NewBest(cluster) => {
                self.handle_new_best(cluster);
            },

            SolverEvent::Finished => {
                self.worker_status = WorkerStatus::Finished;
                self.log("Solver finished.");
            }
        }
    }

    fn handle_new_best(&mut self, cluster: Cluster) {
        let e_new = cluster.energy.unwrap_or(0.0);
        
        // 1. Update Global Record
        let is_global = match &self.current_best {
            Some(curr) => e_new < curr.energy.unwrap_or(f64::MAX),
            None => true
        };

        if is_global {
            self.current_best = Some(cluster.clone());
            self.log(format!(">>> New Global Record: {:.5} eV", e_new));
        }

        // 2. Hall of Fame Deduplication (Isomer Check)
        let mut replaced = false;
        
        let new_hash = cluster.hash_key.as_deref().unwrap_or("INVALID");
        
        if new_hash != "INVALID" && new_hash != "NAN_COORDS" {
            for existing in &mut self.hall_of_fame {
                if let Some(ex_hash) = &existing.hash_key {
                    if ex_hash == new_hash {
                        // Same isomer found.
                        // If new energy is lower (better relaxed), replace it.
                        let e_old = existing.energy.unwrap_or(f64::MAX);
                        if e_new < e_old - 1e-5 {
                            *existing = cluster.clone();
                        }
                        replaced = true;
                        break; 
                    }
                }
            }
        }

        // If not a duplicate (or we allow duplicates due to bad hash), insert and sort
        if !replaced {
            self.hall_of_fame.push(cluster.clone());
            self.hall_of_fame.sort_by(|a, b| 
                a.energy.partial_cmp(&b.energy).unwrap_or(std::cmp::Ordering::Equal)
            );
            if self.hall_of_fame.len() > HOF_CAPACITY {
                self.hall_of_fame.truncate(HOF_CAPACITY);
            }
        }
        
        // Auto-select for visualization
        self.active_cluster = Some(cluster);
    }

    /// Logs a message to the internal buffer.
    fn log(&mut self, msg: impl Into<String>) {
        if self.logs.len() >= LOG_CAPACITY {
            self.logs.pop_front();
        }
        self.logs.push_back(msg.into());
    }

    fn calc_metrics(&mut self) {
        let now = Instant::now();
        if now.duration_since(self.last_ops_check) >= Duration::from_secs(1) {
            if self.ops_per_second == 0.0 && self.ops_counter > 0 {
                self.ops_per_second = self.ops_counter as f64;
            }
            self.ops_counter = 0;
            self.last_ops_check = now;
        }
    }

    // --- Input Handling ---
    
    pub fn on_key(&mut self, key: char) {
        match key {
            'q' => self.should_quit = true,
            '1' => self.mode = AppMode::Dashboard,
            '2' => self.mode = AppMode::Analysis,
            '3' => self.mode = AppMode::HallOfFame,
            '4' => self.mode = AppMode::StructureViewer,
            ' ' => self.toggle_pause(),
            'r' => self.viewport.azimuth = 0.0,
            'j' => self.select_next_hof(),
            'k' => self.select_prev_hof(),
            _ => {}
        }
    }

    fn select_next_hof(&mut self) {
        if self.hall_of_fame.is_empty() { return; }
        let i = match self.hof_state.selected() {
            Some(i) => if i >= self.hall_of_fame.len() - 1 { 0 } else { i + 1 },
            None => 0,
        };
        self.hof_state.select(Some(i));
        self.active_cluster = Some(self.hall_of_fame[i].clone());
    }

    fn select_prev_hof(&mut self) {
        if self.hall_of_fame.is_empty() { return; }
        let i = match self.hof_state.selected() {
            Some(i) => if i == 0 { self.hall_of_fame.len() - 1 } else { i - 1 },
            None => 0,
        };
        self.hof_state.select(Some(i));
        self.active_cluster = Some(self.hall_of_fame[i].clone());
    }

    pub fn toggle_pause(&mut self) {
        match self.worker_status {
            WorkerStatus::Running => {
                self.worker_status = WorkerStatus::Paused;
                self.log("Paused.");
            },
            WorkerStatus::Paused => {
                self.worker_status = WorkerStatus::Running;
                self.log("Resumed.");
            },
            _ => {}
        }
    }
}