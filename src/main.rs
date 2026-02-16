use std::error::Error;
use std::io;
use std::panic;
use std::process::Command;
use std::sync::Arc;
use std::thread;
use std::time::{Duration, Instant};

use anyhow::{anyhow, Context, Result};
use clap::Parser;
use crossbeam_channel::unbounded;
use crossterm::{
    event::{self, DisableMouseCapture, EnableMouseCapture, Event, KeyCode},
    execute,
    terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen},
};
use ratatui::{backend::CrosstermBackend, Terminal};

use klmc_ultimate::core::chemistry::InteractionGrid;
use klmc_ultimate::core::domain::{AlgorithmType, Cluster, Params, Species, SystemDefinition};
use klmc_ultimate::engine::evaluator::Evaluator;
use klmc_ultimate::engine::external::gulp::GulpEvaluator;
use klmc_ultimate::interface::state::AppState;
use klmc_ultimate::interface::ui;
use klmc_ultimate::solvers::bh::BasinHopping;
use klmc_ultimate::solvers::ga::GeneticAlgorithm;

// --- CLI Definitions ---

#[derive(Parser, Debug)]
#[command(author, version, about = "KLMC Ultimate: Knowledge-Led Master Code (Rust Port)", long_about = None)]
struct Args {
    /// Number of worker threads for parallel evaluation
    #[arg(short, long, default_value_t = 4)]
    threads: usize,

    /// Number of atoms in the cluster
    #[arg(short, long, default_value_t = 12)]
    atoms: usize,

    /// Algorithm to run (ga, bh)
    #[arg(short, long, default_value = "ga")]
    algo: String,

    /// Initial box size (Angstroms)
    #[arg(short, long, default_value_t = 6.0)]
    box_size: f64,
}

// --- Terminal Guard (RAII) ---

struct TuiContext {
    terminal: Terminal<CrosstermBackend<io::Stdout>>,
}

impl TuiContext {
    fn new() -> Result<Self> {
        enable_raw_mode().context("Failed to enable raw mode")?;
        let mut stdout = io::stdout();
        execute!(stdout, EnterAlternateScreen, EnableMouseCapture)
            .context("Failed to setup terminal alternate screen")?;
        let backend = CrosstermBackend::new(stdout);
        let terminal = Terminal::new(backend).context("Failed to create terminal backend")?;
        Ok(Self { terminal })
    }
}

impl Drop for TuiContext {
    fn drop(&mut self) {
        // Best-effort restoration of terminal state
        let _ = disable_raw_mode();
        let _ = execute!(
            self.terminal.backend_mut(),
            LeaveAlternateScreen,
            DisableMouseCapture
        );
        let _ = self.terminal.show_cursor();
    }
}

// --- Initialization Helpers ---

fn setup_panic_hook() {
    let original_hook = panic::take_hook();
    panic::set_hook(Box::new(move |panic_info| {
        // Forcefully restore terminal before printing panic
        let _ = disable_raw_mode();
        let _ = execute!(io::stdout(), LeaveAlternateScreen, DisableMouseCapture);
        original_hook(panic_info);
    }));
}

fn create_default_system(args: &Args) -> SystemDefinition {
    // Define MgO system
    // Index 0 = Mg
    let mg = Species {
        symbol: "Mg".to_string(),
        atomic_number: 12,
        mass: 24.305,
        charge: 2.0,
        radius_covalent: 1.30,
        radius_ionic: 0.72,
        color_rgb: (0, 255, 255), // Cyan
    };

    // Index 1 = O
    let o = Species {
        symbol: "O".to_string(),
        atomic_number: 8,
        mass: 15.999,
        charge: -2.0,
        radius_covalent: 0.73,
        radius_ionic: 1.40,
        color_rgb: (255, 0, 0), // Red
    };

    let algo = match args.algo.to_lowercase().as_str() {
        "bh" => AlgorithmType::BasinHopping,
        "scan" => AlgorithmType::ScanBox,
        _ => AlgorithmType::GeneticAlgorithm,
    };

    // Stoichiometry Setup: 50/50 split for MgO
    let n_mg = args.atoms / 2;
    let n_o = args.atoms - n_mg; // Handle odd numbers by giving O one extra
    let atom_counts = vec![n_mg, n_o];

    let params = Params {
        algorithm: algo,
        threads: args.threads,
        atom_count: args.atoms,
        atom_counts, // Explicit stoichiometry
        box_size: args.box_size,
        min_distance: 0.85, // Critical collapse distance
        population_size: 24,
        mutation_rate: 0.2,
        crossover_rate: 0.6,
        elitism_count: 2,
        temperature: 300.0,
        step_size: 0.1,
        max_steps: 1000,
        ..Default::default()
    };

    SystemDefinition {
        species: vec![mg, o],
        params,
    }
}

fn check_dependencies() -> Result<()> {
    // We attempt to run `gulp help`. If GULP is not in PATH, this fails.
    match Command::new("gulp").arg("help").output() {
        Ok(_) => Ok(()),
        Err(_) => Err(anyhow!(
            "Dependency Check Failed: 'gulp' executable not found in PATH.\n\
             KLMC requires GULP to perform energy evaluations.\n\
             Please install GULP or add it to your system PATH."
        )),
    }
}

// --- Main ---

fn main() -> Result<(), Box<dyn Error>> {
    // 1. Safety & Parsing
    setup_panic_hook();
    let args = Args::parse();

    // 2. Pre-flight Checks
    if let Err(e) = check_dependencies() {
        eprintln!("{}", e);
        std::process::exit(1);
    }

    // 3. Initialize System
    let system = create_default_system(&args);

    // 4. Initialize Physics Components
    let grid = Arc::new(InteractionGrid::new(&system.species, 0.75));

    // Evaluator: GULP wrapper
    // We define the Buckingham potential string here
    let species_map = system.species.clone();
    let potentials = r#"
buckingham
Mg core O core 1280.1 0.29969 0.0 0.0 10.0
O core O core 22764.0 0.149 27.88 0.0 10.0
spring
Mg 0.0
O 0.0
    "#
    .trim()
    .to_string();

    let evaluator: Arc<dyn Evaluator> =
        Arc::new(GulpEvaluator::new("gulp", &potentials, species_map));

    // 5. Setup TUI & App State
    let mut tui = TuiContext::new().context("Failed to initialize TUI")?;
    let mut app = AppState::new(system.params.clone());

    // 6. Spawn Solver Thread
    let (tx, rx) = unbounded();
    app.set_channel(rx);

    let params_clone = system.params.clone();
    let grid_clone = grid.clone();
    let eval_clone = evaluator.clone();

    thread::Builder::new()
        .name("Solver-Worker".to_string())
        .spawn(move || {
            // Initialize Rayon global thread pool for parallel evaluations
            let _ = rayon::ThreadPoolBuilder::new()
                .num_threads(params_clone.threads)
                .build_global();

            match params_clone.algorithm {
                AlgorithmType::GeneticAlgorithm => {
                    let solver = GeneticAlgorithm::new(eval_clone, grid_clone, params_clone);
                    solver.solve(tx);
                }
                AlgorithmType::BasinHopping => {
                    let mut rng = rand::thread_rng();
                    // Generate a valid starting cluster with correct stoichiometry
                    let start_cluster = Cluster::new_random(
                        &params_clone.atom_counts, // Pass ref to Vec<usize>
                        params_clone.box_size,
                        &grid_clone,
                        &mut rng,
                    )
                    .unwrap_or_else(|| Cluster::new("Fallback_Empty"));

                    let solver = BasinHopping::new(eval_clone, grid_clone, params_clone);
                    solver.solve(start_cluster, tx);
                }
                _ => {
                    // Stub for other algorithms
                }
            }
        })?;

    // 7. Event Loop
    let tick_rate = Duration::from_millis(50); // 20 FPS
    let mut last_tick = Instant::now();

    while !app.should_quit {
        // Draw
        tui.terminal.draw(|f| ui::draw(f, &mut app))?;

        // Handle Input
        let timeout = tick_rate.saturating_sub(last_tick.elapsed());
        if crossterm::event::poll(timeout)? {
            if let Event::Key(key) = event::read()? {
                if key.kind == event::KeyEventKind::Press {
                    match key.code {
                        KeyCode::Char(c) => app.on_key(c),
                        KeyCode::Esc => app.should_quit = true,
                        _ => {}
                    }
                }
            }
        }

        // Logic Tick
        if last_tick.elapsed() >= tick_rate {
            app.tick();
            last_tick = Instant::now();
        }
    }

    Ok(())
}
