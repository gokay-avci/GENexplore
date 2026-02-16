# KLMC Ultimate: Knowledge-Led Master Code (Rust Port)

**KLMC Ultimate** is a high-performance, Rust-based global optimization framework designed for discovering the lowest-energy structures of atomic clusters and nano-materials. It is a modern port and evolution of the original KLMC software, leveraging Rust's safety, concurrency, and performance.

## 🚀 Key Features

*   **Advanced Algorithms**:
    *   **Genetic Algorithm (GA)**: Evolutionary strategy with tournament selection, elitism, and adaptive mutation rates to explore the potential energy surface efficiently. Includes specific operators like "Cut & Splice" crossover and rotational mutation.
    *   **Basin Hopping (BH)**: A Monte Carlo minimization technique that transforms the energy landscape into a set of basins, effectively finding global minima by hopping between local minima.
*   **Physics Engine Integration**:
    *   Seamlessly integrates with **GULP** (General Utility Lattice Program) for accurate interatomic potential evaluations.
    *   Supports Buckingham, Spring, and other potential models via GULP input generation.
*   **High Performance**:
    *   **Parallel Evaluation**: Utilizes `rayon` for multi-threaded energy calculations, scaling with your CPU cores.
    *   **Efficient Architecture**: Minimizes overhead with a dedicated solver thread and non-blocking TUI updates.
*   **Interactive TUI**:
    *   Built with `ratatui` for a rich terminal user interface.
    *   Real-time visualization of simulation statistics (Best Energy, Diversity, Mutation Rate).
    *   Live control over the simulation.

## 🛠️ Architecture

The project is structured as a library + binary workspace:

*   **`klmc_ultimate` (Library)**: Contains the core logic.
    *   `core`: Domain models (`Cluster`, `Species`, `Atom`), spatial utilities, and chemistry definitions (`InteractionGrid`).
    *   `engine`: Interfaces for physics evaluators (`Evaluator` trait) and mutation/crossover operators.
    *   `solvers`: Implementation of optimization algorithms (`GeneticAlgorithm`, `BasinHopping`).
    *   `interface`: State management for the UI.
    *   `analysis`: Topological analysis and duplicate detection.
*   **`src/main.rs` (Binary)**: The entry point that sets up the CLI, TUI, and spawns the solver thread.

## 📦 Installation & Prerequisites

### Prerequisites
1.  **Rust Toolchain**: Install via [rustup.rs](https://rustup.rs).
2.  **GULP**: The General Utility Lattice Program must be installed and accessible in your system `PATH` as `gulp`.
    *   KLMC generates input files for GULP and parses its standard output.

### Building
```bash
cargo build --release
```

### Running Tests
```bash
cargo test
```

## 🖥️ Usage

Run the program using `cargo run` or the built binary.

```bash
# Run Genetic Algorithm (Default) for MgO cluster with 12 atoms
cargo run --release -- -a ga -n 12

# Run Basin Hopping
cargo run --release -- -a bh -n 12 --threads 8
```

### CLI Options
*   `-a, --algo <ALGO>`: Algorithm to use (`ga`, `bh`). Default: `ga`.
*   `-n, --atoms <N>`: Total number of atoms. Default: `12`.
*   `-t, --threads <N>`: Number of worker threads. Default: `4`.
*   `-b, --box-size <SIZE>`: Initial simulation box size in Ångströms. Default: `6.0`.

## 🧠 How It Works

1.  **Initialization**: Random clusters are generated respecting stoichiometry constraints (e.g., Mg6O6) and checking for atomic overlaps using an `InteractionGrid`.
2.  **Evaluation**: Clusters are sent to the `GulpEvaluator`, which pipes geometry to `gulp`, runs an optimization (relaxation), and parses the final energy and coordinates.
3.  **Evolution (GA)**:
    *   **Selection**: Tournament selection picks parents.
    *   **Crossover**: "Cut and Splice" combines halves of two clusters.
    *   **Mutation**: Rotations, Rattling, Twisting, and Breathing modes perturb structures to escape local minima.
    *   **Reseeding**: If diversity drops or stagnation occurs, the population is partially reseeded (Mass Extinction).
4.  **Exploration (BH)**:
    *   A single walker explores the landscape.
    *   Perturbation -> Local Minimization -> Metropolis Acceptance Criterion.

## 🧪 Testing

The project includes a comprehensive testing suite:
*   **Unit Tests**: Verify core domain logic, spatial operations, and stoichiometry.
*   **Integration Tests**: validating the flow of GA and BH solvers using a `MockEvaluator` (simulating physics without requiring GULP installed).

## 📜 License

This project is open-source.
