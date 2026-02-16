# Future Improvements & Roadmap

## 1. Configuration & Flexibility
- [ ] **Config File Support**: Implement `config.toml` or `YAML` parsing to define species (mass, charge, radii) and potential parameters externally. Currently hardcoded in `create_default_system`.
- [ ] **Dynamic GULP Templates**: Allow users to provide a template file for GULP input generation, enabling support for arbitrary potentials without recompiling.
- [ ] **Extended CLI**: Add command-line arguments for all simulation parameters (e.g., `population_size`, `mutation_rate`, `temperature`).

## 2. Physics Engine
- [ ] **Generic Evaluator Interface**: Abstract the `Evaluator` trait further to support other engines like LAMMPS, VASP, or DFT codes (CP2K).
- [ ] **Robust Error Recovery**: Enhance `GulpEvaluator` to detect specific convergence failures and retry with different minimization algorithms (e.g., `newton` vs `conjugate gradient`).
- [ ] **Native Force Fields**: Implement a simple Lennard-Jones or Buckingham potential directly in Rust for ultra-fast pre-screening before GULP relaxation.

## 3. Algorithm Enhancements
- [ ] **Multi-Objective Optimization**: Support optimizing for properties other than energy (e.g., band gap, bulk modulus).
- [ ] **Advanced Crossover**: Implement plane-cut or sphere-cut crossover operators for better topological mixing.
- [ ] **Meta-Heuristics**: Add support for Simulated Annealing or Particle Swarm Optimization.

## 4. User Interface & Visualization
- [ ] **3D Visualization**: Integrate a lightweight 3D viewer or export `.xyz` trajectory files in real-time.
- [ ] **History Plotting**: Add a line chart widget to the TUI to visualize the energy drop over generations.
- [ ] **Logging**: Implement structured logging to a file for post-run analysis.

## 5. Code Quality & Testing
- [ ] **Benchmark Suite**: Add benchmarks (`criterion`) to optimize the critical paths (e.g., collision detection in `Cluster::new_random`).
- [ ] **Fuzz Testing**: Use `cargo-fuzz` to test the robustness of the input parsers and crossover logic.
- [ ] **Refactor `main.rs`**: Move the TUI event loop into `src/interface/app.rs` to clean up the entry point.
