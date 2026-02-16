use crate::core::domain::Species;
use serde::{Deserialize, Serialize};

/// A flattened 2D matrix storing pre-computed interaction thresholds.
/// Access is O(1) via `index = i * N + j`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InteractionGrid {
    num_species: usize,
    /// Stores (radius_i + radius_j)^2 * limit_factor
    /// We store squared values to avoid sqrt() calls during simulation.
    collision_matrix_sq: Vec<f64>,
}

impl InteractionGrid {
    /// Builds the grid from the system definition.
    /// `covalent_scale`: Multiplier for covalent collision checks (usually ~0.7-0.8).
    pub fn new(species: &[Species], covalent_scale: f64) -> Self {
        let n = species.len();
        let mut grid = vec![0.0; n * n];

        for i in 0..n {
            for j in 0..n {
                // Determine which radius to use.
                // KLMC logic: if same species -> covalent, else -> ionic (heuristic).
                // Here we default to covalent for collision checks to be strict.
                let r_i = species[i].radius_covalent;
                let r_j = species[j].radius_covalent;

                let dist = r_i + r_j;
                let threshold = dist * covalent_scale;

                // Store squared threshold
                grid[i * n + j] = threshold * threshold;
            }
        }

        Self {
            num_species: n,
            collision_matrix_sq: grid,
        }
    }

    /// Returns the squared distance limit below which two atoms are considered colliding.
    #[inline(always)]
    pub fn get_collision_sq(&self, id_a: usize, id_b: usize) -> f64 {
        // Safety check omitted for speed in release builds; ensure IDs are valid upstream.
        self.collision_matrix_sq[id_a * self.num_species + id_b]
    }
}
