use nalgebra::{Point3, Vector3, Matrix3};
use serde::{Deserialize, Serialize};
use uuid::Uuid;
use rand::Rng;
use rand::seq::SliceRandom; // Required for shuffling species

// --- Constants ---
pub const MAX_HISTORY: usize = 50;

// --- Physics Types ---

/// Represents a single chemical element/species properties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Species {
    pub symbol: String,
    pub atomic_number: u8,
    pub mass: f64,             // amu
    pub charge: f64,           // e
    pub radius_covalent: f64,  // Å
    pub radius_ionic: f64,     // Å
    pub color_rgb: (u8, u8, u8), // For TUI visualization
}

impl Default for Species {
    fn default() -> Self {
        Self {
            symbol: "X".to_string(),
            atomic_number: 0,
            mass: 1.0,
            charge: 0.0,
            radius_covalent: 1.0,
            radius_ionic: 1.0,
            color_rgb: (255, 255, 255),
        }
    }
}

/// A single atom instance in a cluster.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element_id: usize, // Index into the SystemConfig.species list
    pub position: Point3<f64>,
    pub velocity: Vector3<f64>,
    pub force: Vector3<f64>,
    pub is_fixed: bool,
}

/// Defines the Periodic Boundary Conditions (if any).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Lattice {
    pub vectors: Matrix3<f64>, // Columns are a, b, c
    pub inverse: Matrix3<f64>, // Precomputed for fractional conversion
}

impl Lattice {
    pub fn new(a: Vector3<f64>, b: Vector3<f64>, c: Vector3<f64>) -> Option<Self> {
        let vectors = Matrix3::from_columns(&[a, b, c]);
        let inverse = vectors.try_inverse()?;
        Some(Self { vectors, inverse })
    }
    
    pub fn to_fractional(&self, p: &Point3<f64>) -> Point3<f64> {
        let v = self.inverse * p.coords;
        Point3::from(v)
    }

    pub fn to_cartesian(&self, p: &Point3<f64>) -> Point3<f64> {
        let v = self.vectors * p.coords;
        Point3::from(v)
    }
}

// --- The Core Entity ---

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ClusterStatus {
    Born,        // Just created
    Valid,       // Checks passed
    Evaluated,   // Physics engine returned energy
    Discarded,   // Failed
    Elite,       // Hall of Fame
}

/// The primary data unit passed between Solver and TUI.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Cluster {
    pub id: Uuid,
    pub generation: u64,
    pub origin: String,
    
    pub atoms: Vec<Atom>,
    pub lattice: Option<Lattice>,
    
    pub energy: Option<f64>,
    pub gradient_norm: Option<f64>,
    pub pmoi: Option<Vector3<f64>>,
    pub hash_key: Option<String>,
    
    pub status: ClusterStatus,
}

impl Cluster {
    pub fn new(origin: &str) -> Self {
        Self {
            id: Uuid::new_v4(),
            generation: 0,
            origin: origin.to_string(),
            atoms: Vec::new(),
            lattice: None,
            energy: None,
            gradient_norm: None,
            pmoi: None,
            hash_key: None,
            status: ClusterStatus::Born,
        }
    }

    /// Tries to generate a random cluster respecting stoichiometry constraints.
    /// 
    /// # Arguments
    /// * `atom_counts`: A slice where index `i` is the count of species `i`.
    ///   Example: `[6, 6]` for 6 Mg and 6 O.
    pub fn new_random<R: Rng + ?Sized>(
        atom_counts: &[usize],
        box_size: f64,
        grid: &crate::core::chemistry::InteractionGrid,
        rng: &mut R,
    ) -> Option<Self> {
        let mut c = Cluster::new("Random");
        
        // 1. Build the exact multiset of element IDs required.
        let mut elements_to_place = Vec::new();
        for (id, &count) in atom_counts.iter().enumerate() {
            for _ in 0..count {
                elements_to_place.push(id);
            }
        }
        
        // 2. Shuffle to randomize initial topology.
        elements_to_place.shuffle(rng);

        // 3. Place atoms (Random Sequential Adsorption)
        for &elem_id in &elements_to_place {
            let mut placed = false;
            
            // Attempt 100 times to place an atom without overlap
            for _ in 0..100 {
                let pos = Point3::new(
                    rng.gen_range(-box_size..box_size),
                    rng.gen_range(-box_size..box_size),
                    rng.gen_range(-box_size..box_size),
                );
                
                // Check overlap with already placed atoms
                let mut clash = false;
                for existing in &c.atoms {
                    let limit_sq = grid.get_collision_sq(elem_id, existing.element_id);
                    // Simple euclidean check for generation (0D logic)
                    // TODO: If 3D PBC generation is needed, wrap logic goes here.
                    let dist_sq = (pos - existing.position).norm_squared();
                    
                    if dist_sq < limit_sq {
                        clash = true;
                        break;
                    }
                }
                
                if !clash {
                    c.atoms.push(Atom {
                        element_id: elem_id,
                        position: pos,
                        velocity: Vector3::zeros(),
                        force: Vector3::zeros(),
                        is_fixed: false,
                    });
                    placed = true;
                    break;
                }
            }
            if !placed { return None; } // Failed to pack
        }
        
        // Center the cluster
        if !c.atoms.is_empty() {
            let mut com = Vector3::zeros();
            for a in &c.atoms { com += a.position.coords; }
            com /= c.atoms.len() as f64;
            for a in &mut c.atoms { a.position -= com; }
        }

        Some(c)
    }

    /// Verifies if the cluster matches the target stoichiometry.
    /// Returns true if counts match exactly.
    pub fn check_stoichiometry(&self, target_counts: &[usize]) -> bool {
        let mut actual_counts = vec![0; target_counts.len()];
        
        for atom in &self.atoms {
            if atom.element_id >= actual_counts.len() {
                return false; // Atom has invalid element_id
            }
            actual_counts[atom.element_id] += 1;
        }
        
        actual_counts == target_counts
    }
}

// --- Configuration Types ---

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub enum AlgorithmType {
    GeneticAlgorithm,
    BasinHopping,
    ScanBox,
    SolidSolution,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    pub algorithm: AlgorithmType,
    pub seed: u64,
    pub threads: usize,
    
    // Physics Constraints
    pub atom_count: usize, // Total atoms
    pub atom_counts: Vec<usize>, // Explicit counts per species (e.g., [6, 6])
    pub box_size: f64,
    pub min_distance: f64,
    
    // GA Specific
    pub population_size: usize,
    pub mutation_rate: f64,
    pub crossover_rate: f64,
    pub elitism_count: usize,
    
    // BH Specific
    pub temperature: f64,
    pub step_size: f64,
    pub bh_steps: usize,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            algorithm: AlgorithmType::GeneticAlgorithm,
            seed: 0,
            threads: 4,
            atom_count: 10,
            atom_counts: vec![], // Must be populated by main
            box_size: 10.0,
            min_distance: 0.8,
            population_size: 32,
            mutation_rate: 0.1,
            crossover_rate: 0.6,
            elitism_count: 2,
            temperature: 300.0,
            step_size: 0.1,
            bh_steps: 100,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SystemDefinition {
    pub species: Vec<Species>,
    pub params: Params,
}