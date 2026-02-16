use crate::core::domain::{Cluster, Atom};
use crate::core::spatial;
use nalgebra::{Vector3, Rotation3, Unit};
use rand::Rng;
use rand::seq::SliceRandom;

/// A composable mutation builder.
#[derive(Clone, Debug)]
pub struct Mutator {
    rotation_intensity: Option<f64>,    // Max angle
    translation_intensity: Option<f64>, // Max displacement
    rattle_intensity: Option<f64>,      // Max atom displacement
    twist_intensity: Option<f64>,       // Max twist factor
    breathing_intensity: Option<f64>,   // Scaling factor range
    swap_count: Option<usize>,          // Pairs to swap
}

impl Mutator {
    pub fn new() -> Self {
        Self {
            rotation_intensity: None,
            translation_intensity: None,
            rattle_intensity: None,
            twist_intensity: None,
            breathing_intensity: None,
            swap_count: None,
        }
    }

    pub fn rotate(mut self, max_angle: f64) -> Self {
        self.rotation_intensity = Some(max_angle);
        self
    }

    pub fn translate(mut self, max_dist: f64) -> Self {
        self.translation_intensity = Some(max_dist);
        self
    }

    pub fn rattle(mut self, max_dist: f64) -> Self {
        self.rattle_intensity = Some(max_dist);
        self
    }

    pub fn twist(mut self, factor: f64) -> Self {
        self.twist_intensity = Some(factor);
        self
    }

    pub fn breathing(mut self, intensity: f64) -> Self {
        self.breathing_intensity = Some(intensity);
        self
    }

    pub fn swap(mut self, count: usize) -> Self {
        self.swap_count = Some(count);
        self
    }

    pub fn apply(&self, cluster: &Cluster, rng: &mut impl Rng) -> Cluster {
        let mut c = cluster.clone();
        c.origin = "Mutation".to_string(); 
        
        spatial::wrap_or_center(&mut c);

        // 1. Breathing (Global Scaling)
        if let Some(mag) = self.breathing_intensity {
            let scale = 1.0 + rng.gen_range(-mag..mag);
            for atom in &mut c.atoms {
                atom.position.coords *= scale;
            }
        }

        // 2. Rotation
        if let Some(mag) = self.rotation_intensity {
            let axis = Unit::new_normalize(Vector3::new(
                rng.gen::<f64>() - 0.5,
                rng.gen::<f64>() - 0.5,
                rng.gen::<f64>() - 0.5,
            ));
            let angle = rng.gen_range(-mag..mag);
            let rot = Rotation3::from_axis_angle(&axis, angle);

            for atom in &mut c.atoms {
                atom.position = rot * atom.position;
            }
        }

        // 3. Twist
        if let Some(mag) = self.twist_intensity {
            let axis_z = 2.0; 
            for atom in &mut c.atoms {
                let theta = atom.position.z * axis_z * mag * (rng.gen::<f64>() - 0.5);
                let (sin, cos) = theta.sin_cos();
                let x = atom.position.x;
                let y = atom.position.y;
                atom.position.x = x * cos - y * sin;
                atom.position.y = x * sin + y * cos;
            }
        }

        // 4. Rattle
        if let Some(mag) = self.rattle_intensity {
            for atom in &mut c.atoms {
                atom.position.x += rng.gen_range(-mag..mag);
                atom.position.y += rng.gen_range(-mag..mag);
                atom.position.z += rng.gen_range(-mag..mag);
            }
        }

        // 5. Swap
        if let Some(count) = self.swap_count {
            let n = c.atoms.len();
            if n >= 2 {
                for _ in 0..count {
                    let i = rng.gen_range(0..n);
                    let j = rng.gen_range(0..n);
                    if i != j {
                        let tmp = c.atoms[i].position;
                        c.atoms[i].position = c.atoms[j].position;
                        c.atoms[j].position = tmp;
                    }
                }
            }
        }

        // 6. Translation
        if let Some(mag) = self.translation_intensity {
            let dx = rng.gen_range(-mag..mag);
            let dy = rng.gen_range(-mag..mag);
            let dz = rng.gen_range(-mag..mag);
            let disp = Vector3::new(dx, dy, dz);

            for atom in &mut c.atoms {
                atom.position += disp;
            }
        }

        spatial::wrap_or_center(&mut c);
        c
    }
}

// --- Helper for Crossover ---

/// Rotates a set of atoms randomly around their geometric center.
/// Defined as a generic function to avoid `dyn Rng` object safety issues.
fn apply_random_rotation<R: Rng + ?Sized>(atoms: &mut [Atom], rng: &mut R) {
    if atoms.is_empty() { return; }

    // Center manually
    let mut com = Vector3::zeros();
    for a in atoms.iter() { com += a.position.coords; }
    com /= atoms.len() as f64;
    for a in atoms.iter_mut() { a.position -= com; }

    // Rotate
    let axis = Unit::new_normalize(Vector3::new(
        rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5, rng.gen::<f64>() - 0.5
    ));
    let angle = rng.gen_range(0.0..6.28);
    let rot = Rotation3::from_axis_angle(&axis, angle);
    
    for a in atoms.iter_mut() { a.position = rot * a.position; }
}

/// Robust "Cut and Splice" Crossover.
pub fn crossover_cut_splice(p1: &Cluster, p2: &Cluster, rng: &mut impl Rng) -> Option<Cluster> {
    if p1.atoms.len() != p2.atoms.len() { return None; }
    let n = p1.atoms.len();
    if n < 2 { return None; }

    // Target Stoichiometry (Source of Truth)
    let max_id = p1.atoms.iter().map(|a| a.element_id).max().unwrap_or(0);
    let mut target_counts = vec![0; max_id + 1];
    for a in &p1.atoms { target_counts[a.element_id] += 1; }

    let mut child = p1.clone();
    child.origin = format!("X({},{})", p1.id.to_string()[0..4].to_string(), p2.id.to_string()[0..4].to_string());
    
    // 1. Prepare Parents (Clone -> Center -> Rotate)
    let mut p1_atoms = p1.atoms.clone();
    let mut p2_atoms = p2.atoms.clone();

    apply_random_rotation(&mut p1_atoms, rng);
    apply_random_rotation(&mut p2_atoms, rng);
    
    // 2. Sort by Z
    p1_atoms.sort_by(|a, b| a.position.z.total_cmp(&b.position.z));
    p2_atoms.sort_by(|a, b| a.position.z.total_cmp(&b.position.z));

    // 3. Cut
    // range 1..n ensures at least 1 atom from P1 and 1 from P2
    let cut_point = rng.gen_range(1..n);
    
    child.atoms.clear();
    child.atoms.extend_from_slice(&p1_atoms[0..cut_point]);
    child.atoms.extend_from_slice(&p2_atoms[cut_point..n]);
    
    // 4. Repair Stoichiometry (Alchemy)
    let mut child_counts = vec![0; max_id + 1];
    for a in &child.atoms { 
        if a.element_id < child_counts.len() { child_counts[a.element_id] += 1; } 
    }

    let mut deficits = Vec::new();
    for (id, &tgt) in target_counts.iter().enumerate() {
        let curr = child_counts.get(id).copied().unwrap_or(0);
        if curr < tgt {
            for _ in 0..(tgt - curr) { deficits.push(id); }
        }
    }
    deficits.shuffle(rng);

    for (id, &tgt) in target_counts.iter().enumerate() {
        let curr = child_counts.get(id).copied().unwrap_or(0);
        if curr > tgt {
            let surplus = curr - tgt;
            let mut indices: Vec<usize> = child.atoms.iter().enumerate()
                .filter(|(_, a)| a.element_id == id)
                .map(|(i, _)| i)
                .collect();
            indices.shuffle(rng);
            
            for i in 0..surplus {
                if let Some(target_idx) = indices.get(i) {
                    if let Some(new_id) = deficits.pop() {
                        child.atoms[*target_idx].element_id = new_id;
                    }
                }
            }
        }
    }

    spatial::wrap_or_center(&mut child);
    Some(child)
}