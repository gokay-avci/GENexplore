use std::process::{Command, Stdio};
use std::io::Write;
use anyhow::{anyhow, Context, Result, bail};

use crate::core::domain::{Cluster, Species};
use crate::engine::evaluator::{Evaluator, EvaluationResult};

/// A high-performance, in-memory wrapper for GULP.
/// Streams input/output via pipes to avoid disk latency where possible.
pub struct GulpEvaluator {
    executable: String,
    potential_parameters: String,
    species_map: Vec<Species>,
}

impl GulpEvaluator {
    /// Creates a new evaluator.
    /// 
    /// # Arguments
    /// * `executable` - Path to GULP binary (e.g., "gulp").
    /// * `potential_parameters` - The potential block (buckingham, spring, etc.).
    /// * `species_map` - Ordered list of species corresponding to element_ids in Clusters.
    pub fn new(executable: &str, potential_parameters: &str, species_map: Vec<Species>) -> Self {
        Self {
            executable: executable.to_string(),
            potential_parameters: potential_parameters.to_string(),
            species_map,
        }
    }

    /// Constructs the GULP input string.
    fn generate_input(&self, cluster: &Cluster) -> Result<String> {
        let mut s = String::with_capacity(1024);

        // 1. Header Keywords
        if cluster.lattice.is_some() {
            s.push_str("opti conv conp properties\n");
        } else {
            s.push_str("opti conv cartesian properties\n");
        }

        // 2. Lattice Vectors (if periodic)
        if let Some(lat) = &cluster.lattice {
            s.push_str("vectors\n");
            let v = lat.vectors;
            // GULP reads vectors as rows
            s.push_str(&format!("{:.9} {:.9} {:.9}\n", v[(0,0)], v[(1,0)], v[(2,0)]));
            s.push_str(&format!("{:.9} {:.9} {:.9}\n", v[(0,1)], v[(1,1)], v[(2,1)]));
            s.push_str(&format!("{:.9} {:.9} {:.9}\n", v[(0,2)], v[(1,2)], v[(2,2)]));
        }

        // 3. Coordinates
        if let Some(lat) = &cluster.lattice {
            s.push_str("fractional\n");
            for atom in &cluster.atoms {
                let spec = self.species_map.get(atom.element_id)
                    .ok_or_else(|| anyhow!("Invalid element_id {}", atom.element_id))?;
                
                let frac = lat.to_fractional(&atom.position);
                s.push_str(&format!("{:<3} core {:.9} {:.9} {:.9}\n", spec.symbol, frac.x, frac.y, frac.z));
            }
        } else {
            s.push_str("cartesian\n");
            for atom in &cluster.atoms {
                let spec = self.species_map.get(atom.element_id)
                    .ok_or_else(|| anyhow!("Invalid element_id {}", atom.element_id))?;
                
                let p = atom.position;
                s.push_str(&format!("{:<3} core {:.9} {:.9} {:.9}\n", spec.symbol, p.x, p.y, p.z));
            }
        }

        // 4. Potentials
        s.push('\n');
        s.push_str(&self.potential_parameters);
        s.push('\n');

        Ok(s)
    }

    /// Executes GULP via stdin/stdout piping.
    fn run_process(&self, input_data: &str) -> Result<String> {
        let mut child = Command::new(&self.executable)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .context("Failed to spawn GULP executable")?;

        if let Some(mut stdin) = child.stdin.take() {
            stdin.write_all(input_data.as_bytes())
                .context("Failed to write to GULP stdin")?;
        }

        let output = child.wait_with_output().context("Failed to read GULP output")?;

        if !output.status.success() {
            let err_msg = String::from_utf8_lossy(&output.stderr);
            bail!("GULP exited with error: {}", err_msg);
        }

        let stdout = String::from_utf8_lossy(&output.stdout).to_string();
        Ok(stdout)
    }

    fn parse_energy(&self, output: &str) -> Result<f64> {
        for line in output.lines() {
            let lower = line.to_ascii_lowercase();
            if lower.contains("final energy") || lower.contains("total lattice energy") {
                if let Some(parts) = line.split('=').nth(1) {
                    let tokens: Vec<&str> = parts.split_whitespace().collect();
                    if let Some(val_str) = tokens.first() {
                        let val = val_str.parse::<f64>()
                            .context("Failed to parse energy float")?;
                        return Ok(val);
                    }
                }
            }
        }
        bail!("Could not find final energy in GULP output");
    }

    fn parse_gnorm(&self, output: &str) -> Option<f64> {
        for line in output.lines() {
            if line.to_ascii_lowercase().contains("final gnorm") {
                if let Some(parts) = line.split('=').nth(1) {
                    if let Ok(val) = parts.trim().parse::<f64>() {
                        return Some(val);
                    }
                }
            }
        }
        None
    }

    fn parse_geometry(&self, output: &str, original: &Cluster) -> Result<Cluster> {
        let mut new_cluster = original.clone();
        let lines: Vec<&str> = output.lines().collect();
        let mut start_idx = None;
        let mut is_fractional = false;

        // Find the LAST occurrence of coordinates
        for (i, line) in lines.iter().enumerate().rev() {
            let lower = line.to_ascii_lowercase();
            if lower.contains("final fractional coordinates") {
                start_idx = Some(i + 5); 
                is_fractional = true;
                break;
            } else if lower.contains("final cartesian coordinates") {
                start_idx = Some(i + 5);
                is_fractional = false;
                break;
            }
        }

        let start = start_idx.ok_or_else(|| anyhow!("No final coordinates found in GULP output"))?;
        let expected_atoms = original.atoms.len();
        let mut count = 0;

        for line in lines.into_iter().skip(start) {
            if count >= expected_atoms { break; }
            
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 6 { continue; } 
            if line.contains("-------") { break; } 

            // Skip shells if present (we only update cores)
            if parts[2].to_lowercase().starts_with('s') { continue; }

            // Parse XYZ
            let x: f64 = parts[3].parse().unwrap_or(f64::NAN);
            let y: f64 = parts[4].parse().unwrap_or(f64::NAN);
            let z: f64 = parts[5].parse().unwrap_or(f64::NAN);

            if x.is_nan() || y.is_nan() || z.is_nan() {
                bail!("Parsed NaN coordinates from GULP output");
            }

            if is_fractional {
                if let Some(lat) = &new_cluster.lattice {
                    let frac = nalgebra::Point3::new(x, y, z);
                    new_cluster.atoms[count].position = lat.to_cartesian(&frac);
                } else {
                    bail!("GULP returned fractional coords but cluster has no lattice");
                }
            } else {
                new_cluster.atoms[count].position = nalgebra::Point3::new(x, y, z);
            }
            
            count += 1;
        }

        // STRICT VALIDATION
        if count != expected_atoms {
            bail!("GULP atom count mismatch: expected {}, got {}. Geometry update aborted.", expected_atoms, count);
        }

        Ok(new_cluster)
    }

    fn check_errors(&self, output: &str) -> Result<()> {
        if output.contains("Conditions for a minimum have not been satisfied") {
            bail!("Convergence failure");
        }
        if output.contains("Interatomic distance too small") {
            bail!("Geometric collapse");
        }
        if output.contains("Dump of error info") {
            bail!("Internal GULP error");
        }
        Ok(())
    }
}

impl Evaluator for GulpEvaluator {
    fn name(&self) -> &str { "GULP (Pipe)" }

    fn evaluate(&self, cluster: &Cluster) -> Result<EvaluationResult> {
        let input_str = self.generate_input(cluster)?;
        let output_str = self.run_process(&input_str)?;

        self.check_errors(&output_str)?;

        let energy = self.parse_energy(&output_str)?;
        let gnorm = self.parse_gnorm(&output_str);
        
        // If geometry parsing fails (e.g. mismatch), we propagate the error 
        // so the solver knows this evaluation is invalid/partial.
        let relaxed_cluster = match self.parse_geometry(&output_str, cluster) {
            Ok(c) => Some(c),
            Err(e) => return Err(anyhow!("Geometry parsing failed: {}", e)),
        };

        Ok(EvaluationResult {
            energy,
            gradient_norm: gnorm,
            relaxed_cluster,
        })
    }
}