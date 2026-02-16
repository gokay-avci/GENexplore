use crate::core::domain::{Cluster, Species};
use crate::engine::evaluator::{Evaluator, EvaluationResult};
use anyhow::{anyhow, Result};
use std::process::Command;
use std::fs;
use uuid::Uuid;

pub struct GulpEvaluator {
    /// Path to the GULP executable.
    executable: String,
    /// The GULP library/species definitions (e.g., buckingham potentials).
    /// This string is appended to the input file.
    potential_parameters: String,
}

impl GulpEvaluator {
    pub fn new(executable: &str, potential_parameters: &str) -> Self {
        Self {
            executable: executable.to_string(),
            potential_parameters: potential_parameters.to_string(),
        }
    }

    /// Generates GULP input content.
    fn generate_input(&self, cluster: &Cluster, species_map: &[Species]) -> String {
        let mut s = String::new();
        
        // Keywords: optimize, constant pressure (if crystal), cartesian coords
        if cluster.lattice.is_some() {
            s.push_str("opti conv conp properties\n");
        } else {
            s.push_str("opti conv cartesian properties\n");
        }

        // Lattice Vectors (if periodic)
        if let Some(lat) = &cluster.lattice {
            s.push_str("vectors\n");
            // GULP expects vectors as rows
            let v = lat.vectors; 
            // col 0 is a, col 1 is b...
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", v[(0,0)], v[(1,0)], v[(2,0)]));
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", v[(0,1)], v[(1,1)], v[(2,1)]));
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", v[(0,2)], v[(1,2)], v[(2,2)]));
        }

        // Coordinates
        if cluster.lattice.is_some() {
             s.push_str("fractional\n");
             if let Some(lat) = &cluster.lattice {
                 for atom in &cluster.atoms {
                     let spec = &species_map[atom.element_id];
                     let frac = lat.to_fractional(&atom.position);
                     s.push_str(&format!("{:<3} core {:.6} {:.6} {:.6}\n", spec.symbol, frac.x, frac.y, frac.z));
                 }
             }
        } else {
            s.push_str("cartesian\n");
            for atom in &cluster.atoms {
                let spec = &species_map[atom.element_id];
                let p = atom.position;
                s.push_str(&format!("{:<3} core {:.6} {:.6} {:.6}\n", spec.symbol, p.x, p.y, p.z));
            }
        }

        // Potentials
        s.push_str("\n");
        s.push_str(&self.potential_parameters);
        
        s
    }

    /// Parses the .gout file to extract energy and new coords.
    fn parse_output(&self, output: &str, original_cluster: &Cluster) -> Result<EvaluationResult> {
        // 1. Check for failure flags
        if output.contains("Conditions for a minimum have not been satisfied") {
             return Err(anyhow!("GULP: Convergence failed"));
        }
        if output.contains("Interatomic distance too small") {
             return Err(anyhow!("GULP: Geometry collapsed"));
        }
        if output.contains("Dump of error info") {
             return Err(anyhow!("GULP: Internal error"));
        }

        // 2. Extract Energy
        let mut energy = None;
        let mut gnorm = None;

        for line in output.lines() {
            if line.contains("Final energy") || line.contains("Total lattice energy") {
                 // Format: "  Final energy =  -123.456 eV"
                 if let Some(val_str) = line.split('=').nth(1) {
                     let clean = val_str.split_whitespace().next().unwrap_or("");
                     if let Ok(v) = clean.parse::<f64>() {
                         energy = Some(v);
                     }
                 }
            }
            if line.contains("Final Gnorm") {
                 if let Some(val_str) = line.split('=').nth(1) {
                     let clean = val_str.split_whitespace().next().unwrap_or("");
                     if let Ok(v) = clean.parse::<f64>() {
                         gnorm = Some(v);
                     }
                 }
            }
        }

        let e = energy.ok_or_else(|| anyhow!("GULP: Could not parse final energy"))?;

        // 3. Extract Geometry (Simplified for this snippet)
        // TODO: Implement full coordinate parsing logic here to update `relaxed_cluster`
        let relaxed_cluster = Some(original_cluster.clone()); 

        Ok(EvaluationResult {
            energy: e,
            gradient_norm: gnorm,
            relaxed_cluster,
        })
    }
}

impl Evaluator for GulpEvaluator {
    fn name(&self) -> &str { "GULP" }

    fn evaluate(&self, cluster: &Cluster) -> Result<EvaluationResult> {
        // Sandboxing
        let job_id = Uuid::new_v4();
        let sandbox_dir = std::env::temp_dir().join(format!("klmc_run_{}", job_id));
        fs::create_dir_all(&sandbox_dir)?;

        let inp_path = sandbox_dir.join("job.gin");

        // TEMPORARY: Panic if we can't map symbols. 
        // Real impl: Store `Vec<Species>` in GulpEvaluator or pass via method.
        // For now, we assume standard MgO if index 0/1 to prevent panics in demo.
        let mg = Species { symbol: "Mg".into(), ..Default::default() };
        let o = Species { symbol: "O".into(), ..Default::default() };
        let species_dummy = vec![mg, o];
        
        let input_content = self.generate_input(cluster, &species_dummy); 

        fs::write(&inp_path, input_content)?;

        // Execute
        let output = Command::new(&self.executable)
            .arg(inp_path.file_name().unwrap()) // GULP usually takes arg or redirection
            .current_dir(&sandbox_dir) // Run inside dir
            .output();

        // Cleanup (always try, ignore errors)
        let _ = fs::remove_dir_all(&sandbox_dir); 

        match output {
            Ok(o) => {
                if !o.status.success() {
                    let err = String::from_utf8_lossy(&o.stderr);
                    return Err(anyhow!("GULP execution failed: {}", err));
                }
                
                let stdout = String::from_utf8_lossy(&o.stdout);
                self.parse_output(&stdout, cluster)
            },
            Err(e) => Err(anyhow!("Failed to spawn GULP: {}", e)),
        }
    }
}