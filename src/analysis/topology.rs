use crate::core::domain::Cluster;
use crate::core::spatial;
use nalgebra::{DMatrix, SymmetricEigen, Matrix3, Vector3, U3};

/// Generates a "Composite Fingerprint" for a cluster.
/// 
/// Definition of a "Gene" (Unique Isomer):
/// 1. Topology (Graph Spectrum): Defines bond connectivity.
/// 2. Geometry (Inertia Tensor): Defines physical shape (Sphere vs Rod vs Disc).
/// 
/// By combining these, we avoid false-positive duplicate detection.
pub fn generate_hash_key(cluster: &Cluster, cutoff_radius: f64) -> String {
    let n = cluster.atoms.len();
    if n == 0 { return "EMPTY".to_string(); }
    if cutoff_radius <= 0.0 { return "INVALID_RADIUS".to_string(); }

    // 1. NaN Check
    for atom in &cluster.atoms {
        if atom.position.coords.iter().any(|c| c.is_nan()) {
            return "NAN_COORDS".to_string();
        }
    }

    // --- Part A: Graph Spectrum (Connectivity) ---
    let mut adjacency = DMatrix::<f64>::zeros(n, n);
    let r_sq = cutoff_radius * cutoff_radius;
    let lattice = cluster.lattice.as_ref();

    for i in 0..n {
        for j in (i + 1)..n {
            let dist_sq = spatial::distance_sq(
                &cluster.atoms[i].position,
                &cluster.atoms[j].position,
                lattice
            );
            if dist_sq < r_sq {
                adjacency[(i, j)] = 1.0;
                adjacency[(j, i)] = 1.0;
            }
        }
    }

    let eigen = SymmetricEigen::new(adjacency);
    let mut evals: Vec<f64> = eigen.eigenvalues.iter().cloned().collect();
    
    // Sort descending for canonical graph representation
    evals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));

    // --- Part B: Principal Moments of Inertia (Shape) ---
    // We assume unit mass for topological comparison to rely purely on geometry
    let pmoi = calculate_pmoi_unit_mass(cluster);

    // --- Part C: Synthesis ---
    // Format: "GS:[e1;e2;...] | PMOI:[i1;i2;i3]"
    // Precision: 3 decimals is usually enough for "Genetic" distinction
    
    let gs_str = evals.iter().map(|e| format!("{:.3}", e)).collect::<Vec<_>>().join(";");
    let pmoi_str = format!("{:.2};{:.2};{:.2}", pmoi[0], pmoi[1], pmoi[2]);

    format!("GS:[{}]|PMOI:[{}]", gs_str, pmoi_str)
}

/// Helper: Calculate PMOI assuming mass=1.0 for all atoms.
/// This provides a purely geometric shape descriptor independent of element types.
fn calculate_pmoi_unit_mass(cluster: &Cluster) -> Vector3<f64> {
    let n = cluster.atoms.len();
    if n < 2 { return Vector3::zeros(); }

    // Center geometry
    let mut center = Vector3::zeros();
    for atom in &cluster.atoms { center += atom.position.coords; }
    center /= n as f64;

    let mut tensor = Matrix3::zeros();
    for atom in &cluster.atoms {
        let r = atom.position.coords - center;
        tensor[(0, 0)] += r.y * r.y + r.z * r.z;
        tensor[(1, 1)] += r.x * r.x + r.z * r.z;
        tensor[(2, 2)] += r.x * r.x + r.y * r.y;
        let i_xy = -r.x * r.y;
        let i_xz = -r.x * r.z;
        let i_yz = -r.y * r.z;
        tensor[(0, 1)] += i_xy; tensor[(1, 0)] += i_xy;
        tensor[(0, 2)] += i_xz; tensor[(2, 0)] += i_xz;
        tensor[(1, 2)] += i_yz; tensor[(2, 1)] += i_yz;
    }

    let eigen: SymmetricEigen<f64, U3> = SymmetricEigen::new(tensor);
    let mut pmoi: Vec<f64> = eigen.eigenvalues.iter().cloned().collect();
    pmoi.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    
    Vector3::new(pmoi[0], pmoi[1], pmoi[2])
}

// Keep the legacy full-mass calculation for other physics needs if present,
// or just use the one above for hashing.
pub fn calculate_pmoi(cluster: &Cluster, atomic_masses: &[f64]) -> Vector3<f64> {
    // ... (Legacy implementation can remain or delegates to similar logic)
    // For this refactor, we just need the hash logic updated.
    // Re-using the logic from previous step for completeness of the file:
    let n = cluster.atoms.len();
    if n < 2 { return Vector3::zeros(); }

    let mut com = Vector3::zeros();
    let mut total_mass = 0.0;
    for atom in &cluster.atoms {
        let m = *atomic_masses.get(atom.element_id).unwrap_or(&1.0);
        com += atom.position.coords * m;
        total_mass += m;
    }
    if total_mass > 0.0 { com /= total_mass; }

    let mut tensor = Matrix3::zeros();
    for atom in &cluster.atoms {
        let m = *atomic_masses.get(atom.element_id).unwrap_or(&1.0);
        let r = atom.position.coords - com;
        tensor[(0, 0)] += m * (r.y * r.y + r.z * r.z);
        tensor[(1, 1)] += m * (r.x * r.x + r.z * r.z);
        tensor[(2, 2)] += m * (r.x * r.x + r.y * r.y);
        let i_xy = -m * r.x * r.y;
        let i_xz = -m * r.x * r.z;
        let i_yz = -m * r.y * r.z;
        tensor[(0, 1)] += i_xy; tensor[(1, 0)] += i_xy;
        tensor[(0, 2)] += i_xz; tensor[(2, 0)] += i_xz;
        tensor[(1, 2)] += i_yz; tensor[(2, 1)] += i_yz;
    }
    let eigen: SymmetricEigen<f64, U3> = SymmetricEigen::new(tensor);
    let mut pmoi: Vec<f64> = eigen.eigenvalues.iter().cloned().collect();
    pmoi.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Vector3::new(pmoi[0], pmoi[1], pmoi[2])
}

pub fn are_duplicates(c1: &Cluster, c2: &Cluster, energy_tol: f64) -> bool {
    match (c1.energy, c2.energy) {
        (Some(e1), Some(e2)) => if (e1 - e2).abs() > energy_tol { return false; },
        _ => return false,
    }
    match (&c1.hash_key, &c2.hash_key) {
        (Some(h1), Some(h2)) => {
            if h1 == "INVALID" || h2 == "INVALID" || h1.contains("NAN") || h2.contains("NAN") { 
                return false; 
            }
            h1 == h2
        },
        _ => false,
    }
}