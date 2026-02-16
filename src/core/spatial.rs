use nalgebra::{Point3, Vector3};
use crate::core::domain::{Cluster, Lattice};
use crate::core::chemistry::InteractionGrid;

/// Calculates the squared distance between two points.
/// If `lattice` is provided, applies Minimum Image Convention (MIC).
#[inline]
pub fn distance_sq(p1: &Point3<f64>, p2: &Point3<f64>, lattice: Option<&Lattice>) -> f64 {
    match lattice {
        Some(lat) => {
            // Periodic: Convert delta to fractional coordinates
            let d_cart = p2 - p1;
            let mut d_frac = lat.inverse * d_cart;

            // Apply MIC: Wrap fractional coordinates to [-0.5, 0.5]
            d_frac.x -= d_frac.x.round();
            d_frac.y -= d_frac.y.round();
            d_frac.z -= d_frac.z.round();

            // Convert back to Cartesian to get real distance
            let d_mic = lat.vectors * d_frac;
            d_mic.norm_squared()
        }
        None => {
            // Euclidean: Standard distance
            nalgebra::distance_squared(p1, p2)
        }
    }
}

/// Checks a cluster for any physical overlaps (hard collisions).
/// Returns `true` if the cluster is valid (no overlaps).
pub fn check_overlap(cluster: &Cluster, grid: &InteractionGrid) -> bool {
    let atoms = &cluster.atoms;
    let n = atoms.len();
    let lattice = cluster.lattice.as_ref();

    for i in 0..n {
        for j in (i + 1)..n {
            let a_i = &atoms[i];
            let a_j = &atoms[j];

            // Get the squared threshold for this pair
            let threshold_sq = grid.get_collision_sq(a_i.element_id, a_j.element_id);
            
            // Calculate actual squared separation
            let dist_sq = distance_sq(&a_i.position, &a_j.position, lattice);

            if dist_sq < threshold_sq {
                return false; // Collision detected
            }
        }
    }
    true
}

/// Moves a point into the primary unit cell (Periodic only) or centers it (0D).
/// 
/// For 3D (Periodic): Wraps atoms into [0, 1) fractional box.
/// For 0D (Cluster): Centers the geometric center to (0,0,0).
/// 
/// **Invariant**: Modifies positions in-place. Does NOT reorder atoms.
pub fn wrap_or_center(cluster: &mut Cluster) {
    if let Some(lat) = &cluster.lattice {
        // 3D: Wrap atoms into [0, 1) fractional box
        for atom in &mut cluster.atoms {
            let mut frac = lat.to_fractional(&atom.position);
            frac.coords.x = frac.coords.x.rem_euclid(1.0);
            frac.coords.y = frac.coords.y.rem_euclid(1.0);
            frac.coords.z = frac.coords.z.rem_euclid(1.0);
            atom.position = lat.to_cartesian(&frac);
        }
    } else {
        // 0D: Center of Geometry to (0,0,0)
        let n = cluster.atoms.len() as f64;
        
        if n == 0.0 { return; }

        let mut center = Vector3::zeros();
        for atom in &cluster.atoms {
            center += atom.position.coords;
        }
        center /= n;

        for atom in &mut cluster.atoms {
            atom.position -= center;
        }
    }
}