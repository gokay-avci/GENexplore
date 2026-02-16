use klmc_ultimate::core::domain::{Cluster, Atom};
use klmc_ultimate::engine::operators::{Mutator, crossover_cut_splice};
use nalgebra::{Point3, Vector3};
use rand::thread_rng;

fn create_dummy_cluster(n: usize) -> Cluster {
    let mut c = Cluster::new("Test");
    for i in 0..n {
        c.atoms.push(Atom {
            element_id: i % 2,
            position: Point3::new(i as f64, 0.0, 0.0),
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            is_fixed: false,
        });
    }
    c
}

#[test]
fn test_mutation_rattle() {
    let c = create_dummy_cluster(2);
    let mut rng = thread_rng();

    // Rattle should change positions relative to each other
    let mutated = Mutator::new().rattle(0.5).apply(&c, &mut rng);

    let dist_orig = (c.atoms[0].position - c.atoms[1].position).norm();
    let dist_new = (mutated.atoms[0].position - mutated.atoms[1].position).norm();

    assert!((dist_orig - dist_new).abs() > 1e-9);
}

#[test]
fn test_crossover() {
    let p1 = create_dummy_cluster(4); // 2 of type 0, 2 of type 1
    let p2 = create_dummy_cluster(4);

    let mut rng = thread_rng();
    let child = crossover_cut_splice(&p1, &p2, &mut rng);

    assert!(child.is_some());
    let c = child.unwrap();
    assert_eq!(c.atoms.len(), 4);

    // Check Stoichiometry (should be preserved)
    let count0 = c.atoms.iter().filter(|a| a.element_id == 0).count();
    let count1 = c.atoms.iter().filter(|a| a.element_id == 1).count();
    assert_eq!(count0, 2);
    assert_eq!(count1, 2);
}
