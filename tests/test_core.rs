use klmc_ultimate::core::chemistry::InteractionGrid;
use klmc_ultimate::core::domain::{Cluster, Species};
use rand::thread_rng;

#[test]
fn test_cluster_creation() {
    let species = vec![
        Species {
            symbol: "A".into(),
            radius_covalent: 1.0,
            ..Default::default()
        },
        Species {
            symbol: "B".into(),
            radius_covalent: 1.0,
            ..Default::default()
        },
    ];
    // Use a small covalent_scale to make packing easy
    let grid = InteractionGrid::new(&species, 0.5);
    let atom_counts = vec![5, 5];

    let mut rng = thread_rng();
    let cluster = Cluster::new_random(&atom_counts, 10.0, &grid, &mut rng);

    assert!(cluster.is_some(), "Cluster creation failed");
    let c = cluster.unwrap();
    assert_eq!(c.atoms.len(), 10);
    assert!(c.check_stoichiometry(&atom_counts));
}

#[test]
fn test_interaction_grid() {
    let species = vec![
        Species {
            radius_covalent: 2.0,
            ..Default::default()
        }, // 0
        Species {
            radius_covalent: 1.0,
            ..Default::default()
        }, // 1
    ];

    let grid = InteractionGrid::new(&species, 1.0);

    // Test 0-0: (2+2)*1 = 4.0 -> sq = 16.0
    assert!((grid.get_collision_sq(0, 0) - 16.0).abs() < 1e-6);

    // Test 0-1: (2+1)*1 = 3.0 -> sq = 9.0
    assert!((grid.get_collision_sq(0, 1) - 9.0).abs() < 1e-6);

    // Test 1-1: (1+1)*1 = 2.0 -> sq = 4.0
    assert!((grid.get_collision_sq(1, 1) - 4.0).abs() < 1e-6);
}
