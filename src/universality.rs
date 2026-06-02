//! Universality classes: systems that share the same fixed point share critical behavior.
//!
//! The key insight of the RG is that microscopic details don't affect critical behavior.
//! Only the fixed point structure matters: same fixed point → same universality class.

use serde::{Deserialize, Serialize};

use crate::critical_exponents::CriticalExponents;

/// A universality class defined by its fixed point and critical exponents.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniversalityClass {
    /// Name of this universality class.
    pub name: String,
    /// Description of the systems that belong to this class.
    pub description: String,
    /// Spatial dimension.
    pub dimension: usize,
    /// Symmetry group (e.g., Z₂, O(n), etc.).
    pub symmetry: SymmetryGroup,
    /// Range of interaction.
    pub interaction_range: InteractionRange,
    /// Critical exponents.
    pub exponents: CriticalExponents,
}

/// Symmetry group of the order parameter.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SymmetryGroup {
    /// Z₂: Ising-like, scalar order parameter.
    Z2,
    /// O(n): n-component vector order parameter.
    On(usize),
    /// Discrete n-state (Potts model).
    Potts(usize),
    /// U(1): XY model.
    U1,
    /// Custom symmetry.
    Custom(String),
}

impl PartialEq for SymmetryGroup {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (SymmetryGroup::Z2, SymmetryGroup::Z2) => true,
            (SymmetryGroup::On(a), SymmetryGroup::On(b)) => a == b,
            (SymmetryGroup::Potts(a), SymmetryGroup::Potts(b)) => a == b,
            (SymmetryGroup::U1, SymmetryGroup::U1) => true,
            (SymmetryGroup::Custom(a), SymmetryGroup::Custom(b)) => a == b,
            _ => false,
        }
    }
}

/// Range of interaction.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum InteractionRange {
    Short,
    LongRange(f64),
}

/// Well-known universality classes.
impl UniversalityClass {
    /// 2D Ising universality class (Z₂ symmetry).
    pub fn ising_2d() -> Self {
        UniversalityClass {
            name: "2D Ising".to_string(),
            description: "Scalar order parameter with Z₂ symmetry in 2D".to_string(),
            dimension: 2,
            symmetry: SymmetryGroup::Z2,
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents::ising_2d(),
        }
    }

    /// 3D Ising universality class (Z₂ symmetry).
    pub fn ising_3d() -> Self {
        UniversalityClass {
            name: "3D Ising".to_string(),
            description: "Scalar order parameter with Z₂ symmetry in 3D".to_string(),
            dimension: 3,
            symmetry: SymmetryGroup::Z2,
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents::ising_3d(),
        }
    }

    /// 3D XY universality class (U(1) symmetry, O(2)).
    pub fn xy_3d() -> Self {
        UniversalityClass {
            name: "3D XY".to_string(),
            description: "Two-component order parameter with O(2) symmetry in 3D".to_string(),
            dimension: 3,
            symmetry: SymmetryGroup::On(2),
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents {
                nu: 0.671,
                beta: 0.348,
                gamma: 1.317,
                delta: 4.779,
                eta: 0.038,
                alpha: -0.015,
                dimension: 3.0,
            },
        }
    }

    /// 3D Heisenberg universality class (O(3) symmetry).
    pub fn heisenberg_3d() -> Self {
        UniversalityClass {
            name: "3D Heisenberg".to_string(),
            description: "Three-component order parameter with O(3) symmetry in 3D".to_string(),
            dimension: 3,
            symmetry: SymmetryGroup::On(3),
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents {
                nu: 0.711,
                beta: 0.366,
                gamma: 1.389,
                delta: 4.803,
                eta: 0.036,
                alpha: -0.133,
                dimension: 3.0,
            },
        }
    }

    /// Mean-field universality class (above upper critical dimension d≥4).
    pub fn mean_field() -> Self {
        UniversalityClass {
            name: "Mean Field".to_string(),
            description: "Mean-field theory, valid above upper critical dimension d_c=4".to_string(),
            dimension: 4,
            symmetry: SymmetryGroup::Z2,
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents::mean_field(),
        }
    }

    /// 2D Potts q=3 universality class.
    pub fn potts_3_2d() -> Self {
        UniversalityClass {
            name: "2D 3-state Potts".to_string(),
            description: "Three-state Potts model in 2D".to_string(),
            dimension: 2,
            symmetry: SymmetryGroup::Potts(3),
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents {
                nu: 5.0 / 6.0,
                beta: 1.0 / 9.0,
                gamma: 13.0 / 9.0,
                delta: 14.0,
                eta: 4.0 / 15.0,
                alpha: 1.0 / 3.0,
                dimension: 2.0,
            },
        }
    }
}

/// Classify whether two systems belong to the same universality class.
///
/// Systems are in the same universality class if they share:
/// - Spatial dimension
/// - Symmetry of the order parameter
/// - Range of interaction
pub fn same_universality_class(a: &UniversalityClass, b: &UniversalityClass) -> bool {
    a.dimension == b.dimension
        && a.symmetry == b.symmetry
        && a.interaction_range == b.interaction_range
}

/// Compare critical exponents of two universality classes.
///
/// Returns the maximum relative difference across all exponents.
pub fn compare_exponents(a: &UniversalityClass, b: &UniversalityClass) -> f64 {
    let ex_a = &a.exponents;
    let ex_b = &b.exponents;

    let diffs = [
        rel_diff(ex_a.nu, ex_b.nu),
        rel_diff(ex_a.beta, ex_b.beta),
        rel_diff(ex_a.gamma, ex_b.gamma),
        rel_diff(ex_a.delta, ex_b.delta),
        rel_diff(ex_a.eta, ex_b.eta),
        rel_diff(ex_a.alpha, ex_b.alpha),
    ];

    diffs.into_iter().fold(0.0_f64, f64::max)
}

fn rel_diff(a: f64, b: f64) -> f64 {
    if a.abs() < 1e-10 && b.abs() < 1e-10 {
        0.0
    } else if a.abs() < 1e-10 || b.abs() < 1e-10 {
        1.0
    } else {
        ((a - b) / a).abs()
    }
}

/// A registry of known universality classes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniversalityRegistry {
    classes: Vec<UniversalityClass>,
}

impl UniversalityRegistry {
    /// Create a registry with standard universality classes.
    pub fn standard() -> Self {
        UniversalityRegistry {
            classes: vec![
                UniversalityClass::ising_2d(),
                UniversalityClass::ising_3d(),
                UniversalityClass::xy_3d(),
                UniversalityClass::heisenberg_3d(),
                UniversalityClass::mean_field(),
                UniversalityClass::potts_3_2d(),
            ],
        }
    }

    /// Add a universality class to the registry.
    pub fn add(&mut self, class: UniversalityClass) {
        self.classes.push(class);
    }

    /// Find universality classes matching the given criteria.
    pub fn find(&self, dimension: usize, symmetry: &SymmetryGroup) -> Vec<&UniversalityClass> {
        self.classes
            .iter()
            .filter(|c| c.dimension == dimension && c.symmetry == *symmetry)
            .collect()
    }

    /// List all registered classes.
    pub fn list(&self) -> &[UniversalityClass] {
        &self.classes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_ising_2d_exponents() {
        let ising = UniversalityClass::ising_2d();
        assert!((ising.exponents.nu - 1.0).abs() < 1e-10);
        assert!((ising.exponents.beta - 0.125).abs() < 1e-10);
    }

    #[test]
    fn test_ising_3d_exponents() {
        let ising = UniversalityClass::ising_3d();
        assert!((ising.exponents.nu - 0.630).abs() < 0.01);
    }

    #[test]
    fn test_same_class_ising() {
        let a = UniversalityClass::ising_3d();
        let b = UniversalityClass::ising_3d();
        assert!(same_universality_class(&a, &b));
    }

    #[test]
    fn test_different_class() {
        let ising = UniversalityClass::ising_3d();
        let heisenberg = UniversalityClass::heisenberg_3d();
        assert!(!same_universality_class(&ising, &heisenberg));
    }

    #[test]
    fn test_registry_standard() {
        let reg = UniversalityRegistry::standard();
        assert!(reg.list().len() >= 5);
    }

    #[test]
    fn test_registry_find() {
        let reg = UniversalityRegistry::standard();
        let found = reg.find(3, &SymmetryGroup::Z2);
        assert!(found.len() >= 1);
        assert_eq!(found[0].name, "3D Ising");
    }

    #[test]
    fn test_compare_exponents_same() {
        let a = UniversalityClass::ising_3d();
        let b = UniversalityClass::ising_3d();
        assert!(compare_exponents(&a, &b) < 1e-10);
    }

    #[test]
    fn test_compare_exponents_different() {
        let ising = UniversalityClass::ising_3d();
        let heisenberg = UniversalityClass::heisenberg_3d();
        assert!(compare_exponents(&ising, &heisenberg) > 0.01);
    }

    #[test]
    fn test_mean_field_class() {
        let mf = UniversalityClass::mean_field();
        assert_eq!(mf.dimension, 4);
        assert!(mf.exponents.verify_rushbrooke().abs() < 1e-10);
    }

    #[test]
    fn test_add_custom_class() {
        let mut reg = UniversalityRegistry::standard();
        let custom = UniversalityClass {
            name: "Custom".to_string(),
            description: "Test".to_string(),
            dimension: 5,
            symmetry: SymmetryGroup::Z2,
            interaction_range: InteractionRange::Short,
            exponents: CriticalExponents::mean_field(),
        };
        reg.add(custom);
        assert_eq!(reg.list().len(), 7);
    }
}
