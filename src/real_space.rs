//! Real-space RG: block spin transformations and Migdal-Kadanoff approximation.
//!
//! Real-space RG operates directly in position space by coarse-graining
//! degrees of freedom (e.g., block spins) and rescaling.

use serde::{Deserialize, Serialize};

use crate::critical_exponents::CriticalExponents;

/// Configuration for a real-space RG transformation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RealSpaceConfig {
    /// Block scaling factor b.
    pub b: f64,
    /// Spatial dimension d.
    pub dimension: usize,
    /// Number of RG iterations.
    pub n_iterations: usize,
}

impl Default for RealSpaceConfig {
    fn default() -> Self {
        RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 20,
        }
    }
}

/// Result of a real-space RG transformation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RealSpaceResult {
    /// Couplings at each RG step.
    pub couplings: Vec<RealSpaceCouplings>,
    /// Estimated critical exponents.
    pub exponents: Option<CriticalExponents>,
    /// Whether the flow converged.
    pub converged: bool,
}

/// Couplings for a lattice model after real-space RG.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RealSpaceCouplings {
    /// Nearest-neighbor coupling K₁ = J/(k_B T).
    pub k1: f64,
    /// Next-nearest-neighbor coupling.
    pub k2: f64,
    /// Four-spin coupling.
    pub k4: f64,
    /// External field h/(k_B T).
    pub h: f64,
}

impl RealSpaceCouplings {
    /// Create with just nearest-neighbor coupling.
    pub fn ising(k1: f64) -> Self {
        RealSpaceCouplings {
            k1,
            k2: 0.0,
            k4: 0.0,
            h: 0.0,
        }
    }
}

/// Majority-rule block spin transformation for 1D Ising model.
///
/// Decimation: sum over every other spin.
/// Result: K' = atanh(tanh(K)²)
pub fn decimation_1d_ising(k: f64) -> f64 {
    let th = k.tanh();
    (th * th).atanh()
}

/// Migdal-Kadanoff bond-moving approximation.
///
/// This is an approximate RG that gives exact results for hierarchical lattices.
/// For d-dimensional Ising model with scale factor b:
///   K' = b^{d-1} * K_eff  (bond moving)
///   then K_eff from solving the 1D decimation chain of length b.
pub fn migdal_kadanoff_step(k: f64, b: f64, d: usize) -> f64 {
    // Bond moving: strengthen by b^{d-1}
    let k_moved = k * b.powi(d as i32 - 1);

    // Decimation chain of length b: apply 1D decimation (b-1) times
    let mut k_eff = k_moved;
    for _ in 1..b as usize {
        k_eff = decimation_1d_ising(k_eff);
    }
    k_eff
}

/// Run the full Migdal-Kadanoff RG flow for the Ising model.
pub fn migdal_kadanoff_flow(config: &RealSpaceConfig, k0: f64) -> RealSpaceResult {
    let mut couplings = Vec::with_capacity(config.n_iterations + 1);
    let mut k = k0;

    couplings.push(RealSpaceCouplings::ising(k));

    let mut converged = false;
    for _ in 0..config.n_iterations {
        k = migdal_kadanoff_step(k, config.b, config.dimension);
        couplings.push(RealSpaceCouplings::ising(k));
    }

    // Check convergence
    if couplings.len() >= 3 {
        let last = couplings.last().unwrap().k1;
        let prev = couplings[couplings.len() - 2].k1;
        if (last - prev).abs() < 1e-10 {
            converged = true;
        }
    }

    RealSpaceResult {
        couplings,
        exponents: None,
        converged,
    }
}

/// Find the critical coupling using Migdal-Kadanoff by bisection.
pub fn find_critical_coupling_mk(config: &RealSpaceConfig, tol: f64) -> f64 {
    let mut k_low = 0.0;
    let mut k_high = 2.0;

    for _ in 0..100 {
        let k_mid = (k_low + k_high) / 2.0;
        let result = migdal_kadanoff_flow(config, k_mid);

        // Check if the coupling flows to 0 (disordered) or infinity (ordered)
        let final_k = result.couplings.last().unwrap().k1;

        if final_k.is_infinite() || final_k > 100.0 || final_k.is_nan() {
            // Ordered phase — k too large
            k_high = k_mid;
        } else if final_k < 0.01 {
            // Disordered phase — k too small
            k_low = k_mid;
        } else {
            // Near critical
            if final_k > k_mid {
                k_low = k_mid;
            } else {
                k_high = k_mid;
            }
        }

        if (k_high - k_low).abs() < tol {
            break;
        }
    }

    (k_low + k_high) / 2.0
}

/// Compute critical exponent ν from the Migdal-Kadanoff RG.
///
/// At the fixed point, ν = ln(b) / ln(∂K'/∂K)|_{K*}.
pub fn nu_from_mk(k_star: f64, b: f64, d: usize) -> f64 {
    let dk = 1e-6;
    let k_plus = k_star + dk;
    let k_minus = k_star - dk;
    let kp_plus = migdal_kadanoff_step(k_plus, b, d);
    let kp_minus = migdal_kadanoff_step(k_minus, b, d);
    let derivative = (kp_plus - kp_minus) / (2.0 * dk);
    let y_t = derivative.ln() / b.ln();
    1.0 / y_t
}

/// Block spin transformation using majority rule.
///
/// For a block of b^d spins, the block spin is +1 if majority are +1.
pub fn majority_rule_block(spins: &[i8]) -> i8 {
    let sum: i32 = spins.iter().map(|&s| s as i32).sum();
    if sum > 0 {
        1
    } else if sum < 0 {
        -1
    } else {
        // Tie: random choice (use +1 as default)
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_decimation_1d_ising_identity() {
        // At K=0, tanh(0)=0, atanh(0)=0
        let k_new = decimation_1d_ising(0.0);
        assert!(k_new.abs() < 1e-10);
    }

    #[test]
    fn test_decimation_1d_ising_reduces() {
        // Decimation always reduces the coupling
        let k_new = decimation_1d_ising(1.0);
        assert!(k_new < 1.0);
    }

    #[test]
    fn test_migdal_kadanoff_step_positive() {
        let k_new = migdal_kadanoff_step(0.5, 2.0, 2);
        assert!(k_new > 0.0);
    }

    #[test]
    fn test_migdal_kadanoff_flow() {
        let config = RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 10,
        };
        let result = migdal_kadanoff_flow(&config, 0.5);
        assert_eq!(result.couplings.len(), 11);
    }

    #[test]
    fn test_migdal_kadanoff_flow_ordered() {
        let config = RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 50,
        };
        let result = migdal_kadanoff_flow(&config, 2.0);
        // Should flow to strong coupling (ordered)
        let final_k = result.couplings.last().unwrap().k1;
        assert!(final_k > 1.0 || final_k.is_infinite());
    }

    #[test]
    fn test_migdal_kadanoff_flow_disordered() {
        let config = RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 50,
        };
        let result = migdal_kadanoff_flow(&config, 0.1);
        // Should flow to zero (disordered)
        let final_k = result.couplings.last().unwrap().k1;
        assert!(final_k < 0.5);
    }

    #[test]
    fn test_find_critical_coupling() {
        let config = RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 50,
        };
        let kc = find_critical_coupling_mk(&config, 0.01);
        // Critical coupling should be positive and finite
        assert!(kc > 0.0);
        assert!(kc < 5.0);
    }

    #[test]
    fn test_majority_rule() {
        assert_eq!(majority_rule_block(&[1, 1, 1, -1, -1]), 1);
        assert_eq!(majority_rule_block(&[-1, -1, -1, 1, 1]), -1);
    }

    #[test]
    fn test_majority_rule_tie() {
        // Tie: defaults to +1
        assert_eq!(majority_rule_block(&[1, -1]), 1);
    }

    #[test]
    fn test_nu_from_mk() {
        let config = RealSpaceConfig {
            b: 2.0,
            dimension: 2,
            n_iterations: 50,
        };
        let kc = find_critical_coupling_mk(&config, 0.001);
        let nu = nu_from_mk(kc, 2.0, 2);
        // ν should be positive and reasonable
        assert!(nu > 0.0);
        assert!(nu < 3.0);
    }
}
