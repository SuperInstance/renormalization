//! Epsilon expansion: systematic d = 4-ε analysis.
//!
//! The epsilon expansion is a perturbative technique where one computes critical
//! exponents as power series in ε = 4 - d, valid near the upper critical dimension.

use serde::{Deserialize, Serialize};

use crate::critical_exponents::CriticalExponents;

/// Epsilon expansion results for φ⁴ theory with O(n) symmetry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EpsilonExpansion {
    /// Number of components n in O(n) model.
    pub n: usize,
    /// Epsilon = 4 - d.
    pub epsilon: f64,
    /// Fixed point coupling g* to order ε².
    pub g_star: Vec<f64>, // [order 0, order 1, order 2, ...]
    /// η to various orders in ε.
    pub eta: Vec<f64>,
    /// ν to various orders in ε.
    pub nu: Vec<f64>,
}

impl EpsilonExpansion {
    /// Compute the Wilson-Fisher fixed point coupling to O(ε).
    ///
    /// For O(n) φ⁴ theory:
    ///   β(g) = -εg + (n+8)/(48π²) g² - (3n+14)/(384π⁴) g³ + ...
    ///   g* = 48π²ε/(n+8) + O(ε²)
    pub fn wilson_fisher_fixed_point_1loop(n: usize, epsilon: f64) -> f64 {
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        48.0 * pi2 * epsilon / (n as f64 + 8.0)
    }

    /// Compute the Wilson-Fisher fixed point coupling to O(ε²).
    pub fn wilson_fisher_fixed_point_2loop(n: usize, epsilon: f64) -> f64 {
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        let g1 = 48.0 * pi2 * epsilon / (n as f64 + 8.0);
        let g2 = g1 * g1 * (3.0 * n as f64 + 14.0) / (384.0 * pi2 * pi2);
        g1 + g2 / (48.0 * pi2 / (n as f64 + 8.0)) * epsilon // approximate
    }

    /// Compute η to O(ε²) for O(n) model.
    ///
    /// η = (n+2)ε² / (2(n+8)²) + O(ε³)
    pub fn eta_o_eps2(n: usize, epsilon: f64) -> f64 {
        let nn = n as f64;
        (nn + 2.0) * epsilon * epsilon / (2.0 * (nn + 8.0) * (nn + 8.0))
    }

    /// Compute ν to O(ε) for O(n) model.
    ///
    /// ν = 1/2 + (n+2)ε / (4(n+8)) + O(ε²)
    pub fn nu_o_eps1(n: usize, epsilon: f64) -> f64 {
        let nn = n as f64;
        0.5 + (nn + 2.0) * epsilon / (4.0 * (nn + 8.0))
    }

    /// Compute ν to O(ε²) for O(n) model.
    pub fn nu_o_eps2(n: usize, epsilon: f64) -> f64 {
        let nn = n as f64;
        let nu1 = 0.5 + (nn + 2.0) * epsilon / (4.0 * (nn + 8.0));
        let nu2 = (nn + 2.0) * (nn + 2.0 + nn + 14.0) * epsilon * epsilon
            / (8.0 * (nn + 8.0) * (nn + 8.0) * 4.0);
        nu1 + nu2
    }

    /// Full epsilon expansion for O(n) model at given ε.
    pub fn compute(n: usize, epsilon: f64) -> Self {
        EpsilonExpansion {
            n,
            epsilon,
            g_star: vec![
                0.0,
                Self::wilson_fisher_fixed_point_1loop(n, epsilon),
            ],
            eta: vec![0.0, Self::eta_o_eps2(n, epsilon)],
            nu: vec![0.5, Self::nu_o_eps1(n, epsilon), Self::nu_o_eps2(n, epsilon)],
        }
    }

    /// Get critical exponents from the epsilon expansion.
    pub fn critical_exponents(&self) -> CriticalExponents {
        let nu = *self.nu.last().unwrap_or(&0.5);
        let eta = *self.eta.last().unwrap_or(&0.0);
        let d = 4.0 - self.epsilon;

        crate::critical_exponents::exponents_from_nu_eta(nu, eta, d)
    }
}

/// Run the epsilon expansion for a range of ε values and extract the fixed point.
pub fn epsilon_scan(n: usize, eps_values: &[f64]) -> Vec<(f64, f64)> {
    eps_values
        .iter()
        .map(|&eps| {
            let g_star = EpsilonExpansion::wilson_fisher_fixed_point_1loop(n, eps);
            (eps, g_star)
        })
        .collect()
}

/// Check convergence of the epsilon expansion by comparing successive orders.
pub fn convergence_check(expansion: &EpsilonExpansion) -> f64 {
    if expansion.nu.len() < 2 {
        return f64::NAN;
    }
    let last = *expansion.nu.last().unwrap();
    let prev = expansion.nu[expansion.nu.len() - 2];
    if prev.abs() < 1e-14 {
        return f64::NAN;
    }
    ((last - prev) / prev).abs()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wilson_fisher_fixed_point_ising() {
        // n=1 (Ising), ε=1 (d=3): g* = 48π²/9
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        let g_star = EpsilonExpansion::wilson_fisher_fixed_point_1loop(1, 1.0);
        let expected = 48.0 * pi2 / 9.0;
        assert!((g_star - expected).abs() < 1e-10);
    }

    #[test]
    fn test_wilson_fisher_vanishes_at_eps0() {
        let g_star = EpsilonExpansion::wilson_fisher_fixed_point_1loop(1, 0.0);
        assert!(g_star.abs() < 1e-10);
    }

    #[test]
    fn test_eta_o_eps2_order() {
        // η ~ O(ε²), so at ε=0 it should vanish
        let eta = EpsilonExpansion::eta_o_eps2(1, 0.0);
        assert!(eta.abs() < 1e-10);
    }

    #[test]
    fn test_nu_o_eps1_mean_field() {
        // At ε=0, ν=0.5 (mean field)
        let nu = EpsilonExpansion::nu_o_eps1(1, 0.0);
        assert!((nu - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_nu_o_eps1_increases_with_eps() {
        let nu0 = EpsilonExpansion::nu_o_eps1(1, 0.0);
        let nu1 = EpsilonExpansion::nu_o_eps1(1, 1.0);
        assert!(nu1 > nu0);
    }

    #[test]
    fn test_full_epsilon_expansion() {
        let exp = EpsilonExpansion::compute(1, 0.5);
        assert_eq!(exp.n, 1);
        assert!((exp.epsilon - 0.5).abs() < 1e-10);
        assert!(exp.g_star[1] > 0.0);
    }

    #[test]
    fn test_epsilon_expansion_critical_exponents() {
        let exp = EpsilonExpansion::compute(1, 1.0);
        let ce = exp.critical_exponents();
        // ν should be close to 0.5 + something
        assert!(ce.nu > 0.5);
        assert!(ce.nu < 1.0);
    }

    #[test]
    fn test_epsilon_scan() {
        let eps_vals: Vec<f64> = (1..=10).map(|i| i as f64 / 10.0).collect();
        let results = epsilon_scan(1, &eps_vals);
        assert_eq!(results.len(), 10);
        // g* should increase with ε
        for i in 1..results.len() {
            assert!(results[i].1 > results[i - 1].1);
        }
    }

    #[test]
    fn test_convergence_check() {
        let exp = EpsilonExpansion::compute(1, 0.1);
        let conv = convergence_check(&exp);
        assert!(conv.is_finite());
        assert!(conv >= 0.0);
    }

    #[test]
    fn test_o_n_symmetry_n_dependence() {
        // Higher n → smaller g* at same ε
        let g1 = EpsilonExpansion::wilson_fisher_fixed_point_1loop(1, 1.0);
        let g3 = EpsilonExpansion::wilson_fisher_fixed_point_1loop(3, 1.0);
        assert!(g1 > g3);
    }
}
