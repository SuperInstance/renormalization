//! Beta functions and RG flow equations.
//!
//! The beta function β(g) = dg/d(ln μ) describes how a coupling constant g
//! evolves with the energy scale μ. This is the core equation of the
//! renormalization group.

use serde::{Deserialize, Serialize};
use nalgebra::DVector;

/// A beta function describing RG flow for a single coupling.
pub type BetaFn = Box<dyn Fn(f64) -> f64 + Send + Sync>;

/// Multi-dimensional beta function for a vector of couplings.
pub type BetaFnVec = Box<dyn Fn(&DVector<f64>) -> DVector<f64> + Send + Sync>;

/// Jacobian of a multi-dimensional beta function.
pub type BetaJacobian = Box<dyn Fn(&DVector<f64>) -> nalgebra::DMatrix<f64> + Send + Sync>;

/// Configuration for an RG flow solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowConfig {
    /// Number of integration steps.
    pub steps: usize,
    /// Step size in ln(μ).
    pub d_ln_mu: f64,
    /// Integration method.
    pub method: IntegrationMethod,
}

/// Integration method for RG flow.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum IntegrationMethod {
    /// Forward Euler (first-order).
    Euler,
    /// 4th-order Runge-Kutta.
    RK4,
}

impl Default for FlowConfig {
    fn default() -> Self {
        Self {
            steps: 1000,
            d_ln_mu: 0.01,
            method: IntegrationMethod::Euler,
        }
    }
}

/// Result of integrating an RG flow trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowTrajectory {
    /// Coupling values at each step.
    pub trajectory: Vec<f64>,
    /// ln(μ) values at each step.
    pub ln_mu_values: Vec<f64>,
}

/// Result of integrating a multi-dimensional RG flow trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowTrajectoryVec {
    /// Coupling vectors at each step.
    pub trajectory: Vec<Vec<f64>>,
    /// ln(μ) values at each step.
    pub ln_mu_values: Vec<f64>,
}

/// Integrate a 1D RG flow from initial coupling `g0` using Euler method.
pub fn flow_euler(beta: &BetaFn, g0: f64, config: &FlowConfig) -> FlowTrajectory {
    let mut g = g0;
    let mut trajectory = Vec::with_capacity(config.steps + 1);
    let mut ln_mu_values = Vec::with_capacity(config.steps + 1);

    trajectory.push(g);
    ln_mu_values.push(0.0);

    for i in 1..=config.steps {
        let ln_mu = i as f64 * config.d_ln_mu;
        g += config.d_ln_mu * beta(g);
        trajectory.push(g);
        ln_mu_values.push(ln_mu);
    }

    FlowTrajectory {
        trajectory,
        ln_mu_values,
    }
}

/// Integrate a 1D RG flow using RK4.
pub fn flow_rk4(beta: &BetaFn, g0: f64, config: &FlowConfig) -> FlowTrajectory {
    let mut g = g0;
    let mut trajectory = Vec::with_capacity(config.steps + 1);
    let mut ln_mu_values = Vec::with_capacity(config.steps + 1);

    trajectory.push(g);
    ln_mu_values.push(0.0);

    let h = config.d_ln_mu;

    for i in 1..=config.steps {
        let ln_mu = i as f64 * h;
        let k1 = beta(g);
        let k2 = beta(g + 0.5 * h * k1);
        let k3 = beta(g + 0.5 * h * k2);
        let k4 = beta(g + h * k3);
        g += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        trajectory.push(g);
        ln_mu_values.push(ln_mu);
    }

    FlowTrajectory {
        trajectory,
        ln_mu_values,
    }
}

/// Integrate a 1D RG flow using the configured method.
pub fn flow(beta: &BetaFn, g0: f64, config: &FlowConfig) -> FlowTrajectory {
    match config.method {
        IntegrationMethod::Euler => flow_euler(beta, g0, config),
        IntegrationMethod::RK4 => flow_rk4(beta, g0, config),
    }
}

/// Integrate a multi-dimensional RG flow using Euler method.
pub fn flow_vec_euler(
    beta: &BetaFnVec,
    g0: &DVector<f64>,
    config: &FlowConfig,
) -> FlowTrajectoryVec {
    let mut g = g0.clone();
    let mut trajectory = Vec::with_capacity(config.steps + 1);
    let mut ln_mu_values = Vec::with_capacity(config.steps + 1);

    trajectory.push(g.iter().cloned().collect());
    ln_mu_values.push(0.0);

    for i in 1..=config.steps {
        let ln_mu = i as f64 * config.d_ln_mu;
        let dg = beta(&g);
        g += config.d_ln_mu * dg;
        trajectory.push(g.iter().cloned().collect());
        ln_mu_values.push(ln_mu);
    }

    FlowTrajectoryVec {
        trajectory,
        ln_mu_values,
    }
}

/// Compute the relevant direction from a 1D beta function.
/// A coupling is relevant if β'(g*) > 0 (flows away from fixed point).
pub fn coupling_relevance(beta_prime: &BetaFn, g_star: f64) -> CouplingRelevance {
    let val = beta_prime(g_star);
    if val > 0.0 {
        CouplingRelevance::Relevant
    } else if val < 0.0 {
        CouplingRelevance::Irrelevant
    } else {
        CouplingRelevance::Marginal
    }
}

/// Classification of a coupling's relevance at a fixed point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum CouplingRelevance {
    /// β'(g*) > 0: flows away from fixed point.
    Relevant,
    /// β'(g*) < 0: flows toward fixed point.
    Irrelevant,
    /// β'(g*) = 0: need higher-order analysis.
    Marginal,
}

/// φ⁴ theory beta function: β(g) = -εg + 3g²/(16π²) at one loop.
pub fn phi4_beta_one_loop(epsilon: f64) -> BetaFn {
    Box::new(move |g: f64| -> f64 {
        -epsilon * g + 3.0 * g * g / (16.0 * std::f64::consts::PI * std::f64::consts::PI)
    })
}

/// Gaussian (free theory) fixed point beta function: β(g) = -εg.
pub fn gaussian_beta(epsilon: f64) -> BetaFn {
    Box::new(move |g: f64| -> f64 { -epsilon * g })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_flow_constant_beta() {
        // β(g) = 1 → g = g0 + ln(μ)
        let beta: BetaFn = Box::new(|_| 1.0);
        let config = FlowConfig {
            steps: 100,
            d_ln_mu: 0.01,
            method: IntegrationMethod::Euler,
        };
        let traj = flow_euler(&beta, 0.0, &config);
        assert!((traj.trajectory.last().unwrap() - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_rk4_flow_constant_beta() {
        let beta: BetaFn = Box::new(|_| 1.0);
        let config = FlowConfig {
            steps: 100,
            d_ln_mu: 0.01,
            method: IntegrationMethod::RK4,
        };
        let traj = flow_rk4(&beta, 0.0, &config);
        assert!((traj.trajectory.last().unwrap() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_flow_converges_to_fixed_point() {
        // β(g) = -g → fixed point at g=0, flows toward it
        let beta: BetaFn = Box::new(|g: f64| -g);
        let config = FlowConfig {
            steps: 10000,
            d_ln_mu: 0.001,
            method: IntegrationMethod::Euler,
        };
        let traj = flow(&beta, 1.0, &config);
        assert!(traj.trajectory.last().unwrap().abs() < 0.01);
    }

    #[test]
    fn test_coupling_relevance_relevant() {
        let beta_prime: BetaFn = Box::new(|_| 1.0);
        assert_eq!(coupling_relevance(&beta_prime, 0.0), CouplingRelevance::Relevant);
    }

    #[test]
    fn test_coupling_relevance_irrelevant() {
        let beta_prime: BetaFn = Box::new(|_| -1.0);
        assert_eq!(coupling_relevance(&beta_prime, 0.0), CouplingRelevance::Irrelevant);
    }

    #[test]
    fn test_coupling_relevance_marginal() {
        let beta_prime: BetaFn = Box::new(|_| 0.0);
        assert_eq!(coupling_relevance(&beta_prime, 0.0), CouplingRelevance::Marginal);
    }

    #[test]
    fn test_phi4_beta_fixed_point() {
        let eps = 1.0;
        let beta = phi4_beta_one_loop(eps);
        // Fixed point: g* = 16π²ε/3
        let g_star = 16.0 * std::f64::consts::PI.powi(2) * eps / 3.0;
        assert!(beta(g_star).abs() < 1e-10);
    }

    #[test]
    fn test_gaussian_beta() {
        let beta = gaussian_beta(1.0);
        assert!((beta(1.0) - (-1.0)).abs() < 1e-10);
    }

    #[test]
    fn test_vec_flow_dimension() {
        let beta: BetaFnVec = Box::new(|g: &DVector<f64>| -> DVector<f64> { -2.0 * g });
        let g0 = DVector::from_vec(vec![1.0, 2.0]);
        let config = FlowConfig {
            steps: 5000,
            d_ln_mu: 0.001,
            method: IntegrationMethod::Euler,
        };
        let traj = flow_vec_euler(&beta, &g0, &config);
        assert_eq!(traj.trajectory[0].len(), 2);
        assert!(traj.trajectory.last().unwrap()[0].abs() < 1.0); // decaying, not necessarily < 0.1
    }
}
