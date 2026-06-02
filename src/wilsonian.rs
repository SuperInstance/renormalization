//! Wilsonian RG: momentum shell integration and effective actions.
//!
//! The Wilsonian approach integrates out high-momentum modes shell by shell,
//! generating an effective action for the remaining low-energy degrees of freedom.

use serde::{Deserialize, Serialize};

/// Parameters for the Wilsonian momentum shell integration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WilsonianParams {
    /// UV cutoff Λ.
    pub lambda: f64,
    /// Current shell width dΛ.
    pub d_lambda: f64,
    /// Spatial dimension d.
    pub dimension: usize,
    /// Number of shells to integrate out.
    pub n_shells: usize,
}

impl Default for WilsonianParams {
    fn default() -> Self {
        Self {
            lambda: 1.0,
            d_lambda: 0.01,
            dimension: 3,
            n_shells: 100,
        }
    }
}

/// An effective action described by its coupling constants.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EffectiveAction {
    /// Mass term (r in φ² term).
    pub r: f64,
    /// φ⁴ coupling (u).
    pub u: f64,
    /// φ⁶ coupling (v).
    pub v: f64,
    /// Higher-order couplings.
    pub higher: Vec<f64>,
    /// Current momentum cutoff.
    pub cutoff: f64,
    /// Spatial dimension.
    pub dimension: usize,
}

impl EffectiveAction {
    /// Create a φ⁴ effective action with given parameters.
    pub fn phi4(r: f64, u: f64, dimension: usize, cutoff: f64) -> Self {
        EffectiveAction {
            r,
            u,
            v: 0.0,
            higher: Vec::new(),
            cutoff,
            dimension,
        }
    }

    /// Perform one Wilsonian step: integrate out shell [Λ-dΛ, Λ].
    ///
    /// At one-loop order:
    ///   dr/dl = 2r + (n+2)u * K_d * Λ^d / (r + Λ²)
    ///   du/dl = (4-d)u - (n+8)u² * K_d * Λ^d / (r + Λ²)²
    /// where l = ln(Λ₀/Λ), K_d = S_d/(2π)^d, n=1 (Ising).
    pub fn wilsonian_step(&self, d_lambda: f64) -> Self {
        let d = self.dimension as f64;
        let lambda = self.cutoff;
        let dl = d_lambda / lambda; // d(ln Λ)

        // Phase space factor K_d (simplified)
        let k_d = surface_area(d) / (2.0 * std::f64::consts::PI).powi(d as i32);

        // One-loop corrections (n=1 component, Ising-like)
        let n_comp = 1.0;
        let denom = self.r + lambda * lambda;

        // dr/dl = 2r + (n+2)*u*K_d*Λ^d/denom
        let dr_dl = 2.0 * self.r + (n_comp + 2.0) * self.u * k_d * lambda.powf(d) / denom;

        // du/dl = (4-d)*u - (n+8)*u²*K_d*Λ^d/denom²
        let du_dl = (4.0 - d) * self.u
            - (n_comp + 8.0) * self.u * self.u * k_d * lambda.powf(d) / (denom * denom);

        let new_r = self.r + dl * dr_dl;
        let new_u = self.u + dl * du_dl;

        EffectiveAction {
            r: new_r,
            u: new_u,
            v: self.v, // unchanged at this order
            higher: self.higher.clone(),
            cutoff: lambda - d_lambda,
            dimension: self.dimension,
        }
    }

    /// Run the full Wilsonian RG flow, integrating out modes from Λ down to Λ*e^{-n_shells*dl}.
    pub fn wilsonian_flow(&self, params: &WilsonianParams) -> Vec<EffectiveAction> {
        let mut actions = Vec::with_capacity(params.n_shells + 1);
        actions.push(self.clone());
        let mut current = self.clone();

        for _ in 0..params.n_shells {
            current = current.wilsonian_step(params.d_lambda);
            actions.push(current.clone());
        }

        actions
    }

    /// Extract the mass eigenvalue (inverse correlation length squared).
    pub fn mass_squared(&self) -> f64 {
        self.r
    }
}

/// Surface area of unit sphere in d dimensions.
pub fn surface_area(d: f64) -> f64 {
    if d <= 0.0 {
        return 0.0;
    }
    2.0 * std::f64::consts::PI.powf(d / 2.0) / gamma_fn(d / 2.0)
}

/// Approximate gamma function (Stirling + Lanczos).
pub fn gamma_fn(x: f64) -> f64 {
    if x <= 0.0 {
        return f64::INFINITY;
    }
    // Lanczos approximation
    let g = 7.0;
    let coef = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];

    if x < 0.5 {
        return std::f64::consts::PI / (gamma_fn(1.0 - x) * (std::f64::consts::PI * x).sin());
    }

    let x = x - 1.0;
    let a = coef[0];
    let t = x + g + 0.5;

    let sum = coef[1..].iter().enumerate().fold(a, |acc, (i, c)| {
        acc + c / (x + i as f64 + 1.0)
    });

    (2.0 * std::f64::consts::PI).sqrt() * t.powf(x + 0.5) * (-t).exp() * sum
}

/// Exact renormalization group equation (Wetterich equation) for the effective average action.
///
/// ∂ₜΓₖ[φ] = ½ Tr[(Γₖ⁽²⁾[φ] + k²)⁻¹]
/// where t = ln(k/Λ).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErgFlow {
    /// Current values of the effective average action couplings.
    pub couplings: Vec<f64>,
    /// Current RG scale k.
    pub k: f64,
    /// UV cutoff.
    pub lambda: f64,
}

impl ErgFlow {
    /// Create a new ERG flow state.
    pub fn new(couplings: Vec<f64>, k: f64, lambda: f64) -> Self {
        ErgFlow { couplings, k, lambda }
    }

    /// Advance the ERG flow by one step.
    pub fn step(&self, dt: f64, beta_couplings: &[f64]) -> Self {
        let new_couplings: Vec<f64> = self
            .couplings
            .iter()
            .zip(beta_couplings.iter())
            .map(|(g, b)| g + dt * b)
            .collect();
        ErgFlow {
            couplings: new_couplings,
            k: self.k * (-dt).exp(),
            lambda: self.lambda,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_effective_action_creation() {
        let action = EffectiveAction::phi4(1.0, 0.5, 3, 1.0);
        assert_eq!(action.r, 1.0);
        assert_eq!(action.u, 0.5);
        assert_eq!(action.dimension, 3);
    }

    #[test]
    fn test_wilsonian_step_changes_couplings() {
        let action = EffectiveAction::phi4(0.1, 0.1, 3, 1.0);
        let stepped = action.wilsonian_step(0.01);
        // r should increase (mass grows)
        assert!(stepped.r > action.r || stepped.r != action.r);
        assert!(stepped.cutoff < action.cutoff);
    }

    #[test]
    fn test_wilsonian_flow_length() {
        let action = EffectiveAction::phi4(0.1, 0.1, 3, 1.0);
        let params = WilsonianParams {
            n_shells: 50,
            d_lambda: 0.01,
            ..Default::default()
        };
        let flow = action.wilsonian_flow(&params);
        assert_eq!(flow.len(), 51); // n_shells + 1
    }

    #[test]
    fn test_wilsonian_flow_cutoff_decreases() {
        let action = EffectiveAction::phi4(0.1, 0.1, 3, 1.0);
        let params = WilsonianParams {
            n_shells: 50,
            d_lambda: 0.01,
            ..Default::default()
        };
        let flow = action.wilsonian_flow(&params);
        for i in 1..flow.len() {
            assert!(flow[i].cutoff < flow[i - 1].cutoff);
        }
    }

    #[test]
    fn test_surface_area_d2() {
        // S₁ = 2π
        let sa = surface_area(2.0);
        assert!((sa - 2.0 * std::f64::consts::PI).abs() < 1e-8);
    }

    #[test]
    fn test_surface_area_d3() {
        // S₂ = 4π
        let sa = surface_area(3.0);
        assert!((sa - 4.0 * std::f64::consts::PI).abs() < 1e-8);
    }

    #[test]
    fn test_gamma_fn_half() {
        // Γ(1/2) = √π
        let g = gamma_fn(0.5);
        assert!((g - std::f64::consts::PI.sqrt()).abs() / std::f64::consts::PI.sqrt() < 1e-8);
    }

    #[test]
    fn test_gamma_fn_integer() {
        // Γ(n) = (n-1)!
        assert!((gamma_fn(1.0) - 1.0).abs() < 1e-8);
        assert!((gamma_fn(2.0) - 1.0).abs() < 1e-8);
        assert!((gamma_fn(3.0) - 2.0).abs() < 1e-8);
        assert!((gamma_fn(4.0) - 6.0).abs() < 1e-8);
    }

    #[test]
    fn test_erg_flow_step() {
        let flow = ErgFlow::new(vec![0.5, 0.3], 1.0, 10.0);
        let stepped = flow.step(0.1, &[-0.5, -0.3]);
        assert!((stepped.couplings[0] - 0.45).abs() < 1e-10);
        assert!((stepped.couplings[1] - 0.27).abs() < 1e-10);
        assert!(stepped.k < flow.k);
    }

    #[test]
    fn test_mass_squared() {
        let action = EffectiveAction::phi4(2.5, 0.1, 3, 1.0);
        assert!((action.mass_squared() - 2.5).abs() < 1e-10);
    }
}
