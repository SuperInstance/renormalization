//! Critical exponents from linearized RG.
//!
//! Near a fixed point g*, the linearized RG gives:
//!   δgᵢ → b^{yᵢ} δgᵢ
//! where yᵢ are the eigenvalues of the stability matrix.
//! Critical exponents are related to these eigenvalues.

use serde::{Deserialize, Serialize};



/// Critical exponents for a phase transition.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CriticalExponents {
    /// Correlation length exponent ν: ξ ∼ |t|^{-ν}.
    /// Related to the thermal eigenvalue by ν = 1/yₜ.
    pub nu: f64,
    /// Order parameter exponent β (not to be confused with beta function).
    /// M ∼ (-t)^β for t < 0.
    pub beta: f64,
    /// Susceptibility exponent γ: χ ∼ |t|^{-γ}.
    pub gamma: f64,
    /// Critical isotherm exponent δ: M ∼ h^{1/δ} at t=0.
    pub delta: f64,
    /// Correlation function exponent η: G(r) ∼ r^{-(d-2+η)} at Tc.
    pub eta: f64,
    /// Specific heat exponent α: C ∼ |t|^{-α}.
    pub alpha: f64,
    /// Spatial dimension.
    pub dimension: f64,
}

impl CriticalExponents {
    /// Verify the Rushbrooke scaling relation: α + 2β + γ = 2.
    pub fn verify_rushbrooke(&self) -> f64 {
        self.alpha + 2.0 * self.beta + self.gamma - 2.0
    }

    /// Verify the Widom scaling relation: γ = β(δ - 1).
    pub fn verify_widom(&self) -> f64 {
        self.gamma - self.beta * (self.delta - 1.0)
    }

    /// Verify the Fisher scaling relation: γ = ν(2 - η).
    pub fn verify_fisher(&self) -> f64 {
        self.gamma - self.nu * (2.0 - self.eta)
    }

    /// Verify the Josephson hyperscaling relation: 2 - α = dν.
    pub fn verify_hyperscaling(&self) -> f64 {
        2.0 - self.alpha - self.dimension * self.nu
    }

    /// Compute ν from the thermal eigenvalue yₜ.
    /// ν = 1/yₜ
    pub fn nu_from_eigenvalue(y_t: f64) -> f64 {
        1.0 / y_t
    }

    /// Mean-field critical exponents (valid above upper critical dimension).
    pub fn mean_field() -> Self {
        CriticalExponents {
            nu: 0.5,
            beta: 0.5,
            gamma: 1.0,
            delta: 3.0,
            eta: 0.0,
            alpha: 0.0,
            dimension: 4.0,
        }
    }

    /// Ising 2D exact exponents.
    pub fn ising_2d() -> Self {
        CriticalExponents {
            nu: 1.0,
            beta: 0.125,
            gamma: 1.75,
            delta: 15.0,
            eta: 0.25,
            alpha: 0.0, // log divergence
            dimension: 2.0,
        }
    }

    /// Ising 3D approximate exponents.
    pub fn ising_3d() -> Self {
        CriticalExponents {
            nu: 0.630,
            beta: 0.326,
            gamma: 1.237,
            delta: 4.789,
            eta: 0.036,
            alpha: 0.110,
            dimension: 3.0,
        }
    }
}

/// Compute critical exponents from a 1D fixed point eigenvalue.
///
/// Given the eigenvalue λ = β'(g*) at a fixed point:
/// - The scaling dimension y = -λ
/// - ν = 1/y (for y > 0, i.e., stable fixed point with t being relevant)
pub fn exponents_from_eigenvalue(eigenvalue: f64, dimension: f64) -> Result<CriticalExponents, String> {
    if eigenvalue >= 0.0 {
        return Err("Eigenvalue must be negative for stable fixed point to compute ν".to_string());
    }
    let y_t = -eigenvalue;
    let nu = 1.0 / y_t;
    let alpha = 2.0 - dimension * nu;

    let gamma = nu * 2.0; // Fisher relation with η=0

    Ok(CriticalExponents {
        nu,
        beta: (2.0 - alpha - gamma) / 2.0,
        gamma,
        delta: 1.0 + gamma / ((2.0 - alpha - gamma) / 2.0),
        eta: 0.0,
        alpha,
        dimension,
    })
}

/// Compute ν from the leading thermal eigenvalue of the stability matrix.
pub fn nu_from_thermal_eigenvalue(y_t: f64) -> f64 {
    1.0 / y_t
}

/// Compute α from ν and dimension d via hyperscaling: α = 2 - dν.
pub fn alpha_from_nu(nu: f64, d: f64) -> f64 {
    2.0 - d * nu
}

/// Compute γ from ν and η via Fisher: γ = ν(2 - η).
pub fn gamma_from_nu_eta(nu: f64, eta: f64) -> f64 {
    nu * (2.0 - eta)
}

/// Compute β_exponent from α and γ via Rushbrooke: β = (2 - α - γ)/2.
pub fn beta_from_alpha_gamma(alpha: f64, gamma: f64) -> f64 {
    (2.0 - alpha - gamma) / 2.0
}

/// Compute δ from β and γ via Widom: δ = 1 + γ/β.
pub fn delta_from_beta_gamma(beta: f64, gamma: f64) -> f64 {
    1.0 + gamma / beta
}

/// Full set of exponents from ν, η, and dimension using scaling relations.
pub fn exponents_from_nu_eta(nu: f64, eta: f64, d: f64) -> CriticalExponents {
    let gamma = gamma_from_nu_eta(nu, eta);
    let alpha = alpha_from_nu(nu, d);
    let beta = beta_from_alpha_gamma(alpha, gamma);
    let delta = delta_from_beta_gamma(beta, gamma);

    CriticalExponents {
        nu,
        beta,
        gamma,
        delta,
        eta,
        alpha,
        dimension: d,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mean_field_rushbrooke() {
        let mf = CriticalExponents::mean_field();
        assert!(mf.verify_rushbrooke().abs() < 1e-10);
    }

    #[test]
    fn test_mean_field_widom() {
        let mf = CriticalExponents::mean_field();
        assert!(mf.verify_widom().abs() < 1e-10);
    }

    #[test]
    fn test_mean_field_fisher() {
        let mf = CriticalExponents::mean_field();
        assert!(mf.verify_fisher().abs() < 1e-10);
    }

    #[test]
    fn test_mean_field_hyperscaling() {
        let mf = CriticalExponents::mean_field();
        assert!(mf.verify_hyperscaling().abs() < 1e-10);
    }

    #[test]
    fn test_ising_2d_scaling_relations() {
        let ising = CriticalExponents::ising_2d();
        // Rushbrooke: α + 2β + γ = 2 (α=0 for log)
        assert!((ising.alpha + 2.0 * ising.beta + ising.gamma - 2.0).abs() < 0.01);
        // Fisher: γ = ν(2-η)
        assert!((ising.gamma - ising.nu * (2.0 - ising.eta)).abs() < 0.01);
    }

    #[test]
    fn test_nu_from_eigenvalue() {
        // Mean field: yₜ = 2 → ν = 0.5
        let nu = nu_from_thermal_eigenvalue(2.0);
        assert!((nu - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_exponents_from_nu_eta_mean_field() {
        let exponents = exponents_from_nu_eta(0.5, 0.0, 4.0);
        let mf = CriticalExponents::mean_field();
        assert!((exponents.alpha - mf.alpha).abs() < 1e-10);
        assert!((exponents.gamma - mf.gamma).abs() < 1e-10);
    }

    #[test]
    fn test_alpha_from_nu() {
        // d=3, ν=0.63 → α = 2-1.89 = 0.11
        let alpha = alpha_from_nu(0.63, 3.0);
        assert!((alpha - 0.11).abs() < 0.01);
    }

    #[test]
    fn test_gamma_from_nu_eta() {
        let gamma = gamma_from_nu_eta(1.0, 0.25);
        assert!((gamma - 1.75).abs() < 1e-10);
    }

    #[test]
    fn test_ising_3d_approximate_relations() {
        let ising = CriticalExponents::ising_3d();
        // These are approximate, so use loose tolerance
        assert!(ising.verify_rushbrooke().abs() < 0.01);
        assert!(ising.verify_fisher().abs() < 0.01);
    }
}
