//! Scaling relations and hyperscaling.
//!
//! Near a critical point, thermodynamic quantities follow power laws.
//! The critical exponents are not independent — they satisfy scaling relations.

use serde::{Deserialize, Serialize};

use crate::critical_exponents::CriticalExponents;

/// Scaling dimension of an operator.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalingDimension {
    /// Name of the operator.
    pub name: String,
    /// Scaling dimension x.
    pub dimension: f64,
    /// Whether it's relevant (x < d), irrelevant (x > d), or marginal (x = d).
    pub relevance: crate::beta::CouplingRelevance,
}

impl ScalingDimension {
    /// Create from a scaling dimension x in spatial dimension d.
    pub fn from_dimension(name: &str, x: f64, d: f64) -> Self {
        let relevance = if x < d - 1e-10 {
            crate::beta::CouplingRelevance::Relevant
        } else if x > d + 1e-10 {
            crate::beta::CouplingRelevance::Irrelevant
        } else {
            crate::beta::CouplingRelevance::Marginal
        };
        ScalingDimension {
            name: name.to_string(),
            dimension: x,
            relevance,
        }
    }
}

/// Verify all standard scaling relations for a set of critical exponents.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScalingVerification {
    /// Rushbrooke: α + 2β + γ = 2.
    pub rushbrooke_error: f64,
    /// Widom: γ = β(δ - 1).
    pub widom_error: f64,
    /// Fisher: γ = ν(2 - η).
    pub fisher_error: f64,
    /// Josephson hyperscaling: 2 - α = dν.
    pub hyperscaling_error: f64,
    /// Whether all relations are satisfied within tolerance.
    pub all_satisfied: bool,
}

impl ScalingVerification {
    /// Verify all scaling relations with given tolerance.
    pub fn verify(exponents: &CriticalExponents, tol: f64) -> Self {
        let r = exponents.verify_rushbrooke();
        let w = exponents.verify_widom();
        let f = exponents.verify_fisher();
        let h = exponents.verify_hyperscaling();

        ScalingVerification {
            rushbrooke_error: r,
            widom_error: w,
            fisher_error: f,
            hyperscaling_error: h,
            all_satisfied: r.abs() < tol && w.abs() < tol && f.abs() < tol && h.abs() < tol,
        }
    }
}

/// Compute the correlation length: ξ ∼ |t|^{-ν}.
pub fn correlation_length(nu: f64, t: f64, xi0: f64) -> f64 {
    xi0 * t.abs().powf(-nu)
}

/// Compute the order parameter: M ∼ (-t)^β for t < 0.
pub fn order_parameter(beta: f64, t: f64, m0: f64) -> f64 {
    if t >= 0.0 {
        0.0
    } else {
        m0 * (-t).powf(beta)
    }
}

/// Compute the susceptibility: χ ∼ |t|^{-γ}.
pub fn susceptibility(gamma: f64, t: f64, chi0: f64) -> f64 {
    chi0 * t.abs().powf(-gamma)
}

/// Compute the specific heat: C ∼ |t|^{-α}.
pub fn specific_heat(alpha: f64, t: f64, c0: f64) -> f64 {
    if alpha.abs() < 1e-10 {
        // Logarithmic divergence
        c0 * t.abs().ln().abs()
    } else {
        c0 * t.abs().powf(-alpha)
    }
}

/// Finite-size scaling: at Tc, the correlation length is cut off by system size L.
/// Observable O scales as O(L) ∼ L^{x_O/ν}.
pub fn finite_size_scaling(exponent: f64, nu: f64, l: f64, amplitude: f64) -> f64 {
    amplitude * l.powf(exponent / nu)
}

/// Scaling collapse: near the critical point, data for different system sizes
/// should collapse when plotted as O*L^{-β/ν} vs t*L^{1/ν}.
pub fn scaling_collapse_transform(
    t: f64,
    observable: f64,
    l: f64,
    nu: f64,
    beta_exp: f64,
) -> (f64, f64) {
    let x = t * l.powf(1.0 / nu);
    let y = observable * l.powf(-beta_exp / nu);
    (x, y)
}

/// Corrections to scaling: leading correction exponent ω.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrectionsToScaling {
    /// Leading correction-to-scaling exponent ω.
    pub omega: f64,
    /// Subleading correction exponent.
    pub omega_sub: f64,
}

impl CorrectionsToScaling {
    /// For the 3D Ising model: ω ≈ 0.832.
    pub fn ising_3d() -> Self {
        CorrectionsToScaling {
            omega: 0.832,
            omega_sub: 1.67,
        }
    }

    /// Apply correction to a scaling observable.
    /// O(t) = O₀ |t|^{-x} (1 + a₁|t|^ω + a₂|t|^(ω_sub) + ...)
    pub fn apply_correction(&self, base_value: f64, t: f64, a1: f64, a2: f64) -> f64 {
        base_value * (1.0 + a1 * t.abs().powf(self.omega) + a2 * t.abs().powf(self.omega_sub))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_correlation_length_divergence() {
        // ξ → ∞ as t → 0
        let xi1 = correlation_length(0.63, 0.1, 1.0);
        let xi2 = correlation_length(0.63, 0.01, 1.0);
        assert!(xi2 > xi1);
    }

    #[test]
    fn test_correlation_length_power_law() {
        // ξ ∼ |t|^{-ν}
        let xi1 = correlation_length(0.5, 0.1, 1.0);
        let xi2 = correlation_length(0.5, 0.01, 1.0);
        // Ratio should be (0.1/0.01)^{0.5} = 10^{0.5} = √10
        let expected_ratio = (0.1_f64 / 0.01_f64).powf(0.5);
        assert!((xi2 / xi1 - expected_ratio).abs() < 0.01);
    }

    #[test]
    fn test_order_parameter_zero_above_tc() {
        assert_eq!(order_parameter(0.33, 0.1, 1.0), 0.0);
    }

    #[test]
    fn test_order_parameter_nonzero_below_tc() {
        let m = order_parameter(0.33, -0.1, 1.0);
        assert!(m > 0.0);
    }

    #[test]
    fn test_susceptibility_divergence() {
        let chi1 = susceptibility(1.24, 0.1, 1.0);
        let chi2 = susceptibility(1.24, 0.01, 1.0);
        assert!(chi2 > chi1);
    }

    #[test]
    fn test_specific_heat_log_divergence() {
        // α=0 → logarithmic divergence
        let c = specific_heat(0.0, 0.01, 1.0);
        assert!(c > 0.0);
    }

    #[test]
    fn test_scaling_verification_mean_field() {
        let mf = CriticalExponents::mean_field();
        let ver = ScalingVerification::verify(&mf, 1e-10);
        assert!(ver.all_satisfied);
    }

    #[test]
    fn test_finite_size_scaling() {
        // At Tc, M(L) ∼ L^{-β/ν} — exponent is negative so M decreases with L
        let m_l1 = finite_size_scaling(-0.326, 0.63, 10.0, 1.0);
        let m_l2 = finite_size_scaling(-0.326, 0.63, 100.0, 1.0);
        assert!(m_l2 < m_l1); // order parameter decreases with L at critical point
    }

    #[test]
    fn test_scaling_collapse_transform() {
        let (x, y) = scaling_collapse_transform(0.01, 1.5, 100.0, 0.63, 0.326);
        assert!(x.is_finite());
        assert!(y.is_finite());
    }

    #[test]
    fn test_corrections_to_scaling() {
        let corr = CorrectionsToScaling::ising_3d();
        let base = 1.0;
        let corrected = corr.apply_correction(base, 0.01, 0.1, 0.01);
        assert!(corrected > base);
    }

    #[test]
    fn test_scaling_dimension_relevant() {
        let sd = ScalingDimension::from_dimension("thermal", 1.0, 3.0);
        assert_eq!(sd.relevance, crate::beta::CouplingRelevance::Relevant);
    }

    #[test]
    fn test_scaling_dimension_irrelevant() {
        let sd = ScalingDimension::from_dimension("higher", 5.0, 3.0);
        assert_eq!(sd.relevance, crate::beta::CouplingRelevance::Irrelevant);
    }

    #[test]
    fn test_scaling_dimension_marginal() {
        let sd = ScalingDimension::from_dimension("marginal", 3.0, 3.0);
        assert_eq!(sd.relevance, crate::beta::CouplingRelevance::Marginal);
    }
}
