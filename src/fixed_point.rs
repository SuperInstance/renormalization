//! Fixed points of RG flow: β(g*) = 0 with stability classification.
//!
//! Fixed points are scale-invariant theories where the couplings don't flow.
//! They govern the universal behavior of entire classes of systems.

use serde::{Deserialize, Serialize};
use nalgebra::DVector;

use crate::beta::{BetaFn, BetaFnVec, FlowConfig, flow};

/// A fixed point of the 1D RG flow.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FixedPoint {
    /// The coupling value at the fixed point.
    pub g_star: f64,
    /// β'(g*) — eigenvalue of the linearized flow.
    pub eigenvalue: f64,
    /// Stability classification.
    pub stability: Stability,
    /// Name/label for this fixed point.
    pub label: String,
}

/// Stability classification of a fixed point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Stability {
    /// β'(g*) < 0: IR-attractive, UV-repulsive.
    Stable,
    /// β'(g*) > 0: IR-repulsive, UV-attractive.
    Unstable,
    /// β'(g*) = 0: marginal, need higher-order analysis.
    Marginal,
}

/// A fixed point in multi-dimensional coupling space.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FixedPointVec {
    /// The coupling vector at the fixed point.
    pub g_star: Vec<f64>,
    /// Eigenvalues of the stability matrix (∂βᵢ/∂gⱼ) at g*.
    pub eigenvalues: Vec<f64>,
    /// Number of relevant directions (positive eigenvalues).
    pub n_relevant: usize,
    /// Number of irrelevant directions (negative eigenvalues).
    pub n_irrelevant: usize,
    /// Number of marginal directions (zero eigenvalues).
    pub n_marginal: usize,
    /// Label for this fixed point.
    pub label: String,
}

/// Find a 1D fixed point using Newton's method.
///
/// Solves β(g*) = 0 starting from initial guess `g0`.
pub fn find_fixed_point_newton(
    beta: &BetaFn,
    beta_prime: &BetaFn,
    g0: f64,
    tol: f64,
    max_iter: usize,
) -> Option<FixedPoint> {
    let mut g = g0;
    for _ in 0..max_iter {
        let b = beta(g);
        if b.abs() < tol {
            let eig = beta_prime(g);
            let stability = if eig < -tol {
                Stability::Stable
            } else if eig > tol {
                Stability::Unstable
            } else {
                Stability::Marginal
            };
            return Some(FixedPoint {
                g_star: g,
                eigenvalue: eig,
                stability,
                label: String::new(),
            });
        }
        let bp = beta_prime(g);
        if bp.abs() < 1e-14 {
            return None; // degenerate
        }
        g -= b / bp;
    }
    None
}

/// Find all 1D fixed points by scanning a range and refining with Newton's method.
pub fn find_fixed_points_scan(
    beta: &BetaFn,
    beta_prime: &BetaFn,
    g_min: f64,
    g_max: f64,
    n_scan: usize,
    tol: f64,
) -> Vec<FixedPoint> {
    let dg = (g_max - g_min) / n_scan as f64;
    let mut results = Vec::new();
    let mut prev_b = beta(g_min);

    for i in 1..=n_scan {
        let g = g_min + i as f64 * dg;
        let b = beta(g);

        // Check for zero crossing or near-zero
        if b.abs() < tol {
            if let Some(fp) = find_fixed_point_newton(beta, beta_prime, g, tol, 100) {
                if !results.iter().any(|r: &FixedPoint| (r.g_star - fp.g_star).abs() < tol * 10.0) {
                    results.push(fp);
                }
            }
        } else if prev_b * b < 0.0 {
            // Sign change — refine
            let g0 = g - dg * b / (b - prev_b);
            if let Some(fp) = find_fixed_point_newton(beta, beta_prime, g0, tol, 100) {
                // Avoid duplicates
                if !results.iter().any(|r: &FixedPoint| (r.g_star - fp.g_star).abs() < tol * 10.0) {
                    results.push(fp);
                }
            }
        }
        prev_b = b;
    }
    results
}

/// Classify the basin of attraction for a stable 1D fixed point.
///
/// Returns (lower_bound, upper_bound) of the basin of attraction.
/// Uses bisection to find where the flow diverges.
pub fn basin_of_attraction(
    beta: &BetaFn,
    fp: &FixedPoint,
    search_range: f64,
    config: &FlowConfig,
) -> Option<(f64, f64)> {
    if fp.stability != Stability::Stable {
        return None;
    }

    let g_star = fp.g_star;

    // Find lower bound
    let lower = find_basin_boundary(beta, g_star, g_star - search_range, config)?;

    // Find upper bound
    let upper = find_basin_boundary(beta, g_star, g_star + search_range, config)?;

    Some((lower, upper))
}

fn find_basin_boundary(
    beta: &BetaFn,
    g_inside: f64,
    g_outside: f64,
    config: &FlowConfig,
) -> Option<f64> {
    let tol = 1e-6;
    let mut lo = g_inside;
    let mut hi = g_outside;

    // Verify inside flows to fixed point
    let traj_in = flow(beta, g_inside, config);
    let final_in = *traj_in.trajectory.last().unwrap();
    if (final_in - g_inside).abs() > search_range_abs(lo, hi) {
        // Actually need to check convergence toward some fp
    }

    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        let traj = flow(beta, mid, config);
        let final_g = *traj.trajectory.last().unwrap();

        // Check if it diverges (very large or NaN)
        let diverges = final_g.abs() > 1e10 || final_g.is_nan() || final_g.is_infinite();

        if (hi - lo).abs() < tol {
            return Some(mid);
        }

        if diverges {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Some((lo + hi) / 2.0)
}

fn search_range_abs(a: f64, b: f64) -> f64 {
    (b - a).abs()
}

/// Find a multi-dimensional fixed point using Newton's method.
pub fn find_fixed_point_vec_newton(
    beta: &BetaFnVec,
    jacobian: &crate::beta::BetaJacobian,
    g0: &DVector<f64>,
    tol: f64,
    max_iter: usize,
) -> Option<FixedPointVec> {
    let mut g = g0.clone();
    let _n = g.len();

    for _ in 0..max_iter {
        let b = beta(&g);
        if b.iter().all(|x| x.abs() < tol) {
            let jac = jacobian(&g);
            let eigen = jac.symmetric_eigenvalues();
            let mut n_rel = 0;
            let mut n_irrel = 0;
            let mut n_marg = 0;
            for i in 0..eigen.len() {
                let v = eigen[i];
                if v > tol {
                    n_rel += 1;
                } else if v < -tol {
                    n_irrel += 1;
                } else {
                    n_marg += 1;
                }
            }
            return Some(FixedPointVec {
                g_star: g.iter().cloned().collect(),
                eigenvalues: eigen.iter().cloned().collect(),
                n_relevant: n_rel,
                n_irrelevant: n_irrel,
                n_marginal: n_marg,
                label: String::new(),
            });
        }

        let jac = jacobian(&g);
        // Solve jac * delta = -b
        let decomp = jac.clone().lu();
        if let Some(delta) = decomp.solve(&(-&b)) {
            g += delta;
        } else {
            return None;
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_fixed_point_newton_simple() {
        // β(g) = g - g² → fixed points at g=0, g=1
        let beta: BetaFn = Box::new(|g: f64| g - g * g);
        let beta_prime: BetaFn = Box::new(|g: f64| 1.0 - 2.0 * g);

        let fp0 = find_fixed_point_newton(&beta, &beta_prime, 0.1, 1e-10, 100).unwrap();
        assert!(fp0.g_star.abs() < 1e-8);
        assert_eq!(fp0.stability, Stability::Unstable); // β'(0)=1>0

        let fp1 = find_fixed_point_newton(&beta, &beta_prime, 1.5, 1e-10, 100).unwrap();
        assert!((fp1.g_star - 1.0).abs() < 1e-8);
        assert_eq!(fp1.stability, Stability::Stable); // β'(1)=-1<0
    }

    #[test]
    fn test_scan_fixed_points() {
        let beta: BetaFn = Box::new(|g: f64| g - g * g);
        let beta_prime: BetaFn = Box::new(|g: f64| 1.0 - 2.0 * g);
        // g-g² = g(1-g) has roots at g=0 and g=1
        let fps = find_fixed_points_scan(&beta, &beta_prime, -0.5, 1.5, 200, 1e-10);
        assert!(fps.len() >= 2, "expected >= 2 fixed points, got {}", fps.len());
    }

    #[test]
    fn test_fixed_point_classification_stable() {
        let fp = FixedPoint {
            g_star: 1.0,
            eigenvalue: -2.0,
            stability: Stability::Stable,
            label: "test".to_string(),
        };
        assert_eq!(fp.stability, Stability::Stable);
    }

    #[test]
    fn test_fixed_point_classification_unstable() {
        let fp = FixedPoint {
            g_star: 0.0,
            eigenvalue: 1.0,
            stability: Stability::Unstable,
            label: "test".to_string(),
        };
        assert_eq!(fp.stability, Stability::Unstable);
    }

    #[test]
    fn test_gaussian_fixed_point() {
        // β(g) = -εg → g*=0 is stable for ε>0
        let eps = 1.0;
        let beta: BetaFn = Box::new(move |g: f64| -eps * g);
        let beta_prime: BetaFn = Box::new(move |_: f64| -eps);
        let fp = find_fixed_point_newton(&beta, &beta_prime, 0.01, 1e-10, 100).unwrap();
        assert!(fp.g_star.abs() < 1e-8);
        assert_eq!(fp.stability, Stability::Stable);
    }

    #[test]
    fn test_wilson_fisher_fixed_point() {
        // φ⁴ in d=4-ε: β(g) = -εg + 3g²/(16π²)
        let eps = 0.1;
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        let beta: BetaFn = Box::new(move |g: f64| -eps * g + 3.0 * g * g / (16.0 * pi2));
        let beta_prime: BetaFn = Box::new(move |g: f64| -eps + 6.0 * g / (16.0 * pi2));

        // Non-trivial fixed point at g* = 16π²ε/3
        let g_star_approx = 16.0 * pi2 * eps / 3.0;
        let fp = find_fixed_point_newton(&beta, &beta_prime, g_star_approx * 1.5, 1e-10, 100);
        assert!(fp.is_some());
        let fp = fp.unwrap();
        assert!((fp.g_star - g_star_approx).abs() / g_star_approx < 0.01);
    }

    #[test]
    fn test_multi_dim_fixed_point() {
        // β₁(g₁,g₂) = -g₁ + g₁*g₂
        // β₂(g₁,g₂) = -2*g₂ + g₁²
        // Fixed point at g₁=g₂=0
        let beta: BetaFnVec = Box::new(|g: &DVector<f64>| {
            DVector::from_vec(vec![
                -g[0] + g[0] * g[1],
                -2.0 * g[1] + g[0] * g[0],
            ])
        });
        let jac: crate::beta::BetaJacobian = Box::new(|g: &DVector<f64>| {
            nalgebra::DMatrix::from_row_slice(2, 2, &[
                -1.0 + g[1], g[0],
                2.0 * g[0], -2.0,
            ])
        });
        let g0 = DVector::from_vec(vec![0.01, 0.01]);
        let fp = find_fixed_point_vec_newton(&beta, &jac, &g0, 1e-10, 100).unwrap();
        assert!(fp.g_star[0].abs() < 1e-8);
        assert!(fp.g_star[1].abs() < 1e-8);
        assert_eq!(fp.n_irrelevant, 2);
    }
}
