# renormalization

**Renormalization in Rust. From lattice models to universal behavior.**

[![Rust](https://img.shields.io/badge/rust-2021-orange.svg)](https://www.rust-lang.org/)
[![License: MIT OR Apache-2.0](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE)

79 tests · 2,100+ lines of Rust · 8 modules

---

## What This Does

The renormalization group explains why wildly different physical systems behave identically near critical points. Water boiling, magnets losing magnetization, binary alloys ordering — they all follow the same mathematical law. The RG tells you which details matter and which don't — coarse-grain the microscopic description, and most details wash out, leaving only a few "relevant" parameters that determine universal behavior.

You get:
- **Beta functions** — the RG flow equations β(g) = dg/d(ln μ) with Euler and RK4 integrators
- **Fixed point analysis** — Newton's method and scanning for β(g*)=0, stability classification
- **Critical exponents** — ν, α, β, γ, δ, η with scaling relation verification (Rushbrooke, Widom, Fisher, Josephson)
- **Wilsonian RG** — momentum shell integration, effective actions, one-loop flow equations
- **Epsilon expansion** — systematic d=4−ε analysis for O(n) models to O(ε²)
- **Universality classes** — Ising 2D/3D, XY, Heisenberg, Potts, mean-field with exponent registries
- **Scaling relations** — correlation length, order parameter, susceptibility, finite-size scaling
- **Real-space RG** — decimation, Migdal-Kadanoff approximation, critical coupling search

---

## Key Idea

The renormalization group is a mathematical microscope that zooms out from microscopic details:

| RG Concept | Physical Meaning |
|---|---|
| Beta function β(g) | How coupling constants change with scale |
| Fixed point β(g*)=0 | Scale-invariant theory (phase transition) |
| Relevant coupling | Grows under RG — controls the physics |
| Irrelevant coupling | Shrinks under RG — microscopic detail that doesn't matter |
| Universality class | All systems flowing to the same fixed point |
| Critical exponents | Eigenvalues of the linearized RG at the fixed point |
| Epsilon expansion | Perturbative expansion in d = 4−ε |
| Wilsonian RG | Integrate out high-momentum modes shell by shell |

The deep insight: **microscopic details don't matter**. Systems with different Hamiltonians but the same symmetry and dimensionality flow to the same fixed point and share identical critical exponents.

---

## Install

```toml
[dependencies]
renormalization = "0.1"
```

Requires Rust 2021 edition. Dependencies: `serde`, `nalgebra`.

---

## Quick Start

```rust
use renormalization::*;

// 1. RG flow: integrate a beta function
let beta: BetaFn = Box::new(|g: f64| -g + g * g);
let config = FlowConfig::default();
let traj = flow(&beta, 0.5, &config);

// 2. Find fixed points
let beta_prime: BetaFn = Box::new(|g: f64| 1.0 - 2.0 * g);
let fp0 = find_fixed_point_newton(&beta, &beta_prime, 0.1, 1e-10, 100).unwrap();
assert!(fp0.g_star.abs() < 1e-8);    // g* = 0
assert_eq!(fp0.stability, Stability::Unstable);

let fp1 = find_fixed_point_newton(&beta, &beta_prime, 1.5, 1e-10, 100).unwrap();
assert!((fp1.g_star - 1.0).abs() < 1e-8); // g* = 1
assert_eq!(fp1.stability, Stability::Stable);

// 3. φ⁴ theory beta function
let phi4_beta = phi4_beta_one_loop(1.0); // ε=1 (d=3)
let g_star = 16.0 * PI * PI / 3.0;       // Wilson-Fisher fixed point
assert!(phi4_beta(g_star).abs() < 1e-10);

// 4. Critical exponents for Ising 3D
let ising = CriticalExponents::ising_3d();
let mf = CriticalExponents::mean_field();
let ising2d = CriticalExponents::ising_2d();

// Verify scaling relations
assert!(ising.verify_rushbrooke().abs() < 0.01); // α + 2β + γ = 2
assert!(ising.verify_widom().abs() < 0.01);      // γ = β(δ-1)
assert!(ising.verify_fisher().abs() < 0.01);     // γ = ν(2-η)

// 5. Epsilon expansion for O(n) model
let exp = EpsilonExpansion::compute(1, 1.0); // Ising (n=1) in d=3 (ε=1)
let ce = exp.critical_exponents();
assert!(ce.nu > 0.5); // ν > mean-field value

// 6. Universality classes
let ising_cls = UniversalityClass::ising_3d();
let heisenberg_cls = UniversalityClass::heisenberg_3d();
assert!(!same_universality_class(&ising_cls, &heisenberg_cls));

let registry = UniversalityRegistry::standard();
let found = registry.find(3, &SymmetryGroup::Z2); // finds 3D Ising
```

---

## API Reference

### Beta Functions and Flow Integration

```rust
// One-loop φ⁴ beta function: β(g) = -εg + 3g²/(16π²)
let beta = phi4_beta_one_loop(1.0);

// Gaussian beta: β(g) = -εg
let gauss = gaussian_beta(1.0);

// Custom beta function
let beta: BetaFn = Box::new(|g: f64| -g + g*g);

// Integrate the flow
let config = FlowConfig {
    steps: 1000,
    d_ln_mu: 0.01,
    method: IntegrationMethod::RK4,
};
let traj = flow(&beta, 0.5, &config);

// Classify coupling relevance at a fixed point
let relevance = coupling_relevance(&beta_prime, g_star);
// Relevant (flows away), Irrelevant (flows toward), Marginal (need higher order)
```

Multi-dimensional flows via `flow_vec_euler` with `BetaFnVec` and `BetaJacobian`.

### Fixed Point Analysis

```rust
// Newton's method to find β(g*)=0
let fp = find_fixed_point_newton(&beta, &beta_prime, 1.0, 1e-10, 100);

// Scan a range for all fixed points
let fps = find_fixed_points_scan(&beta, &beta_prime, -1.0, 2.0, 200, 1e-10);

// Basin of attraction for a stable fixed point
let basin = basin_of_attraction(&beta, &fp, 5.0, &config);

// Multi-dimensional fixed points
let fp_vec = find_fixed_point_vec_newton(&beta_vec, &jacobian, &g0, 1e-10, 100);
// Returns eigenvalues, n_relevant, n_irrelevant, n_marginal
```

### Critical Exponents

```rust
// Pre-built exponent sets
let mf = CriticalExponents::mean_field();     // ν=0.5, α=0, β=0.5, γ=1, δ=3, η=0
let ising2d = CriticalExponents::ising_2d();   // ν=1, α=0 (log), β=1/8, γ=7/4, η=1/4
let ising3d = CriticalExponents::ising_3d();   // ν≈0.630, β≈0.326, γ≈1.237

// Verify scaling relations
ising3d.verify_rushbrooke();  // α + 2β + γ = 2
ising3d.verify_widom();       // γ = β(δ-1)
ising3d.verify_fisher();      // γ = ν(2-η)
ising3d.verify_hyperscaling(); // 2-α = dν

// Compute from eigenvalue at fixed point
let ce = exponents_from_eigenvalue(-2.0, 3.0); // ν=0.5 in d=3

// Compute from ν and η
let ce = exponents_from_nu_eta(0.63, 0.036, 3.0);
```

### Wilsonian RG

```rust
// φ⁴ effective action
let action = EffectiveAction::phi4(0.1, 0.1, 3, 1.0); // (r, u, d, Λ)

// One Wilsonian step: integrate out shell [Λ-dΛ, Λ]
let stepped = action.wilsonian_step(0.01);

// Full flow
let params = WilsonianParams { n_shells: 100, d_lambda: 0.01, ..Default::default() };
let flow = action.wilsonian_flow(&params);

// Exact renormalization group (Wetterich equation)
let erg = ErgFlow::new(vec![0.5, 0.3], 1.0, 10.0);
let stepped = erg.step(0.1, &[-0.5, -0.3]);
```

### Epsilon Expansion

```rust
// Wilson-Fisher fixed point to O(ε) and O(ε²)
let g1 = EpsilonExpansion::wilson_fisher_fixed_point_1loop(1, 1.0); // Ising
let g2 = EpsilonExpansion::wilson_fisher_fixed_point_2loop(1, 1.0);

// Anomalous dimension η to O(ε²)
let eta = EpsilonExpansion::eta_o_eps2(1, 1.0); // η ~ (n+2)ε²/(2(n+8)²)

// Correlation length exponent ν to O(ε) and O(ε²)
let nu1 = EpsilonExpansion::nu_o_eps1(1, 1.0); // ν = 1/2 + (n+2)ε/(4(n+8))
let nu2 = EpsilonExpansion::nu_o_eps2(1, 1.0);

// Full expansion
let exp = EpsilonExpansion::compute(1, 1.0); // n=1, ε=1
let ce = exp.critical_exponents();

// Scan ε values
let results = epsilon_scan(1, &[0.1, 0.2, 0.5, 1.0]);
```

### Universality Classes

```rust
// Pre-built classes
let classes = [
    UniversalityClass::ising_2d(),     // Z₂, d=2, exact exponents
    UniversalityClass::ising_3d(),     // Z₂, d=3, numerical exponents
    UniversalityClass::xy_3d(),        // O(2), d=3
    UniversalityClass::heisenberg_3d(),// O(3), d=3
    UniversalityClass::mean_field(),   // above d_c=4
    UniversalityClass::potts_3_2d(),   // 3-state Potts in d=2
];

// Check if two systems share a universality class
same_universality_class(&ising, &heisenberg); // false (different symmetry)

// Compare critical exponents quantitatively
let max_diff = compare_exponents(&ising, &heisenberg);

// Registry of known classes
let reg = UniversalityRegistry::standard();
let found = reg.find(3, &SymmetryGroup::Z2); // → 3D Ising
```

### Scaling Relations

```rust
// Physical observables near criticality
let xi = correlation_length(0.63, 0.01, 1.0);    // ξ ∼ |t|^{-ν}
let m = order_parameter(0.326, -0.01, 1.0);       // M ∼ (-t)^β
let chi = susceptibility(1.237, 0.01, 1.0);        // χ ∼ |t|^{-γ}
let c = specific_heat(0.110, 0.01, 1.0);           // C ∼ |t|^{-α}

// Finite-size scaling
let scaled = finite_size_scaling(-0.326, 0.63, 100.0, 1.0);

// Scaling collapse for data analysis
let (x, y) = scaling_collapse_transform(t, observable, L, nu, beta);

// Corrections to scaling (ω ≈ 0.832 for 3D Ising)
let corr = CorrectionsToScaling::ising_3d();
let corrected = corr.apply_correction(base, t, a1, a2);
```

### Real-Space RG

```rust
// 1D Ising decimation: K' = atanh(tanh(K)²)
let k_new = decimation_1d_ising(1.0);

// Migdal-Kadanoff approximation
let k_next = migdal_kadanoff_step(0.5, 2.0, 2); // (K, b=2, d=2)

// Full MK flow
let config = RealSpaceConfig { b: 2.0, dimension: 2, n_iterations: 50 };
let result = migdal_kadanoff_flow(&config, 0.5);

// Find critical coupling by bisection
let kc = find_critical_coupling_mk(&config, 0.001);

// Extract ν from the fixed point
let nu = nu_from_mk(kc, 2.0, 2);

// Majority rule block spin
let block_spin = majority_rule_block(&[1, 1, 1, -1, -1]); // → 1
```

---

## How It Works

The crate implements the RG as a layered analysis pipeline:

```
Layer 1: Flow Equations       beta (Euler/RK4 integrators, coupling relevance)
              │
              ▼
Layer 2: Fixed Points         fixed_point (Newton's method, stability, basins)
              │
              ▼
Layer 3: Critical Exponents   critical_exponents (ν,α,β,γ,δ,η + scaling relations)
              │
              ├──▶ 4a: Epsilon Expansion   epsilon_expansion (d=4-ε perturbation theory)
              │
              ├──▶ 4b: Wilsonian RG        wilsonian (momentum shell, effective actions)
              │
              └──▶ 4c: Real-Space RG       real_space (decimation, Migdal-Kadanoff)
              │
              ▼
Layer 5: Universality         universality (Ising, XY, Heisenberg, Potts, mean-field)
              │
              ▼
Layer 6: Scaling              scaling (ξ, M, χ, C, finite-size, corrections)
```

---

## The Math

### Beta Functions and RG Flow

The **beta function** β(g) = dg/d(ln μ) describes how a coupling constant g evolves with energy scale μ. In φ⁴ theory at one loop:

**β(g) = −εg + (n+8)g²/(48π²)**

where ε = 4−d and n is the number of field components. Fixed points β(g*)=0 correspond to scale-invariant theories.

### Critical Exponents

Near a fixed point, the linearized RG gives scaling dimensions yᵢ (eigenvalues of the stability matrix). The critical exponents are:

- **ν = 1/yₜ**: correlation length ξ ∼ |t|^{−ν}
- **α = 2−dν**: specific heat C ∼ |t|^{−α} (Josephson hyperscaling)
- **β**: order parameter M ∼ (−t)^β
- **γ = ν(2−η)**: susceptibility χ ∼ |t|^{−γ} (Fisher relation)
- **δ**: critical isotherm M ∼ h^{1/δ}
- **η**: correlation function G(r) ∼ r^{−(d−2+η)}

These satisfy **Rushbrooke** (α + 2β + γ = 2) and **Widom** (γ = β(δ−1)) as exact inequalities (thermodynamic) or equalities (scaling hypothesis).

### Epsilon Expansion

Wilson & Fisher's systematic expansion in ε = 4−d:

- Wilson-Fisher fixed point: g* = 48π²ε/(n+8) + O(ε²)
- ν = ½ + (n+2)ε/(4(n+8)) + O(ε²)
- η = (n+2)ε²/(2(n+8)²) + O(ε³)

At ε=1 (d=3), these give surprisingly accurate results despite being asymptotic series.

### Universality

Systems sharing the same spatial dimension, order parameter symmetry, and interaction range flow to the same fixed point and share identical critical exponents. Water-liquid transition and uniaxial magnets are both in the 3D Ising universality class (Z₂ symmetry, d=3, short-range interactions).

### Wilsonian RG

Integrate out high-momentum modes shell by shell. The effective action evolves:

**dr/dl = 2r + (n+2)u·K_d·Λ^d/(r+Λ²)**
**du/dl = (4−d)u − (n+8)u²·K_d·Λ^d/(r+Λ²)²**

The first equation says the mass grows (relevant), the second shows the coupling is marginal at d=4 and relevant for d<4.

### Migdal-Kadanoff

A real-space approximation: move bonds to create a 1D chain, then decimate. For the Ising model with scale factor b in d dimensions:

1. Bond moving: K → b^{d−1}·K
2. Decimation chain: K' = atanh(tanh(K)²) applied (b−1) times

This gives approximate critical exponents — exact on hierarchical lattices.

---

## Module Overview

| Module | Tests | Key Types | Purpose |
|--------|-------|-----------|---------|
| `beta` | 9 | `BetaFn`, `FlowConfig`, `FlowTrajectory` | RG flow equations and integrators |
| `fixed_point` | 6 | `FixedPoint`, `Stability`, `FixedPointVec` | Fixed point finding and classification |
| `critical_exponents` | 9 | `CriticalExponents` | Exponents with scaling relation verification |
| `wilsonian` | 9 | `EffectiveAction`, `WilsonianParams`, `ErgFlow` | Momentum shell RG |
| `universality` | 10 | `UniversalityClass`, `UniversalityRegistry`, `SymmetryGroup` | Universality class catalog |
| `epsilon_expansion` | 10 | `EpsilonExpansion` | d=4−ε perturbation theory |
| `scaling` | 11 | `ScalingVerification`, `CorrectionsToScaling` | Scaling relations and finite-size effects |
| `real_space` | 10 | `RealSpaceConfig`, `RealSpaceResult` | Decimation and Migdal-Kadanoff |

---

## References

- **Renormalization Group:** Wilson & Kogut, *The Renormalization Group and the ε Expansion* (1974)
- **Phase Transitions:** Goldenfeld, *Lectures on Phase Transitions and the Renormalization Group* (1992)
- **Critical Phenomena:** Cardy, *Scaling and Renormalization in Statistical Physics* (1996)
- **Epsilon Expansion:** Zinn-Justin, *Quantum Field Theory and Critical Phenomena* (2002)
- **Real-Space RG:** Nelson & Fisher, *Soluble Renormalization Groups and Scaling Theory* (1975)

---

## License

MIT OR Apache-2.0
