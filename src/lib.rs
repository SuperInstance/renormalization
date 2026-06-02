#![deny(unsafe_code)]
//! # renormalization
//!
//! Renormalization group theory — the mathematical framework for understanding
//! scale-dependent behavior in physical and complex systems.
//!
//! ## Core Concepts
//!
//! - **RG Flow Equations**: Beta functions β(g) = dg/d(ln μ) governing coupling evolution
//! - **Fixed Points**: Solutions to β(g*)=0 with stability classification
//! - **Critical Exponents**: Eigenvalues of β' at fixed points via linearized RG
//! - **Wilsonian RG**: Momentum shell integration and effective actions
//! - **Universality Classes**: Systems sharing the same fixed point
//! - **Epsilon Expansion**: Systematic d = 4-ε analysis
//! - **Scaling Relations**: ξ ∼ |t|^{-ν}, hyperscaling, and more
//! - **Real-Space RG**: Block spin transformations, Migdal-Kadanoff

pub mod beta;
pub mod fixed_point;
pub mod critical_exponents;
pub mod wilsonian;
pub mod universality;
pub mod epsilon_expansion;
pub mod scaling;
pub mod real_space;

pub use beta::*;
pub use fixed_point::*;
pub use critical_exponents::*;
pub use wilsonian::*;
pub use universality::*;
pub use epsilon_expansion::*;
pub use scaling::*;
pub use real_space::*;
