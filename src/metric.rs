//! Schwarzschild metric quantities on the equatorial slice (M = 1, sin θ = 1).
//!
//! Every null geodesic in Schwarzschild lies in a plane through the origin,
//! so each ray is integrated in its own orbital plane with coordinates
//! (r, φ) and the equatorial-slice Christoffel symbols. The t-equation
//! decouples through the conserved energy E = (1 − 2/r) dt/dλ, leaving a
//! first-order system on the state [r, φ, dr/dλ, dφ/dλ].

/// Event horizon radius (coordinate singularity — terminate before it).
pub const R_HORIZON: f64 = 2.0;
/// Unstable circular photon orbit.
pub const R_PHOTON_SPHERE: f64 = 3.0;
/// Innermost stable circular orbit; inner edge of the accretion disk.
pub const R_ISCO: f64 = 6.0;
/// Critical impact parameter 3√3: inbound rays below this are captured.
/// This is the apparent radius of the black hole shadow.
pub const B_CRIT: f64 = 5.196152422706632;
/// Terminate integration this far above the horizon. Below ~1e-3 the two
/// 1/(r−2) terms in the RHS (which cancel only on-shell) start to lose
/// precision, and the emitted light is redshifted to black anyway.
pub const HORIZON_EPS: f64 = 1e-3;

/// Geodesic state in the ray's orbital plane: [r, φ, dr/dλ, dφ/dλ].
pub type State = [f64; 4];

/// Right-hand side of the null-geodesic system, with e_sq = E².
///
/// From the geodesic equation with the equatorial-slice Christoffel symbols
///   Γ^r_tt = (1−2/r)/r²,  Γ^r_rr = −1/(r(r−2)),  Γ^r_φφ = −(r−2),  Γ^φ_rφ = 1/r
/// and dt/dλ = E/(1−2/r):
#[inline]
pub fn geodesic_rhs(y: State, e_sq: f64) -> State {
    let [r, _phi, p_r, p_phi] = y;
    let f = 1.0 - 2.0 / r;
    [
        p_r,
        p_phi,
        -e_sq / (r * r * f) + p_r * p_r / (r * (r - 2.0)) + (r - 2.0) * p_phi * p_phi,
        -2.0 * p_r * p_phi / r,
    ]
}

/// Null-condition residual g_μν p^μ p^ν = −E²/f + (dr/dλ)²/f + r²(dφ/dλ)².
/// Zero on-shell; its drift measures integration error (E enters analytically,
/// so this doubles as the energy-conservation check).
pub fn null_residual(y: State, e_sq: f64) -> f64 {
    let [r, _phi, p_r, p_phi] = y;
    let f = 1.0 - 2.0 / r;
    (p_r * p_r - e_sq) / f + r * r * p_phi * p_phi
}
