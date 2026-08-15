//! RK4 integration of null geodesics with an adaptive step h ∝ (r − 2).
//!
//! The step shrinks geometrically approaching the horizon (never overshooting
//! the coordinate singularity), is small near the photon sphere where the
//! trajectory curvature peaks, and grows ∝ r in the weak field — so a ray can
//! cheaply be followed out to very large radii for accurate escape directions.

use crate::metric::{HORIZON_EPS, R_HORIZON, State, geodesic_rhs};

pub struct TraceParams {
    /// Conserved energy E of the ray (sets dt/dλ; enters the RHS as E²).
    pub e: f64,
    /// Step size factor: h = step_scale · (r − 2).
    pub step_scale: f64,
    /// Safety cap on steps; only rays within ~1e-6 of the critical impact
    /// parameter (orbiting the photon sphere indefinitely) hit this.
    pub max_steps: u32,
    /// Outward-moving rays beyond this radius terminate as escaped.
    pub r_far: f64,
}

pub enum RayOutcome {
    /// Fell through the horizon (or the integrator stepped inside it).
    Horizon,
    /// Reached r_far moving outward; state gives the asymptotic direction.
    Escaped { state: State },
    /// Step cap hit — treated as captured (photon-sphere limbo).
    MaxSteps,
}

#[inline]
fn axpy(y: State, k: State, h: f64) -> State {
    [
        y[0] + h * k[0],
        y[1] + h * k[1],
        y[2] + h * k[2],
        y[3] + h * k[3],
    ]
}

/// One classic RK4 step of size h.
#[inline]
pub fn rk4_step(y: State, h: f64, e_sq: f64) -> State {
    let k1 = geodesic_rhs(y, e_sq);
    let k2 = geodesic_rhs(axpy(y, k1, h * 0.5), e_sq);
    let k3 = geodesic_rhs(axpy(y, k2, h * 0.5), e_sq);
    let k4 = geodesic_rhs(axpy(y, k3, h), e_sq);
    [
        y[0] + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        y[1] + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        y[2] + h / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        y[3] + h / 6.0 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]),
    ]
}

/// Integrate a ray from y0 until it falls in, escapes, or exhausts max_steps.
pub fn trace(y0: State, p: &TraceParams) -> RayOutcome {
    let e_sq = p.e * p.e;
    let mut y = y0;
    for _ in 0..p.max_steps {
        let h = p.step_scale * (y[0] - R_HORIZON);
        let y_new = rk4_step(y, h, e_sq);
        // Negated comparison so NaN (from a step landing on the singularity)
        // also terminates as Horizon instead of looping forever.
        if !(y_new[0] > R_HORIZON + HORIZON_EPS) {
            return RayOutcome::Horizon;
        }
        if y_new[0] > p.r_far && y_new[2] > 0.0 {
            return RayOutcome::Escaped { state: y_new };
        }
        y = y_new;
    }
    RayOutcome::MaxSteps
}
