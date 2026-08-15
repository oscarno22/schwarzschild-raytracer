//! RK4 integration of null geodesics with an adaptive step h ∝ (r − 2).
//!
//! The step shrinks geometrically approaching the horizon (never overshooting
//! the coordinate singularity), is small near the photon sphere where the
//! trajectory curvature peaks, and grows ∝ r in the weak field — so a ray can
//! cheaply be followed out to very large radii for accurate escape directions.

use crate::metric::{HORIZON_EPS, R_HORIZON, State, geodesic_rhs};

pub struct TraceParams {
    /// Conserved energy E of the ray (sets dt/dλ = E/(1 − 2/r)).
    pub e: f64,
    /// Step size factor: h = step_scale · (r − 2).
    pub step_scale: f64,
    /// Safety cap on steps; only rays within ~1e-6 of the critical impact
    /// parameter (orbiting the photon sphere indefinitely) hit this.
    pub max_steps: u32,
    /// Outward-moving rays beyond this radius terminate as escaped.
    pub r_far: f64,
    /// z-components of the ray-plane basis (e1.z, e2.z): the global height
    /// of the ray is z = r(a·cos φ + b·sin φ) — the only piece of 3D
    /// geometry the hot loop needs.
    pub plane_az: f64,
    pub plane_bz: f64,
    /// Accretion disk annulus (inner, outer) in the global equatorial plane,
    /// or None to trace through it (shadow-only renders).
    pub disk: Option<(f64, f64)>,
}

pub enum RayOutcome {
    /// Fell through the horizon (or the integrator stepped inside it).
    Horizon,
    /// Crossed the equatorial plane inside the disk annulus; state is
    /// interpolated to the crossing (state[0] is the hit radius).
    Disk { state: State },
    /// Reached r_far moving outward; state gives the asymptotic direction.
    Escaped { state: State },
    /// Step cap hit — treated as captured (photon-sphere limbo).
    MaxSteps,
}

#[inline]
fn axpy(y: State, k: State, h: f64) -> State {
    std::array::from_fn(|i| y[i] + h * k[i])
}

/// One classic RK4 step of size h.
#[inline]
pub fn rk4_step(y: State, h: f64, e: f64) -> State {
    let k1 = geodesic_rhs(y, e);
    let k2 = geodesic_rhs(axpy(y, k1, h * 0.5), e);
    let k3 = geodesic_rhs(axpy(y, k2, h * 0.5), e);
    let k4 = geodesic_rhs(axpy(y, k3, h), e);
    std::array::from_fn(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]))
}

/// Global height above the equatorial plane, from in-plane coordinates.
#[inline]
fn global_z(y: State, p: &TraceParams) -> f64 {
    let (s, c) = y[1].sin_cos();
    y[0] * (p.plane_az * c + p.plane_bz * s)
}

/// Integrate a ray from y0 until it falls in, hits the disk, escapes, or
/// exhausts max_steps.
pub fn trace(y0: State, p: &TraceParams) -> RayOutcome {
    let mut y = y0;
    let mut z = global_z(y0, p);
    for _ in 0..p.max_steps {
        let h = p.step_scale * (y[0] - R_HORIZON);
        let y_new = rk4_step(y, h, p.e);
        // Negated comparison so NaN (from a step landing on the singularity)
        // also terminates as Horizon instead of looping forever.
        if !(y_new[0] > R_HORIZON + HORIZON_EPS) {
            return RayOutcome::Horizon;
        }
        let z_new = global_z(y_new, p);
        if let Some((inner, outer)) = p.disk
            && z * z_new < 0.0
        {
            // Crossed the equatorial plane inside this step; refine by
            // linear interpolation in λ (steps in the disk region are ~0.1,
            // far below a pixel's footprint).
            let s = z / (z - z_new);
            let hit: State = std::array::from_fn(|k| y[k] + s * (y_new[k] - y[k]));
            // A crossing outside the annulus (the inner gap or beyond the
            // rim) continues: later crossings image the disk's far side.
            if hit[0] >= inner && hit[0] <= outer {
                return RayOutcome::Disk { state: hit };
            }
        }
        if y_new[0] > p.r_far && y_new[2] > 0.0 {
            return RayOutcome::Escaped { state: y_new };
        }
        y = y_new;
        z = z_new;
    }
    RayOutcome::MaxSteps
}
