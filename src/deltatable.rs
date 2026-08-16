//! Precomputed δ-ray table for real-time rendering.
//!
//! At fixed camera radius the entire geodesic trace is a function of the
//! single screen angle ε (measured off the toward-black-hole direction):
//! every ray lies in its own orbital plane with the same in-plane shape.
//! We tabulate, for a geometrically refined grid of ε values, the full
//! trajectory polyline (φ_p, r, t) plus the asymptotic escape direction.
//! A pixel then needs a couple of binary searches instead of ~500 RK4
//! steps, which is what makes the interactive viewer real-time on CPU.
//!
//! Conventions: the polyline's φ_p is strictly increasing (the camera
//! construction gives every ray L ≥ 0), while r is *not* monotonic
//! (periapsis) — all interpolation is therefore done in φ_p, never in r.

use rayon::prelude::*;

use crate::integrator::rk4_step;
use crate::metric::{B_CRIT, HORIZON_EPS, R_HORIZON, State};

pub struct TableParams {
    /// Number of ε rows.
    pub n_rows: usize,
    /// Offset of the innermost row above ε_crit, radians.
    pub w_min: f64,
    /// Angular coverage; must exceed the half-diagonal of the widest FOV
    /// the renderer will be asked for (zoom must never leave the table).
    pub eps_max: f64,
    /// RK4 step factor, as in the offline renderer.
    pub step_scale: f64,
    /// Step cap; only near-critical rows (winding many times) approach it.
    pub max_steps: u32,
    /// Polyline decimation: store a point when Δφ since the last stored
    /// point reaches this…
    pub dphi_max: f64,
    /// …or when |Δr| does (densifies periapsis and the steep radial legs).
    pub dr_max: f64,
    /// Escape radius. MUST equal `render::R_FAR`: the asymptotic direction
    /// keeps bending measurably beyond r = 50, and a mismatch shifts the
    /// starfield relative to the offline renderer.
    pub r_far: f64,
}

impl Default for TableParams {
    fn default() -> Self {
        // Coverage: half-diagonal of a 100° horizontal FOV at 16:9 under
        // the gnomonic projection, +2% margin. The viewer clamps its zoom
        // inside this, so FOV changes never rebuild the table.
        let fov_max = 100.0f64.to_radians();
        let eps_max = ((fov_max / 2.0).tan() * (1.0 + (9.0f64 / 16.0).powi(2)).sqrt())
            .atan()
            * 1.02;
        Self {
            n_rows: 4096,
            w_min: 1e-6,
            eps_max,
            step_scale: 0.02,
            max_steps: 300_000,
            dphi_max: 0.01,
            dr_max: 0.1,
            r_far: crate::render::R_FAR,
        }
    }
}

/// One tabulated ray: the trajectory polyline in the ray's orbital plane.
pub struct DeltaRow {
    pub eps: f64,
    /// Strictly increasing, phi[0] = 0 at the camera.
    pub phi: Vec<f64>,
    pub r: Vec<f64>,
    /// Coordinate time along the (backward) trace, t[0] = 0.
    pub t: Vec<f64>,
    /// Unwrapped in-plane asymptotic angle: the escape direction is
    /// cos ψ·e1 + sin ψ·e2. Stored as φ_end + atan2(r·p_φ, p_r), which is
    /// continuous across rows (never a wrapped atan2 — lerping across a 2π
    /// seam is the classic bug). NaN marks a row that exhausted max_steps
    /// (photon-sphere limbo) and is treated as captured.
    pub psi_inf: f64,
}

pub struct DeltaTable {
    pub r_cam: f64,
    pub sqrt_f: f64,
    /// Capture boundary: an inbound ray is captured iff b < 3√3, and
    /// b = r_cam·sin ε/√f, so ε_crit = asin(B_CRIT·√f/r_cam) — analytic;
    /// validated against trace-bisection in tests/table.rs.
    pub eps_crit: f64,
    pub params: TableParams,
    /// Ascending ε, all supercritical: ε_j = ε_crit + w_min·ρ^j.
    pub rows: Vec<DeltaRow>,
}

impl DeltaTable {
    pub fn build(r_cam: f64, params: TableParams) -> Self {
        assert!(r_cam >= 6.0, "table requires r_cam ≥ 6 (shadow < 90°)");
        let sqrt_f = (1.0 - 2.0 / r_cam).sqrt();
        let eps_crit = (B_CRIT * sqrt_f / r_cam).asin();
        // Geometric refinement: row spacing shrinks toward ε_crit where
        // deflection diverges ~ −ln(ε − ε_crit), and grows to ~0.002 rad
        // (still sub-pixel at typical FOVs) at the far end.
        let span = params.eps_max - eps_crit;
        let rho = (span / params.w_min).powf(1.0 / (params.n_rows - 1) as f64);
        let rows: Vec<DeltaRow> = (0..params.n_rows)
            .into_par_iter()
            .map(|j| {
                let eps = eps_crit + params.w_min * rho.powi(j as i32);
                trace_row(r_cam, sqrt_f, eps, &params)
            })
            .collect();
        Self {
            r_cam,
            sqrt_f,
            eps_crit,
            params,
            rows,
        }
    }

    /// Bracketing row index j (rows[j].eps ≤ eps, rows[j+1].eps ≥ eps when
    /// possible) and the linear blend weight toward row j+1.
    pub fn row_at(&self, eps: f64) -> (usize, f64) {
        let n = self.rows.len();
        let hi = self.rows.partition_point(|row| row.eps < eps);
        let j = hi.clamp(1, n - 1) - 1;
        let (e0, e1) = (self.rows[j].eps, self.rows[j + 1].eps);
        let w = ((eps - e0) / (e1 - e0)).clamp(0.0, 1.0);
        (j, w)
    }
}

/// Integrate one table ray. Same numerics as `integrator::trace` but with
/// no disk detection (disk crossings are resolved per pixel at lookup time)
/// and with polyline recording.
fn trace_row(r_cam: f64, sqrt_f: f64, eps: f64, p: &TableParams) -> DeltaRow {
    let e = sqrt_f;
    let mut y: State = [r_cam, 0.0, -eps.cos() * sqrt_f, eps.sin() / r_cam, 0.0];
    let mut phi = vec![y[1]];
    let mut r = vec![y[0]];
    let mut t = vec![y[4]];
    let mut psi_inf = f64::NAN;
    for _ in 0..p.max_steps {
        let h = p.step_scale * (y[0] - R_HORIZON);
        let y_new = rk4_step(y, h, e);
        if !(y_new[0] > R_HORIZON + HORIZON_EPS) {
            // Supercritical rows never reach the horizon in exact math;
            // treat a numerical graze as limbo (rendered captured).
            break;
        }
        y = y_new;
        if y[1] - phi.last().unwrap() >= p.dphi_max
            || (y[0] - r.last().unwrap()).abs() >= p.dr_max
        {
            phi.push(y[1]);
            r.push(y[0]);
            t.push(y[4]);
        }
        if y[0] > p.r_far && y[2] > 0.0 {
            psi_inf = y[1] + (y[0] * y[3]).atan2(y[2]);
            break;
        }
    }
    if *phi.last().unwrap() < y[1] {
        phi.push(y[1]);
        r.push(y[0]);
        t.push(y[4]);
    }
    DeltaRow {
        eps,
        phi,
        r,
        t,
        psi_inf,
    }
}
