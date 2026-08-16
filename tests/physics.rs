//! Physics sanity tests from the project spec: photon-sphere orbit, weak-field
//! deflection, critical impact parameter, and conserved-quantity drift.

use schwarzschild_raytracer::integrator::{RayOutcome, TraceParams, rk4_step, trace};
use schwarzschild_raytracer::metric::{B_CRIT, State, null_residual};

/// Initial state for a ray at radius r0 with impact parameter b, moving
/// inbound, under the local normalization E = √f(r0).
fn inbound_ray(r0: f64, b: f64) -> (State, f64) {
    let f = 1.0 - 2.0 / r0;
    let e = f.sqrt();
    let l = b * e;
    let p_phi = l / (r0 * r0);
    let p_r = -(e * e - f * l * l / (r0 * r0)).sqrt();
    assert!(p_r.is_finite(), "b too large to reach r0");
    ([r0, 0.0, p_r, p_phi, 0.0], e)
}

/// In-plane Cartesian velocity d/dλ of (r cos φ, r sin φ).
fn plane_velocity(y: State) -> (f64, f64) {
    let [r, phi, p_r, p_phi, _t] = y;
    let (s, c) = phi.sin_cos();
    (p_r * c - r * p_phi * s, p_r * s + r * p_phi * c)
}

/// (a) A photon launched tangentially at r = 3 should stay near-circular over
/// a full orbit. The orbit is an unstable equilibrium, so we only bound the
/// deviation over 2π — long-term drift is physically correct.
#[test]
fn photon_sphere_orbit_stays_near_r3() {
    // p_r = 0, p_φ = 1/3 → E² = f·r²·p_φ² = 1/3 satisfies the null condition.
    let e: f64 = (1.0f64 / 3.0).sqrt();
    let mut y: State = [3.0, 0.0, 0.0, 1.0 / 3.0, 0.0];
    let mut max_dev: f64 = 0.0;
    while y[1] < 2.0 * std::f64::consts::PI {
        y = rk4_step(y, 1e-3, e);
        max_dev = max_dev.max((y[0] - 3.0).abs());
    }
    assert!(max_dev < 0.01, "max |r - 3| over one orbit = {max_dev:.2e}");
}

/// (b) Weak-field deflection: a ray with b = 100 should bend by ≈ 4/b.
#[test]
fn weak_deflection_matches_4_over_b() {
    let b = 100.0;
    let r0 = 2000.0;
    let (y0, e) = inbound_ray(r0, b);
    let params = TraceParams {
        e,
        step_scale: 0.02,
        max_steps: 200_000,
        r_far: r0,
        plane_az: 0.0,
        plane_bz: 0.0,
        disk: None,
    };
    let RayOutcome::Escaped { state } = trace(y0, &params) else {
        panic!("b = 100 ray must escape");
    };
    let (vx0, vy0) = plane_velocity(y0);
    let (vx1, vy1) = plane_velocity(state);
    let cos_ang = (vx0 * vx1 + vy0 * vy1)
        / ((vx0 * vx0 + vy0 * vy0).sqrt() * (vx1 * vx1 + vy1 * vy1).sqrt());
    let deflection = cos_ang.clamp(-1.0, 1.0).acos();
    let expected = 4.0 / b;
    assert!(
        (deflection - expected).abs() / expected < 0.05,
        "deflection {deflection:.5} vs expected {expected:.5}"
    );
}

/// (c) b_crit = 3√3 separates capture from escape: 0.1% above escapes,
/// 0.1% below falls through the horizon.
#[test]
fn critical_impact_parameter_separates_capture_and_escape() {
    let r0 = 1000.0;
    for (factor, expect_escape) in [(1.001, true), (0.999, false)] {
        let (y0, e) = inbound_ray(r0, B_CRIT * factor);
        let params = TraceParams {
            e,
            step_scale: 0.02,
            max_steps: 200_000,
            r_far: r0,
            plane_az: 0.0,
            plane_bz: 0.0,
            disk: None,
        };
        match (trace(y0, &params), expect_escape) {
            (RayOutcome::Escaped { .. }, true) | (RayOutcome::Horizon, false) => {}
            (_, true) => panic!("b = {factor}·b_crit should escape"),
            (_, false) => panic!("b = {factor}·b_crit should be captured"),
        }
    }
}

/// (d) Along a strongly-bent ray (b = 6, periapsis ≈ 4.2), the conserved
/// angular momentum L = r²·dφ/dλ and the null-condition residual (the
/// energy-conservation proxy — E enters the RHS analytically) stay tiny.
#[test]
fn conserved_quantities_drift_stays_small() {
    let (y0, e) = inbound_ray(30.0, 6.0);
    let e_sq = e * e;
    let l0 = y0[0] * y0[0] * y0[3];
    let mut y = y0;
    let mut max_l_drift: f64 = 0.0;
    let mut max_null: f64 = 0.0;
    let mut escaped = false;
    for _ in 0..200_000 {
        let h = 0.02 * (y[0] - 2.0);
        y = rk4_step(y, h, e);
        max_l_drift = max_l_drift.max((y[0] * y[0] * y[3] - l0).abs() / l0);
        max_null = max_null.max(null_residual(y, e_sq).abs() / e_sq);
        if y[0] > 100.0 && y[2] > 0.0 {
            escaped = true;
            break;
        }
    }
    assert!(escaped, "b = 6 ray from r = 30 must escape past r = 100");
    assert!(max_l_drift < 1e-6, "L drift = {max_l_drift:.2e}");
    assert!(max_null < 1e-6, "null residual = {max_null:.2e}");
}

/// (e) Coordinate time along a radial ray has the closed form
///   Δt = Δr + 2·ln((r₁ − 2)/(r₀ − 2))
/// (the tortoise coordinate r* = r + 2 ln(r/2 − 1); dt = dr/f on a radial
/// null ray). The integrated 5th state component must match it.
#[test]
fn radial_ray_time_matches_tortoise_coordinate() {
    let e = 1.0;
    // Outbound radial null ray: p_r = +E, p_φ = 0.
    let y0: State = [5.0, 0.0, e, 0.0, 0.0];
    let params = TraceParams {
        e,
        step_scale: 0.02,
        max_steps: 10_000,
        r_far: 500.0,
        plane_az: 0.0,
        plane_bz: 0.0,
        disk: None,
    };
    let RayOutcome::Escaped { state } = trace(y0, &params) else {
        panic!("outbound radial ray must escape");
    };
    let [r1, _, _, _, t] = state;
    let expected = (r1 - 5.0) + 2.0 * ((r1 - 2.0) / (5.0 - 2.0)).ln();
    assert!(
        ((t - expected) / expected).abs() < 1e-6,
        "Δt = {t:.8}, tortoise-coordinate value {expected:.8}"
    );
}

/// (f) Light from the disk's far side (imaged above the shadow) travels
/// longer than light from the near side (below the shadow): scanning the
/// center column, the topmost disk hit must carry a larger Δt than the
/// bottommost one.
#[test]
fn far_side_light_arrives_later_than_near_side() {
    use schwarzschild_raytracer::camera::Camera;
    let cam = Camera::new(30.0, 80.0, 0.0, 75.0, 100, 100);
    let mut hits: Vec<(u32, f64)> = vec![]; // (j, Δt at disk crossing)
    for j in 0..100 {
        let ray = cam.make_ray(cam.pixel_dir(50, j, 0.5, 0.5));
        let params = TraceParams {
            e: ray.e,
            step_scale: 0.02,
            max_steps: 60_000,
            r_far: 2000.0,
            plane_az: ray.e1.z,
            plane_bz: ray.e2.z,
            disk: Some((6.0, 20.0)),
        };
        if let RayOutcome::Disk { state } = trace(ray.y0, &params) {
            hits.push((j, state[4]));
        }
    }
    let (j_top, t_far) = hits[0];
    let (j_bot, t_near) = *hits.last().unwrap();
    assert!(j_top < 50 && j_bot > 50, "expected hits both above and below center");
    assert!(
        t_far > t_near,
        "far-side Δt = {t_far:.3} must exceed near-side Δt = {t_near:.3}"
    );
}
