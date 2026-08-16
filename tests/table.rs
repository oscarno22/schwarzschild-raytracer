//! δ-table verification: the analytic capture boundary, classification
//! parity against the full RK4 trace, and interpolation accuracy of the
//! tabulated escape directions.

use schwarzschild_raytracer::camera::Camera;
use schwarzschild_raytracer::deltatable::{DeltaTable, Lookup, TableParams};
use schwarzschild_raytracer::integrator::{RayOutcome, TraceParams, trace};
use schwarzschild_raytracer::metric::State;
use schwarzschild_raytracer::scene::{DISK_INNER, DISK_OUTER};

const R_CAM: f64 = 30.0;

fn small_params() -> TableParams {
    TableParams {
        n_rows: 1024,
        w_min: 1e-4,
        max_steps: 120_000,
        ..TableParams::default()
    }
}

/// Full RK4 trace of a screen-angle-ε ray, no disk.
fn trace_eps(eps: f64, r_far: f64) -> RayOutcome {
    let sqrt_f = (1.0 - 2.0 / R_CAM).sqrt();
    let y0: State = [R_CAM, 0.0, -eps.cos() * sqrt_f, eps.sin() / R_CAM, 0.0];
    trace(
        y0,
        &TraceParams {
            e: sqrt_f,
            step_scale: 0.02,
            max_steps: 300_000,
            r_far,
            plane_az: 0.0,
            plane_bz: 0.0,
            disk: None,
        },
    )
}

/// (1) The analytic ε_crit = asin(b_crit·√f/r_cam) must agree with a
/// capture/escape bisection through the integrator (as in tests/shadow.rs).
#[test]
fn eps_crit_matches_bisection() {
    let table = DeltaTable::build(R_CAM, small_params());
    let captured = |eps: f64| !matches!(trace_eps(eps, 2000.0), RayOutcome::Escaped { .. });
    let (mut lo, mut hi) = (0.05f64, 0.4f64);
    assert!(captured(lo) && !captured(hi));
    for _ in 0..45 {
        let mid = 0.5 * (lo + hi);
        if captured(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let eps_bisect = 0.5 * (lo + hi);
    assert!(
        (eps_bisect - table.eps_crit).abs() < 1e-6,
        "bisection {eps_bisect:.8} vs analytic {:.8}",
        table.eps_crit
    );
}

/// (2) Table classification (captured vs escaped) matches the full trace
/// for ~200 ε values; disagreement is permitted only within 2 local row
/// spacings of ε_crit.
#[test]
fn classification_matches_trace() {
    let table = DeltaTable::build(R_CAM, small_params());
    for k in 0..200 {
        let eps = 0.5 * table.eps_crit
            + (table.params.eps_max * 0.98 - 0.5 * table.eps_crit) * k as f64 / 199.0;
        let (j, _) = table.row_at(eps);
        let table_captured = eps <= table.eps_crit
            || table.rows[j].psi_inf.is_nan()
            || table.rows[j + 1].psi_inf.is_nan();
        let full_captured = !matches!(trace_eps(eps, 2000.0), RayOutcome::Escaped { .. });
        if table_captured != full_captured {
            let spacing = table.rows[j + 1].eps - table.rows[j].eps;
            assert!(
                (eps - table.eps_crit).abs() < 2.0 * spacing,
                "classification mismatch at ε = {eps:.6} (ε_crit {:.6}, spacing {spacing:.2e})",
                table.eps_crit
            );
        }
    }
}

/// (3) Lerped ψ_inf at row midpoints matches a fresh trace's asymptotic
/// in-plane angle to < 1e-3 rad away from the critical ring (the star PSF
/// floor is 6e-4 rad, so stars stay within ~1 PSF width).
#[test]
fn escape_direction_interpolation() {
    let table = DeltaTable::build(R_CAM, small_params());
    let mut checked = 0;
    for j in (0..table.rows.len() - 1).step_by(23) {
        let (r0, r1) = (&table.rows[j], &table.rows[j + 1]);
        let eps = 0.5 * (r0.eps + r1.eps);
        if eps < table.eps_crit + 5e-3 || r0.psi_inf.is_nan() || r1.psi_inf.is_nan() {
            continue;
        }
        let RayOutcome::Escaped { state } = trace_eps(eps, table.params.r_far) else {
            panic!("midpoint ε = {eps} must escape");
        };
        let psi_full = state[1] + (state[0] * state[3]).atan2(state[2]);
        let psi_lerp = 0.5 * (r0.psi_inf + r1.psi_inf);
        assert!(
            (psi_full - psi_lerp).abs() < 1e-3,
            "ψ mismatch at ε = {eps:.6}: full {psi_full:.6} vs lerp {psi_lerp:.6}"
        );
        checked += 1;
    }
    assert!(checked > 15, "too few midpoints checked ({checked})");
}

/// (5) Disk-hit parity through the camera path: for pixels both paths
/// classify as disk, interpolated hit radius and retarded time must match
/// the full trace to 0.05 M (spot azimuth error Ω·Δt ≈ 2e-3 rad —
/// invisible). Classification itself must agree away from the ring and the
/// annulus edges.
#[test]
fn disk_hits_match_full_trace() {
    let table = DeltaTable::build(R_CAM, small_params());
    let cam = Camera::new(R_CAM, 80.0, 0.0, 75.0, 60, 34);
    let sqrt_f = (1.0 - 2.0 / R_CAM).sqrt();
    let (mut compared, mut mismatched, mut disk_pixels) = (0, 0, 0);
    for j in 0..34 {
        for i in 0..60 {
            let ray = cam.make_ray(cam.pixel_dir(i, j, 0.5, 0.5));
            let eps = (ray.l / R_CAM).atan2(-ray.y0[2] / sqrt_f);
            let full = trace(
                ray.y0,
                &TraceParams {
                    e: ray.e,
                    step_scale: 0.02,
                    max_steps: 300_000,
                    r_far: table.params.r_far,
                    plane_az: ray.e1.z,
                    plane_bz: ray.e2.z,
                    disk: Some((DISK_INNER, DISK_OUTER)),
                },
            );
            let tab = table.sample(eps, ray.e1.z, ray.e2.z);
            match (&full, &tab) {
                (RayOutcome::Disk { state }, Lookup::Disk { r, t, .. }) => {
                    disk_pixels += 1;
                    // Annulus-edge hits can land on the other side of the
                    // 6/20 boundary under interpolation; skip the rim.
                    if state[0] - DISK_INNER > 0.2 && DISK_OUTER - state[0] > 0.2 {
                        assert!(
                            (state[0] - r).abs() < 0.05,
                            "hit radius {} vs table {r} at pixel ({i},{j})",
                            state[0]
                        );
                        assert!(
                            (state[4] - t).abs() < 0.05,
                            "hit time {} vs table {t} at pixel ({i},{j})",
                            state[4]
                        );
                        compared += 1;
                    }
                }
                (RayOutcome::Disk { .. }, _) | (_, Lookup::Disk { .. }) => {
                    disk_pixels += 1;
                    // Classification mismatch: tolerated only near the ring
                    // or the annulus rim (checked in aggregate below).
                    if (eps - table.eps_crit).abs() > 5e-3 {
                        mismatched += 1;
                    }
                }
                _ => {}
            }
        }
    }
    assert!(compared > 50, "too few disk pixels compared ({compared})");
    assert!(
        (mismatched as f64) < 0.05 * disk_pixels as f64,
        "{mismatched}/{disk_pixels} disk classification mismatches"
    );
}
