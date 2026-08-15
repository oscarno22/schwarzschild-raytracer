//! Shadow-size validation: the sharp edge of the black hole's silhouette must
//! sit at the critical impact parameter b_crit = 3√3, both in impact-parameter
//! space (bisection through the camera code) and in rendered pixels.

use schwarzschild_raytracer::camera::Camera;
use schwarzschild_raytracer::integrator::{RayOutcome, TraceParams, trace};
use schwarzschild_raytracer::metric::B_CRIT;
use schwarzschild_raytracer::render::{R_FAR, render};
use schwarzschild_raytracer::{Config, vec3::Vec3};

/// (e) Bisect the screen angle between captured and escaped rays; the
/// recovered impact parameter must converge to b_crit. This exercises the
/// full camera initial-condition path, so it catches tetrad-normalization
/// bugs that pure-integrator tests cannot.
#[test]
fn shadow_edge_converges_to_b_crit() {
    let r_cam = 30.0;
    let cam = Camera::new(r_cam, 90.0, 75.0, 100, 100);
    let toward_bh = -(cam.pos.normalize());
    let perp = Vec3::new(0.0, 0.0, 1.0); // perpendicular to the radial line
    let captured = |eps: f64| {
        let d = toward_bh * eps.cos() + perp * eps.sin();
        let ray = cam.make_ray(d);
        let params = TraceParams {
            e: ray.e,
            step_scale: 0.02,
            max_steps: 200_000,
            r_far: R_FAR,
        };
        // MaxSteps means photon-sphere limbo at the exact boundary: captured.
        !matches!(trace(ray.y0, &params), RayOutcome::Escaped { .. })
    };
    let (mut lo, mut hi) = (0.05f64, 0.4f64); // radians off the BH direction
    assert!(captured(lo) && !captured(hi), "bisection bracket invalid");
    for _ in 0..45 {
        let mid = 0.5 * (lo + hi);
        if captured(mid) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let eps_edge = 0.5 * (lo + hi);
    let b_edge = r_cam * eps_edge.sin() / (1.0 - 2.0 / r_cam).sqrt();
    assert!(
        (b_edge - B_CRIT).abs() < 1e-3,
        "shadow edge at b = {b_edge:.6}, expected {B_CRIT:.6}"
    );
}

/// The rendered silhouette's pixel radius must match the analytic prediction
///   α = asin(b_crit·√(1 − 2/r_cam)/r_cam)   (angular radius)
///   px = (W/2)·tan α / tan(fov/2)           (gnomonic screen projection)
/// At r_cam = 30, fov = 75°, W = 160 that is 17.7 px — note the naive
/// angle-proportional mapping would predict 20.9 px and must fail.
#[test]
fn rendered_shadow_radius_matches_prediction() {
    let cfg = Config {
        width: 160,
        height: 90,
        samples: 1,
        serial: true,
        ..Config::default()
    };
    let buf = render(&cfg);
    // Scan the middle row for the black run around screen center.
    let row = cfg.height / 2;
    let is_black = |i: u32| {
        let idx = ((row * cfg.width + i) * 3) as usize;
        buf[idx] == 0 && buf[idx + 1] == 0 && buf[idx + 2] == 0
    };
    let center = cfg.width / 2;
    assert!(is_black(center), "screen center must be inside the shadow");
    let mut left = center;
    while left > 0 && is_black(left - 1) {
        left -= 1;
    }
    let mut right = center;
    while right < cfg.width - 1 && is_black(right + 1) {
        right += 1;
    }
    let measured = (right - left + 1) as f64 / 2.0;
    let sqrt_f = (1.0 - 2.0 / cfg.r_cam).sqrt();
    let alpha = (B_CRIT * sqrt_f / cfg.r_cam).asin();
    let predicted =
        cfg.width as f64 / 2.0 * alpha.tan() / (cfg.fov_deg.to_radians() / 2.0).tan();
    assert!(
        (measured - predicted).abs() < 1.5,
        "shadow radius {measured:.1} px, predicted {predicted:.1} px"
    );
}
