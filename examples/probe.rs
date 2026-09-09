//! Debug probe: dump per-sample trace outcomes for chosen pixels.
//! Usage: cargo run --release --example probe -- <i> <j> [width height]

use schwarzschild_raytracer::Config;
use schwarzschild_raytracer::camera::Camera;
use schwarzschild_raytracer::integrator::{RayOutcome, TraceParams, trace};
use schwarzschild_raytracer::render::R_FAR;
use schwarzschild_raytracer::scene::{DISK_INNER, DISK_OUTER};

fn main() {
    let args: Vec<f64> = std::env::args()
        .skip(1)
        .map(|s| s.parse().expect("numeric args"))
        .collect();
    let (pi, pj) = (args[0] as u32, args[1] as u32);
    let cfg = Config {
        width: *args.get(2).unwrap_or(&800.0) as u32,
        height: *args.get(3).unwrap_or(&450.0) as u32,
        samples: 2,
        ..Config::default()
    };
    let cam = Camera::new(
        cfg.r_cam,
        cfg.inclination_deg,
        cfg.azimuth_deg,
        cfg.fov_deg,
        cfg.width,
        cfg.height,
    );
    let n = cfg.samples;
    let mut dirs = vec![];
    for a in 0..n {
        for b in 0..n {
            let sx = (a as f64 + 0.5) / n as f64;
            let sy = (b as f64 + 0.5) / n as f64;
            let ray = cam.make_ray(cam.pixel_dir(pi, pj, sx, sy));
            let params = TraceParams {
                e: ray.e,
                step_scale: cfg.step_scale,
                max_steps: cfg.max_steps,
                r_far: R_FAR,
                plane_az: ray.e1.z,
                plane_bz: ray.e2.z,
                disk: Some((DISK_INNER, DISK_OUTER)),
            };
            match trace(ray.y0, &params) {
                RayOutcome::Horizon => println!("({a},{b}) horizon"),
                RayOutcome::MaxSteps => println!("({a},{b}) max-steps"),
                RayOutcome::Disk { state } => {
                    println!("({a},{b}) disk r={:.3} phi={:.3}", state[0], state[1])
                }
                RayOutcome::Escaped { state } => {
                    let [r, phi, p_r, p_phi, _t] = state;
                    let (s, c) = phi.sin_cos();
                    let v = (ray.e1 * (p_r * c - r * p_phi * s)
                        + ray.e2 * (p_r * s + r * p_phi * c))
                        .normalize();
                    println!(
                        "({a},{b}) sky dir=({:+.4},{:+.4},{:+.4}) wind={:.2}pi",
                        v.x,
                        v.y,
                        v.z,
                        phi / std::f64::consts::PI
                    );
                    dirs.push(v);
                }
            }
        }
    }
    if dirs.len() > 1 {
        let cloud = dirs[1..]
            .iter()
            .map(|d| d.dot(dirs[0]).clamp(-1.0, 1.0).acos())
            .fold(0.0f64, f64::max);
        let px_ang = cfg.fov_deg.to_radians() / cfg.width as f64;
        println!("cloud={cloud:.5} rad = {:.1} px", cloud / px_ang);
    }
}
