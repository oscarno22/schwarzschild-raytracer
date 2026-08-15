//! Per-pixel rendering. Each pixel is a pure function of (i, j, config), so
//! the parallel and serial paths are byte-identical by construction.

use crate::Config;
use crate::camera::Camera;
use crate::integrator::{RayOutcome, TraceParams, trace};
use crate::scene::{DISK_INNER, DISK_OUTER, disk_gradient};

/// Escaping rays are integrated out to this radius before termination.
/// Classification is already certain at r ≈ 50 (no turning points outside the
/// photon sphere for outward rays), but the asymptotic direction — needed for
/// the starfield — keeps bending by over a degree beyond that.
pub const R_FAR: f64 = 2000.0;

pub fn render(cfg: &Config) -> Vec<u8> {
    let cam = Camera::new(
        cfg.r_cam,
        cfg.inclination_deg,
        cfg.fov_deg,
        cfg.width,
        cfg.height,
    );
    let mut buf = vec![0u8; cfg.width as usize * cfg.height as usize * 3];
    let row_len = cfg.width as usize * 3;
    for (j, row) in buf.chunks_mut(row_len).enumerate() {
        render_row(&cam, cfg, j as u32, row);
    }
    buf
}

fn render_row(cam: &Camera, cfg: &Config, j: u32, row: &mut [u8]) {
    for i in 0..cfg.width {
        let px = shade_pixel(cam, cfg, i, j);
        row[i as usize * 3..i as usize * 3 + 3].copy_from_slice(&px);
    }
}

fn shade_pixel(cam: &Camera, cfg: &Config, i: u32, j: u32) -> [u8; 3] {
    let n = cfg.samples;
    let mut acc = [0.0f64; 3];
    for a in 0..n {
        for b in 0..n {
            let sx = (a as f64 + 0.5) / n as f64;
            let sy = (b as f64 + 0.5) / n as f64;
            let ray = cam.make_ray(cam.pixel_dir(i, j, sx, sy));
            let params = TraceParams {
                e: ray.e,
                step_scale: cfg.step_scale,
                max_steps: cfg.max_steps,
                r_far: R_FAR,
                plane_az: ray.e1.z,
                plane_bz: ray.e2.z,
                disk: Some((DISK_INNER, DISK_OUTER)),
            };
            let c = match trace(ray.y0, &params) {
                RayOutcome::Horizon | RayOutcome::MaxSteps => [0.0; 3],
                RayOutcome::Disk { state } => disk_gradient(state[0]),
                // Flat gray placeholder background; starfield lands later.
                RayOutcome::Escaped { .. } => [0.2; 3],
            };
            acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
        }
    }
    let inv = 1.0 / (n * n) as f64;
    std::array::from_fn(|k| ((acc[k] * inv).clamp(0.0, 1.0) * 255.0).round() as u8)
}
