//! Per-pixel rendering. Each pixel is a pure function of (i, j, config), so
//! the parallel and serial paths are byte-identical by construction.

use rayon::prelude::*;

use crate::Config;
use crate::camera::Camera;
use crate::color::tonemap;
use crate::integrator::{RayOutcome, TraceParams, trace};
use crate::scene::{DISK_INNER, DISK_OUTER, Scene};

/// Escaping rays are integrated out to this radius before termination.
/// Classification is already certain at r ≈ 50 (no turning points outside the
/// photon sphere for outward rays), but the asymptotic direction — needed for
/// the starfield — keeps bending by over a degree beyond that.
pub const R_FAR: f64 = 2000.0;

pub fn render(cfg: &Config) -> Vec<u8> {
    let cam = Camera::new(
        cfg.r_cam,
        cfg.inclination_deg,
        cfg.azimuth_deg,
        cfg.fov_deg,
        cfg.width,
        cfg.height,
    );
    // Star PSF width: 1.5 pixel angular widths, floored at an intrinsic
    // ~0.0006 rad so extreme resolutions stay sane.
    let star_sigma = (1.5 * cfg.fov_deg.to_radians() / cfg.width as f64).max(6e-4);
    let scene = Scene::new(cfg.r_cam, cfg.exposure, star_sigma);
    let mut buf = vec![0u8; cfg.width as usize * cfg.height as usize * 3];
    let row_len = cfg.width as usize * 3;
    if cfg.serial {
        for (j, row) in buf.chunks_mut(row_len).enumerate() {
            render_row(&cam, &scene, cfg, j as u32, row);
        }
    } else {
        // Rows are disjoint slices and each pixel is a pure function of
        // (i, j, cfg), so the output is byte-identical to the serial path.
        buf.par_chunks_mut(row_len)
            .enumerate()
            .for_each(|(j, row)| render_row(&cam, &scene, cfg, j as u32, row));
    }
    buf
}

fn render_row(cam: &Camera, scene: &Scene, cfg: &Config, j: u32, row: &mut [u8]) {
    for i in 0..cfg.width {
        let px = shade_pixel(cam, scene, cfg, i, j);
        row[i as usize * 3..i as usize * 3 + 3].copy_from_slice(&px);
    }
}

fn shade_pixel(cam: &Camera, scene: &Scene, cfg: &Config, i: u32, j: u32) -> [u8; 3] {
    let n = cfg.samples;
    let mut acc = [0.0f64; 3];
    // Escaped samples are shaded after all traces: the angular scatter of a
    // pixel's own escape directions measures the local lensing magnification,
    // which sets how wide the star PSF must be to avoid aliasing into
    // speckle near the photon ring.
    let mut sky_dirs: Vec<crate::vec3::Vec3> = Vec::with_capacity((n * n) as usize);
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
            match trace(ray.y0, &params) {
                RayOutcome::Horizon | RayOutcome::MaxSteps => {}
                RayOutcome::Disk { state } => {
                    // z-axis angular momentum is conserved and computable at
                    // ray setup: L_z = L·(e1×e2)·ẑ.
                    let b_z = ray.l / ray.e * ray.e1.cross(ray.e2).z;
                    let c = scene.disk_radiance(state[0], b_z);
                    acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
                }
                RayOutcome::Escaped { state } => {
                    // Asymptotic direction reconstructed from the in-plane
                    // final state: dx/dλ = p_r·r̂ + r·p_φ·φ̂ in the (e1,e2)
                    // plane. Euclidean normalization is fine at r = 2000.
                    let [r, phi, p_r, p_phi] = state;
                    let (s, c) = phi.sin_cos();
                    let v = ray.e1 * (p_r * c - r * p_phi * s)
                        + ray.e2 * (p_r * s + r * p_phi * c);
                    sky_dirs.push(v.normalize());
                }
            }
        }
    }
    if !sky_dirs.is_empty() {
        // Angular radius of the sample cloud on the sky ≈ magnification ×
        // pixel size; the PSF must bridge the gap between adjacent samples
        // (cloud size / n). The starfield's cell search radius scales with
        // the widened PSF, so the cap here is just a cost guard.
        let cloud = sky_dirs[1..]
            .iter()
            .map(|d| d.dot(sky_dirs[0]).clamp(-1.0, 1.0).acos())
            .fold(0.0f64, f64::max);
        let spread = (cloud / n as f64 / scene.star_sigma).clamp(1.0, 25.0);
        for d in &sky_dirs {
            let c = scene.starfield(*d, spread);
            acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
        }
    }
    let inv = 1.0 / (n * n) as f64;
    tonemap([acc[0] * inv, acc[1] * inv, acc[2] * inv])
}
