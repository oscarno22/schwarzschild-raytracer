//! Per-pixel rendering. Each pixel is a pure function of (i, j, config), so
//! the parallel and serial paths are byte-identical by construction.
//!
//! Tracing and shading are split: `CachedRender` traces every sample once and
//! stores the per-pixel hit data (disk hits with retarded emission time, and
//! already-shaded starfield colors — the sky is time-independent). Re-shading
//! a frame then only re-evaluates disk radiance, which is what makes
//! fixed-camera animation (`--frames`) ~100× cheaper than re-tracing.

use rayon::prelude::*;

use crate::Config;
use crate::camera::Camera;
use crate::color::tonemap;
use crate::integrator::{RayOutcome, TraceParams, trace};
use crate::scene::{DISK_INNER, DISK_OUTER, Scene, Spot};

/// Escaping rays are integrated out to this radius before termination.
/// Classification is already certain at r ≈ 50 (no turning points outside the
/// photon sphere for outward rays), but the asymptotic direction — needed for
/// the starfield — keeps bending by over a degree beyond that.
pub const R_FAR: f64 = 2000.0;

/// One disk-hit sample: everything shading needs, camera geometry already
/// resolved to global disk coordinates.
pub struct DiskHit {
    /// Hit radius in the equatorial plane.
    pub r: f64,
    /// Global azimuth of the hit point (atan2 of the 3D position).
    pub phi: f64,
    /// z-axis impact parameter L_z/E of the traced ray (sets the Doppler
    /// shift sign: positive = aimed at the receding side).
    pub b_z: f64,
    /// Emission coordinate time relative to the camera's t = 0 (negative:
    /// the light left the disk before it arrived).
    pub t_emit: f64,
}

/// Per-pixel trace results. Contributions are stored in the exact order the
/// straight-through renderer accumulates them (disk hits in sample order,
/// then shaded sky colors), so re-shading reproduces its bytes exactly.
struct PixelCache {
    hits: Vec<DiskHit>,
    sky: Vec<[f64; 3]>,
}

/// A fully traced image that can be re-shaded at any frame time.
pub struct CachedRender {
    cfg: Config,
    scene: Scene,
    pixels: Vec<PixelCache>,
}

impl CachedRender {
    pub fn new(cfg: &Config) -> Self {
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
        let spot = (cfg.spot_amp != 0.0).then(|| Spot::new(cfg.spot_r, cfg.spot_amp));
        let scene = Scene::new(cfg.r_cam, cfg.exposure, star_sigma, spot, cfg.profile);
        let n_px = cfg.width as usize * cfg.height as usize;
        let mut pixels = Vec::with_capacity(n_px);
        let trace_row = |j: u32, row: &mut [PixelCache]| {
            for (i, px) in row.iter_mut().enumerate() {
                *px = trace_pixel(&cam, &scene, cfg, i as u32, j);
            }
        };
        // Placeholder-fill then trace by disjoint rows; each pixel is a pure
        // function of (i, j, cfg), so serial and parallel agree exactly.
        pixels.resize_with(n_px, || PixelCache {
            hits: Vec::new(),
            sky: Vec::new(),
        });
        let row_len = cfg.width as usize;
        if cfg.serial {
            for (j, row) in pixels.chunks_mut(row_len).enumerate() {
                trace_row(j as u32, row);
            }
        } else {
            pixels
                .par_chunks_mut(row_len)
                .enumerate()
                .for_each(|(j, row)| trace_row(j as u32, row));
        }
        Self {
            cfg: cfg.clone(),
            scene,
            pixels,
        }
    }

    /// Shade every pixel at frame coordinate time `t_frame`.
    pub fn shade(&self, t_frame: f64) -> Vec<u8> {
        let cfg = &self.cfg;
        let mut buf = vec![0u8; cfg.width as usize * cfg.height as usize * 3];
        let row_len = cfg.width as usize * 3;
        let shade_row = |j: usize, row: &mut [u8]| {
            for i in 0..cfg.width as usize {
                let px = shade_cached(
                    &self.scene,
                    &self.pixels[j * cfg.width as usize + i],
                    t_frame,
                    cfg.samples,
                );
                row[i * 3..i * 3 + 3].copy_from_slice(&px);
            }
        };
        if cfg.serial {
            for (j, row) in buf.chunks_mut(row_len).enumerate() {
                shade_row(j, row);
            }
        } else {
            buf.par_chunks_mut(row_len)
                .enumerate()
                .for_each(|(j, row)| shade_row(j, row));
        }
        buf
    }
}

/// Single-frame render at cfg.time (the frames mode shares this path, so a
/// frame emitted by `CachedRender` is identical to a one-shot invocation).
pub fn render(cfg: &Config) -> Vec<u8> {
    CachedRender::new(cfg).shade(cfg.time)
}

fn trace_pixel(cam: &Camera, scene: &Scene, cfg: &Config, i: u32, j: u32) -> PixelCache {
    let n = cfg.samples;
    let mut hits: Vec<DiskHit> = Vec::new();
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
                    // Global position of the hit: r(cos φ_p·e1 + sin φ_p·e2).
                    let (s, c) = state[1].sin_cos();
                    let x = ray.e1 * c + ray.e2 * s;
                    hits.push(DiskHit {
                        r: state[0],
                        phi: x.y.atan2(x.x),
                        b_z,
                        // The trace runs backward with t increasing, so the
                        // photon left the disk state[4] before camera time.
                        t_emit: -state[4],
                    });
                }
                RayOutcome::Escaped { state } => {
                    // Asymptotic direction reconstructed from the in-plane
                    // final state: dx/dλ = p_r·r̂ + r·p_φ·φ̂ in the (e1,e2)
                    // plane. Euclidean normalization is fine at r = 2000.
                    let [r, phi, p_r, p_phi, _t] = state;
                    let (s, c) = phi.sin_cos();
                    let v = ray.e1 * (p_r * c - r * p_phi * s)
                        + ray.e2 * (p_r * s + r * p_phi * c);
                    sky_dirs.push(v.normalize());
                }
            }
        }
    }
    let mut sky: Vec<[f64; 3]> = Vec::with_capacity(sky_dirs.len());
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
            sky.push(scene.starfield(*d, spread));
        }
    }
    PixelCache { hits, sky }
}

fn shade_cached(scene: &Scene, px: &PixelCache, t_frame: f64, samples: u32) -> [u8; 3] {
    let mut acc = [0.0f64; 3];
    for hit in &px.hits {
        let boost = match &scene.spot {
            // Retarded time: the spot is where its orbit had it when the
            // light left, t_emit (< 0) before this frame's camera time.
            Some(spot) => spot.temp_boost(hit.r, hit.phi, t_frame + hit.t_emit),
            None => 1.0,
        };
        let c = scene.disk_radiance(hit.r, hit.b_z, boost);
        acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
    }
    for c in &px.sky {
        acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
    }
    let inv = 1.0 / (samples * samples) as f64;
    tonemap([acc[0] * inv, acc[1] * inv, acc[2] * inv])
}
