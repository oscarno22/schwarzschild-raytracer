//! Image rendering through the δ-table: the fast approximate path used by
//! the real-time viewer, sharing all shading code (Scene, Spot, tonemap)
//! with the offline renderer so the two agree pixel-for-pixel away from the
//! photon ring (enforced in tests/table.rs).

use rayon::prelude::*;

use crate::Config;
use crate::camera::Camera;
use crate::color::tonemap;
use crate::deltatable::{DeltaTable, Lookup};
use crate::scene::{Scene, Spot};

/// How sky pixels pick their star-PSF widening.
pub enum SkySpread {
    /// spread = 1.0 exactly — what the offline renderer computes at
    /// samples = 1 (its sample-cloud measurement is empty). Parity tests
    /// use this so sky differences are purely direction error.
    Reference,
    /// spread from the tabulated magnification |Δψ/Δε|: the table's stand-in
    /// for the offline sample-cloud measurement, used by the live viewer.
    Magnification,
}

/// A view + scene bound to a prebuilt δ-table, shadeable at any frame time.
pub struct TableFrame<'t> {
    table: &'t DeltaTable,
    cam: Camera,
    scene: Scene,
    width: u32,
    height: u32,
    samples: u32,
    serial: bool,
    /// Precomputed pixel angular size / star_sigma factor for Magnification.
    mag_scale: f64,
    spread: SkySpread,
}

impl<'t> TableFrame<'t> {
    pub fn new(cfg: &Config, table: &'t DeltaTable, spread: SkySpread) -> Self {
        debug_assert!(
            (cfg.r_cam - table.r_cam).abs() < 1e-12,
            "table built for a different r_cam"
        );
        let cam = Camera::new(
            cfg.r_cam,
            cfg.inclination_deg,
            cfg.azimuth_deg,
            cfg.fov_deg,
            cfg.width,
            cfg.height,
        );
        // Same star PSF floor as CachedRender::new.
        let star_sigma = (1.5 * cfg.fov_deg.to_radians() / cfg.width as f64).max(6e-4);
        let spot = (cfg.spot_amp != 0.0).then(|| Spot::new(cfg.spot_r, cfg.spot_amp));
        let scene = Scene::new(cfg.r_cam, cfg.exposure, star_sigma, spot, cfg.profile);
        let pixel_ang = cfg.fov_deg.to_radians() / cfg.width as f64;
        Self {
            table,
            cam,
            scene,
            width: cfg.width,
            height: cfg.height,
            samples: cfg.samples,
            serial: cfg.serial,
            mag_scale: 0.5 * pixel_ang / star_sigma,
            spread,
        }
    }

    /// Shade every pixel at frame time `t_frame` into an RGB8 buffer.
    pub fn shade_rgb(&self, t_frame: f64) -> Vec<u8> {
        let mut buf = vec![0u8; self.width as usize * self.height as usize * 3];
        let row_len = self.width as usize * 3;
        let shade_row = |j: usize, row: &mut [u8]| {
            for i in 0..self.width as usize {
                let px = self.shade_pixel(i as u32, j as u32, t_frame);
                row[i * 3..i * 3 + 3].copy_from_slice(&px);
            }
        };
        if self.serial {
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

    /// Shade into a packed 0RGB u32 buffer (softbuffer's little-endian
    /// format: r << 16 | g << 8 | b), parallel over rows.
    pub fn shade_0rgb(&self, t_frame: f64, out: &mut [u32]) {
        assert_eq!(out.len(), self.width as usize * self.height as usize);
        out.par_chunks_mut(self.width as usize)
            .enumerate()
            .for_each(|(j, row)| {
                for (i, px) in row.iter_mut().enumerate() {
                    let [r, g, b] = self.shade_pixel(i as u32, j as u32, t_frame);
                    *px = (r as u32) << 16 | (g as u32) << 8 | b as u32;
                }
            });
    }

    fn shade_pixel(&self, i: u32, j: u32, t_frame: f64) -> [u8; 3] {
        let n = self.samples;
        let mut acc = [0.0f64; 3];
        for a in 0..n {
            for b in 0..n {
                let sx = (a as f64 + 0.5) / n as f64;
                let sy = (b as f64 + 0.5) / n as f64;
                let ray = self.cam.make_ray(self.cam.pixel_dir(i, j, sx, sy));
                let eps = (ray.l / self.table.r_cam).atan2(-ray.y0[2] / self.table.sqrt_f);
                let c = match self.table.sample(eps, ray.e1.z, ray.e2.z) {
                    Lookup::Captured => [0.0; 3],
                    Lookup::Disk { r, phi_p, t } => {
                        let b_z = ray.l / ray.e * ray.e1.cross(ray.e2).z;
                        let (s, c) = phi_p.sin_cos();
                        let x = ray.e1 * c + ray.e2 * s;
                        let phi_hit = x.y.atan2(x.x);
                        let boost = match &self.scene.spot {
                            Some(spot) => spot.temp_boost(r, phi_hit, t_frame - t),
                            None => 1.0,
                        };
                        self.scene.disk_radiance(r, b_z, boost)
                    }
                    Lookup::Sky { psi_inf, dpsi_deps } => {
                        let (s, c) = psi_inf.sin_cos();
                        let dir = ray.e1 * c + ray.e2 * s;
                        let spread = match self.spread {
                            SkySpread::Reference => 1.0,
                            SkySpread::Magnification => {
                                (dpsi_deps * self.mag_scale).clamp(1.0, 25.0)
                            }
                        };
                        self.scene.starfield(dir, spread)
                    }
                };
                acc = [acc[0] + c[0], acc[1] + c[1], acc[2] + c[2]];
            }
        }
        let inv = 1.0 / (n * n) as f64;
        tonemap([acc[0] * inv, acc[1] * inv, acc[2] * inv])
    }
}
