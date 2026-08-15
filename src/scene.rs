//! Scene content: accretion disk emission with the full relativistic
//! treatment — thin-disk temperature profile, gravitational + Doppler
//! frequency shift, g⁴ beaming — plus a deterministic procedural starfield.

use crate::color::BlackbodyLut;
use crate::metric::R_ISCO;
use crate::vec3::Vec3;

/// Disk annulus: ISCO to the spec's outer edge, in units of M.
pub const DISK_INNER: f64 = R_ISCO;
pub const DISK_OUTER: f64 = 20.0;

/// Disk temperature at the inner edge. Sets the color scale; the shape is
/// the standard thin-disk scaling T ∝ r^(−3/4). 7000 K puts the outer disk
/// at a deep orange ~2800 K and lets Doppler beaming push the approaching
/// inner edge to a white-blue ~15 000 K.
pub const T_ISCO: f64 = 7000.0;

pub struct Scene {
    /// √(1 − 2/r_cam) — gravitational blueshift factor of the static camera.
    pub sqrt_f_cam: f64,
    /// Overall intensity scale applied to disk emission before tonemapping.
    pub exposure: f64,
    /// Angular radius (radians) of a star's Gaussian point-spread falloff.
    /// Floored at ~1.5 pixel widths by the caller so stars neither shimmer
    /// nor drop out between neighboring pixels.
    pub star_sigma: f64,
    lut: BlackbodyLut,
}

impl Scene {
    pub fn new(r_cam: f64, exposure: f64, star_sigma: f64) -> Self {
        Self {
            sqrt_f_cam: (1.0 - 2.0 / r_cam).sqrt(),
            exposure,
            star_sigma,
            lut: BlackbodyLut::new(),
        }
    }

    /// Thin-disk temperature profile, hottest at the ISCO.
    /// (The Novikov–Thorne profile T ∝ [r⁻³(1 − √(6/r))]^(1/4) is the
    /// drop-in alternative; it vanishes at the ISCO and peaks near r ≈ 8.2.)
    pub fn disk_temperature(r: f64) -> f64 {
        T_ISCO * (r / R_ISCO).powf(-0.75)
    }

    /// Combined gravitational + Doppler shift g = ν_obs/ν_em for an emitter
    /// on a circular equatorial geodesic orbit (Ω = r^(−3/2), u^t =
    /// 1/√(1 − 3/r)) observed by the static camera:
    ///
    ///   g = √(1 − 3/r) / (√f_cam · (1 + Ω·b_z))
    ///
    /// where b_z = L_z/E of the *traced* (backward) ray. The physical photon
    /// runs the path in reverse, so its conserved L_z is the negative of the
    /// traced one — hence `+` here where the textbook formula (in terms of
    /// the photon's own angular momentum) has `−`. A traced ray aimed at the
    /// receding side of the disk has b_z > 0 and must come out redshifted.
    pub fn redshift(&self, r: f64, b_z_traced: f64) -> f64 {
        (1.0 - 3.0 / r).sqrt() / (self.sqrt_f_cam * (1.0 + r.powf(-1.5) * b_z_traced))
    }

    /// Observed linear-sRGB radiance of the disk at hit radius r, for a
    /// traced ray with z-axis impact parameter b_z.
    ///
    /// A blackbody at T_em seen with shift g is exactly a blackbody at
    /// g·T_em, and the frequency-integrated intensity transforms as g⁴
    /// (I/ν³ invariance) — with I_em ∝ T_em⁴ both effects reduce to
    /// intensity ∝ T_obs⁴. Beaming and color stay mutually consistent.
    pub fn disk_radiance(&self, r: f64, b_z_traced: f64) -> [f64; 3] {
        let g = self.redshift(r, b_z_traced);
        let t_obs = g * Self::disk_temperature(r);
        let intensity = (t_obs / T_ISCO).powi(4) * self.exposure;
        self.lut.sample(t_obs).map(|c| c * intensity)
    }

    /// Deterministic hash-based starfield, sampled at the escaped ray's
    /// asymptotic direction. Stars live in cells of a 256×512 lat-long grid;
    /// presence probability ∝ sin θ cancels the grid's pole clustering, so
    /// the sky is statistically uniform. Everything derives from splitmix64
    /// of the cell id — no RNG state, no data files, resolution-independent.
    ///
    /// `spread` (≥ 1) widens the Gaussian PSF for strongly-lensed rays,
    /// whose pixel footprint on the sky is sheared far beyond one pixel —
    /// without it, the region around the photon ring aliases into speckle.
    /// The widening conserves energy (peak ∝ 1/spread²), and the cell
    /// search radius grows with the PSF so wide kernels see every star
    /// that contributes.
    pub fn starfield(&self, dir: Vec3, spread: f64) -> [f64; 3] {
        const BANDS: i64 = 256;
        const COLS: i64 = 512;
        const DENSITY: f64 = 0.12;
        const SEED: u64 = 0x0B5E12FE_D51C_E5EE;
        let theta = dir.z.clamp(-1.0, 1.0).acos();
        let phi = dir.y.atan2(dir.x).rem_euclid(std::f64::consts::TAU);
        let band = ((theta / std::f64::consts::PI) * BANDS as f64) as i64;
        let col = ((phi / std::f64::consts::TAU) * COLS as f64) as i64;
        let sigma = self.star_sigma * spread;
        let band_height = std::f64::consts::PI / BANDS as f64;
        let col_width = std::f64::consts::TAU / COLS as f64 * theta.sin().max(0.15);
        let rb = ((2.5 * sigma / band_height).ceil() as i64).clamp(1, 6);
        let rc = ((2.5 * sigma / col_width).ceil() as i64).clamp(1, 12);
        // The deep-space floor keeps escaped pixels strictly nonzero, so
        // pure black remains an exact marker for captured rays.
        let mut acc = [1.2e-3, 1.2e-3, 1.8e-3];
        for db in -rb..=rb {
            let b = band + db;
            if !(0..BANDS).contains(&b) {
                continue;
            }
            for dc in -rc..=rc {
                let c = (col + dc).rem_euclid(COLS);
                let h0 = splitmix64((b * COLS + c) as u64 ^ SEED);
                let theta_c = (b as f64 + 0.5) / BANDS as f64 * std::f64::consts::PI;
                if hash01(h0) > DENSITY * theta_c.sin() {
                    continue;
                }
                let h1 = splitmix64(h0);
                let h2 = splitmix64(h1);
                let h3 = splitmix64(h2);
                let h4 = splitmix64(h3);
                let theta_s = (b as f64 + hash01(h1)) / BANDS as f64 * std::f64::consts::PI;
                let phi_s = (c as f64 + hash01(h2)) / COLS as f64 * std::f64::consts::TAU;
                let star = Vec3::new(
                    theta_s.sin() * phi_s.cos(),
                    theta_s.sin() * phi_s.sin(),
                    theta_s.cos(),
                );
                let psi = dir.dot(star).clamp(-1.0, 1.0).acos();
                // Energy-conserving widening: peak drops as the PSF spreads.
                let falloff = (-(psi / sigma).powi(2)).exp() / (spread * spread);
                if falloff < 1e-6 {
                    continue;
                }
                // Power-law brightness (few bright, many faint) and a
                // temperature distribution skewed toward cool stars.
                let brightness = (0.08 * hash01(h3).powf(-0.9)).min(8.0);
                let t_star = 2200.0 * 10.0f64.powf(1.3 * hash01(h4).powi(2));
                let rgb = self.lut.sample(t_star);
                for k in 0..3 {
                    acc[k] += rgb[k] * brightness * falloff;
                }
            }
        }
        acc
    }
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D049BB133111EB);
    x ^ (x >> 31)
}

/// Uniform f64 in [0, 1) from the high 53 bits of a hash.
fn hash01(x: u64) -> f64 {
    (x >> 11) as f64 / (1u64 << 53) as f64
}
