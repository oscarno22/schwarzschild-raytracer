//! Scene content: accretion disk emission with the full relativistic
//! treatment — thin-disk temperature profile, gravitational + Doppler
//! frequency shift, and g⁴ beaming.

use crate::color::BlackbodyLut;
use crate::metric::R_ISCO;

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
    lut: BlackbodyLut,
}

impl Scene {
    pub fn new(r_cam: f64, exposure: f64) -> Self {
        Self {
            sqrt_f_cam: (1.0 - 2.0 / r_cam).sqrt(),
            exposure,
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
}
