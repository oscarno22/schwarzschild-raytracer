//! Pinhole camera: pixel → view direction → geodesic initial conditions.
//!
//! Screen angles are angles in the local orthonormal frame of a static
//! observer at the camera (tetrad e_r̂ = √f ∂_r, e_φ̂ = (1/r) ∂_φ). The
//! overall null-vector normalization is an affine convention (we pick the
//! locally-measured frequency = 1), but the √f factor relating radial and
//! transverse components is physical — dropping it inflates the shadow by
//! √(1/f) ≈ 3.4% at r = 30 while looking entirely plausible.

use crate::metric::State;
use crate::vec3::Vec3;

pub struct Camera {
    pub pos: Vec3,
    pub r_cam: f64,
    /// √(1 − 2/r_cam); also the conserved E of every ray under our
    /// normalization.
    pub sqrt_f: f64,
    forward: Vec3,
    right: Vec3,
    up: Vec3,
    tan_half_fov: f64,
    width: f64,
    height: f64,
}

/// A pixel's ray: geodesic initial state plus the orbital-plane basis
/// needed to reconstruct global 3D positions and directions.
pub struct Ray {
    /// [r, φ_plane, dr/dλ, dφ/dλ] at the camera.
    pub y0: State,
    /// Conserved energy (= camera √f under local normalization p^t̂ = 1).
    pub e: f64,
    /// Conserved in-plane angular momentum L = r²·dφ/dλ.
    pub l: f64,
    /// Radial basis vector of the ray's orbital plane (unit, at the camera).
    pub e1: Vec3,
    /// Transverse basis vector; global position is r(cos φ·e1 + sin φ·e2).
    pub e2: Vec3,
}

impl Camera {
    pub fn new(
        r_cam: f64,
        inclination_deg: f64,
        azimuth_deg: f64,
        fov_deg: f64,
        width: u32,
        height: u32,
    ) -> Self {
        let inc = inclination_deg.to_radians();
        let az = azimuth_deg.to_radians();
        let pos = Vec3::new(inc.sin() * az.cos(), inc.sin() * az.sin(), inc.cos()) * r_cam;
        let forward = (-pos).normalize();
        let cross = forward.cross(Vec3::new(0.0, 0.0, 1.0));
        // Looking straight down the polar axis leaves "up" unconstrained.
        let right = if cross.length() < 1e-9 {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            cross.normalize()
        };
        let up = right.cross(forward);
        Self {
            pos,
            r_cam,
            sqrt_f: (1.0 - 2.0 / r_cam).sqrt(),
            forward,
            right,
            up,
            tan_half_fov: (fov_deg.to_radians() / 2.0).tan(),
            width: width as f64,
            height: height as f64,
        }
    }

    /// View direction for pixel (i, j) with intra-pixel offset (sx, sy) ∈ [0,1).
    /// Gnomonic projection: horizontal FOV spans the width, square pixels.
    pub fn pixel_dir(&self, i: u32, j: u32, sx: f64, sy: f64) -> Vec3 {
        let u = (2.0 * (i as f64 + sx) / self.width - 1.0) * self.tan_half_fov;
        let v = (1.0 - 2.0 * (j as f64 + sy) / self.height)
            * self.tan_half_fov
            * (self.height / self.width);
        (self.forward + self.right * u + self.up * v).normalize()
    }

    /// Geodesic initial conditions for unit view direction d.
    ///
    /// With δ the angle between d and the outward radial direction:
    ///   dr/dλ = cos δ·√f,  dφ/dλ = sin δ / r,  E = √f,  L = r·sin δ
    /// which satisfies the null condition exactly and gives b = L/E =
    /// r sin δ/√f.
    pub fn make_ray(&self, d: Vec3) -> Ray {
        let e1 = self.pos * (1.0 / self.r_cam);
        let cos_d = d.dot(e1);
        let d_perp = d - e1 * cos_d;
        let sin_d = d_perp.length();
        // A ray dead-on radial has no defined orbital plane; any plane
        // through the radial line works (the geodesic never leaves it).
        let e2 = if sin_d < 1e-12 {
            e1.any_perpendicular()
        } else {
            d_perp * (1.0 / sin_d)
        };
        Ray {
            y0: [
                self.r_cam,
                0.0,
                cos_d * self.sqrt_f,
                sin_d / self.r_cam,
            ],
            e: self.sqrt_f,
            l: self.r_cam * sin_d,
            e1,
            e2,
        }
    }
}
