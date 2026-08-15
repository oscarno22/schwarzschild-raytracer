//! Minimal f64 3-vector. Hand-rolled on purpose — the only geometry this
//! project needs is dot/cross/normalize for camera setup and ray-plane bases.

use std::ops::{Add, Mul, Neg, Sub};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vec3 {
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    pub fn cross(self, o: Self) -> Self {
        Self::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    pub fn length(self) -> f64 {
        self.dot(self).sqrt()
    }

    pub fn normalize(self) -> Self {
        self * (1.0 / self.length())
    }

    /// An arbitrary unit vector perpendicular to `self` (assumed unit length).
    /// Crosses with whichever global axis is least parallel, so it never
    /// degenerates.
    pub fn any_perpendicular(self) -> Self {
        let axis = if self.x.abs() < 0.9 {
            Self::new(1.0, 0.0, 0.0)
        } else {
            Self::new(0.0, 1.0, 0.0)
        };
        self.cross(axis).normalize()
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.x, -self.y, -self.z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cross_is_right_handed() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.cross(y), Vec3::new(0.0, 0.0, 1.0));
    }

    #[test]
    fn normalize_gives_unit_length() {
        let v = Vec3::new(3.0, -4.0, 12.0).normalize();
        assert!((v.length() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn any_perpendicular_is_perpendicular() {
        for v in [
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.6, -0.48, 0.64),
        ] {
            let p = v.any_perpendicular();
            assert!(v.dot(p).abs() < 1e-15);
            assert!((p.length() - 1.0).abs() < 1e-15);
        }
    }
}
