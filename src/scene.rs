//! Scene content: the accretion disk (and, later, the starfield background).

use crate::metric::R_ISCO;

/// Disk annulus: ISCO to the spec's outer edge, in units of M.
pub const DISK_INNER: f64 = R_ISCO;
pub const DISK_OUTER: f64 = 20.0;

/// Placeholder radius-gradient shading: hot near the ISCO fading to a dim
/// red rim. Replaced by the blackbody + redshift pipeline in the polish
/// stage.
pub fn disk_gradient(r: f64) -> [f64; 3] {
    let t = ((r - DISK_INNER) / (DISK_OUTER - DISK_INNER)).clamp(0.0, 1.0);
    let hot = [1.0, 0.93, 0.75];
    let cold = [0.45, 0.08, 0.01];
    std::array::from_fn(|k| hot[k] + t * (cold[k] - hot[k]))
}
