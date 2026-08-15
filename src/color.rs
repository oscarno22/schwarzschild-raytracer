//! Blackbody temperature → linear sRGB with no data files: the Planck
//! spectrum is integrated against analytic fits to the CIE 1931 color
//! matching functions, normalized so temperature controls chroma only
//! (intensity is applied separately from the T⁴ law).

/// Piecewise Gaussian with separate widths left/right of the mean.
#[inline]
fn pw_gauss(x: f64, mu: f64, s_left: f64, s_right: f64) -> f64 {
    let s = if x < mu { s_left } else { s_right };
    let t = (x - mu) / s;
    (-0.5 * t * t).exp()
}

/// Wyman–Sloan–Shirley (JCGT 2013) multi-lobe fits to the CIE 1931
/// 2° color matching functions x̄, ȳ, z̄; λ in nm.
fn cie_xyz_bar(l: f64) -> [f64; 3] {
    [
        1.056 * pw_gauss(l, 599.8, 37.9, 31.0) + 0.362 * pw_gauss(l, 442.0, 16.0, 26.7)
            - 0.065 * pw_gauss(l, 501.1, 20.4, 26.2),
        0.821 * pw_gauss(l, 568.8, 46.9, 40.5) + 0.286 * pw_gauss(l, 530.9, 16.3, 31.1),
        1.217 * pw_gauss(l, 437.0, 11.8, 36.0) + 0.681 * pw_gauss(l, 459.0, 26.0, 13.8),
    ]
}

/// Planck spectral radiance at wavelength λ (nm) and temperature T (K),
/// arbitrary overall scale (we normalize by luminance afterwards).
/// c₂ = hc/k_B = 1.4388e7 nm·K.
fn planck(l_nm: f64, t: f64) -> f64 {
    1.0 / (l_nm.powi(5) * ((1.4388e7 / (l_nm * t)).exp_m1()))
}

/// Chromaticity of a blackbody at temperature T as linear sRGB, normalized
/// to luminance Y = 1 (negative out-of-gamut components clamped).
fn blackbody_rgb(t: f64) -> [f64; 3] {
    let mut xyz = [0.0f64; 3];
    let mut l = 380.0;
    while l <= 780.0 {
        let b = planck(l, t);
        let bar = cie_xyz_bar(l);
        for k in 0..3 {
            xyz[k] += b * bar[k];
        }
        l += 5.0;
    }
    let inv_y = 1.0 / xyz[1];
    let [x, y, z] = xyz.map(|c| c * inv_y);
    // XYZ → linear sRGB, D65.
    [
        (3.2406 * x - 1.5372 * y - 0.4986 * z).max(0.0),
        (-0.9689 * x + 1.8758 * y + 0.0415 * z).max(0.0),
        (0.0557 * x - 0.2040 * y + 1.0570 * z).max(0.0),
    ]
}

const LUT_SIZE: usize = 256;
const LUT_T_MIN: f64 = 500.0;
const LUT_T_MAX: f64 = 50_000.0;

/// Log-spaced in-memory lookup table over [500 K, 50 000 K], built once at
/// startup. Out-of-range temperatures clamp (chromaticity saturates toward
/// deep red / blue well inside the range anyway).
pub struct BlackbodyLut {
    entries: Vec<[f64; 3]>,
}

impl BlackbodyLut {
    pub fn new() -> Self {
        let entries = (0..LUT_SIZE)
            .map(|i| {
                let t = LUT_T_MIN
                    * (LUT_T_MAX / LUT_T_MIN).powf(i as f64 / (LUT_SIZE - 1) as f64);
                blackbody_rgb(t)
            })
            .collect();
        Self { entries }
    }

    pub fn sample(&self, t: f64) -> [f64; 3] {
        let x = (t / LUT_T_MIN).ln() / (LUT_T_MAX / LUT_T_MIN).ln()
            * (LUT_SIZE - 1) as f64;
        let x = x.clamp(0.0, (LUT_SIZE - 1) as f64);
        let i = (x as usize).min(LUT_SIZE - 2);
        let frac = x - i as f64;
        let (a, b) = (self.entries[i], self.entries[i + 1]);
        std::array::from_fn(|k| a[k] + frac * (b[k] - a[k]))
    }
}

impl Default for BlackbodyLut {
    fn default() -> Self {
        Self::new()
    }
}

/// HDR linear radiance → display u8: Reinhard x/(1+x) per channel, then the
/// sRGB transfer function.
pub fn tonemap(c: [f64; 3]) -> [u8; 3] {
    c.map(|v| {
        let v = (v.max(0.0)) / (1.0 + v.max(0.0));
        let s = if v <= 0.003_130_8 {
            12.92 * v
        } else {
            1.055 * v.powf(1.0 / 2.4) - 0.055
        };
        (s * 255.0).round().clamp(0.0, 255.0) as u8
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (g) Color sanity: 6500 K is near-neutral, 3000 K is reddish,
    /// 20 000 K is bluish.
    #[test]
    fn blackbody_chromaticity_sanity() {
        let lut = BlackbodyLut::new();
        let neutral = lut.sample(6500.0);
        for k in 1..3 {
            assert!(
                (neutral[k] / neutral[0] - 1.0).abs() < 0.15,
                "6500 K should be near-neutral, got {neutral:?}"
            );
        }
        let warm = lut.sample(3000.0);
        assert!(warm[0] > warm[1] && warm[1] > warm[2], "3000 K: {warm:?}");
        let hot = lut.sample(20_000.0);
        assert!(hot[2] > hot[0], "20 000 K: {hot:?}");
    }

    #[test]
    fn lut_matches_direct_evaluation() {
        let lut = BlackbodyLut::new();
        for t in [800.0, 4500.0, 10_000.0, 42_000.0] {
            let direct = blackbody_rgb(t);
            let sampled = lut.sample(t);
            for k in 0..3 {
                assert!((direct[k] - sampled[k]).abs() < 0.02);
            }
        }
    }
}
