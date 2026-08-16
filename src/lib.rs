//! Schwarzschild black hole raytracer: backward ray tracing through null
//! geodesics of the Schwarzschild metric (M = 1 geometric units).

pub mod camera;
pub mod color;
pub mod deltatable;
pub mod integrator;
pub mod metric;
pub mod render;
pub mod scene;
pub mod vec3;

/// Render configuration, populated from CLI flags. All distances in units of M.
#[derive(Debug, Clone)]
pub struct Config {
    pub width: u32,
    pub height: u32,
    pub output: String,
    pub r_cam: f64,
    pub inclination_deg: f64,
    pub azimuth_deg: f64,
    pub fov_deg: f64,
    pub samples: u32,
    pub serial: bool,
    pub max_steps: u32,
    pub step_scale: f64,
    pub exposure: f64,
    /// Coordinate time of the frame (camera clock sits at t = 0 + this).
    pub time: f64,
    /// Hot-spot temperature amplitude; 0 disables the spot entirely.
    pub spot_amp: f64,
    /// Orbital radius of the hot spot (circular geodesic, Ω = r^(−3/2)).
    pub spot_r: f64,
    /// Number of frames to emit (fixed camera: trace once, re-shade N times).
    pub frames: u32,
    /// Coordinate-time step between frames.
    pub frame_dt: f64,
    /// Disk temperature law: thin (T ∝ r^(−3/4)) or Novikov–Thorne.
    pub profile: scene::DiskProfile,
    /// Color (physical shading) or Echo (retarded-time false color).
    pub render_mode: RenderMode,
}

/// What each pixel displays.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenderMode {
    /// Physical shading: blackbody disk + starfield.
    Color,
    /// False color of the disk hit's light-travel delay Δt — the light-echo
    /// geometry made directly visible (far side and wrapped images are
    /// measurably "older"). Sky renders dark navy, captured rays black.
    Echo,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            width: 1920,
            height: 1080,
            output: "render.png".to_string(),
            r_cam: 30.0,
            inclination_deg: 80.0,
            azimuth_deg: 0.0,
            fov_deg: 75.0,
            samples: 2,
            serial: false,
            max_steps: 60_000,
            step_scale: 0.02,
            exposure: 4.0,
            time: 0.0,
            spot_amp: 0.0,
            spot_r: 7.0,
            frames: 1,
            frame_dt: 1.0,
            profile: scene::DiskProfile::Thin,
            render_mode: RenderMode::Color,
        }
    }
}

pub const USAGE: &str = "\
schwarzschild-raytracer [OPTIONS]
  --width N          image width in px          (default 1920)
  --height N         image height in px         (default 1080)
  --output PATH      output PNG path            (default render.png)
  --r-cam R          camera radius in M         (default 30.0)
  --inclination DEG  polar angle from +z axis   (default 80.0; 90 = in disk plane)
  --azimuth DEG      camera azimuth around +z   (default 0.0)
  --fov DEG          horizontal field of view   (default 75.0)
  --samples N        NxN supersampling          (default 2)
  --serial           render single-threaded (benchmark comparison)
  --max-steps N      per-ray integration cap    (default 60000)
  --step-scale Q     RK4 step h = Q*(r-2)       (default 0.02)
  --exposure X       disk intensity scale       (default 4.0)
  --time T           frame coordinate time in M (default 0.0)
  --spot-amp A       hot-spot amplitude, 0=off  (default 0.0)
  --spot-r R         hot-spot orbit radius in M (default 7.0)
  --frames N         frames to emit (fixed cam) (default 1)
  --frame-dt DT      time step between frames   (default 1.0)
  --profile P        disk temperature law: thin | nt (default thin)
  --render-mode M    color | echo (light-delay false color; default color)";

impl Config {
    pub fn parse(args: impl Iterator<Item = String>) -> Result<Self, String> {
        let mut cfg = Self::default();
        let mut args = args.peekable();
        while let Some(flag) = args.next() {
            let mut value = |name: &str| {
                args.next()
                    .ok_or_else(|| format!("{name} requires a value"))
            };
            match flag.as_str() {
                "--width" => cfg.width = parse_num(&value("--width")?)?,
                "--height" => cfg.height = parse_num(&value("--height")?)?,
                "--output" => cfg.output = value("--output")?,
                "--r-cam" => cfg.r_cam = parse_num(&value("--r-cam")?)?,
                "--inclination" => cfg.inclination_deg = parse_num(&value("--inclination")?)?,
                "--azimuth" => cfg.azimuth_deg = parse_num(&value("--azimuth")?)?,
                "--fov" => cfg.fov_deg = parse_num(&value("--fov")?)?,
                "--samples" => cfg.samples = parse_num(&value("--samples")?)?,
                "--serial" => cfg.serial = true,
                "--max-steps" => cfg.max_steps = parse_num(&value("--max-steps")?)?,
                "--step-scale" => cfg.step_scale = parse_num(&value("--step-scale")?)?,
                "--exposure" => cfg.exposure = parse_num(&value("--exposure")?)?,
                "--time" => cfg.time = parse_num(&value("--time")?)?,
                "--spot-amp" => cfg.spot_amp = parse_num(&value("--spot-amp")?)?,
                "--spot-r" => cfg.spot_r = parse_num(&value("--spot-r")?)?,
                "--frames" => cfg.frames = parse_num(&value("--frames")?)?,
                "--frame-dt" => cfg.frame_dt = parse_num(&value("--frame-dt")?)?,
                "--profile" => {
                    cfg.profile = match value("--profile")?.as_str() {
                        "thin" => scene::DiskProfile::Thin,
                        "nt" | "novikov-thorne" => scene::DiskProfile::NovikovThorne,
                        other => return Err(format!("unknown profile: {other} (thin | nt)")),
                    }
                }
                "--render-mode" => {
                    cfg.render_mode = match value("--render-mode")?.as_str() {
                        "color" => RenderMode::Color,
                        "echo" => RenderMode::Echo,
                        other => return Err(format!("unknown render mode: {other} (color | echo)")),
                    }
                }
                "--help" | "-h" => return Err(String::new()),
                other => return Err(format!("unknown flag: {other}")),
            }
        }
        if cfg.r_cam <= 3.0 {
            return Err("--r-cam must be > 3 (outside the photon sphere)".into());
        }
        if cfg.samples == 0 || cfg.width == 0 || cfg.height == 0 {
            return Err("--width, --height, --samples must be nonzero".into());
        }
        if !(cfg.fov_deg > 0.0 && cfg.fov_deg < 180.0) {
            return Err("--fov must be in (0, 180)".into());
        }
        if cfg.frames == 0 {
            return Err("--frames must be nonzero".into());
        }
        if cfg.spot_amp != 0.0 {
            if cfg.spot_amp <= -1.0 {
                return Err("--spot-amp must be > -1 (temperature stays positive)".into());
            }
            if !(cfg.spot_r > 6.0 && cfg.spot_r < 20.0) {
                return Err("--spot-r must lie inside the disk annulus (6, 20)".into());
            }
        }
        Ok(cfg)
    }
}

fn parse_num<T: std::str::FromStr>(s: &str) -> Result<T, String> {
    s.parse().map_err(|_| format!("invalid number: {s}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_defaults_and_flags() {
        let cfg = Config::parse(std::iter::empty()).unwrap();
        assert_eq!(cfg.width, 1920);
        let cfg = Config::parse(
            ["--width", "640", "--serial", "--fov", "60"]
                .iter()
                .map(|s| s.to_string()),
        )
        .unwrap();
        assert_eq!((cfg.width, cfg.serial, cfg.fov_deg), (640, true, 60.0));
        assert!(Config::parse(["--bogus"].iter().map(|s| s.to_string())).is_err());
    }
}
