//! Schwarzschild black hole raytracer: backward ray tracing through null
//! geodesics of the Schwarzschild metric (M = 1 geometric units).

pub mod integrator;
pub mod metric;
pub mod vec3;

/// Render configuration, populated from CLI flags. All distances in units of M.
#[derive(Debug, Clone)]
pub struct Config {
    pub width: u32,
    pub height: u32,
    pub output: String,
    pub r_cam: f64,
    pub inclination_deg: f64,
    pub fov_deg: f64,
    pub samples: u32,
    pub serial: bool,
    pub max_steps: u32,
    pub step_scale: f64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            width: 1920,
            height: 1080,
            output: "render.png".to_string(),
            r_cam: 30.0,
            inclination_deg: 80.0,
            fov_deg: 75.0,
            samples: 2,
            serial: false,
            max_steps: 60_000,
            step_scale: 0.02,
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
  --fov DEG          horizontal field of view   (default 75.0)
  --samples N        NxN supersampling          (default 2)
  --serial           render single-threaded (benchmark comparison)
  --max-steps N      per-ray integration cap    (default 60000)
  --step-scale Q     RK4 step h = Q*(r-2)       (default 0.02)";

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
                "--fov" => cfg.fov_deg = parse_num(&value("--fov")?)?,
                "--samples" => cfg.samples = parse_num(&value("--samples")?)?,
                "--serial" => cfg.serial = true,
                "--max-steps" => cfg.max_steps = parse_num(&value("--max-steps")?)?,
                "--step-scale" => cfg.step_scale = parse_num(&value("--step-scale")?)?,
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
