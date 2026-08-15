use std::time::Instant;

use schwarzschild_raytracer::{Config, USAGE, render::render};

fn main() {
    let cfg = match Config::parse(std::env::args().skip(1)) {
        Ok(cfg) => cfg,
        Err(msg) => {
            if !msg.is_empty() {
                eprintln!("error: {msg}\n");
            }
            eprintln!("{USAGE}");
            std::process::exit(if msg.is_empty() { 0 } else { 1 });
        }
    };
    println!(
        "rendering {}x{} ({}x{} samples/px), camera r={}M inclination={}° fov={}°{}",
        cfg.width,
        cfg.height,
        cfg.samples,
        cfg.samples,
        cfg.r_cam,
        cfg.inclination_deg,
        cfg.fov_deg,
        if cfg.serial { " [serial]" } else { "" },
    );
    let t0 = Instant::now();
    let buf = render(&cfg);
    let elapsed = t0.elapsed();
    let rays = cfg.width as u64 * cfg.height as u64 * (cfg.samples as u64).pow(2);
    println!(
        "traced {rays} rays in {elapsed:.2?} ({:.2} Mrays/s)",
        rays as f64 / elapsed.as_secs_f64() / 1e6
    );
    let img = image::RgbImage::from_raw(cfg.width, cfg.height, buf)
        .expect("buffer size matches dimensions");
    if let Err(e) = img.save(&cfg.output) {
        eprintln!("error: failed to write {}: {e}", cfg.output);
        std::process::exit(1);
    }
    println!("wrote {}", cfg.output);
}
