use std::time::Instant;

use schwarzschild_raytracer::deltatable::{DeltaTable, TableParams};
use schwarzschild_raytracer::render::CachedRender;
use schwarzschild_raytracer::tablerender::{SkySpread, TableFrame};
use schwarzschild_raytracer::{Config, USAGE};

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
        "rendering {}x{} ({}x{} samples/px), camera r={}M inclination={}° azimuth={}° fov={}°{}",
        cfg.width,
        cfg.height,
        cfg.samples,
        cfg.samples,
        cfg.r_cam,
        cfg.inclination_deg,
        cfg.azimuth_deg,
        cfg.fov_deg,
        if cfg.serial { " [serial]" } else { "" },
    );
    let table;
    let cached;
    let shade: Box<dyn Fn(f64) -> Vec<u8>> = if cfg.table {
        let t0 = Instant::now();
        table = DeltaTable::build(cfg.r_cam, TableParams::default());
        let points: usize = table.rows.iter().map(|r| r.phi.len()).sum();
        println!(
            "built δ-table in {:.2?}: {} rows, {points} points, ε_crit = {:.6}",
            t0.elapsed(),
            table.rows.len(),
            table.eps_crit
        );
        let frame = TableFrame::new(&cfg, &table, SkySpread::Magnification);
        Box::new(move |t| frame.shade_rgb(t))
    } else {
        let t0 = Instant::now();
        cached = CachedRender::new(&cfg);
        let elapsed = t0.elapsed();
        let rays = cfg.width as u64 * cfg.height as u64 * (cfg.samples as u64).pow(2);
        println!(
            "traced {rays} rays in {elapsed:.2?} ({:.2} Mrays/s)",
            rays as f64 / elapsed.as_secs_f64() / 1e6
        );
        Box::new(move |t| cached.shade(t))
    };
    if cfg.frames == 1 {
        save(shade(cfg.time), &cfg, &cfg.output);
    } else {
        // Fixed camera: geometry is frame-independent, so emit every frame
        // from the single trace/table above by re-shading only.
        let (stem, ext) = match cfg.output.rsplit_once('.') {
            Some((stem, ext)) => (stem, ext),
            None => (cfg.output.as_str(), "png"),
        };
        let t0 = Instant::now();
        for f in 0..cfg.frames {
            let t = cfg.time + f as f64 * cfg.frame_dt;
            save(shade(t), &cfg, &format!("{stem}_{f:04}.{ext}"));
        }
        println!(
            "shaded {} frames in {:.2?} (dt = {} M)",
            cfg.frames,
            t0.elapsed(),
            cfg.frame_dt
        );
    }
}

fn save(buf: Vec<u8>, cfg: &Config, path: &str) {
    let img = image::RgbImage::from_raw(cfg.width, cfg.height, buf)
        .expect("buffer size matches dimensions");
    if let Err(e) = img.save(path) {
        eprintln!("error: failed to write {path}: {e}");
        std::process::exit(1);
    }
    println!("wrote {path}");
}
