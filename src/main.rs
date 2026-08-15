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
        "rendering {}x{} ({}x supersampling), camera at r={}M, inclination {}°",
        cfg.width, cfg.height, cfg.samples, cfg.r_cam, cfg.inclination_deg
    );
    // Render pipeline lands in later stages.
}
