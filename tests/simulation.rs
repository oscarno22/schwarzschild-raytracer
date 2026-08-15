//! Gates for the hot-spot simulation: a disabled spot must leave the image
//! untouched, an enabled spot must move with time, and frames mode must be
//! identical to equivalent single-frame invocations.

use schwarzschild_raytracer::render::{CachedRender, render};
use schwarzschild_raytracer::Config;

fn small() -> Config {
    Config {
        width: 96,
        height: 54,
        samples: 1,
        ..Config::default()
    }
}

/// --spot-amp 0 disables the spot entirely: neither --time nor --spot-r may
/// leak into the output.
#[test]
fn spot_off_render_matches_plain() {
    let plain = render(&small());
    let with_knobs = render(&Config {
        spot_amp: 0.0,
        spot_r: 12.0,
        time: 99.0,
        ..small()
    });
    assert_eq!(plain, with_knobs);
}

/// An enabled spot must change the image, and must move as frame time
/// advances (half an orbital period at r = 7 swings it to the far side).
#[test]
fn spot_renders_and_orbits() {
    let base = render(&small());
    let t0 = render(&Config {
        spot_amp: 0.6,
        ..small()
    });
    let t_half = render(&Config {
        spot_amp: 0.6,
        time: 58.0,
        ..small()
    });
    assert_ne!(base, t0, "spot must be visible");
    assert_ne!(t0, t_half, "spot must move with frame time");
}

/// A frame emitted by the trace-once cache must be byte-identical to a
/// one-shot render invoked directly at that frame's time.
#[test]
fn frames_mode_matches_single_frame() {
    let cfg = Config {
        spot_amp: 0.6,
        time: 3.0,
        frame_dt: 20.0,
        frames: 3,
        ..small()
    };
    let cached = CachedRender::new(&cfg);
    for f in 0..cfg.frames {
        let t = cfg.time + f as f64 * cfg.frame_dt;
        let single = render(&Config {
            time: t,
            frames: 1,
            ..cfg.clone()
        });
        assert_eq!(cached.shade(t), single, "frame {f} diverges");
    }
}
