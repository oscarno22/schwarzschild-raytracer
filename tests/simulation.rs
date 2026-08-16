//! Gates for the hot-spot simulation: a disabled spot must leave the image
//! untouched, an enabled spot must move with time, and frames mode must be
//! identical to equivalent single-frame invocations.

use schwarzschild_raytracer::render::{CachedRender, render};
use schwarzschild_raytracer::scene::{DiskProfile, NT_PEAK_R, Scene, T_ISCO};
use schwarzschild_raytracer::{Config, RenderMode};

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

/// --profile thin must be byte-identical to the default; Novikov–Thorne
/// must differ, vanish at the ISCO, and peak at r = 49/6 with value T_ISCO.
#[test]
fn novikov_thorne_profile() {
    let plain = render(&small());
    let thin = render(&Config {
        profile: DiskProfile::Thin,
        ..small()
    });
    let nt = render(&Config {
        profile: DiskProfile::NovikovThorne,
        ..small()
    });
    assert_eq!(plain, thin);
    assert_ne!(plain, nt);

    let scene = Scene::new(30.0, 4.0, 1e-3, None, DiskProfile::NovikovThorne);
    assert_eq!(scene.disk_temperature(6.0), 0.0, "NT must vanish at ISCO");
    let t_peak = scene.disk_temperature(NT_PEAK_R);
    assert!((t_peak - T_ISCO).abs() < 1e-6, "NT peak {t_peak} != {T_ISCO}");
    assert!(scene.disk_temperature(NT_PEAK_R - 0.5) < t_peak);
    assert!(scene.disk_temperature(NT_PEAK_R + 0.5) < t_peak);
}

/// Echo mode maps light-travel delay to a colormap whose red channel is
/// strictly increasing: scanning the center column, the topmost disk pixel
/// (far side) must read redder — i.e. older light — than the bottommost
/// (near side).
#[test]
fn echo_mode_far_side_reads_older() {
    let cfg = Config {
        render_mode: RenderMode::Echo,
        ..small()
    };
    let buf = render(&cfg);
    let px = |i: u32, j: u32| {
        let idx = ((j * cfg.width + i) * 3) as usize;
        [buf[idx], buf[idx + 1], buf[idx + 2]]
    };
    let is_disk = |c: [u8; 3]| c != [0, 0, 0] && c != [8, 8, 24];
    let col = cfg.width / 2;
    let top = (0..cfg.height).find(|&j| is_disk(px(col, j))).unwrap();
    let bottom = (0..cfg.height).rev().find(|&j| is_disk(px(col, j))).unwrap();
    assert!(top < cfg.height / 2 && bottom > cfg.height / 2);
    assert!(
        px(col, top)[0] > px(col, bottom)[0],
        "far side {:?} must map later (redder) than near side {:?}",
        px(col, top),
        px(col, bottom)
    );
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
