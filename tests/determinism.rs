//! The parallel render must be byte-identical to the serial one.

use schwarzschild_raytracer::{Config, render::render};

#[test]
fn parallel_render_matches_serial() {
    let base = Config {
        width: 64,
        height: 64,
        samples: 1,
        ..Config::default()
    };
    let serial = render(&Config {
        serial: true,
        ..base.clone()
    });
    let parallel = render(&Config {
        serial: false,
        ..base
    });
    assert_eq!(serial, parallel);
}
