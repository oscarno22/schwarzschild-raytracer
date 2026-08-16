//! Real-time interactive viewer: δ-table rendering into a softbuffer
//! window. Build with `cargo run --release --features viewer --bin viewer`
//! (accepts the same CLI flags as the offline renderer).
//!
//! Controls:
//!   ← / →      orbit azimuth
//!   ↑ / ↓      inclination (toward / away from the pole)
//!   Z / X      zoom in / out (FOV, clamped inside the table's coverage)
//!   Q / E      move camera in / out (rebuilds the δ-table)
//!   Space      pause / resume the simulation clock
//!   , / .      scrub time backward / forward 5 M
//!   H          toggle the hot spot
//!   Esc        quit

use std::num::NonZeroU32;
use std::rc::Rc;
use std::time::Instant;

use schwarzschild_raytracer::deltatable::{DeltaTable, TableParams};
use schwarzschild_raytracer::tablerender::{SkySpread, TableFrame};
use schwarzschild_raytracer::{Config, USAGE};
use winit::application::ApplicationHandler;
use winit::dpi::LogicalSize;
use winit::event::{ElementState, KeyEvent, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{Key, NamedKey};
use winit::window::{Window, WindowId};

/// Simulation clock speed: coordinate time per wall-clock second. One spot
/// orbit at r = 7 (period ≈ 116 M) takes ~29 s at this rate.
const TIME_SCALE: f64 = 4.0;

struct App {
    cfg: Config,
    table: DeltaTable,
    t_sim: f64,
    paused: bool,
    window: Option<Rc<Window>>,
    surface: Option<softbuffer::Surface<Rc<Window>, Rc<Window>>>,
    last_frame: Instant,
    fps_window: (Instant, u32),
}

impl App {
    fn new(mut cfg: Config) -> Self {
        // The offline default (1920×1080×2×2 samples) is far too heavy for
        // a live window; start at 960×540×1 unless the user asked otherwise.
        if (cfg.width, cfg.height) == (1920, 1080) {
            (cfg.width, cfg.height) = (960, 540);
        }
        cfg.samples = 1;
        let t_sim = cfg.time;
        let table = build_table(cfg.r_cam);
        Self {
            cfg,
            table,
            t_sim,
            paused: false,
            window: None,
            surface: None,
            last_frame: Instant::now(),
            fps_window: (Instant::now(), 0),
        }
    }

    fn handle_key(&mut self, event_loop: &ActiveEventLoop, key: &Key) {
        match key {
            Key::Named(NamedKey::Escape) => event_loop.exit(),
            Key::Named(NamedKey::ArrowLeft) => self.cfg.azimuth_deg -= 2.0,
            Key::Named(NamedKey::ArrowRight) => self.cfg.azimuth_deg += 2.0,
            Key::Named(NamedKey::ArrowUp) => {
                self.cfg.inclination_deg = (self.cfg.inclination_deg - 2.0).max(1.0);
            }
            Key::Named(NamedKey::ArrowDown) => {
                self.cfg.inclination_deg = (self.cfg.inclination_deg + 2.0).min(179.0);
            }
            Key::Named(NamedKey::Space) => self.paused = !self.paused,
            Key::Character(c) => match c.as_str() {
                // FOV stays inside the table's angular coverage, so zooming
                // never needs a rebuild.
                "z" => self.cfg.fov_deg = (self.cfg.fov_deg - 5.0).max(25.0),
                "x" => self.cfg.fov_deg = (self.cfg.fov_deg + 5.0).min(100.0),
                "q" | "e" => {
                    let dr = if c.as_str() == "q" { -2.0 } else { 2.0 };
                    let r = (self.cfg.r_cam + dr).clamp(6.0, 200.0);
                    if r != self.cfg.r_cam {
                        self.cfg.r_cam = r;
                        println!("rebuilding δ-table for r_cam = {r}");
                        self.table = build_table(r);
                    }
                }
                "," => self.t_sim -= 5.0,
                "." => self.t_sim += 5.0,
                "h" => {
                    self.cfg.spot_amp = if self.cfg.spot_amp == 0.0 { 0.6 } else { 0.0 };
                }
                _ => {}
            },
            _ => {}
        }
    }

    fn redraw(&mut self) {
        let (Some(window), Some(surface)) = (&self.window, &mut self.surface) else {
            return;
        };
        let size = window.inner_size();
        let (Some(w), Some(h)) = (NonZeroU32::new(size.width), NonZeroU32::new(size.height))
        else {
            return;
        };
        surface.resize(w, h).expect("surface resize");
        let dt = self.last_frame.elapsed().as_secs_f64();
        self.last_frame = Instant::now();
        if !self.paused {
            self.t_sim += TIME_SCALE * dt.min(0.25);
        }
        self.cfg.width = size.width;
        self.cfg.height = size.height;
        let frame = TableFrame::new(&self.cfg, &self.table, SkySpread::Magnification);
        let mut buffer = surface.buffer_mut().expect("surface buffer");
        frame.shade_0rgb(self.t_sim, &mut buffer);
        buffer.present().expect("present");

        self.fps_window.1 += 1;
        if self.fps_window.0.elapsed().as_secs_f64() >= 1.0 {
            let fps = self.fps_window.1 as f64 / self.fps_window.0.elapsed().as_secs_f64();
            window.set_title(&format!(
                "Schwarzschild viewer — r={} incl={}° az={}° fov={}° t={:.0}M{} — {fps:.0} fps",
                self.cfg.r_cam,
                self.cfg.inclination_deg,
                self.cfg.azimuth_deg,
                self.cfg.fov_deg,
                self.t_sim,
                if self.paused { " [paused]" } else { "" },
            ));
            self.fps_window = (Instant::now(), 0);
        }
    }
}

fn build_table(r_cam: f64) -> DeltaTable {
    let t0 = Instant::now();
    let table = DeltaTable::build(r_cam, TableParams::default());
    println!("δ-table for r_cam = {r_cam}: {:.2?}", t0.elapsed());
    table
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        let attrs = Window::default_attributes()
            .with_title("Schwarzschild viewer")
            .with_inner_size(LogicalSize::new(self.cfg.width, self.cfg.height));
        let window = Rc::new(event_loop.create_window(attrs).expect("create window"));
        let context = softbuffer::Context::new(window.clone()).expect("softbuffer context");
        self.surface =
            Some(softbuffer::Surface::new(&context, window.clone()).expect("surface"));
        self.window = Some(window);
        self.last_frame = Instant::now();
    }

    fn window_event(
        &mut self,
        event_loop: &ActiveEventLoop,
        _id: WindowId,
        event: WindowEvent,
    ) {
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::KeyboardInput {
                event:
                    KeyEvent {
                        logical_key,
                        state: ElementState::Pressed,
                        ..
                    },
                ..
            } => self.handle_key(event_loop, &logical_key),
            WindowEvent::RedrawRequested => self.redraw(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _event_loop: &ActiveEventLoop) {
        // Continuous animation: request the next frame immediately.
        if let Some(window) = &self.window {
            window.request_redraw();
        }
    }
}

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
    let event_loop = EventLoop::new().expect("event loop (needs a display server)");
    event_loop.set_control_flow(ControlFlow::Poll);
    let mut app = App::new(cfg);
    event_loop.run_app(&mut app).expect("event loop run");
}
