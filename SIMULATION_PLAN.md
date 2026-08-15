# Simulation Plan — Animation, Orbiting Hot Spot, Real-Time Viewer

Follow-up to the completed static renderer (see README / CLAUDE.md). Goal:
turn the renderer into a running simulation in three phases, in this order.
Phases A and B are the committed scope; C is an optional stretch. Keep the
project ethos: physics first, minimal dependencies, staged commits with
verification at each gate.

Approved sequencing (2026-08-15): orbit animation → hot spot with retarded
time → (optional) δ-table real-time viewer. GPU port explicitly rejected.

---

## Phase A — Orbit animation (offline, ffmpeg)

**A1. `--azimuth DEG` flag** (default 0). Camera position becomes
`r_cam·(sin i·cos α, sin i·sin α, cos i)`. The camera basis
(`forward = −pos/r`, `right = normalize(forward × ẑ)`, `up = right × forward`)
and all downstream math (per-ray plane basis, b_z via e1×e2, starfield in
global directions) are already position-generic — nothing else assumes
azimuth 0. The scene is axisymmetric, so orbiting visibly pans the starfield
and shifts the Doppler asymmetry correctly. Verify: `--azimuth 90` render has
the same disk morphology with rotated stars; determinism test still passes.

**A2. `scripts/orbit.sh`** — shell loop over frames calling the binary with
varying `--azimuth` (and optionally a slow inclination drift), writing
`frames/frame_%04d.png`, then:
`ffmpeg -framerate 24 -i frames/frame_%04d.png -c:v libx264 -pix_fmt yuv420p orbit.mp4`.
Suggested defaults: 720 frames over 360°, 1280×720, `--samples 2`
(~2 s/frame ≈ 25 min total). Keep dimensions even for yuv420p.

Commit A: `animate: camera azimuth + orbit script`.

---

## Phase B — Orbiting hot spot with retarded time (the real simulation)

The scene becomes time-dependent: a bright spot orbits the disk at
r_spot with Ω = r_spot^(−3/2), and every disk hit is shaded at the photon's
**retarded** emission time, so light-travel delays are physically correct.
Payoff: the secondary (wrapped) image of the spot visibly lags the primary,
light echoes sweep the ring, and the spot pulses through its Doppler cycle.

**B1. Integrate coordinate time along rays.** Extend
`metric::State` to `[r, φ, p_r, p_φ, t]` with `dt/dλ = E/(1 − 2/r)` in the
RHS (t feeds back into nothing; RK4 carries it as a 5th component).
Update the `[r, _phi, p_r, p_phi]` destructurings in metric/integrator/tests.
~20% more RHS work — accept it unconditionally rather than keeping two paths.
Sign convention: the trace runs backward with dt/dλ > 0, so emission happens
at `t_emit = −Δt_trace` relative to `t_cam = 0`. Divergence of 1/f near the
horizon is irrelevant: Δt is only consumed at disk hits (r ≥ 6).

- Test: for a radial ray, Δt matches the tortoise-coordinate analytic result
  `Δt = Δr + 2·ln((r₁−2)/(r₂−2))` to ~1e-6 relative.
- Test: Δt to the disk's far side exceeds Δt to the near side.

**B2. Hot spot model.** At the hit, compute the global azimuth
`φ_hit = atan2(x.y, x.x)` from `x = r(cos φ_p·e1 + sin φ_p·e2)` (e1/e2 are
already in scope in `shade_pixel`). Spot center at
`φ_spot(t) = φ₀ + Ω_spot·(τ_frame + t_emit)` (t_emit is negative), wrap with
rem_euclid. In-disk distance
`d² = (r_hit − r_spot)² + (r_hit·wrap(φ_hit − φ_spot))²`; temperature
multiplier `1 + amp·exp(−d²/(2σ_spot²))` so the spot is hotter (bluer) as
well as brighter through the existing T_obs⁴ pipeline. Defaults:
r_spot = 7 (orbital period 2π·7^1.5 ≈ 116 M), σ_spot = 0.7, amp ≈ 0.6.

CLI: `--time T` (frame time in M, default 0), `--spot-amp A` (default 0 =
off; amp 0 must render byte-identically to pre-spot output — test this),
`--spot-r R`. Animation clock is coordinate time; note in README that camera
proper time is the constant factor √(1−2/30) ≈ 0.966 off.

**B3. Fixed-camera animation mode: `--frames N --frame-dt DT`.** For a fixed
camera the ray geometry is identical across frames — only shading changes.
Trace the image once, cache per-sample hit data
(hit kind, r_hit, φ_hit, b_z, t_emit; ~30 B/sample ≈ 120 MB at 720p×2×2),
then emit N frames by re-shading only: ~100× faster than re-tracing, which
makes 30-second spot clips a couple of minutes of work. Orbit + spot combos
fall back to per-frame tracing via the shell script. Suggested frame-dt:
~1 M per frame at 24 fps puts the spot period at ~5 s of video.

Commit B in two: `physics: coordinate-time integration along rays`, then
`simulation: orbiting hot spot with retarded-time shading + frames mode`.

README additions: a clip (GIF or linked mp4), the retarded-time formula, and
a note on what the echo/lag demonstrates.

---

## Phase C (stretch, decide later) — δ-table real-time viewer

Key symmetry: at fixed r_cam the entire trace is a function of the single
screen angle δ. Precompute, for ~4–8k δ values (geometrically refined near
the critical angle), the trajectory polyline (φ_p, r, t) until
escape/capture. Per pixel at runtime: disk crossings solve
`a·cos φ_p + b·sin φ_p = 0` analytically → interpolate r (and t) from the
polyline → shade; sky rays look up the final direction. A few lookups per
pixel instead of ~500 RK4 steps → real-time orbit/look-around on CPU.
Table rebuild only when r_cam changes.

Scope: separate binary (`src/bin/viewer.rs`) behind a cargo feature, since a
window needs one extra dependency (`minifb` or `softbuffer`+`winit`) — a
deliberate, contained exception to the minimal-deps rule. Do not start this
without deciding that tradeoff explicitly.

---

## Order of work / verification gates

| Step | Gate |
|---|---|
| A1 azimuth | tests green; rotated-starfield render sanity check |
| A2 script | orbit.mp4 plays, no seams at the 0/360 wrap |
| B1 time | tortoise-coordinate test passes |
| B2 spot | `--spot-amp 0` byte-identical to base; secondary-image lag visible in a two-frame diff |
| B3 frames mode | frames-mode output identical to equivalent single-frame invocations |
| C viewer | only after an explicit deps decision |
