# Schwarzschild Black Hole Raytracer — Project Spec

## Overview

A CPU-parallel raytracer, written in Rust, that renders a Schwarzschild black
hole by numerically integrating null geodesics through curved spacetime —
not a shader trick, actual general relativity. Output: a rendered PNG showing
gravitational lensing, a warped accretion disk, the photon sphere, and the
event horizon "shadow."

**Goal for the weekend:** a working, physically correct renderer with a
lensed accretion disk and starfield background. Stretch goals below if time
allows.

---

## Physics Background (so Claude Code and I stay on the same page)

We work in **Schwarzschild coordinates** `(t, r, θ, φ)` with the metric
(using geometric units, `G = c = 1`, and setting `M = 1` for the black hole
mass so all distances are in units of `M`):

```
ds² = -(1 - 2/r) dt² + (1 - 2/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²
```

Key radii (in units of M):
- Event horizon: `r = 2`
- Photon sphere (unstable circular photon orbit): `r = 3`
- ISCO (innermost stable circular orbit, useful for accretion disk inner edge): `r = 6`

### Null geodesics

Light rays follow null geodesics. Because the Schwarzschild metric is static
and spherically symmetric, motion is confined to a plane (we can always
orient coordinates so `θ = π/2` for a given ray), which reduces the problem
to 2D: `(r, φ)`.

Two conserved quantities along a geodesic:
- Energy: `E = (1 - 2/r) dt/dλ`
- Angular momentum: `L = r² dφ/dλ`

Define the impact parameter `b = L/E`. The radial "orbit equation" governing
photon trajectories is:

```
(dr/dφ)² = r⁴/b² - r²(1 - 2/r)
```

For ray tracing it's more numerically stable to integrate in an affine
parameter λ using the full second-order geodesic equations (via Christoffel
symbols) rather than the reduced orbit equation directly — this avoids
coordinate singularities and handles rays that plunge into the horizon
cleanly. Use a first-order system of 4 ODEs per ray:

```
d²r/dλ²     = -(Γ^r_tt)(dt/dλ)² - (Γ^r_rr)(dr/dλ)² - (Γ^r_φφ)(dφ/dλ)²
d²φ/dλ²     = -2(Γ^φ_rφ)(dr/dλ)(dφ/dλ)
```

with the relevant nonzero Christoffel symbols for the equatorial slice:

```
Γ^r_tt  = (1 - 2/r)/r²
Γ^r_rr  = -1/(r(r - 2))
Γ^r_φφ  = -(r - 2)
Γ^φ_rφ  = 1/r
```

Integrate this system with **RK4**, stepping λ backward from the camera
(backward ray tracing — start at the eye, trace where each pixel's ray
*came from*).

### Termination conditions per ray

- `r < 2 + ε` → ray fell into the horizon → pixel is black
- Ray crosses the accretion disk plane (`θ = π/2`) within `r ∈ [6, 20]` (ISCO
  to outer disk edge) → sample disk color/temperature at that `(r, φ)`
- `r > R_max` (e.g. 50) and still receding → ray escaped to infinity → sample
  background starfield using the ray's asymptotic direction

---

## Architecture

```
src/
  main.rs           - CLI entry, image setup, output
  metric.rs         - Christoffel symbols, geodesic RHS function
  integrator.rs      - RK4 stepper, adaptive step size near horizon
  camera.rs          - Pixel -> initial ray (position, direction) in Schwarzschild coords
  scene.rs           - Disk intersection + shading, starfield sampling
  render.rs           - Parallel loop over pixels (rayon), assembles image buffer
  color.rs           - Blackbody temperature -> RGB, redshift/Doppler adjustment
```

### Dependencies

- `rayon` — parallel iterator over pixels
- `image` — PNG output
- `nalgebra` (optional) — vector/matrix convenience, not required
- `glam` (optional alternative to nalgebra, lighter weight)

No async, no GPU, no network — deliberately simple dependency footprint so
the physics/numerics is the whole story.

---

## Milestones

### 1. Metric + integrator (get this exactly right first)
- Implement Christoffel symbols and the geodesic RHS as a pure function
- Implement RK4 integrator with fixed step size
- **Sanity test:** a photon fired exactly at the photon sphere radius
  (`r = 3`, tangential) should orbit near-circularly for many steps before
  eventually escaping or falling in (numerically unstable equilibrium —
  it *should* drift, that's physically correct)
- **Sanity test:** a photon with large impact parameter (`b >> 3√3`) should
  pass by with negligible deflection, converging to the flat-space straight
  line
- **Sanity test:** a photon with `b` near the critical value `b_crit = 3√3 M`
  should show strong deflection (this is the lensing "ring" boundary)

### 2. Camera + basic scene
- Place camera at some `r_cam` (e.g. 30M), pointed at the black hole
- For each pixel, compute initial `(r, φ)` and `(dr/dλ, dφ/dλ)` from pixel
  screen coordinates (this is the fiddly part — get the local orthonormal
  tetrad at the camera right so "straight ahead" maps to the correct initial
  ray direction)
- Render just the horizon shadow against a flat gray background first, to
  confirm the silhouette size and shape (should be ~5.2M angular radius,
  bigger than the horizon itself due to lensing) matches known results

### 3. Accretion disk
- Add disk intersection check during integration (crossing `θ=π/2` within
  `[6, 20] M`)
- Simple shading first: solid color or radius-based gradient
- This alone should already produce the iconic "warped ring" image — the
  disk appears both in front of and behind/above the black hole because
  lensed light from the far side of the disk bends around and reaches the
  camera

### 4. Parallelize
- Swap the pixel loop to `rayon`'s `par_iter`, confirm identical output to
  serial version (determinism check), measure speedup

### 5. Polish / stretch goals (pick what's fun)
- Blackbody temperature gradient on the disk (hotter near ISCO, ~10,000K+,
  cooler further out) mapped to RGB
- Relativistic Doppler beaming — disk material orbiting near light-speed on
  the side moving toward the camera should appear brighter/blueshifted, the
  receding side dimmer/redshifted
- Gravitational redshift near the horizon
- Background starfield (procedural or a real star catalog) so lensing of
  background stars is visible as an Einstein-ring-like distortion around the
  shadow
- Adaptive step size (smaller steps near the photon sphere where curvature
  is extreme) for cleaner images with fewer artifacts
- Multiple camera positions / an orbit animation (frames -> video via ffmpeg)

---

## Numerical Notes / Gotchas

- Coordinate singularity at `r = 2` (event horizon) is a coordinate
  artifact, not physical, but it *will* cause `Γ^r_rr` to blow up
  numerically — always check `r < 2 + ε` and terminate before the
  integrator steps into the singularity.
- Step size matters a lot near the photon sphere (`r ≈ 3`); use a smaller
  fixed step there or implement adaptive RK4/RKF45 if renders show ringing
  artifacts around the shadow edge.
- Validate against the known critical impact parameter `b_crit = 3√3 M ≈
  5.196 M` — this is the sharp edge of the black hole's shadow. If your
  shadow radius doesn't converge to this, something's off in the camera or
  integrator.
- Keep `M = 1` throughout and treat all output distances as "in units of M"
  — avoids floating point scale issues and matches how this is discussed in
  the literature, so it's easy to cross-check against papers/textbooks
  (Chandrasekhar's *Mathematical Theory of Black Holes*, or the visualization
  papers behind *Interstellar*'s Gargantua, are good references if you want
  to compare).

---

## Definition of Done (for the weekend)

- [ ] Renders a PNG showing a lensed accretion disk and correctly-sized
      shadow against a background
- [ ] Physics sanity tests above pass (documented in a `tests/` module or at
      least verified manually and noted in README)
- [ ] Parallelized with `rayon`, with a noted speedup vs. serial
- [ ] README with a couple of example renders and a short explanation of
      what's physically happening (good for the portfolio — this is the part
      that makes it "impressive" rather than just "a cool image")