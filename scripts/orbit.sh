#!/usr/bin/env bash
# Render an orbit animation: sweep --azimuth over a full revolution and
# assemble the frames with ffmpeg. All knobs via environment variables:
#
#   FRAMES=720 WIDTH=1280 HEIGHT=720 SAMPLES=2 INCLINATION=80 ./scripts/orbit.sh
#
# Frames land in frames/frame_%04d.png; the video in orbit.mp4. Keep WIDTH
# and HEIGHT even — yuv420p requires it. Extra CLI flags for the renderer
# can be passed as script arguments.
#
# TIME_PER_FRAME > 0 advances coordinate time --time per frame (starting at
# TIME_START), so the scene evolves while the camera orbits — combine with
# a hot spot for the full simulation:
#
#   TIME_PER_FRAME=1 FRAMES=360 ./scripts/orbit.sh --spot-amp 1.2

set -euo pipefail
cd "$(dirname "$0")/.."

FRAMES=${FRAMES:-720}
WIDTH=${WIDTH:-1280}
HEIGHT=${HEIGHT:-720}
SAMPLES=${SAMPLES:-2}
INCLINATION=${INCLINATION:-80}
FPS=${FPS:-24}
OUT=${OUT:-orbit.mp4}
TIME_START=${TIME_START:-0}
TIME_PER_FRAME=${TIME_PER_FRAME:-0}

cargo build --release
mkdir -p frames

for ((f = 0; f < FRAMES; f++)); do
    az=$(awk -v f="$f" -v n="$FRAMES" 'BEGIN { printf "%.6f", 360.0 * f / n }')
    t=$(awk -v f="$f" -v t0="$TIME_START" -v dt="$TIME_PER_FRAME" \
        'BEGIN { printf "%.6f", t0 + f * dt }')
    printf 'frame %d/%d azimuth=%s time=%s\n' "$((f + 1))" "$FRAMES" "$az" "$t"
    ./target/release/schwarzschild-raytracer \
        --width "$WIDTH" --height "$HEIGHT" --samples "$SAMPLES" \
        --inclination "$INCLINATION" --azimuth "$az" --time "$t" \
        --output "$(printf 'frames/frame_%04d.png' "$f")" "$@"
done

ffmpeg -y -framerate "$FPS" -i frames/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p "$OUT"
echo "wrote $OUT"
