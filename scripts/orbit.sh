#!/usr/bin/env bash
# Render an orbit animation: sweep --azimuth over a full revolution and
# assemble the frames with ffmpeg. All knobs via environment variables:
#
#   FRAMES=720 WIDTH=1280 HEIGHT=720 SAMPLES=2 INCLINATION=80 ./scripts/orbit.sh
#
# Frames land in frames/frame_%04d.png; the video in orbit.mp4. Keep WIDTH
# and HEIGHT even — yuv420p requires it. Extra CLI flags for the renderer
# (e.g. --spot-amp 0.6 --time 0) can be passed as script arguments.

set -euo pipefail
cd "$(dirname "$0")/.."

FRAMES=${FRAMES:-720}
WIDTH=${WIDTH:-1280}
HEIGHT=${HEIGHT:-720}
SAMPLES=${SAMPLES:-2}
INCLINATION=${INCLINATION:-80}
FPS=${FPS:-24}
OUT=${OUT:-orbit.mp4}

cargo build --release
mkdir -p frames

for ((f = 0; f < FRAMES; f++)); do
    az=$(awk -v f="$f" -v n="$FRAMES" 'BEGIN { printf "%.6f", 360.0 * f / n }')
    printf 'frame %d/%d azimuth=%s\n' "$((f + 1))" "$FRAMES" "$az"
    ./target/release/schwarzschild-raytracer \
        --width "$WIDTH" --height "$HEIGHT" --samples "$SAMPLES" \
        --inclination "$INCLINATION" --azimuth "$az" \
        --output "$(printf 'frames/frame_%04d.png' "$f")" "$@"
done

ffmpeg -y -framerate "$FPS" -i frames/frame_%04d.png \
    -c:v libx264 -pix_fmt yuv420p "$OUT"
echo "wrote $OUT"
