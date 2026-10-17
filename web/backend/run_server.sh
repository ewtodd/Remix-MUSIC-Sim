#!/usr/bin/env bash
# Run the Remix-MUSIC-Sim web backend (dev). Requires the backend venv.
set -e
cd /home/e-work/Remix-MUSIC-Sim
# Prefer a freshly `make`-built musicsim in the repo root (has the latest sim).
REPO_BIN=/home/e-work/Remix-MUSIC-Sim/musicsim
if [ -x "$REPO_BIN" ]; then
  export MUSICSIM_BIN="${MUSICSIM_BIN:-$REPO_BIN}"
else
  export MUSICSIM_BIN="${MUSICSIM_BIN:-/home/e-work/Remix-MUSIC-Sim/result/bin/musicsim}"
fi
exec web/backend/.venv/bin/python -m uvicorn web.backend.app:app \
  --host "${MUSIC_HOST:-127.0.0.1}" --port "${MUSIC_PORT:-8141}"
