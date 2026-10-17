"""FastAPI app for the Remix-MUSIC-Sim web interface.

Endpoints (from the plan):
  GET  /api/reactions                         list reactions + controls
  GET  /api/control/{reaction}/{name}         a control's toml + curated fields
  POST /api/run                               start a musicsim run
  GET  /api/run/{run_id}                      run status (+ log tail)
  GET  /api/run/{run_id}/result.root          serve the produced .root to JSROOT
  GET  /api/run/{run_id}/events               JSON projection of the MC tree
  GET  /api/health                            liveness

The server is thin: it spawns musicsim and serves files/JSON. Everything is
kept out of process-global mutable state where possible; runs live in a
per-process dict keyed by run_id (fine for a single-user local tool).
"""

from __future__ import annotations

import logging
import os
import traceback
import uuid
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from . import control as control_mod
from . import reactions as reactions_mod
from . import runner as runner_mod

# ---------------------------------------------------------------------------
# Logging. uvicorn logs to stderr; we also mirror the app's own log lines to
# a file (web/backend/server.log) so a 500 traceback is captured even when the
# process is started in the background without a visible stderr.
# ---------------------------------------------------------------------------

_LOG_PATH = Path(__file__).resolve().parent / "server.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(_LOG_PATH),
    ],
)
logger = logging.getLogger("musicsim.web")

# ---------------------------------------------------------------------------
# Configuration (overridable via env)
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[2]  # .../Remix-MUSIC-Sim
CONTROLS_ROOT = Path(os.environ.get("MUSIC_CONTROLS_ROOT", REPO_ROOT / "ControlExamples"))
RUNS_ROOT = Path(os.environ.get("MUSIC_RUNS_ROOT", REPO_ROOT / "web" / "runs"))
FRONTEND_DIST = Path(os.environ.get("MUSIC_FRONTEND_DIST", REPO_ROOT / "web" / "frontend" / "dist"))
MUSICSIM_BIN = runner_mod.find_musicsim()
env = os.environ.copy()

# ---------------------------------------------------------------------------
# App
# ---------------------------------------------------------------------------

app = FastAPI(title="Remix-MUSIC-Sim web", version="0.1.0")
logger.info("Remix-MUSIC-Sim web starting; controls_root=%s runs_root=%s musicsim=%s",
            CONTROLS_ROOT, RUNS_ROOT, MUSICSIM_BIN)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

_runs: dict[str, runner_mod.Run] = {}


# ---------------------------------------------------------------------------
# Global exception handler: log the full traceback (to server.log) so a 500 is
# diagnosable, and return a minimal JSON body to the client.
# ---------------------------------------------------------------------------


@app.exception_handler(Exception)
async def unhandled_exception_handler(request, exc: Exception):
    logger.error("Unhandled exception on %s: %s", request.url.path, exc)
    logger.error(traceback.format_exc())
    return JSONResponse(status_code=500, content={"detail": str(exc)})


# ---------------------------------------------------------------------------
# Request/response models
# ---------------------------------------------------------------------------


class RunRequest(BaseModel):
    """A run request.

    Either `control_path` (an existing control file) OR curated `fields`
    (optionally with a raw `toml` overlay that's merged on top) is provided.
    `normalize` toggles the per-strip beam normalization feature.
    """

    control_path: str | None = None
    fields: dict[str, Any] | None = None
    toml: str | None = None
    normalize: bool = False
    events_per_strip: int | None = None
    output: str | None = None


def _resolve_run_dir(run_id: str) -> Path:
    d = RUNS_ROOT / run_id
    d.mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Endpoints
# ---------------------------------------------------------------------------


@app.get("/api/health")
def health() -> dict[str, Any]:
    return {"status": "ok", "musicsim": MUSICSIM_BIN, "controls_root": str(CONTROLS_ROOT)}


@app.get("/api/reactions")
def list_reactions() -> dict[str, Any]:
    return {"reactions": reactions_mod.scan_reactions(str(CONTROLS_ROOT))}


@app.get("/api/control/{reaction}/{name}")
def get_control(reaction: str, name: str) -> dict[str, Any]:
    # Guard against path traversal: resolve within CONTROLS_ROOT.
    base = (CONTROLS_ROOT / reaction).resolve()
    path = (base / name).resolve()
    if not path.is_file() or base not in path.parents:
        raise HTTPException(404, "control not found")
    try:
        return reactions_mod.read_control(str(path))
    except Exception as exc:  # invalid TOML etc.
        raise HTTPException(400, f"could not read control: {exc}") from exc


@app.post("/api/run")
def start_run(req: RunRequest) -> dict[str, Any]:
    run_id = uuid.uuid4().hex[:12]
    workdir = _resolve_run_dir(run_id)
    output_name = req.output or "result.root"
    output_path = str(workdir / output_name)
    logger.info("start_run id=%s output=%s normalize=%s fields=%s toml=%s",
                run_id, output_name, req.normalize,
                bool(req.fields), bool(req.toml))

    if req.control_path:
        # Load the existing control file.
        with open(req.control_path, "rb") as fh:
            control_text = fh.read().decode("utf-8")
        control = control_mod.parse_control_text(control_text)
    else:
        # Build from curated fields + optional raw TOML overlay.
        if not req.fields and not req.toml:
            raise HTTPException(400, "provide control_path or fields/toml")
        base = req.fields or {}
        # The output path is required by the schema; inject it before building.
        base.setdefault("run", {})["output"] = output_name
        if req.toml:
            overlay = control_mod.parse_control_text(req.toml)
            merged = control_mod.merge_control(base, overlay)
        else:
            merged = control_mod.build_control(base)
        control = merged
        control_text = control_mod.dump_control(control)

    # Normalize the output path to this run's workdir, and apply run overrides.
    control.setdefault("run", {})
    control["run"]["output"] = output_path
    if req.events_per_strip is not None:
        control["run"]["events_per_strip"] = req.events_per_strip

    # Normalize toggle: if on, write a reference beam control (strip = -1) and
    # mark the run to use it. Full wiring lands in Phase 5; for now we record
    # the flag on the run and leave the beam-reference generation for Phase 5.
    control_text = control_mod.dump_control(control)

    run = runner_mod.run_control(
        run_id,
        workdir,
        control,
        control_text,
        MUSICSIM_BIN,
        env,
    )
    _runs[run_id] = run
    return {"run_id": run_id, "status": run.status}


@app.get("/api/run/{run_id}")
def run_status(run_id: str) -> dict[str, Any]:
    run = _runs.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")
    return {
        "run_id": run_id,
        "status": run.status,
        "phase": run.phase,
        "error": run.error,
        "log_tail": run.log[-4000:],
    }


@app.get("/api/run/{run_id}/result.root")
def run_result(run_id: str) -> FileResponse:
    run = _runs.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")
    out = run.output_path
    if not out.is_file():
        raise HTTPException(404, f"output not found: {out}")
    return FileResponse(str(out), media_type="application/x-root", filename="result.root")


@app.get("/api/run/{run_id}/events")
def run_events(run_id: str) -> dict[str, Any]:
    """Placeholder for the MC-tree JSON projection."""
    run = _runs.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")
    return {
        "run_id": run_id,
        "status": run.status,
        "projection": "pending",  # Phase 3/4 wires the jsroot client projection
    }


# ---------------------------------------------------------------------------
# Serve the built frontend if present
# ---------------------------------------------------------------------------

if FRONTEND_DIST.is_dir():
    app.mount("/", StaticFiles(directory=str(FRONTEND_DIST), html=True), name="static")
