"""Subprocess wrapper around `musicsim` (+ optional cache pre-passes).

The sim resolves paths via compiled-in CATIMA_PREFIX, with TALYS_BIN /
TALYS_ATOMKI_V2_DIR env overrides. Following the plan, `runner.py`:

- writes the generated control into a per-run workdir,
- if the control needs it (per-strip TALYS class tables, or a stopping model
  that requires SRIM), runs the corresponding cache pre-pass first,
- spawns `musicsim <control.toml>` in that workdir (run.output is made
  absolute to the workdir so the tables/ output lands there),
- tracks status + collects a log tail for the UI.

This module is deliberately free of any FastAPI/uvicorn dependency so it can be
unit-tested standalone.
"""

from __future__ import annotations

import os
import subprocess
import threading
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Path resolution
# ---------------------------------------------------------------------------

_MUSICSIM_NAME = "musicsim"
_SRIM_CACHE_NAME = "srim-cache"
_TALYS_CACHE_NAME = "talys-cache"


def find_musicsim(explicit: str | None = None) -> str:
    """Locate the musicsim binary.

    Priority: explicit arg > env MUSICSIM_BIN > `musicsim` on PATH.
    """
    if explicit:
        return explicit
    env = os.environ.get("MUSICSIM_BIN")
    if env:
        return env
    return _MUSICSIM_NAME


def _needs_cache(control: dict[str, Any]) -> dict[str, bool]:
    """Work out whether a control requires the srim-cache / talys-cache pre-pass.

    Conservative, schema-informed heuristic:
    - TALYS per-strip class tables => `[talys]` table present (or a step uses
      residue_excitation = "talys" with level tables). We treat the presence
      of a `[talys]` table as "needs the talys-cache pre-pass".
    - stopping = "srim" (or a low_energy/other stopping model that uses SRIM
      tables) => needs the srim-cache pre-pass.
    """
    needs = {"srim": False, "talys": False}

    physics = control.get("physics", {})
    stopping = physics.get("stopping", "catima")
    if stopping == "srim":
        needs["srim"] = True

    if "talys" in control:
        needs["talys"] = True

    # per-step residue_excitation = "talys" uses level tables; that path needs
    # the class tables only if `[talys]` is present, but flag it conservatively.
    for step in control.get("reaction", {}).get("step", []):
        if step.get("residue_excitation") == "talys":
            needs["talys"] = True
    return needs


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------


class Run:
    """State of a single run: path, subprocess, status, log."""

    def __init__(self, run_id: str, workdir: Path, control_path: Path,
                 musicsim_bin: str, env: dict[str, str] | None = None,
                 output_path: Path | None = None):
        self.run_id = run_id
        self.workdir = workdir
        self.control_path = control_path
        self.musicsim_bin = musicsim_bin
        self.env = env or {}
        self.output_path = output_path or (workdir / "result.root")
        self.status = "pending"  # pending -> running -> done | error
        self.phase = "spawning"  # spawning | cache | running | done
        self.log = ""
        self.log_path = workdir / "run.log"
        self.proc: subprocess.Popen[str] | None = None
        self.lock = threading.Lock()
        self.error = ""

    def _update(self, **kw: Any) -> None:
        with self.lock:
            for k, v in kw.items():
                setattr(self, k, v)

    def start_cache(self, cache: str) -> None:
        """Run a cache pre-pass in the workdir."""
        self._update(status="running", phase="cache")
        if cache == "srim":
            cached = [self._find_bin(_SRIM_CACHE_NAME)]
        else:
            cached = [self._find_bin(_TALYS_CACHE_NAME)]
        cached.append(str(self.control_path))
        proc = subprocess.Popen(
            cached,
            cwd=str(self.workdir),
            env=self._env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        tail, _ = proc.communicate()
        self._append_log(tail)
        if proc.returncode != 0:
            self._update(status="error", error=f"{cache}-cache pre-pass failed", phase="cache")

    def start_musicsim(self) -> None:
        """Spawn musicsim in the workdir (must be called after cache pre-pass)."""
        self._update(status="running", phase="running")
        cmd = [self.musicsim_bin, str(self.control_path)]
        proc = subprocess.Popen(
            cmd,
            cwd=str(self.workdir),
            env=self._env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        self._update(proc=proc)

        def _pump() -> None:
            assert proc.stdout is not None
            for line in proc.stdout:
                self._append_log(line)
            rc = proc.wait()
            if rc == 0:
                self._update(status="done", phase="done")
            else:
                self._update(status="error", error=f"musicsim exited {rc}", phase="done")

        t = threading.Thread(target=_pump, daemon=True)
        t.start()

    def validate(self) -> tuple[bool, str]:
        """Run `musicsim --check control.toml`. Returns (ok, output)."""
        check = [self.musicsim_bin, "--check", str(self.control_path)]
        proc = subprocess.run(
            check,
            cwd=str(self.workdir),
            env=self._env(),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        return proc.returncode == 0, proc.stdout

    def _append_log(self, chunk: str) -> None:
        with self.lock:
            self.log += chunk
            if len(self.log) > 200_000:
                self.log = self.log[-200_000:]
        with open(self.log_path, "a") as fh:
            fh.write(chunk)

    def _env(self) -> dict[str, str]:
        env = dict(os.environ)
        env.update(self.env)
        return env

    def _find_bin(self, name: str) -> str:
        env = os.environ.get(name.upper().replace("-", "_") + "_BIN")
        if env:
            return env
        return name


def run_control(
    run_id: str,
    workdir: Path,
    control: dict[str, Any],
    control_text: str,
    musicsim_bin: str,
    env: dict[str, str] | None = None,
) -> Run:
    """Write the control to the workdir and kick off a run.

    Returns a Run object; callers poll `run.status`. The control's run.output
    is rewritten to be absolute inside the workdir.
    """
    workdir.mkdir(parents=True, exist_ok=True)
    control_path = workdir / "control.toml"
    control_path.write_text(control_text)

    output_path = Path(control.get("run", {}).get("output", workdir / "result.root"))
    run = Run(run_id, workdir, control_path, musicsim_bin, env, output_path=output_path)
    needs = _needs_cache(control)
    if needs["srim"]:
        run.start_cache("srim")
    if run.status == "error":
        return run
    if needs["talys"]:
        run.start_cache("talys")
    if run.status == "error":
        return run
    run.start_musicsim()
    return run
