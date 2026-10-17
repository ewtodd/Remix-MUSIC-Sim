"""Build and parse Remix-MUSIC-Sim control TOML files.

There are two sources of truth here:

1. The README "Strict TOML control schema" (required vs optional keys).
2. The actual GENERATED control must validate with `musicsim --check`.

The tricky bit is that the loader is strict: unknown keys, wrong types, and
incomplete reaction definitions are errors. So this module only ever emits
keys it knows the loader accepts, and the unit tests assert the result passes
`musicsim --check`.
"""

from __future__ import annotations

import tomllib
from typing import Any

import tomli_w

# ---------------------------------------------------------------------------
# Defaults (from README "Required and optional values"). These are only filled
# in when the caller doesn't supply them; they ensure a valid minimal control.
# ---------------------------------------------------------------------------

_RUN_DEFAULTS: dict[str, Any] = {
    "events_per_strip": 100,  # web-app default; the sim's own README example has 10000
    "threads": 4,
    "wait": 0,
    "update": 0,
    "max_time": 2000.0,
    "sim_step": 0.001,
    "method": 0,
    "file_opt": "recreate",
    "print_opt": 0,
}

_BEAM_DEFAULTS: dict[str, Any] = {
    "energy_fwhm": 0.0,
    "dedx_scale": 1.0,
}

# Default noise: no resolution broadening, no noise. The loader requires the
# "at most one" rule to hold, so we omit noise keys entirely by default rather
# than emitting a disabled -1 table.
_DETECTOR_BASE = {
    "eloss_bins": 500,
    "max_eloss": 25.0,
}


def _normalize_strip(detector: dict[str, Any]) -> None:
    """Enforce the strip selector forms: `strip` OR `strip_first`+`strip_last`."""
    has_single = "strip" in detector
    has_range = "strip_first" in detector or "strip_last" in detector
    if has_single and has_range:
        raise ValueError(
            "detector cannot combine `strip` with `strip_first`/`strip_last`"
        )
    if not has_single and not has_range:
        # Default to an unreacted-beam run; the loader treats strip == -1 as
        # having no reaction.
        detector["strip"] = -1


def build_control(fields: dict[str, Any]) -> dict[str, Any]:
    """Assemble a control dict from curated fields.

    `fields` is a partial/overlay dict using the SAME structure as the TOML
    (nested dicts). Missing optional keys are filled with defaults so the
    result is a valid minimal control. Reaction steps are accepted as
    `reaction_steps` (a list of dicts) or `reaction` -> `step` list.
    """
    out: dict[str, Any] = {}

    # gas: required species/pressure/temperature
    gas = dict(fields.get("gas", {}))
    for required in ("species", "pressure", "temperature"):
        if required not in gas:
            raise ValueError(f"[gas] missing required key '{required}'")
    out["gas"] = gas

    # beam: required species/energy
    beam = dict(fields.get("beam", {}))
    for required in ("species", "energy"):
        if required not in beam:
            raise ValueError(f"[beam] missing required key '{required}'")
    for k, v in _BEAM_DEFAULTS.items():
        beam.setdefault(k, v)
    out["beam"] = beam

    # target: required species/compound
    target = dict(fields.get("target", {}))
    for required in ("species", "compound"):
        if required not in target:
            raise ValueError(f"[target] missing required key '{required}'")
    out["target"] = target

    # windows: optional; pass through if present
    if "windows" in fields:
        out["windows"] = fields["windows"]

    # physics: optional; pass through if present
    if "physics" in fields:
        out["physics"] = fields["physics"]

    # detector: base + strip selector + optional noise
    detector = dict(_DETECTOR_BASE)
    detector.update(fields.get("detector", {}))
    _normalize_strip(detector)
    out["detector"] = detector

    # reaction steps: array of tables. A reacted run (strip >= 0) REQUIRES at
    # least one step. An unreacted beam run (strip == -1) MAY carry steps; the
    # sim accepts them and simply never fires them (beam-only). So preserve the
    # steps either way, matching `musicsim --check` (ControlFile.cpp:888-891).
    steps = fields.get("reaction_steps")
    if steps is None:
        steps = fields.get("reaction", {}).get("step")
    if steps:
        out["reaction"] = {"step": steps}

    # run: required output/seed + defaults (events_per_strip default 100, threads 4)
    run = dict(_RUN_DEFAULTS)
    run.update(fields.get("run", {}))
    for required in ("output", "seed"):
        if required not in run:
            raise ValueError(f"[run] missing required key '{required}'")
    out["run"] = run

    return out


def dump_control(control: dict[str, Any]) -> str:
    """Serialize a control dict to TOML text (tomli_w)."""
    return tomli_w.dumps(control)


def load_control(path: str) -> dict[str, Any]:
    """Load a control TOML file into a dict."""
    with open(path, "rb") as fh:
        return tomllib.load(fh)


def parse_control_text(text: str) -> dict[str, Any]:
    """Parse raw TOML text into a dict (used for the raw editor overlay)."""
    import io

    return tomllib.loads(text)


def merge_control(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    """Deep-merge `overlay` on top of `base`, then re-validate shape.

    Used for the "curated fields + raw TOML overlay" path. The overlay may
    supply any subset; missing required keys still come from `base`.
    """
    merged = _deep_merge(copy(base), copy(overlay))
    return build_control(merged)


def _deep_merge(base: Any, overlay: Any) -> Any:
    import copy

    if isinstance(base, dict) and isinstance(overlay, dict):
        result = copy.deepcopy(base)
        for k, v in overlay.items():
            if k in result and isinstance(result[k], dict) and isinstance(v, dict):
                result[k] = _deep_merge(result[k], v)
            else:
                result[k] = copy.deepcopy(v)
        return result
    return copy.deepcopy(overlay)


def copy(obj: Any) -> Any:
    import copy

    return copy.deepcopy(obj)
