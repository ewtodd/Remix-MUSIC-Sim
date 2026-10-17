"""Scan a folder of .toml controls into a list of "reactions".

Layout convention (matches `ControlExamples/`):
  <controls_root>/<reaction>/<control>.toml
Each `<reaction>` directory is a reaction type; each `.toml` inside is one
control (e.g. a bulk run, an interactive run, a beam run).

We also expose the curated fields of a chosen control so the frontend can
prefill the form. The "curated" extraction is a best-effort projection of the
usual fields; everything else stays editable via the raw TOML textarea.
"""

from __future__ import annotations

import os
from typing import Any

from . import control as control_mod


def scan_reactions(controls_root: str) -> list[dict[str, Any]]:
    """Return a list of reactions, each with its controls.

    Returns:
        [{"id": "<reaction>", "label": "<reaction>",
          "controls": [{"name": "<control>", "label": "<control>",
                        "path": "<abspath>", "default": True}, ...]}, ...]
    """
    if not os.path.isdir(controls_root):
        return []

    reactions: list[dict[str, Any]] = []
    for entry in sorted(os.listdir(controls_root)):
        full = os.path.join(controls_root, entry)
        if not os.path.isdir(full):
            continue
        # The first .toml in a reaction dir is the "default" control.
        tomls = sorted(
            f
            for f in os.listdir(full)
            if f.endswith(".toml") and not f.startswith(".")
        )
        if not tomls:
            continue
        controls = []
        for i, name in enumerate(tomls):
            controls.append(
                {
                    "name": name,
                    "label": _friendly(name),
                    "path": os.path.join(full, name),
                    "default": i == 0,
                }
            )
        reactions.append(
            {
                "id": entry,
                "label": _friendly(entry),
                "controls": controls,
            }
        )
    return reactions


def read_control(path: str) -> dict[str, Any]:
    """Read a control file, returning raw TOML text + parsed dict + curated fields."""
    with open(path, "rb") as fh:
        text = fh.read().decode("utf-8")
    fields = control_mod.parse_control_text(text)
    curated = extract_curated(fields)
    return {
        "path": path,
        "toml": text,
        "curated": curated,
    }


def extract_curated(fields: dict[str, Any]) -> dict[str, Any]:
    """Project parsed fields into the curated fields the frontend form edits.

    This is a subset of the schema. The frontend uses these to prefill the
    form, and the raw TOML textarea covers the rest.
    """
    curated: dict[str, Any] = {}
    for key in ("gas", "beam", "target", "detector", "physics", "windows"):
        if key in fields:
            curated[key] = fields[key]
    if "run" in fields:
        run = dict(fields["run"])
        curated["run"] = run
    if "reaction" in fields and "step" in fields["reaction"]:
        # Copy the array-of-tables verbatim.
        curated["reaction_steps"] = fields["reaction"]["step"]
    return curated


def _friendly(name: str) -> str:
    """Turn a file/dir name into a human label: '37Cl_alpha_n' -> '37Cl alpha n'."""
    base = name
    if base.endswith(".toml"):
        base = base[: -len(".toml")]
    base = base.replace("_", " ")
    return " ".join(base.split())
