"""Tests for control.py generation against `musicsim --check` as the oracle.

These verify the strict-schema correctness of generated controls. The plan
explicitly uses `musicsim --check` as the authority for what a valid control
shapes like, so the crucial assertion on every generated control is that it
passes `--check`.
"""

from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import pytest

from web.backend import control as control_mod


_REPO_BIN = Path("/home/e-work/Remix-MUSIC-Sim/musicsim")
_FALLBACK = Path("/home/e-work/Remix-MUSIC-Sim/result/bin/musicsim")
# Prefer a freshly `make`-built musicsim (has the latest sim changes).
MUSICSIM = str(_REPO_BIN if _REPO_BIN.exists() else _FALLBACK)


def write_and_check(control: dict) -> tuple[bool, str]:
    """Write a control dict to a temp .toml and run `musicsim --check`."""
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "c.toml"
        p.write_text(control_mod.dump_control(control))
        proc = subprocess.run(
            [MUSICSIM, "--check", str(p)],
            capture_output=True, text=True, cwd=td,
        )
        return proc.returncode == 0, proc.stdout + proc.stderr


# ---------------------------------------------------------------------------
# A valid based control (matches the reference 37Cl_alpha_n_bulk.toml shape)
# ---------------------------------------------------------------------------


def _base_control() -> dict:
    return {
        "gas": {"species": "4He", "pressure": 220.0, "temperature": 293.0},
        "beam": {"species": "37Cl", "energy": 92.0, "energy_fwhm": 0.9},
        "target": {"species": "4He", "compound": "41K"},
        "detector": {"strip_first": 3, "strip_last": 13, "eloss_bins": 555, "max_eloss": 10.0},
        "run": {"events_per_strip": 50, "seed": 12345, "output": "out.root", "threads": 1},
        "reaction_steps": [
            {
                "residue_excitation": "ground",
                "angular_distribution": "isotropic",
                "evap": {"name": "n", "color": 4},
                "res": {"name": "40K", "color": 3, "dedx_scale": 1.0},
            }
        ],
    }


def test_build_reacted_control_passes_check():
    ok, out = write_and_check(control_mod.build_control(_base_control()))
    assert ok, out


def test_build_unreacted_beam_control_passes_check():
    fields = _base_control()
    fields["detector"] = {"strip": -1, "eloss_bins": 555, "max_eloss": 10.0}
    fields.pop("reaction_steps")  # unreacted beam has no steps
    ok, out = write_and_check(control_mod.build_control(fields))
    assert ok, out


def test_build_normalize_control_passes_check():
    """normalize_per_strip + strip_gains are preserved and validate."""
    fields = _base_control()
    fields["detector"]["normalize_per_strip"] = True
    fields["detector"]["strip_gains"] = [1.0] * 18
    fields["detector"]["strip_gains"][3] = 0.434
    control = control_mod.build_control(fields)
    assert control["detector"]["normalize_per_strip"] is True
    assert control["detector"]["strip_gains"][3] == 0.434
    ok, out = write_and_check(control)
    assert ok, out


def test_build_normalize_rejects_bad_gain_length():
    """A non-18-length strip_gains is rejected by the schema (not the builder)."""
    fields = _base_control()
    fields["detector"]["normalize_per_strip"] = True
    fields["detector"]["strip_gains"] = [1.0, 2.0]  # wrong length
    control = control_mod.build_control(fields)
    # The builder passes it through; the schema check is what rejects it.
    assert control["detector"]["strip_gains"] == [1.0, 2.0]
    ok, out = write_and_check(control)
    assert not ok  # musicsim rejects the wrong-length gains


def test_build_defaults_events_and_threads():
    """A run with no events_per_strip/threads defaults to 100 / 4 (web default)."""
    fields = _base_control()
    del fields["run"]["events_per_strip"]
    del fields["run"]["threads"]
    control = control_mod.build_control(fields)
    assert control["run"]["events_per_strip"] == 100
    assert control["run"]["threads"] == 4


def test_build_strip_beam_with_steps_passes_check():
    """strip = -1 (unreacted beam) carrying reaction steps is valid (sim accepts)."""
    # Mirrors ControlExamples/*_beam.toml / *_alpha_alpha.toml that have both.
    fields = _base_control()
    fields["detector"] = {"strip": -1, "eloss_bins": 555, "max_eloss": 10.0}
    # keep the reaction steps (strip == -1 with steps is allowed by the sim)
    control = control_mod.build_control(fields)
    assert control["detector"]["strip"] == -1
    assert "reaction" in control  # steps preserved, not dropped
    ok, out = write_and_check(control)
    assert ok, out


def test_missing_required_field_raises():
    fields = _base_control()
    del fields["gas"]["species"]
    with pytest.raises(ValueError):
        control_mod.build_control(fields)


def test_strip_and_range_conflict_raises():
    fields = _base_control()
    fields["detector"]["strip"] = 3
    with pytest.raises(ValueError):
        control_mod.build_control(fields)


def test_merge_overlay_valid():
    base = control_mod.build_control(_base_control())
    overlay = {"beam": {"energy": 100.0}}
    merged = control_mod.merge_control(base, overlay)
    assert merged["beam"]["energy"] == 100.0
    ok, out = write_and_check(merged)
    assert ok, out


def test_roundtrip_parse_dump_stable():
    text = control_mod.dump_control(control_mod.build_control(_base_control()))
    parsed = control_mod.parse_control_text(text)
    assert parsed["gas"]["species"] == "4He"
    assert parsed["detector"]["strip_first"] == 3
