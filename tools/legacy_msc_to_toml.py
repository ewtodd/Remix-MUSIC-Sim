#!/usr/bin/env python3
"""Convert legacy upstream-musicsim ``.msc`` control files to Remix-MUSIC-Sim
``.toml``.

Use this when you want to run an existing configuration from the original ANL
musicsim (https://gitlab.phy.anl.gov/music/sim, ``.msc`` "PARAMETER VALUE
COMMENT" rows) against this fork's TOML loader.

Usage:
    tools/legacy_msc_to_toml.py path/to/file.msc                # → path/to/file.toml
    tools/legacy_msc_to_toml.py path/to/file.msc -o out.toml    # → out.toml
    tools/legacy_msc_to_toml.py /path/to/musicsim/Examples/**/*.msc

Caveats — the upstream schema doesn't map 1:1 onto Remix's:

  * ``Kb`` was the beam KE *at the gas surface*. Remix's ``beam.energy`` is the
    KE *at the accelerator* (before the entrance window / optional degrader).
    The converter writes ``Kb`` into ``beam.energy`` and emits a TOML comment
    noting the unit change — adjust by hand if you want to model the windows.
  * SRIM-table keys (``SRIMdir``, ``SRIMbeam``, ``SRIM{evap,res}N``) and the
    ``AnodeGeom`` key are silently dropped: catima handles dE/dx and the anode
    geometry is hardcoded in this fork.
  * Gas composition, entrance/exit windows, degrader, and ``Threads`` aren't
    in the upstream schema; their TOML sections are simply omitted (Remix
    uses its compiled-in defaults).

Per-key comments from the .msc file are preserved as trailing TOML comments.
"""

from __future__ import annotations
import argparse
import glob
import sys
from collections import OrderedDict
from pathlib import Path

# section, toml-key, kind  ;  kind in {"int","float","str"}
SCALAR_MAP = {
    "Gas":              ("gas",      "species",     "str"),
    "Pressure":         ("gas",      "pressure",    "float"),
    "pressure":         ("gas",      "pressure",    "float"),
    "Temperature":      ("gas",      "temperature", "float"),
    "temperature":      ("gas",      "temperature", "float"),

    "beam":             ("beam",     "species",     "str"),
    "BeamEnergy":       ("beam",     "energy",      "float"),
    "KbFWHM":           ("beam",     "energy_fwhm", "float"),
    "BeamEnergyFWHM":   ("beam",     "energy_fwhm", "float"),
    "EbeamFWHM":        ("beam",     "energy_fwhm", "float"),
    "dEdxScaleBeam":    ("beam",     "dedx_scale",  "float"),
    # Legacy aliases mapped (with a warning comment, see SEMANTIC_NOTES)
    "Kb":               ("beam",     "energy",      "float"),
    "Ebeam":            ("beam",     "energy",      "float"),

    "target":           ("target",   "species",     "str"),
    "compound":         ("target",   "compound",    "str"),

    "ELossBins":        ("detector", "eloss_bins",  "int"),
    "MaxELoss":         ("detector", "max_eloss",   "float"),
    "strip":            ("detector", "strip",       "int"),
    "stripFirst":       ("detector", "strip_first", "int"),
    "stripLast":        ("detector", "strip_last",  "int"),
    "Eres":             ("detector", "eres",        "float"),

    "NEvents":          ("run",      "n_events",    "int"),
    "Threads":          ("run",      "threads",     "int"),
    "Wait":             ("run",      "wait",        "int"),
    "Update":           ("run",      "update",      "int"),
    "MaxTime":          ("run",      "max_time",    "float"),
    "SimStep":          ("run",      "sim_step",    "float"),
    "Method":           ("run",      "method",      "int"),
    "FileName":         ("run",      "output",      "str"),
    "FileOpt":          ("run",      "file_opt",    "str"),
    "PrintOpt":         ("run",      "print_opt",   "int"),
    "reacClass":        ("run",      "reac_class",  "int"),
}

# Inline tables under [windows]: msc keys -> (sub-table-name, field, kind)
WINDOW_MAP = {
    "EntranceMaterial":  ("entrance", "material",  "str"),
    "EntranceThickness": ("entrance", "thickness", "float"),
    "ExitMaterial":      ("exit",     "material",  "str"),
    "ExitThickness":     ("exit",     "thickness", "float"),
    "DegraderMaterial":  ("degrader", "material",  "str"),
    "DegraderLength":    ("degrader", "length",    "float"),
}

# Reaction-step keys: msc base name -> (which-subtable, field, kind)
REACTION_FIELDS = {
    "Name":    ("evap_or_res", "name",       "str"),
    "Color":   ("evap_or_res", "color",      "int"),
    "dEdx":    ("evap_or_res", "dedx_scale", "float"),
}

# Keys we silently drop (handled implicitly or no longer supported)
DROP_KEYS = {
    "NumEvapPart", "AnodeGeom", "SRIMdir", "CSVfile", "SRIMbeam",
    "SRIMevap0", "SRIMevap1", "SRIMevap2", "SRIMevap3",
    "SRIMres0",  "SRIMres1",  "SRIMres2",  "SRIMres3",
    "Control", "PARAMETER", "Parameter",  # spurious header tokens
}

# Keys whose semantics shifted: convert + append an extra warning comment so
# the user has to look at it instead of silently getting wrong physics.
SEMANTIC_NOTES = {
    "Kb":    "upstream Kb was beam KE at the gas surface; Remix beam.energy "
             "is at the accelerator (pre-window). Adjust if you want window losses applied.",
    "Ebeam": "upstream Ebeam was beam KE at the gas surface; Remix beam.energy "
             "is at the accelerator (pre-window). Adjust if you want window losses applied.",
}


def parse_msc(path: Path):
    """Return list of (key, value, comment) tuples preserving file order."""
    rows = []
    with path.open() as f:
        lines = f.readlines()
    # Skip the two-line header (a title line and a column header)
    for raw in lines[2:]:
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split(None, 2)
        if len(parts) < 2:
            continue
        key, val = parts[0], parts[1]
        comment = parts[2].strip() if len(parts) == 3 else ""
        rows.append((key, val, comment))
    return rows


def cast(val: str, kind: str):
    if kind == "int":
        return int(val)
    if kind == "float":
        return float(val)
    return val


def toml_scalar(val, kind: str) -> str:
    if kind == "str":
        # Drop unhelpful sentinel values
        return f'"{val}"'
    if kind == "int":
        return str(int(val))
    if kind == "float":
        # Keep the user's literal but ensure it parses as a TOML float
        try:
            f = float(val)
        except ValueError:
            return f'"{val}"'
        if f == int(f) and "." not in str(val) and "e" not in str(val).lower():
            return f"{int(f)}.0"
        return str(f)
    return str(val)


def to_toml(rows):
    """Translate (key,value,comment) rows into a TOML string."""
    sections = OrderedDict()
    sections.setdefault("gas", [])
    sections.setdefault("beam", [])
    sections.setdefault("target", [])
    sections.setdefault("windows", [])
    sections.setdefault("detector", [])

    window_subs = {"entrance": {}, "exit": {}, "degrader": {}}
    window_comments = {"entrance": {}, "exit": {}, "degrader": {}}

    # reaction_steps[i] = {"evap": {...}, "res": {...}, "comments": {...}}
    reaction_steps = {}

    unknown = []

    for key, val, comment in rows:
        if key in DROP_KEYS:
            continue

        if key in SCALAR_MAP:
            section, tkey, kind = SCALAR_MAP[key]
            full_comment = comment
            if key in SEMANTIC_NOTES:
                note = f"NOTE ({key}→{tkey}): {SEMANTIC_NOTES[key]}"
                full_comment = f"{comment} -- {note}" if comment else note
            sections.setdefault(section, []).append(
                (tkey, toml_scalar(val, kind), full_comment)
            )
            continue

        if key in WINDOW_MAP:
            sub, field, kind = WINDOW_MAP[key]
            window_subs[sub][field] = toml_scalar(val, kind)
            if comment:
                window_comments[sub][field] = comment
            continue

        # Reaction-step keys: evapNName, evapNColor, dEdxScaleEvapN,
        #                    resNName,  resNColor,  dEdxScaleResN
        matched = False
        for base, side in (("evap", "evap"), ("res", "res")):
            for suffix, field, kind in (
                ("Name",  "name",       "str"),
                ("Color", "color",      "int"),
            ):
                prefix = f"{base}"
                if key.startswith(prefix) and key.endswith(suffix):
                    middle = key[len(prefix):-len(suffix)]
                    if middle.isdigit():
                        idx = int(middle)
                        step = reaction_steps.setdefault(idx, {"evap": {}, "res": {}, "c": {}})
                        step[side][field] = toml_scalar(val, kind)
                        if comment:
                            step["c"][f"{side}.{field}"] = comment
                        matched = True
                        break
            if matched:
                break
        if matched:
            continue
        # dEdxScaleEvapN / dEdxScaleResN
        for prefix, side in (("dEdxScaleEvap", "evap"), ("dEdxScaleRes", "res")):
            if key.startswith(prefix):
                middle = key[len(prefix):]
                if middle.isdigit():
                    idx = int(middle)
                    step = reaction_steps.setdefault(idx, {"evap": {}, "res": {}, "c": {}})
                    step[side]["dedx_scale"] = toml_scalar(val, "float")
                    if comment:
                        step["c"][f"{side}.dedx_scale"] = comment
                    matched = True
                    break
        if matched:
            continue

        unknown.append((key, val, comment))

    # Emit
    out = []
    section_order = ["gas", "beam", "target", "windows", "detector", "reaction", "run"]

    for sec in section_order:
        if sec == "windows":
            entries = [(s, window_subs[s]) for s in ("entrance", "exit", "degrader") if window_subs[s]]
            if not entries:
                continue
            out.append(f"[{sec}]")
            for name, fields in entries:
                ordered_keys = ["material", "thickness", "length"]
                pairs = []
                for k in ordered_keys:
                    if k in fields:
                        pairs.append(f"{k} = {fields[k]}")
                comment_bits = []
                for k in ordered_keys:
                    c = window_comments[name].get(k)
                    if c:
                        comment_bits.append(c)
                trailing = f"  # {' / '.join(comment_bits)}" if comment_bits else ""
                out.append(f"{name} = {{ {', '.join(pairs)} }}{trailing}")
            out.append("")
            continue

        if sec == "reaction":
            if not reaction_steps:
                continue
            for idx in sorted(reaction_steps):
                step = reaction_steps[idx]
                if not step["evap"] and not step["res"]:
                    continue
                out.append("[[reaction.step]]")
                for side in ("evap", "res"):
                    fields = step[side]
                    if not fields:
                        continue
                    pairs = []
                    for k in ("name", "color", "dedx_scale"):
                        if k in fields:
                            pairs.append(f"{k} = {fields[k]}")
                    cmt = step["c"].get(f"{side}.name") or step["c"].get(f"{side}.color")
                    trailing = f"  # {cmt}" if cmt else ""
                    out.append(f"{side} = {{ {', '.join(pairs)} }}{trailing}")
                out.append("")
            continue

        rows_for_sec = sections.get(sec) or []
        if not rows_for_sec:
            continue
        out.append(f"[{sec}]")
        # Align: pad keys to longest in section
        keylen = max(len(k) for k, _, _ in rows_for_sec)
        for k, v, c in rows_for_sec:
            line = f"{k.ljust(keylen)} = {v}"
            if c:
                line += f"  # {c}"
            out.append(line)
        out.append("")

    if unknown:
        out.append("# Unrecognized parameters from the .msc file (left as comments for review):")
        for k, v, c in unknown:
            cmt = f"  # {c}" if c else ""
            out.append(f"# {k} = {v}{cmt}")
        out.append("")

    return "\n".join(out).rstrip() + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("inputs", nargs="+", help=".msc files or globs")
    ap.add_argument("-o", "--out", help="Output file (only when a single input is given)")
    ap.add_argument("--stdout", action="store_true", help="Print to stdout instead of writing files")
    ap.add_argument("--keep", action="store_true", help="Keep the .msc file after conversion (default: leave it)")
    args = ap.parse_args()

    # Expand globs the shell didn't
    paths: list[Path] = []
    for inp in args.inputs:
        matches = glob.glob(inp, recursive=True)
        if not matches and Path(inp).exists():
            matches = [inp]
        for m in matches:
            paths.append(Path(m))

    if args.out and len(paths) != 1:
        print("--out only valid with a single input", file=sys.stderr)
        return 2

    for src in paths:
        rows = parse_msc(src)
        toml_text = to_toml(rows)
        if args.stdout:
            print(f"# === {src} ===")
            print(toml_text)
            continue
        dst = Path(args.out) if args.out else src.with_suffix(".toml")
        dst.write_text(toml_text)
        print(f"wrote {dst}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
