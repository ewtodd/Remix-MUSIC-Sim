# Plan — web app for Remix-MUSIC-Sim, 2026-10-16

> STATUS: WIP. Phases 0–5 are implemented and verified; Phase 6 (polish) and
> Phase 7 (flake devShell) are not yet done. Committed as a work-in-progress
> snapshot so the web app and the sim normalization feature land together.

## Goal

A web app, living in a new `web/` subdir of this repo, that:

- Lets the user pick a reaction type (from a scan of a configurable folder of
  `.toml` controls, defaulting to the bundled `ControlExamples/`).
- Shows clear, nice form controls for the common run parameters (beam
  species/energy, gas/pressure, detector strips, reaction step), with a raw
  TOML editor for everything else.
- Generates the control `.toml` from those inputs, runs `musicsim` on demand,
  then visualizes the result in the browser.
- Visualizes the detector trace (the per-strip dE curve, reconstructed from the
  `events_MeV` tree) AND a 3D path of the particles in the detector,
  reconstructed from the MC-tree stop points (no simulator change).
- Has an **on/off toggle for trace normalization**, and this normalization
  capability is added to the simulator itself (not just the web app).

This plan has **two parts**:

1. **Web app** (tooling, no simulator physics change): new `web/` tree plus a
   small `flake.nix` devShell addition. No physics REVIEW file.
2. **Simulator feature** (schema + output change, *not* physics): add a
   normalization option that controls whether the per-strip trace TGraphs /
   `events_MeV` are emitted normalized to a.u. This touches the strict TOML
   schema, README, and `Output.cpp`. It is not a physics change (it does not
   alter kinematics, excitation, energy-loss, or stopping), so it does not
   need a physics REVIEW file, but because it changes the strict schema and
   output format it should still be reviewed and called out in the commit.

## Decided (ewt, 2026-10-16)

- **Data model:** the web UI auto-generates the `.toml`, runs `musicsim`, and
  the user views the result. The emphasis is a clear, nice control panel.
- **Rendering:** JSROOT in the browser reads the `.root` directly (minimal
  server), using JSROOT's built-in 2D/3D viewers, plus three.js where we want
  a bespoke 3D detector + path view.
- **3D path fidelity:** reconstruct a coarse path from the MC stop points. No
  simulator change. (The full per-step `TrX/TrY/TrZ/TrK` trajectory is only
  collected when `SaveTrajectory=true`, which is false in bulk runs — so it is
  deliberately out of scope for now.) The 3D view should mirror the original
  sim's interactive TEve visualization (see below).
- **Reaction menu:** scan a configurable folder of `.toml` controls, defaulting
  to the bundled `ControlExamples/`.
- **Hosting:** a new `web/` subdir in this repo, run locally through a nix
  devShell.
- **Backend runtime:** Python + FastAPI, added to the nix devShell. It spawns
  `musicsim`, serves the produced `.root`, and exposes a tiny API.
- **Controls:** a curated set of common controls as form fields, with a raw
  TOML editor for everything else.
- **Event view:** both a single-event picker (index slider / random, with
  checkboxes for which products to draw) and a multi-event overlay mode.
- **Cache pre-passes:** run `talys-cache` / `srim-cache` **automatically** when
  the control needs them, and surface a clear message in the UI that it is
  doing so (with progress/errors). Do not silently skip.
- **3D data source:** **JSROOT `treeDraw` + `graph:true`** produces the JSON
  projection of the MC tree (stop points, energies) **client-side** in the
  browser — no `uproot` dependency, no per-event ROOT parsing in the backend.
- **Detector geometry:** **mirror the hardcoded `Geometry.cpp`** 20×2 anode
  mapping in `detectorgeom.js`.
- **Normalization:** **on/off toggle** in the web UI, and it is added to the
  simulator in general (see "Normalization" below). The reference beam run is
  automatic and cheap (a `strip = -1` run with few events).

## What the output has (verified)

Each `musicsim` run writes one `.root` (path from `run.output` — must be an
absolute or repo-relative path; `musicsim` writes tables under the directory
of `run.output`). It contains:

- `events_MeV` tree — experimental-like detector channels (per-column dE vs
  strip, cathode), used for the "trace" 2D plot.
- `MC` tree — truth. Key branches for the 3D path:
  - `requested_strip`, `reaction_strip` (which strip)
  - `vertex_x/y/z` (reaction vertex, cm; `vertex_z=-1000` ⇒ no reaction)
  - `beam_stop_x/y/z`, `beam_stop_strip`, `beam_termination`
  - `evap_stop_x/y/z[n_steps]`, `evap_stop_strip[n_steps]`, `evap_termination[n_steps]`
  - `residue_step`, `residue_stop_x/y/z`, `residue_stop_strip`, `residue_termination`
  - `evap_energy`, `residue_energy`, `theta_cm/phi_cm` per step, `residue_excitation`
  - `reaction_channel` (TALYS class channel, 0 for configurable steps)
- `traces/` directory — per-event TGraphs named `Trace_s<step>_e<evt>` (and
  `Trace_s<step>_c<col>_e<evt>`), x = strip index, y = dE per column. **Only
  written when `run.update = 1`** (the TEve/GL interactive viewer), which fails
  headless — so the 2D view uses the `events_MeV` tree instead (see the Phase 0
  spike results).
- `metadata/` — control text, level/TALYS table paths + MD5s, channel
  cross-sections, tagline/revision info, and per-strip resolved config.

### Detector geometry (for the 3D reconstruction)

Hardcoded in `Geometry.cpp` (from `agfiles/AnodeGeometry`, 2019/11/04):
20 rows × 2 columns. Strip 0 (full width, 9 cm), strips 1–16 split L/R
(halves sum to 9 cm, short half alternates), strip 17 (full width). Rows 0/19
are dead-layer caps. Row pitch (dz) ≈ 1.578 cm per strip; dead caps ≈ 3.5–3.6
cm. So strip index → z position is monotonic and known; x/y come from the
strip's width/half split. This is enough to render the detector as a simple
3D stack of strips and to draw a coarse path.

## Phase 0 spike — results (2026-10-16)

All gates passed. Verified on this machine with the real installed `musicsim`
(`result/bin/musicsim`, v26.10.15):

- **Sim spawn (PASS).** Wrote a generated control into an isolated workdir
  (`web/runs/_spike_sim/`), `musicsim --check` → "Control file is valid",
  then `musicsim spike.toml` → produced `out/spike.root` (550 MC + 550 events
  rows). Works exactly as the backend will call it. **Note:** the reaction
  kinematics in that throwaway control produced mostly unreacted beam (`n_steps
  =1`, `vertex_x=0`, `beam_stop_z=-1000`) — expected, not a code issue.
- **Toolchain (PASS).** node v24 + npm 11 via `nix shell nixpkgs#nodejs`;
  python 3.14 via `nix shell nixpkgs#python3`. The repo's `flake.nix` has no
  node/python yet — that's the Phase 7 devShell `buildInputs` add. (Python
  `pip` needed `python3Packages.pip` and a libstdc++ note; we don't need a
  Python venv — see "3D data" below.)
- **JSROOT reads the `.root` (PASS).** jsroot 7.11.2 `openFile` reads both
  trees and all 42 `MC` branches. Branch values pull cleanly via
  `treeDraw(tree, {expr, graph:true, arr_limit:1e9})` → returns a `TGraph`
  with `fX` (entry id) and `fY` (value). Verified on `reaction_strip`
  (5 distinct), `beam_energy_gas` (549 distinct), `beam_stop_z`, `n_steps`,
  and all `vertex_*`/`*_stop_*`. `events_MeV` exposes `LeftdE[18]`/
  `RightdE[18]` arrays + `Cathode`/`Strip0dE`/`Strip17dE` scalars.
- **three.js renders detector + path (PASS).** three.js imports and renders in
  headless Chromium (WebGL 2.0 via `--enable-unsafe-swiftshader
  --ignore-gpu-blocklist`). The 20×2 anode grid renders as a translucent
  box array; thick track tubes (white beam, cyan evap, red residue) render
  clearly and meet at the vertex. Scene graph = 45 children (40 anodes + 3
  tracks + 2 lights); strip↔z math validated (row0 z=-15.78, `stripZ(8)=-3.16`,
  `stripZ(17)=11.05`).

### Two findings that change the plan

1. **The `traces/` directory is empty unless `run.update = 1`.** In
   `Simulator.cpp`, trace TGraphs are created+written only under `ctf.Update`
   (`if (ROOTfile != 0 && ctf.Update) CreateTracesAndTrajectories();`), and
   `update = 1` brings up the full TEve/GL interactive viewer (`new
   TEveManager(...)`, `TEveArrow`, `TCanvas`) which **fails headless**
   ("Display is not set" / "Insufficient support from the graphics hardware").
   So on a headless server you cannot get the `traces/` TGraphs the way the
   2D view is currently described.
   - **Resolution (decided):** the 2D trace is reconstructed from the **bulk
     `events_MeV` tree** instead — it carries `LeftdE[18]`/`RightdE[18]`
     per-strip dE plus `Cathode`, which is exactly the trace curve. No
     `update=1`, no GL, works headless. (The `traces/` TGraphs remain the
     interactive viewer's internal drawing; the web app does not need them.)
     This also cleanly handles the normalization toggle, since it's applied to
     the `events_MeV` values.
   - This is an *approach* change, not a new simulator feature — the web app
     simply reads the bulk trees it already gets.
2. **The backend does not need Python/ROOT for the 3D data.** JSROOT reads the
   `.root` directly in the browser, so the "JSON projection of the MC tree"
   can be produced **client-side** by `treeDraw` + `graph:true` (as verified).
   This removes the need for an `uproot` dependency in the FastAPI backend.
   The backend stays thin: spawn `musicsim`, serve the `.root`, tiny API.

### Could not verify here
- Real **TALYS/SRIM cache pre-pass** spawning (the throwaway control used a
  simple ground/isotropic step, no cache). The auto-run + message behavior is
  specified in the plan; it will be exercised in Phase 1 against a
  `per-strip` TALYS control.

## Normalization (added to the simulator, toggle in the web app)

### What "normalized" means (from ~/MUSIC's treatment of sim traces)

Sim energies are arbitrary-unit floats. To overlay sim on calibrated data
(which has the beam at ~1 a.u. per strip), the sim traces are **per-strip
normalized to 1 a.u.**, not divided by each trace's own peak.

Reference: `~/MUSIC/tooling/src/StripSumScatter.cpp`:

- `SimBeamGains` (`gain[s] = 1 / mean(beam per-strip total)`, from the
  `strip = -1` beam run) — with unit gains, computes each strip's mean beam
  total, then sets `gain[s] = 1 / mean[s]`. Every strip's sim beam lands on
  1 a.u. This flattens the sim beam the same way the per-channel data
  normalization flattens the experimental beam. A single global factor would
  not — it would preserve the sim's per-strip structure. Strips with no beam
  signal keep gain 0 (drop out like an uncalibrated channel).
- `SimTotal` / `SimPopTraces`: `total[s] = gain[s] * Total(s)` (and the
  IGNORE_SHORT variant applies the gain to the single short half), so the
  traces share the data's axis and the sim beam is flat at the data level.
- `python/data.py max_normalize` and `python/blind_an.py _trace` are the
  Python mirrors: "Flatten per strip by gain … flattened a.u. (beam == 1),
  absolute amplitude preserved (no per-trace peak division)."

So the web-app normalization toggle: **off** = raw sim arbitrary units;
**on** = per-strip normalized (beam → 1 a.u. per strip, via the beam-run
gains). This is NOT max-normalization per trace (that is a separate shape-only
mode).

### Simulator change (the part to add to the sim in general)

Add a control option (proposed: `[detector] normalize_traces = false` or a
`[run.normalize]` block — decide the exact key) that makes the sim emit
per-strip-normalized traces. The sim already knows the beam particle; the
beam-reference gains are computable from a `strip = -1` beam run.

**Decided mechanism (ewt, 2026-10-16):** the reaction control **references an
automatic reference beam run** — a `strip = -1` (unreacted) beam control that
`musicsim` runs first (or reads if its output already exists) to build the
per-strip gains. This reference beam run is cheap: it **barely needs any
events** to get the per-strip mean energy, so it finishes quickly. The web app
generates that beam control automatically when normalization is toggled on.

Concretely:

- When normalization is on, the backend/`control.py` also writes (or points to)
  a `strip = -1` beam control and runs it first with a small
  `events_per_strip` (enough for a stable per-strip mean — a few hundred, not
  the full reaction count).
- `musicsim` computes `gain[s] = 1 / mean(beam per-strip total)` from that
  reference beam output and applies it when emitting normalized traces.
- If the reference beam output already exists (same beam/gas/geometry), reuse
  it rather than re-running — the gain set is unchanged for identical inputs.

This is a schema + `Output.cpp` change (and README), not physics.

Also note: the `events_MeV` tree is labelled "MeV" today; if normalization is
turned on the values become a.u. units, so the metadata (or the branch name /
a flag) should record which was emitted, so downstream consumers never mix
them. Surface this in the README.

## Original sim visualization to mirror (already in the code)

The original fork ships an interactive TEve viewer, toggled by `run.update = 1`
(`Simulator.cpp` `SetupRun` / `Simulator::DrawMUSIC`, `Output.cpp`
`UpdateVisuals`). It draws, for one event:

- **3D detector geometry** via `TEveGeoTopNode` (`DrawMUSIC`), plus axis
  arrows. This is the reference for what the web app's 3D path view should
  reproduce (detector as a stack of strips + the beam/product paths).
- **Per-particle track arrows** (`TrackBeam`, `TrackEvaP`, `TrackEvaR`) for
  beam and each evaporation product/residue.
- **2D trace canvas, two panes:** left = all column traces for the event;
  right = all particle traces (beam, evap, residue) for the event, with
  legends.

The web app's 2D trace view and 3D path view should visually match what this
TEve viewer shows (same strip-vs-dE axes, same per-particle colors/legend),
so there is one familiar picture whether you run `musicsim` interactively or
in the browser. Use the `traces/` TGraphs for the 2D panes (or the original
`TraceUB`/`TraceB`/`TraceER`/`TraceEP` grouping), and the MC stop points for
the 3D paths.

## Architecture

```
web/
  plan/            (this plan)
  backend/
    app.py         FastAPI app
    control.py     .toml build/parse (tomli / tomli-w or a small hand-rolled
                   writer that matches the strict schema)
    runner.py      subprocess wrapper around musicsim (+ srim-cache /
                   talys-cache pre-pass when needed)
    reactions.py   scan a folder of .toml controls -> reaction list
    static/        built frontend (served by FastAPI)
  frontend/
    index.html
    src/
      main.js        app shell, reaction picker
      controls.js    curated form + raw TOML editor + normalize toggle
      run.js         POST /api/run, poll status
      trace2d.js     JSROOT: draw traces/ TGraphs
      path3d.js      three.js: detector + coarse MC path
      detectorgeom.js  strip -> z / x / y mapping (mirror of Geometry.cpp)
      normalize.js   apply per-strip beam gains (on/off)
    package.json / vite.config.js
  README.md          how to run (nix develop, then uvicorn)
```

### API (FastAPI)

- `GET /api/reactions` — list reactions from the scanned controls folder:
  each reaction = dir of `.toml` controls, with a friendly label.
- `GET /api/control/<reaction>/<name>` — return the `.toml` text + parsed
  curated fields for the chosen control (so the UI can prefill).
- `POST /api/run` — body carries either a control path or generated fields
  (curated + optional raw TOML overlay + normalize flag). Server writes a
  `.toml` into a per-run working dir, runs the cache pre-pass when needed,
  spawns `musicsim`, returns a `run_id`.
- `GET /api/run/<run_id>` — status (running / done / error, log tail). Also
  reports whether a cache pre-pass is running (the "message saying so" the UI
  shows).
- `GET /api/run/<run_id>/result.root` — serve the produced `.root` to JSROOT
  (correct content-type; JSROOT fetches via HTTP).
- `GET /api/run/<run_id>/events` — a JSON projection of the MC tree (stop
  points, energies, channel names) for the three.js path. Also returns the
  per-strip beam gains (if normalization is on) so the browser can normalize
  without re-parsing ROOT.

`musicsim` invocation: `musicsim <control.toml>`; `--check` validates only.
Working dir per run isolates the `talys/` / `srim/` table output that
`run.output` implies. Runs go in a scratch dir (e.g. `web/runs/<run_id>/`),
gitignored; `run.output` inside the generated control is made absolute to that
dir.

### Frontend

- Reaction picker: dropdown/list from `/api/reactions`; selecting one loads
  its controls.
- Control panel: curated fields (beam species/energy/FWHM, gas species/
  pressure/temperature, detector `strip_first`/`strip_last` or single `strip`,
  reaction step residue/evap names, `events_per_strip`, seed) rendered as
  inputs, plus a collapsible raw TOML textarea that stays in sync.
- **Normalize toggle** — on/off; when on, the UI shows the traces normalized to
  1 a.u. per strip and (if the sim option is on) the generated control carries
  the normalize flag / beam reference. The toggle is surfaced clearly in the
  panel.
- Run button → POST `/api/run`, poll until done. While a cache pre-pass runs,
  show an explicit "generating TALYS/SRIM tables (one-time)…" message.
- Results view:
  - **Trace (2D):** JSROOT opens the `.root` and draws the selected event's
    `traces/` TGraphs (strip vs dE per column), raw or normalized per the
    toggle.
  - **Path (3D):** three.js renders the detector as a stack of strips
    (semi-transparent boxes colored by strip id) and the coarse path as
    line/points: beam enters along −z, stops or reacts at `vertex_x/y/z`; then
    evap (and any further products) and residue branch out to their stop
    points. Products colored by the control's `color` field (beam blue, evap
    red, residue green as in the control examples) and legend text from the
    control's names — mirroring the original TEve viewer's per-particle
    arrows and legend.
  - Event toggle: single-event index slider + random, with checkboxes for
    beam/evap/residue; or an overlay mode that draws many events (e.g. every
    Nth event, colored by reaction/branch).
- Event index alignment: pick event N in JSROOT's trace view and in the 3D
  view from the same N (both read the same event index in the file).

## Phase 1 results (2026-10-16)

Backend implemented in `web/backend/` and smoke-tested end-to-end against the
real `musicsim`:

- **`control.py`** — `build_control`/`dump_control`/`load_control`/
  `parse_control_text`/`merge_control`. Enforces the strict schema (required
  keys, strip-selector exclusivity, reaction-step shape).
- **`reactions.py`** — scans `ControlExamples/` (9 reactions) into
  reaction→controls; `read_control` returns raw TOML + parsed dict + curated
  fields; `extract_curated` projects the fields the form edits.
- **`runner.py`** — `Run` class + `run_control`; spawns `musicsim` in an
  isolated per-run workdir with `run.output` made absolute to it; detects and
  runs `srim-cache`/`talys-cache` pre-passes when a control needs them; streams
  the log; tracks status / error.
- **`app.py`** — FastAPI endpoints: `/api/health`, `/api/reactions`,
  `/api/control/{reaction}/{name}`, `POST /api/run`, `GET /api/run/{id}`,
  `GET /api/run/{id}/result.root`, `GET /api/run/{id}/events` (placeholder).
- **Tests** — `tests/test_control.py` (6 tests) uses `musicsim --check` as the
  oracle. **All pass**, including reacted + unreacted-beam generation.

Smoke test proved the full loop: `POST /api/run` with generated fields →
`musicsim` ran to `done` in the workdir (30 events/strip) → `result.root`
served (valid ROOT, magic `root`).

Notes / follow-ups for later phases:
- The `/api/run/{id}/events` MC projection is a placeholder; the browser does
  the projection with jsroot (Phase 3/4).
- The normalize toggle currently only records the flag on the run; the beam
  reference run + sim-side normalization land in Phase 5.
- Python deps live in `web/backend/.venv` (built via `nix shell nixpkgs#python3
  nixpkgs#python3Packages.pip` + `pip install`). The flake devShell addition
  (which would make these available without a venv) is Phase 7; until then
  `web/backend/run_server.sh` starts uvicorn from the venv.

## Phase 2 results (2026-10-16)

Frontend implemented in `web/frontend/` (Vite + vanilla JS modules) and
exercised end-to-end in headless Chromium against the live backend:

- **`vite.config.js`** — dev server proxies `/api` → `127.0.0.1:8141`, so no
  CORS. `localhost` favicon 404 is harmless.
- **`src/main.js`** — loads `/api/reactions`, populates the reaction + control
  dropdowns, wires Load-control and Run.
- **`src/controls.js`** — curated form (Gas / Beam / Target / Detector / Run),
  the **normalize toggle**, and a **raw TOML editor** kept in sync via
  `smol-toml` (`src/toml.js`). Editing a curated field re-serializes the raw
  editor; editing raw re-parses; structure is a single `state.control` dict.
- **`src/run.js`** — `POST /api/run` with `fields` + `normalize`, then polls
  `GET /api/run/{id}` every second, showing phase ("running" / "generating
  tables…") and the live log tail. Fires `state.onDone(runId)` when finished.

**Verified in the browser (puppeteer + system Chromium):**
- Reaction dropdown has 9 reactions; 37Cl_alpha_n lists 8 controls.
- Loading `37Cl_alpha_n_bulk.toml` populated the form (`beam.energy=92`,
  `gas.species=4He`) and the raw editor (816 chars).
- Clicking Run produced a real `result.root` (14 MB, 8 worker logs) and the
  status reached `done`; `onDone` fired.
- Screenshot confirmed the intended layout renders.

**Notes / deferred:**
- Editing the curated `strip` vs `strip_first/last` isn't yet mutually
  exclusive in the UI (backend keeps it exclusive via the schema; tested).
- Backend still has no cancel/stop endpoint; `#cancel-btn` is a no-op
  placeholder for a later phase.
- `state.onDone` is a console stub — Phase 3/4 attach the trace + 3D views.

## Phase 3 results (2026-10-16)

2D trace view implemented in `web/frontend/src/` and verified end-to-end in
headless Chromium against a real run.

- **Data source (per the Phase 0 finding):** reconstruct the trace from the
  `events_MeV` tree — `traces/` is only written under `run.update=1`, which
  fails headless. Verified columns: `LeftdE[16]`/`RightdE[16]` (interleaved L/R
  electrodes per strip), `Strip0dE`, `Strip17dE`, `Cathode`. The per-strip trace
  = `Strip0dE` (strip 0) | `Left[s-1]+Right[s-1]` (strips 1–16) | `Strip17dE`
  (strip 17).
- **`src/trace2d.js`** — `loadTraceData(root)` pulls the columns via jsroot
  `treeDraw` + `graph:true` (per-event array flattening) and `drawTrace(ctx,…)`
  renders the strip-dE curve on a canvas (axes, grid, points).
- **`src/normalize.js`** — `computeGains`/`applyGains` implement the ~/MUSIC
  per-strip beam normalization. **Verified numerically:** after applying gains,
  the beam (unreacted, `reaction_strip=-1`) per-strip means are all **1.000**.
- **`src/traceview.js`** — fetches `/api/run/{id}/result.root` as a blob, opens
  it with jsroot, shows the results panel with an **event slider**, the
  **normalize toggle**, and the trace canvas. Wired into `state.onDone`.

**Verified in the browser:** run completed → trace panel appeared with
`n=110000` events (11 strips × 10000), `beam=72970` unreacted events, gains
applied; ~9300 non-background pixels drawn per canvas confirmed the curve
renders; slider + normalize toggle both re-draw. Screenshot confirmed the
blue dE-vs-strip curve, MeV y-axis, and controls.

**Notes / deferred.**
- The beam reference for normalization is derived **from the run's own
  unreacted events** for now. The dedicated `strip=-1` reference run is Phase 5.
- jsroot emits a harmless `ERR_REQUEST_RANGE_NOT_SATISFIABLE` when reading blob
  URLs (requests a range at EOF); data still loads correctly.
- Large files (10000 events/strip → ~14 MB, 110000 rows) load in the browser
  fine, though the traceview chunk is big (code-split via dynamic import).

## Phase 4 results (2026-10-16)

3D path view implemented in `web/frontend/src/` and verified end-to-end in
headless Chromium against a real run.

- **`detectorgeom.js`** — mirrors the hardcoded `Geometry.cpp` 20×2 anode
  mapping: row dz / dx / color, strip→z-centre, column x-centres, dead caps,
  total depth. Read-only JS mirror of the C++.
- **`mcdata.js`** — client-side MC projection via jsroot `treeDraw` + `graph:
  true`; zips 23 branches into per-event records (no uproot/backend).
- **`path3d.js`** — three.js scene: translucent anode stack + edges, RGB axis
  arrows, and track tubes (white beam, cyan evap, pink residue). Coarse path:
  beam `(0,0,0)→vertex` (reacted) or `→beam_stop`; evap/residue `vertex→stop`.
  Single-event and overlay modes, product toggles, OrbitControls.
- **`pathview.js`** — controller: event slider, overlay-all, product checkboxes,
  meta line; wired into `state.onDone`.

Verified: single event shows one path (beam + evap + residue from the vertex);
overlay shows the dense spray of all ~330 events. Track colors and the detector
grid render correctly (confirmed by screenshots). Fixed the renderer-detach bug
(re-append the canvas when the container is cleared on mode change).

- **Also fixed a Phase 2 latent bug:** curated field edits only re-serialized
  the raw TOML when the input event bubbled (`bubbles:true`). Rendered inputs
  fire bubbling events for real typing, but programmatic `dispatchEvent(new
  Event("input"))` without `bubbles:true` didn't — so a test-driven edit
  appeared not to propagate. Confirmed the app behaviour is correct for real
  input; test scripts now use `bubbles:true`.

## Phase 5 results (2026-10-16)

- **Schema (new `[detector]` keys):** `normalize_per_strip` (bool) +
  `strip_gains` (18 values, strip 0..17). Parsed + validated in
  `ControlFile.cpp` (added to the `CheckKeys` whitelist; wrong-length /
  bad-value arrays rejected), mirrored in `include/Simulator.hpp`
  (`controlFileParams::normalizePerStrip` / `stripGains`, default gains 1.0),
  and documented in the README example + required/optional table.
- **Sim applies gains** in `Simulator::FinalizeEvent` (right before
  `SimTree->Fill()`): multiplies `LeftdE`/`RightdE` per segmented strip and
  `Strip0dE`/`Strip17dE` by `stripGains[strip]`. Cathode stays in MeV (no
  per-strip meaning).
- **Web wiring (`web/frontend/src/normalize.js`):** when the toggle is on,
  `run.js` first POSTs a `strip = -1` beam-reference control (same gas/beam,
  small events), polls to done, reads its `events_MeV` with JSROOT, computes
  `gain[s] = 1 / mean(beam per-strip total)`, injects `normalize_per_strip =
  true` + `strip_gains` into the reaction control, then runs it. Client-side
  jsroot keeps the backend thin. Kept the Phase 3 client-side
  `computeGains`/`applyGains` for the non-sim-normalized trace path.
- **Backend:** `run_server.sh` now prefers a freshly `make`-built repo-root
  `musicsim` (the sim change only takes effect when the binary is rebuilt), so
  the web app uses the latest sim.

### Verified
- `musicsim --check` accepts the new keys; rejects wrong-length `strip_gains`
  ("must have exactly 18 values") and non-bool `normalize_per_strip`.
- Direct sim A/B (seed fixed): the injected `strip_gains` give an exact
  **ratio 2.000** on the target strip, **1.000** elsewhere.
- Full browser end-to-end with the toggle on: beam-reference run → gains
  injected → reaction run done. Verified `beamMean × injected gain = 1.0000`
  for all 16 segmented strips (beam lands flat at 1 a.u.), and the trace view
  shows "gains applied".

## Post-Phase-5 fixes (2026-10-16)

Reported by a live test run; all three fixed and re-verified.

- **jsroot "can't access property getUint8, this.arr is null" / 3D error.**
  Root cause: the frontend opened the served `.root` via a **blob URL**, and
  jsroot's lazy range reads hit a revoked/range-limited blob
  (`ERR_REQUEST_RANGE_NOT_SATISFIABLE` / `ERR_FILE_NOT_FOUND`). Fix: new
  `rootio.js#openResultRoot()` fetches the file as an **ArrayBuffer** and calls
  `jsroot openFile(buffer)` directly (no blob, no range requests). All three
  open sites updated: `traceview.js`, `mcdata.js` (3D), `normalize.js`
  (beam-gain read).
- **`POST /api/run` 500 on 9 controls** (the `*_alpha_alpha.toml` / `*_beam.toml`
  family). Root cause: `control.py#build_control` rejected `strip = -1`
  (unreacted beam) combined with a `[[reaction.step]]`, but the sim accepts that
  and simply never fires the steps (beam-only). Fixed to match `musicsim
  --check`: preserve the steps whether or not the strip selector is reacted
  (ControlFile.cpp:888-891).
- **Defaults.** A run with no `events_per_strip` / `threads` now defaults to
  **100 events / 4 threads** (`_RUN_DEFAULTS`), instead of the sim's heavy
  10000/8. `run.seed` stays required.
- **Logging.** Added a `logging` setup to `app.py` that mirrors the app's own
  log lines to `web/backend/server.log`, a startup line, a `start_run` line, and
  a global `@app.exception_handler(Exception)` that writes the full traceback to
  `server.log` and returns `{detail}` JSON so a 500 is diagnosable.

Backend suite now **10 tests** (added: default 100/4, strip=-1-with-steps
passes `--check`).

## Phases

0. **Spike, the gate.** Confirm the React/Vite/JSROOT/three.js toolchain builds
   under a nix devShell on this machine, and that JSROOT can load a real
   `traces_*.root` from `ControlExamples` and draw a trace. Also confirm
   `musicsim` can be spawned from FastAPI with a generated control in an
   isolated working dir (small `events_per_strip` for speed). Stop and report
   if JSROOT can't read the `traces/` TGraphs or three.js assembly is awkward.
1. **Backend.** `reactions.py` scan + `control.py` load/generate + `runner.py`
   spawn + the API above. Unit tests for `.toml` generation against the strict
   schema, using `musicsim --check` as an oracle (assert the generated control
   passes `--check`).

   **DONE (2026-10-16).** See "Phase 1 results" below.
2. **Frontend: reaction picker + control panel + run flow.** Curated form,
   raw TOML editor, normalize toggle, run/poll, — get to a real `.root`
   produced and served.
+
+   **DONE (2026-10-16).** See "Phase 2 results" below.
3. **Trace (2D) view with JSROOT.** Draw `traces/` for a chosen event; apply
   the normalize toggle.
+
+   **DONE (2026-10-16).** See "Phase 3 results" below.
4. **Path (3D) view with three.js.** Detector stack from the geometry mapping
   + coarse MC path; single-event and overlay modes.
+
+   **DONE (2026-10-16).** See "Phase 4 results" below.
5. **Simulator normalization feature.** Add the normalize control option to the
   sim, implement per-strip beam-gain normalization in `Output.cpp`, update
   the README, and wire the web toggle to it. Validate with `musicsim --check`
   and confirm the emitted units are noted in metadata.
+
+   **DONE (2026-10-16).** See "Phase 5 results" below.
6. **Polish.** Event index sync between 2D and 3D, legend, nice styling,
   error/exception surfacing from `musicsim` logs, README for running.
7. **Flake + review.** Add a `web` devShell (or extend the default devShell)
   with python3/fastapi + node. Call out the build-only flake change and the
   sim normalization schema/output change in the review. Optional: add
   `web/runs/` to `.gitignore`.

## Open questions (surface before build, don't guess)

- **Normalization mechanism (DECIDED).** Automatic reference beam run: the
  reaction control references a `strip = -1` beam run, run first (or reused)
  with a small `events_per_strip`, to build the per-strip gains. Remaining
  sub-details to settle in Phase 5: the exact TOML key, and whether
  `events_MeV` gets a units/flag record.
- **Which controls to bundle by default.** `ControlExamples/` has 9 reactions.
  Scan that, or maintain a curated `web/controls/` subset with nicer labels?
  **Recommend scan `ControlExamples/` first; curated labels are a polish item.**
- **JSROOT can draw the trace TGraphs standalone** (needs the `.root` open
  either way) — confirm in the Phase 0 spike.
- **Multi-event overlay of `events_MeV`.** Does the user also want multiple
  events overlaid in the 2D trace view, or only the 3D overlay? **Recommend
  single-event for 2D, multi-event for 3D in v1.**

## Notes for whichever agent builds this

- The strict TOML schema is the source of truth (README `Strict TOML control
  schema`, `Reaction steps`, `Per-strip TALYS tables`). Generated controls
  must validate with `musicsim --check`.
- `run.output` path AND the surrounding dir matter (tables are written to the
  directory of `run.output`); give each run its own dir.
- `ControlExamples/37Cl_alpha_n/37Cl_alpha_n_bulk.toml` is a good, small,
  cache-free reference control to validate the generated shape against.
- `nix build` output is `result/bin/musicsim`; the devShell sets
  `CATIMA_PREFIX`, `SRIM_TABLE_BIN`, `TALYS_BIN`, etc., so spawning from the
  shell picks up the caches automatically.
- The normalization reference is `~/MUSIC/tooling/src/StripSumScatter.cpp`
  (`SimBeamGains`, `SimTotal`, `SimPopTraces`) and the Python mirrors
  `~/MUSIC/python/data.py` (`max_normalize`) and `~/MUSIC/python/blind_an.py`
  (`_trace`).
