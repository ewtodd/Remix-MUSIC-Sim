# Review checklist — window/degrader dedx_scale + SRIM-nix table generator, 2026-09-03
<!---->
Per the AI-disclosure section of the README, changes to simulation physics
(kinematics, energy-loss chain, straggling sampling) must be human-reviewed.
Line numbers are as of HEAD (`5bf30f2`).
No code is reproduced here; open the
file at the line and review.
Checks indicate review by ewt.
<!---->
Claude-authored commits in the last few days (identified by
`Co-Authored-By: Claude` trailer):
<!---->
| Commit | Title | Model | Date |
| --- | --- | --- | --- |
| `5bf30f2` | Read dedx_scale on the window and degrader layers | Fable 5.1 | 2026-09-03 |
| `a260dff` | Take make-srim-table from the SRIM-nix flake | Fable 5.1 | 2026-09-03 |
<!---->
`5bf30f2` touches the beam energy-loss chain: the window and degrader layers
finally read the `dedx_scale = 1.08367` the 87Rb controls have carried since
the window loss was calibrated to 46.9 MeV.
`a260dff` is tooling: which generator the SRIM-table cache (`srim-cache`)
drives. It changes how `stopping = "srim"/"mean"` tables are produced, not
the simulation physics.
<!---->
## Core physics — review these lines
<!---->
### Layer stopping-power scale, applied as equivalent thickness (`5bf30f2`)
<!---->
- [x] `src/ControlFile.cpp:133-139` — the key is read per layer (defaults to
      1.0 when absent, via the member initializers at
      `include/Simulator.hpp:394-399`) and validated `> 0` →
      `exit(EXIT_FAILURE)`, matching the file's error convention (cf.
      `src/ControlFile.cpp:152-156`).
      Check that rejecting zero/negative (rather than treating them as
      "disable the layer") is intended.
NOTE by ewt: Yes, disabling the layer is done by giving a thickness < 0. So this is intended.
- [x] `src/Materials.cpp:138-143` — the comment asserting the equivalence:
      scaling dE/dx by s and scaling the length by s are the same integral,
      and the one difference is straggling, which grows as the square root of
      the thickness.
      Check the claim against `EnergyThroughWithStraggling`
      (`src/Materials.cpp:199-218`): `sigma_E` is computed from
      `catima::calculate` on the *scaled* material, so σ ∝ √(s·L) — a layer
      that stops harder is treated as a thicker one, straggling included.
      This is the key physics decision to sign off on: the layer scale
      changes **both** the mean and the straggling, unlike the beam /
      reaction-step `dedx_scale`, which scales only the mean
      (`src/EnergyLoss.cpp:248-249`).
NOTE by ewt: agreed in principle; will revise the approach if it proves problematic. so far so good for 37Cl/87Rb. straggling generally seems to contribute little in this sim anyway. 
- [x] `src/Materials.cpp:145-149` — the physics itself: `buildLayer`
      multiplies the thickness by `scale` in both branches (mg/cm² →
      `BuildSolidMaterial`, μm → `BuildBulkMaterial`).
      Check: scaling the areal density and scaling the linear length are the
      same scaling of the material (the conversion between the two is catima's
      constant density per material), so either sub-table form gives a
      consistent layer.
- [x] `src/Materials.cpp:150-153, 158-161` — entrance and exit scales
      applied when the materials are built.
- [x] `src/Materials.cpp:173-179` — degrader scale applied in both branches
      of `BuildDegrader`.
- [x] `src/EventLoop.cpp:152-155` — the per-event chain: degrader and
      entrance window both go through `EnergyThroughWithStraggling`, so
      their scales widen the per-event beam-energy Gaussian by √scale
      (×1.041 at 1.08367) on top of the mean shift.
- [x] `src/Propagation.cpp:317-321` — the exit window is asymmetric
      (pre-existing design, now carrying a scale): exit traversal is
      `EnergyOutOfMaterial`, mean only, no straggling sampling — so the exit
      scale shifts the `Kbeam_exit` / `Kh_exit` / `Kl_exit` means but adds no
      straggling there.
      Check that the exit window is deliberately straggling-free.
NOTE by ewt: there's no point in calculating straggling on the exit window because we never measure after it. only case would be if there was Si diode data after, but I do not have that. it's an easy change if we ever need to.
- [x] `src/Simulator.cpp:235-242` — the startup `Kb_at_gas` estimate goes
      through the scaled materials (mean-only path), so the printed beam
      chain reflects the scale; it feeds `SetInitialKinematics` at
      `src/EventLoop.cpp:158`.
      Check the estimate path is the only place the scaled material is used
      at startup.
<!---->
### The calibration itself
<!---->
- [x] 46.9 / 1.08367 = 43.28 MeV — recompute the factor against the 87Rb
      window-loss calibration (nominal loss 43.2 MeV per the commit message)
      and confirm the constant, including that the accompanying √1.08367 ≈ 1.041
      broadening of the window/degrader
      straggling is consistent with what was measured.
      Per the commit message, every run that set the key needs re-running.
NOTE by ewt: re-run in progress
<!---->
## Tooling — SRIM table generator from the flake (`a260dff`)
<!---->
- [x] `flake.nix:10-13` — `srim-nix` as a flake input with its own nixpkgs
      (`flake.lock:72-88`: nixpkgs-unstable pinned 2025-11-02, rev
      `7241bcb`), left alone because the wine wrapper is pinned to what that
      flake was tested with.
      Check that keeping the input's nixpkgs un-pinned-by-us is intended, and
      that `srim-nix` is pinned by the lock (rev `2a4b65e`,
      `flake.lock:96-113`) rather than by a rev in the URL.
- [x] `flake.nix:28, 49-52` — `make-srim-table` is selected from srim-nix and
      its store path is passed to make, so only `srim-cache` (not `musicsim`)
      carries it.
      Check that `musicsim` itself needs no generator at runtime (it only
      reads tables).
- [x] `flake.nix:70-74, 78` — the dev shell puts `make-srim-table` on PATH
      *and* exports `SRIM_TABLE_BIN` pointing at the same store path; both
      agree, env wins.
- [x] `Makefile:33-38` — bakes the store path into `srim-cache` via
      `-DMUSICSIM_SRIM_TABLE_BIN`; empty means PATH lookup at run time.
- [x] `src/srim-cache.cpp:38-45` — precedence: `SRIM_TABLE_BIN` env >
      build-baked path > bare `make-srim-table` on PATH (pre-commit behavior
      preserved).
      Check that env-first is the intended override mechanism (e.g. testing a
      new SRIM build).
- [x] `src/srim-cache.cpp:145-149` — the generator path is single-quoted in
      the `system()` command; a Nix store path is quote-safe.
      Note: there is no existence check on the baked-in path — on a machine
      where that store path is absent, generation fails with the new message
      (`src/srim-cache.cpp:151-156`) rather than falling back to PATH.
      Check that's acceptable for `nix build` output run elsewhere.
<!---->
## Checked and excluded (not core physics)
<!---->
- `README.md:241-246` — new per-layer `dedx_scale` documentation; consistent
  with the behavior, including the separation from the `[beam]` /
  reaction-step scales (gas only).
- `src/ControlFile.cpp:8-12` — header comment now lists `dedx_scale` under
  the window sub-tables.
- `src/Materials.cpp:165-168` — verbose "Window stopping scale … (applied as
  thickness)" line.
- `include/Simulator.hpp:391-393` — comment on the layer fields.
- `src/srim-cache.cpp:1-11, 151-156` — header comment and failure-message
  text.
<!---->
Nits:
- `src/Materials.cpp:185-186` — the `EnergyOutOfMaterial` comment still says
  "Used for thin entrance/exit windows"; the entrance window goes through the
  straggling path (and now carries a scale), so only the *exit* window uses
  this function.
- `src/Simulator.cpp:244-253` — the verbose beam chain prints the *nominal*
  window/degrader thickness; the scale is visible only in the separate
  "Window stopping scale" line, which doesn't cover the degrader at all.
<!---->
