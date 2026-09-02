# Core-physics review checklist — Claude commits, 2026-09-01
<!---->
Per the AI-disclosure section of the README, changes to simulation physics
(kinematics, energy-loss chain, straggling sampling) must be human-reviewed.
Line numbers are as of HEAD (`5c9e238`).
No code is reproduced here; open the
file at the line and review.
Checks indicate review by ewt.
<!---->
Claude-authored commits in the last few days (identified by
`Co-Authored-By: Claude` trailer):
<!---->
| Commit | Title | Model | Date |
| --- | --- | --- | --- |
| `52a0428` | Set the SRIM globals on the master thread only | Opus 5 | 2026-09-01 |
| `5c9e238` | Add stopping = "mean": average catima and SRIM per energy point | Fable 5.1 | 2026-09-01 |
<!---->
Neither commit touches `Kinematics.cpp` or the Vavilov sampler; all physical
change is in the mean energy-loss chain and the plumbing around it.
<!---->
## Core physics — review these lines
<!---->
### Energy-loss chain: new `stopping = "mean"` model (`5c9e238`)
<!---->
- [x] `src/EnergyLoss.cpp:176-177` — the physics itself: per energy-grid
      point, mean dE/dx becomes `0.5 * (catima mean + SRIM table value)` when
      model 2 is selected.
      Check: equal-weight arithmetic mean is what ApJ
      983:142 sec 2.2 describes; the averaging happens on raw model values
      *before* the user `dedx_scale` is applied (scale is multiplied in at
      `src/EnergyLoss.cpp:215` and `:253`), which is the sane ordering; the
      `std::max(0.0, ...)` clamp is preserved.
- [x] `src/EnergyLoss.cpp:141-142` — gate inverted from
      `gStoppingModel != 1` to `== 0`, so `"mean"` now also enters the SRIM
      table path.
      Consequence: a missing table is fatal (`exit(1)` at
      `src/EnergyLoss.cpp:157`) under `"mean"` too, not just `"srim"`.
      Check
      that hard-failing (rather than silently falling back to catima) is the
      intended behavior for the mean model.
- [x] `src/EnergyLoss.cpp:136-140` — comment asserting `"mean"`/`"srim"`
      modify only the *mean* dE/dx, and that straggling stays catima because
      SRIM tables carry no variance.
      Check the claim against the code: the
      variance table `sigma2_per_cm_` is filled from catima at
      `src/EnergyLoss.cpp:133` and is not touched anywhere below.
- [x] `src/EnergyLoss.cpp:159-177` — the loop over the log-energy grid that
      interpolates the SRIM table and blends it in.
      Check the clamp/interpolate
      logic (`:163-173`) is unchanged from the `"srim"` path and that the blend
      cannot run on a partially populated `eloss_per_cm_`.
- [x] `src/ControlFile.cpp:308-309` — new accepted value:
      `physics.stopping = "mean"` maps to `stoppingModel = 2`.
      Check the
      validation branch (`src/ControlFile.cpp:302-316`) rejects everything
      else as before.
<!---->
### Straggling — verify untouched (no lines changed)
<!---->
- [x] `src/EnergyLoss.cpp:128-133` — sanity read: the straggling variance
      still comes from catima's `gStragglingConfig` only. Neither commit
      modifies `sigma2_per_cm_`, `src/VavilovSampler.*`, or the Gaussian
      window/degrader smearing — confirm by eye that nothing downstream of
      line 142 in `BuildTables` writes to `sigma2_per_cm_`.
<!---->
### Threading around the physics model selector (`52a0428`)
<!---->
- [x] `src/ControlFile.cpp:415-424` — new `if (workerId_ == 0)` guard around
      the `gStoppingModel` / `gSrimDir` / `gSrimGasTag` assignments.
      `gStoppingModel` *is* the gas stopping-power model selector, so this
      changes who writes a physics knob. Check: master (worker id 0, default
      at `include/Simulator.hpp:48`) loads the control file before any worker
      is spawned (`src/Simulator.cpp:395-415` spawns ids 1..N, each reloading
      the *same* control file), so workers only ever read. Caveat: if a
      per-worker control file with a different `[physics] stopping` were ever
      introduced, workers would silently run the master's model.
<!---->
## Checked and excluded (not core physics)
<!---->
- `src/EnergyLoss.cpp:186-189` — startup log now prints "mean(catima, SRIM)"
  vs "SRIM".
  (Nit: the missing-table error at `src/EnergyLoss.cpp:154` still
  says `physics.stopping = srim` even when the active model is `"mean"`.)
- `src/Particle.cpp:256-258` — per-event log header labels which stopping
  model was used.
- `src/ControlFile.cpp:294-301, 311-313` — comments and error-message text
  around the `stopping` key.
- `src/srim-cache.cpp:75-76` — the cache tool now also generates tables for
  `"mean"`.
  Tooling, not simulator physics, but a `"mean"` run hard-fails
  without the tables it produces.
<!---->
<!---->
