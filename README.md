# MUSIC simulator
<!---->
Fork of the ANL MUSIC simulator (originally by Daniel Santiago-Gonzalez).
Stopping powers are computed using [catima](https://github.com/hrosiak/catima).
<!---->
## Build
<!---->
Dependencies are managed by Nix.
<!---->
```bash
nix build
./result/bin/musicsim ControlExamples/37Cl_alpha_n/37Cl_alpha_n_bulk.toml
```
<!---->
For iterative development:
<!---->
```bash
nix develop
make
./musicsim ControlExamples/37Cl_alpha_n/37Cl_alpha_n_bulk.toml
```
<!---->
## Control file (TOML)
<!---->
Control files are TOML.
Sections and keys are case-sensitive.
Every section is optional.
Anything you don't write will be its compiled-in default.
Example:
<!---->
```toml
[gas]
species     = "4He"
pressure    = 220.0     # Torr
temperature = 293.0     # K
#
[beam]
species     = "37Cl"
energy      = 92.0      # MeV at the accelerator (before windows)
energy_fwhm = 0.90      # MeV FWHM
dedx_scale  = 1.0
#
[target]
species  = "4He"
compound = "41K"
#
[windows.entrance]
material         = "Ti"
thickness_mg_cm2 = 0.9       # alternatively: thickness_um = X
#
[windows.exit]
material         = "Ti"
thickness_mg_cm2 = 1.3
#
[windows.degrader]
material     = "Mylar"
thickness_um = 6.0           # alternatively: thickness_mg_cm2 = X
#
[detector]
eloss_bins  = 555
max_eloss   = 10.0      # MeV
strip_first = 3         # alternatively: strip = -1 for unreacted-beam runs
strip_last  = 13
eres        = 0.05      # MeV; -1 disables
#
[[reaction.step]]                            # add more [[reaction.step]] for decay chains
[reaction.step.evap]
name       = "4He"
color      = 2
dedx_scale = 1.0
[reaction.step.res]
name       = "37Cl"
color      = 4
dedx_scale = 1.0
#
[run]
n_events  = 10000
threads   = 8
wait      = 0
update    = 0
max_time  = 2000.0      # ns
sim_step  = 0.001       # cm
method    = 0
output    = "traces_37Cl_aa_bulk.root"
file_opt  = "recreate"
print_opt = 0
```
<!---->
The per-event beam-energy chain is:
`beam.energy ± beam.energy_fwhm` (accelerator) → degrader (if any) → entrance
window → `Kbi` (gas surface). At the degrader and each window the beam
traverses a single thick step, which puts the Vavilov parameter κ comfortably
into the Bohr regime (κ ≫ 10) where Vavilov has already collapsed to a
Gaussian — so straggling there is sampled as a single Gaussian draw from
catima's `sigma_E`.
The accelerator energy spread, the degrader straggling,
and the entrance-window straggling are independent Gaussians and add in
quadrature naturally.
<!---->
Particles that exit the downstream face of the gas are propagated through the
exit window before their energy is recorded as `Kbeam_exit` / `Kh_exit[]` /
`Kl_exit[]` on the `MC` tree.
<!---->
### Per-step gas straggling (Vavilov)
<!---->
Inside the active gas the propagator takes ~10 μm steps by default, which puts the
per-step Vavilov parameter κ ≈ 0.05–0.5.
In that band the symmetric Bohr
Gaussian that catima exposes via `sigma_E` is the wrong distribution shape; using it would induce warnings about unphysical energy loss.
The real per-step energy-loss spectrum is asymmetric, with a long high-loss
tail from δ-electron knock-ons and no negative-loss left tail.
Each step
samples instead from the full Vavilov distribution Φ(λ_V; κ, β²) at the
step's κ and β², standardised against catima's mean dE and σ_E so the first
two moments still match catima exactly.
Per-step samples are clamped to
[0, Ki] (energy can only decrease via dE/dx, never increase).
<!---->
The sampler uses the convolution decomposition of Yi & Han[^vavmc] —
Φ(λ_V; κ, β²) factors as Φ(λ_V; (1−β²)κ, 0) ⋆ Φ(λ_V; β²κ, 1), so a draw
at any (κ, β²) is the sum of two independent draws from 1-D κ-indexed
tables.
The tables (`src/VavilovSampler.cpp`) are pre-filled once at startup
from ROOT's `Math::VavilovAccurate` (~3 ms of one-time work) and looked up
via bilinear interpolation in the hot path — no per-step VavilovAccurate
construction. Outside the Vavilov band (κ < 10⁻³ or κ > 10) the sampler
falls back to a Gaussian, matching the Landau and Bohr limits respectively.
<!---->
[^vavmc]: Chul-Young Yi, Hyon-Soo Han, *A Monte Carlo algorithm for the
Vavilov distribution*, Nucl. Instr. and Meth. in Phys. Res. B **149**(3),
263–271 (1999). [doi:10.1016/S0168-583X(98)00803-9](https://doi.org/10.1016/S0168-583X(98)00803-9).
<!---->
### Reaction strip selection
<!---->
Set **one** of: `detector.strip = N`, or both `detector.strip_first` and
`detector.strip_last`. `strip = -1` runs unreacted-beam (no reaction vertex).
<!---->
### Stopping-power models
<!---->
catima exposes a couple of physics-model switches that can shift dE/dx
predictions noticeably for heavy ions at low velocity. They're set through
an optional `[physics]` section:
<!---->
```toml
[physics]
z_effective = "atima14"   # default
low_energy  = "srim_95"   # default
```
<!---->
- `z_effective` controls the projectile-charge-state model.
Valid:
  `none`, `pierce_blann`, `anthony_landorf`, `hubert`, `winger`,
  `schiwietz`, `global`, `atima14`.
  catima's compiled default is
  `pierce_blann`; `atima14` is used here.
- `low_energy` controls the table used below the Bethe-Bloch regime
  (where 37Cl in Ti sits at ~0.3 MeV/u, for instance). Valid: `srim_85`
  (catima's compiled default) or `srim_95` (newer SRIM/Ziegler curves;
  default here).
<!---->
### Disabling layers
<!---->
Any of the four physical layers can be turned off independently — the
simulator prints a warning at startup and skips both the mean dE/dx and the
straggling for that layer.
<!---->
| Disable                         | How                                                                            |
| ---                             | ---                                                                            |
| Gas energy loss + straggling    | `gas.pressure = 0` (or any non-positive value)                                |
| Entrance window                 | Set `thickness_mg_cm2` (or `thickness_um`) to `-1` under `[windows.entrance]` |
| Exit window                     | Same, under `[windows.exit]`                                                  |
| Degrader                        | Omit the `[windows.degrader]` block, or set its thickness to `-1`             |
<!---->
Each of the three layer sub-tables (`[windows.entrance]`, `[windows.exit]`,
`[windows.degrader]`) accepts **exactly one** of `thickness_mg_cm2` (areal
density, mg/cm²) or `thickness_um` (linear length along the beam axis, μm).
Setting both in the same sub-table is an error.
Pick whichever matches the
spec sheet for the physical layer; the loader dispatches to catima's areal
or linear thickness path accordingly.
<!---->
### Multi-threading
<!---->
`run.threads = N` fans the event loop out across `N` worker threads using
`std::async`.
Each worker is a separate `MUSIC_Simulator` instance with its
own RNG seed and per-worker ROOT output; the master pre-warms catima's
internal cache and then merges the per-worker output files with `TFileMerger`
into the configured `run.output` at the end of the run.
Workers force
`update=0`/`wait=0`, so visualization is only available in single-threaded
mode (the default, `threads = 1`).
<!---->
`reaction.step` is an array of tables — add more `[[reaction.step]]` blocks
for multi-step decay chains; up to `MaxNumEvapPart` (10) are supported.
`evap.dedx_scale` and `res.dedx_scale` are optional and default to `1.0`.
Particle names use `AEl` notation (e.g. `37Cl`, `4He`, `1H`, `n`).
<!---->
### Bringing over upstream `.msc` files
<!---->
The pre-TOML `.msc` format used by the original ANL musicsim isn't accepted by
this fork.
A python converter is included:
<!---->
```bash
tools/legacy_msc_to_toml.py path/to/upstream/file.msc
# writes path/to/upstream/file.toml; per-key comments are preserved.
```
<!---->
Glob-friendly: `tools/legacy_msc_to_toml.py /path/to/musicsim/Examples/**/*.msc`.
<!---->
`Kb` is mapped to `beam.energy` with an inline warning — the meaning shifted
from "gas-surface KE" (upstream) to "accelerator KE before windows" (Remix),
so adjust the value by hand if you want the windows modelled.
SRIM-table keys
and `AnodeGeom` are dropped (catima handles dE/dx; geometry is hardcoded).
<!---->
## Output tree layout
<!---->
Each run produces one ROOT file with two trees, mirroring the format produced by
the upstream `EventBuilderNearestGrid` analysis used [here](https://github.com/ewtodd/MUSIC/tree/main/ProductionMode_37Cl/macros):
<!---->
- **`events_MeV`** — detector-level branches.
`LeftdE[18]`, `RightdE[18]`, and
  `TotaldE[18]` (Float, MeV).
  Energies are MeV truth (hence the `_MeV` suffix
  vs the upstream ADC-valued `event` tree), so analysis macros should
  calibrate data to MeV rather than rescaling sim to ADC.
- **`MC`** — truth-only branches, friended to `events_MeV`: reaction strip, beam
  kinematics (`BeamEnergyAccel` at the accelerator, `Kbi` at the gas surface,
  `Kbr` at the reaction point, `Kbeam_exit` after the exit window),
  light/heavy four-vector components, vertex and exit positions, per-evap
  exit energies (`Kl_exit[]`, `Kh_exit[]`), and dead-layer energy deposits
  `DeadUS_dE` / `DeadDS_dE` (MeV) for the unread gas regions upstream and
  downstream of the readout strips.
  Sentinel values on exit-energy branches:
  `-1` = particle stopped in the gas, `-2` = no such particle for this event
  (e.g. `Kbeam_exit` on a reacted event).
  A literal `0` means the particle
  reached the exit window but stopped inside it.
<!---->
`Hits[k]` triggers when the corresponding channel energy exceeds
`max(0.02 MeV, 3·Eres)` so that pure-noise cells are not counted.
<!---->
## A note on AI-assisted development
<!---->
Parts of this fork — the catima migration, std::async multi-threading,
window / degrader plumbing, and most boilerplate around them — were
written with the help of [Claude Code](https://claude.ai/claude-code).
All changes to the simulation physics (kinematics, energy-loss chain,
straggling sampling) are human-reviewed and approved before being
committed.
