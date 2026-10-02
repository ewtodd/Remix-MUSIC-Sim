# MUSIC simulator

```
 ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩
   ██████╗  ███████╗ ███╗   ███╗ ██╗ ██╗  ██╗
   ██╔══██╗ ██╔════╝ ████╗ ████║ ██║ ╚██╗██╔╝
   ██████╔╝ █████╗   ██╔████╔██║ ██║  ╚███╔╝
   ██╔══██╗ ██╔══╝   ██║╚██╔╝██║ ██║  ██╔██╗
   ██║  ██║ ███████╗ ██║ ╚═╝ ██║ ██║ ██╔╝ ██╗
   ╚═╝  ╚═╝ ╚══════╝ ╚═╝     ╚═╝ ╚═╝ ╚═╝  ╚═╝
 ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩ ♪ ♫ ♬ ♩
```

Remix-MUSIC-Sim is a fork of the ANL MUSIC simulator originally written by
Daniel Santiago-Gonzalez. Stopping powers are computed with
[catima](https://github.com/hrosiak/catima).

## AI-assisted development disclosure

Parts of this fork—including the catima migration, multi-threading, window and
degrader plumbing, and supporting code—were written with help from
[Claude Code](https://claude.ai/claude-code),
[son-of-anton](https://github.com/ewtodd/son-of-anton), and Codex. Simulation
physics changes must be reviewed by a human before they are committed. The
unchecked review files under `checks/` identify the exact source lines that
need sign-off.

## Build and verification

Dependencies are managed by Nix:

```bash
nix build
./result/bin/musicsim ControlExamples/37Cl_alpha_n/37Cl_alpha_n_bulk.toml
```

`nix build` also ships the regression-suite binaries
(`legacy-msc-to-toml`, `root_checks`, `vavilov_statistics`, `nuclide_loader`)
in `result/bin`. To run the suite — schema, converter, nuclide-loader,
distribution, conservation, kinematics, boundary, reproducibility, and ROOT
update tests — build first, then:

```bash
bash tests/run-tests.sh
```

It finds everything under `result/bin` by default; the `MUSICSIM_BIN`,
`CONVERTER_BIN`, and `*_CHECK_BIN` environment variables override each path.

For working on the sources inside the dev shell (`nix develop`), the Makefile
targets are still there: `make -j2` builds by hand, `make -j2 check` runs the
suite against the local build, `make strict` treats the repository's extended
GCC warning set as errors, `make -j2 sanitize-check` repeats the suite under
AddressSanitizer and UndefinedBehaviorSanitizer, and `make format-check`
checks C/C++ with clang-format and TOML with Taplo.

To validate a control without running events:

```bash
./musicsim --check path/to/control.toml
```

## Strict TOML control schema

Names and keys are case-sensitive. Unknown keys, wrong types, non-finite
numbers, invalid ranges, and incomplete reaction definitions are errors. A
failed load is transactional: no partial configuration is used.

The required top-level tables are `[gas]`, `[beam]`, `[target]`, `[detector]`,
and `[run]`. `[windows]` and `[physics]` are optional. A reacted run also needs
one to ten complete `[[reaction.step]]` entries. The loader validates every
nuclide, beam-plus-target compound A/Z, per-step A/Z, and decay-chain
continuity.

This is a complete representative control:

```toml
[gas]
species     = "4He"
pressure    = 220.0 # Torr; 0 disables gas energy loss
temperature = 293.0 # K

[beam]
species     = "37Cl"
energy      = 92.0 # MeV at the accelerator, before any degrader/window
energy_fwhm = 0.90 # MeV FWHM; optional, default 0
dedx_scale  = 1.0  # optional gas stopping multiplier

[target]
species  = "4He"
compound = "41K"

[windows.entrance]
material         = "Ti"
thickness_mg_cm2 = 0.9 # use exactly one thickness representation
dedx_scale       = 1.0

[windows.exit]
material         = "Ti"
thickness_mg_cm2 = 1.3

[windows.degrader]
material     = "Mylar"
thickness_um = 6.0

[detector]
eloss_bins  = 555
max_eloss   = 10.0 # MeV
strip_first = 3    # or use one `strip`; do not combine the forms
strip_last  = 13

# Choose at most one detector-noise model. A scalar applies to all 34 anode
# electrodes but not the cathode. A channel table may use Cathode, S0, S17,
# and L1..L16/R1..R16. Missing channels and -1 are disabled.
resolution_fwhm_percent = 5.0
# noise_sigma_mev = 0.05

[[reaction.step]]
residue_excitation  = "ground"
angular_distribution = "rutherford"
theta_cm_min_deg     = 1.0

[reaction.step.evap]
name       = "4He"
color      = 2
dedx_scale = 1.0

[reaction.step.res]
name       = "37Cl"
color      = 4
dedx_scale = 1.0

[physics]
z_effective            = "atima14"
low_energy             = "srim_95"
straggling_z_effective = "pierce_blann"
straggling             = true
stopping               = "catima"

[run]
events_per_strip = 10000
threads          = 8
wait             = 0
update           = 0
max_time         = 2000.0 # ns
sim_step         = 0.001  # cm, finite and > 0
method           = 0
output           = "traces_37Cl_aa_bulk.root"
file_opt         = "recreate"
print_opt        = 0
seed             = 12345
```

### Required and optional values

| Table | Required | Optional defaults and constraints |
| --- | --- | --- |
| `gas` | `species`, `pressure`, `temperature` | Pressure may be zero to disable gas loss; temperature must be positive. |
| `beam` | `species`, `energy` | `energy_fwhm = 0`, `dedx_scale = 1`; energy and scale must be positive. |
| `target` | `species`, `compound` | Both are required even for an unreacted-beam run. |
| `detector` | `eloss_bins`, `max_eloss`, and one strip selector | `strip = -1` is unreacted; otherwise use `strip = 0..17` or a complete `strip_first`/`strip_last` range. |
| `run` | `events_per_strip`, `output`, `seed` | `threads = 1`, `wait = 0`, `update = 0`, `max_time = 2000`, `sim_step = 0.001`, `method = 0`, `file_opt = "recreate"`, `print_opt = 0`. |

`events_per_strip` is literal: selecting strips 3 through 13 with 10,000 events
writes 110,000 entries. `seed` must be in `[1, 4294967295]`. The same
`(seed, requested strip, event_index)` tuple produces the same event regardless
of worker count or scheduling.

`run.update = 1` enables interactive ROOT visualization; `run.wait = 1` pauses
between visualized events. These options are unavailable with multiple
threads. `run.file_opt = "update"` instead means append to an existing ROOT
file; it is supported only with one thread.

### Detector noise

The two explicit and mutually exclusive noise keys are:

- `noise_sigma_mev`: absolute Gaussian sigma in MeV.
- `resolution_fwhm_percent`: Gaussian FWHM as a percentage of that event's
  channel deposit.

Each key accepts either a scalar or a table. For example:

```toml
[detector.noise_sigma_mev]
Cathode = 0.05
S0      = 0.04
L1      = 0.06
R1      = 0.05
S17     = -1
```

### Reaction steps

Every reacted step requires `evap.name`, `res.name`, `residue_excitation`, and
`angular_distribution`. Product `color` and `dedx_scale` are optional. Particle
names use `AEl` notation, such as `37Cl`, `4He`, `1H`, and `n`.

The excitation model is selected independently for each step:

- `ground`: the residue is left in its ground state.
- `uniform`: excitation is uniform from zero through the energy available at
  that step.
- `forced`: an intermediate residue receives at least the excitation needed to
  permit the remaining configured decay chain. It is rejected on the final
  step.

The angular model is also per step:

- `isotropic`: uniform in solid angle.
- `rutherford`: inverse-CDF Rutherford sampling above
  `theta_cm_min_deg` (default 1 degree).

An elastic first step—products equal to target plus beam—is required to use
`ground` excitation and `rutherford` angles.

### Layers and stopping-power models

Each present window/degrader table must contain exactly one of
`thickness_mg_cm2` or `thickness_um`. A zero or negative thickness disables the
layer; a positive layer also requires `material`. Supported materials are Ti,
Al, Havar, Kapton, and Mylar. `dedx_scale` defaults to 1 and is implemented as
equivalent thickness, so it changes both mean loss and straggling.

The beam-energy path is accelerator energy plus independent Gaussian beam
spread, then degrader, entrance window, and gas. Window/degrader output is
clamped to `[0, incoming energy]`. Only downstream gas exits traverse the exit
window, with planar thickness scaled by `1 / |cos(theta)|`. Exit-window
straggling is intentionally omitted because no downstream detector is modeled.

`physics.stopping` selects `catima`, `srim`, or the arithmetic `mean` of the two
mean stopping powers. SRIM-backed modes require tables produced in advance:

```bash
./srim-cache path/to/control.toml
```

The remaining physics defaults are shown in the full example. Valid effective
charge names are `none`, `pierce_blann`, `anthony_landorf`, `hubert`, `winger`,
`schiwietz`, `global`, and `atima14`. `low_energy` accepts `srim_85` or
`srim_95`. The straggling charge model cannot be `atima14`, because that catima
path does not produce a variance. `physics.straggling = false` disables both
gas and entrance-layer/degrader energy-loss fluctuations; detector noise is
independent.

### Gas straggling and transport

Gas transport ends a step exactly at the first cell edge, detector surface,
requested vertex, time limit, or stopping point it reaches. Every traversed
fraction contributes its energy deposit. Passive matter cannot add energy.

For Vavilov-band gas steps, a `(kappa, beta^2, probability)` quantile table is
built from ROOT's `Math::VavilovAccurate` and interpolated during transport.
The table spans `kappa = 1e-3..10`; thinner steps use the lowest supported,
Landau-like Vavilov shape, while thicker steps use the Gaussian limit. The
distribution is standardized with ROOT's analytic moments and scaled to
catima's step variance. Interpolation and physical loss truncation mean sampled
moments are close to, but not mathematically identical to, the catima moments.
Kappa uses the same effective-charge prescription as the configured
straggling variance.

### Multi-threading and output updates

Multi-threaded runs create isolated simulator instances and temporary ROOT
files. Catima calls are serialized because its cache is process-global. The
driver checks every worker result and entry count, merges into a temporary
file, validates it, and atomically publishes the final output. Per-worker logs
use distinct names.

Multi-threaded `file_opt = "update"` is rejected. A single-threaded update is
also transactional: it copies the destination to `<output>.tmp`, requires an
exact tree schema and byte-identical embedded control TOML, checks unique and
contiguous event identities, appends, validates, and then atomically replaces
the destination. A stale `.tmp` file is never overwritten automatically.

## Legacy `.msc` conversion

The converter is a dependency-free C++17 executable built by `make`. The old
Python path remains as a compatibility launcher:

```bash
make legacy-msc-to-toml
./legacy-msc-to-toml path/to/input.msc
tools/legacy_msc_to_toml.py path/to/input.msc
```

By default it writes a sibling `.toml`, preserves the source, refuses to
overwrite a destination, and fails for missing/unmatched inputs. Use `--force`
to overwrite, `--remove-source` to delete the source only after successful
publication, `--stdout` to inspect output, or `-o FILE` for one explicit
destination. `--keep` states the default preservation policy and cannot be
combined with `--remove-source`.

The converter escapes TOML strings and emits the strict schema. It marks
inferred excitation/angular models with `REVIEW` comments; those choices must
be checked before running. Legacy `Kb` is treated as accelerator energy, while
upstream files used gas-surface energy, so window-enabled conversions also need
manual review. Unsupported trace-database and obsolete SRIM/geometry settings
are rejected or reported.

## ROOT output

Each successful run publishes one ROOT file containing a single cycle of two
equal-length trees. `events_MeV` has `MC` as a friend.

`events_MeV` contains Float_t energies in MeV:

| Branch | Meaning |
| --- | --- |
| `LeftdE[16]`, `RightdE[16]` | Segmented strip 1..16 deposits, indexed by strip minus one. |
| `Strip0dE`, `Strip17dE` | Full-width end-strip deposits. |
| `Cathode` | Sum seen by the cathode, with its independently configured noise. |

`MC` contains:

| Branches | Meaning |
| --- | --- |
| `n_steps`, `event_index` | Configured reaction-step count and deterministic per-strip event identity. |
| `requested_strip`, `reaction_strip` | Configured vertex strip and actual successful reaction strip; the latter is `-1` when no reaction occurred. |
| `beam_energy_accel`, `beam_energy_gas`, `beam_energy_reaction`, `beam_energy_exit` | Beam kinetic-energy chain in MeV. |
| `vertex_x`, `vertex_y`, `vertex_z` | Successful reaction vertex in cm. |
| `beam_stop_x/y/z`, `beam_stop_strip`, `beam_termination` | Final beam state and termination. |
| `evap_energy[n_steps]`, `residue_energy[n_steps]` | Per-step product kinetic energies at creation. |
| `theta_cm[n_steps]`, `phi_cm[n_steps]` | Per-step sampled center-of-mass angles in degrees. |
| `evap_theta/phi[n_steps]`, `residue_theta/phi[n_steps]` | Per-step lab angles in degrees. |
| `evap_energy_exit[n_steps]`, `evap_stop_x/y/z[n_steps]`, `evap_stop_strip[n_steps]`, `evap_termination[n_steps]` | Light-product final states. |
| `residue_energy_exit[n_steps]` | Exit-energy slots; only the surviving residue can be populated. |
| `residue_step`, `residue_stop_x/y/z`, `residue_stop_strip`, `residue_termination` | Surviving-residue identity and final state. |
| `DeadUS_dE`, `DeadDS_dE` | Deposits in unread upstream/downstream gas regions, in MeV. |

Exit-energy sentinels are `-2` for not applicable/not propagated and `-1` for a
transported particle that did not reach the downstream exit. A literal zero
can mean that the particle reached but stopped in the exit foil. Stop-strip
sentinels are `-2` for not applicable and `-1` for an exit or a stop outside a
readout strip.

Termination codes are stored in `metadata/termination_codes` and currently
mean: 0 not propagated, 1 reached the requested vertex, 2 stopped, 3 downstream
exit, 4 upstream exit, 5 side exit, 6 timeout, and 7 invalid state.

The `metadata` directory also stores the software version, control path, exact
control text, resolved configuration (schema version 3), and straggling model.
These records are part of update-file compatibility checks.

## Physics scope and limitations

- A requested reaction vertex is sampled uniformly within each selected strip.
  The simulator conditions on that choice; it does not calculate an absolute
  reaction yield, cross section, beam attenuation, or target-density weighting.
- Transport is straight-line. Multiple Coulomb scattering, charge exchange,
  angular straggling, and additional nuclear interactions are not modeled.
- The detector is the fixed 20-row, two-column MUSIC geometry encoded in the
  source.
- Gas stopping/variance tables cover total kinetic energy from `1e-4` to
  `1e4` MeV and clamp outside that range.
- `forced` and `uniform` excitation are continuum models, not discrete nuclear
  level schemes with spin/parity, widths, branching ratios, or lifetimes.
- Detector response includes only the selected independent Gaussian noise
  model; thresholds, saturation, gain drift, correlations, digitization, and
  non-Gaussian electronics response are absent.

## C++ API compatibility

The command-line TOML workflow is the supported interface. Unsafe or no-op
legacy public APIs were removed, including angular-grid simulation,
trace-database generation, excitation-range/legacy energetics helpers,
`WriteTraces`, and `Particle::GetEnergyLoss`. Downstream C++ code that called
those routines must migrate; the project does not currently promise source ABI
compatibility for its internal simulation classes.
