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

`nix build` also ships `srim-cache` and `talys-cache` and the regression-suite
binaries (`legacy-msc-to-toml`, `root_checks`, `vavilov_statistics`,
`nuclide_loader`, `excitation_table`, `talys_tables`) in `result/bin`. To run
the suite — schema, converter, nuclide-loader, distribution, level-table,
conservation, kinematics, branch-weighting, boundary, reproducibility, and
ROOT update tests — build first, then:

```bash
bash tests/run-tests.sh
```

It finds everything under `result/bin` by default; the `MUSICSIM_BIN`,
`CONVERTER_BIN`, `TALYSCACHE_BIN`, `SRIMCACHE_BIN`, and `*_CHECK_BIN`
environment variables override each path. The suite runs `talys-cache` and
`srim-cache` against stand-ins for TALYS and SRIM, so it needs neither.

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

### Binary cache

`nix build` includes TALYS for `talys-cache`, with its 8 GB nuclear-structure
database, and builds catima and TALYS from source. TALYS's source comes from
nds.iaea.org, which sometimes refuses automated downloads. e-desktop serves
its Nix store as a binary cache. Another machine can fetch these from it
instead, along with this simulator at any revision built there:

- URL: `https://cache.ethanwtodd.com`
- Signing key: `e-desktop:35K0AY3HcDOSHVQ/lklmbvmrXjIspM/LYf7yek5lyVA=`

Add both to that machine's Nix configuration, for example in `nix.conf`:

```ini
extra-substituters = https://cache.ethanwtodd.com
extra-trusted-public-keys = e-desktop:35K0AY3HcDOSHVQ/lklmbvmrXjIspM/LYf7yek5lyVA=
```

On NixOS, use `nix.settings.substituters` and `nix.settings.trusted-public-keys`.
Do not make the flake's inputs follow another nixpkgs (for example
`talys-nix.inputs.nixpkgs.follows`). That changes the derivation hashes, and
the cache then has nothing to serve.

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

# Optional per-strip normalization to 1 a.u. (beam-reference gains). When
# normalize_per_strip = true, every event's per-strip dE in the events_MeV tree
# is multiplied by the strip's gain from strip_gains (18 values, strip 0..17).
# A gain of 0 drops that strip's contribution; the Cathode full-plate sum is
# left in MeV (it has no per-strip meaning). Omit both keys to disable.
# normalize_per_strip = true
# strip_gains = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
#                1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

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
z_effective            = "winger"
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
| `detector` | `eloss_bins`, `max_eloss`, and one strip selector | `strip = -1` is unreacted; otherwise use `strip = 0..17` or a complete `strip_first`/`strip_last` range. `normalize_per_strip = false`, `strip_gains = [1.0 × 18]` (per-strip beam-gain multipliers). |
| `run` | `events_per_strip`, `output`, `seed` | `threads = 1`, `wait = 0`, `update = 0`, `max_time = 2000`, `sim_step = 0.001`, `method = 0`, `file_opt = "recreate"`, `print_opt = 0`. |

`events_per_strip` is literal: selecting strips 3 through 13 with 10,000 events
writes 110,000 entries. `seed` must be in `[1, 4294967295]`. The same
`(seed, requested strip, event_index)` tuple produces the same event regardless
of worker count or scheduling, and different seeds give independent events.

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
- `talys`: excitation is drawn from the level table named by `level_table`,
  a path resolved relative to the control file's directory. For TALYS's own
  channels and tables per strip, use a [TALYS class](#talys-classes) instead.

The angular model is also per step:

- `isotropic`: uniform in solid angle.
- `rutherford`: inverse-CDF Rutherford sampling above
  `theta_cm_min_deg` (default 1 degree).

A step may instead mix models through weighted `[[reaction.step.branch]]`
tables. Each branch requires `residue_excitation`, `angular_distribution`,
and `weight`, and may set `theta_cm_min_deg` (default: the step's value) and
`level_table`. A branched step must not also set `residue_excitation` or
`angular_distribution`; `weight` is only valid inside a branch. Weights are
relative, must be finite and `>= 0`, need at least one `> 0`, and are
normalized at load. A step has at most 8 branches, and `forced` is not
allowed in a branch. One branch is drawn per reacted event before its
excitation and angles; a step without branch tables draws nothing extra, so
its event stream is unchanged.

```toml
[[reaction.step]]
[reaction.step.evap]
name = "4He"

[reaction.step.res]
name = "37Cl"

[[reaction.step.branch]] # elastic
residue_excitation   = "ground"
angular_distribution = "rutherford"
weight               = 0.99

[[reaction.step.branch]] # inelastic
residue_excitation   = "talys"
level_table          = "talys_levels_37Cl_aap.tab"
angular_distribution = "isotropic"
weight               = 0.01
```

Branch weights are user-supplied cross-section ratios within the simulated
acceptance (for example, Rutherford integrated above `theta_cm_min_deg`
against the Talys inelastic cross section); a TALYS class, below, takes them
from TALYS instead. A `talys` branch's weight is its cross section with the
whole table open.

A level table is whitespace-separated text with `#` comments. `ground 0 w`
and `level E_ex w` rows are point masses; `bin lo hi w` is a uniform
continuum over `[lo, hi]` MeV with total weight `w`. Weights are relative.
Unknown tags, `hi <= lo`, negative values, non-numeric fields, and a table
with no positive weight are errors naming the file and line.

At each vertex, only the part of a table at or below the energy available is
open; a `bin` straddling it counts pro rata. A branch is drawn with
probability proportional to its weight times its open fraction (always 1 for
`ground` and `uniform`), still with one random draw, and its `E_ex` is drawn
from the open part, so it never exceeds the energy available. As levels
close, the inelastic share shrinks and the open branches take it up. The
reaction is rejected, like a below-threshold step, only when every branch is
closed, for example a step-level `talys` table with nothing open.

An elastic first step—products equal to target plus beam—must list the
target-like product as `evap` and the beam-like product as `res`, because
`theta_cm` is the `res` direction; the swapped labelling is rejected. It is
required to use `ground` excitation and `rutherford` angles. With branch
tables, it must instead include at least one `ground` + `rutherford` branch
with weight `> 0`; its other branches may be inelastic. `rutherford` is valid
only on step 0, in a branch or not: later steps are decays of a moving
residue, not scattering off the beam.

### TALYS classes

MUSIC sorts reactions by the residue's charge, so it cannot tell (a,n) from
(a,2n): after 87Rb + a, both leave yttrium. A TALYS class simulates that
mixture. Instead of `[[reaction.step]]` tables, `[reaction]` names the final
residue's element:

```toml
[reaction]
talys_class        = "Y"  # every TALYS channel ending in yttrium
theta_cm_min_deg   = 30.0 # the elastic cut; only in the beam's class
residue_dedx_scale = 1.0  # every residue; default 1
evap_dedx_scale    = 1.0  # every light product; default 1
```

The class is every exclusive channel TALYS gives that ends in that element.
For 87Rb at E_cm 14 MeV, "Y" is (a,g), (a,n), (a,2n), and (a,gn), which emits
a gamma first; "Sr" is (a,p), (a,np), (a,pn), and (a,d); "Rb" is mainly (a,a')
and elastic. A channel is a chain of steps in TALYS's emission order, so (a,np)
and (a,pn) are separate channels. A gamma step emits a photon, transported
like a neutron. The target must be `4He`, the run must be reacted, and
`run.update` must be 0.

Each strip has one TALYS table, computed at the CM energy of the mean beam at
the strip's centre: no straggling, integrated in `sim_step` slices as the
transport does. It gives each channel's exclusive cross section and, for each
step, the excitation of the residue that step leaves, in mb, from TALYS's full
decay cascade. Each event draws one channel, with probability proportional to
its cross section times the open fraction of its first step's table at the
vertex, as for weighted branches, and runs that channel's chain. Each step's
excitation is drawn from its table between a floor and the energy available.
The floor is the lowest excitation from which the rest of the chain can
happen: the next emission's separation energy plus the lowest state its own
table can reach. The open fraction is taken above that floor and normalised,
so a fully open channel weighs its cross section.

In the beam's class, elastic scattering above `theta_cm_min_deg` (default 1°)
is a channel too: TALYS's shape plus compound elastic cross section, with θ_cm
drawn from TALYS's distribution in one random draw. TALYS tabulates it as a
ratio to Rutherford every 2°. Within each 2° cell the ratio is held at its mean
and the angle is exact Rutherford, so the weight is exactly the integral of
what is sampled. Near and above the Coulomb barrier this matters at large
angles: for 37Cl at E_cm 8.65 MeV, elastic above 60° is 0.39 of Rutherford and
above 120° 0.17; above 10° it is within 2%. Every other step is isotropic.

`[talys] min_channel_fraction` (default 0.001) and `min_channel_mb` (default
1e-4 mb, 100 nb) leave out channels too small to matter: a channel is kept
only if, at one selected strip, its cross section reaches both that fraction of
the class's, elastic excluded, and that size. Each strip's left-out cross
section is recorded.

Strips where a class cannot be simulated as TALYS describes are warned about.
`musicsim` prints each as a `musicsim WARNING:` line at startup and keeps it in
`metadata/warnings`:

- TALYS has no channel ending in the element open at the strip's energy, so
  its events do not react, or, in the beam's class, scatter only elastically;
- TALYS gives no elastic angular distribution far below the barrier (all zeros,
  from E_cm of about 0.1 MeV down), so elastic there is pure Rutherford, which
  it is in practice;
- a strip whose centre lies beyond the mean beam's range gets no table. Beam
  that straggles that far scatters elastically with Rutherford angles in the
  beam's class, and does not react in any other. Its resolved-configuration
  line says so instead of giving cross sections;
- a channel whose tables leave no state from which its chain can finish is left
  out of that strip;
- the left-out channels add up to more than `min_channel_fraction` of the
  class.

A beam meant to stop in the detector is the usual way to meet the first three.

`musicsim` only reads these tables, and stops with a message naming the
missing one. Generate them first; like `srim-cache`, it skips tables that
already exist:

```bash
./result/bin/srim-cache path/to/control.toml # when physics.stopping reads SRIM
./result/bin/talys-cache path/to/control.toml
```

`talys-cache` sets up the beam and detector exactly as `musicsim` does. It runs
TALYS for alpha on the beam nucleus at rest, at the alpha energy with the
strip's E_cm, for each strip's table not yet cached. Tables go in `talys/`
under the directory of `run.output`, beside the SRIM tables. A table's file
name ends in the MD5 of its exact TALYS input, so a change to the beam, gas,
geometry, strips, or TALYS settings gives new tables, while identical inputs
share one: the an, ap, and aa controls of one beam use the same tables. Each
table also records its input, and `musicsim` refuses a table whose recorded
input is not the strip's. Run `srim-cache` first when `physics.stopping` reads
SRIM tables, since the strip energies depend on them. Only the TALYS tables
name the class's products, so `talys-cache` then makes their SRIM tables too.

A table holds every exclusive channel TALYS gives at that energy: a
`channel <emitted> <mb>` row, then per step a `step <k> <Z> <A>` row and its
`level` and `bin` rows, then `elastic θ ratio` rows. It comes from TALYS's
`channels`, `outpopulation`, and `outdecay` output. Each channel's cross
section must match TALYS's exclusive cross sections (its section 6a) to 0.1%,
and each step's rows must sum to it. TALYS's shape elastic must match
Rutherford at the strip's E_cm to 1%, which catches a wrong frame or energy.
The decay output this needs is about 2 GB for 87Rb; it is streamed, then
deleted.

The optional `[talys]` table configures TALYS and the class, and is an error
without `talys_class`:

```toml
[talys]
alpha_omp            = 9           # default: 9 when an Atomki-V2 potential exists for the beam nucleus, else 2
keywords             = ["bins 40"] # appended to the defaults; TALYS keeps a keyword's last value
min_channel_fraction = 0.001
min_channel_mb       = 1e-4
```

`alpha_omp` 2 is McFadden-Satchler and 9 is Atomki-V2, which covers Z 26–83.
The defaults are the MUSIC analysis's shared TALYS settings: `ldmodel 1`,
`maxlevelstar 30`, `transpower 20`, `ecisstep 0.02`, `xseps`, `transeps`, and
`popeps 1.e-30`, `elow 1.e-06`, `widthfluc 20`, `widthmode 1`,
`preequilibrium n`, `bins 80`, `equidistant y`, `segment 4`. A keyword may not
set `projectile`, `element`, `mass`, `energy`, `alphaomp`, `outpopulation`,
`outangle`, `fileelastic`, `channels`, or `outdecay`. TALYS is `TALYS_BIN`,
else the one the flake built in; the potentials are in `TALYS_ATOMKI_V2_DIR`,
else the flake's. Tables are keyed by TALYS input, not by TALYS build, so
delete `talys/` after changing TALYS.

### Layers and stopping-power models

Each present window/degrader table must contain exactly one of
`thickness_mg_cm2` or `thickness_um`. A zero or negative thickness disables the
layer; a positive layer also requires `material`. Supported materials are Ti,
Al, Havar, Kapton, and Mylar. `dedx_scale` defaults to 1 and is implemented as
equivalent thickness, so it changes both mean loss and straggling. The beam and
reaction-product `dedx_scale` values differ: they multiply only that particle's
mean stopping in the gas, and its gas straggling keeps catima's variance. Gas
pressure is measured, so those scales correct the stopping model rather than
the amount of matter.

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
and `schiwietz`. `low_energy` accepts `srim_85` or `srim_95`.
`physics.straggling = false` disables both gas and entrance-layer/degrader
energy-loss fluctuations; detector noise is independent.

The effective charge only affects catima's mean stopping above 10 MeV/u. Below
that, catima uses an empirical stopping fit that does not use it. Straggling
variance uses it at every energy.

`atima14` and `global` are rejected. Both need catima's optional GLOBAL
library, which this build does not include. Without it, catima does not fail;
it silently returns an effective charge of zero (`atima14`) or a wrong one
(`global`). For Z <= 28 below 1000 MeV/u, and for any Z below 30 MeV/u, catima's
`atima14` is Winger, so the default `winger` reproduces it for every MUSIC beam.

Before 2026-10-14 the default was `atima14`. Mean stopping was correct below
10 MeV/u and too low above it: -11% at 12 MeV/u, -20% at 13.75 MeV/u, -32% at
16 MeV/u, and -53% at 20 MeV/u. Runs with any particle above 10 MeV/u in the
gas should be redone. That includes fast light reaction products, not only the
beam.

### Gas straggling and transport

Gas transport ends a step exactly at the first cell edge, detector surface,
requested vertex, time limit, or stopping point it reaches. Every traversed
fraction contributes its energy deposit. Passive matter cannot add energy.

For Vavilov-band gas steps, a `(kappa, beta^2, probability)` quantile table is
built from ROOT's `Math::VavilovAccurate` and interpolated during transport.
The table spans `kappa = 1e-3..10`; thinner steps use the lowest supported,
Landau-like Vavilov shape, while thicker steps use the Gaussian limit. The
distribution is standardized with ROOT's analytic moments and scaled to
catima's step variance. Kappa uses the same effective-charge prescription as
the configured straggling variance.

How closely the sampled moments match catima's depends on kappa:

- At `kappa >= 0.01`, sampled draws have zero mean and their variance is within
  1% of catima's. Heavy ions at the default 10 um step fall in this range (37Cl
  and 87Rb in 400-555 Torr He are at about 0.03-1).
- At `kappa = 1e-3`, the sampled variance is about 7% low, because the
  Landau-like tail beyond the last table point is cut off.
- Below `kappa = 1e-3`, the shape is frozen. A draw that would give the ion
  energy is clamped to zero loss, so the mean loss is biased high: about +50%
  for 5 MeV protons and +5% for 20 MeV alphas in 10 um steps. Light reaction
  products at the default step are in this range.

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
be checked before running. An elastic first step becomes `ground` plus
`rutherford`, any other final step `uniform`, and an intermediate step
`forced`. Legacy `Eres` becomes `noise_sigma_mev`: every `.msc` file predates
the percent-FWHM convention (June 2026), so its value is an absolute sigma in
MeV. Legacy `Kb` is treated as accelerator energy, while upstream files used
gas-surface energy, so window-enabled conversions also need manual review. Unsupported trace-database and obsolete SRIM/geometry settings
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
| `n_steps`, `event_index` | Configured reaction-step count, or a TALYS class's longest chain, and deterministic per-strip event identity. |
| `requested_strip`, `reaction_strip` | Configured vertex strip and actual successful reaction strip; the latter is `-1` when no reaction occurred. |
| `reaction_channel` | Index into `metadata/channels`: the TALYS class channel the event ran, or `0` for configured steps; `-1` when no reaction occurred. |
| `beam_energy_accel`, `beam_energy_gas`, `beam_energy_reaction`, `beam_energy_exit` | Beam kinetic-energy chain in MeV. |
| `vertex_x`, `vertex_y`, `vertex_z` | Successful reaction vertex in cm. |
| `beam_stop_x/y/z`, `beam_stop_strip`, `beam_termination` | Final beam state and termination. |
| `evap_energy[n_steps]`, `residue_energy[n_steps]` | Per-step product kinetic energies at creation. |
| `theta_cm[n_steps]`, `phi_cm[n_steps]` | Per-step sampled center-of-mass angles in degrees. |
| `residue_excitation[n_steps]`, `reaction_branch[n_steps]` | Per-step residue excitation at creation in MeV (`-2` when no reaction), and the index of the drawn branch (`0` for a step without branch tables or of a TALYS class, `-1` when no reaction). |
| `evap_theta/phi[n_steps]`, `residue_theta/phi[n_steps]` | Per-step lab angles in degrees. |
| `evap_energy_exit[n_steps]`, `evap_stop_x/y/z[n_steps]`, `evap_stop_strip[n_steps]`, `evap_termination[n_steps]` | Light-product final states. |
| `residue_energy_exit[n_steps]` | Exit-energy slots; only the surviving residue can be populated. |
| `residue_step`, `residue_stop_x/y/z`, `residue_stop_strip`, `residue_termination` | Surviving-residue identity and final state. |
| `DeadUS_dE`, `DeadDS_dE` | Deposits in unread upstream/downstream gas regions, in MeV. |

### Sentinels

Every `MC` branch is reset before each event and keeps its sentinel unless the
event produced that quantity. Two conventions run through the table: `-2`
means not applicable (the particle or quantity does not exist in this event),
and `-1` means the particle was transported but did not leave through the
downstream face. Angles use `-1000`, which lies outside every valid angle.

| Branch | Sentinel | Otherwise |
| --- | --- | --- |
| `requested_strip` | `-1`: unreacted run (`strip = -1`). | The configured vertex strip, whether or not a reaction happened. |
| `reaction_strip` | `-1`: no reaction (unreacted run, beam ended before the vertex, or chain rejected). | Equal to `requested_strip`. |
| `reaction_channel` | `-1`: no reaction. | The channel the event ran, as listed in `metadata/channels`. |
| `beam_energy_accel` | None. | This event's accelerator KE, including the `energy_fwhm` draw. The nominal value is in `metadata`. |
| `beam_energy_gas` | None. `0` means the beam stopped in the degrader or entrance window. | KE entering the gas. |
| `beam_energy_reaction` | `-2`: no reaction. | Beam KE at the vertex. |
| `beam_energy_exit` | `-2`: beam consumed by a reaction. `-1`: transported but did not exit downstream. | KE after the exit window; `0` means it stopped in the window. |
| `vertex_x/y/z` | `vertex_z = -1000`: no reaction. `vertex_x/y` are then `0`, which is also a valid coordinate, so test `vertex_z` or `reaction_strip`. | Reaction vertex in cm. |
| `beam_stop_x/y/z` | `(0, 0, -1000)`: beam consumed by a reaction. | Final beam position; the crossing point for exits. |
| `beam_stop_strip` | `-2`: beam consumed by a reaction. `-1`: stopped in a dead layer, or ended any way other than stopping. | Readout strip 0..17 where the beam stopped. |
| `beam_termination` | None; always a termination code. Reacted events report `1`. | |
| `evap_energy`, `residue_energy` | `-2`: the step did not happen. A rejected chain leaves every step at `-2`, as do the steps past a TALYS class event's chain. | Lab KE at creation. |
| `theta_cm`, `evap_theta`, `residue_theta` | `-1000`: the step did not happen. | Angle in degrees, in `[0, 180]`. |
| `phi_cm`, `evap_phi`, `residue_phi` | `-1000`: the step did not happen. | Angle in degrees, in `[-180, 180]`. |
| `evap_energy_exit` | `-2`: not transported. `-1`: transported but did not exit downstream. | KE after the exit window; neutrals keep their KE. |
| `residue_energy_exit` | `-2` in every slot except `residue_step`; intermediate residues decay at the vertex. | As `evap_energy_exit`, for the surviving residue. |
| `evap_stop_x/y/z`, `residue_stop_x/y/z` | `(0, 0, -1000)`: not transported. | Final position; the crossing point for exits. |
| `evap_stop_strip`, `residue_stop_strip` | `-2`: not transported. `-1`: stopped in a dead layer, or ended any way other than stopping. | Readout strip 0..17 where the particle stopped. |
| `evap_termination`, `residue_termination` | `0`: not transported. | Termination code. |
| `residue_step` | `-1`: no surviving residue (no reaction). | Index of the final step, whose residue was transported. |
| `DeadUS_dE`, `DeadDS_dE` | None; `0` is a real zero deposit. | Deposits in MeV. |

How the beam fields combine for each kind of event:

| Event | `beam_termination` | `beam_energy_exit` | `beam_stop_strip` | `vertex_z` / `reaction_strip` | Product fields |
| --- | --- | --- | --- | --- | --- |
| Reacted | `1` | `-2` | `-2` | vertex / strip | Filled; only the surviving residue has an exit state. |
| Reached the vertex, chain rejected | Final code of the resumed beam | KE after window if `3`, else `-1` | Strip if `2`, else `-1` | `-1000` / `-1` | All sentinels. |
| Unreacted run | Final code | KE after window if `3`, else `-1` | Strip if `2`, else `-1` | `-1000` / `-1` | All sentinels. |
| Beam stopped before the vertex | `2` | `-1` | Strip where it stopped (`-1` in a dead layer) | `-1000` / `-1` | All sentinels. |
| Beam escaped, timed out, or invalid | `4`/`5`, `6`, `7` | `-1` | `-1` | `-1000` / `-1` | All sentinels. |

Termination codes are stored in `metadata/termination_codes` and currently
mean: 0 not propagated, 1 reached the requested vertex, 2 stopped, 3 downstream
exit, 4 upstream exit, 5 side exit, 6 timeout, and 7 invalid state.

The `metadata` directory also stores the software version, the control path,
the exact control text as read at startup, each level table and TALYS table as
read at startup (`level_table_N`, with its path in `level_table_N_path`), the
resolved configuration (schema version 6), what each `reaction_channel` value
means (`channels`), the straggling model, and any startup warnings
(`warnings`, written only when there are some). Worker
threads use these copies, so a file edited mid-run cannot reach part of the
run.

The resolved configuration records every physics default actually applied
(`z_effective`, `low_energy`, `straggling_z_effective`, `sim_step_cm`,
`max_time_ns`, `beam_energy_fwhm_mev`, ...), each reaction step's branches
with their models, weights, and level-table paths and MD5s, a TALYS class's
settings and channels (`reaction_channel_N`), each strip's beam energy, E_cm,
channel cross sections in channel order, left-out cross section, and TALYS
table path and MD5 (`talys_strip_N`), and
`git_revision` (a `-dirty` suffix means uncommitted changes),
`catima_revision`, and `catima_global`. Defaults can change between builds, so
a control file alone does not determine the physics.

An update run appends only if both the control text and the resolved
configuration (except `total_entries`) match the file. A file is therefore
never extended by a different build, with different defaults, or with an
edited level table. Schema 6 added the `reaction_channel` truth branch, so an
update run also refuses any older file; regenerate it instead.

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
- `forced` and `uniform` excitation are continuum models. `talys` samples
  discrete level energies and continuum bins from a precomputed table, but
  carries no spin/parity, widths, gamma decay, or lifetimes. A `level_table`
  is computed at one beam energy, and a TALYS class's tables at each strip's
  centre, not per vertex. Branch weights are fixed user inputs.
- A TALYS class follows TALYS's cascade only through its tables. A step's
  table is the residue's population right after that emission. If the residue
  gamma-decays before its next emission, that step starts from the pre-gamma
  excitation, so the gamma's energy goes into the next emission's kinematics.
  Only the first gamma of a gamma-first channel is emitted as a particle. Each
  step's excitation is drawn from its own table, cut by the energy available
  and the floor, and is otherwise uncorrelated with the previous state. TALYS's
  spins and parities are used only to build the tables, and its non-elastic
  angles are replaced by isotropic ones.
- Within a strip, a class channel's energy dependence comes only from its
  states closing as the vertex energy falls below the strip centre, as for
  weighted branches. Near a channel's threshold this matters: for 87Rb at
  E_cm 3.9 MeV, (a,n) is closed over the downstream third of the strip, and
  (a,g) runs in 35% of the yttrium class's events against 18% of its cross
  section at the centre.
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
