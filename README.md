# MUSIC simulator

Fork of the ANL MUSIC simulator (originally by Daniel Santiago-Gonzalez).
Stopping powers are computed on the fly via [catima](https://github.com/hrosiak/catima);
the anode geometry is hardcoded (20 strips × 2 columns).

## Build

Dependencies are managed by Nix.

```bash
nix build
./result/bin/musicsim Examples/37Cl_alpha_n/37Cl_alpha_n.msc
```

For iterative development:

```bash
nix develop
make
./musicsim Examples/37Cl_alpha_n/37Cl_alpha_n.msc
```

## Control-file keys

Stopping-power inputs are now specified per-run rather than per-particle:

| Key                 | Meaning                                                          |
| ---                 | ---                                                              |
| `Gas`               | Fill gas: `4He`, `3He`, `Ar`, `CF4`, `CH4`, `P10`, `iC4H10`      |
| `Pressure`          | Gas pressure in Torr                                             |
| `Temperature`       | Gas temperature in K                                             |
| `BeamEnergy`        | Beam KE in MeV at the accelerator (before the entrance window)   |
| `KbFWHM`            | Beam energy FWHM in MeV at the accelerator                       |
| `EntranceMaterial`  | Entrance window: `Ti` (default), `Havar`, `Kapton`, or `Mylar`   |
| `EntranceThickness` | Entrance window areal density in mg/cm²                          |
| `ExitMaterial`      | Exit window material (same vocabulary, plus `Al`/`Aluminum`)     |
| `ExitThickness`     | Exit window areal density in mg/cm²                              |
| `DegraderMaterial`  | Optional bulk degrader between the accelerator and the entrance window. Empty (default) = no degrader. |
| `DegraderLength`    | Degrader length along the beam axis in μm                        |

The per-event beam-energy chain is:
`BeamEnergy ± KbFWHM` (accelerator) → degrader (if any) → entrance window →
`Kbi` (gas surface). Energy straggling is sampled per event from catima's
`sigma_E` at both the degrader and the entrance window — the accelerator
spread (`KbFWHM`) and straggling add in quadrature naturally because each is
applied as an independent Gaussian draw.

Particles that exit the downstream face of the gas are propagated through the
exit window before their energy is recorded as `Kbeam_exit` / `Kh_exit[]` /
`Kl_exit[]` on the `MC` tree.

The old `Kb` key (beam KE at the gas surface) is no longer accepted — use
`BeamEnergy` (accelerator KE) and the entrance-window parameters instead.

### Disabling layers

Any of the four physical layers can be turned off independently — the
simulator prints a warning at startup and skips both the mean dE/dx and the
straggling for that layer. Other layers still apply normally.

| Disable                         | How                                       |
| ---                             | ---                                       |
| Gas energy loss + straggling    | `Pressure 0` (or any non-positive value)  |
| Entrance window                 | `EntranceThickness -1`                    |
| Exit window                     | `ExitThickness -1`                        |
| Degrader                        | Omit `DegraderMaterial`, or `DegraderLength -1` |

Useful for sanity checks (e.g. "what is the beam energy at the gas surface
with no degrader?") and for debugging specific contributions in isolation.

### Multi-threading

`Threads N` in the ctrl file fans the event loop out across `N` worker threads
using `std::async`. Each worker is a separate `MUSIC_Simulator` instance with
its own RNG seed and per-worker ROOT output; the master pre-warms catima's
internal cache and then merges the per-worker output files with `TFileMerger`
into the configured `FileName` at the end of the run. Workers force
`Update=0`/`Wait=0`, so visualization is only available in single-threaded
mode (the default, `Threads 1`).

Per-particle `SRIMbeam`, `SRIMres*`, `SRIMevap*`, `AnodeGeom`, and `SRIMdir` keys
are ignored (with a warning) — catima computes dE/dx from the gas composition,
and the anode geometry is no longer file-driven. The per-particle
`dEdxScale*` knobs still apply, multiplying the catima energy loss.

Supported particle naming (`AEl`, e.g. `37Cl`) and the rest of the control-file
schema are unchanged.

## Output tree layout

Each run produces one ROOT file with two trees, mirroring the format produced by
the upstream `EventBuilderNearestGrid` analysis pipeline:

- **`event`** — detector-level branches. `LeftdE[18]`, `RightdE[18]`,
  `TotaldE[18]` (Float, MeV); `AllTimestamps[36]` (ULong64), `AllFlags[36]`
  (UInt, always 0 in sim), `Hits[36]` (Int, 0/1 above noise threshold);
  `Cathode`, `Grid` (Float, MeV); `IsComplete` (Bool, same heuristic as the
  upstream event builder). Energies are MeV truth, so analysis macros should
  calibrate data to MeV rather than rescaling sim to ADC.
- **`MC`** — truth-only branches, friended to `event`: reaction strip, beam
  kinematics (`BeamEnergyAccel` at the accelerator, `Kbi` at the gas surface,
  `Kbr` at the reaction point, `Kbeam_exit` after the exit window),
  light/heavy four-vector components, vertex and exit positions, and per-evap
  exit energies (`Kl_exit[]`, `Kh_exit[]`). Sentinel values on exit-energy
  branches: `-1` = particle stopped in the gas, `-2` = no such particle for
  this event (e.g. `Kbeam_exit` on a reacted event). A literal `0` means the
  particle reached the exit window but stopped inside it.

`Hits[k]` triggers when the corresponding channel energy exceeds
`max(0.02 MeV, 3·Eres)` so that pure-noise cells are not counted.

## A note on AI-assisted development

Parts of this fork — the catima migration, std::async multi-threading,
window / degrader plumbing, and most boilerplate around them — were
written with the help of [Claude Code](https://claude.ai/claude-code).
All changes to the simulation physics (kinematics, energy-loss chain,
straggling sampling) are human-reviewed and approved before being
committed.
