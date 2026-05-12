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
| `ExitMaterial`      | Exit window material (same vocabulary)                           |
| `ExitThickness`     | Exit window areal density in mg/cm²                              |

The simulator propagates `BeamEnergy` through the entrance window with catima at
startup to obtain `Kbi` (KE at the gas surface). Particles that exit the
downstream face of the gas are propagated through the exit window before their
energy is recorded as `Kbeam_exit` / `Kh_exit[]` / `Kl_exit[]` on the `MC` tree.

The old `Kb` key (beam KE at the gas surface) is no longer accepted — use
`BeamEnergy` (accelerator KE) and the entrance-window parameters instead.

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
