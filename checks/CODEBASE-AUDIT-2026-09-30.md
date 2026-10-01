# Remix-MUSIC-Sim codebase audit — 2026-09-30

## Executive summary

**Disposition: not ready for scientifically trusted production runs.** The
project builds and a small multi-threaded smoke run completes, but the audit
found defects that can silently produce impossible events, biased energy
deposits, incorrect truth metadata, and configuration-dependent undefined
behaviour. Multi-threaded runs additionally have C++ data races in shared
physics configuration and potentially in catima's process-wide cache.

The most consequential findings are:

1. An energetically forbidden primary reaction can be reported as accepted;
   the incident beam then disappears at the sampled vertex.
2. Multi-threaded workers mutate process-global stopping and straggling
   configuration while other workers can be building physics tables.
3. The default residue-excitation policy is explicitly described in the source
   as wrong for single-step reactions, yet every shipped reacting example
   inherits it.
4. The propagator drops the final energy-loss step when a particle stops and
   drops boundary-crossing fractions, breaking event-level energy accounting.
5. A reacted configuration with no reaction steps dereferences uninitialized
   pointers and crashes under ASan/UBSan.
6. Multi-step `theta_cm`/`phi_cm` truth arrays store the last sampled angle in
   every element rather than the angle used for each step.

Finding count:

| Severity | Count | Meaning |
| --- | ---: | --- |
| Critical | 3 | Can silently invalidate ordinary scientific output or invoke undefined behaviour in a supported production mode. |
| High | 15 | Material correctness, reliability, or operability failure in a plausible run. |
| Medium | 18 | Conditional defect, unsafe edge case, misleading interface, or important validation gap. |
| Low | 4 | Maintainability/style debt with limited immediate runtime impact. |

The immediate recommendation is to stop producing or publishing reacted and
multi-threaded samples until F-01 through F-18 are addressed and regression
tests are in place. Existing samples should be triaged using the
[affected-results guidance](#affected-results-and-rerun-guidance).

## Audit scope and evidence

- Revision: `c9cc7e8435881ab48955685145f3b629a80782ac`.
- Audit date: 2026-09-30.
- Scope: all tracked C++ sources and headers, build/packaging files, the TOML
  loader, all 29 shipped `ControlExamples`, `basic.toml`, the legacy converter,
  README, and existing review notes.
- A pre-existing local edit to `flake.nix` was present. It was inspected as part
  of the working tree but was not modified or attributed to this audit.
- No simulator source, configuration, or generated data was changed by the
  audit. Temporary controls, binaries, and ROOT outputs were kept under `/tmp`.

Evidence labels used below:

| Label | Meaning |
| --- | --- |
| **Reproduced** | Observed by executing the audited revision or its sanitizer build. |
| **Static** | Follows directly from source/control flow, compiler diagnostics, or static analysis. |
| **Model** | Physics-model inconsistency or an assumption that needs domain validation. |
| **Documentation** | Shipped documentation/configuration disagrees with implemented behaviour. |

This was a code audit, not an experimental validation campaign. It checked
conservation logic, limiting behaviour, units, transport accounting, and model
consistency, but did not compare large Monte Carlo ensembles against measured
MUSIC data, SRIM output, LISE++, Geant4, or an independent Vavilov
implementation.

## Severity definitions

- **Critical:** silent corruption of normal scientific results, impossible
  events accepted as real, or concurrency-related undefined behaviour in a
  documented production mode.
- **High:** crash/hang, significant physics bias, data loss, or materially wrong
  shipped configuration under plausible use.
- **Medium:** conditional correctness issue, missing validation, misleading
  output/documentation, or maintainability problem likely to cause a future
  defect.
- **Low:** localized style, dead-code, or robustness issue with little current
  effect.

## Findings register

| ID | Severity | Evidence | Area | Finding |
| --- | --- | --- | --- | --- |
| F-01 | **Critical** | Reproduced | Kinematics | Energetically forbidden primary reactions can return “allowed,” consume the beam, and be written with the requested reaction strip. |
| F-02 | **Critical** | Static | Concurrency/physics | Workers race while writing process-global catima and straggling configuration; cache pre-warming is not sufficient for every supported configuration. |
| F-03 | **Critical** | Model + static | Kinematics/defaults | The default `forced` residue-excitation model locks away about 83% of available energy and is known to be wrong for single-step reactions; all shipped examples inherit it. |
| F-04 | **High** | Reproduced | Kinematics/configuration | A reacted run with zero reaction steps dereferences uninitialized `EvaP[0]`/`EvaR[0]`; ASan/UBSan reports invalid access and the process segfaults. |
| F-05 | **High** | Reproduced | Transport | The stopping step is discarded, so deposited energy plus residual state is not energy-conserving. |
| F-06 | **High** | Static | Transport/kinematics | The caller does not verify that the beam reached the sampled reaction vertex; a stopped/escaped/timed-out beam can react at an unreachable downstream point. |
| F-07 | **High** | Static | Decay chains | Failure after an earlier decay step is not atomic: the code can resume the original beam while retaining truth for an earlier ejectile. |
| F-08 | **High** | Reproduced | Output/kinematics | Every multi-step CM-angle slot receives the final step's angle. |
| F-09 | **High** | Static | Transport | A step that crosses a vertex or detector boundary is discarded rather than shortened to the boundary. |
| F-10 | **High** | Static | Configuration/transport | `sim_step = 0` causes an exact infinite loop; a negative value reverses the normal event transport. |
| F-11 | **High** | Model + static | Straggling | The low-κ Vavilov limit is replaced by a Gaussian even though the low-κ limit is Landau-like. |
| F-12 | **High** | Model + documentation | Shipped examples | The shipped elastic-scattering controls inherit both forced excitation and isotropic scattering, so they do not model elastic Rutherford scattering. |
| F-13 | **High** | Documentation + static | Shipped examples/noise | 24 scalar-noise examples and all three per-channel tables describe absolute MeV σ while the implementation interprets values as percent FWHM. |
| F-14 | **High** | Documentation + static | Shipped examples/beam | 19 of 29 shipped examples omit beam energy and silently run at the compiled 100 MeV default. |
| F-15 | **High** | Reproduced | CLI/error handling | Control-file parse failure does not stop execution, normal batch success returns status 1, and an argument-count error returns status 0. |
| F-16 | **High** | Static | Configuration/physics | Reaction steps are not checked for completeness, nuclide validity, chain continuity, or charge/mass-number conservation. |
| F-17 | **High** | Static | Multi-threaded output | Worker return values are ignored and the final merge always uses `RECREATE`, defeating `file_opt = "update"` and risking partial output. |
| F-18 | **High** | Static | Resources | Pervasive raw ownership and unconditional 50,000-point particle traces cause large, cumulative leaks; the memory guard is disabled. |
| F-19 | **Medium** | Model + static | Decay thresholds | `minEx` rejects the equality case of an exothermic ground-state secondary decay even though the decay is energetically open. |
| F-20 | **Medium** | Static | Window/degrader physics | Gaussian layer straggling is floored at zero but not capped at the incoming energy, so a window or degrader can add energy. |
| F-21 | **Medium** | Model + static | Exit window | All downstream particles traverse the exit foil at normal thickness, without the planar `1/|cos(theta)|` path correction. |
| F-22 | **Medium** | Static | Output semantics | Side/upstream escape and max-time termination are all recorded as “stopped in gas.” |
| F-23 | **Medium** | Model + static | Straggling | The quantile table clamps both probability tails and estimates moments from point samples; exact first-two-moment preservation is not achieved. |
| F-24 | **Medium** | Static | Straggling | Convolution components below the κ table minimum are dropped; some valid total κ values therefore produce exactly zero fluctuation. |
| F-25 | **Medium** | Model | Straggling | Vavilov κ uses bare nuclear `Z²` while the width comes from a separately selected effective-charge model. |
| F-26 | **Medium** | Static | Configuration | Most numeric domains, finiteness, enums, ranges, unknown keys, and wrong value types are not validated; a layer may specify neither documented thickness key. |
| F-27 | **Medium** | Static | ROOT I/O | File, directory, tree, and update-schema results are dereferenced without checking for zombie/null/incompatible objects. |
| F-28 | **Medium** | Static analyzer | Nuclide loaders | Optional AME/NUBASE loaders have unchecked indexed writes, a non-standard VLA, short-line hazards, leaks, and an uninitialized spin flag. |
| F-29 | **Medium** | Static | Particle state | `GetEexc()` resamples an excitation state on every read, so repeated `GetKE()` calls can describe different states for one particle. |
| F-30 | **Medium** | Static | Reproducibility | Single-thread runs auto-seed nondeterministically; seeds/config provenance are not persisted, and MT streams depend on worker count. |
| F-31 | **Medium** | Reproduced | Shipped examples | `17F_alpha_alpha.toml` specifies both strip-selection forms and is rejected by the parser. |
| F-32 | **Medium** | Documentation | README/schema | Required inputs, physics switches, layer behaviour, and output branch names are incomplete or stale. |
| F-33 | **Medium** | Static + documentation | Run semantics | `n_events` is applied per selected strip, not to the whole run, without being documented as such. |
| F-34 | **Medium** | Static | Tool security | `srim-cache` builds a shell command from configuration/environment strings and passes it to `system()`. |
| F-35 | **Medium** | Reproduced + static | Converter/tooling | The documented converter cannot run in the Nix dev shell and has silent/unsafe conversion behaviours. |
| F-36 | **Medium** | Static | Quality controls | There is no automated test suite or CI, despite physics-critical code and recently changing output/schema behaviour. |
| F-37 | **Low** | Static analyzer | Maintainability | Detector-response scratch arrays are fixed at length 3 and depend implicitly on `AnodeCols == 2`. |
| F-38 | **Low** | Static | Dormant APIs | Public legacy routines contain null/out-of-bounds and divide-by-zero paths if called. |
| F-39 | **Low** | Static | Dead code | Parsed/no-op/unused interfaces obscure the active execution path. |
| F-40 | **Low** | Compiler/formatters | Style | The strict build emits 188 warnings, one C++ formatting site fails, and `basic.toml` is not Taplo-formatted. |

## Critical findings

### F-01 — Forbidden primary reaction is accepted

**Evidence:** Reproduced.

When `EneAvail < 0`, the reaction loop marks the current and subsequent
products `DoNotPropagate` and breaks, but never clears `ReactionAllowed`
([`Kinematics.cpp:108–121`](../src/Kinematics.cpp#L108-L121)). The event loop
therefore enters the “reaction allowed” branch
([`EventLoop.cpp:193–244`](../src/EventLoop.cpp#L193-L244)) instead of resuming
the incident beam. Detector response also records the requested strip
unconditionally ([`Propagation.cpp:211–218`](../src/Propagation.cpp#L211-L218)).

A deterministic one-event test used gas/window loss disabled and
`4He + 4He -> 1H + 7Li` at 1 MeV. The output contained:

```text
reaction_strip       = 0
beam_energy_reaction = 1 MeV
beam_energy_exit     = -1
evap_energy[0]       = -2
residue_energy[0]    = -2
```

The products were not realized, but the beam was not transported beyond the
vertex and was classified as stopped. This creates an impossible event and can
silently distort both traces and reaction yields.

**Remediation:** Return a structured reaction result. A failed primary step
must set `allowed = false`, clear all product state, leave the beam state at the
vertex, and let the caller continue it. Record `reaction_strip = -1` (or a
separate attempted-vertex field) when no reaction occurred. Add an invariant
test for every failure branch.

### F-02 — Multi-threaded physics state has data races

**Evidence:** Static; a two-worker smoke run succeeding does not disprove the
race.

Every worker re-runs `loadCtrlFile()`. It writes
`catima::default_config` ([`ControlFile.cpp:58–65`](../src/ControlFile.cpp#L58-L65),
[`ControlFile.cpp:247–264`](../src/ControlFile.cpp#L247-L264)),
`music::gStragglingConfig`, and `music::gStragglingEnabled`
([`ControlFile.cpp:266–300`](../src/ControlFile.cpp#L266-L300)). Only the SRIM
selector/path globals are guarded to the master
([`ControlFile.cpp:412–432`](../src/ControlFile.cpp#L412-L432)). Once an early
worker reaches `SetupRun()`, it can call catima and build `EnergyLoss` tables
while later workers are still mutating the globals
([`Simulator.cpp:395–415`](../src/Simulator.cpp#L395-L415)). This is a C++ data
race even when each worker happens to assign the same value.

`PreWarmCatima()` acknowledges that catima's cache is process-wide and unsafe
for concurrent writes ([`Simulator.cpp:128–176`](../src/Simulator.cpp#L128-L176)).
The installed catima revision has only 60 cache entries. The supported maximum
can require as many as 23 ions (beam, target, compound, and two particles for
each of ten steps), multiple materials, and two configurations. Pre-warming can
therefore evict earlier entries; a later cache miss can race in the workers.

Every worker also truncates/writes the same `musicsim.log` and
`energetics.log` ([`Simulator.cpp:101–120`](../src/Simulator.cpp#L101-L120),
[`Simulator.cpp:296–312`](../src/Simulator.cpp#L296-L312)). Those streams can
lose or interleave diagnostics.

**Remediation:** Parse once into an immutable, instance-owned physics config;
pass explicit `catima::Config` values into every calculation; build all
per-particle tables serially before launching workers or give each worker a
private/thread-safe cache. Use per-worker log files or synchronized logging.
Run ThreadSanitizer where ROOT/catima permit it, plus a high-ion-count stress
test that forces cache churn.

### F-03 — The default reaction excitation model is knowingly inappropriate

**Evidence:** Model/static.

The compiled default is `residueExc = 0`, “forced.” The source itself says this
is wrong for a single-step reaction and explains the consequence
([`Simulator.hpp:351–365`](../include/Simulator.hpp#L351-L365)). The
implementation samples residue excitation uniformly over
`[2 Eavail / 3, Eavail]`
([`Kinematics.cpp:128–140`](../src/Kinematics.cpp#L128-L140)); its mean removes
`5/6`, or about 83%, of the energy available to outgoing kinetic energy.

None of the shipped controls sets `reaction.residue_excitation`, and the README
does not document the option. Thus every shipped reacting example silently
uses this policy, including single-step channels.

**Impact:** Product energies, ranges, stopping strips, and detector traces are
systematically shifted. This is not merely a sampling choice: for a nominally
elastic channel it turns an elastic final nucleus into a highly excited one.

**Remediation:** Require an explicit excitation model for reacted runs, or make
`ground` the safe single-step default. For decay chains, accept a per-step
excitation distribution or level scheme rather than using one global heuristic.
Record the chosen model in output metadata and regenerate affected samples.

## High-severity findings

### F-04 — Zero-step reacted configurations invoke undefined behaviour

The parser deliberately sets `NumEvapPart = 0` when no `[[reaction.step]]`
array exists ([`ControlFile.cpp:372–400`](../src/ControlFile.cpp#L372-L400)). A
reacted run nevertheless constructs its reaction label from `EvaP[0]` and
`EvaR[0]` ([`Kinematics.cpp:74–84`](../src/Kinematics.cpp#L74-L84)). Those array
slots were allocated but never initialized.

An ASan/UBSan run reported misaligned accesses through `0xbebebebe...` at
`Kinematics.cpp:82–83` and terminated with a segmentation violation (exit
129). A non-sanitized run happened to survive once, which is typical of
undefined behaviour and not evidence of safety.

Require at least one complete step when `strip >= 0`; explicitly allow zero
steps only for `strip = -1`. Initialize pointer arrays to null regardless.

### F-05 — Terminal stopping energy is not deposited

The propagator calculates `Kf`, then immediately breaks when `Kf < 1 keV`
before adding `|Ki-Kf|` to the detector response or updating energy/time/state
([`Propagation.cpp:131–172`](../src/Propagation.cpp#L131-L172)). After the loop,
it stores the tentative new position but the previous step's energy and
momentum ([`Propagation.cpp:193–205`](../src/Propagation.cpp#L193-L205)).

In a deterministic, straggling-disabled 2 MeV alpha test:

```text
input KE                            2.000000000 MeV
DeadUS_dE                           1.441781163 MeV
sum of readout deposits (Cathode)   0.553954899 MeV
accounted total                     1.995736062 MeV
missing                             0.004263938 MeV
```

The stored stop position advanced to the rejected step endpoint while its
momentum remained from the preceding state. Deposit the remaining kinetic
energy in the traversed cell (or solve the exact stopping fraction), set the
particle to zero kinetic energy at the physical stop point, and return an
explicit `Stopped` status.

### F-06 — Reactions can occur beyond the beam range

`PropagateParticle()` returns the same success value after stopping, escaping,
timing out, or reaching `endZ` ([`Propagation.cpp:60–93`](../src/Propagation.cpp#L60-L93),
[`Propagation.cpp:146–208`](../src/Propagation.cpp#L146-L208)). The event loop
does not inspect the reached position or a termination reason. It always reads
the residual beam energy and calls `SetReactionKinematics()` at the originally
sampled vertex ([`EventLoop.cpp:169–193`](../src/EventLoop.cpp#L169-L193)), which
rewrites the beam position to that vertex.

A low-energy beam that stops before a downstream selected strip can therefore
react there using the stale pre-stop kinetic energy. Return a termination enum
and final state, and only construct a reaction if the propagator reached the
vertex within tolerance.

### F-07 — Later-step failure does not roll back an earlier reaction

For a later step with `EneAvail <= minEx`, the code can first mark the previous
residue non-propagating, set `ReactionAllowed = 0`, and park only the current
products ([`Kinematics.cpp:122–142`](../src/Kinematics.cpp#L122-L142),
[`Kinematics.cpp:232–241`](../src/Kinematics.cpp#L232-L241)). Earlier ejectiles
remain marked as real. The caller sees a single false result and resumes the
original beam instead of transporting any products
([`EventLoop.cpp:216–255`](../src/EventLoop.cpp#L216-L255)). Later output logic
can still treat the earlier ejectile as transported truth.

This mixes mutually exclusive histories: an intact beam and an already emitted
particle. Decide whether a chain is all-or-nothing or may terminate at the last
allowed residue, then implement that policy explicitly. A per-step result is
preferable to one chain-wide boolean.

### F-08 — Multi-step CM truth angles are overwritten

Angles are sampled/mutated inside the reaction loop
([`Kinematics.cpp:146–165`](../src/Kinematics.cpp#L146-L165)), but after the loop
the current final values are copied to every array element
([`Kinematics.cpp:244–261`](../src/Kinematics.cpp#L244-L261)).

In a deterministic two-neutron test, the event log showed independent step
angles of approximately 115.52° and 96.08°, while the ROOT tree stored
`theta_cm = [96.08495, 96.08495]`; `phi_cm` is affected identically. Save each
step's angles at the point they are drawn and do not refill them after the loop.

### F-09 — Boundary-crossing path length and loss are discarded

The next full endpoint is calculated before boundary checking. If it lies past
the reaction vertex, downstream face, upstream face, or a side wall, the loop
breaks without transporting the fractional distance to the surface
([`Propagation.cpp:70–93`](../src/Propagation.cpp#L70-L93)). The post-loop state
is then placed at that overshot endpoint. Consequences include:

- reaction energy and TOF correspond to a point up to one `sim_step` before the
  vertex while kinematics are placed at the vertex;
- exit gas loss is low by up to one full step;
- exit/stop coordinates lie outside the physical boundary; and
- results have a first-order dependence on step size near every surface.

Intersect the segment with the first crossed surface/endZ, transport exactly
that fractional length, and store the surface point.

### F-10 — Nonpositive `sim_step` is unsafe in the normal run path

`sim_step` is read without validation
([`ControlFile.cpp:402–410`](../src/ControlFile.cpp#L402-L410)). With zero step,
velocity can be positive but `Dt`, position, time, and energy all remain
unchanged, so the loop at [`Propagation.cpp:60–78`](../src/Propagation.cpp#L60-L78)
never terminates. A negative step takes the inverse-energy branch and moves
ordinary event particles backward; negative steps are intended only by a
legacy angular-scan caller ([`EventLoop.cpp:427–430`](../src/EventLoop.cpp#L427-L430)).

Require finite `sim_step > 0` in user controls and make propagation direction
an explicit API value rather than encoding it in the sign of a length.

### F-11 — The low-κ Vavilov fallback has the wrong limiting shape

`SampleStandardized()` uses a standard Gaussian for both
`κ <= 10^-3` and `κ >= 10`
([`VavilovSampler.cpp:123–131`](../src/VavilovSampler.cpp#L123-L131)). ROOT's
Vavilov documentation identifies the small-κ limit as Landau-like (commonly
`κ < 0.01`) and the large-κ limit as Gaussian. The source/README claim that a
Gaussian matches both limits
([`VavilovSampler.hpp:23–27`](../include/VavilovSampler.hpp#L23-L27),
[`README.md:146–169`](../README.md#L146-L169)).

This replaces the long high-loss tail with a symmetric distribution precisely
in the thin-absorber regime. Use ROOT's Landau/Vavilov path below the supported
table band, or extend and validate the table into that regime.

### F-12 — Shipped “elastic” runs are not configured as elastic scattering

The 37Cl alpha channel declares `4He + 37Cl` products
([`37Cl_alpha_alpha.toml:73–81`](../ControlExamples/37Cl_alpha_n/37Cl_alpha_alpha.toml#L73-L81))
but does not select either ground-state residue excitation or Rutherford
angles. It therefore inherits forced excitation and isotropic CM sampling. The
source notes that isotropic sampling can make the heavy product lose roughly
19% where a forward-peaked Rutherford event loses roughly 2%
([`Simulator.hpp:366–379`](../include/Simulator.hpp#L366-L379)).

Make elastic controls explicit with `residue_excitation = "ground"` and the
validated angular law. Also note that the current single global angular option
is applied to secondary decay steps too
([`Kinematics.cpp:146–164`](../src/Kinematics.cpp#L146-L164)); distributions
should be per step so Rutherford scattering is not reused for evaporation.

### F-13 — Shipped detector-resolution values use the wrong units

The implementation treats every `eres` value as **relative percent FWHM** and
converts it to an event-dependent Gaussian σ
([`Propagation.cpp:240–247`](../src/Propagation.cpp#L240-L247)). Of 29 shipped
examples, 26 use scalar `eres`; 24 of those comments describe it as an absolute
noise in MeV and/or σ (for example
[`17F_alpha_alpha.toml:9–15`](../ControlExamples/17F_alpha/17F_alpha_alpha.toml#L9-L15)).
All three per-channel tables also say `σ [MeV]`. The active `0.05` values in
[`37Cl_beam.toml:34–71`](../ControlExamples/37Cl_alpha_n/37Cl_beam.toml#L34-L71)
are consequently interpreted as 0.05% FWHM; for a 1 MeV deposit this is
approximately 0.000212 MeV σ, not 0.05 MeV σ.

Choose one schema and migrate every control. Prefer names carrying the unit,
such as `resolution_fwhm_percent` or `noise_sigma_mev`, and reject the legacy
ambiguous key after a transition period.

### F-14 — Most shipped controls silently use a 100 MeV beam

Nineteen of the 29 tracked `ControlExamples` omit `beam.energy`. The loader
retains the compiled default of 100 MeV
([`ControlFile.cpp:475–495`](../src/ControlFile.cpp#L475-L495)); representative
older controls specify only species and spread
([`17F_alpha_alpha.toml:1–4`](../ControlExamples/17F_alpha/17F_alpha_alpha.toml#L1-L4)).
Beam energy is not a harmless UI default: it sets reaction thresholds, product
kinematics, ranges, and detector response.

Require beam energy for every run and add it to every example. If 100 MeV was
actually intended, make that value explicit so the control is auditable.

### F-15 — CLI errors and process status are inverted

`main()` logs a failed `loadCtrlFile()` but continues into `run()`
([`main.cpp:41–54`](../src/main.cpp#L41-L54)). A missing-file test therefore
continued with placeholder/default configuration and, in the headless audit
environment, segfaulted while creating the default GUI (exit 129). Independently
of that environment-specific crash, continuing after parse failure is wrong.

Normal batch completion returns 1, while too many arguments return 0
([`main.cpp:36–39`](../src/main.cpp#L36-L39),
[`main.cpp:53–62`](../src/main.cpp#L53-L62)). Shell scripts and schedulers thus
see success as failure and a usage error as success.

Return `EXIT_FAILURE` immediately on argument/load/run failures and
`EXIT_SUCCESS` after a successful batch. Make interactive-loop entry a distinct
mode rather than overloading `Simulator::run()`'s integer return.

### F-16 — Reaction definitions are not structurally or physically validated

Each array element increments the step count even if `evap` or `res` is absent,
and missing/wrong-type fields retain placeholders
([`ControlFile.cpp:372–400`](../src/ControlFile.cpp#L372-L400)). More than ten
steps are silently truncated after a warning. There is no check that:

- beam + target matches the declared compound;
- each step conserves mass number and charge;
- a step's parent equals the previous residue;
- all named nuclides exist before setup; or
- the chain has a valid number of complete steps for the selected run mode.

Kinematics then trusts these masses and charges
([`Kinematics.cpp:82–92`](../src/Kinematics.cpp#L82-L92)). Invalid nuclear
reactions can therefore be simulated if their masses happen to pass the energy
test. Validate the whole reaction graph before allocating physics objects and
print a canonical reaction/Q-value summary for user confirmation.

### F-17 — Multi-threaded completion and update semantics are unsafe

Worker lambdas reduce all successful `run()` results to 1, but the parent
ignores every returned value from `future::get()`
([`Simulator.cpp:395–421`](../src/Simulator.cpp#L395-L421)). It then ignores
individual `AddFile()` results and attempts a merge. The final output is always
opened with `RECREATE`
([`Simulator.cpp:423–435`](../src/Simulator.cpp#L423-L435)), regardless of the
requested `file_opt`.

Thus a failed worker can yield a partial/failed merge without a precise cause,
and `file_opt = "update"` can overwrite the existing final file. Propagate
structured worker status, abort before merge on any failure, validate expected
entry counts, and either implement atomic update semantics or reject `update`
when `threads > 1`.

### F-18 — Ownership leaks make supported runs memory-heavy

`Simulator`, `Particle`, and `NuclideFinder` allocate extensively but declare no
destructors. Every `Particle` unconditionally allocates five arrays of 50,000
floats—1,000,000 bytes (0.954 MiB)—even when trajectory saving is false
([`Particle.cpp:3–40`](../src/Particle.cpp#L3-L40),
[`Particle.hpp:64–96`](../include/Particle.hpp#L64-L96)). Setup can create 23
particles at the ten-step limit, and every strip simulation leaks two more
particles ([`EventLoop.cpp:34–44`](../src/EventLoop.cpp#L34-L44)). With 32
workers and a multi-strip run, trace buffers alone can reach hundreds of MB or
more than 1 GiB.

The intended memory guard returns immediately because the member `gSystem` is
initialized to null and normal CLI setup never assigns it
([`Simulator.cpp:74–93`](../src/Simulator.cpp#L74-L93)). `Chroma` is another
analyzer-confirmed leak ([`Output.cpp:153–188`](../src/Output.cpp#L153-L188)).

Use RAII (`std::vector`, `std::array`, `std::unique_ptr` where ROOT ownership
allows it), allocate traces only when requested, remove the unused `BeamInit`,
and make memory checks functional or remove the false safeguard.

## Medium-severity findings

### Physics and transport

#### F-19 — `minEx` rejects an open exothermic equality case

`SetupRun()` sets `minEx = -Q0` only when `Q0 < 0`
([`Simulator.cpp:273–293`](../src/Simulator.cpp#L273-L293)), then kinematics
requires `EneAvail > minEx`, not `>=`
([`Kinematics.cpp:128–142`](../src/Kinematics.cpp#L128-L142)). For an exothermic
secondary decay of a ground-state parent, `EneAvail == -Q0`; the decay is open
but this gate rejects it. The direct invariant-mass availability test already
contains the physical threshold. Remove the redundant gate or derive a single
threshold consistently, with a documented numerical tolerance.

#### F-20 — Gaussian window/degrader smearing can increase energy

The sampled material output is `mean Eout + Gaussian(0, sigma_E)` and is only
floored at zero ([`Materials.cpp:197–216`](../src/Materials.cpp#L197-L216)). A
positive fluctuation can make `Eout > Ein`, so the beam gains energy across a
passive layer. Bound the sampled loss to `[0, Ein]`, or sample a positive-loss
distribution whose support is physical, and quantify the bias introduced by
truncation.

#### F-21 — Exit-window loss ignores trajectory angle

Any particle with `z >= AnodeDepth` is passed through the same material object
([`Propagation.cpp:317–340`](../src/Propagation.cpp#L317-L340)). That material's
thickness is the normal beam-axis thickness, while an oblique crossing of a
planar foil has path length `t/|cos(theta)|`. This underestimates exit loss for
non-forward products. Scale areal/linear thickness by the actual surface
incidence angle and define behaviour for edge exits.

#### F-22 — Escape and timeout are mislabeled as a gas stop

The transport loop exits for downstream, upstream, or lateral boundary crossing
using one path, and can also end on `MaxTime`
([`Propagation.cpp:60–93`](../src/Propagation.cpp#L60-L93)). Output recognizes
only `z >= AnodeDepth` as an exit; everything else returns `-1`, documented as
“stopped inside the gas”
([`Propagation.cpp:313–342`](../src/Propagation.cpp#L313-L342)). Store an
explicit termination reason (`stopped`, `downstream_exit`, `side_exit`,
`upstream_exit`, `timeout`, `not_propagated`) and do not overload an energy
sentinel.

#### F-23 — Vavilov tails and moments are discretization-biased

The CDF grid covers only `u = 1/258 ... 257/258`
([`VavilovSampler.cpp:53–69`](../src/VavilovSampler.cpp#L53-L69)). Uniform draws
outside that interval are clamped to the endpoint quantiles
([`VavilovSampler.cpp:21–32`](../src/VavilovSampler.cpp#L21-L32),
[`VavilovSampler.cpp:98–112`](../src/VavilovSampler.cpp#L98-L112)), creating
point masses in both tails. Table “moments” are the unweighted mean of 257
quantile samples rather than an integral over the represented distribution
([`VavilovSampler.cpp:71–95`](../src/VavilovSampler.cpp#L71-L95)). Bilinear
interpolation and the final physical clamp further change them. The README's
claim that the first two moments match catima “exactly” is therefore too strong
([`README.md:146–159`](../README.md#L146-L159)).

Validate sampled mean, variance, quantiles, and tail probabilities against
direct `ROOT::Math::VavilovAccurate` draws over the full κ/β² domain; document
measured interpolation error rather than exactness.

#### F-24 — Small convolution components can vanish

The total κ may be inside the supported band while both decomposed terms
`(1-beta²)κ` and `beta²κ` are at or below `10^-3`. Both are then skipped and the
routine returns exactly zero
([`VavilovSampler.cpp:134–158`](../src/VavilovSampler.cpp#L134-L158)). Around
`κ = 0.001–0.002` and intermediate β² this creates a discontinuous region with
no straggling. Handle each component with its correct limiting distribution or
sample the two-dimensional distribution directly.

#### F-25 — κ and σ use inconsistent charge models

The width comes from catima using configurable
`gStragglingConfig`, specifically intended to provide an effective-charge model
([`EnergyLoss.cpp:120–133`](../src/EnergyLoss.cpp#L120-L133)). The Vavilov shape
parameter independently uses bare nuclear `Z²`
([`EnergyLoss.cpp:255–278`](../src/EnergyLoss.cpp#L255-L278)). For partially
stripped low-energy heavy ions, this combination needs justification and
validation. Derive κ from a charge prescription consistent with the variance,
or document and benchmark why bare Z is correct for the chosen formalism.

### Configuration, I/O, and tooling

#### F-26 — The TOML schema is permissive where physics requires strictness

The generic getters silently ignore missing keys and wrong types
([`ControlFile.cpp:77–109`](../src/ControlFile.cpp#L77-L109)). Unknown keys are
also ignored. There are no systematic finite/range checks for pressure,
temperature, beam energy/spread, stopping scales, histogram dimensions,
resolution values, strip bounds, event/thread counts, max time, method, or
output mode. NaN/Inf and negative values can reach square roots, density
calculations, loops, and ROOT constructors.

Layer documentation requires exactly one thickness representation, but the
loader rejects only the “both” case; “neither” silently retains a compiled
default ([`ControlFile.cpp:126–170`](../src/ControlFile.cpp#L126-L170),
[`README.md:220–246`](../README.md#L220-L246)). Introduce a strict schema phase
that accumulates all diagnostics, rejects unknown/wrong-type/non-finite values,
and distinguishes required fields from explicit optional defaults.

#### F-27 — ROOT output objects are not validated

Fresh/update file opening, `Get("events_MeV")`, `Get("MC")`, `Get("traces")`,
tree schema, branch-address results, and directory creation are used without
zombie/null/error checks
([`Output.cpp:9–54`](../src/Output.cpp#L9-L54),
[`Simulator.cpp:314–329`](../src/Simulator.cpp#L314-L329)). A missing tree or
incompatible prior schema can become a null dereference or corrupt update.
Validate every ROOT operation, compare the full schema before append, write to
a temporary file, and atomically rename only after success.

#### F-28 — External nuclide loaders are unsafe

Clang's analyzer found possible writes beyond both dynamically sized loader
arrays ([`NuclideFinder.cpp:240–329`](../src/NuclideFinder.cpp#L240-L329),
[`NuclideFinder.cpp:372–388`](../src/NuclideFinder.cpp#L372-L388)). The NUBASE
loader also uses a non-standard stack VLA, performs several fixed-column reads
without first validating line length, and writes `halfLifeUncMeas` when parsing
spin provenance, leaving `spinUnc` uninitialized
([`NuclideFinder.cpp:583–635`](../src/NuclideFinder.cpp#L583-L635)); that flag is
later read to decide whether to accept a spin
([`NuclideFinder.cpp:647–660`](../src/NuclideFinder.cpp#L647-L660)).
`GetMass(index)` itself lacks the bounds checks used by adjacent accessors
([`NuclideFinder.cpp:60–80`](../src/NuclideFinder.cpp#L60-L80)).

Replace raw arrays/VLA with vectors, parse one validated record at a time, use
checked indexing, initialize every field, and add malformed/truncated-file
tests.

#### F-29 — Reading excitation energy mutates particle state

For multiple excitation states, `GetEexc()` draws from its private RNG on every
call ([`Particle.cpp:84–101`](../src/Particle.cpp#L84-L101)). `GetKE()` calls it
implicitly ([`Particle.cpp:144–151`](../src/Particle.cpp#L144-L151)), so logging,
branch filling, or repeated calculations can silently select different states
for the same particle. Sample once when the particle/event state is created;
make accessors const and deterministic.

#### F-30 — Runs are not reproducible from their outputs

Single-thread execution constructs `TRandom3(0)`, which auto-seeds
nondeterministically ([`Simulator.cpp:64–66`](../src/Simulator.cpp#L64-L66)). MT
uses hardcoded per-worker seeds, so changing thread count changes event streams
([`Simulator.cpp:400–413`](../src/Simulator.cpp#L400-L413)). No user seed,
resolved configuration, software revision, physics-model version, or seed is
written to the ROOT output. Add a master seed to the schema and metadata; derive
counter-based/per-event streams so results do not depend on scheduling or
thread count.

#### F-31 — A shipped example is rejected

[`17F_alpha_alpha.toml:9–15`](../ControlExamples/17F_alpha/17F_alpha_alpha.toml#L9-L15)
sets `strip` and `strip_first`/`strip_last` together. The parser correctly
rejects that combination
([`ControlFile.cpp:437–459`](../src/ControlFile.cpp#L437-L459)); executing the
example produced that error. Fix the example and add a smoke test that loads
every shipped control.

#### F-32 — README and implemented schema/output disagree

Notable discrepancies:

- “Every section is optional” is false: detector strip selection is mandatory,
  reacted runs need complete particles, and placeholder defaults are not usable
  ([`README.md:47–53`](../README.md#L47-L53),
  [`ControlFile.cpp:437–459`](../src/ControlFile.cpp#L437-L459)).
- `physics.stopping`, `reaction.residue_excitation`, and
  `reaction.angular_distribution` are implemented but absent from the control
  reference ([`ControlFile.cpp:302–370`](../src/ControlFile.cpp#L302-L370),
  [`README.md:180–218`](../README.md#L180-L218)).
- The README names `BeamEnergyAccel`, `Kbi`, `Kbr`, `Kbeam_exit`, `Kl_exit`,
  `Kh_exit`, and light/heavy four-vectors
  ([`README.md:284–312`](../README.md#L284-L312)). Actual branches use names such
  as `beam_energy_accel`, `beam_energy_gas`, `beam_energy_reaction`,
  `beam_energy_exit`, `evap_energy_exit`, and `residue_energy_exit`, and no
  four-vector branches are created
  ([`Output.cpp:62–105`](../src/Output.cpp#L62-L105)).
- The README says κ below `10^-3` falls back in a way that matches the Landau
  limit, but the implementation draws a Gaussian (F-11).

Generate schema/output documentation from a tested definition, or at minimum
make a checked example and branch-list test authoritative.

#### F-33 — `n_events` means events per strip

The run loops through every selected strip and calls `Simulate(..., NEvents, ...)`
for each ([`Simulator.cpp:331–343`](../src/Simulator.cpp#L331-L343)). A range of
11 strips with `n_events = 100000` writes 1.1 million entries, not 100,000. This
may be intentional, but neither the key name nor README says “per strip.” Rename
or document it and print the resolved total before launch.

#### F-34 — `srim-cache` permits shell interpretation of inputs

The tool interpolates the generator path, gas, pressure, temperature, and output
path into a fixed buffer, then calls `std::system()`
([`srim-cache.cpp:133–156`](../src/srim-cache.cpp#L133-L156)). Ion names are
restricted, but gas is unquoted and single quotes in the generator/output path
are not escaped. A crafted control or environment variable can alter the shell
command. Invoke the generator directly with an argument vector (`exec`/spawn),
validate all tokens, and avoid the 1024-byte truncation risk.

#### F-35 — The legacy converter is not dependable in the documented shell

The executable uses `#!/usr/bin/env python3`, but Python is absent from the Nix
dev shell ([`flake.nix:68–81`](../flake.nix#L68-L81)); running its documented
help command there exits 127 with `env: 'python3': No such file or directory`.
Additionally:

- `--keep` is parsed but never used
  ([`legacy_msc_to_toml.py:326–358`](../tools/legacy_msc_to_toml.py#L326-L358));
- an unmatched input glob produces no output and returns success;
- destinations are overwritten without confirmation; and
- TOML strings are quoted without escaping quotes, backslashes, or control
  characters ([`legacy_msc_to_toml.py:143–157`](../tools/legacy_msc_to_toml.py#L143-L157)).

Add Python to the shell/package, validate nonempty inputs, use a TOML writer or
correct escaping, make overwrite explicit, and test conversion fixtures.

#### F-36 — No automated regression suite or CI exists

The Makefile has only build and clean targets
([`Makefile:28–49`](../Makefile#L28-L49)); the flake has no `checks`/test phase
([`flake.nix:36–67`](../flake.nix#L36-L67)). The pre-commit hook formats staged
files but runs no compile, schema, invariant, or physics tests. The dynamically
reproduced defects are compact enough to become deterministic regression tests.

## Low-severity maintainability and style findings

### F-37 — Detector scratch storage assumes two columns

`baseDE[3]` and `noisedDE[3]` are indexed through `AnodeCols + 1`
([`Propagation.cpp:219–240`](../src/Propagation.cpp#L219-L240)). The hardcoded
geometry currently fixes `AnodeCols = 2`
([`Geometry.cpp:3–9`](../src/Geometry.cpp#L3-L9)), so current execution is safe,
but a geometry change would overflow these arrays. Size them from the geometry
or make the two-column invariant compile-time and explicit.

### F-38 — Dormant public routines contain immediate hazards

- `CalculateExcEnergyRange()` writes through `SegEexcRange`, which is initialized
  to null and never allocated
  ([`Simulator.cpp:9–15`](../src/Simulator.cpp#L9-L15),
  [`Kinematics.cpp:401–444`](../src/Kinematics.cpp#L401-L444)).
- `PrintEnergetics()` loops through `AnodeRows + 1` and indexes
  `DeltaEB[AnodeRows]` on the last iteration
  ([`Kinematics.cpp:447–519`](../src/Kinematics.cpp#L447-L519)).
- The angular-grid API divides by `ThSteps - 1` and `PhiSteps - 1` without
  validating that either count exceeds one
  ([`EventLoop.cpp:392–396`](../src/EventLoop.cpp#L392-L396)).

These paths are not used by the primary TOML run today. Remove them or restore
them behind tests before advertising them as public methods.

### F-39 — Dead and misleading interfaces remain

`run.reac_class` is parsed but never consumed
([`ControlFile.cpp:433–435`](../src/ControlFile.cpp#L433-L435));
`WriteTraces()` and `Particle::GetEnergyLoss()` are public no-ops
([`Output.cpp:340–343`](../src/Output.cpp#L340-L343),
[`Particle.cpp:109–113`](../src/Particle.cpp#L109-L113)); `BeamInit` is allocated
and copied but unused (F-18). Removing or deprecating these APIs would make the
supported execution path substantially easier to reason about.

### F-40 — Warning and formatter debt is substantial

An isolated GCC build with the project's flags plus
`-Wall -Wextra -Wpedantic -Wconversion -Wshadow -Wnon-virtual-dtor`
`-Wold-style-cast -Woverloaded-virtual -Wnull-dereference -Wdouble-promotion`
`-Wformat=2` completed with 188 warnings:

| Warning class | Count |
| --- | ---: |
| float conversion | 90 |
| shadowing | 34 |
| old-style cast | 30 |
| double promotion | 17 |
| other conversion | 10 |
| deprecated copy | 4 |
| VLA | 1 |
| unused variable | 1 |
| unused parameter | 1 |

Many conversions are at ROOT visualization/output boundaries and are benign,
but the volume hides meaningful warnings such as the VLA and RNG seed
narrowing. The repository's hook-equivalent clang-format check fails only at
[`EventLoop.cpp:418–419`](../src/EventLoop.cpp#L418-L419). `taplo fmt --check`
fails only for [`basic.toml`](../basic.toml). Establish a warning baseline,
eliminate it by category, and enforce non-mutating checks in CI.

## Physics assumptions and validation gaps

These are not necessarily implementation defects, but they bound what results
can mean and should be explicit in output metadata and scientific validation:

1. **Conditioned reaction placement.** A reaction strip is chosen by the user
   and the vertex is uniform in z inside that strip
   ([`EventLoop.cpp:163–170`](../src/EventLoop.cpp#L163-L170)). There is no
   energy-dependent reaction cross section, beam attenuation, target-density
   weighting, or absolute yield model.
2. **Straight-line transport.** Direction is fixed during each propagation;
   there is no angular straggling, multiple Coulomb scattering, charge exchange,
   or nuclear interaction in gas/windows. Neutrals traverse the gas without
   interaction or deposit.
3. **Hardcoded detector.** Geometry is a fixed 20-row, two-column layout
   ([`Geometry.cpp:3–60`](../src/Geometry.cpp#L3-L60)); configuration cannot
   represent a different chamber or misalignment.
4. **Stopping-table interpolation.** Gas stopping and variance are
   pre-tabulated over total KE `10^-4` to `10^4` MeV and clamped outside that
   range ([`EnergyLoss.cpp:95–133`](../src/EnergyLoss.cpp#L95-L133),
   [`EnergyLoss.cpp:193–207`](../src/EnergyLoss.cpp#L193-L207)). Runs outside the
   table range receive endpoint physics without a warning.
5. **Mixed model for mean and width.** Mean dE/dx may be catima, SRIM, or their
   arithmetic mean, while variance always comes from catima under a distinct
   charge model ([`EnergyLoss.cpp:120–178`](../src/EnergyLoss.cpp#L120-L178)).
   This is a deliberate construction but needs validation for every ion/energy
   regime used.
6. **Excitation continuum.** `forced` and `uniform` draw a continuous excitation
   energy, without discrete levels, widths, spin/parity selection, branching
   ratios, or lifetime/displaced-decay treatment.
7. **Detector response.** Noise is a relative Gaussian width proportional to
   each event's deposit. Correlations, thresholds, saturation, gain variation,
   non-Gaussian electronics response, and digitization are not modeled.

Two choices were explicitly accepted in the repository's prior human review
and are therefore recorded as assumptions, not reopened here as defects:

- Exit-window straggling is intentionally omitted because no downstream
  detector is presently modeled
  ([existing review](REVIEW-window-dedx-scale-srim-table-2026-09-03.md)).
- A layer `dedx_scale` is implemented as equivalent thickness, so it scales
  both mean loss and straggling. This was provisionally accepted against the
  37Cl/87Rb calibration; it should still be revalidated for other materials and
  ions.

## Affected results and rerun guidance

Treat the following existing outputs as suspect:

| Output/run type | Why | Action after fixes |
| --- | --- | --- |
| Any reacted run that omitted `reaction.residue_excitation` | Inherits F-03. | Select a physically justified model and regenerate. |
| 37Cl + alpha “alpha”/elastic samples | F-03 and F-12; likely wrong excitation and angular law. | Reconfigure as ground-state plus validated Rutherford/experimental distribution, then regenerate. |
| Any multi-step chain using `theta_cm[]`/`phi_cm[]` | F-08 corrupts truth angles; F-07 can create inconsistent failure histories. | Regenerate; do not repair angles post hoc. |
| Events in which a particle stops or reaches a boundary/vertex | F-05, F-06, F-09, and F-22. | Regenerate after fractional-step and termination-state fixes. |
| Any `threads > 1` production | F-02 makes physics initialization undefined; F-17 weakens completion guarantees. | Regenerate after race fixes and compare statistically with one-thread reference. |
| Runs based on shipped `eres` comments | F-13 changes noise by units and scaling. | Establish intended units, migrate controls, regenerate. |
| The 19 examples without explicit beam energy | They ran at 100 MeV regardless of likely experiment intent. | Recover intended beam energy and regenerate. |
| Update-mode MT runs | Final merge recreates the destination. | Verify backups/source data; do not use MT update until fixed. |

## Verification matrix

| Check | Result | Notes |
| --- | --- | --- |
| `nix build --no-link` | **Pass** | Full packaged build completed. |
| Isolated strict GCC build | **Pass with 188 warnings** | Warning breakdown is in F-40. |
| Clang static analysis (`clang-analyzer-*`, selected `bugprone-*`) | **Completed with findings** | Confirmed possible null parser input, fixed-size detector scratch coupling, `Chroma` leak, and AME/NUBASE bounds hazards. Third-party header noise was suppressed when interpreting results. |
| Hook-style clang-format check | **Fail** | Only `src/EventLoop.cpp:418–419`. |
| `taplo check` over 30 tracked control TOMLs | **Pass** | Syntax only; semantic invalidity remains in F-31. |
| `taplo fmt --check` | **Fail** | `basic.toml` only. |
| Missing control file | **Fail as expected, then erroneous continuation** | Parse error was printed; process continued and segfaulted in headless GUI setup, exit 129 (F-15). |
| Shipped `17F_alpha_alpha.toml` | **Rejected** | Conflicting strip forms (F-31). |
| Forbidden primary-reaction invariant test | **Fail** | Reproduced impossible accepted event described in F-01. |
| Stopping-energy conservation test | **Fail** | 2 MeV in, 1.995736062 MeV recorded; about 4.264 keV omitted (F-05). |
| Two-step CM-angle truth test | **Fail** | Two sampled angles, identical stored array elements (F-08). |
| Zero-step reacted run under ASan/UBSan | **Fail** | Misaligned uninitialized-pointer access and segmentation violation (F-04). |
| Two-worker/eight-event smoke run | **Pass, limited** | Both trees had 8 entries and `events_MeV` had one friend; this does not exercise away F-02 or prove stress safety. |
| Converter inside `nix develop` | **Fail** | `python3` missing, exit 127 (F-35). |

No unit, integration, golden-output, statistical-distribution, sanitizer, or
cross-platform test suite exists in the repository, so the checks above were
performed manually.

## Prioritized remediation plan

### P0 — Blockers before new production data

1. Make control loading transactional and fatal on failure. Validate a strict
   typed schema, require a complete reaction for reacted runs, check A/Z
   conservation and chain continuity, and reject non-finite/out-of-range values.
2. Replace the chain-wide reaction boolean with explicit per-step state. Fix the
   negative-availability return, define partial-chain semantics, and store each
   sampled CM angle immediately.
3. Redesign propagation around an explicit result containing final state and
   termination reason. Shorten the last step to a vertex/boundary/stop point and
   account for all energy before returning.
4. Eliminate mutable global physics configuration from worker setup. Guarantee
   read-only/thread-local catima operation, validate cache capacity/eviction,
   isolate logs, and fail the whole job on any worker error.
5. Remove the unsafe excitation default. Require an explicit per-channel/per-step
   physical model and correct the elastic examples.
6. Correct and statistically validate the Vavilov implementation, particularly
   the Landau limit, component thresholds, tails, moments, and heavy-ion charge
   treatment.

### P1 — Reliability and scientific traceability

1. Migrate all example controls: explicit beam energies, unambiguous detector
   resolution units, explicit reaction models, and one valid strip selector.
2. Add event-level invariant tests: four-momentum conservation at each reaction,
   A/Z conservation at configuration load, and energy accounting across gas and
   windows when stochastic effects are disabled.
3. Add propagation convergence tests over step sizes and exact tests for every
   termination surface/reason, including unreachable vertices and oblique exit
   windows.
4. Add distribution tests against direct ROOT Vavilov evaluation across κ/β²,
   with tolerances for mean, variance, central quantiles, and tails.
5. Make ROOT writes checked and atomic; define or reject multi-threaded update;
   store resolved controls, revision, model selections, and RNG seed in every
   output.
6. Convert core ownership to RAII and remove unconditional trace buffers from
   non-visual runs.

### P2 — Maintenance controls

1. Add CI jobs for Nix build, strict compile, clang-format, Taplo syntax/format,
   sanitizer tests, example-load tests, and the converter.
2. Remove or repair dormant APIs and unused settings; split large legacy classes
   into configuration, physics, transport, and output components.
3. Pay down warnings by category and make new warnings fail CI.
4. Replace `system()` in `srim-cache` and package/test Python for the converter.

## Suggested acceptance tests

A fix should not be considered complete until at least these checks pass:

- forbidden primary and later-chain steps produce a single coherent history;
- no-step reacted controls fail validation without entering setup;
- with straggling/noise/windows disabled, initial KE equals all gas deposits plus
  surviving exit KE within a stated floating-point tolerance;
- varying `sim_step` across at least an order of magnitude converges at vertices,
  stopping points, and boundaries;
- a vertex beyond the beam range is rejected and labeled distinctly;
- each multi-step truth angle matches the angle used to build that step's
  four-vectors;
- four-momentum is conserved for every allowed step before transport, and A/Z
  are conserved by validated controls;
- sampled Vavilov moments and quantiles meet documented tolerances across the
  supported κ/β² grid and approach Landau/Gaussian limits correctly;
- one-thread and many-thread runs are race-free and statistically consistent;
- failed workers cannot produce a “successful” merged file, and update mode
  cannot overwrite existing data unexpectedly;
- every shipped TOML loads in a dry-run schema check and carries explicit beam,
  noise-unit, and reaction-model choices; and
- CLI success/error statuses follow normal Unix conventions.

## Positive observations

- The Nix package build is reproducible enough to build the complete simulator
  and helper against pinned ROOT/catima inputs.
- Relativistic beam momentum and two-body CM momentum are used, rather than a
  nonrelativistic shortcut, in event kinematics
  ([`Kinematics.cpp:32–67`](../src/Kinematics.cpp#L32-L67),
  [`Kinematics.cpp:175–206`](../src/Kinematics.cpp#L175-L206)).
- Isotropic sampling correctly draws uniform `cos(theta)`, and the Rutherford
  inverse-CDF implementation is clearly documented
  ([`Kinematics.cpp:146–164`](../src/Kinematics.cpp#L146-L164)).
- Gas composition/density handling is explicit, including the P10 distinction
  between mixture molar mass and catima's stoichiometric mass
  ([`Materials.cpp:23–74`](../src/Materials.cpp#L23-L74)).
- The per-event accelerator → degrader → entrance-window → gas chain is easy to
  follow and keeps mean setup values separate from event sampling
  ([`EventLoop.cpp:142–158`](../src/EventLoop.cpp#L142-L158)).
- Gas-step sampled loss is clamped to physical support, detector reset/sentinel
  conventions are documented in code, and the two-tree friend layout worked in
  a fresh-file smoke test.
- The parser already rejects conflicting strip forms, unknown named physics
  models, and an invalid Rutherford cutoff. These checks provide a good pattern
  for the missing schema validation.
- Existing review notes clearly record human decisions around stopping-power
  averaging, exit-window straggling, and layer `dedx_scale`; preserving that
  audit trail is good practice.

## Audit limitations

- No ThreadSanitizer run was completed; the concurrency findings follow from
  direct unsynchronized accesses and dependency cache design.
- GUI behaviour was not validated with a display server. The missing-control
  continuation was tested headlessly; the continuation bug is independent of
  the eventual GUI crash.
- SRIM table generation under Wine was not executed, and no external SRIM table
  corpus was available for numerical comparison.
- Dynamic tests were deliberately small and deterministic. They prove the named
  invariants fail but do not quantify bias across production ensembles.
- Optional legacy/angular-scan/nuclide-loader paths received static inspection,
  not complete end-to-end execution.

