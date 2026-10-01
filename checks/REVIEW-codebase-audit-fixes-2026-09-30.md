# Review checklist — codebase-audit fixes, 2026-09-30
<!---->
This checklist covers the uncommitted remediation of every finding in
`checks/CODEBASE-AUDIT-2026-09-30.md`. Line numbers refer to the formatted
working tree on 2026-09-30. No checkbox is pre-approved: checks indicate human
review by the maintainer.
<!---->
The dedicated physics checklist is
`checks/REVIEW-core-physics-audit-fixes-2026-09-30.md`.
<!---->
## Blocking compatibility and release decisions
<!---->
- [ ] **Accept the control-schema break.** `src/ControlFile.cpp:162-180,
      230-308, 365-489` now requires explicit beam energy, event count, output,
      seed, and per-step physics choices; renames `n_events` to
      `events_per_strip`; and replaces ambiguous `eres` with
      `noise_sigma_mev` or `resolution_fwhm_percent`. Confirm downstream control
      generators have been migrated. The new schema is documented at
      `README.md:61-210`.
- [ ] **Accept the ROOT-schema break.** `src/Output.cpp:18-263` uses the
      `events_MeV` plus `MC` layout, schema version 3, explicit requested/actual
      reaction strips, stable event indices, stop positions, and termination
      codes. Confirm existing analysis code is updated before production runs.
- [ ] **Accept the C++ source-API break.** The supported interface at
      `include/Simulator.hpp:47-162` removes unsafe angular-grid, trace-database,
      excitation-range, and legacy-energetics entry points; `include/Particle.hpp:17-72`
      removes the no-op particle energy-loss method and makes state reads const.
      The compatibility impact is called out at `README.md:360-367`.
- [ ] **Approve the rerun boundary.** Results made with the audited revision may
      contain forbidden reactions, lost terminal/boundary energy, repeated CM
      angles, forced excitation, or mis-scaled detector noise. Decide which
      existing samples must be regenerated after physics review.
<!---->
## Critical findings
<!---->
- [ ] **F-01 — forbidden primary reactions.** Review the invariant-mass
      threshold and false return at `src/Kinematics.cpp:112-119`, full rollback
      at `src/Kinematics.cpp:242-258`, and the event-level distinction between
      attempted and realized reactions at `src/EventLoop.cpp:164-213` and
      `src/Propagation.cpp:251-259`. Confirm the regression assertions at
      `tests/root_checks.cpp:113-169` leave no vertex, product truth, or
      `reaction_strip` after rejection.
- [ ] **F-02 — shared physics races.** Confirm physics settings are
      instance-owned at `include/EnergyLoss.hpp:17-32`, all catima cache access
      is serialized at `src/EnergyLoss.cpp:118-135, 356-371` and
      `src/Materials.cpp:187-226`, prewarming remains an optimization only at
      `src/Simulator.cpp:244-274`, and worker logs are isolated at
      `src/Simulator.cpp:228-231, 458-461`.
- [ ] **F-03 — unsafe excitation default.** Confirm each reaction step must
      choose its excitation and angular model at `src/ControlFile.cpp:420-455,
      556-570`; the final step cannot be `forced`; and the implementation at
      `src/Kinematics.cpp:121-145` gives forced intermediate residues only the
      excitation range needed by the remaining chain. Review every migrated
      `ControlExamples/**/*.toml` choice as a physics decision, not merely a
      schema migration.
<!---->
## High-severity findings
<!---->
- [ ] **F-04 — zero-step reacted run.** Confirm reacted controls without at
      least one complete step are rejected at `src/ControlFile.cpp:365-459,
      556-570`, particle-pointer arrays are null-initialized at
      `src/Simulator.cpp:32-53`, and `tests/run-tests.sh:24-36` exercises the
      rejection under normal and sanitizer builds.
- [ ] **F-05 — discarded stopping energy.** Review the stopped-step distance,
      deposit, and final-energy accounting at `src/EnergyLoss.cpp:251-324` and
      its application before termination at `src/Propagation.cpp:198-244`.
      Confirm the event-level conservation assertion at
      `tests/root_checks.cpp:85-110` uses gas-entry KE versus all gas deposits.
- [ ] **F-06 — unreachable reaction vertices.** Confirm the explicit transport
      result at `include/Simulator.hpp:114-141` and that kinematics are invoked
      only after `ReachedVertex` at `src/EventLoop.cpp:168-213`; stopped,
      escaped, timed-out, and invalid beams must remain unreacted.
- [ ] **F-07 — partial decay-chain rollback.** Review transactional product
      construction and truth clearing at `src/Kinematics.cpp:76-84, 242-258`.
      Confirm `tests/controls/later-step-forbidden.toml:26-47` reaches a valid
      first step followed by a forbidden second step, and
      `tests/root_checks.cpp:113-169` rejects every partial truth field.
- [ ] **F-08 — overwritten multi-step angles.** Confirm each step stores its CM
      angles immediately at `src/Kinematics.cpp:149-172` and the two-step test
      at `tests/root_checks.cpp:172-193` rejects repeated final-step values.
- [ ] **F-09 — discarded boundary fractions.** Review first-intersection
      selection at `src/Propagation.cpp:117-196`, fractional transport and
      deposit at `src/Propagation.cpp:198-240`, and coarse/fine boundary checks
      at `tests/root_checks.cpp:196-218` and `tests/run-tests.sh:82-89`.
- [ ] **F-10 — unsafe nonpositive step.** Confirm strict finite-positive
      validation at `src/ControlFile.cpp:515-530` and the defensive transport
      guard at `src/Propagation.cpp:28-38`. Verify the invalid fixture is
      rejected by `tests/run-tests.sh:24-36`.
- [ ] **F-11 — wrong low-kappa limit.** Confirm the sampler retains the
      lowest-supported Vavilov shape below kappa `1e-3` and uses a Gaussian only
      above kappa 10 at `src/VavilovSampler.cpp:107-117`; review the quantile and
      limiting checks at `tests/vavilov_statistics.cpp:11-70`.
- [ ] **F-12 — non-elastic elastic examples.** Confirm the parser enforces
      ground-state plus Rutherford settings for elastic identities at
      `src/ControlFile.cpp:598-624`, angular models are per step at
      `src/Kinematics.cpp:149-170`, and the representative 37Cl elastic control
      is explicit at `ControlExamples/37Cl_alpha_n/37Cl_alpha_alpha.toml:73-83`.
- [ ] **F-13 — ambiguous detector-resolution units.** Review the mutually
      exclusive schema at `src/ControlFile.cpp:230-308`, the absolute-sigma and
      relative-percent calculations at `src/Propagation.cpp:266-353`, and all
      migrated examples. Confirm the values—not only their key names—match the
      intended detector model.
- [ ] **F-14 — implicit 100 MeV beam.** Confirm `beam.energy` is required at
      `src/ControlFile.cpp:168-174` and finite-positive at
      `src/ControlFile.cpp:515-527`; audit the explicit energies in every
      `ControlExamples/**/*.toml`, which `tests/run-tests.sh:29-36` loads.
- [ ] **F-15 — CLI/status inversion.** Confirm argument, parse, check-only,
      simulation, and success exits at `src/main.cpp:39-69` follow Unix status
      conventions; the negative CLI cases are exercised at
      `tests/run-tests.sh:24-28`.
- [ ] **F-16 — unvalidated reactions.** Review complete per-step parsing at
      `src/ControlFile.cpp:365-459`, required model checks at
      `src/ControlFile.cpp:556-570`, and nuclide plus compound/per-step A/Z
      validation at `src/ControlFile.cpp:583-624`.
- [ ] **F-17 — unsafe MT completion/update.** Confirm multi-threaded update is
      rejected at `src/ControlFile.cpp:547-554`; worker failures and entry-count
      mismatches abort at `src/Simulator.cpp:646-742`; and merge output is
      validated and atomically renamed at `src/Simulator.cpp:744-803`.
- [ ] **F-18 — leaks and unconditional trace allocation.** Review simulator
      ownership cleanup at `src/Simulator.cpp:100-188`, lazy particle trace
      allocation at `src/Particle.cpp:26-42`, RAII energy-loss ownership at
      `include/Particle.hpp:98-108`, visualization-only trace creation at
      `src/Simulator.cpp:582-587`, and the restored memory guard at
      `src/Simulator.cpp:190-206`.
<!---->
## Medium-severity findings
<!---->
- [ ] **F-19 — exothermic equality threshold.** Confirm the single
      invariant-mass tolerance at `src/Kinematics.cpp:112-145` permits a
      zero-available-energy state, while recursive forced-chain thresholds are
      computed consistently at `src/Simulator.cpp:392-417`.
- [ ] **F-20 — passive layers adding energy.** Confirm window/degrader
      straggling output is clamped to `[0, Ein]` at
      `src/Materials.cpp:201-226`.
- [ ] **F-21 — oblique exit foil.** Confirm downstream charged particles use
      thickness divided by `abs(cos(theta))` at
      `src/Propagation.cpp:372-400`, while side/upstream exits do not traverse
      the downstream foil.
- [ ] **F-22 — conflated stop/escape/timeout.** Review termination codes at
      `include/Simulator.hpp:114-131`, propagation assignments at
      `src/Propagation.cpp:66-115, 146-195, 232-248`, and persisted branches at
      `src/Output.cpp:59-97, 202-252`.
- [ ] **F-23 — Vavilov tail/moment bias.** Review midpoint probability nodes and
      tail extrapolation at `src/VavilovSampler.cpp:41-104`, analytic
      standardization at `src/VavilovSampler.cpp:56-69`, documented limitations
      at `README.md:250-258`, and numeric tolerances at
      `tests/vavilov_statistics.cpp:11-70`.
- [ ] **F-24 — vanished small convolution terms.** Confirm the old split
      convolution is gone and one direct three-dimensional quantile table is
      constructed/sampled at `include/VavilovSampler.hpp:4-43` and
      `src/VavilovSampler.cpp:41-117`.
- [ ] **F-25 — inconsistent charge models.** Check the mass-fraction sum of
      `z_eff^2 * Z/A` for every gas component at
      `src/EnergyLoss.cpp:267-307`; it must use the same configured charge model
      as the catima variance and must not multiply separate mixture averages.
- [ ] **F-26 — permissive TOML.** Review unknown-key/type accumulation at
      `src/ControlFile.cpp:40-160`, exact layer/noise schemas at
      `src/ControlFile.cpp:181-308`, and all finite/range/enum constraints at
      `src/ControlFile.cpp:464-581`. Confirm assignment is transactional only
      after an empty error set at `src/ControlFile.cpp:627-647`.
- [ ] **F-27 — unchecked ROOT I/O.** Review schema and branch validation at
      `src/Output.cpp:24-193`, checked fills/writes and staged-file validation at
      `src/Simulator.cpp:465-639`, and update/duplicate-cycle assertions at
      `tests/root_checks.cpp:55-82, 270-301`.
- [ ] **F-28 — unsafe AME/NUBASE loaders.** Confirm checked fixed-width parsing,
      temporary vectors, and transactional replacement at
      `src/NuclideFinder.cpp:336-535`; review bounds-checked accessors at
      `src/NuclideFinder.cpp:129-255` and malformed-file tests at
      `tests/nuclide_loader.cpp:33-75`.
- [ ] **F-29 — excitation resampled on read.** Confirm `GetEexc()` and `GetKE()`
      are deterministic const reads at `src/Particle.cpp:93-105, 147-155`, and
      random selection is an explicit state transition at
      `src/Particle.cpp:295-342`.
- [ ] **F-30 — irreproducible RNG/provenance.** Review required seed parsing at
      `src/ControlFile.cpp:479-490`, per-event SplitMix64 derivation at
      `src/Simulator.cpp:421-435`, event reseeding at
      `src/EventLoop.cpp:108-110`, metadata at `src/Output.cpp:305-372`, and the
      one-thread/many-thread exact comparison at
      `tests/root_checks.cpp:220-265`.
- [ ] **F-31 — rejected shipped example.** Confirm
      `ControlExamples/17F_alpha/17F_alpha_alpha.toml:21-29` uses only one strip
      selector and that every shipped TOML is checked by
      `tests/run-tests.sh:29-36`.
- [ ] **F-32 — stale README/schema/output documentation.** Review the strict
      schema and physics switches at `README.md:61-258`, output branch and
      sentinel contract at `README.md:295-342`, and model limitations at
      `README.md:343-358` against `src/ControlFile.cpp` and `src/Output.cpp`.
- [ ] **F-33 — ambiguous event count.** Confirm the key is now
      `events_per_strip` at `src/ControlFile.cpp:464-470`, its multiplicative
      meaning is explicit at `README.md:160-165`, and the strip loop at
      `src/Simulator.cpp:588-592` matches that contract.
- [ ] **F-34 — shell injection in `srim-cache`.** Confirm generator arguments
      remain distinct strings and are launched directly with `posix_spawnp` at
      `src/srim-cache.cpp:79-121`; no configuration or environment value is
      evaluated by a shell.
- [ ] **F-35 — unsafe/unavailable converter.** Review dependency-free C++
      conversion, TOML escaping, strict value parsing, collision refusal,
      atomic writes, and explicit source removal at
      `src/legacy-msc-to-toml.cpp:49-130, 486-608`; confirm its fixtures at
      `tests/run-tests.sh:38-54`.
- [ ] **F-36 — missing tests/CI.** Review build, strict-warning, sanitizer, and
      formatter targets at `Makefile:67-166`, the eight-part suite at
      `tests/run-tests.sh:24-102`, and CI execution at
      `.github/workflows/ci.yml:1-25`.
<!---->
## Low-severity findings
<!---->
- [ ] **F-37 — fixed-length detector scratch arrays.** Confirm response scratch
      storage is sized from `AnodeCols` at `src/Propagation.cpp:266-299`, while
      the fixed geometry is declared at `src/Geometry.cpp:3-28`.
- [ ] **F-38 — hazardous dormant APIs.** Confirm the public surface at
      `include/Simulator.hpp:47-162` no longer exposes excitation-range,
      legacy-energetics, or angular-grid routines with null, bounds, or
      divide-by-zero paths; accept the source-compatibility break.
- [ ] **F-39 — dead/no-op interfaces.** Confirm `run.reac_class`,
      `WriteTraces`, `Particle::GetEnergyLoss`, `BeamInit`, and the obsolete
      trace-database path are absent from the active API/schema; see
      `README.md:360-367` for the migration notice.
- [ ] **F-40 — warning/format debt.** Confirm strict warnings are errors at
      `Makefile:48-50, 89-94`, formatter coverage is enforced at
      `Makefile:45-47, 139-155`, and CI invokes both at
      `.github/workflows/ci.yml:18-25`.
<!---->
## Independent verification
<!---->
- [ ] Run `make -j2 check` and inspect all eight reported test groups.
- [ ] Run `make strict`; require zero diagnostics under the configured warning
      set.
- [ ] Run `make -j4 sanitize-check`; require no AddressSanitizer or
      UndefinedBehaviorSanitizer report.
- [ ] Run `make format-check` and `git diff --check`.
- [ ] From a clean build tree, run `nix build path:. --no-link --print-build-logs`.
- [ ] Confirm `sha256sum flake.nix` is
      `d61be43e0a8814961e20ccffd443931623d01f27b3f24ef1b72a5e2b425d7a6a`.
- [ ] Inspect `git status --short`; do not stage, commit, or push until this
      review and the physics checklist are discussed.
<!---->
