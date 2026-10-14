# Core-physics review checklist — codebase-audit fixes, 2026-09-30
<!---->
Per the AI-disclosure section of the README, all simulation-physics changes
must be human-reviewed. Line numbers are exact for the formatted, uncommitted
working tree on 2026-09-30. No code is reproduced here; open each referenced
range and review it. All boxes are intentionally unchecked.
<!---->
The originating findings are F-01 through F-14 and F-19 through F-25 in
`checks/CODEBASE-AUDIT-2026-09-30.md`. Configuration, output, reproducibility,
and tests are included where they affect scientific interpretation.
<!---->
## Physics-policy decisions requiring explicit acceptance
<!---->
- [x] **Accept the schema/API break before using these changes.** Reacted steps
      now require explicit excitation and angular models
      (`src/ControlFile.cpp:365-455, 556-570`), run controls require a seed and
      `events_per_strip` (`src/ControlFile.cpp:464-490`), ambiguous `eres` is
      rejected (`src/ControlFile.cpp:230-308`), and truth output is schema 3
      (`src/Output.cpp:195-255`, `src/Output.cpp:329-363`). Update control
      generators and analysis consumers first.
      > 2026-10-14 ewt: accepted. MUSIC ROOT readers unaffected. Follow-up:
      > migrate the 53 controls in `~/MUSIC/analysis/*/sim_control`
      > (`[detector.eres]` → `resolution_fwhm_percent`, not `noise_sigma_mev`).
- [x] **Accept the low-kappa approximation.** For kappa below ROOT's accurate
      range, the implementation freezes the shape at kappa `1e-3` rather than
      attempting an asymptotic Landau parameter conversion
      (`src/VavilovSampler.cpp:74-117`). Confirm this is adequate for the
      experiment's step sizes and ions.
      > 2026-10-14 ewt: accepted. Correction found at item EnergyLoss.cpp:251-324:
      > with the Eloss >= 0 clamp, the frozen shape biases the MEAN loss of light
      > particles in 10 um steps (~+50% for 5 MeV p, ~+5% for 20 MeV alpha; heavy
      > ions unaffected). Negligible in strip sums; Landau follow-up logged.
- [x] **Accept the excitation models.** `ground`, `uniform`, and `forced` remain
      continuum choices rather than a nuclear level/width/branching model
      (`src/Kinematics.cpp:121-145`; limitations at `README.md:343-358`). Review
      every shipped channel choice independently.
      > 2026-10-14 ewt: accepted, including the ControlExamples choices. Old
      > `forced` was U(2/3·Eavail, Eavail) on every step: all 87Rb sim results
      > and the 37Cl (α,α) (`uniform`) are biased and need rerunning. MUSIC
      > migration: (α,α) → ground/rutherford until TALYS branches land;
      > (α,n)/(α,p) → uniform.
- [x] **Accept layer-scale semantics.** A layer `dedx_scale` is implemented as
      equivalent thickness, so it scales both mean loss and straggling
      (`src/Materials.cpp:138-179`). Confirm this remains consistent with the
      calibrations that introduced those factors.
      > 2026-10-14 ewt: accepted; unchanged since the 2026-09-03 sign-off.
<!---->
## Reaction graph and kinematics
<!---->
- [x] `src/ControlFile.cpp:365-459` — every reaction entry is a complete table
      with named ejectile/residue, explicit per-step excitation and angular
      policies, bounded Rutherford cutoff, and positive stopping scales.
      > 2026-10-14 ewt: accepted. Style follow-up logged below.
- [x] `src/ControlFile.cpp:583-624` — every named nuclide exists; beam plus target
      matches the compound; every parent decomposes into products with conserved
      A and Z; identity elastic scattering is constrained to ground-state plus
      Rutherford sampling.
      > 2026-10-14 ewt: accepted; fixed: a swapped elastic labelling (evap = beam,
      > res = target) is now rejected, since Rutherford theta_cm is the res
      > direction. Uncommitted, untested.
- [x] `src/Simulator.cpp:392-417` — recursive `minEx` is the minimum excitation
      needed to feed all subsequently forced steps. Check signs and indexing for
      exothermic and endothermic chains.
      > 2026-10-14 ewt: accepted. (Old code had the minEx sign inverted and never
      > used it for forced sampling.)
- [x] `src/Kinematics.cpp:3-12, 34-69` — beam momentum, total four-momentum, CM
      boost, and compound excitation use relativistic energy-momentum relations
      and the intended atomic-mass convention.
      > 2026-10-14 ewt: accepted (pre-existing, LISE++-validated).
- [x] `src/Kinematics.cpp:95-145` — invariant-mass availability is authoritative,
      the `1e-9 MeV` tolerance accepts threshold equality, `ground` gives zero
      excitation, `uniform` samples `[0,Eavail]`, and `forced` is limited to
      intermediate residues with enough energy for the remaining chain.
      > 2026-10-14 ewt: accepted.
- [x] `src/Kinematics.cpp:149-172` — isotropic sampling is uniform in
      `cos(theta)` and azimuth; Rutherford inverse-CDF sampling and its minimum
      angle are correct; each step's sampled truth angle is stored immediately.
      > 2026-10-14 ewt: accepted; fixed: 'rutherford' is now rejected after step 0
      > (src/ControlFile.cpp, required-key loop). Uncommitted, untested.
- [x] `src/Kinematics.cpp:182-233` — the two-body momentum radicand uses
      `(mh + Ex)` consistently, products are opposite in the CM frame, the lab
      boost sign is correct, and each heavy product becomes the next decay
      parent's four-vector.
      > 2026-10-14 ewt: accepted.
- [x] `src/Kinematics.cpp:76-84, 242-274` — reaction construction is atomic:
      products remain disabled until the full chain succeeds, every product and
      truth field is cleared on failure, and only the final residue is eligible
      for transport.
      > 2026-10-14 ewt: accepted.
- [x] `src/EventLoop.cpp:164-213` — `beam_energy_reaction` and the vertex are
      recorded only for a realized reaction; an attempted but forbidden chain
      keeps the N/A sentinels.
      > 2026-10-14 ewt: accepted.
- [x] `tests/root_checks.cpp:113-169` and
      `tests/controls/later-step-forbidden.toml:26-47` — independently verify
      that both primary and later-step failures retain no product energies,
      angles, exit states, or residue identity while the beam continues.
      > 2026-10-14 ewt: accepted; adding an assertion that step 0 succeeded and
      > step 1 failed (in progress).
- [x] `tests/root_checks.cpp:172-193` — strengthen or accept the current
      multi-step-angle test. It proves independent stored values, but a deeper
      validation would reconstruct the four-vectors and verify that each stored
      angle is the one actually used for that step.
      > 2026-10-14 ewt: strengthen: adding a four-vector reconstruction check of
      > each stored theta_cm/phi_cm (in progress).
<!---->
## Gas transport and detector boundaries
<!---->
- [x] `include/Simulator.hpp:114-141` — the termination enum and result object
      distinguish not-propagated, reached-vertex, stopped, downstream exit,
      upstream exit, side exit, timeout, and invalid state without overloading
      an energy sentinel.
      > 2026-10-14 ewt: accepted.
- [x] `src/Propagation.cpp:28-95` — null/disabled particles and invalid time/step
      inputs terminate safely; initial energy, rest energy, direction, and
      outside-volume classification are physically consistent.
      > 2026-10-14 ewt: accepted.
- [x] `src/Propagation.cpp:97-196` — each iteration selects the first of the time,
      requested-vertex, detector-face, row, and column boundaries. Check all
      signs for forward, backward, and transverse tracks, including a start
      exactly on a boundary.
      > 2026-10-14 ewt: accepted (near-transverse row-boundary InvalidState edge,
      > ~1e-7 of light particles, not logged).
- [x] `src/Propagation.cpp:198-240` — the chosen fractional distance is passed to
      stopping/straggling, its complete deposit is assigned to the traversed
      cell, position/time advance by the same distance, and a stopping step is
      accounted before returning.
      > 2026-10-14 ewt: accepted.
- [x] `src/EnergyLoss.cpp:251-324` — a sampled loss is never negative or above
      the incoming KE; a stopping step deposits all remaining KE and shortens
      its path proportionally. Review the linear-in-step approximation used to
      locate that stop.
      > 2026-10-14 ewt: accepted; Landau follow-up logged (see low-kappa note).
- [x] `src/EventLoop.cpp:168-213, 277-312` — a reaction is attempted only after
      exact vertex arrival; a rejected reaction resumes the same beam from that
      state without resetting already accumulated deposits.
      > 2026-10-14 ewt: accepted.
- [x] `src/Propagation.cpp:372-438` — only a downstream exit crosses the exit
      foil, oblique tracks use `thickness / abs(cos(theta))`, neutrals retain
      their KE, and the final residue is selected from reaction state rather
      than inferred from transport success.
      > 2026-10-14 ewt: accepted.
- [x] `src/Propagation.cpp:440-453` and `src/Output.cpp:266-303` — termination,
      stop-position, residue-step, and energy sentinels are internally
      consistent for reacted, rejected, stopped, escaped, timed-out, and invalid
      events.
      > 2026-10-14 ewt: accepted.
- [x] `tests/root_checks.cpp:85-110` — with stochastic effects and windows off,
      gas-entry KE equals cathode plus upstream/downstream dead-layer deposits
      within `2e-5 MeV` for a stopping event.
      > 2026-10-14 ewt: accepted.
- [x] `tests/root_checks.cpp:196-218` — coarse and fine step runs hit the exact
      downstream surface and converge in exit energy. Decide whether additional
      upstream, side, timeout, stopped-before-vertex, and oblique-exit fixtures
      are required before production.
      > 2026-10-14 ewt: accepted; adding stopped-before-vertex and product side-exit
      > fixtures (in progress). README sentinel documentation requested.
      > Done in README "Sentinels"; angle sentinel changed -1 → -1000
      > (Simulator::kAngleUnset; -1 is a legal azimuth); output schema 3 → 4.
      > beam_energy_accel now records the per-event accelerator KE (spread
      > included). Fixtures stopped-before-vertex and side-exit added.
<!---->
## Materials, windows, and detector response
<!---->
- [x] `src/Materials.cpp:3-75` — pressure/temperature units and ideal-gas density
      are correct; P10 uses its mole-fraction molar mass for density but the
      stoichiometric mass for catima component weight fractions.
      > 2026-10-14 ewt: accepted.
- [x] `src/Materials.cpp:77-136` — elemental/compound definitions, densities,
      Havar mass fractions, and conversions between mg/cm2 and micrometres are
      correct for the physical foils in use.
      > 2026-10-14 ewt: accepted.
- [x] `src/Materials.cpp:138-183` — disabled layers, positive scale factors, and
      equivalent-thickness handling agree for entrance, exit, and degrader
      layers.
      > 2026-10-14 ewt: accepted.
- [x] `src/Materials.cpp:185-226` — catima receives kinetic energy per nucleon,
      returns energy/sigma per nucleon, and sampled passive-layer output is
      bounded to `[0, Ein]`.
      > 2026-10-14 ewt: accepted.
- [x] `src/EventLoop.cpp:143-162` — the event energy chain is accelerator spread,
      then degrader, then entrance window, then gas entry; the same event RNG is
      used for each stochastic term.
      > 2026-10-14 ewt: accepted.
- [x] `src/Propagation.cpp:266-353` — detector noise uses exactly one selected
      convention: absolute Gaussian sigma in MeV or relative percent FWHM. Check
      that anode channels are independent, the cathode receives one independent
      fluctuation after physical summing, and dead layers remain noiseless.
      > 2026-10-14 ewt: accepted. No pedestal noise on empty channels:
      > it would be cut in the real analysis anyway.
- [x] `ControlExamples/37Cl_alpha_n/37Cl_alpha_alpha.toml:6-10, 29-83` — verify
      the explicit beam energy, absolute noise table, and elastic reaction
      choices against experiment-specific intent; repeat this review for every
      migrated control under `ControlExamples/`.
      > 2026-10-14 ewt: accepted. All 19 restored beam energies match the legacy
      > `.msc` Kb (window and degrader off, so accelerator = gas entry); before
      > the audit they silently ran at the 100 MeV default. Noise follows each
      > file's own comment: 24 `[MeV]` files -> `noise_sigma_mev` (upstream
      > absolute sigma), 16C x2 and basic -> `resolution_fwhm_percent`.
      > Correction to the item-1 discussion: the examples did have noise; from
      > 0908577 (2026-06-10) to the audit the MeV-commented ones ran with
      > near-zero noise because their values were read as percent FWHM.
<!---->
## Stopping tables, effective charge, and Vavilov sampling
<!---->
- [x] `include/EnergyLoss.hpp:17-32` — mean and straggling catima configurations
      are immutable per simulator instance, and the process-global catima cache
      is protected by one mutex.
      > 2026-10-14 ewt: accepted. The only unlocked catima call on the hot path
      > is `z_effective`, which is pure (no cache or static state).
      > Pre-existing bug found here and fixed: catima is built without GLOBAL,
      > so the default mean model `atima14` returned z_eff = 0. That zeroed the
      > Bethe term and gave mean stopping that was too low above 10 MeV/u (-11%
      > at 12, -20% for 16C at 13.75, -32% at 16, -53% at 20 MeV/u), with no
      > change below 10 MeV/u. It is also why atima14 gave sigma_E = 0 (8dc0080).
      > The default is now `winger`, which is what atima14 computes below
      > 30 MeV/u. `atima14` and `global` are rejected for both keys.
      > Building with GLOBAL was considered and rejected: no numerical change in
      > the MUSIC energy range, it needs unlicensed F77 source plus our own glue,
      > and `z_eff_global` passes MeV/u where it expects beta below 30 MeV/u
      > anyway. Converter and README updated. Tests: fast-alpha default ==
      > explicit winger (fails on the old binary); invalid-atima14 and
      > invalid-global-straggling.
- [x] `src/EnergyLoss.cpp:93-193` — mean stopping and variance tables use the
      intended catima configurations, units, reference thickness, interpolation
      grid, SRIM replacement/mean policy, and user-scale ordering.
      > 2026-10-14 ewt: accepted, with the table raised from 256 to 1024 points.
      > The 1e-5 cm thin-layer rate matches catima dE/dx to 2e-9. Below 1 keV/u
      > catima reports the particle stopped, and energy is still conserved. SRIM
      > units match the real tables (0.5 keV-3.7 GeV). Worst interpolation error
      > went from 0.9% to 0.1% at the end of a track and is under 0.03% above
      > 0.25 MeV/u. Suite runtime unchanged.
- [x] `src/EnergyLoss.cpp:225-243` — gas `Z/A` is a mass-fraction sum and remains
      valid for all supported elemental and compound gases.
      > 2026-10-14 ewt: accepted. It reduces to sum(stn*Z)/sum(stn*A) for all
      > seven gases (CF4 0.477, P10 0.458). Only used as the xi fallback.
- [x] `src/EnergyLoss.cpp:259-308` — catima variance scales as path length, while
      gas `dedx_scale` changes the mean only. Confirm that separation is the
      intended calibrated model.
      > 2026-10-14 ewt: accepted. Beam and product scales correct the stopping
      > model, since pressure is measured. Window scales are equivalent thickness
      > (calibrated foil loss). All MUSIC gas scales are 1.0. README now states
      > that the gas scale only affects the mean.
- [x] `src/EnergyLoss.cpp:267-307` — verify dimensions and constants in xi,
      maximum transferable electron energy, and kappa. For mixtures, confirm the
      additive `weight_fraction * z_eff^2 * Z/A` calculation and that
      `physics_.straggling` is the same effective-charge prescription used by
      catima's variance.
      > 2026-10-14 ewt: accepted. K, m_e and u match PDG/CODATA 2018. The mixture
      > sum uses per-component z_eff with the same config that built catima's
      > variance. 87Rb at 4.6 MeV/u, 1 mm, 400 Torr He gives kappa of about 4.3.
- [x] `src/VavilovSampler.cpp:41-72` — the `(log kappa, beta2, probability)` grid
      covers the declared domain and standardizes every ROOT quantile with
      ROOT's analytic mean and variance.
      > 2026-10-14 ewt: accepted. Table is 64x17x513, about 0.5 s per process.
- [x] `src/VavilovSampler.cpp:13-32, 74-105` — interpolation clamps kappa/beta2,
      extrapolates only the probability tails to avoid endpoint point masses,
      and indexes all eight trilinear corners correctly.
      > 2026-10-14 ewt: accepted.
- [x] `src/VavilovSampler.cpp:107-117` — low-kappa requests use the kappa `1e-3`
      Vavilov shape and high-kappa requests use a standard Gaussian; no split
      convolution component can silently vanish.
      > 2026-10-14 ewt: accepted. Note: the switch to a Gaussian at kappa = 10 drops
      > skewness from about 0.16 to 0; mean and sigma stay continuous. Negligible
      > in the thick-step regime. Low-kappa freeze: see item 2 and the follow-up.
- [x] `tests/vavilov_statistics.cpp:11-70` — review the kappa/beta2/probability
      grid and tolerances against direct `ROOT::Math::VavilovAccurate`
      quantiles, plus Gaussian-limit moments and deterministic low-limit
      behavior. Decide whether production acceptance also needs larger sampled
      moment/tail tests for the experiment's exact ion/pressure/step regime.
      > 2026-10-14 ewt: accepted with fixes. Added off-node beta2 (0.01, 0.27) and
      > kappa (0.0137, 0.6, 2.2); before, beta2 interpolation was never
      > exercised. Tightened tolerances from 0.12/0.35 to 0.005/0.06 sigma
      > (measured worst error 0.0007/0.024). Added sampled moments at kappa
      > 0.0137/0.137/2.2, beta2 = 0.01 (|mean| < 0.01, |var-1| < 0.05; worst of 20
      > seeds 0.0043/0.025). It catches swapped beta2/kappa corners, which the
      > old test passed. Production regime: the sampler depends only on
      > (kappa, beta2), and heavy ions are at 0.03-1, so no per-ion test is needed.
- [x] `README.md:244-258` — documentation accurately states boundary accounting,
      low/high-kappa behavior, analytic standardization, physical truncation,
      and the fact that sampled moments are close to—not exactly—the catima
      moments.
      > 2026-10-14 ewt: accepted with a rewrite. "Close to catima" is now broken out
      > by kappa: within 1% at kappa >= 0.01, 7% low at 1e-3, and below that the
      > frozen shape plus clamp biases mean loss (+50% 5 MeV p, +5% 20 MeV alpha).
<!---->
## Determinism, concurrency, and scientific provenance
<!---->
- [x] `src/ControlFile.cpp:310-360, 633-647` — physics selections are parsed into
      a temporary instance-owned configuration and installed only after the
      whole control validates; no worker mutates catima defaults.
      > 2026-10-14 ewt: accepted. Fixed a provenance bug found here: the TOML was
      > re-read from disk by every worker, and again by the metadata writer AFTER
      > the run, so an edit mid-run changed the recorded text (and possibly a
      > late worker's physics). Now `loadCtrlFile` reads it once into
      > `ctrlFileText_`; workers parse that text (`loadCtrlText`), and both the
      > metadata and the update check use it. Test: a threaded run from a
      > read-once pipe must match repro-multi, with stored text identical to the
      > fixture. Mutants that re-read in workers or in metadata both fail.
      > Side effect: a control can now be read from a pipe (toml++ parse_file
      > could not).
- [x] `src/EnergyLoss.cpp:118-135, 356-371`, `src/Materials.cpp:187-226`, and
      `src/Simulator.cpp:244-274` — every catima calculation/cache access is
      serialized, including cache misses after prewarming.
      > 2026-10-14 ewt: accepted. All 9 catima call sites in src/ hold
      > CatimaMutex except the pure `z_effective` (see item 33).
- [x] `src/Simulator.cpp:421-435` and `src/EventLoop.cpp:108-110` — each event RNG
      stream is a deterministic function of master seed, requested strip, and
      per-strip event index; TRandom3 seed zero is avoided.
      > 2026-10-14 ewt: accepted. RNG is per instance (no gRandom). Limit: TRandom3
      > takes a 32-bit seed, so about 1% of 10k-event strips have one duplicated
      > event pair. Negligible.
- [x] `src/Simulator.cpp:682-705` — worker slices receive disjoint stable event
      offsets, so scheduling and thread count cannot alter an event's random
      stream.
      > 2026-10-14 ewt: accepted. 8 threads x 4 events (4 idle workers) gives the
      > same events as 2 threads.
- [x] `src/Output.cpp:305-372` — output metadata records source control TOML,
      resolved physics choices, master seed, version, counts, and termination
      semantics needed to reproduce and interpret a run.
      > 2026-10-14 ewt: accepted with fixes. Resolved config now records the applied
      > defaults (z_effective, low_energy, straggling_z_effective, sim_step,
      > max_time, energy_fwhm), git_revision (flake self.shortRev or
      > dirtyShortRev; Makefile asks git, with a stamp so Output.o rebuilds),
      > catima_revision, and catima_global. Update runs now also require the
      > resolved configuration to match, so a file never mixes builds or
      > defaults. Tests: `standard` requires the keys and known revisions;
      > alter-resolved plus update must be refused.
- [x] `tests/root_checks.cpp:220-265` and `tests/run-tests.sh:91-97` — one-thread
      and multi-thread ROOT rows compare exactly by `(requested_strip,
      event_index)` and merged entry counts are checked.
      > 2026-10-14 ewt: accepted. Exact per-leaf comparison; duplicate keys fail;
      > production merge enforces strips x events.
<!---->
## Physics acceptance run
<!---->
- [x] Run `make -j2 check`; inspect the stopping, primary/later forbidden,
      multi-step-angle, boundary-convergence, and reproducibility fixtures.
      > 2026-10-14 ewt: run as `nix build && bash tests/run-tests.sh` (same script)
      > after every change in this review, 8/8 each time.
- [x] Run `make -j4 sanitize-check` and confirm the same physics fixtures pass
      with AddressSanitizer and UndefinedBehaviorSanitizer enabled.
      > 2026-10-14 ewt: `make -j8 sanitize-check` passes 8/8 with zero ASan/UBSan
      > reports, including all new fixtures, under STRICT_FLAGS (-Werror).
- [x] Compare representative production channels against an independent
      kinematics calculator and independently generated stopping tables.
      > 2026-10-14 ewt: kinematics are a permanent test (two-body-ground: 200
      > 87Rb(a,p)90Sr events vs closed-form E3(theta, T_beam), worst 8% of the
      > rounding tolerance; a 10 keV Q error fails). Stopping (one-off, catima
      > winger vs your SRIM tables, 400/555 Torr He): catima/SRIM is 0.89-0.95
      > above ~2.5 MeV/u and 1.03-1.15 below ~1 MeV/u for 87Rb, 90Sr, 37Cl,
      > 40K, 4He and p. Consistent with the 87Rb barometer note; `mean` sits
      > between. The data calibration remains the arbiter.
- [x] Compare sampled energy-loss distributions against direct ROOT Vavilov (or
      an experimental reference) at the actual kappa/beta2 grid used in
      production, including truncated tails and stopped steps.
      > 2026-10-14 ewt: KS check added to vavilov_statistics at kappa 1e-3, 0.03,
      > 0.3 and 1 (beta2 0.006/0.01), 100k draws each, D = 0.002-0.004 against the
      > 0.1% bound of 0.0062. Stopped steps are covered by the conservation and
      > stopped-before-vertex fixtures.
- [ ] After accepting all model decisions above, regenerate any result affected
      by the audit's rerun guidance before treating it as scientifically trusted.
<!---->
## Follow-ups from the 2026-10-14 review
<!---->
- [ ] Migrate the 53 controls in `~/MUSIC/analysis/*/sim_control` to schema 3:
      `[detector.eres]` → `resolution_fwhm_percent`, `n_events` →
      `events_per_strip`, add `seed`; (α,α) → ground/rutherford, (α,n)/(α,p) →
      uniform. Rerun all 87Rb and the 37Cl (α,α) simulations.
- [ ] Rerun anything with a particle above 10 MeV/u in the gas: the old
      `atima14` default underestimated mean stopping there. That covers the 16C
      examples (13.75 MeV/u beam) and 87Rb runs with fast recoil alphas or
      protons; the 87Rb and 37Cl beams themselves are unaffected. The ignored
      `37Cl_beam.toml` in the repo root still sets `atima14` and now fails to load.
- [ ] Sample a Tmax-truncated Landau for kappa < 1e-3 instead of freezing the
      kappa = 1e-3 Vavilov shape (removes the light-particle mean-loss bias).
- [ ] Replace the ungated per-event "reaction energetically not allowed" stdout
      warning (src/EventLoop.cpp) with a per-run count and one summary line.
- [ ] Style cleanup of the audit commit (0035262) against the AGENTS.md C++
      rules: `auto` 58 → 120, range-for 8 → 22, long lambdas (e.g. `loadProduct`),
      7 smart pointers (unique_ptr/make_unique; previously none). Remove the
      uncalled EnergyLoss::GetEnergyLoss and Particle::GetInitialEnergy /
      GetPathLength (plus their EnergyLoss targets).
      Also: one named amu→MeV constant (now hand-written in 8 places, two CODATA
      values). Also: `ResetBranches()` sets `beam_energy_reaction = 0` but the
      event loop sentinel is -2.
<!---->
