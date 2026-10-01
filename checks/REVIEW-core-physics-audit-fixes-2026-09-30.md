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
- [ ] **Accept the schema/API break before using these changes.** Reacted steps
      now require explicit excitation and angular models
      (`src/ControlFile.cpp:365-455, 556-570`), run controls require a seed and
      `events_per_strip` (`src/ControlFile.cpp:464-490`), ambiguous `eres` is
      rejected (`src/ControlFile.cpp:230-308`), and truth output is schema 3
      (`src/Output.cpp:195-255`, `src/Output.cpp:329-363`). Update control
      generators and analysis consumers first.
- [ ] **Accept the low-kappa approximation.** For kappa below ROOT's accurate
      range, the implementation freezes the shape at kappa `1e-3` rather than
      attempting an asymptotic Landau parameter conversion
      (`src/VavilovSampler.cpp:74-117`). Confirm this is adequate for the
      experiment's step sizes and ions.
- [ ] **Accept the excitation models.** `ground`, `uniform`, and `forced` remain
      continuum choices rather than a nuclear level/width/branching model
      (`src/Kinematics.cpp:121-145`; limitations at `README.md:343-358`). Review
      every shipped channel choice independently.
- [ ] **Accept layer-scale semantics.** A layer `dedx_scale` is implemented as
      equivalent thickness, so it scales both mean loss and straggling
      (`src/Materials.cpp:138-179`). Confirm this remains consistent with the
      calibrations that introduced those factors.
<!---->
## Reaction graph and kinematics
<!---->
- [ ] `src/ControlFile.cpp:365-459` — every reaction entry is a complete table
      with named ejectile/residue, explicit per-step excitation and angular
      policies, bounded Rutherford cutoff, and positive stopping scales.
- [ ] `src/ControlFile.cpp:583-624` — every named nuclide exists; beam plus target
      matches the compound; every parent decomposes into products with conserved
      A and Z; identity elastic scattering is constrained to ground-state plus
      Rutherford sampling.
- [ ] `src/Simulator.cpp:392-417` — recursive `minEx` is the minimum excitation
      needed to feed all subsequently forced steps. Check signs and indexing for
      exothermic and endothermic chains.
- [ ] `src/Kinematics.cpp:3-12, 34-69` — beam momentum, total four-momentum, CM
      boost, and compound excitation use relativistic energy-momentum relations
      and the intended atomic-mass convention.
- [ ] `src/Kinematics.cpp:95-145` — invariant-mass availability is authoritative,
      the `1e-9 MeV` tolerance accepts threshold equality, `ground` gives zero
      excitation, `uniform` samples `[0,Eavail]`, and `forced` is limited to
      intermediate residues with enough energy for the remaining chain.
- [ ] `src/Kinematics.cpp:149-172` — isotropic sampling is uniform in
      `cos(theta)` and azimuth; Rutherford inverse-CDF sampling and its minimum
      angle are correct; each step's sampled truth angle is stored immediately.
- [ ] `src/Kinematics.cpp:182-233` — the two-body momentum radicand uses
      `(mh + Ex)` consistently, products are opposite in the CM frame, the lab
      boost sign is correct, and each heavy product becomes the next decay
      parent's four-vector.
- [ ] `src/Kinematics.cpp:76-84, 242-274` — reaction construction is atomic:
      products remain disabled until the full chain succeeds, every product and
      truth field is cleared on failure, and only the final residue is eligible
      for transport.
- [ ] `src/EventLoop.cpp:164-213` — `beam_energy_reaction` and the vertex are
      recorded only for a realized reaction; an attempted but forbidden chain
      keeps the N/A sentinels.
- [ ] `tests/root_checks.cpp:113-169` and
      `tests/controls/later-step-forbidden.toml:26-47` — independently verify
      that both primary and later-step failures retain no product energies,
      angles, exit states, or residue identity while the beam continues.
- [ ] `tests/root_checks.cpp:172-193` — strengthen or accept the current
      multi-step-angle test. It proves independent stored values, but a deeper
      validation would reconstruct the four-vectors and verify that each stored
      angle is the one actually used for that step.
<!---->
## Gas transport and detector boundaries
<!---->
- [ ] `include/Simulator.hpp:114-141` — the termination enum and result object
      distinguish not-propagated, reached-vertex, stopped, downstream exit,
      upstream exit, side exit, timeout, and invalid state without overloading
      an energy sentinel.
- [ ] `src/Propagation.cpp:28-95` — null/disabled particles and invalid time/step
      inputs terminate safely; initial energy, rest energy, direction, and
      outside-volume classification are physically consistent.
- [ ] `src/Propagation.cpp:97-196` — each iteration selects the first of the time,
      requested-vertex, detector-face, row, and column boundaries. Check all
      signs for forward, backward, and transverse tracks, including a start
      exactly on a boundary.
- [ ] `src/Propagation.cpp:198-240` — the chosen fractional distance is passed to
      stopping/straggling, its complete deposit is assigned to the traversed
      cell, position/time advance by the same distance, and a stopping step is
      accounted before returning.
- [ ] `src/EnergyLoss.cpp:251-324` — a sampled loss is never negative or above
      the incoming KE; a stopping step deposits all remaining KE and shortens
      its path proportionally. Review the linear-in-step approximation used to
      locate that stop.
- [ ] `src/EventLoop.cpp:168-213, 277-312` — a reaction is attempted only after
      exact vertex arrival; a rejected reaction resumes the same beam from that
      state without resetting already accumulated deposits.
- [ ] `src/Propagation.cpp:372-438` — only a downstream exit crosses the exit
      foil, oblique tracks use `thickness / abs(cos(theta))`, neutrals retain
      their KE, and the final residue is selected from reaction state rather
      than inferred from transport success.
- [ ] `src/Propagation.cpp:440-453` and `src/Output.cpp:266-303` — termination,
      stop-position, residue-step, and energy sentinels are internally
      consistent for reacted, rejected, stopped, escaped, timed-out, and invalid
      events.
- [ ] `tests/root_checks.cpp:85-110` — with stochastic effects and windows off,
      gas-entry KE equals cathode plus upstream/downstream dead-layer deposits
      within `2e-5 MeV` for a stopping event.
- [ ] `tests/root_checks.cpp:196-218` — coarse and fine step runs hit the exact
      downstream surface and converge in exit energy. Decide whether additional
      upstream, side, timeout, stopped-before-vertex, and oblique-exit fixtures
      are required before production.
<!---->
## Materials, windows, and detector response
<!---->
- [ ] `src/Materials.cpp:3-75` — pressure/temperature units and ideal-gas density
      are correct; P10 uses its mole-fraction molar mass for density but the
      stoichiometric mass for catima component weight fractions.
- [ ] `src/Materials.cpp:77-136` — elemental/compound definitions, densities,
      Havar mass fractions, and conversions between mg/cm2 and micrometres are
      correct for the physical foils in use.
- [ ] `src/Materials.cpp:138-183` — disabled layers, positive scale factors, and
      equivalent-thickness handling agree for entrance, exit, and degrader
      layers.
- [ ] `src/Materials.cpp:185-226` — catima receives kinetic energy per nucleon,
      returns energy/sigma per nucleon, and sampled passive-layer output is
      bounded to `[0, Ein]`.
- [ ] `src/EventLoop.cpp:143-162` — the event energy chain is accelerator spread,
      then degrader, then entrance window, then gas entry; the same event RNG is
      used for each stochastic term.
- [ ] `src/Propagation.cpp:266-353` — detector noise uses exactly one selected
      convention: absolute Gaussian sigma in MeV or relative percent FWHM. Check
      that anode channels are independent, the cathode receives one independent
      fluctuation after physical summing, and dead layers remain noiseless.
- [ ] `ControlExamples/37Cl_alpha_n/37Cl_alpha_alpha.toml:6-10, 29-83` — verify
      the explicit beam energy, absolute noise table, and elastic reaction
      choices against experiment-specific intent; repeat this review for every
      migrated control under `ControlExamples/`.
<!---->
## Stopping tables, effective charge, and Vavilov sampling
<!---->
- [ ] `include/EnergyLoss.hpp:17-32` — mean and straggling catima configurations
      are immutable per simulator instance, and the process-global catima cache
      is protected by one mutex.
- [ ] `src/EnergyLoss.cpp:93-193` — mean stopping and variance tables use the
      intended catima configurations, units, reference thickness, interpolation
      grid, SRIM replacement/mean policy, and user-scale ordering.
- [ ] `src/EnergyLoss.cpp:225-243` — gas `Z/A` is a mass-fraction sum and remains
      valid for all supported elemental and compound gases.
- [ ] `src/EnergyLoss.cpp:259-308` — catima variance scales as path length, while
      gas `dedx_scale` changes the mean only. Confirm that separation is the
      intended calibrated model.
- [ ] `src/EnergyLoss.cpp:267-307` — verify dimensions and constants in xi,
      maximum transferable electron energy, and kappa. For mixtures, confirm the
      additive `weight_fraction * z_eff^2 * Z/A` calculation and that
      `physics_.straggling` is the same effective-charge prescription used by
      catima's variance.
- [ ] `src/VavilovSampler.cpp:41-72` — the `(log kappa, beta2, probability)` grid
      covers the declared domain and standardizes every ROOT quantile with
      ROOT's analytic mean and variance.
- [ ] `src/VavilovSampler.cpp:13-32, 74-105` — interpolation clamps kappa/beta2,
      extrapolates only the probability tails to avoid endpoint point masses,
      and indexes all eight trilinear corners correctly.
- [ ] `src/VavilovSampler.cpp:107-117` — low-kappa requests use the kappa `1e-3`
      Vavilov shape and high-kappa requests use a standard Gaussian; no split
      convolution component can silently vanish.
- [ ] `tests/vavilov_statistics.cpp:11-70` — review the kappa/beta2/probability
      grid and tolerances against direct `ROOT::Math::VavilovAccurate`
      quantiles, plus Gaussian-limit moments and deterministic low-limit
      behavior. Decide whether production acceptance also needs larger sampled
      moment/tail tests for the experiment's exact ion/pressure/step regime.
- [ ] `README.md:244-258` — documentation accurately states boundary accounting,
      low/high-kappa behavior, analytic standardization, physical truncation,
      and the fact that sampled moments are close to—not exactly—the catima
      moments.
<!---->
## Determinism, concurrency, and scientific provenance
<!---->
- [ ] `src/ControlFile.cpp:310-360, 633-647` — physics selections are parsed into
      a temporary instance-owned configuration and installed only after the
      whole control validates; no worker mutates catima defaults.
- [ ] `src/EnergyLoss.cpp:118-135, 356-371`, `src/Materials.cpp:187-226`, and
      `src/Simulator.cpp:244-274` — every catima calculation/cache access is
      serialized, including cache misses after prewarming.
- [ ] `src/Simulator.cpp:421-435` and `src/EventLoop.cpp:108-110` — each event RNG
      stream is a deterministic function of master seed, requested strip, and
      per-strip event index; TRandom3 seed zero is avoided.
- [ ] `src/Simulator.cpp:682-705` — worker slices receive disjoint stable event
      offsets, so scheduling and thread count cannot alter an event's random
      stream.
- [ ] `src/Output.cpp:305-372` — output metadata records source control TOML,
      resolved physics choices, master seed, version, counts, and termination
      semantics needed to reproduce and interpret a run.
- [ ] `tests/root_checks.cpp:220-265` and `tests/run-tests.sh:91-97` — one-thread
      and multi-thread ROOT rows compare exactly by `(requested_strip,
      event_index)` and merged entry counts are checked.
<!---->
## Physics acceptance run
<!---->
- [ ] Run `make -j2 check`; inspect the stopping, primary/later forbidden,
      multi-step-angle, boundary-convergence, and reproducibility fixtures.
- [ ] Run `make -j4 sanitize-check` and confirm the same physics fixtures pass
      with AddressSanitizer and UndefinedBehaviorSanitizer enabled.
- [ ] Compare representative production channels against an independent
      kinematics calculator and independently generated stopping tables.
- [ ] Compare sampled energy-loss distributions against direct ROOT Vavilov (or
      an experimental reference) at the actual kappa/beta2 grid used in
      production, including truncated tails and stopped steps.
- [ ] After accepting all model decisions above, regenerate any result affected
      by the audit's rerun guidance before treating it as scientifically trusted.
<!---->
