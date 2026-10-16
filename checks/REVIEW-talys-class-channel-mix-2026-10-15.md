# Review checklist — TALYS classes: each MUSIC class as TALYS's channel mix, 2026-10-15
<!---->
Per the AI-disclosure section of the README, changes to simulation physics
(excitation sampling, branch weighting, beam energies) must be human-reviewed.
This covers the **uncommitted** change on top of ee3e0d2, planned in
`PLAN-talys-channel-mix-2026-10-15.md`. A MUSIC class (an, ap, aa) becomes
every TALYS exclusive channel ending in one element, weighted per strip by
TALYS's cross sections, each event running one channel's decay chain. On the
TALYS side, `talys-cache` now rebuilds TALYS's full decay cascade; on the
simulator side, `[reaction] talys_class` replaces the two-branch (a,a') form.
It also fixes the per-event seeding, so every run's events change. Line
numbers are for the working tree on 2026-10-15.
<!---->
Written by Claude Code (Opus 5.5). Verified so far: `nix build && bash
tests/run-tests.sh` passes all 10 stages, as do `make check` and `make
sanitize-check` (no ASan/UBSan reports); `make strict`, `make format-check`,
and `git diff --check` are clean. The two reference controls that pin the
pre-branch random draws were regenerated from fea3f4b itself with only the
new seeding patched in, and the working tree reproduces them, so controls
without `talys_class` still draw exactly as fea3f4b does.
<!---->
## Decisions made in conversation (2026-10-15) — confirm
<!---->
- [x] **A class is the final residue's element.** It is every exclusive
      channel ending there, since MUSIC sees only the residue's charge, and
      each strip weights them by TALYS's exclusive cross sections at its
      centre energy (`src/StripTables.cpp:209-218`).
      > 2026-10-15 ewt: accepted.
- [x] **Each event records its channel.** `reaction_channel` in the MC tree,
      and `metadata/channels` says what each value is
      (`src/Output.cpp:505-516`).
      > 2026-10-15 ewt: accepted.
- [x] **Small channels are left out by two thresholds**
      (`include/TalysTables.hpp:17-20`, `src/StripTables.cpp:209-218`): a
      channel is kept if, at one selected strip, it reaches both
      `min_channel_fraction` of the class's cross section (default 0.001) and
      `min_channel_mb` (default 1e-4 mb).
      > 2026-10-15 ewt: accepted.
- [x] **The full TALYS cascade decides each channel**, not each state's own
      gamma fraction. Phase 0 on 87Rb at E_cm 14 MeV: (a,n) 158.0441 mb from
      the cascade against TALYS's 158.046; the state's own gamma fraction
      gave +0.15%, and cutting at S_n gave −53%.
      > 2026-10-15 ewt: accepted.
- [x] **`talys_class` replaces the two-branch (a,a') form**, which is now an
      error pointing to it (`src/ControlFile.cpp:107-134`).
      > 2026-10-15 ewt: accepted.
- [x] **Both emission orders are kept** as separate channels, (a,np) and
      (a,pn), each with TALYS's own flows.
      > 2026-10-15 ewt: accepted.
- [x] **One TALYS run per strip, at its centre**, as before. Within a strip,
      a channel's energy dependence comes only from its states closing as the
      vertex energy falls below the centre. Near a threshold this is a large
      effect: at 87Rb strip 16 (E_cm 3.93 MeV), (a,n) is closed below E_cm
      3.75, so (a,g) runs in 35% of the Y class's events against 18.5% of the
      centre cross section, both about 1e-9 mb there.
      > 2026-10-15 ewt: accepted.
- [x] **Per-event seeds are fixed** (`src/Simulator.cpp:423-449`). The old
      `EventSeed` XORed the master seed with the event index before one mix,
      so seeds differing in low bits gave the same streams a few indices
      apart: seeds 104 and 105 replayed 784 of 800 events of
      `tests/controls/branches.toml`. The seed, strip, and event index are
      now each mixed in turn. Every run's events change.
      > 2026-10-15 ewt: accepted.
- [x] **Disk and time.** TALYS's decay output is about 2.1 GB for 87Rb per
      strip, streamed and deleted; a TALYS run takes about 17 s.
      > 2026-10-15 ewt: accepted.
<!---->
## Physics choices I made — confirm or overturn
<!---->
- [x] **`min_channel_mb` defaults to 1e-6 mb (1 nb)**, far below what MUSIC
      can measure, so it removes only channels that matter nowhere: 87Rb's aa
      class no longer keeps (a,ga), which had passed the relative threshold
      only at strip 17, where the whole class is 3e-23 mb. Both thresholds
      must hold at the same strip. A class with no channel reaching 1 nb at
      any selected strip is an error naming both thresholds
      (`src/StripTables.cpp:293-299`). Check: the default value.
      > 2026-10-15 ewt: make it 100 nb. Done: the default is now 1e-4 mb
      > (`include/TalysTables.hpp:20`, README). The real 87Rb and 37Cl runs
      > keep the same channels. At 87Rb 330 MeV the Sr class falls under
      > 100 nb from strip 10 (62 nb), where a single-strip ap control would
      > now be refused. The MUSIC controls' 401 MeV beam puts their strips at
      > higher E_cm, so they should stay above it (not run).
- [x] **Floors** (`src/StripTables.cpp:355-381`). A step's excitation is
      drawn no lower than the excitation from which the rest of its chain can
      happen: the next emission's separation energy plus the lowest state the
      next step's table can reach, recursively. Without it, a TALYS continuum
      bin straddling a threshold puts events below it, and their next step
      fails. Check: in the suite's (a,2n), the lowest step-0 excitation drawn
      is 7.79956 MeV, 40K's neutron separation energy, and no event goes
      unreacted.
      > 2026-10-15 ewt: accepted.
- [x] **A bin straddling a floor counts pro rata**
      (`src/ExcitationTable.cpp:22-31`), as a bin straddling the energy
      available already did. The open fraction is taken above the floor and
      divided by the table's weight above it (`src/StripTables.cpp:371-373,
      397-433`), so a fully open channel weighs its cross section.
      > 2026-10-15 ewt: accepted.
- [x] **One draw picks the channel**, by cross section times step 0's open
      fraction at the vertex (`src/StripTables.cpp:397-433`). Later steps'
      openness is not in the weight; with the floors, a later step can fail
      only on an exact edge, and the event is then unreacted, like any
      rejected chain.
      > 2026-10-15 ewt: accepted.
- [x] **Gamma steps** (`src/StripTables.cpp:303-317`). A channel's leading
      gamma ((a,g), and the gamma-first (a,gn), (a,gp), ...) is emitted as a
      massless, chargeless particle, transported like a neutron. Gammas within
      a residue are not emitted (next item). An (a,g) residue keeps the
      excitation of the state the first gamma feeds.
      > 2026-10-15 ewt: accepted.
- [x] **A residue's own gamma cascade is not simulated.** Each step's table
      is the residue's population right after its emission. If it then
      gamma-decays before the next emission, the next step starts from the
      pre-gamma excitation, so that gamma's energy goes into the next
      particle. Each step's excitation is otherwise uncorrelated with the
      previous step's state, and TALYS's spins and parities are used only to
      build the tables. These are in the README's limitations.
      > 2026-10-15 ewt: accepted.
- [x] **Non-elastic angles are isotropic**; elastic θ_cm comes from TALYS's
      distribution as before (`src/Kinematics.cpp:197-216`).
      > 2026-10-15 ewt: accepted.
- [x] **Elastic is only in the beam's class**, and is excluded from the class
      cross section the relative threshold compares against
      (`src/StripTables.cpp:231-246`). `theta_cm_min_deg` is refused in any
      other class (`src/ControlFile.cpp:1018-1031`).
      > 2026-10-15 ewt: accepted.
- [x] **Stopped strips** (`src/StripTables.cpp:168-178, 238-241`). A strip
      past the mean beam's range has no table: elastic with pure Rutherford
      angles in the beam's class, no reaction in any other. Check: the suite's
      30 MeV beam with 6 MeV FWHM reacts elastically in strips 15 and 16.
      > 2026-10-15 ewt: accepted.
- [x] **Left-out cross section** (`src/StripTables.cpp:276-292`) is recorded
      per strip as `dropped_mb`, and warned about when, summed, it would have
      passed both thresholds. 87Rb's Sr class at strip 0 drops 0.0009 mb of
      0.73 mb, (a,np) and the like, and is warned about; the 3.7e-26 mb (a,ga)
      leaves at strip 17 of the aa class is not.
      > 2026-10-15 ewt: accepted.
- [x] **A channel whose chain cannot finish at a strip is left out there**,
      with a warning (`src/StripTables.cpp:374-380`).
      > 2026-10-15 ewt: accepted.
- [x] **One stopping scale per kind of product**
      (`src/ControlFile.cpp:628-654`): `residue_dedx_scale` for every
      residue, `evap_dedx_scale` for every light product, both default 1. The
      beam-like residue of elastic and (a,a') uses `residue_dedx_scale`, not
      `beam.dedx_scale`.
      > 2026-10-15 ewt: accepted.
<!---->
## Found during the work
<!---->
- [x] **TALYS exits with status 0 after its own errors**
      (`src/talys-cache.cpp:159-168`). With alphaomp 9 and no potential file
      it prints `TALYS-error: ...` and exits 0. `talys-cache` now requires
      TALYS's "successful calculation" line and reports the error line.
      > 2026-10-15 ewt: accepted.
- [x] **`talys-cache` makes the products' SRIM tables**
      (`src/talys-cache.cpp:297-315`). Only the TALYS tables name a class's
      products, and the MUSIC controls use `stopping = "mean"`, which needs a
      SRIM table per ion. `srim-cache`'s generator moved to
      `include/SrimGenerator.hpp` so both tools share it.
      > 2026-10-15 ewt: accepted.
<!---->
## Code — review these lines
<!---->
- [x] `src/TalysCascade.cpp:50-75` — each state's fate, bottom-up: the
      probability of each ordered emission sequence before it stops.
      > 2026-10-15 ewt: accepted.
- [x] `src/TalysCascade.cpp:111-233` — streams TALYS's population file and
      handles its three traps: population blocks include feeding from other
      nuclei (`:214-231`); the per-bin table continues headerless between
      states, so only `entries` rows are read (`:166-169`); decay columns are
      in the source's spin convention (`:192-205`).
      > 2026-10-15 ewt: accepted.
- [x] `src/TalysCascade.cpp:235-320` — every ordered channel with its cross
      section, and per step the residue's population right after that
      emission, carried down each residue's gamma cascade into the next
      emission.
      > 2026-10-15 ewt: accepted.
- [x] `src/TalysCascade.cpp:322-414` — reads TALYS's section 6a and checks
      each channel against it to 0.1%, and each step's sum against the
      channel. Check: the 5e-6 mb allowance is TALYS printing the gamma-first
      total `f12.5`.
      > 2026-10-15 ewt: accepted.
- [x] `src/TalysCascade.cpp:416-512` — writes and reads a table's `channel`
      and `step` rows: levels are bins up to NL, continuum bins span Ex ±
      half the bin size.
      > 2026-10-15 ewt: accepted.
- [x] `src/TalysTables.cpp:55-60, 100-104` — `channels y` and `outdecay y`
      added to every TALYS input, and refused as `[talys]` keywords.
      > 2026-10-15 ewt: accepted.
- [x] `src/talys-cache.cpp:152-223` — reads TALYS's output, requires its
      sign-off, builds and checks the channels, writes the table.
      > 2026-10-15 ewt: accepted.
- [x] `src/ExcitationTable.cpp:22-31, 138-208` — floored allowed weight and
      sampling; the unfloored forms call them with floor 0, giving identical
      draws. `LowestAbove` finds a floor's lowest reachable state.
      > 2026-10-15 ewt: accepted.
- [x] `src/ControlFile.cpp:69-78, 107-134, 238-291, 623-661, 852, 893-902,
      1018-1031` — the element's Z from the mass table; `level_table`
      required again; `min_channel_fraction` and `min_channel_mb`; the
      `[reaction]` class keys; class runs need no steps, a reacted run, no
      event display, a `4He` target, and an element at or below the
      compound's.
      > 2026-10-15 ewt: accepted.
- [x] `src/StripTables.cpp:147-301` — reads each strip's table once (the
      master's copy in workers), checks its input, selects the class's
      channels, and builds elastic.
      > 2026-10-15 ewt: accepted.
- [x] `src/StripTables.cpp:303-395` — each channel's own particles per step
      and the idle particles for slots past a chain. Check: a channel's
      steps never share a particle object, so (a,2n)'s two neutrons are
      transported separately.
      > 2026-10-15 ewt: accepted.
- [x] `src/Kinematics.cpp:87-107` — a class event draws its channel and
      points the step slots at its chain, idle particles after it.
      > 2026-10-15 ewt: accepted.
- [x] `src/Kinematics.cpp:122, 153-190` — the loop runs the chain's steps; a
      class step samples its own table above its floor, an elastic step is
      ground state; other controls take the unchanged branch path.
      > 2026-10-15 ewt: accepted.
- [x] `src/Kinematics.cpp:311-320` — slots past the chain stay
      non-propagating, only the last residue is transported, and
      `reaction_channel` is set.
      > 2026-10-15 ewt: accepted.
- [x] `src/Simulator.cpp:103-104, 396-397` and `src/ParticleSetup.cpp:60, 77`
      — particles are owned by one list, so a slot pointing at a channel's
      particle is never freed twice; `forced`'s minimum excitations are
      skipped for a class.
      > 2026-10-15 ewt: accepted.
- [x] `src/Simulator.cpp:423-449` — the seed, strip, and event index each
      mixed in turn by SplitMix64's finaliser.
      > 2026-10-15 ewt: accepted.
- [x] `src/Output.cpp:69, 167, 222, 294, 380` — the `reaction_channel`
      branch, also checked and bound on update runs; schema 6.
      > 2026-10-15 ewt: accepted.
- [x] `src/Output.cpp:436-477` — the class's resolved configuration: its
      settings, each channel (`reaction_channel_N`), and per strip the beam
      energy, E_cm, channel cross sections, left-out cross section, and table
      MD5 (`talys_strip_N`). An update run refuses any change to them.
      > 2026-10-15 ewt: accepted.
- [x] `src/Output.cpp:479-516` — channel descriptions, e.g.
      "87Rb(a,2n)89Y: 1n + 90Y, then 1n + 89Y", and the legend.
      > 2026-10-15 ewt: accepted. (The function's comment wrote the neutron
      > as "n"; it now says "1n", as the output does.)
<!---->
## Schema, build, and docs
<!---->
- [x] `include/Simulator.hpp:63-73, 233-244, 363-364, 484-490` — the channel
      struct, step-slot ownership, the new branch, and the class settings.
      > 2026-10-15 ewt: accepted.
- [x] `include/SrimGenerator.hpp:1-149`, `src/srim-cache.cpp:1-90` — the
      SRIM generator, moved unchanged in behaviour.
      > 2026-10-15 ewt: accepted.
- [x] `Makefile:71, 88, 95-98, 127-132, 153, 203-208, 213-215, 219-221,
      240` — TalysCascade in the TALYS checks, the generator for
      `talys-cache`, a sanitized `srim-cache`, and `SRIMCACHE_BIN` for the
      suite.
      > 2026-10-15 ewt: accepted.
- [x] `README.md:48-51, 190-193, 235, 276, 304-427, 562, 587, 596, 624-644,
      658-677` — seeds, the "TALYS classes" section (replacing "Per-strip
      TALYS tables"), the truth branch, the metadata and update paragraphs,
      and two limitations.
      > 2026-10-15 ewt: accepted.
- [x] `flake.nix:31` — version stays `26.10.15` if committed today (a
      same-day commit keeps the date, as 042674f did).
      > 2026-10-15 ewt: accepted.
<!---->
## Tests
<!---->
- [x] `tests/fake-talys.sh` writes a small cascade in TALYS's formats with
      all three traps: (a,g), (a,gn), (a,n), (a,2n) above E_cm 12 MeV, (a,p),
      (a,a'). Its (a,g) state sits at the strip centre's excitation, so it is
      open only upstream of each centre; its gamma-then-neutron state is above
      41K's neutron threshold. A `fakeerror` keyword fails as TALYS does.
      > 2026-10-15 ewt: accepted.
- [x] `tests/talys_tables.cpp:127-210` — the fake's channels and per-step
      tables against hand-worked values, 6a mismatches, and table round
      trips and refusals.
      > 2026-10-15 ewt: accepted.
- [x] `tests/excitation_table.cpp:138-171` — floored weights, sampling within
      the floor and its uniformity, `LowestAbove`, and floor 0 drawing
      exactly as before.
      > 2026-10-15 ewt: accepted.
- [x] `tests/fixtures/legacy-draws-chain.angles`,
      `tests/fixtures/legacy-draws-elastic.angles` — regenerated from fea3f4b
      with only the new `EventSeed`, as their headers say
      (`tests/run-tests.sh:190-200`).
      > 2026-10-15 ewt: accepted.
- [x] `tests/run-tests.sh:49-60` — ten invalid class and `[talys]` controls,
      each rejected for its own reason.
      > 2026-10-15 ewt: accepted.
- [x] `tests/run-tests.sh:204-330` with `tests/root_checks.cpp:946-1240`:
      - aa class (37Cl, 92 MeV): elastic and (a,a') weights against each
        stored table, draw fractions, elastic angles, vertex energies, the
        closed strip warned about.
      - aa class with a stopping beam: warnings, and beam that straggles past
        the mean range reacting elastically (`class-channels unreacted-ok`).
      - an class (37Cl, 135.5 MeV, 12,000 events), `class-channels`: every
        event reacts through a listed channel; chain lengths, sentinels past
        each chain, the transported residue, and every excitation above its
        next emission's threshold; draw fractions within 5σ of the recorded
        cross sections per strip. (a,g) and (a,2n) dip below the 1.2%
        threshold at some strips and are kept; (a,2n) is unreachable at
        strip 5 and warned about.
      - thresholds: 0.6 mb drops (a,g), 0.5 mb, with no warning; at 6.5%,
        (a,g) and (a,gn) are each left out but together warned about.
      - a failed TALYS run is refused.
      - SRIM stopping (`tests/fake-srim.sh`): `srim-cache` makes the beam and
        target tables, `talys-cache` the products' (39K, 40K, 41K), and
        `musicsim` reads them; a 2% threshold drops (a,g) and records its
        0.5 mb.
      > 2026-10-15 ewt: accepted.
- [x] `tests/run-tests.sh:339-347` with `tests/root_checks.cpp:1256-1285` —
      seeds 104 and 105 share at most 1% of their events.
      > 2026-10-15 ewt: accepted.
- [x] Mutation-checked in a scratch copy, each mutant built with `nix build`
      and run through the suite on its `result/bin`. Each of these 24 fails
      at least one check:
      - the channel draw: no open fraction; no normalisation above the floor;
        another strip's cross sections; every channel recorded as 0;
        unreacted events recorded as channel 0;
      - floors: none; sampling below the floor; a straddling bin not cut by
        the floor; the bin weight's denominator;
      - thresholds: no relative threshold; no absolute threshold; the class
        total over all channels; kept only if passing at every strip; the
        left-out cross section not counted; a left-out sum warned about below
        `min_channel_mb`; no left-out warning;
      - transport: idle slots transported; the wrong residue transported;
      - elastic: weight at 1° instead of the cut; pure Rutherford angles; a
        stopped strip's weight 0;
      - tools: no products' SRIM tables; no TALYS sign-off check;
      - seeding: the old one, caught by the pinned draws and by the
        independence check alone (784 of 800 events replayed).

      Not covered: a class chosen by the first residue's element instead of
      the last, since no fake channel changes element mid-chain (the real runs
      below do); the two thresholds met at different strips. The TALYS side
      was mutation-checked in Phase 1 (9 mutants, in the plan), with the fake
      as it was before its (a,g) and gamma-then-neutron states moved.
      > 2026-10-15 ewt: accepted.
- [x] Real TALYS, by hand (not in the suite), 2000 events per strip:
      - 87Rb at 330 MeV, strips 0–17: Y class (a,g), (a,n), (a,gn), (a,2n);
        at strip 0, (a,2n) is 58.7% of the class and 57.9% of events. Sr
        class (a,p) alone at the default thresholds; with none, ten channels,
        including (a,np) as n + 90Y then p + 89Sr, and the three-step
        (a,gnp). Rb class elastic and (a,a').
      - 37Cl at 92 MeV: K class (a,g), (a,n); Ar class (a,p); Cl class
        elastic, (a,a'), (a,ga).
      > 2026-10-15 ewt: accepted.
