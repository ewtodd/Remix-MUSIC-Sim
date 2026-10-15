# Review checklist — per-strip TALYS tables and computed weights, 2026-10-15
<!---->
Per the AI-disclosure section of the README, changes to simulation physics
(excitation sampling, branch weighting, beam energies) must be human-reviewed.
This covers the **uncommitted** change on top of 042674f that replaces the
plan's MUSIC-repo producer (plan §6): TALYS is a flake input here, a
`talys-cache` pre-pass builds one table per strip, and a step with such a
table gets its branch weights per strip. Line numbers are for the working tree
on 2026-10-15. Checks indicate review by ewt.
<!---->
Written by Claude Code (Opus 5.5). Verified so far: `nix build && bash
tests/run-tests.sh` passes all 10 stages, as do `make check` and `make
sanitize-check` (no ASan/UBSan reports); `make strict`,
`make format-check`, and `git diff --check` are clean. The two reference
controls that pin the pre-branch random draws (legacy-draws-chain and
legacy-draws-elastic) still reproduce fea3f4b.
<!---->
## Decisions made in conversation (2026-10-15) — confirm
<!---->
- [x] **TALYS runs here, like SRIM.** `talys-nix` is a flake input locked
      at MUSIC's revision a7f8ba3 with its own nixpkgs, as `srim-nix` is
      (`flake.nix:16, 35-36`). A
      `talys-cache <control>` pre-pass builds missing tables and `musicsim`
      only reads them (`src/talys-cache.cpp:192-236`), so parallel batch runs
      never race to build one.
      > 2026-10-15 ewt: accepted.
- [x] **One table per strip**, at the strip centre's mean beam energy
      (`src/StripTables.cpp:45-96`).
      > 2026-10-15 ewt: accepted; sub-strip bins only if per-vertex accuracy matters later.
- [x] **Weights computed per strip** (`src/StripTables.cpp:126-144`): elastic
      = Rutherford above the branch's `theta_cm_min_deg`, inelastic = the
      strip table's TALYS cross section, both in mb.
      > 2026-10-15 ewt: accepted.
- [x] **`[talys]` block**, recorded through each table's TALYS input
      (`src/ControlFile.cpp:259-292`). `alpha_omp` defaults to 9 (Atomki-V2)
      when a potential exists for the beam nucleus, else 2 (McFadden-Satchler)
      (`src/StripTables.cpp:57`). Atomki-V2 has 87Rb, not 37Cl. Extra
      `keywords` follow MUSIC's shared defaults (`src/TalysTables.cpp:65-72`).
      > 2026-10-15 ewt: accepted.
<!---->
## Physics choices I made — confirm or overturn
<!---->
- [x] **The elastic weight is pure Rutherford** above `theta_cm_min_deg`.
      Measured at the most upstream strip, true elastic over Rutherford above
      θ_min is:

      | θ_min | 37Cl, E_cm 8.65 | 87Rb, E_cm 14 |
      |---|---|---|
      | 10° | 0.98 | 0.99 |
      | 30° | 0.79 | 0.92 |
      | 60° | 0.39 | 0.46 |
      | 120° | 0.17 | 0.09 |

      > 2026-10-15 ewt: overturned — do it right. Elastic now comes from
      > TALYS (next item).
- [x] **Elastic from TALYS** (`src/TalysTables.cpp:274-445`,
      `src/StripTables.cpp:150-171`, `src/Kinematics.cpp:148-151`).
      `talys-cache` adds `outangle y` and `fileelastic y` and reads
      `aa*ang.L00`, the (a,a) file (`src/talys-cache.cpp:152`); a plain
      `*ang.L00` would also match the (a,g) file.
      - **What is stored:** `elastic θ ratio` rows on TALYS's 2° grid, where
        ratio = TALYS's shape-elastic ratio + compound elastic ÷ our
        Rutherford (`:339-352`). TALYS's own "xs/Rutherford" column is shape
        elastic only. Compound elastic (1.9 mb/sr, isotropic, for 37Cl at
        8.65 MeV) also leaves the residue in its ground state, so it is
        elastic to the detector.
      - **Frame/energy check** (`:322-336`): TALYS's shape elastic ÷ its ratio
        must equal our Rutherford to 1% between 30° and 150°. On real output
        it is a constant 1.0019 (37Cl) and 1.0034 (87Rb) at every angle —
        TALYS's own constants or E_cm rounding.
      - **Weight and sampling:** the ratio is held at its cell mean, with
        exact Rutherford inside each cell. The weight is therefore exactly
        the integral of what is sampled, and θ still takes one draw.
      - **Inelastic angles** stay isotropic.
      > 2026-10-15 ewt: accepted.
- [x] **Inelastic = the residue's population right after binary alpha
      emission** (TALYS `outpopulation`, first block of the beam nucleus's
      section, `src/TalysTables.cpp:146-272`). Discrete levels become `level`
      rows, continuum bins `bin` rows. The bin centre is TALYS's Ex, and the
      width comes from the header. The ground-state row is dropped
      (`:230`): it is gamma-then-alpha feeding, not compound elastic. The sum
      must match TALYS's binary alpha cross section to 0.1% (`:250`). In the
      spike runs it matched to ~5e-6 from 14 mb down to 6e-19 mb, the
      residue being small gamma-then-alpha feeding of excited levels.
      > 2026-10-15 ewt: accepted.
- [x] **A closed strip gets an empty table and weight 0**
      (`src/TalysTables.cpp:162`, `src/ExcitationTable.cpp:124-132`). Below
      the first excited level TALYS reports a binary alpha cross section of 0.
      > 2026-10-15 ewt: accepted, provided the user is told, for example when
      > a beam is meant to stop. Now warned about (next item).
- [x] **Warnings** (`src/StripTables.cpp:113-176`, `src/Output.cpp:444-451,
      494-500`, `src/talys-cache.cpp:170-180, 228-266`). Three cases are
      printed by `talys-cache` as it runs, printed by `musicsim` at startup
      as `musicsim WARNING:` lines (master only), and kept in
      `metadata/warnings`:
      - a strip with no (a,a') channel open;
      - a strip where TALYS gives no elastic distribution. Real TALYS writes
        all zeros from E_cm ≈ 0.09 MeV down; elastic is then pure Rutherford,
        which it physically is there.
      - a strip whose centre is beyond the mean beam's range. It gets no
        TALYS run (TALYS cannot run at E = 0), inelastic weight 0, and pure
        Rutherford elastic, and its resolved-configuration line says "beam
        stops before the strip centre: elastic only". `talys-cache` reports
        "N beyond the beam's range".

      Real 37Cl at 30 MeV over strips 10–17 gives closed channels at 10–14,
      Rutherford-only at 12–14, and stopped strips at 15–17.
      > 2026-10-15 ewt: accepted.
- [x] **Strict step shape** (`src/ControlFile.cpp:222-243, 861-891,
      1001-1007`). Step 0 only, since TALYS needs the beam energy. The target
      is `4He` (TALYS projectile `a`), labelled as (a,a'): evap = target,
      res = beam. Exactly one `ground` + `rutherford` branch and one
      `talys` + `isotropic` branch, and no `weight` keys. Isotropic, because
      the TALYS weight is angle-integrated. Under audit F-12, an elastic-
      labelled step needs a weighted elastic branch, so this cannot be a
      step-level (inelastic-only) table.
      > 2026-10-15 ewt: accepted.
- [x] **TALYS kinematics** (`src/StripTables.cpp:82-89`): E_cm is the
      relativistic invariant at the strip centre's beam energy. TALYS gets
      E_α = E_cm (m_α + m_A)/m_A, which its non-relativistic conversion turns
      back into that E_cm.
      > 2026-10-15 ewt: accepted for now; I need to think about this more
      > later.
- [x] **Cache identity** (`src/StripTables.cpp:90-97, 140-145`). A table's
      file name ends in the MD5 of its exact TALYS input, which includes the
      Atomki-V2 potential's MD5. `musicsim` also refuses a table whose
      recorded input is not the strip's. TALYS's own build is not part of the
      key, as with SRIM tables, so delete `talys/` after changing TALYS.
      > 2026-10-15 ewt: accepted.
<!---->
## Found during the work
<!---->
- [x] **Pre-existing: `EnergyLoss::GetFinalEnergy` is one linear
      step** (`src/EnergyLoss.cpp:212-223`): `Eloss = dE/dx(E₀) × L`. Over a
      long path it ignores the beam slowing down. The first per-strip
      energies used it and were up to 4.2 MeV of beam too high at strip 17
      (37Cl, 92 MeV, 220 Torr). The fix integrates in `sim_step` slices
      (`src/StripTables.cpp:18-27`), and each strip's mean vertex energy now
      matches within ±0.07 MeV. The same single long step is still in
      log-only values: `Kb_max`/`Kb_min` in the per-strip log
      (`src/EventLoop.cpp:52-55`) and `CMEMin`, "Ecom(final)" in
      `energetics.log` (`src/Kinematics.cpp:405`). `energetics.log`'s
      per-strip table steps once per strip, which is closer. No event uses
      these. Check: fix them separately, or leave them.
      > 2026-10-15 ewt: fix. Done: those calls now go through
      > `MeanBeamEnergyAfter` (`src/StripTables.cpp:18-27`), and the time of
      > flight through the integrating `GetTimeOfFlight(medium, E, L)`.
      > - **`energetics.log`** (37Cl at 92 MeV): each strip centre now
      >   matches the integrated value to 0.004 MeV, where it was up to 0.29
      >   MeV high. "Ecom(final)" goes from 4.86 to 4.31 MeV.
      > - **`musicsim.log`, strip 15:** Kb was 59.6 → 57.7 MeV and is now
      >   56.5 → 54.2. tof was 0 → 0 ns (stale) and is now 13.9 → 14.8 ns.
      >
      > Events are unchanged: the suite and the old-draw reference controls
      > pass.
<!---->
## Code — review these lines
<!---->
- [x] `src/Kinematics.cpp:148-151, 323-362` — `SelectBranch` takes the strip's weight
      in a computed step, and `BranchTable` the strip's table. Still one draw
      for a multi-branch step and none for a single branch.
      > 2026-10-15 ewt: accepted.
- [x] `src/EventLoop.cpp:50-57` — the strip's z range now comes from
      `StripZRange` (the old loop, moved), and the strip being simulated is
      recorded for the kinematics.
      > 2026-10-15 ewt: accepted.
- [x] `src/StripTables.cpp:100-178` — tables are read once by the master
      (workers use its copies), checked against the strip's input, and stored
      in `levelTableText_`, so metadata and update checks cover them.
      > 2026-10-15 ewt: accepted.
- [x] `src/Output.cpp:417-500` — computed branches print `weight = computed
      per strip`, and each strip gets a `reaction_step_0_strip_N` line: beam
      energy, E_cm, weights, table path and MD5. An update run therefore
      refuses a changed table, weight, or strip energy.
      > 2026-10-15 ewt: accepted.
- [x] `src/talys-cache.cpp:37-214` — forks TALYS in `<table>.work`. The
      table is written to `.tmp` and renamed. The work directory is removed
      on success and kept on failure.
      > 2026-10-15 ewt: accepted, with a fix. Work and `.tmp` names were
      > fixed per table, so two concurrent `talys-cache` runs needing the same
      > table could clash. They now carry the PID, and the atomic rename
      > makes the last (identical) writer win. Checked by running four at
      > once on the test control: all succeed, 8 tables, no leftover
      > files, and `computed-weights` passes on the result.
- [x] `src/ControlFile.cpp:150-165, 205-255` — a branch `weight` is now
      optional while parsing; a step without a talys-cache branch still
      requires one, with the old message.
      > 2026-10-15 ewt: accepted.
<!---->
## Schema, build, and docs
<!---->
- [x] `README.md:280-380`: the "Per-strip TALYS tables" section (TALYS
      elastic, warnings, the beam-stops case), the metadata paragraph, and
      the limitations.
      > 2026-10-15 ewt: accepted.
- [x] `flake.nix:61-62, 78-84, 105-106` and the `Makefile`. `talys-cache`
      gets TALYS and the potentials compiled in. `musicsim` gets the
      potentials directory, since it picks the default model. The dev shell
      exports both.
      > 2026-10-15 ewt: accepted, keeping it in the default package (the
      > closure grows by TALYS's 8 GB database). The build should go to the
      > binary cache, and the README should say so.
- [x] `flake.nix:31`: version, was `26.10.14`.
      > 2026-10-15 ewt: bump. Done: `26.10.15`, the date convention of
      > 372f249. The build is served by e-desktop's binary cache
      > (`https://cache.ethanwtodd.com`), which the README's new "Binary
      > cache" section documents with its signing key. It serves whatever
      > e-desktop has built, so a revision is there once it is built there.
<!---->
## Tests
<!---->
- [x] `tests/talys_tables.cpp`:
      - Rutherford against a direct integral at 1, 10 and 120°;
      - keyword refusals and deck order;
      - the level parser on real TALYS 2.2 output (37Cl at E_cm 8.65 MeV,
        trimmed), behind three decoys: a synthetic 37Ar isobar, a 41K
        section, and a second parity block. Also covered: a closed channel,
        a total mismatch, a missing residue, an unfinished run, shifted bins;
      - `CheckElasticRows` (`:185`) on real TALYS elastic output: shape plus
        compound at 180°, refusal at another E_cm, the all-zero fallback, a
        short grid;
      - `CheckElasticTable` (`:256`): an empty table reproduces the old
        Rutherford draws exactly, the integral matches numerical integration
        to 1e-6, the sampled angle bands are within 5σ of the table, and
        sampling takes exactly one draw.
      > 2026-10-15 ewt: accepted.
- [x] `tests/run-tests.sh:192-243` with `tests/fake-talys.sh`:
      - `musicsim` refuses to run without tables; `talys-cache` builds 8 and
        then reuses 8.
      - `root_checks computed-weights` (`tests/root_checks.cpp:937-1076`)
        checks each strip:
        - the elastic weight is the stored table's elastic integral, which
          `ElasticIntegral` (`:902`) recomputes independently;
        - the inelastic weight is the table's sum, and the MD5 matches;
        - the inelastic draw fraction is within 5σ;
        - the mean vertex energy is within 0.15 MeV of the recorded centre
          energy;
        - pooled over strips, the elastic angles beyond 150° are within 5σ.
          The fake cuts elastic to a quarter there.
      - The closed strip 9 is reported, warned about, and kept in metadata.
      - `talys-cache-stops` (37Cl at 30 MeV, strips 13–16, two threads) has
        two closed strips, one Rutherford-only, and two beyond range. All of
        it is printed, recorded, and on the strip lines.
      - A table copied under another strip's name is refused.
      > 2026-10-15 ewt: accepted.
- [x] Eight invalid controls, each rejected with its own message
      (`tests/run-tests.sh:48-55`).
      > 2026-10-15 ewt: accepted.
- [x] Mutation-checked in a scratch copy. Each of these 18 fails at least
      one check:
      - weights from the wrong strip;
      - the single-step centre energy;
      - normalized instead of mb weights;
      - a dropped cot²;
      - the ground row kept;
      - a section matched without its Z;
      - no table-provenance check;
      - `talys-cache` rebuilding cached tables;
      - a closed table refused;
      - the elastic weight as pure Rutherford;
      - the elastic angles as pure Rutherford;
      - compound elastic dropped;
      - the frame check on total elastic instead of shape elastic;
      - any `*ang.L00` file accepted;
      - the cell ratio taken at its start;
      - no beam-stops warning;
      - warnings not recorded;
      - no all-zero fallback.

      Dropping only one of the two guards that stop the parser at the first
      population block survives, because the other still stops it.
      > 2026-10-15 ewt: accepted.
- [x] Real TALYS, by hand (not in the suite):
      - 37Cl at 92 MeV, θ_min 60°: strip 0's elastic weight is 345 mb, where
        pure Rutherford gave 787. The inelastic share rises from 1.5% to
        3.4%.
      - 87Rb at 330 MeV, θ_min 30° (`result/bin`, no environment variables):
        Atomki-V2 is picked by default. Elastic is 0.93 of Rutherford at
        strip 0 (E_cm 13.1), 0.997 at strip 6, and 1.000 at strip 12.
        Inelastic falls from 5.7 mb to 3e-23 mb.
      - 37Cl at 30 MeV, strips 10–17: the warnings listed above.
      > 2026-10-15 ewt: accepted.
