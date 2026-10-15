# Review checklist — weighted branches + Talys level-table excitation, 2026-10-14
<!---->
Per the AI-disclosure section of the README, changes to simulation physics
(excitation sampling, kinematics, branch weighting) must be human-reviewed.
This covers the **uncommitted** implementation of
`PLAN-talys-weighted-branch-excitation-2026-10-14.md` (sim side only; the
MUSIC-repo Talys producer, plan §6, is not started). Line numbers are for the
working tree on 2026-10-14. No code is reproduced here; open the file at the
line and review. Checks indicate review by ewt.
<!---->
Written by Claude Code (Opus 5.5). Verified so far:
`nix build && bash tests/run-tests.sh` passes all 10 stages, and the four
reacted legacy controls (elastic, uniform, multistep, later-step-forbidden)
reproduce HEAD's events bit for bit on every shared leaf.
<!---->
## Merge onto fea3f4b (2026-10-14, before review)
<!---->
The work was rebased from eecdc8a onto the reviewed audit commit. Line numbers
below predate the merge. Changes made while merging:
<!---->
- Output schema is **5**, not 4: fea3f4b already used 4 for the angle
  sentinel, `beam_energy_accel`, and the expanded resolved configuration.
- The elastic rule keeps fea3f4b's labelled order (beam-like product is
  `res`; swapped labels are rejected) instead of "either order". The
  weighted-elastic-branch rule applies to the labelled form.
- `rutherford` after step 0 is now also rejected inside branches.
- The per-branch step description moved into fea3f4b's `RunConfiguration()`,
  so the update-run configuration check covers branches and weights.
- `level_table` resolution uses the read-once `sourcePath`.
- README: the elastic paragraph documents the label order and the step-0
  Rutherford rule.
- 10/10 stages pass; `format-check` and `git diff --check` are clean.
  `make strict` only looked clean: the dev shell puts `$PWD/include` on
  `CPLUS_INCLUDE_PATH`, which GCC treats as system headers, so it hid the
  `-Wdeprecated-copy` note below. That note applied to fea3f4b too; it is
  fixed separately (Makefile `unexport`, defaulted `FourVector` copy).
<!---->
## Decisions the plan left open or got wrong — confirm or overturn
<!---->
- [x] **Branch tables replace step-level models outright.** The plan's §1 TOML
      example puts step-level `ground`/`rutherford` + `weight` next to a
      `[[reaction.step.branch]]`, but its "implement exactly" rules (and Task
      1.2's ambiguity test) forbid that. I followed the rules
      (`src/ControlFile.cpp:180-184`) and made a step-level `weight` an error
      instead of a silently ignored key (`src/ControlFile.cpp:540-543`).
      Check: the two-branch form in `README.md:228-246` is the schema you want.
      > 2026-10-14 ewt: accepted.
- [x] **`forced` is refused inside a branch** (`src/ControlFile.cpp:109-110,
      127-130`). `minEx` (`src/Simulator.cpp:407`) assumes one model per step,
      and a mixed forced/non-forced step would make it ambiguous.
      Check: no planned control needs a forced branch.
      > 2026-10-14 ewt: accepted.
- [x] **`level_table` resolves relative to the control file**
      (`src/ControlFile.cpp:82-107, 230-231`), as the plan recommends. Note that
      SRIM tables resolve relative to the *output* directory
      (`src/ControlFile.cpp:843`), so the two conventions differ.
      > 2026-10-14 ewt: accepted. Level tables are inputs that travel with the control; SRIM tables are generated caches next to the output.
- [x] **Over-`EneAvail` tail: truncate and renormalize, one draw**
      (`src/ExcitationTable.cpp:127-170`), instead of the plan's "re-draw up to
      N, then clamp". Clamping would pile probability at `E_ex = EneAvail`, a
      zero-energy ejectile. A bin straddling `EneAvail` counts pro rata.
      Check: conditioning on `E_ex <= EneAvail` is the physics you want.
      > 2026-10-14 ewt: accepted within the branch, but the closed part now also lowers the branch's draw weight (see the next item), so a partly closed table no longer keeps its full share.
- [x] **A closed inelastic channel rejects the reaction.** If the talys
      branch is drawn and no tabulated weight lies below `EneAvail`, the event
      goes unreacted like a below-threshold step (`src/Kinematics.cpp:128-131`,
      `src/ExcitationTable.cpp:139-142`, per plan §3 step 4). It does not fall
      back to the elastic branch. Check: acceptable, or renormalize over the
      open branches instead.
      > 2026-10-14 ewt: overturned. `SelectBranch` draws by weight × open fraction at the vertex (`ExcitationTable::AllowedWeight` for talys, 1 otherwise), still one draw; a closed branch is never drawn, and the reaction is rejected only when every branch is closed. New runs branches-partial (expect 1/7) and branches-closed (expect 0, no rejected events); both fail under the fixed-weight rule. README updated.
- [x] **Two new MC truth branches**, `residue_excitation[n_steps]` and
      `reaction_branch[n_steps]` (`src/Output.cpp:239-242, 300-301`). The plan
      did not list them, but its tests and manual validation need per-event
      `E_ex`/branch truth, which did not exist. Consequence: `file_opt =
      "update"` now refuses a schema-3 file as an incompatible tree schema
      (verified; the file is left untouched).
      > 2026-10-14 ewt: accepted; schema 5 after the merge. Their rejected-chain sentinels (-2, -1) are now checked on the forbidden and stopped runs.
- [x] **Level tables are read once and recorded** (found during review).
      Each worker re-read its tables from disk and the metadata kept only the
      path, so an edited table could reach part of a run, and an update run
      could append to a file made with another table.
      > 2026-10-14 ewt: fix accepted. The master reads each table once (`levelTableText_`) and workers parse its copies; metadata stores `level_table_N` and `level_table_N_path`; each talys branch line in the resolved configuration carries the table's MD5, so an update run refuses a changed table. Tests: stage 8 metadata check, stage 9 FIFO read-once run, stage 10 edited-table update refusal; each fails against a build with the old behaviour.
<!---->
## Core physics — review these lines
<!---->
- [x] `src/Kinematics.cpp:124-135` — the branch is drawn before any
      excitation or angle draw, and the drawn branch's model feeds excitation
      and the angular code (`:142-146`).
      > 2026-10-14 ewt: accepted. Now `src/Kinematics.cpp:131-145`; a step with every branch closed is rejected there (D5).
- [x] `src/Kinematics.cpp:313-335` — `SelectBranch`: a one-branch step returns
      before drawing (`:318-319`), so legacy controls add **zero** RNG draws.
      A multi-branch step costs exactly one draw; zero-weight branches can
      never be drawn. Check: this is the only extra draw. Note that nothing in
      the suite pins this permanently: the repro tests are unreacted beam runs.
      > 2026-10-14 ewt: accepted, with weight x open fraction (D5); a closed branch is never drawn either. Now pinned: legacy-draws-chain (forced, uniform, isotropic) and legacy-draws-elastic (ground, Rutherford) must reproduce CM angles from a fea3f4b build (tests/fixtures/*.angles, stage 8). Both fail when a one-branch step burns a draw; all 2608 leaf values both builds share were bit-identical when the reference was made.
- [x] `src/Kinematics.cpp:337-365` — `SampleExcitation`: cases 0/1/2 move the
      old switch verbatim (same draws in the same order); case 3 calls the
      table (`:358-360`).
      > 2026-10-14 ewt: accepted. Cases 0/1/2 match fea3f4b's switch draw for draw.
- [x] `src/Kinematics.cpp:243-244` — a rejected chain also clears the new
      truth fields.
      > 2026-10-14 ewt: accepted; the sentinels are now checked on the forbidden and stopped runs.
- [x] `src/ExcitationTable.cpp:137-170` — inverse CDF over point masses, then
      bins, restricted to `E_ex <= EneAvail`. Ground is always allowed
      (`EneAvail` is clamped `>= 0`, and `ground` rows must sit at 0, `:80-81`).
      Row order does not matter.
      > 2026-10-14 ewt: accepted.
- [x] `src/ExcitationTable.cpp:37-125` — loader: unknown tag, `hi <= lo`,
      negative values, non-numeric or extra fields, and no positive weight are
      hard errors with `file:line`; weights normalize to 1.
      > 2026-10-14 ewt: accepted. The loader is now Parse(text, source) plus Load(path), so the master's single read can be shared.
- [x] `src/ControlFile.cpp:809-829` — **REVIEW-GATED relaxation of audit
      F-12.** Without branches, the strict rule is unchanged. With branches,
      it needs at least one `ground` + `rutherford` branch with **weight > 0**,
      so a zero-weighted elastic branch cannot hide an all-inelastic channel.
      Both labelling orders are covered (`:804-810`, matching main's split).
      > 2026-10-14 ewt: accepted. After the merge only the labelled order is elastic; swapped labels are rejected.
- [x] `src/ControlFile.cpp:197-208` — weights are user-authored relative cross
      sections, normalized at load. Before production, confirm each channel's
      values against Talys (inelastic σ) and the Rutherford integral above
      `theta_cm_min_deg` (elastic σ). No `ControlExamples/` file was migrated
      to branches.
      > 2026-10-14 ewt: acknowledged; check every branch weight against Talys and the Rutherford integral before production.
<!---->
## Schema, output, and docs
<!---->
- [x] `src/ControlFile.cpp:168-209, 605-611` — branch parsing: required keys,
      the 8-branch cap, a per-branch `theta_cm_min_deg` defaulting to the
      step's value.
      > 2026-10-14 ewt: accepted. Now `src/ControlFile.cpp:125-223`.
- [x] `src/ControlFile.cpp:716-755` — required-key checks skip branched steps.
      Every other step gets one implicit branch, and level tables load at
      parse time, so `--check` catches a bad table.
      > 2026-10-14 ewt: accepted without a `--check` test for a malformed table; the unit tests cover the loader directly.
- [x] `src/Output.cpp:343, 361-384` — schema version 4; a branched step
      records one `reaction_step_N_branch_M` metadata line per branch.
      > 2026-10-14 ewt: accepted; schema 5 after the merge, and talys branch lines now carry the table MD5.
- [x] `README.md:194-265` (schema, table format, elastic rule), `:375` (truth
      branches), `:393-396` (schema 4, update refusal), `:409-413` (no longer
      claims continuum-only; states the table's limits).
      > 2026-10-14 ewt: accepted, including the merge and review updates (open-fraction rule, stored tables, update refusal).
- [ ] `flake.nix:27` — version `26.10.14`, bumped the same way as `372f249`.
<!---->
## Tests
<!---->
- [x] `tests/excitation_table.cpp:74-196` — shape, truncation, a level exactly
      at `EneAvail`, no draw on rejection, exactly one draw per sample,
      unsorted rows, malformed tables. Mutation-checked: the test catches
      unclipped bins, a two-draw sampler, and dropping the per-level
      `EneAvail` guard.
      > 2026-10-14 ewt: accepted. Re-verified: all three mutants fail with their own messages.
- [x] `tests/root_checks.cpp:209-292` — `branches`: drawn fraction within 5σ
      of 0.25, Q = 0 energy balance to 1e-3 MeV, angular model per branch;
      `talys`: a step-level table on every event. Mutation-checked:
      weight-blind selection and angles read from the legacy arrays both fail.
      > 2026-10-14 ewt: accepted; the D5 partial and closed runs extend it.
- [x] `tests/run-tests.sh:29-45` — each invalid control must fail **with its
      own message**, not just fail; `:79-80, 111-117` are the new stages.
      > 2026-10-14 ewt: accepted.
- [x] Pre-existing, untouched: `make strict` (and so the sanitizer build)
      stops on `-Wdeprecated-copy` in `src/FourVector.cpp`. The changed files
      pass the strict flags apart from that inherited warning.
      > 2026-10-14 ewt: fixed by the separate strict-masking change (Makefile unexport, defaulted FourVector copy).
