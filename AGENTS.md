# Remix-MUSIC-Sim — notes for agents

README.md is the reference for building, the strict TOML control schema,
per-strip TALYS tables, and the test suite. Read the section you need.

## Physics changes need a review checklist before commit

- Simulation physics (kinematics, excitation and branch sampling, the
  energy-loss/straggling chain, beam energies, stopping models) is signed off
  by ewt before it is committed (README, "AI-assisted development disclosure").
- For such a change, write `checks/REVIEW-<slug>-<YYYY-MM-DD>.md` against the
  uncommitted working tree, copying the structure of the newest file in
  `checks/`: what it covers and that line numbers are the working tree's on
  that date; "Written by <tool> (<model>). Verified so far: …"; decisions and
  physics choices to confirm or overturn; then code, schema/build/docs, tests.
- Items are `- [ ] path:line-range — what the lines do`. Add a `Check:` only
  for a specific claim ewt must confirm (an invariant, value, ordering), not a
  restatement. Line numbers come from the working tree (`grep -n`), never a
  diff; within a file, ascending order. No code in the checklist.
- Leave boxes unchecked. ewt ticks them and replies inline
  (`> YYYY-MM-DD ewt: …`); the change and its checklist are committed together.
- Choosing an excitation model, angular model, or branch weight for a reaction
  channel is a physics decision: surface each one, never pick across a batch.

## Controls

- `legacy-msc-to-toml` converts only old whitespace `.msc` controls. Nothing
  converts an older TOML schema to the current one; migrate by hand and
  validate with `musicsim --check <control.toml>`.
- ~/MUSIC (ewtodd/MUSIC) keeps its sim controls in
  `analysis/<dataset>/sim_control/*.toml`. They are gitignored there, so a
  schema change here never reaches them. MUSIC's `RemixSim` reads only each
  control's `output = ".../traces_<iso>_<tag>.root"` line and loads the result
  from `analysis/<dataset>/sim_root_files/`.

## Commits

- A sentence-case imperative summary line with no prefix ("Add weighted
  reaction branches and Talys level-table excitation"), a blank line, then a
  body that explains what changed and why.
- Agent commits: `--author="son-of-anton-bot <307402699+son-of-anton-bot@users.noreply.github.com>"`.
