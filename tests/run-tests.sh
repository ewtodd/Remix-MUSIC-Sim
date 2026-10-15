#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
# Defaults target the `nix build` workflow: run `nix build` in the repo, then
# `bash tests/run-tests.sh` with the result/ symlink in place. Override with
# MUSICSIM_BIN / CONVERTER_BIN / *_CHECK_BIN for a make-in-dev-shell build.
simulator=${MUSICSIM_BIN:-"$repo_dir/result/bin/musicsim"}
converter=${CONVERTER_BIN:-"$repo_dir/result/bin/legacy-msc-to-toml"}
root_checks=${ROOT_CHECK_BIN:-"$repo_dir/result/bin/root_checks"}
vavilov_checks=${VAVILOV_CHECK_BIN:-"$repo_dir/result/bin/vavilov_statistics"}
nuclide_checks=${NUCLIDE_CHECK_BIN:-"$repo_dir/result/bin/nuclide_loader"}
excitation_checks=${EXCITATION_CHECK_BIN:-"$repo_dir/result/bin/excitation_table"}
temp_dir=$(mktemp -d "${TMPDIR:-/tmp}/remix-music-tests.XXXXXX")
trap 'test -n "${temp_dir:-}" && rm -rf -- "$temp_dir"' EXIT

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

expect_failure() {
  if "$@" >"$temp_dir/expected-failure.log" 2>&1; then
    fail "command unexpectedly succeeded: $*"
  fi
}

# Each invalid control must be rejected, and for the reason it was written for.
expect_check_error() {
  expect_failure "$simulator" --check "$repo_dir/tests/controls/$1.toml"
  grep -F -- "$2" "$temp_dir/expected-failure.log" >/dev/null ||
    fail "$1 was not rejected with: $2"
}

printf '%s\n' '[1/10] validating CLI and strict schema'
expect_failure "$simulator"
expect_failure "$simulator" --check "$temp_dir/missing.toml"
expect_check_error invalid-zero-step 'at least one complete reaction.step'
expect_check_error invalid-branch-ambig 'alongside [[reaction.step.branch]]'
expect_check_error invalid-zero-weight 'needs at least one weight > 0'
expect_check_error invalid-too-many-branches 'has more than 8 branches'
expect_check_error invalid-branch-forced \
  "must be 'ground', 'uniform', or 'talys'"
expect_check_error invalid-elastic-unweighted \
  'elastic step 0 requires a branch'
expect_check_error invalid-sim-step "'run.sim_step' must be finite and > 0"
expect_check_error invalid-atima14 "requires catima's GLOBAL library"
expect_check_error invalid-global-straggling "requires catima's GLOBAL library"
# The README's full example must stay a valid control.
awk '/^```toml/ { blocks++; if (blocks == 1) { inside = 1; next } }
     /^```/ { inside = 0 } inside' "$repo_dir/README.md" \
  >"$temp_dir/readme-example.toml"
"$simulator" --check "$temp_dir/readme-example.toml" \
  >"$temp_dir/check.log" 2>&1 || {
  sed -n '1,80p' "$temp_dir/check.log" >&2
  fail "README's full TOML example is not a valid control"
}
while IFS= read -r control; do
  "$simulator" --check "$control" >"$temp_dir/check.log" 2>&1 || {
    sed -n '1,160p' "$temp_dir/check.log" >&2
    fail "schema validation failed for $control"
  }
done < <(printf '%s\n' "$repo_dir/basic.toml"; \
         find "$repo_dir/ControlExamples" -type f -name '*.toml' | sort; \
         find "$repo_dir/tests/controls" -type f -name '*.toml' ! -name 'invalid-*' | sort)

printf '%s\n' '[2/10] testing legacy converter safety and escaping'
expect_failure "$converter" "$temp_dir/no-match-*.msc"
expect_failure "$converter" --keep --remove-source \
  "$repo_dir/tests/fixtures/legacy.msc"
"$converter" --stdout "$repo_dir/tests/fixtures/legacy.msc" \
  >"$temp_dir/legacy-stdout.toml"
"$simulator" --check "$temp_dir/legacy-stdout.toml" >"$temp_dir/check.log" 2>&1
grep -F 'output' "$temp_dir/legacy-stdout.toml" | grep -F '\\' | grep -F '\"' >/dev/null
cp "$repo_dir/tests/fixtures/legacy.msc" "$temp_dir/legacy.msc"
"$converter" "$temp_dir/legacy.msc" >/dev/null
test -f "$temp_dir/legacy.toml" || fail 'converter did not create output'
test -f "$temp_dir/legacy.msc" || fail 'converter removed source by default'
expect_failure "$converter" "$temp_dir/legacy.msc"
"$converter" --force "$temp_dir/legacy.msc" >/dev/null
cp "$repo_dir/tests/fixtures/legacy.msc" "$temp_dir/remove.msc"
"$converter" --remove-source "$temp_dir/remove.msc" >/dev/null
test ! -e "$temp_dir/remove.msc" || fail '--remove-source kept its input'
# Every .msc predates the June 2026 percent-FWHM eres, so Eres is an absolute
# sigma in MeV. A non-elastic final step is inferred as uniform excitation.
grep -Fx 'noise_sigma_mev = 0.050000000000000003 # legacy absolute sigma [MeV]' \
  "$temp_dir/legacy-stdout.toml" >/dev/null ||
  fail 'converter did not map legacy Eres to an absolute MeV sigma'
"$converter" --stdout "$repo_dir/tests/fixtures/legacy-reaction.msc" \
  >"$temp_dir/legacy-reaction.toml"
"$simulator" --check "$temp_dir/legacy-reaction.toml" \
  >"$temp_dir/check.log" 2>&1 ||
  fail 'converted reaction control failed schema validation'
grep -F 'residue_excitation = "uniform"' "$temp_dir/legacy-reaction.toml" \
  >/dev/null || fail 'converter did not infer uniform for a final (a,p) step'

printf '%s\n' '[3/10] testing nuclide loaders (ERROR lines below are expected: malformed-input assertions)'
"$nuclide_checks" "$temp_dir"

printf '%s\n' '[4/10] validating Vavilov interpolation and limits'
"$vavilov_checks"

printf '%s\n' '[5/10] testing Talys level-table loading and sampling'
"$excitation_checks" "$repo_dir/tests/fixtures/talys-levels.tab" "$temp_dir"

run_control() {
  local control=$1
  (cd "$temp_dir" && "$simulator" "$repo_dir/tests/controls/$control.toml" \
    >"$temp_dir/$control.log" 2>&1) || {
    sed -n '1,200p' "$temp_dir/$control.log" >&2
    fail "simulation failed for $control"
  }
}

printf '%s\n' '[6/10] testing stopping-energy conservation, fast-ion stopping, and forbidden reactions'
run_control stopping
"$root_checks" standard "$temp_dir/stopping.root" 1
"$root_checks" conservation "$temp_dir/stopping.root"
# A 16 MeV/u alpha is above catima's 10 MeV/u switch to Bethe stopping, where
# the charge model matters. The default must be the explicit winger result.
run_control fast-alpha-default
run_control fast-alpha-winger
"$root_checks" standard "$temp_dir/fast-alpha-default.root" 1
"$root_checks" compare "$temp_dir/fast-alpha-default.root" \
  "$temp_dir/fast-alpha-winger.root"
run_control forbidden-reaction
"$root_checks" standard "$temp_dir/forbidden.root" 1
"$root_checks" forbidden "$temp_dir/forbidden.root"
run_control later-step-forbidden
"$root_checks" standard "$temp_dir/later-step-forbidden.root" 1
"$root_checks" forbidden "$temp_dir/later-step-forbidden.root"
# Only the log proves step 0 succeeded and step 1 failed; the fixture sets
# print_opt = 1 so the kinematics write it.
grep -F 'Reaction step 1 is below threshold.' "$temp_dir/musicsim.log" \
  >/dev/null || fail 'later-step fixture did not reject step 1 after step 0'

printf '%s\n' '[7/10] testing multi-step truth angles, boundaries, and stopped or side-exiting particles'
run_control multistep
"$root_checks" standard "$temp_dir/multistep.root" 1
"$root_checks" multistep "$temp_dir/multistep.root"
"$root_checks" reconstruct "$temp_dir/multistep.root"
run_control boundary-coarse
run_control boundary-fine
"$root_checks" boundary "$temp_dir/boundary-coarse.root" \
  "$temp_dir/boundary-fine.root"
run_control stopped-before-vertex
"$root_checks" standard "$temp_dir/stopped-before-vertex.root" 1
"$root_checks" stopped "$temp_dir/stopped-before-vertex.root" 15
run_control side-exit
"$root_checks" standard "$temp_dir/side-exit.root" 1
"$root_checks" side-exit "$temp_dir/side-exit.root"
run_control two-body-ground
"$root_checks" standard "$temp_dir/two-body-ground.root" 200
"$root_checks" two-body "$temp_dir/two-body-ground.root"

printf '%s\n' '[8/10] testing weighted branches and Talys excitation'
run_control branches
"$root_checks" standard "$temp_dir/branches.root" 800
"$root_checks" branches "$temp_dir/branches.root" 0.25
"$root_checks" level-table "$temp_dir/branches.root" \
  "$repo_dir/tests/controls/../fixtures/talys-levels.tab" \
  "$repo_dir/tests/fixtures/talys-levels.tab"
# Effective weight is weight x open fraction: 1 x 0.5 against 3, then 0.
run_control branches-partial
"$root_checks" branches "$temp_dir/branches-partial.root" 0.142857142857
run_control branches-closed
"$root_checks" branches "$temp_dir/branches-closed.root" 0
# Step-level models keep the reaction draws of fea3f4b, the last build before
# weighted branches: their CM angles must match that build's reference.
for control in legacy-draws-chain legacy-draws-elastic; do
  run_control "$control"
  "$root_checks" cm-angles "$temp_dir/$control.root" \
    >"$temp_dir/$control.angles"
  grep -v '^#' "$repo_dir/tests/fixtures/$control.angles" |
    diff - "$temp_dir/$control.angles" >&2 ||
    fail "$control no longer reproduces the pre-branch reaction draws"
done
run_control talys-step
"$root_checks" standard "$temp_dir/talys-step.root" 40
"$root_checks" talys "$temp_dir/talys-step.root"

printf '%s\n' '[9/10] testing thread-count-independent events and merge output'
run_control repro-single
run_control repro-multi
"$root_checks" standard "$temp_dir/repro-single.root" 4
"$root_checks" standard "$temp_dir/repro-multi.root" 4
"$root_checks" compare "$temp_dir/repro-single.root" \
  "$temp_dir/repro-multi.root"
# A pipe can be read only once, so this fails if workers or the metadata
# writer re-read the control file instead of using the text loaded at start.
(cd "$temp_dir" && "$simulator" \
  <(cat "$repo_dir/tests/controls/control-read-once.toml") \
  >"$temp_dir/control-read-once.log" 2>&1) || {
  sed -n '1,200p' "$temp_dir/control-read-once.log" >&2
  fail 'run from a read-once control pipe failed'
}
"$root_checks" compare "$temp_dir/repro-multi.root" \
  "$temp_dir/control-read-once.root"
"$root_checks" control-text "$temp_dir/control-read-once.root" \
  "$repo_dir/tests/controls/control-read-once.toml"
# Level tables are read once too: a worker re-reading this FIFO would block
# until the timeout, since nothing writes to it a second time.
mkfifo "$temp_dir/levels.fifo"
sed -e "s|^level_table .*|level_table          = \"$temp_dir/levels.fifo\"|" \
  -e 's|^events_per_strip = 100|events_per_strip = 5|' \
  -e 's|"branches.root"|"table-read-once.root"|' \
  "$repo_dir/tests/controls/branches.toml" >"$temp_dir/table-read-once.toml"
cat "$repo_dir/tests/fixtures/talys-levels.tab" >"$temp_dir/levels.fifo" &
writer=$!
(cd "$temp_dir" && timeout 120 "$simulator" "$temp_dir/table-read-once.toml" \
  >"$temp_dir/table-read-once.log" 2>&1) || {
  kill "$writer" 2>/dev/null || true
  sed -n '1,200p' "$temp_dir/table-read-once.log" >&2
  fail 'run from a read-once level-table pipe failed'
}
wait "$writer"
"$root_checks" level-table "$temp_dir/table-read-once.root" \
  "$temp_dir/levels.fifo" "$repo_dir/tests/fixtures/talys-levels.tab"
# A failed threaded run removes the ROOT files it made, so the next run is not
# blocked by its .tmp, and keeps files it did not make. Worker 1 meets a stale
# temporary file; the publish rename meets a directory.
touch "$temp_dir/mt-worker-fail_w1.root.tmp"
mkdir "$temp_dir/mt-publish-fail.root"
for control in mt-worker-fail mt-publish-fail; do
  (cd "$temp_dir" && expect_failure "$simulator" \
    "$repo_dir/tests/controls/$control.toml")
  for leftover in "${control}_w0.root" "${control}_w1.root" \
    "$control.root.tmp"; do
    test ! -e "$temp_dir/$leftover" || fail "failed $control run left $leftover"
  done
done
test -e "$temp_dir/mt-worker-fail_w1.root.tmp" ||
  fail 'threaded-run cleanup removed a temporary file it did not create'
test -d "$temp_dir/mt-publish-fail.root" ||
  fail 'threaded-run cleanup removed the blocking directory'

printf '%s\n' '[10/10] testing checked, atomic update output'
run_control update
run_control update
"$root_checks" standard "$temp_dir/update.root" 4
"$root_checks" unique "$temp_dir/update.root"
# Same control text, but the file claims another build's default: refuse.
"$root_checks" alter-resolved "$temp_dir/update.root"
(cd "$temp_dir" && expect_failure "$simulator" \
  "$repo_dir/tests/controls/update.toml")
grep -F 'resolved defaults, or code revision' "$temp_dir/expected-failure.log" \
  >/dev/null || fail 'update appended to a file from another configuration'
# Same control text, but its level table changed after the file was written.
cp "$repo_dir/tests/fixtures/talys-levels.tab" "$temp_dir/update-levels.tab"
sed -e 's|^level_table .*|level_table          = "update-levels.tab"|' \
  -e 's|^events_per_strip = 100|events_per_strip = 2|' \
  -e 's|^threads          = 2|threads          = 1|' \
  -e 's|^file_opt         = "recreate"|file_opt         = "update"|' \
  -e 's|"branches.root"|"table-update.root"|' \
  "$repo_dir/tests/controls/branches.toml" >"$temp_dir/table-update.toml"
(cd "$temp_dir" && "$simulator" "$temp_dir/table-update.toml" \
  >"$temp_dir/table-update.log" 2>&1) || {
  sed -n '1,200p' "$temp_dir/table-update.log" >&2
  fail 'level-table update run failed'
}
sed 's|^level 1.727 0.40|level 1.727 0.50|' \
  "$repo_dir/tests/fixtures/talys-levels.tab" >"$temp_dir/update-levels.tab"
grep -Fx 'level 1.727 0.50' "$temp_dir/update-levels.tab" >/dev/null ||
  fail 'could not edit the copied level table'
(cd "$temp_dir" && expect_failure "$simulator" "$temp_dir/table-update.toml")
grep -F 'level tables, resolved defaults' "$temp_dir/expected-failure.log" \
  >/dev/null || fail 'update appended to a file made with another level table'

printf '%s\n' 'All tests passed.'
