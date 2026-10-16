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
talys_checks=${TALYS_CHECK_BIN:-"$repo_dir/result/bin/talys_tables"}
talys_cache=${TALYSCACHE_BIN:-"$repo_dir/result/bin/talys-cache"}
srim_cache=${SRIMCACHE_BIN:-"$repo_dir/result/bin/srim-cache"}
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
expect_check_error invalid-talys-steplevel 'use reaction.talys_class'
expect_check_error invalid-talys-target "needs target.species = '4He'"
expect_check_error invalid-talys-block \
  '[talys] is only used with reaction.talys_class'
expect_check_error invalid-class-steps \
  'cannot be combined with reaction.talys_class'
expect_check_error invalid-class-element 'is not an element at or below'
expect_check_error invalid-class-theta 'elastic is only in the beam'
expect_check_error invalid-class-fraction "min_channel_fraction' must be a number"
expect_check_error invalid-class-update 'does not support run.update = 1'
expect_check_error invalid-class-key 'is only used with reaction.talys_class'
expect_check_error invalid-talys-keyword "sets 'energy', which talys-cache"
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

printf '%s\n' '[5/10] testing level tables, TALYS output parsing, and the Rutherford integral'
"$excitation_checks" "$repo_dir/tests/fixtures/talys-levels.tab" "$temp_dir"
# The fake TALYS's cascade for 37Cl at E_cm = 13 MeV, where every channel of
# its toy model is open.
mkdir "$temp_dir/fake-cascade"
(cd "$temp_dir/fake-cascade" &&
  printf '%s\n' 'projectile a' 'element 17' 'mass 37' 'energy 14.405405' |
  "$repo_dir/tests/fake-talys.sh")
"$talys_checks" "$repo_dir/tests/fixtures" "$temp_dir/fake-cascade"

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
# weighted branches, run with the current event seeding: their CM angles must
# match that build's reference.
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
# Per-strip TALYS tables: musicsim will not run without them, talys-cache
# builds each once (with a fake TALYS here), and the aa class's elastic and
# (a,a') weights follow strip by strip.
(cd "$temp_dir" && expect_failure "$simulator" \
  "$repo_dir/tests/controls/talys-cache.toml")
grep -F 'run: talys-cache' "$temp_dir/expected-failure.log" >/dev/null ||
  fail 'musicsim ran without its per-strip TALYS tables'
build_talys_tables() {
  (cd "$temp_dir" && TALYS_BIN="$repo_dir/tests/fake-talys.sh" \
    SRIM_TABLE_BIN="$repo_dir/tests/fake-srim.sh" \
    "$talys_cache" "$repo_dir/tests/controls/$1.toml" \
    >"$temp_dir/$1-build.log" 2>&1) || {
    sed -n '1,80p' "$temp_dir/$1-build.log" >&2
    fail "talys-cache failed for $1"
  }
  grep -F "$2" "$temp_dir/$1-build.log" >/dev/null ||
    fail "talys-cache did not report for $1: $2"
}
build_talys_tables talys-cache '8 generated, 0 already cached, 0 failed'
build_talys_tables talys-cache '0 generated, 8 already cached, 0 failed'
run_control talys-cache
"$root_checks" standard "$temp_dir/talys-cache.root" 3200
"$root_checks" computed-weights "$temp_dir/talys-cache.root" 34 120
# The fake closes (a,a') below E_cm = 7 MeV, in strip 9: musicsim warns and
# records it.
grep -F 'musicsim WARNING: strip 9 (E_cm' "$temp_dir/talys-cache.log" |
  grep -F 'no channel ending in Cl open' >/dev/null ||
  fail 'musicsim did not warn about the closed strip'
"$root_checks" metadata-contains "$temp_dir/talys-cache.root" warnings \
  'no channel ending in Cl open'
# A table under another strip's name was made from another TALYS input.
tables=("$temp_dir"/talys/*.tab)
cp "${tables[0]}" "${tables[1]}"
(cd "$temp_dir" && expect_failure "$simulator" \
  "$repo_dir/tests/controls/talys-cache.toml")
grep -F 'was not made from strip' "$temp_dir/expected-failure.log" \
  >/dev/null || fail "musicsim used a table made from another strip's input"
# A beam that stops in strip 14: strips past it get no table, and both those
# and the closed strips are reported, printed, and recorded. Beam that
# straggles past the mean range still reacts there, elastically.
build_talys_tables talys-cache-stops \
  "2 generated, 0 already cached, 0 failed, 2 beyond the beam's range"
run_control talys-cache-stops
"$root_checks" standard "$temp_dir/talys-cache-stops.root" 800
"$root_checks" class-channels "$temp_dir/talys-cache-stops.root" unreacted-ok
for warning in 'strip 15: the mean beam stops before the strip centre' \
  'strip 16: the mean beam stops before the strip centre' \
  'strip 13 (E_cm' 'gave no elastic angular distribution'; do
  grep -F 'musicsim WARNING:' "$temp_dir/talys-cache-stops.log" |
    grep -F -- "$warning" >/dev/null || fail "musicsim did not warn: $warning"
  "$root_checks" metadata-contains "$temp_dir/talys-cache-stops.root" \
    warnings "$warning"
done
"$root_checks" metadata-contains "$temp_dir/talys-cache-stops.root" \
  resolved_configuration 'beam stops before the strip centre: no TALYS table'
# The an class: (a,g), (a,gn), (a,n), and (a,2n), drawn per strip by cross
# section; (a,g) only upstream of each strip centre. (a,g) and (a,2n) dip below
# the 1.2% threshold at some strips, but a channel passing at any is kept.
# (a,2n) must leave 40K above its neutron threshold, which cuts the fake's
# bins at each strip and leaves none at strip 5.
build_talys_tables talys-class '6 generated, 0 already cached, 0 failed'
run_control talys-class
"$root_checks" standard "$temp_dir/talys-class.root" 12000
"$root_checks" class-channels "$temp_dir/talys-class.root"
for channel in '0 = 37Cl(a,g)41K' '3 = 37Cl(a,2n)39K'; do
  "$root_checks" metadata-contains "$temp_dir/talys-class.root" channels \
    "$channel"
done
grep -F 'musicsim WARNING: strip 5 (E_cm' "$temp_dir/talys-class.log" |
  grep -F "channel 'nn' has no state from which its chain can finish" \
    >/dev/null || fail 'musicsim did not warn that (a,2n) cannot finish'
# A 0.6 mb minimum drops (a,g), 0.5 mb, though it passes the 1.2% threshold.
sed -e '/^min_channel_fraction/a min_channel_mb = 0.6' \
  -e 's/^events_per_strip .*/events_per_strip = 20/' \
  -e 's/"talys-class.root"/"talys-class-mb.root"/' \
  "$repo_dir/tests/controls/talys-class.toml" >"$temp_dir/talys-class-mb.toml"
(cd "$temp_dir" && "$simulator" "$temp_dir/talys-class-mb.toml" \
  >"$temp_dir/talys-class-mb.log" 2>&1) ||
  fail 'simulation failed for talys-class-mb'
"$root_checks" metadata-contains "$temp_dir/talys-class-mb.root" channels \
  '0 = 37Cl(a,n)40K'
"$root_checks" metadata-contains "$temp_dir/talys-class-mb.root" \
  resolved_configuration 'dropped_mb = 0.5;'
if grep -F 'left out by the thresholds' "$temp_dir/talys-class-mb.log" \
  >/dev/null; then
  fail 'musicsim warned about left-out channels below min_channel_mb'
fi
# At 6.5%, (a,g) and (a,gn) each fail at every strip, but together they pass
# at strips 4 and 5, which is warned about.
sed -e 's/^min_channel_fraction .*/min_channel_fraction = 0.065/' \
  -e 's/^events_per_strip .*/events_per_strip = 20/' \
  -e 's/"talys-class.root"/"talys-class-sum.root"/' \
  "$repo_dir/tests/controls/talys-class.toml" >"$temp_dir/talys-class-sum.toml"
(cd "$temp_dir" && "$simulator" "$temp_dir/talys-class-sum.toml" \
  >"$temp_dir/talys-class-sum.log" 2>&1) ||
  fail 'simulation failed for talys-class-sum'
grep -F 'musicsim WARNING: strip 5 (E_cm' "$temp_dir/talys-class-sum.log" |
  grep -F 'left out by the thresholds add up to 2.5 mb' >/dev/null ||
  fail 'musicsim did not warn about channels left out together'
# TALYS exits with status 0 after its own errors; talys-cache must notice.
sed -e 's/^keywords .*/keywords = ["bins 40", "fakeerror y"]/' \
  "$repo_dir/tests/controls/talys-class.toml" >"$temp_dir/talys-error.toml"
(cd "$temp_dir" && export TALYS_BIN="$repo_dir/tests/fake-talys.sh" &&
  expect_failure "$talys_cache" "$temp_dir/talys-error.toml")
grep -F 'TALYS-error: fake failure' "$temp_dir/expected-failure.log" \
  >/dev/null || fail 'talys-cache used the output of a failed TALYS run'
# With SRIM stopping, talys-cache also makes the stopping tables of the
# class's products, which only its TALYS tables name.
(cd "$temp_dir" && SRIM_TABLE_BIN="$repo_dir/tests/fake-srim.sh" \
  "$srim_cache" "$repo_dir/tests/controls/talys-class-srim.toml" \
  >"$temp_dir/srim-cache.log" 2>&1) || fail 'srim-cache failed'
grep -F '2 generated, 0 already cached, 0 failed' "$temp_dir/srim-cache.log" \
  >/dev/null || fail 'srim-cache did not make the beam and target tables'
build_talys_tables talys-class-srim \
  "products: 3 generated, 0 already cached, 0 failed"
build_talys_tables talys-class-srim \
  "products: 0 generated, 3 already cached, 0 failed"
run_control talys-class-srim
# Its 2% threshold drops (a,g), 1% of the class, and records its 0.5 mb.
"$root_checks" metadata-contains "$temp_dir/talys-class-srim.root" channels \
  '0 = 37Cl(a,n)40K'
"$root_checks" metadata-contains "$temp_dir/talys-class-srim.root" \
  resolved_configuration 'dropped_mb = 0.5;'
for ion in 39K 40K 41K; do
  grep -F "SRIM table for $ion:" "$temp_dir/talys-class-srim.log" \
    >/dev/null || fail "musicsim did not read the $ion stopping table"
done

printf '%s\n' '[9/10] testing thread-count-independent events and merge output'
run_control repro-single
run_control repro-multi
"$root_checks" standard "$temp_dir/repro-single.root" 4
"$root_checks" standard "$temp_dir/repro-multi.root" 4
"$root_checks" compare "$temp_dir/repro-single.root" \
  "$temp_dir/repro-multi.root"
# A neighbouring master seed gives new events, not the same ones a few
# indices along.
sed -e "s|^level_table .*|level_table          = \"$repo_dir/tests/fixtures/talys-levels.tab\"|" \
  -e 's|^seed .*|seed = 104|' -e 's|"branches.root"|"branches-seed.root"|' \
  "$repo_dir/tests/controls/branches.toml" >"$temp_dir/branches-seed.toml"
(cd "$temp_dir" && "$simulator" "$temp_dir/branches-seed.toml" \
  >"$temp_dir/branches-seed.log" 2>&1) || fail 'simulation failed for branches-seed'
"$root_checks" independent "$temp_dir/branches.root" \
  "$temp_dir/branches-seed.root"
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
