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

printf '%s\n' '[1/8] validating CLI and strict schema'
expect_failure "$simulator"
expect_failure "$simulator" --check "$temp_dir/missing.toml"
expect_failure "$simulator" --check \
  "$repo_dir/tests/controls/invalid-zero-step.toml"
while IFS= read -r control; do
  "$simulator" --check "$control" >"$temp_dir/check.log" 2>&1 || {
    sed -n '1,160p' "$temp_dir/check.log" >&2
    fail "schema validation failed for $control"
  }
done < <(printf '%s\n' "$repo_dir/basic.toml"; \
         find "$repo_dir/ControlExamples" -type f -name '*.toml' | sort; \
         find "$repo_dir/tests/controls" -type f -name '*.toml' ! -name 'invalid-*' | sort)

printf '%s\n' '[2/8] testing legacy converter safety and escaping'
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

printf '%s\n' '[3/8] testing nuclide loaders (ERROR lines below are expected: malformed-input assertions)'
"$nuclide_checks" "$temp_dir"

printf '%s\n' '[4/8] validating Vavilov interpolation and limits'
"$vavilov_checks"

run_control() {
  local control=$1
  (cd "$temp_dir" && "$simulator" "$repo_dir/tests/controls/$control.toml" \
    >"$temp_dir/$control.log" 2>&1) || {
    sed -n '1,200p' "$temp_dir/$control.log" >&2
    fail "simulation failed for $control"
  }
}

printf '%s\n' '[5/8] testing stopping-energy conservation and forbidden reactions'
run_control stopping
"$root_checks" standard "$temp_dir/stopping.root" 1
"$root_checks" conservation "$temp_dir/stopping.root"
run_control forbidden-reaction
"$root_checks" standard "$temp_dir/forbidden.root" 1
"$root_checks" forbidden "$temp_dir/forbidden.root"
run_control later-step-forbidden
"$root_checks" standard "$temp_dir/later-step-forbidden.root" 1
"$root_checks" forbidden "$temp_dir/later-step-forbidden.root"

printf '%s\n' '[6/8] testing multi-step truth angles and boundary convergence'
run_control multistep
"$root_checks" standard "$temp_dir/multistep.root" 1
"$root_checks" multistep "$temp_dir/multistep.root"
run_control boundary-coarse
run_control boundary-fine
"$root_checks" boundary "$temp_dir/boundary-coarse.root" \
  "$temp_dir/boundary-fine.root"

printf '%s\n' '[7/8] testing thread-count-independent events and merge output'
run_control repro-single
run_control repro-multi
"$root_checks" standard "$temp_dir/repro-single.root" 4
"$root_checks" standard "$temp_dir/repro-multi.root" 4
"$root_checks" compare "$temp_dir/repro-single.root" \
  "$temp_dir/repro-multi.root"

printf '%s\n' '[8/8] testing checked, atomic update output'
run_control update
run_control update
"$root_checks" standard "$temp_dir/update.root" 4
"$root_checks" unique "$temp_dir/update.root"

printf '%s\n' 'All tests passed.'
