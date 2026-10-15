#!/usr/bin/env bash
# Stands in for TALYS in the suite: (a,a') falls linearly to zero at E_cm =
# 7 MeV, so a strip range spans open and closed channels.
set -euo pipefail

input=$(cat)
read -r energy Z A < <(awk '$1 == "energy" { e = $2 } $1 == "element" { z = $2 }
  $1 == "mass" { a = $2 } END { print e, z, a }' <<<"$input")
sigma=$(awk -v e="$energy" -v a="$A" 'BEGIN {
  s = 40.0 * (e * a / (a + 4) - 7.0); if (s < 0) s = 0; printf "%.6E", s }')
share() { awk -v s="$sigma" -v f="$1" 'BEGIN { printf "%.6E", s * f }'; }

{
  printf '%s\n' '    TALYS-fake (regression-suite stand-in)' '' \
    ' 2. Binary non-elastic cross sections (non-exclusive)' '' \
    " alpha   = $sigma" '' \
    ' The TALYS team congratulates you with this successful calculation.'
} >talys.out

# Levels at 1 and 2 MeV, then one continuum bin centred at 2.25 MeV, as
# TALYS lays them out: the first bin starts at the last level.
{
  printf '%s\n' '# header:' '#   title: fake population' '# reaction:' \
    '#   type: (a,x)' '#   residual:' "#     Z: $Z" "#     A: $A" \
    '#   parameters:' '#     number of discrete levels (NL): 2' \
    '#     continuum bin size [MeV]:  5.000000E-01' '#   quantity:' \
    '#     type: population' '#   datablock:'
  printf '%6d %15s %15s\n' 0 0.000000E+00 1.000000E-03 \
    1 1.000000E+00 "$(share 0.5)" 2 2.000000E+00 "$(share 0.3)" \
    3 2.250000E+00 "$(share 0.2)"
} >"populationE$(printf '%08.3f' "$energy").out"

# Elastic: Rutherford to 150 deg, a quarter beyond; all zeros below E_cm 0.05
# MeV, as real TALYS. The (a,g) file is a decoy talys-cache must not read.
awk -v e="$energy" -v z="$Z" -v a="$A" 'BEGIN {
  pi = 3.14159265358979; ecm = e * a / (a + 4); k = 2 * z * 1.439964548 / ecm
  print "# header:"; print "#   title: fake elastic angular distribution"
  for (t = 0; t <= 180; t += 2) {
    ratio = ecm < 0.05 ? 0 : t < 150 ? 1.0 : 0.25
    s = sin(t * pi / 360); direct = t > 0 ? ratio * 10 * (k / 4)^2 / s^4 : ratio * 1e30
    printf "%15.6E%15.6E%15.6E%15.6E%15.6E%15.6E\n", t, direct, direct, 0, ratio, 0
  } }' >"aa$(printf '%08.3f' "$energy")ang.L00"
printf '%s\n' '# decoy' '   0.0 -1 -1 -1 -1' >"ag$(printf '%08.3f' "$energy")ang.L00"
