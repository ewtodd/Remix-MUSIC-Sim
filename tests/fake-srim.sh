#!/usr/bin/env bash
# Stands in for make-srim-table in the suite: a flat 0.05 MeV/mm stopping table
# in SRIM's layout for any ion, written where --out says.
set -euo pipefail

ion='' out=''
while (($#)); do
  case $1 in
  --ion) ion=$2 ;;
  --out) out=$2 ;;
  esac
  shift 2
done
test -n "$ion" && test -n "$out"
{
  printf '%s\n' " Ion = $ion (fake)" ' Stopping Units =  MeV/mm'
  for energy in 0.001 0.01 0.1 1 10 100 1000 10000; do
    printf '%s MeV   5.000E-02  0.000E+00\n' "$energy"
  done
} >"$out"
