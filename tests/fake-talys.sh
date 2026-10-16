#!/usr/bin/env bash
# Stands in for TALYS in the suite with a small cascade in TALYS's formats;
# its exclusive cross sections follow by hand from the numbers below.
set -euo pipefail

input=$(cat)
# A `fakeerror` keyword fails as TALYS does: an error message, then status 0.
if grep -q '^fakeerror' <<<"$input"; then
  echo ' TALYS-error: fake failure requested by the input'
  exit 0
fi
read -r energy Z A < <(awk '$1 == "energy" { e = $2 } $1 == "element" { z = $2 }
  $1 == "mass" { a = $2 } END { print e, z, a }' <<<"$input")
tag=$(printf '%08.3f' "$energy")

# Cross sections [mb] at this E_cm: (a,a') falls linearly to zero at 7 MeV,
# so a strip range spans open and closed; (a,2n) opens above 12 MeV.
read -r ecm sAA f2 < <(awk -v e="$energy" -v a="$A" 'BEGIN {
  ecm = e * a / (a + 4); s = 40 * (ecm - 7); f = 15 * (ecm - 12)
  printf "%.6f %.6E %.6E\n", ecm, (s > 0 ? s : 0), (f > 0 ? f : 0) }')

# Nuclides: compound C (fed by its first gamma 2 mb at 12 MeV, above its
# neutron threshold, which emits a neutron, and 0.5 mb at the excitation 37Cl
# + a reaches at E_cm, so (a,g) is open only upstream of a strip centre), C-n
# (10 mb ground, 20 mb at 1 MeV, f2 in a bin above its neutron threshold and
# 0.4 f2 just above that, which gamma-decays into it; that bin keeps a quarter
# by gamma), C-2n, C-p (5 mb ground), and the beam nucleus ((a,a') at 1, 2,
# 2.25 MeV).
awk -v Z="$Z" -v A="$A" -v E="$energy" -v ecm="$ecm" -v sAA="$sAA" \
  -v f2="$f2" -f /dev/stdin >"populationE$tag.out" <<'AWK'
function nuclide(id, z, a, nl, size) {
  nz[id] = z; na[id] = a; nnl[id] = nl; nsize[id] = size; nbins[id] = 0
}
function bin(id, b, ex, feed) {
  bex[id, b] = ex; bfeed[id, b] = feed
  if (b + 1 > nbins[id]) nbins[id] = b + 1
}
# A state of nuclide id at bin b decays by type t into bin db of nuclide did.
function flow(id, b, t, did, db, mb) {
  nf[id, b]++; ft[id, b, nf[id, b]] = t; fd[id, b, nf[id, b]] = did
  fb[id, b, nf[id, b]] = db; fmb[id, b, nf[id, b]] = mb; dec[id, b, t] += mb
}
function labels(a,    s) {
  s = a % 2 ? "JP= 0.5          JP= 1.5" : "JP= 0.0          JP= 1.0"
  return "##      bin            Ex           " s
}
function emit(id,    b, p, t, k, d, did, slots) {
  print "# header:"
  print "#   title: fake population of nuclide " id " at " E " MeV"
  print "# reaction:"; print "#   type: (a,x)"; print "#   residual:"
  print "#     Z: " nz[id]; print "#     A: " na[id]; print "#   parameters:"
  print "#     number of discrete levels (NL): " nnl[id]
  printf "#     continuum bin size [MeV]:  %.6E\n", nsize[id]
  # Every state is at the second spin slot and positive parity.
  for (p = -1; p <= 1; p += 2) {
    print "#   quantity:"; print "#     type: population"
    print "#     parity: " p; print "#   datablock:"
    print "#     columns: 6"; print "#     entries: " nbins[id]
    print labels(na[id])
    for (b = 0; b < nbins[id]; b++)
      printf "%6d %15.6E %15.6E %15.6E %15.6E %15.6E\n", b, bex[id, b],
        bfeed[id, b], (p > 0 ? bfeed[id, b] : 0), 0, (p > 0 ? bfeed[id, b] : 0)
  }
  # TALYS starts this table here and continues it between the states.
  print "#   quantity:"; print "#     type: population per bin before decay"
  print "#   datablock:"; print "#     columns: 6"
  print "#     entries: " 2 * nbins[id]
  for (b = nbins[id] - 1; b >= 0; b--) {
    for (p = -1; p <= 1; p += 2)
      printf "%6d %6d %15.6E %15.6E %15.6E %15.6E\n", b, p, bex[id, b],
        pop[id, b], (p > 0 ? pop[id, b] : 0), (p > 0 ? pop[id, b] : 0)
    if (!(pop[id, b] > 0)) continue
    print "#   parameters:"; print "#     bin: " b
    printf "#     excitation energy [MeV]:  %.6E\n", bex[id, b]
    printf "#     bin width [MeV]:  %.6E\n", nsize[id]
    printf "#     spin:  %.6E\n", (na[id] % 2 ? 1.5 : 1.0)
    print "#     parity: 1"
    printf "#     population [mb]:  %.6E\n", pop[id, b]
    print "#   quantity:"; print "#     type: particle decay from this bin"
    print "#   datablock:"; print "#     columns: 2"; print "#     entries: 7"
    for (t = 0; t < 7; t++) printf "%-10s %15.6E\n", tname[t], dec[id, b, t]
    for (t = 0; t < 7; t++) {
      if (!(dec[id, b, t] > 0)) continue
      for (p = -1; p <= 1; p += 2) {
        print "#   parameters:"
        printf "#     spin:  %.6E\n", (na[id] % 2 ? 1.5 : 1.0)
        print "#     parity: 1"
        printf "#     population [mb]:  %.6E\n", pop[id, b]
        printf "#     %s decay [mb]:  %.6E\n", tname[t], dec[id, b, t]
        printf "#     %s decay to P=%s [mb]:  %.6E\n", tname[t],
          (p < 0 ? "-1" : " 1"), (p > 0 ? dec[id, b, t] : 0)
        print "#   quantity:"
        print "#     type: " tname[t] " decay from J,P bin to residual J',P' bin"
        did = -1
        for (k = 1; k <= nf[id, b]; k++) if (ft[id, b, k] == t) did = fd[id, b, k]
        print "#   datablock:"; print "#     columns: 4"
        print "#     entries: " nbins[did]
        # Labelled in this (source) nucleus's spin convention.
        print labels(na[id])
        for (d = 0; d < nbins[did]; d++) {
          slots = 0
          for (k = 1; k <= nf[id, b]; k++)
            if (ft[id, b, k] == t && fb[id, b, k] == d && p > 0)
              slots += fmb[id, b, k]
          printf "%6d %15.6E %15.6E %15.6E\n", d, bex[did, d], 0, slots
        }
      }
    }
  }
}
BEGIN {
  split("gamma neutron proton deuteron triton helium-3 alpha", names, " ")
  for (t = 0; t < 7; t++) tname[t] = names[t + 1]
  zc = Z + 2; ac = A + 4
  nuclide("C", zc, ac, 2, 0.5); bin("C", 0, 0.0, 0.0); bin("C", 1, 12.0, 2.0)
  bin("C", 2, ecm + 6.2229, 0.5)
  nuclide("R1", zc, ac - 1, 1, 0.5)
  bin("R1", 0, 0.0, 10.0); bin("R1", 1, 1.0, 22.0)
  if (f2 > 0) { bin("R1", 2, ecm - 5.0, f2); bin("R1", 3, ecm - 4.5, 0.4 * f2) }
  nuclide("R2", zc, ac - 2, 0, 0.5); bin("R2", 0, 0.0, 1.05 * f2)
  nuclide("Rp", zc - 1, ac - 1, 0, 0.5); bin("Rp", 0, 0.0, 5.0)
  nuclide("B", Z, A, 2, 0.5)
  bin("B", 0, 0.0, 0.0); bin("B", 1, 1.0, 0.5 * sAA)
  bin("B", 2, 2.0, 0.3 * sAA); bin("B", 3, 2.25, 0.2 * sAA)
  # Populations before decay and the flows that carry them down.
  pop["C", 1] = 2.0; flow("C", 1, 1, "R1", 1, 2.0); pop["C", 2] = 0.5
  pop["R1", 3] = 0.4 * f2; pop["R1", 2] = 1.4 * f2
  if (f2 > 0) flow("R1", 3, 0, "R1", 2, 0.4 * f2)
  if (f2 > 0) { flow("R1", 2, 0, "R1", 1, 0.35 * f2); flow("R1", 2, 1, "R2", 0, 1.05 * f2) }
  pop["R1", 1] = 22.0 + 0.35 * f2; flow("R1", 1, 0, "R1", 0, pop["R1", 1])
  pop["R1", 0] = 10.0 + pop["R1", 1]; pop["Rp", 0] = 5.0
  pop["B", 3] = 0.2 * sAA; pop["B", 2] = 0.5 * sAA; pop["B", 1] = sAA
  pop["B", 0] = sAA
  if (sAA > 0) {
    flow("B", 3, 0, "B", 2, pop["B", 3]); flow("B", 2, 0, "B", 1, pop["B", 2])
    flow("B", 1, 0, "B", 0, pop["B", 1])
  }
  emit("C"); emit("R1"); emit("R2"); emit("Rp"); emit("B")
}
AWK

# Exclusive cross sections as TALYS's section 6a writes them: (a,n) is the
# 30 mb fed below the threshold plus the quarter of 1.4 f2 kept by gamma.
awk -v sAA="$sAA" -v f2="$f2" 'BEGIN {
  print "    TALYS-fake (regression-suite stand-in)"; print ""
  print " 6. Exclusive cross sections"; print ""
  print " 6a. Total exclusive cross sections "; print ""
  print "     Emitted particles     cross section reaction"
  print "    n   p   d   t   h   a"
  printf "    0   0   0   0   0   0    %.5E  (a,g)\n", 0.5
  printf "    1   0   0   0   0   0    %.5E  (a,n)\n", 30.0 + 0.35 * f2
  printf "    0   1   0   0   0   0    %.5E  (a,p)\n", 5.0
  if (sAA > 0) printf "    0   0   0   0   0   1    %.5E  (a,a%c)\n", sAA, 39
  if (f2 > 0) printf "    2   0   0   0   0   0    %.5E  (a,2n)\n", 1.05 * f2
  print ""; printf " (n,gn) + (n,gp) +...(n,ga) cross sections:     %.5f\n", 2.0
  print ""; print " The TALYS team congratulates you with this successful calculation."
}' >talys.out

# Elastic: Rutherford to 150 deg, a quarter beyond; all zeros below E_cm 0.05
# MeV, as real TALYS. The (a,g) file is a decoy talys-cache must not read.
awk -v e="$energy" -v z="$Z" -v a="$A" 'BEGIN {
  pi = 3.14159265358979; ecm = e * a / (a + 4); k = 2 * z * 1.439964548 / ecm
  print "# header:"; print "#   title: fake elastic angular distribution"
  for (t = 0; t <= 180; t += 2) {
    ratio = ecm < 0.05 ? 0 : t < 150 ? 1.0 : 0.25
    s = sin(t * pi / 360); direct = t > 0 ? ratio * 10 * (k / 4)^2 / s^4 : ratio * 1e30
    printf "%15.6E%15.6E%15.6E%15.6E%15.6E%15.6E\n", t, direct, direct, 0, ratio, 0
  } }' >"aa${tag}ang.L00"
printf '%s\n' '# decoy' '   0.0 -1 -1 -1 -1' >"ag${tag}ang.L00"
