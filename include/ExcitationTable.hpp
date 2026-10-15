#ifndef EXCITATIONTABLE_HPP
#define EXCITATIONTABLE_HPP

#include <array>
#include <string>
#include <utility>
#include <vector>

#include <Rtypes.h>
#include <TRandom.h>

namespace music {
// Residue excitation distribution read from a Talys-derived level table:
// point masses (`ground`/`level` rows) plus uniform continuum bins (`bin`).
struct ExcitationTable {
  std::vector<std::pair<Double_t, Double_t>> levels; // {E_ex, weight}
  std::vector<std::array<Double_t, 3>> bins;         // {lo, hi, weight}
  Double_t total = 0.0; // sum of all weights; 1 after a successful Parse

  // Validates and normalizes table text read from `source`. On failure returns
  // kFALSE, leaves the table empty, and names source and line in `error`.
  Bool_t Parse(const std::string &text, const std::string &source,
               std::string &error);
  // Reads the file at `path`, then parses it as above.
  Bool_t Load(const std::string &path, std::string &error);
  // Table weight at or below EneAvail; a bin straddling it counts pro rata.
  Double_t AllowedWeight(Double_t EneAvail) const;
  // Samples E_ex from the table truncated to [0, EneAvail] and renormalized,
  // using exactly one RNG draw. With no weight there, ok = kFALSE, no draw.
  Double_t Sample(Double_t EneAvail, TRandom *rng, Bool_t &ok) const;
};
} // namespace music

#endif
