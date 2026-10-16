#ifndef TALYSCASCADE_HPP
#define TALYSCASCADE_HPP

#include <istream>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <Rtypes.h>

#include "ExcitationTable.hpp"

namespace music {

// Emission types in TALYS's order; "gnpdtha" names them in channel strings.
const Int_t kEmissionTypes = 7;
const char kEmissionLetters[] = "gnpdtha";

// One residue state: a TALYS bin at a spin (2J) and parity.
struct CascadeKey {
  Int_t Z = 0, N = 0, bin = 0, twoJ = 0, parity = 0;
  bool operator<(const CascadeKey &other) const;
};

struct CascadeFlow {
  Int_t type = 0; // index into kEmissionLetters
  CascadeKey dest;
  Double_t mb = 0.0;
};

struct CascadeState {
  Double_t pop = 0.0;                  // population before decay [mb]
  Double_t decay[kEmissionTypes] = {}; // its decay by each type [mb]
  std::vector<CascadeFlow> flows;      // where that decay goes
};

// A residue's bin layout: levels are bins 0..levels, then continuum bins of
// width binSize centred at their energies.
struct CascadeNuclide {
  Int_t levels = -1;
  Double_t binSize = 0.0;
  std::map<Int_t, Double_t> binEx;
};

// TALYS's decay cascade from a population file written with outpopulation
// and outdecay, and each state's feeding by the compound's first emission.
struct Cascade {
  Int_t compoundZ = 0, compoundN = 0;
  std::map<CascadeKey, CascadeState> states;
  std::map<CascadeKey, Double_t> feeding;
  std::map<std::pair<Int_t, Int_t>, CascadeNuclide> nuclides;
};

// Streams the population file; the compound is TALYS's projectile plus target.
Bool_t ReadCascade(std::istream &input, Int_t compoundZ, Int_t compoundN,
                   Cascade &cascade, std::string &error);

// An exclusive channel, as the particles emitted in order ("nn", "np", "g").
// Step k's table maps the bins of the residue its k-th emission leaves to mb.
struct CascadeChannel {
  std::string emitted;
  Double_t sigmaMb = 0.0;
  std::vector<std::pair<Int_t, Int_t>> residues; // (Z, N) per step
  std::vector<std::map<Int_t, Double_t>> steps;
};

// Every ordered channel with a positive cross section.
std::vector<CascadeChannel> CascadeChannels(const Cascade &cascade);

// TALYS's section 6a: exclusive cross sections keyed by the sorted letters of
// the particles emitted ("nn", "np", "" for (a,g)); gammaFirst is TALYS's
// separate total of channels that start with a gamma.
Bool_t ExclusiveCrossSections(const std::string &talysOut,
                              std::map<std::string, Double_t> &sigma,
                              Double_t &gammaFirst, std::string &error);

// Checks the channels against section 6a to 0.1%, or names the first miss.
Bool_t CheckChannels(const std::vector<CascadeChannel> &channels,
                     const std::map<std::string, Double_t> &sigma,
                     Double_t gammaFirst, std::string &error);

// Channel-table text: a `channel` line per channel, then per step a `step`
// line and its level and bin rows (mb).
std::string ChannelTableRows(const Cascade &cascade,
                             const std::vector<CascadeChannel> &channels);

// A channel read back from a talys-cache table; residues are (Z, A).
struct TalysChannel {
  std::string emitted;
  Double_t sigmaMb = 0.0;
  std::vector<std::pair<Int_t, Int_t>> residues;
  std::vector<ExcitationTable> steps;
};

// Reads a table's channels; comments and `elastic` rows are skipped.
Bool_t ParseChannelTable(const std::string &text, const std::string &source,
                         std::vector<TalysChannel> &channels,
                         std::string &error);

} // namespace music

#endif
