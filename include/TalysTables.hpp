#ifndef TALYSTABLES_HPP
#define TALYSTABLES_HPP

#include <string>
#include <vector>

#include <Rtypes.h>
#include <TRandom.h>

namespace music {

// The control's [talys] block. TALYS runs alpha on the beam nucleus at rest,
// at the alpha energy that gives the strip's CM energy.
struct TalysSettings {
  Int_t alphaOmp = 0; // 0 = 9 (Atomki-V2) if a potential exists, else 2
  std::vector<std::string> keywords; // appended after the defaults
};

// One per-strip table a run needs, named by the MD5 of its exact TALYS input.
struct TalysRequest {
  Int_t strip = -1;
  Bool_t stopped = kFALSE; // the mean beam stops before the strip centre
  Int_t Z = 0; // the beam nucleus: TALYS's target and the residue read back
  Int_t A = 0;
  Double_t beamEnergy = 0.0; // mean beam KE at the strip centre [MeV]
  Double_t ecm = 0.0;        // CM kinetic energy there [MeV]
  std::string deck;          // TALYS input
  std::string potential;     // linked as alphaomp9real.gnu, or empty
  std::string path;          // cached table
};

// Keywords every run starts from: the MUSIC analysis's shared TALYS settings.
const std::vector<std::string> &TalysDefaultKeywords();

// Why a [talys] keyword cannot be used, or empty if it can.
std::string TalysKeywordProblem(const std::string &keyword);

// TALYS_ATOMKI_V2_DIR if set, else the directory the build compiled in.
std::string AtomkiV2Dir();

// The Atomki-V2 real-potential file for (Z, A), or empty if there is none.
std::string AtomkiV2Potential(const std::string &dir, Int_t Z, Int_t A);

// The exact TALYS input. potentialMd5 is recorded when alphaOmp is 9.
std::string TalysDeck(Int_t Z, Int_t A, Double_t alphaLabMeV, Int_t alphaOmp,
                      const std::string &potentialMd5,
                      const std::vector<std::string> &keywords);

// The comment block a cached table carries, so a reader can confirm the table
// was made from this input.
std::string TalysDeckComment(const std::string &deck);

// Level-table rows (mb) for the residue (Z, A) after binary alpha emission;
// a closed channel gives no rows and sigmaMb = 0.
Bool_t TalysLevelRows(const std::string &talysOut,
                      const std::string &population, Int_t Z, Int_t A,
                      std::string &rows, Double_t &sigmaMb, std::string &error);

// Rutherford cross section above theta_cm = thetaMinDeg [mb]: the integral of
// what the rutherford angular model samples.
Double_t RutherfordAboveMb(Int_t Z1, Int_t Z2, Double_t ecm,
                           Double_t thetaMinDeg);

// Rutherford dsigma/dOmega at theta_cm = thetaDeg [mb/sr].
Double_t RutherfordDsigmaMb(Int_t Z1, Int_t Z2, Double_t ecm,
                            Double_t thetaDeg);

// `elastic theta ratio` rows: TALYS's shape plus compound elastic over our
// Rutherford. All zeros (far below the barrier): no rows, rutherfordOnly.
Bool_t TalysElasticRows(const std::string &angular, Int_t Z1, Int_t Z2,
                        Double_t ecm, std::string &rows, Bool_t &rutherfordOnly,
                        std::string &error);

// Elastic as Rutherford times a ratio on a CM-angle grid: the ratio is its
// cell mean, Rutherford exact within a cell. No rows means ratio 1.
struct ElasticTable {
  std::vector<Double_t> theta; // deg, from 0 to 180, increasing
  std::vector<Double_t> ratio;

  // Reads the `elastic` rows of a talys-cache table; other rows are ignored.
  Bool_t Parse(const std::string &text, const std::string &source,
               std::string &error);
  // Cross section above thetaMinDeg [mb].
  Double_t AboveMb(Int_t Z1, Int_t Z2, Double_t ecm,
                   Double_t thetaMinDeg) const;
  // Samples theta_cm [rad] above thetaMinDeg with exactly one RNG draw.
  Double_t Sample(Double_t thetaMinDeg, TRandom *rng) const;

private:
  size_t CellCount() const { return theta.empty() ? 1 : theta.size() - 1; }
  // Cell i's part above thetaMinDeg: its Rutherford weight 1/uLow - 1/uHigh,
  // with u = sin^2(theta/2), times its mean ratio; 0 when wholly below.
  Double_t CellWeight(size_t i, Double_t thetaMinDeg, Double_t &uLow,
                      Double_t &uHigh) const;
};

// Hex MD5 of the bytes.
std::string Md5Hex(const std::string &text);

// Reads a whole file as bytes; kFALSE if it cannot be opened.
Bool_t ReadFileText(const std::string &path, std::string &text);

} // namespace music

#endif
