#include "TalysTables.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <sstream>

#include <TMD5.h>
#include <TMath.h>

#ifndef MUSICSIM_TALYS_ATOMKI_V2_DIR
#define MUSICSIM_TALYS_ATOMKI_V2_DIR ""
#endif

namespace {

// Parses one whole token as a finite number.
Bool_t ParseNumber(const std::string &token, Double_t &value) {
  size_t consumed = 0;
  try {
    value = std::stod(token, &consumed);
  } catch (const std::exception &) {
    return kFALSE;
  }
  return consumed == token.size() && std::isfinite(value);
}

} // namespace

namespace music {

const std::vector<std::string> &TalysDefaultKeywords() {
  static const std::vector<std::string> keywords = {
      "ldmodel 1",     "maxlevelstar 30", "transpower 20",    "ecisstep 0.02",
      "xseps 1.e-30",  "transeps 1.e-30", "popeps 1.e-30",    "elow 1.e-06",
      "widthfluc 20",  "widthmode 1",     "preequilibrium n", "bins 80",
      "equidistant y", "segment 4"};
  return keywords;
}

std::string TalysKeywordProblem(const std::string &keyword) {
  for (size_t i = 0; i < keyword.size(); ++i)
    if (keyword[i] == '\n' || keyword[i] == '\r' || keyword[i] == '#')
      return "must be one keyword line without '#'";
  std::istringstream fields(keyword);
  std::string name;
  if (!(fields >> name))
    return "is empty";
  for (size_t i = 0; i < name.size(); ++i)
    name[i] =
        static_cast<char>(std::tolower(static_cast<unsigned char>(name[i])));
  if (name == "alphaomp")
    return "sets alphaomp; use talys.alpha_omp";
  if (name == "projectile" || name == "element" || name == "mass" ||
      name == "energy" || name == "outpopulation" || name == "outangle" ||
      name == "fileelastic" || name == "channels" || name == "outdecay")
    return "sets '" + name + "', which talys-cache controls";
  return "";
}

std::string AtomkiV2Dir() {
  const char *env = std::getenv("TALYS_ATOMKI_V2_DIR");
  if (env && env[0] != '\0')
    return env;
  return MUSICSIM_TALYS_ATOMKI_V2_DIR;
}

std::string AtomkiV2Potential(const std::string &dir, Int_t Z, Int_t A) {
  if (dir.empty())
    return "";
  char name[64];
  std::snprintf(name, sizeof(name), "/z%03da%03da_talys_alphaomp9real.gnu", Z,
                A);
  const std::string path = dir + name;
  std::ifstream file(path);
  return file ? path : "";
}

std::string TalysDeck(Int_t Z, Int_t A, Double_t alphaLabMeV, Int_t alphaOmp,
                      const std::string &potentialMd5,
                      const std::vector<std::string> &keywords) {
  std::ostringstream deck;
  if (alphaOmp == 9)
    deck << "# alphaomp9real.gnu md5 = " << potentialMd5 << std::endl;
  char energy[64];
  std::snprintf(energy, sizeof(energy), "energy %.6f", alphaLabMeV);
  deck << "projectile a" << std::endl
       << "element " << Z << std::endl
       << "mass " << A << std::endl
       << energy << std::endl;
  const std::vector<std::string> &defaults = TalysDefaultKeywords();
  for (size_t i = 0; i < defaults.size(); ++i)
    deck << defaults[i] << std::endl;
  deck << "alphaomp " << alphaOmp << std::endl;
  for (size_t i = 0; i < keywords.size(); ++i)
    deck << keywords[i] << std::endl;
  deck << "outpopulation y" << std::endl
       << "outangle y" << std::endl
       << "fileelastic y" << std::endl
       << "channels y" << std::endl
       << "outdecay y" << std::endl;
  return deck.str();
}

std::string TalysDeckComment(const std::string &deck) {
  std::ostringstream comment;
  comment << "# TALYS input:" << std::endl;
  std::istringstream lines(deck);
  std::string line;
  while (std::getline(lines, line))
    comment << "#   " << line << std::endl;
  return comment.str();
}

Double_t RutherfordAboveMb(Int_t Z1, Int_t Z2, Double_t ecm,
                           Double_t thetaMinDeg) {
  const Double_t e2 = 1.439964548; // MeV fm
  const Double_t k = Z1 * Z2 * e2 / ecm;
  const Double_t half = 0.5 * thetaMinDeg * TMath::Pi() / 180.0;
  const Double_t cot = std::cos(half) / std::sin(half);
  // (pi/4) k^2 cot^2(theta_min/2) fm^2; 1 fm^2 = 10 mb.
  return 10.0 * 0.25 * TMath::Pi() * k * k * cot * cot;
}

Double_t RutherfordDsigmaMb(Int_t Z1, Int_t Z2, Double_t ecm,
                            Double_t thetaDeg) {
  const Double_t k = Z1 * Z2 * 1.439964548 / ecm;
  const Double_t s = std::sin(0.5 * thetaDeg * TMath::Pi() / 180.0);
  return 10.0 * (k / 4.0) * (k / 4.0) / (s * s * s * s);
}

Bool_t TalysElasticRows(const std::string &angular, Int_t Z1, Int_t Z2,
                        Double_t ecm, std::string &rows, Bool_t &rutherfordOnly,
                        std::string &error) {
  rows.clear();
  rutherfordOnly = kFALSE;
  std::vector<Double_t> theta, direct, compound, ratio;
  std::istringstream input(angular);
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '#')
      continue;
    std::istringstream fields(line);
    std::string tokens[5];
    Double_t values[5];
    for (Int_t i = 0; i < 5; ++i)
      if (!(fields >> tokens[i]) || !ParseNumber(tokens[i], values[i]) ||
          (i > 0 && values[i] < 0.0)) {
        error = "malformed elastic angular row '" + line + "'";
        return kFALSE;
      }
    if (!theta.empty() && !(values[0] > theta.back())) {
      error = "elastic angles are not increasing at '" + line + "'";
      return kFALSE;
    }
    theta.push_back(values[0]);
    direct.push_back(values[2]);
    compound.push_back(values[3]);
    ratio.push_back(values[4]);
  }
  if (theta.size() < 2 || theta.front() != 0.0 || theta.back() != 180.0) {
    error = "the elastic angular distribution does not run from 0 to 180 deg";
    return kFALSE;
  }

  Bool_t anyRatio = kFALSE;
  for (size_t i = 0; i < theta.size(); ++i)
    anyRatio = anyRatio || ratio[i] > 0.0 || compound[i] > 0.0;
  if (!anyRatio) {
    rutherfordOnly = kTRUE;
    return kTRUE;
  }
  // TALYS's ratio is shape elastic over Rutherford, so shape elastic over it
  // must be our Rutherford: this catches a lab-frame angle or another E_cm.
  for (size_t i = 0; i < theta.size(); ++i) {
    if (theta[i] < 30.0 || theta[i] > 150.0 || !(ratio[i] > 0.0))
      continue;
    const Double_t talys = direct[i] / ratio[i];
    const Double_t ours = RutherfordDsigmaMb(Z1, Z2, ecm, theta[i]);
    if (std::fabs(talys - ours) > 0.01 * ours) {
      char numbers[192];
      std::snprintf(numbers, sizeof(numbers),
                    "TALYS's shape elastic at %g deg implies "
                    "Rutherford %.6e mb/sr, not %.6e",
                    theta[i], talys, ours);
      error = numbers;
      return kFALSE;
    }
  }
  // Compound elastic also leaves the residue in its ground state, so it counts
  // as elastic (nothing at 0 deg, where Rutherford diverges).
  std::ostringstream text;
  char formatted[96];
  for (size_t i = 0; i < theta.size(); ++i) {
    const Double_t total =
        theta[i] > 0.0
            ? ratio[i] + compound[i] / RutherfordDsigmaMb(Z1, Z2, ecm, theta[i])
            : ratio[i];
    std::snprintf(formatted, sizeof(formatted), "elastic %.4f %.9e", theta[i],
                  total);
    text << formatted << std::endl;
  }
  rows = text.str();
  return kTRUE;
}

Bool_t ElasticTable::Parse(const std::string &text, const std::string &source,
                           std::string &error) {
  theta.clear();
  ratio.clear();
  std::istringstream input(text);
  std::string line;
  Int_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    std::istringstream fields(line);
    std::string tag, angle, value, extra;
    if (!(fields >> tag) || tag != "elastic")
      continue;
    const std::string where = source + ":" + std::to_string(lineNumber);
    Double_t a = 0.0, r = 0.0;
    if (!(fields >> angle >> value) || (fields >> extra) ||
        !ParseNumber(angle, a) || !ParseNumber(value, r) || r < 0.0 ||
        (!theta.empty() && !(a > theta.back()))) {
      error = where + ": malformed or out-of-order 'elastic' row";
      theta.clear();
      ratio.clear();
      return kFALSE;
    }
    theta.push_back(a);
    ratio.push_back(r);
  }
  if (!theta.empty() &&
      (theta.size() < 2 || theta.front() != 0.0 || theta.back() != 180.0)) {
    error = source + ": 'elastic' rows must run from 0 to 180 deg";
    theta.clear();
    ratio.clear();
    return kFALSE;
  }
  return kTRUE;
}

Double_t ElasticTable::CellWeight(size_t i, Double_t thetaMinDeg,
                                  Double_t &uLow, Double_t &uHigh) const {
  Double_t a = theta.empty() ? 0.0 : theta[i];
  const Double_t b = theta.empty() ? 180.0 : theta[i + 1];
  Double_t ra = ratio.empty() ? 1.0 : ratio[i];
  const Double_t rb = ratio.empty() ? 1.0 : ratio[i + 1];
  if (b <= thetaMinDeg)
    return 0.0;
  if (a < thetaMinDeg) {
    ra += (rb - ra) * (thetaMinDeg - a) / (b - a);
    a = thetaMinDeg;
  }
  const Double_t deg = TMath::Pi() / 180.0;
  uLow = std::pow(std::sin(0.5 * a * deg), 2);
  uHigh = std::pow(std::sin(0.5 * b * deg), 2);
  return (1.0 / uLow - 1.0 / uHigh) * 0.5 * (ra + rb);
}

Double_t ElasticTable::AboveMb(Int_t Z1, Int_t Z2, Double_t ecm,
                               Double_t thetaMinDeg) const {
  Double_t sum = 0.0, uLow = 0.0, uHigh = 0.0;
  for (size_t i = 0; i < CellCount(); ++i)
    sum += CellWeight(i, thetaMinDeg, uLow, uHigh);
  // Rutherford over a cell is 4pi (k/4)^2 (1/uLow - 1/uHigh) fm^2.
  const Double_t k = Z1 * Z2 * 1.439964548 / ecm;
  return 10.0 * 0.25 * TMath::Pi() * k * k * sum;
}

Double_t ElasticTable::Sample(Double_t thetaMinDeg, TRandom *rng) const {
  Double_t total = 0.0, uLow = 0.0, uHigh = 0.0;
  for (size_t i = 0; i < CellCount(); ++i)
    total += CellWeight(i, thetaMinDeg, uLow, uHigh);
  Double_t pick = rng->Uniform(0.0, total);
  if (!(total > 0.0))
    return thetaMinDeg * TMath::Pi() / 180.0;
  // Rounding can leave pick past every cell: then the last cell's top edge.
  Double_t lastLow = 0.0, lastHigh = 0.0, fraction = 1.0;
  for (size_t i = 0; i < CellCount(); ++i) {
    const Double_t weight = CellWeight(i, thetaMinDeg, uLow, uHigh);
    if (weight <= 0.0)
      continue;
    lastLow = uLow;
    lastHigh = uHigh;
    if (pick < weight) {
      fraction = pick / weight;
      break;
    }
    pick -= weight;
  }
  // Exact Rutherford within the cell: 1/u is uniform between its edges.
  const Double_t u =
      1.0 / (1.0 / lastLow - fraction * (1.0 / lastLow - 1.0 / lastHigh));
  return 2.0 * std::asin(std::sqrt(std::min(1.0, u)));
}

std::string Md5Hex(const std::string &text) {
  TMD5 md5;
  md5.Update(reinterpret_cast<const UChar_t *>(text.data()),
             static_cast<UInt_t>(text.size()));
  md5.Final();
  return md5.AsString();
}

Bool_t ReadFileText(const std::string &path, std::string &text) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    return kFALSE;
  std::ostringstream buffer;
  buffer << input.rdbuf();
  text = buffer.str();
  return kTRUE;
}

} // namespace music
