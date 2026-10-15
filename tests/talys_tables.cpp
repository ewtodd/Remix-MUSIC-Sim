#include "ExcitationTable.hpp"
#include "TalysTables.hpp"

#include <TMath.h>
#include <TRandom3.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

Bool_t Fail(const std::string &message) {
  std::cerr << "talys_tables: " << message << std::endl;
  return kFALSE;
}

Bool_t Near(Double_t actual, Double_t expected, Double_t relative) {
  return std::fabs(actual - expected) <= relative * std::fabs(expected);
}

// The text with its first `from` replaced; empty if `from` is absent.
std::string Replace(std::string text, const std::string &from,
                    const std::string &to) {
  const size_t at = text.find(from);
  if (at == std::string::npos)
    return "";
  return text.replace(at, from.size(), to);
}

// The closed form against a direct integral of (k/4E)^2 / sin^4(theta/2) over
// theta_cm in [thetaMin, pi], taken in ln(theta) where the integrand is steep.
Bool_t CheckRutherford() {
  const Double_t thetas[3] = {1.0, 10.0, 120.0};
  for (Int_t t = 0; t < 3; ++t) {
    const Double_t ecm = 8.0;
    const Double_t k = 17 * 2 * 1.439964548 / ecm;
    const Double_t lo = std::log(thetas[t] * TMath::Pi() / 180.0);
    const Double_t hi = std::log(TMath::Pi());
    const Int_t steps = 200000;
    const Double_t h = (hi - lo) / steps;
    Double_t integral = 0.0;
    for (Int_t i = 0; i <= steps; ++i) {
      const Double_t theta = std::exp(lo + i * h);
      const Double_t s = std::sin(0.5 * theta);
      const Double_t f = 2.0 * TMath::Pi() * std::sin(theta) * theta *
                         std::pow(k / 4.0, 2) / std::pow(s, 4);
      const Double_t simpson =
          (i == 0 || i == steps) ? 1.0 : (i % 2 ? 4.0 : 2.0);
      integral += simpson * f;
    }
    integral *= 10.0 * h / 3.0; // fm^2 -> mb
    const Double_t closed = music::RutherfordAboveMb(17, 2, ecm, thetas[t]);
    if (!Near(closed, integral, 1e-7))
      return Fail("Rutherford integral above " + std::to_string(thetas[t]) +
                  " deg is " + std::to_string(closed) + " mb, not " +
                  std::to_string(integral));
  }
  return kTRUE;
}

Bool_t CheckKeywords() {
  if (!music::TalysKeywordProblem("ldmodel 2").empty())
    return Fail("an ordinary keyword was refused");
  const std::vector<std::string> refused = {
      "energy 5",   "  Projectile a", "alphaomp 2", "outpopulation n",
      "mass 37",    "element 17",     "",           "bins 60 # coarse",
      "outangle n", "fileelastic n"};
  for (size_t i = 0; i < refused.size(); ++i)
    if (music::TalysKeywordProblem(refused[i]).empty())
      return Fail("keyword '" + refused[i] + "' was not refused");
  return kTRUE;
}

// The deck's lines in order, without line terminators.
std::vector<std::string> Lines(const std::string &text) {
  std::vector<std::string> lines;
  std::istringstream input(text);
  std::string line;
  while (std::getline(input, line))
    lines.push_back(line);
  return lines;
}

Bool_t CheckDeck() {
  const std::vector<std::string> deck =
      Lines(music::TalysDeck(17, 37, 9.396226, 2, "", {"bins 60"}));
  const std::vector<std::string> defaults = music::TalysDefaultKeywords();
  const size_t expected = 4 + defaults.size() + 5;
  if (deck.size() != expected || deck[0] != "projectile a" ||
      deck[1] != "element 17" || deck[2] != "mass 37" ||
      deck[3] != "energy 9.396226" || deck[4] != defaults[0] ||
      deck[expected - 5] != "alphaomp 2" || deck[expected - 4] != "bins 60" ||
      deck[expected - 3] != "outpopulation y" ||
      deck[expected - 2] != "outangle y" ||
      deck[expected - 1] != "fileelastic y")
    return Fail("the TALYS input has the wrong lines or order");
  // A control's keywords follow the defaults, so TALYS keeps the control's.
  const std::vector<std::string> atomki =
      Lines(music::TalysDeck(37, 87, 14.6448, 9, "0123abcd", {}));
  if (atomki[0] != "# alphaomp9real.gnu md5 = 0123abcd")
    return Fail("an Atomki-V2 input does not record its potential's MD5");
  const std::vector<std::string> comment =
      Lines(music::TalysDeckComment(music::TalysDeck(17, 37, 9.0, 2, "", {})));
  if (comment.size() != expected || comment[0] != "# TALYS input:" ||
      comment[1] != "#   projectile a")
    return Fail("the input comment block is malformed");
  return kTRUE;
}

Bool_t ExpectRowsError(const std::string &talysOut,
                       const std::string &population, Int_t Z, Int_t A,
                       const std::string &message, const std::string &what) {
  std::string rows, error;
  Double_t sigma = 0.0;
  if (talysOut.empty())
    return Fail(what + ": the fixture edit did not apply");
  if (music::TalysLevelRows(talysOut, population, Z, A, rows, sigma, error))
    return Fail(what + " was accepted");
  if (error.find(message) == std::string::npos)
    return Fail(what + " failed with '" + error + "', not '" + message + "'");
  return kTRUE;
}

// Real TALYS 2.2 output for alpha + 37Cl at E_cm = 8.65 MeV, trimmed, behind
// a synthetic 37Ar isobar; it and the 41K and second parity blocks are decoys.
Bool_t CheckRows(const std::string &fixtures) {
  std::string talysOut, population;
  if (!music::ReadFileText(fixtures + "/talys-cl37.out", talysOut) ||
      !music::ReadFileText(fixtures + "/talys-cl37-population.out", population))
    return Fail("cannot read the TALYS fixtures in " + fixtures);

  std::string rows, error;
  Double_t sigma = 0.0;
  if (!music::TalysLevelRows(talysOut, population, 17, 37, rows, sigma, error))
    return Fail("the 37Cl fixture was refused: " + error);
  if (!Near(sigma, 1.435315186e+01, 1e-9))
    return Fail("the 37Cl cross section is " + std::to_string(sigma));
  const std::vector<std::string> lines = Lines(rows);
  if (lines.size() != 100 || lines[0] != "level 1.726580 3.319079000e+00")
    return Fail("the 37Cl rows are wrong");
  music::ExcitationTable table;
  if (!table.Parse(rows, "fixture rows", error))
    return Fail("the rows do not parse as a level table: " + error);
  if (table.levels.size() != 30 || table.bins.size() != 70 ||
      !Near(table.sum, sigma, 1e-9) || table.bins[0][0] != 5.2287)
    return Fail("the rows parse to the wrong table");

  std::string closedRows;
  const std::string closed =
      Replace(talysOut, "alpha   = 1.43531E+01", "alpha   = 0.00000E+00");
  if (closed.empty() ||
      !music::TalysLevelRows(closed, population, 17, 37, closedRows, sigma,
                             error) ||
      !closedRows.empty() || sigma != 0.0)
    return Fail("a closed channel did not give an empty table");
  music::ExcitationTable empty;
  Bool_t ok = kTRUE;
  TRandom3 rng(1);
  if (!empty.Parse(closedRows, "closed", error, kTRUE) || empty.sum != 0.0 ||
      empty.AllowedWeight(100.0) != 0.0 ||
      empty.Sample(100.0, &rng, ok) != 0.0 || ok)
    return Fail("an empty talys-cache table is not a closed channel");

  return ExpectRowsError(Replace(talysOut, "alpha   = 1.43531E+01",
                                 "alpha   = 1.53531E+01"),
                         population, 17, 37, "does not match",
                         "a population that misses the binary alpha total") &&
         ExpectRowsError(talysOut, population, 18, 40, "no population of",
                         "a residue absent from the population file") &&
         ExpectRowsError(Replace(talysOut, "congratulates", "regrets"),
                         population, 17, 37, "TALYS did not finish",
                         "an unfinished TALYS run") &&
         ExpectRowsError(talysOut,
                         Replace(population, "bin size [MeV]:  4.273605E-02",
                                 "bin size [MeV]:  5.273605E-02"),
                         17, 37, "do not start at the last discrete level",
                         "bins that are not centred where TALYS puts them");
}

// Real TALYS 2.2 elastic angular distribution for alpha + 37Cl at E_cm =
// 8.6476 MeV: shape elastic matches Rutherford, compound elastic is added.
Bool_t CheckElasticRows(const std::string &fixtures) {
  std::string angular;
  if (!music::ReadFileText(fixtures + "/talys-cl37-elastic.ang", angular))
    return Fail("cannot read the elastic fixture in " + fixtures);
  const Double_t ecm = 8.647571;
  std::string rows, error;
  Bool_t rutherfordOnly = kTRUE;
  if (!music::TalysElasticRows(angular, 17, 2, ecm, rows, rutherfordOnly,
                               error) ||
      rutherfordOnly)
    return Fail("the 37Cl elastic fixture was refused: " + error);
  music::ElasticTable table;
  if (!table.Parse(rows, "elastic rows", error) || table.theta.size() != 91 ||
      table.ratio[0] != 1.0)
    return Fail("the elastic rows parse to the wrong table: " + error);
  // At 180 deg: TALYS's shape ratio 0.150617 plus compound 1.897182 mb/sr.
  const Double_t expected =
      0.150617 + 1.897182 / music::RutherfordDsigmaMb(17, 2, ecm, 180.0);
  if (!Near(table.ratio.back(), expected, 1e-9))
    return Fail("the 180 deg ratio is " + std::to_string(table.ratio.back()) +
                ", not shape plus compound " + std::to_string(expected));

  // At another E_cm, TALYS's shape elastic is not our Rutherford.
  if (music::TalysElasticRows(angular, 17, 2, 1.05 * ecm, rows, rutherfordOnly,
                              error) ||
      error.find("implies Rutherford") == std::string::npos)
    return Fail("an elastic file for another E_cm was accepted: " + error);
  std::ostringstream zeros;
  zeros << "0 0 0 0 0" << std::endl
        << "90 0 0 0 0" << std::endl
        << "180 0 0 0 0" << std::endl;
  if (!music::TalysElasticRows(zeros.str(), 17, 2, ecm, rows, rutherfordOnly,
                               error) ||
      !rutherfordOnly || !rows.empty())
    return Fail("an all-zero distribution did not fall back to Rutherford");
  if (music::TalysElasticRows("0 1 1 0 1", 17, 2, ecm, rows, rutherfordOnly,
                              error))
    return Fail("a distribution that stops short of 180 deg was accepted");
  return kTRUE;
}

// Rutherford x ratio over thetaMin..180, integrated numerically with the
// ratio held at its cell mean, as the table defines it [mb].
Double_t NumericElastic(const music::ElasticTable &table, Double_t thetaMin) {
  const Double_t deg = TMath::Pi() / 180.0;
  Double_t sum = 0.0;
  for (size_t i = 0; i + 1 < table.theta.size(); ++i) {
    Double_t a = table.theta[i];
    const Double_t b = table.theta[i + 1];
    if (b <= thetaMin)
      continue;
    Double_t ra = table.ratio[i];
    if (a < thetaMin) {
      ra += (table.ratio[i + 1] - ra) * (thetaMin - a) / (b - a);
      a = thetaMin;
    }
    const Int_t steps = 2000;
    const Double_t h = (b - a) / steps;
    Double_t cell = 0.0;
    for (Int_t s = 0; s <= steps; ++s) {
      const Double_t theta = a + s * h;
      const Double_t simpson =
          (s == 0 || s == steps) ? 1.0 : (s % 2 ? 4.0 : 2.0);
      cell += simpson * music::RutherfordDsigmaMb(17, 2, 8.0, theta) * 2.0 *
              TMath::Pi() * std::sin(theta * deg);
    }
    sum += cell * h * deg / 3.0 * 0.5 * (ra + table.ratio[i + 1]);
  }
  return sum;
}

Bool_t CheckElasticTable(const std::string &fixtures) {
  // No rows is ratio 1: pure Rutherford, drawn exactly as before.
  music::ElasticTable rutherford;
  if (!Near(rutherford.AboveMb(17, 2, 8.0, 7.0),
            music::RutherfordAboveMb(17, 2, 8.0, 7.0), 1e-12))
    return Fail("an empty elastic table is not pure Rutherford");
  TRandom3 table(5), direct(5);
  for (Int_t i = 0; i < 1000; ++i) {
    const Double_t umin =
        std::pow(std::sin(0.5 * 7.0 * TMath::Pi() / 180.0), 2);
    const Double_t u =
        1.0 / (1.0 / umin - direct.Uniform(0.0, 1.0) * (1.0 / umin - 1.0));
    const Double_t expected = 2.0 * std::asin(std::sqrt(std::min(1.0, u)));
    if (!Near(rutherford.Sample(7.0, &table), expected, 1e-9))
      return Fail("an empty elastic table samples differently from Rutherford");
  }

  std::string angular, rows, error;
  Bool_t rutherfordOnly = kFALSE;
  music::ElasticTable real;
  if (!music::ReadFileText(fixtures + "/talys-cl37-elastic.ang", angular) ||
      !music::TalysElasticRows(angular, 17, 2, 8.647571, rows, rutherfordOnly,
                               error) ||
      !real.Parse(rows, "elastic rows", error))
    return Fail("cannot build the 37Cl elastic table: " + error);
  const Double_t thetas[3] = {1.0, 10.0, 63.0};
  for (Int_t t = 0; t < 3; ++t)
    if (!Near(real.AboveMb(17, 2, 8.0, thetas[t]),
              NumericElastic(real, thetas[t]), 1e-6))
      return Fail("the elastic cross section above " +
                  std::to_string(thetas[t]) + " deg is not its integral");

  // Sampled angles follow the table, with one draw each.
  const Double_t edges[4] = {20.0, 40.0, 90.0, 180.0};
  Long64_t counts[3] = {0, 0, 0};
  const Long64_t samples = 400000;
  TRandom3 rng(11), reference(11);
  for (Long64_t i = 0; i < samples; ++i) {
    const Double_t theta = real.Sample(20.0, &rng) * 180.0 / TMath::Pi();
    for (Int_t band = 0; band < 3; ++band)
      if (theta >= edges[band] && theta <= edges[band + 1])
        ++counts[band];
    reference.Rndm();
  }
  if (rng.Rndm() != reference.Rndm())
    return Fail("elastic sampling does not take exactly one draw");
  const Double_t total = real.AboveMb(17, 2, 8.0, 20.0);
  for (Int_t band = 0; band < 3; ++band) {
    const Double_t p =
        (real.AboveMb(17, 2, 8.0, edges[band]) -
         (band < 2 ? real.AboveMb(17, 2, 8.0, edges[band + 1]) : 0.0)) /
        total;
    const Double_t fraction = static_cast<Double_t>(counts[band]) / samples;
    if (std::fabs(fraction - p) > 5.0 * std::sqrt(p * (1.0 - p) / samples))
      return Fail("elastic angles between " + std::to_string(edges[band]) +
                  " and " + std::to_string(edges[band + 1]) + " deg are " +
                  std::to_string(fraction) + " of draws, expected " +
                  std::to_string(p));
  }
  return kTRUE;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "usage: talys_tables <fixtures-dir>" << std::endl;
    return 2;
  }
  if (!CheckRutherford() || !CheckKeywords() || !CheckDeck() ||
      !CheckRows(argv[1]) || !CheckElasticRows(argv[1]) ||
      !CheckElasticTable(argv[1]))
    return 1;
  std::cout << "talys_tables: all checks passed" << std::endl;
  return 0;
}
