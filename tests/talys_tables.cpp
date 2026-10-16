#include "ExcitationTable.hpp"
#include "TalysCascade.hpp"
#include "TalysTables.hpp"

#include <TMath.h>
#include <TRandom3.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
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
      "outangle n", "fileelastic n",  "channels n", "outdecay n"};
  for (size_t i = 0; i < refused.size(); ++i)
    if (music::TalysKeywordProblem(refused[i]).empty())
      return Fail("keyword '" + refused[i] + "' was not refused");
  return kTRUE;
}

// The deck's lines in order, without line terminators.
// Lines joined into one text, each ended as a file's would be.
std::string Text(const std::vector<std::string> &lines) {
  std::ostringstream text;
  for (size_t i = 0; i < lines.size(); ++i)
    text << lines[i] << std::endl;
  return text.str();
}

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
  const size_t expected = 4 + defaults.size() + 7;
  if (deck.size() != expected || deck[0] != "projectile a" ||
      deck[1] != "element 17" || deck[2] != "mass 37" ||
      deck[3] != "energy 9.396226" || deck[4] != defaults[0] ||
      deck[expected - 7] != "alphaomp 2" || deck[expected - 6] != "bins 60" ||
      deck[expected - 5] != "outpopulation y" ||
      deck[expected - 4] != "outangle y" ||
      deck[expected - 3] != "fileelastic y" ||
      deck[expected - 2] != "channels y" || deck[expected - 1] != "outdecay y")
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

// One step's table as {bin: mb} in a readable form for messages.
std::string Bins(const std::map<Int_t, Double_t> &table) {
  std::ostringstream out;
  for (std::map<Int_t, Double_t>::const_iterator b = table.begin();
       b != table.end(); ++b)
    out << " " << b->first << ":" << b->second;
  return out.str();
}

// The fake TALYS cascade for 37Cl at E_cm = 13 MeV: its channels are worked
// by hand in tests/fake-talys.sh, and it has all three TALYS format traps.
Bool_t CheckCascade(const std::string &dir) {
  std::string talysOut, error;
  std::ifstream population(dir + "/populationE0014.405.out");
  music::Cascade cascade;
  if (!population || !music::ReadFileText(dir + "/talys.out", talysOut) ||
      !music::ReadCascade(population, 19, 22, cascade, error))
    return Fail("cannot read the fake cascade in " + dir + ": " + error);
  // Trap 1: 22 mb sit in 40K's 1 MeV population, 2 of them via 41K's gamma.
  const music::CascadeKey level{19, 21, 1, 2, 1};
  if (!cascade.feeding.count(level) ||
      !Near(cascade.feeding.at(level), 20.0, 1e-9))
    return Fail("feeding from another nucleus was not removed");

  struct Expected {
    const char *emitted;
    Double_t sigma;
    std::vector<std::map<Int_t, Double_t>> steps;
  };
  const std::vector<Expected> expected = {
      {"a", 240.0, {{{1, 120.0}, {2, 72.0}, {3, 48.0}}}},
      {"g", 0.5, {{{2, 0.5}}}},
      {"gn", 2.0, {{{1, 2.0}}, {{1, 2.0}}}},
      {"n", 35.25, {{{0, 10.0}, {1, 20.0}, {2, 3.75}, {3, 1.5}}}},
      {"nn", 15.75, {{{2, 11.25}, {3, 4.5}}, {{0, 15.75}}}},
      {"p", 5.0, {{{0, 5.0}}}}};
  const std::vector<music::CascadeChannel> channels =
      music::CascadeChannels(cascade);
  if (channels.size() != expected.size())
    return Fail("the fake cascade gives " + std::to_string(channels.size()) +
                " channels, not 6");
  for (size_t c = 0; c < channels.size(); ++c) {
    const music::CascadeChannel &got = channels[c];
    const Expected &want = expected[c];
    Bool_t same = got.emitted == want.emitted &&
                  Near(got.sigmaMb, want.sigma, 1e-6) &&
                  got.steps.size() == want.steps.size();
    for (size_t k = 0; same && k < got.steps.size(); ++k) {
      same = got.steps[k].size() == want.steps[k].size();
      for (std::map<Int_t, Double_t>::const_iterator b = want.steps[k].begin();
           same && b != want.steps[k].end(); ++b)
        same = got.steps[k].count(b->first) &&
               Near(got.steps[k].at(b->first), b->second, 1e-6);
    }
    if (!same)
      return Fail("channel " + got.emitted + " is " +
                  std::to_string(got.sigmaMb) + " mb with step 0" +
                  Bins(got.steps[0]) + "; expected " + want.emitted);
  }

  std::map<std::string, Double_t> sigma;
  Double_t gammaFirst = 0.0;
  if (!music::ExclusiveCrossSections(talysOut, sigma, gammaFirst, error) ||
      !music::CheckChannels(channels, sigma, gammaFirst, error))
    return Fail("the fake cascade fails TALYS's 6a check: " + error);
  std::map<std::string, Double_t> wrong = sigma;
  wrong["n"] += 1.0;
  if (music::CheckChannels(channels, wrong, gammaFirst, error) ||
      error.find("(a,n)") == std::string::npos)
    return Fail("a cascade that misses 6a's (a,n) was accepted");
  if (music::CheckChannels(channels, sigma, gammaFirst + 1.0, error))
    return Fail("a cascade that misses 6a's gamma-first total was accepted");

  // The table text reads back to the same channels and bins.
  std::vector<music::TalysChannel> read;
  const std::string rows = music::ChannelTableRows(cascade, channels);
  if (!music::ParseChannelTable(rows, "fake table", read, error) ||
      read.size() != channels.size())
    return Fail("the channel table does not read back: " + error);
  const music::TalysChannel &n = read[3];
  if (n.emitted != "n" || !Near(n.sigmaMb, 35.25, 1e-6) ||
      n.residues[0] != std::make_pair(19, 40) ||
      n.steps[0].levels.size() != 2 || n.steps[0].bins.size() != 2 ||
      n.steps[0].bins[0][0] != 7.75 || n.steps[0].bins[0][1] != 8.25)
    return Fail("the (a,n) channel reads back wrong");
  const std::vector<std::string> broken = {
      Text({"channel n 1.0", "level 0 1.0"}),
      Text({"channel nn 1.0", "step 0 19 40", "level 0 1"}),
      Text({"channel n 1.0", "step 1 19 40", "level 0 1"}),
      Text({"channel x 1.0", "step 0 1 1"})};
  for (size_t i = 0; i < broken.size(); ++i)
    if (music::ParseChannelTable(broken[i], "broken", read, error))
      return Fail("a malformed channel table was accepted: " + broken[i]);
  return kTRUE;
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
  if (argc != 3) {
    std::cerr << "usage: talys_tables <fixtures-dir> <fake-cascade-dir>"
              << std::endl;
    return 2;
  }
  if (!CheckRutherford() || !CheckKeywords() || !CheckDeck() ||
      !CheckCascade(argv[2]) || !CheckElasticRows(argv[1]) ||
      !CheckElasticTable(argv[1]))
    return 1;
  std::cout << "talys_tables: all checks passed" << std::endl;
  return 0;
}
