#include "ExcitationTable.hpp"

#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

Bool_t Fail(const std::string &message) {
  std::cerr << "excitation_table: " << message << std::endl;
  return kFALSE;
}

Bool_t Near(Double_t actual, Double_t expected, Double_t tolerance) {
  return std::fabs(actual - expected) <= tolerance;
}

// Loads the fixture and checks its parsed, normalized contents.
Bool_t CheckLoad(const std::string &fixture, music::ExcitationTable &table) {
  std::string error;
  if (!table.Load(fixture, error))
    return Fail("fixture did not load: " + error);
  if (table.levels.size() != 3 || table.bins.size() != 1)
    return Fail("fixture has the wrong number of levels or bins");
  if (!Near(table.total, 1.0, 1e-12) ||
      !Near(table.levels[0].second, 0.40, 1e-12) ||
      !Near(table.bins[0][2], 0.25, 1e-12))
    return Fail("fixture weights were not normalized");
  return kTRUE;
}

// Sample counts per level and in the bin, drawn at one EneAvail.
struct SampleStats {
  Long64_t atLevel[3] = {0, 0, 0};
  Long64_t inBin = 0;
  Double_t binSum = 0.0;
  Double_t maxEx = 0.0;
  Long64_t offSupport = 0;
};

SampleStats Draw(const music::ExcitationTable &table, Double_t EneAvail,
                 Long64_t count, UInt_t seed) {
  TRandom3 rng(seed);
  SampleStats stats;
  for (Long64_t i = 0; i < count; ++i) {
    Bool_t ok = kFALSE;
    const Double_t ex = table.Sample(EneAvail, &rng, ok);
    stats.maxEx = std::max(stats.maxEx, ex);
    Bool_t matched = kFALSE;
    for (size_t level = 0; level < table.levels.size(); ++level) {
      if (ex == table.levels[level].first) {
        ++stats.atLevel[level];
        matched = kTRUE;
      }
    }
    if (!matched && ex >= table.bins[0][0] && ex <= table.bins[0][1]) {
      ++stats.inBin;
      stats.binSum += ex;
      matched = kTRUE;
    }
    if (!matched || !ok)
      ++stats.offSupport;
  }
  return stats;
}

// Full range: the sampled shape follows the table weights.
Bool_t CheckFullRange(const music::ExcitationTable &table) {
  const Long64_t count = 200000;
  const SampleStats stats = Draw(table, 20.0, count, 11);
  const Double_t n = static_cast<Double_t>(count);
  if (stats.offSupport != 0)
    return Fail("a full-range sample fell off the table support");
  if (!Near(static_cast<Double_t>(stats.atLevel[0]) / n, 0.40, 0.01) ||
      !Near(static_cast<Double_t>(stats.atLevel[1]) / n, 0.25, 0.01) ||
      !Near(static_cast<Double_t>(stats.atLevel[2]) / n, 0.10, 0.01) ||
      !Near(static_cast<Double_t>(stats.inBin) / n, 0.25, 0.01))
    return Fail("full-range sample fractions do not follow the weights");
  const Double_t binMean = stats.binSum / static_cast<Double_t>(stats.inBin);
  if (!Near(binMean, 5.0, 0.02))
    return Fail("continuum-bin samples are not uniform over [4, 6]");
  return kTRUE;
}

// Over-EneAvail tail: the table is truncated at EneAvail and renormalized,
// so nothing lands above it and the allowed shape keeps its ratios.
Bool_t CheckTruncation(const music::ExcitationTable &table) {
  if (!Near(table.AllowedWeight(5.0), 0.875, 1e-12))
    return Fail("bin straddling EneAvail is not weighted pro rata");
  const Long64_t count = 200000;
  const SampleStats stats = Draw(table, 5.0, count, 12);
  const Double_t n = static_cast<Double_t>(count);
  if (stats.offSupport != 0 || stats.maxEx > 5.0)
    return Fail("a truncated sample exceeded EneAvail");
  if (!Near(static_cast<Double_t>(stats.atLevel[0]) / n, 0.40 / 0.875, 0.01) ||
      !Near(static_cast<Double_t>(stats.inBin) / n, 0.125 / 0.875, 0.01))
    return Fail("truncated sample was not renormalized over [0, EneAvail]");

  // A level exactly at EneAvail is allowed; one just above it is not.
  if (!Near(table.AllowedWeight(3.086), 0.65, 1e-12))
    return Fail("a level at exactly EneAvail was excluded");
  return kTRUE;
}

// Below the lowest state: no weight, so the step is rejected without a draw.
Bool_t CheckNoAllowedWeight(const music::ExcitationTable &table) {
  TRandom3 used(13);
  TRandom3 fresh(13);
  Bool_t ok = kTRUE;
  table.Sample(1.0, &used, ok);
  if (ok)
    return Fail("sampling below the lowest state did not report failure");
  if (used.Rndm() != fresh.Rndm())
    return Fail("a rejected sample consumed an RNG draw");
  return kTRUE;
}

// Each successful sample uses exactly one draw, so event streams stay aligned.
Bool_t CheckOneDraw(const music::ExcitationTable &table) {
  TRandom3 used(14);
  TRandom3 fresh(14);
  for (Int_t i = 0; i < 5; ++i) {
    Bool_t ok = kFALSE;
    table.Sample(20.0, &used, ok);
    fresh.Rndm();
  }
  if (used.Rndm() != fresh.Rndm())
    return Fail("a sample did not consume exactly one RNG draw");
  return kTRUE;
}

// A floor cuts the table from below as EneAvail does from above: a level at
// the floor is kept, a straddling bin counts pro rata, and floor 0 draws as
// the unfloored form does.
Bool_t CheckFloor(const music::ExcitationTable &table) {
  if (!Near(table.AllowedWeight(3.086, 20.0), 0.60, 1e-12) ||
      !Near(table.AllowedWeight(5.0, 20.0), 0.125, 1e-12) ||
      !Near(table.AllowedWeight(3.0, 5.0), 0.475, 1e-12))
    return Fail("floored allowed weights are wrong");
  Double_t lowest = -1.0;
  if (!table.LowestAbove(0.0, lowest) || lowest != 1.727 ||
      !table.LowestAbove(1.8, lowest) || lowest != 3.086 ||
      !table.LowestAbove(4.5, lowest) || lowest != 4.5 ||
      table.LowestAbove(6.0, lowest))
    return Fail("LowestAbove does not find the lowest state above a floor");
  TRandom3 floored(15);
  TRandom3 plain(15);
  TRandom3 rng(16);
  Double_t binSum = 0.0;
  const Int_t count = 100000;
  for (Int_t i = 0; i < count; ++i) {
    Bool_t okFloored = kFALSE, okPlain = kFALSE, ok = kFALSE;
    if (table.Sample(0.0, 5.0, &floored, okFloored) !=
            table.Sample(5.0, &plain, okPlain) ||
        okFloored != okPlain)
      return Fail("a zero floor changed a draw");
    const Double_t ex = table.Sample(5.0, 20.0, &rng, ok);
    if (!ok || ex < 5.0 || ex > 6.0)
      return Fail("a floored sample left [floor, the bin's top]");
    binSum += ex;
  }
  if (!Near(binSum / count, 5.5, 0.01))
    return Fail("floored bin samples are not uniform over [5, 6]");
  return kTRUE;
}

void WriteTable(const std::filesystem::path &path,
                const std::vector<std::string> &lines) {
  std::ofstream output(path);
  for (size_t i = 0; i < lines.size(); ++i)
    output << lines[i] << std::endl;
}

// Rows need not be sorted: a level above EneAvail listed first is never drawn.
Bool_t CheckUnsortedRows(const std::filesystem::path &directory) {
  const std::filesystem::path path = directory / "unsorted.tab";
  WriteTable(path, {"level 3.0 0.5", "level 1.0 0.5"});
  music::ExcitationTable table;
  std::string error;
  if (!table.Load(path.string(), error))
    return Fail("unsorted table did not load: " + error);
  TRandom3 rng(15);
  for (Int_t i = 0; i < 1000; ++i) {
    Bool_t ok = kFALSE;
    if (table.Sample(2.0, &rng, ok) != 1.0 || !ok)
      return Fail("an unsorted table returned a level above EneAvail");
  }
  return kTRUE;
}

// Malformed tables are hard errors that name the file and line.
Bool_t CheckRejects(const std::filesystem::path &directory) {
  const std::vector<std::vector<std::string>> badTables = {
      {"level 1.0 0.5", "foo 2.0 0.5"},
      {"bin 3.0 2.0 0.5"},
      {"level 1.0 -0.5"},
      {"level one 0.5"},
      {"ground 1.0 0.5"},
      {"level 1.0"},
      {"bin 2.0 2.0 0.5"},
      {"level 1.0 0.5 7"}};
  for (size_t i = 0; i < badTables.size(); ++i) {
    const std::filesystem::path path =
        directory / ("bad-" + std::to_string(i) + ".tab");
    WriteTable(path, badTables[i]);
    music::ExcitationTable table;
    std::string error;
    if (table.Load(path.string(), error))
      return Fail("accepted malformed table " + path.string());
    if (error.find(path.string() + ":") != 0)
      return Fail("error does not name the file and line: " + error);
    if (!table.levels.empty() || !table.bins.empty())
      return Fail("a rejected table kept partial rows");
  }

  const std::filesystem::path zero = directory / "zero-weight.tab";
  WriteTable(zero, {"level 1.0 0.0"});
  music::ExcitationTable table;
  std::string error;
  if (table.Load(zero.string(), error))
    return Fail("accepted a table with no positive weight");
  if (table.Load((directory / "missing.tab").string(), error))
    return Fail("accepted a missing table");
  return kTRUE;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: excitation_table FIXTURE TEMP_DIR" << std::endl;
    return 2;
  }
  music::ExcitationTable table;
  const Bool_t ok = CheckLoad(argv[1], table) && CheckFullRange(table) &&
                    CheckTruncation(table) && CheckNoAllowedWeight(table) &&
                    CheckOneDraw(table) && CheckFloor(table) &&
                    CheckUnsortedRows(argv[2]) && CheckRejects(argv[2]);
  return ok ? 0 : 1;
}
