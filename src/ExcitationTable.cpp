#include "ExcitationTable.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <fstream>
#include <sstream>

namespace {

// Parses one whole token as a finite number; "1.0x" and "nan" are rejected.
Bool_t ParseNumber(const std::string &token, Double_t &value) {
  size_t consumed = 0;
  try {
    value = std::stod(token, &consumed);
  } catch (const std::exception &) {
    return kFALSE;
  }
  return consumed == token.size() && std::isfinite(value);
}

// Weight of the part of `bin` that lies in [floor, EneAvail].
Double_t ClippedBinWeight(const std::array<Double_t, 3> &bin, Double_t floor,
                          Double_t EneAvail) {
  const Double_t lo = std::max(bin[0], floor);
  const Double_t hi = bin[1];
  if (lo >= EneAvail || lo >= hi)
    return 0.0;
  const Double_t top = std::min(hi, EneAvail);
  return bin[2] * (top - lo) / (hi - bin[0]);
}

} // namespace

namespace music {

Bool_t ExcitationTable::Load(const std::string &path, std::string &error) {
  std::ifstream file(path, std::ios::binary);
  if (!file) {
    levels.clear();
    bins.clear();
    total = 0.0;
    error = "cannot open level table '" + path + "'";
    return kFALSE;
  }
  std::ostringstream text;
  text << file.rdbuf();
  return Parse(text.str(), path, error);
}

Bool_t ExcitationTable::Parse(const std::string &text,
                              const std::string &source, std::string &error) {
  levels.clear();
  bins.clear();
  total = 0.0;
  std::istringstream input(text);

  std::string line;
  Int_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::string where = source + ":" + std::to_string(lineNumber);
    const size_t comment = line.find('#');
    if (comment != std::string::npos)
      line.erase(comment);
    std::istringstream fields(line);
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token)
      tokens.push_back(token);
    if (tokens.empty())
      continue;

    const std::string &tag = tokens[0];
    std::vector<Double_t> values(tokens.size() - 1, 0.0);
    for (size_t i = 1; i < tokens.size(); ++i) {
      if (!ParseNumber(tokens[i], values[i - 1])) {
        error = where + ": '" + tokens[i] + "' is not a finite number";
        levels.clear();
        bins.clear();
        return kFALSE;
      }
    }

    std::string problem;
    if (tag == "ground" || tag == "level") {
      if (values.size() != 2)
        problem = "'" + tag + "' rows take E_ex and weight";
      else if (values[0] < 0.0)
        problem = "E_ex must be >= 0";
      else if (tag == "ground" && values[0] != 0.0)
        problem = "a 'ground' row must have E_ex = 0";
      else if (values[1] < 0.0)
        problem = "weight must be >= 0";
      else
        levels.emplace_back(values[0], values[1]);
    } else if (tag == "bin") {
      if (values.size() != 3)
        problem = "'bin' rows take lo, hi, and weight";
      else if (values[0] < 0.0)
        problem = "bin lo must be >= 0";
      else if (values[1] <= values[0])
        problem = "bin hi must be > lo";
      else if (values[2] < 0.0)
        problem = "weight must be >= 0";
      else
        bins.push_back({values[0], values[1], values[2]});
    } else {
      problem = "unknown row tag '" + tag + "'";
    }
    if (!problem.empty()) {
      error = where + ": " + problem;
      levels.clear();
      bins.clear();
      return kFALSE;
    }
  }

  Double_t sum = 0.0;
  for (size_t i = 0; i < levels.size(); ++i)
    sum += levels[i].second;
  for (size_t i = 0; i < bins.size(); ++i)
    sum += bins[i][2];
  if (!(sum > 0.0)) {
    error = "level table '" + source + "' has no positive weight";
    levels.clear();
    bins.clear();
    return kFALSE;
  }
  for (size_t i = 0; i < levels.size(); ++i)
    levels[i].second /= sum;
  for (size_t i = 0; i < bins.size(); ++i)
    bins[i][2] /= sum;
  total = 1.0;
  return kTRUE;
}

Double_t ExcitationTable::AllowedWeight(Double_t EneAvail) const {
  return AllowedWeight(0.0, EneAvail);
}

Double_t ExcitationTable::AllowedWeight(Double_t floor,
                                        Double_t EneAvail) const {
  Double_t allowed = 0.0;
  for (size_t i = 0; i < levels.size(); ++i)
    if (levels[i].first >= floor && levels[i].first <= EneAvail)
      allowed += levels[i].second;
  for (size_t i = 0; i < bins.size(); ++i)
    allowed += ClippedBinWeight(bins[i], floor, EneAvail);
  return allowed;
}

Double_t ExcitationTable::Sample(Double_t EneAvail, TRandom *rng,
                                 Bool_t &ok) const {
  return Sample(0.0, EneAvail, rng, ok);
}

Double_t ExcitationTable::Sample(Double_t floor, Double_t EneAvail,
                                 TRandom *rng, Bool_t &ok) const {
  const Double_t allowed = AllowedWeight(floor, EneAvail);
  ok = allowed > 0.0;
  if (!ok)
    return 0.0;

  // Walk levels, then bins, as one CDF over the allowed weight.
  Double_t u = rng->Uniform(0.0, allowed);
  Double_t lastAllowed = 0.0;
  for (size_t i = 0; i < levels.size(); ++i) {
    const Double_t energy = levels[i].first;
    const Double_t weight = levels[i].second;
    if (energy < floor || energy > EneAvail || weight <= 0.0)
      continue;
    if (u < weight)
      return energy;
    u -= weight;
    lastAllowed = energy;
  }
  for (size_t i = 0; i < bins.size(); ++i) {
    const Double_t weight = ClippedBinWeight(bins[i], floor, EneAvail);
    if (weight <= 0.0)
      continue;
    const Double_t lo = std::max(bins[i][0], floor);
    const Double_t top = std::min(bins[i][1], EneAvail);
    if (u < weight)
      return lo + (top - lo) * (u / weight);
    u -= weight;
    lastAllowed = top;
  }
  // Rounding can leave u at the very top of the CDF.
  return lastAllowed;
}

Bool_t ExcitationTable::LowestAbove(Double_t floor, Double_t &Ex) const {
  Bool_t found = kFALSE;
  for (size_t i = 0; i < levels.size(); ++i)
    if (levels[i].second > 0.0 && levels[i].first >= floor &&
        (!found || levels[i].first < Ex)) {
      Ex = levels[i].first;
      found = kTRUE;
    }
  for (size_t i = 0; i < bins.size(); ++i)
    if (bins[i][2] > 0.0 && bins[i][1] > floor &&
        (!found || std::max(bins[i][0], floor) < Ex)) {
      Ex = std::max(bins[i][0], floor);
      found = kTRUE;
    }
  return found;
}

} // namespace music
