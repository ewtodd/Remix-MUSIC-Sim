#include "TalysCascade.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <sstream>

namespace {

// Change in (Z, N) for each emission type, in kEmissionLetters order.
const Int_t kDZ[music::kEmissionTypes] = {0, 0, 1, 1, 1, 2, 2};
const Int_t kDN[music::kEmissionTypes] = {0, 1, 0, 1, 2, 1, 2};
const char *const kTypeNames[music::kEmissionTypes] = {
    "gamma", "neutron", "proton", "deuteron", "triton", "helium-3", "alpha"};

Int_t TypeIndex(const std::string &name) {
  for (Int_t t = 0; t < music::kEmissionTypes; ++t)
    if (name == kTypeNames[t])
      return t;
  return -1;
}

Int_t LetterIndex(char letter) {
  for (Int_t t = 0; t < music::kEmissionTypes; ++t)
    if (music::kEmissionLetters[t] == letter)
      return t;
  return -1;
}

Double_t ValueAfterColon(const std::string &line) {
  return std::strtod(line.c_str() + line.find_last_of(':') + 1, nullptr);
}

Bool_t StartsWith(const std::string &line, const char *prefix) {
  return line.rfind(prefix, 0) == 0;
}

// What each state goes on to emit, as ordered letters, with probabilities.
typedef std::map<std::string, Double_t> Fate;

const Fate &FateOf(const std::map<music::CascadeKey, Fate> &fates,
                   const music::CascadeKey &key) {
  static const Fate stays = {{"", 1.0}};
  const std::map<music::CascadeKey, Fate>::const_iterator found =
      fates.find(key);
  return found == fates.end() ? stays : found->second;
}

// Bottom-up: every destination sorts before its source (lower bin, or a
// nuclide with lower Z or N), so it is final when its source is reached.
std::map<music::CascadeKey, Fate> Fates(const music::Cascade &cascade) {
  std::map<music::CascadeKey, Fate> fates;
  for (std::map<music::CascadeKey, music::CascadeState>::const_iterator it =
           cascade.states.begin();
       it != cascade.states.end(); ++it) {
    const music::CascadeState &state = it->second;
    Fate fate;
    if (!(state.pop > 0.0) || state.flows.empty()) {
      fate[""] = 1.0;
    } else {
      for (size_t f = 0; f < state.flows.size(); ++f) {
        const music::CascadeFlow &flow = state.flows[f];
        const std::string emitted =
            flow.type == 0 ? ""
                           : std::string(1, music::kEmissionLetters[flow.type]);
        const Fate &rest = FateOf(fates, flow.dest);
        for (Fate::const_iterator r = rest.begin(); r != rest.end(); ++r)
          fate[emitted + r->first] += flow.mb / state.pop * r->second;
      }
    }
    fates[it->first] = fate;
  }
  return fates;
}

Double_t Probability(const std::map<music::CascadeKey, Fate> &fates,
                     const music::CascadeKey &key, const std::string &rest) {
  const Fate &fate = FateOf(fates, key);
  const Fate::const_iterator found = fate.find(rest);
  return found == fate.end() ? 0.0 : found->second;
}

std::string Sorted(std::string letters) {
  std::sort(letters.begin(), letters.end());
  return letters;
}

void AppendLine(std::string &text, const std::string &line) {
  std::ostringstream out;
  out << line << std::endl;
  text += out.str();
}

} // namespace

namespace music {

bool CascadeKey::operator<(const CascadeKey &other) const {
  if (Z != other.Z)
    return Z < other.Z;
  if (N != other.N)
    return N < other.N;
  if (bin != other.bin)
    return bin < other.bin;
  if (twoJ != other.twoJ)
    return twoJ < other.twoJ;
  return parity < other.parity;
}

Bool_t ReadCascade(std::istream &input, Int_t compoundZ, Int_t compoundN,
                   Cascade &cascade, std::string &error) {
  cascade = Cascade();
  cascade.compoundZ = compoundZ;
  cascade.compoundN = compoundN;
  std::map<CascadeKey, Double_t> population; // all feeding from outside
  std::string line, quantity;
  Bool_t inResidual = kFALSE;
  Int_t Z = -1, N = -1, A = -1, blockParity = 0, flowType = -1;
  Int_t stateBin = -1, stateTwoJ = 0, stateParity = 0, flowParity = 0;
  Long64_t remaining = 0; // rows left in the current block (its "entries")
  while (std::getline(input, line)) {
    if (StartsWith(line, "# header:")) {
      inResidual = kFALSE;
      quantity.clear();
      stateBin = -1;
      remaining = 0;
    } else if (StartsWith(line, "#   residual:")) {
      inResidual = kTRUE;
    } else if (inResidual && StartsWith(line, "#     Z:")) {
      Z = static_cast<Int_t>(ValueAfterColon(line));
    } else if (inResidual && StartsWith(line, "#     A:")) {
      A = static_cast<Int_t>(ValueAfterColon(line));
      N = A - Z;
      inResidual = kFALSE;
    } else if (StartsWith(line, "#     number of discrete levels (NL):")) {
      cascade.nuclides[std::make_pair(Z, N)].levels =
          static_cast<Int_t>(ValueAfterColon(line));
    } else if (StartsWith(line, "#     continuum bin size [MeV]:")) {
      cascade.nuclides[std::make_pair(Z, N)].binSize = ValueAfterColon(line);
    } else if (StartsWith(line, "#     entries:")) {
      remaining = static_cast<Long64_t>(ValueAfterColon(line));
    } else if (StartsWith(line, "#     type:")) {
      quantity = line.substr(line.find(':') + 1);
      quantity.erase(0, quantity.find_first_not_of(' '));
      quantity.erase(quantity.find_last_not_of(' ') + 1);
    } else if (StartsWith(line, "#     parity:")) {
      blockParity = stateParity = static_cast<Int_t>(ValueAfterColon(line));
    } else if (StartsWith(line, "#     bin:")) {
      stateBin = static_cast<Int_t>(ValueAfterColon(line));
    } else if (StartsWith(line, "#     excitation energy [MeV]:")) {
      cascade.nuclides[std::make_pair(Z, N)].binEx[stateBin] =
          ValueAfterColon(line);
    } else if (StartsWith(line, "#     spin:")) {
      stateTwoJ = static_cast<Int_t>(std::lround(2.0 * ValueAfterColon(line)));
    } else if (StartsWith(line, "#     population [mb]:") && stateBin >= 0) {
      cascade.states[CascadeKey{Z, N, stateBin, stateTwoJ, stateParity}].pop =
          ValueAfterColon(line);
    } else if (StartsWith(line, "#     ") &&
               line.find(" decay to P=") != std::string::npos) {
      std::istringstream fields(line.substr(6));
      std::string name;
      fields >> name;
      flowType = TypeIndex(name);
      flowParity = line.find("P=-1") != std::string::npos ? -1 : 1;
    } else if (!line.empty() && line[0] != '#' && remaining > 0) {
      // A block's rows; past its declared count, rows of a table TALYS
      // continues between states follow, which are not read.
      --remaining;
      std::istringstream fields(line);
      if (quantity == "population" && stateBin < 0) {
        // bin, Ex, population, population at this parity, then per spin
        // slot: J = slot, plus 1/2 in an odd-A nuclide.
        Double_t bin = 0, ex = 0, total = 0, partial = 0, value = 0;
        fields >> bin >> ex >> total >> partial;
        cascade.nuclides[std::make_pair(Z, N)].binEx[static_cast<Int_t>(bin)] =
            ex;
        for (Int_t slot = 0; fields >> value; ++slot)
          if (value > 0.0)
            population[CascadeKey{Z, N, static_cast<Int_t>(bin),
                                  2 * slot + A % 2, blockParity}] += value;
      } else if (quantity == "particle decay from this bin" && stateBin >= 0) {
        std::string name;
        Double_t value = 0;
        const Int_t type = (fields >> name >> value) ? TypeIndex(name) : -1;
        if (type >= 0)
          cascade.states[CascadeKey{Z, N, stateBin, stateTwoJ, stateParity}]
              .decay[type] = value;
      } else if (quantity.find("decay from J,P bin to residual J',P' bin") !=
                     std::string::npos &&
                 stateBin >= 0 && flowType >= 0) {
        // Columns are labelled in the source's spin convention; the slot is
        // the destination's J (plus 1/2 if its A is odd).
        const Int_t destZ = Z - kDZ[flowType], destN = N - kDN[flowType];
        Double_t bin = 0, ex = 0, value = 0;
        fields >> bin >> ex;
        CascadeState &state =
            cascade.states[CascadeKey{Z, N, stateBin, stateTwoJ, stateParity}];
        for (Int_t slot = 0; fields >> value; ++slot)
          if (value > 0.0)
            state.flows.push_back(CascadeFlow{
                flowType,
                CascadeKey{destZ, destN, static_cast<Int_t>(bin),
                           2 * slot + (destZ + destN) % 2, flowParity},
                value});
      }
    }
  }
  if (cascade.states.empty() && population.empty()) {
    error = "no residue populations in the TALYS population file";
    return kFALSE;
  }

  // The population blocks also hold feeding from other nuclei's decays;
  // what the compound's first emission fed is the rest.
  std::map<CascadeKey, Double_t> inflow;
  for (std::map<CascadeKey, CascadeState>::const_iterator it =
           cascade.states.begin();
       it != cascade.states.end(); ++it)
    for (size_t f = 0; f < it->second.flows.size(); ++f)
      if (it->second.flows[f].type != 0)
        inflow[it->second.flows[f].dest] += it->second.flows[f].mb;
  for (std::map<CascadeKey, Double_t>::const_iterator it = population.begin();
       it != population.end(); ++it) {
    const std::map<CascadeKey, Double_t>::const_iterator in =
        inflow.find(it->first);
    const Double_t fed = it->second - (in == inflow.end() ? 0.0 : in->second);
    // Both are printed to 7 digits, so a difference that small is zero.
    if (fed > 1e-5 * it->second)
      cascade.feeding[it->first] = fed;
  }
  return kTRUE;
}

std::vector<CascadeChannel> CascadeChannels(const Cascade &cascade) {
  const std::map<CascadeKey, Fate> fates = Fates(cascade);

  // Ordered channels and their cross sections, from what each first
  // emission fed and where that state goes next.
  std::map<std::string, Double_t> sigma;
  for (std::map<CascadeKey, Double_t>::const_iterator it =
           cascade.feeding.begin();
       it != cascade.feeding.end(); ++it) {
    for (Int_t t = 0; t < kEmissionTypes; ++t) {
      if (it->first.Z != cascade.compoundZ - kDZ[t] ||
          it->first.N != cascade.compoundN - kDN[t])
        continue;
      const Fate &fate = FateOf(fates, it->first);
      for (Fate::const_iterator r = fate.begin(); r != fate.end(); ++r)
        sigma[std::string(1, kEmissionLetters[t]) + r->first] +=
            it->second * r->second;
    }
  }

  std::vector<CascadeChannel> channels;
  for (std::map<std::string, Double_t>::const_iterator c = sigma.begin();
       c != sigma.end(); ++c) {
    if (!(c->second > 0.0))
      continue;
    CascadeChannel channel;
    channel.emitted = c->first;
    channel.sigmaMb = c->second;
    // Population that has made exactly the emissions so far, entering each
    // state of the current residue.
    Int_t first = LetterIndex(channel.emitted[0]);
    Int_t Z = cascade.compoundZ - kDZ[first],
          N = cascade.compoundN - kDN[first];
    std::map<CascadeKey, Double_t> entering;
    for (std::map<CascadeKey, Double_t>::const_iterator it =
             cascade.feeding.begin();
         it != cascade.feeding.end(); ++it)
      if (it->first.Z == Z && it->first.N == N)
        entering[it->first] = it->second;
    for (size_t step = 0; step < channel.emitted.size(); ++step) {
      const std::string rest = channel.emitted.substr(step + 1);
      std::map<Int_t, Double_t> table;
      for (std::map<CascadeKey, Double_t>::const_iterator it = entering.begin();
           it != entering.end(); ++it) {
        const Double_t mb = it->second * Probability(fates, it->first, rest);
        if (mb > 0.0)
          table[it->first.bin] += mb;
      }
      channel.residues.push_back(std::make_pair(Z, N));
      channel.steps.push_back(table);
      if (rest.empty())
        break;
      // Down this residue's gamma cascade, highest bins first, then through
      // the next emission into the next residue.
      const Int_t next = LetterIndex(rest[0]);
      std::map<CascadeKey, Double_t> following;
      std::map<CascadeKey, Double_t> present = entering;
      const CascadeKey low{Z, N, -1, 0, 0};
      const CascadeKey high{Z, N + 1, -1, 0, 0};
      std::map<CascadeKey, CascadeState>::const_iterator begin =
          cascade.states.lower_bound(low);
      std::map<CascadeKey, CascadeState>::const_iterator it =
          cascade.states.lower_bound(high);
      while (it != begin) {
        --it;
        const std::map<CascadeKey, Double_t>::const_iterator here =
            present.find(it->first);
        if (here == present.end() || !(it->second.pop > 0.0))
          continue;
        for (size_t f = 0; f < it->second.flows.size(); ++f) {
          const CascadeFlow &flow = it->second.flows[f];
          const Double_t mb = here->second * flow.mb / it->second.pop;
          if (flow.type == 0)
            present[flow.dest] += mb;
          else if (flow.type == next)
            following[flow.dest] += mb;
        }
      }
      entering.swap(following);
      Z -= kDZ[next];
      N -= kDN[next];
    }
    channels.push_back(channel);
  }
  return channels;
}

Bool_t ExclusiveCrossSections(const std::string &talysOut,
                              std::map<std::string, Double_t> &sigma,
                              Double_t &gammaFirst, std::string &error) {
  sigma.clear();
  gammaFirst = 0.0;
  const size_t start = talysOut.find("6a. Total exclusive cross sections");
  if (start == std::string::npos) {
    error = "no exclusive cross sections (6a) in the TALYS output";
    return kFALSE;
  }
  std::istringstream input(talysOut.substr(start));
  std::string line;
  while (std::getline(input, line)) {
    if (line.find("(n,gn) + (n,gp)") != std::string::npos) {
      gammaFirst = ValueAfterColon(line);
      return kTRUE;
    }
    // "   n   p   d   t   h   a    sigma  (a,xx)": six multiplicities first.
    std::istringstream fields(line);
    Int_t count[6];
    Double_t value = 0;
    std::string label;
    if (!(fields >> count[0] >> count[1] >> count[2] >> count[3] >> count[4] >>
          count[5] >> value >> label) ||
        label[0] != '(')
      continue;
    std::string key;
    for (Int_t t = 0; t < 6; ++t)
      key +=
          std::string(static_cast<size_t>(count[t]), kEmissionLetters[t + 1]);
    sigma[Sorted(key)] += value;
  }
  error = "the TALYS exclusive cross sections (6a) end unexpectedly";
  return kFALSE;
}

Bool_t CheckChannels(const std::vector<CascadeChannel> &channels,
                     const std::map<std::string, Double_t> &sigma,
                     Double_t gammaFirst, std::string &error) {
  std::map<std::string, Double_t> ours;
  Double_t oursGammaFirst = 0.0, total = gammaFirst;
  for (size_t c = 0; c < channels.size(); ++c) {
    const CascadeChannel &channel = channels[c];
    if (channel.emitted[0] == 'g' && channel.emitted.size() > 1)
      oursGammaFirst += channel.sigmaMb;
    else
      ours[Sorted(channel.emitted == "g" ? "" : channel.emitted)] +=
          channel.sigmaMb;
    for (size_t step = 0; step < channel.steps.size(); ++step) {
      Double_t sum = 0.0;
      for (std::map<Int_t, Double_t>::const_iterator b =
               channel.steps[step].begin();
           b != channel.steps[step].end(); ++b)
        sum += b->second;
      if (std::fabs(sum - channel.sigmaMb) > 1e-6 * channel.sigmaMb) {
        error = "channel " + channel.emitted + " step " + std::to_string(step) +
                " does not carry the channel's cross section";
        return kFALSE;
      }
    }
  }
  for (std::map<std::string, Double_t>::const_iterator s = sigma.begin();
       s != sigma.end(); ++s)
    total += s->second;
  std::map<std::string, Double_t> keys = ours;
  keys.insert(sigma.begin(), sigma.end());
  keys["(gamma first)"] = 0.0;
  for (std::map<std::string, Double_t>::const_iterator k = keys.begin();
       k != keys.end(); ++k) {
    // 6a prints 6 significant digits, the gamma-first total 5 decimals [mb];
    // 1e-9 of the total covers channels too small to print.
    Double_t mine = 0.0, theirs = 0.0, rounding = 0.0;
    if (k->first == "(gamma first)") {
      mine = oursGammaFirst;
      theirs = gammaFirst;
      rounding = 5e-6;
    } else {
      mine = ours.count(k->first) ? ours.at(k->first) : 0.0;
      theirs = sigma.count(k->first) ? sigma.at(k->first) : 0.0;
    }
    if (std::fabs(mine - theirs) >
        1e-3 * std::max(mine, theirs) + 1e-9 * total + rounding) {
      char numbers[200];
      std::snprintf(numbers, sizeof(numbers),
                    "the cascade gives %.6e mb for (a,%s), TALYS's exclusive "
                    "cross section is %.6e mb",
                    mine, k->first.c_str(), theirs);
      error = numbers;
      return kFALSE;
    }
  }
  return kTRUE;
}

std::string ChannelTableRows(const Cascade &cascade,
                             const std::vector<CascadeChannel> &channels) {
  std::ostringstream text;
  char row[160];
  for (size_t c = 0; c < channels.size(); ++c) {
    const CascadeChannel &channel = channels[c];
    std::snprintf(row, sizeof(row), "channel %s %.9e", channel.emitted.c_str(),
                  channel.sigmaMb);
    text << row << std::endl;
    for (size_t step = 0; step < channel.steps.size(); ++step) {
      const std::pair<Int_t, Int_t> &residue = channel.residues[step];
      text << "step " << step << " " << residue.first << " "
           << residue.first + residue.second << std::endl;
      const std::map<std::pair<Int_t, Int_t>, CascadeNuclide>::const_iterator
          nuclide = cascade.nuclides.find(residue);
      for (std::map<Int_t, Double_t>::const_iterator b =
               channel.steps[step].begin();
           b != channel.steps[step].end(); ++b) {
        if (!(b->second > 0.0))
          continue;
        const Double_t ex = nuclide == cascade.nuclides.end() ||
                                    !nuclide->second.binEx.count(b->first)
                                ? 0.0
                                : nuclide->second.binEx.at(b->first);
        if (nuclide == cascade.nuclides.end() ||
            b->first <= nuclide->second.levels)
          std::snprintf(row, sizeof(row), "level %.6f %.9e", ex, b->second);
        else
          std::snprintf(row, sizeof(row), "bin %.6f %.6f %.9e",
                        ex - 0.5 * nuclide->second.binSize,
                        ex + 0.5 * nuclide->second.binSize, b->second);
        text << row << std::endl;
      }
    }
  }
  return text.str();
}

Bool_t ParseChannelTable(const std::string &text, const std::string &source,
                         std::vector<TalysChannel> &channels,
                         std::string &error) {
  channels.clear();
  std::vector<std::vector<std::string>> rows; // per channel: each step's rows
  std::istringstream input(text);
  std::string line;
  Int_t lineNumber = 0;
  while (std::getline(input, line)) {
    ++lineNumber;
    const std::string where = source + ":" + std::to_string(lineNumber);
    std::istringstream fields(line);
    std::string tag;
    if (!(fields >> tag) || tag[0] == '#' || tag == "elastic")
      continue;
    if (tag == "channel") {
      TalysChannel channel;
      if (!(fields >> channel.emitted >> channel.sigmaMb) ||
          !std::isfinite(channel.sigmaMb) || channel.sigmaMb < 0.0 ||
          channel.emitted.find_first_not_of(kEmissionLetters) !=
              std::string::npos) {
        error = where + ": malformed 'channel' line";
        return kFALSE;
      }
      channels.push_back(channel);
      rows.push_back(std::vector<std::string>());
    } else if (tag == "step") {
      Int_t step = -1, Z = 0, A = 0;
      if (channels.empty() || !(fields >> step >> Z >> A) ||
          step != static_cast<Int_t>(rows.back().size())) {
        error = where + ": 'step' out of place or out of order";
        return kFALSE;
      }
      channels.back().residues.push_back(std::make_pair(Z, A));
      rows.back().push_back("");
    } else if (channels.empty() || rows.back().empty()) {
      error = where + ": '" + tag + "' row outside a channel step";
      return kFALSE;
    } else {
      AppendLine(rows.back().back(), line);
    }
  }
  for (size_t c = 0; c < channels.size(); ++c) {
    TalysChannel &channel = channels[c];
    if (rows[c].size() != channel.emitted.size()) {
      error = source + ": channel " + channel.emitted + " has " +
              std::to_string(rows[c].size()) + " steps";
      return kFALSE;
    }
    channel.steps.resize(rows[c].size());
    for (size_t step = 0; step < rows[c].size(); ++step)
      if (!channel.steps[step].Parse(rows[c][step],
                                     source + " (channel " + channel.emitted +
                                         ", step " + std::to_string(step) + ")",
                                     error))
        return kFALSE;
  }
  return kTRUE;
}

} // namespace music
