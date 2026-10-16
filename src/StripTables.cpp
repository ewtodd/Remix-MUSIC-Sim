#include "Simulator.hpp"
#include "TalysCascade.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <set>

// The z range of a readout strip, as the event loop samples vertices in it.
void Simulator::StripZRange(Int_t StpID, Double_t &MinZ, Double_t &MaxZ) const {
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (Int_t stp = 1; stp < AnodeRows; stp++) {
    MinZ += AnodeDZ[stp - 1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0] == StpID)
      break;
  }
}

// Mean beam KE after `length` cm of gas, integrated in sim_step slices as the
// transport steps it; a single GetFinalEnergy call is one linear step.
Double_t Simulator::MeanBeamEnergyAfter(Double_t energy, Double_t length) {
  const Long64_t steps = static_cast<Long64_t>(std::ceil(length / ctf.SimStep));
  const Double_t dx = steps > 0 ? length / static_cast<Double_t>(steps) : 0.0;
  for (Long64_t i = 0; i < steps && energy > 0.0; ++i)
    energy = Beam->GetFinalEnergy(0, energy, dx);
  return energy;
}

Bool_t Simulator::PlanTalysTables(std::vector<music::TalysRequest> &requests) {
  DisableVisualization();
  planTalysOnly_ = kTRUE;
  const Bool_t ok = SetupRun();
  planTalysOnly_ = kFALSE;
  requests.clear();
  for (size_t strip = 0; strip < talysRequests_.size(); ++strip)
    if (talysRequests_[strip].strip >= 0)
      requests.push_back(talysRequests_[strip]);
  return ok;
}

namespace {

// The (Z, A) a TALYS emission letter stands for; a gamma is (0, 0).
std::pair<Int_t, Int_t> EmittedNucleus(char letter) {
  switch (letter) {
  case 'n':
    return std::make_pair(0, 1);
  case 'p':
    return std::make_pair(1, 1);
  case 'd':
    return std::make_pair(1, 2);
  case 't':
    return std::make_pair(1, 3);
  case 'h':
    return std::make_pair(2, 3);
  case 'a':
    return std::make_pair(2, 4);
  }
  return std::make_pair(0, 0);
}

// Legend order: fewer emissions first, then by name.
bool ChannelBefore(const std::string &a, const std::string &b) {
  return a.size() != b.size() ? a.size() < b.size() : a < b;
}

std::string StripWhere(const music::TalysRequest &request) {
  char where[96];
  std::snprintf(where, sizeof(where), "strip %d (E_cm = %.4g MeV)",
                request.strip, request.ecm);
  return where;
}

} // namespace

// A talys_class run's TALYS inputs, one per strip at its centre energy;
// unless planning, reads the tables and builds the class's channels.
Bool_t Simulator::PrepareStripTables() {
  talysRequests_.clear();
  warnings_.clear();
  channels_.clear();
  if (ctf.talysClass.IsNull() || ctf.stripFirst < 0)
    return kTRUE;

  const Int_t beamZ = NuF->GetZ(ctf.beamName.Data());
  const Int_t beamA = NuF->GetA(ctf.beamName.Data());
  const std::string potential =
      music::AtomkiV2Potential(music::AtomkiV2Dir(), beamZ, beamA);
  Int_t alphaOmp = ctf.talys.alphaOmp;
  if (alphaOmp == 0)
    alphaOmp = potential.empty() ? 2 : 9;
  std::string potentialMd5;
  if (alphaOmp == 9) {
    std::string bytes;
    if (potential.empty() || !music::ReadFileText(potential, bytes)) {
      std::cerr << "musicsim ERROR: talys.alpha_omp = 9 needs an Atomki-V2 "
                   "potential for "
                << ctf.beamName << " in TALYS_ATOMKI_V2_DIR ('"
                << music::AtomkiV2Dir() << "')." << std::endl;
      return kFALSE;
    }
    potentialMd5 = music::Md5Hex(bytes);
  }

  const Double_t mb = Beam->Mass;
  const Double_t mt = Target->Mass;
  talysRequests_.assign(kReadoutStrips, music::TalysRequest());
  for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
    Double_t MinZ = 0;
    Double_t MaxZ = 0;
    StripZRange(strip, MinZ, MaxZ);
    music::TalysRequest &request = talysRequests_[strip];
    request.strip = strip;
    request.Z = beamZ;
    request.A = beamA;
    request.beamEnergy = MeanBeamEnergyAfter(Kb_at_gas, 0.5 * (MinZ + MaxZ));
    request.ecm =
        std::sqrt((mb + mt) * (mb + mt) + 2.0 * mt * request.beamEnergy) - mb -
        mt;
    request.stopped = !(request.beamEnergy > 0.0);
    if (request.stopped)
      continue;
    // Normal kinematics: the alpha lab energy with this CM energy.
    const Double_t alphaLab = request.ecm * (mt + mb) / mb;
    request.deck = music::TalysDeck(beamZ, beamA, alphaLab, alphaOmp,
                                    potentialMd5, ctf.talys.keywords);
    request.potential = alphaOmp == 9 ? potential : "";
    char name[192];
    std::snprintf(name, sizeof(name), "/talys/talys_%s_%s_Ecm%.4f_%s.tab",
                  ctf.beamName.Data(), ctf.target.Data(), request.ecm,
                  music::Md5Hex(request.deck).substr(0, 12).c_str());
    request.path = physics_.srimDir + name;
  }
  if (planTalysOnly_)
    return kTRUE;
  if (!SelectChannels() || !BuildChannelParticles())
    return kFALSE;
  for (size_t w = 0; workerId_ == 0 && w < warnings_.size(); ++w)
    std::cerr << "musicsim WARNING: " << warnings_[w] << std::endl;
  return kTRUE;
}

// Reads each strip's table and keeps the channels ending in the class element
// that reach min_channel_fraction of the class and min_channel_mb at one strip,
// plus elastic when the class is the beam's element.
Bool_t Simulator::SelectChannels() {
  channels_.clear();
  const Int_t beamZ = NuF->GetZ(ctf.beamName.Data());
  const Int_t beamA = NuF->GetA(ctf.beamName.Data());
  const Int_t targetZ = NuF->GetZ(ctf.target.Data());
  const Int_t targetA = NuF->GetA(ctf.target.Data());
  const Bool_t elastic = ctf.classZ == beamZ;
  const Double_t fraction = ctf.talys.minChannelFraction;
  const Double_t minMb = ctf.talys.minChannelMb;
  const std::string element(ctf.talysClass.Data());
  std::vector<std::vector<music::TalysChannel>> read(kReadoutStrips);
  std::vector<Double_t> classMb(kReadoutStrips, 0.0);
  std::map<std::string, Bool_t> kept;
  stripElastic_.assign(kReadoutStrips, music::ElasticTable());
  stripDroppedMb_.assign(kReadoutStrips, 0.0);
  for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
    const music::TalysRequest &request = talysRequests_[strip];
    const std::string where = StripWhere(request);
    if (request.stopped) {
      warnings_.push_back(
          "strip " + std::to_string(strip) +
          ": the mean beam stops before the strip centre, so it has no "
          "TALYS table; " +
          (elastic ? std::string("any reaction there is elastic, with pure "
                                 "Rutherford angles")
                   : "no channel ending in " + element +
                         " is simulated there, so its events do not react"));
      continue;
    }
    std::string text;
    Bool_t found = kFALSE;
    if (masterTables_) {
      const std::map<std::string, std::string>::const_iterator table =
          masterTables_->find(request.path);
      found = table != masterTables_->end();
      if (found)
        text = table->second;
    } else {
      found = music::ReadFileText(request.path, text);
    }
    if (!found) {
      std::cerr << "musicsim ERROR: missing TALYS table '" << request.path
                << "' for " << where << "; run: talys-cache " << ctrlFilePath_
                << std::endl;
      return kFALSE;
    }
    if (text.find(music::TalysDeckComment(request.deck)) == std::string::npos) {
      std::cerr << "musicsim ERROR: TALYS table '" << request.path
                << "' was not made from strip " << strip
                << "'s TALYS input; delete it and rerun talys-cache."
                << std::endl;
      return kFALSE;
    }
    std::string error;
    if (!music::ParseChannelTable(text, request.path, read[strip], error) ||
        (elastic && !stripElastic_[strip].Parse(text, request.path, error))) {
      std::cerr << "musicsim ERROR: " << error << std::endl;
      return kFALSE;
    }
    const std::vector<music::TalysChannel> &channels = read[strip];
    for (size_t ch = 0; ch < channels.size(); ++ch)
      if (channels[ch].residues.back().first == ctf.classZ)
        classMb[strip] += channels[ch].sigmaMb;
    for (size_t ch = 0; ch < channels.size(); ++ch)
      if (channels[ch].residues.back().first == ctf.classZ)
        kept[channels[ch].emitted] =
            kept[channels[ch].emitted] ||
            (channels[ch].sigmaMb >= fraction * classMb[strip] &&
             channels[ch].sigmaMb >= minMb);
    if (!(classMb[strip] > 0.0))
      warnings_.push_back(
          where + ": TALYS has no channel ending in " + element + " open, " +
          (elastic ? "so every reaction in this strip is elastic"
                   : "so its events do not react"));
    if (elastic && stripElastic_[strip].theta.empty())
      warnings_.push_back(where + ": TALYS gave no elastic angular "
                                  "distribution, so elastic is pure "
                                  "Rutherford there");
    levelTableText_[request.path] = text;
  }

  if (elastic) {
    ReactionChannel channel;
    channel.evap.push_back(std::make_pair(targetZ, targetA));
    channel.res.push_back(std::make_pair(beamZ, beamA));
    channel.stripSigma.assign(kReadoutStrips, 0.0);
    for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
      const music::TalysRequest &request = talysRequests_[strip];
      // A stopped strip has nothing else, so its weight is only nonzero.
      channel.stripSigma[strip] =
          request.stopped
              ? 1.0
              : stripElastic_[strip].AboveMb(beamZ, targetZ, request.ecm,
                                             ctf.classThetaCmMinDeg);
    }
    channels_.push_back(channel);
  }
  std::vector<std::string> names;
  for (std::map<std::string, Bool_t>::const_iterator it = kept.begin();
       it != kept.end(); ++it)
    if (it->second)
      names.push_back(it->first);
  std::sort(names.begin(), names.end(), ChannelBefore);
  for (size_t n = 0; n < names.size(); ++n) {
    ReactionChannel channel;
    channel.emitted = names[n];
    channel.stripSigma.assign(kReadoutStrips, 0.0);
    channel.stripSteps.assign(kReadoutStrips,
                              std::vector<music::ExcitationTable>());
    for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
      for (size_t ch = 0; ch < read[strip].size(); ++ch) {
        const music::TalysChannel &talys = read[strip][ch];
        if (talys.emitted != channel.emitted)
          continue;
        channel.stripSigma[strip] = talys.sigmaMb;
        channel.stripSteps[strip] = talys.steps;
        if (!channel.res.empty())
          continue;
        for (size_t k = 0; k < talys.residues.size(); ++k) {
          channel.evap.push_back(EmittedNucleus(talys.emitted[k]));
          channel.res.push_back(talys.residues[k]);
        }
      }
    }
    channels_.push_back(channel);
  }
  for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
    for (size_t ch = 0; ch < read[strip].size(); ++ch)
      if (read[strip][ch].residues.back().first == ctf.classZ &&
          !kept[read[strip][ch].emitted])
        stripDroppedMb_[strip] += read[strip][ch].sigmaMb;
    // Warned about when, taken together, they would have passed the thresholds.
    if (!(stripDroppedMb_[strip] > fraction * classMb[strip] &&
          stripDroppedMb_[strip] >= minMb))
      continue;
    char sums[160];
    std::snprintf(sums, sizeof(sums),
                  "%.3g mb, more than %.3g%% of the class's %.3g mb",
                  stripDroppedMb_[strip], 100.0 * fraction, classMb[strip]);
    warnings_.push_back(StripWhere(talysRequests_[strip]) +
                        ": the channels ending in " + element +
                        " left out by the thresholds add up to " + sums);
  }
  if (channels_.empty()) {
    std::cerr << "musicsim ERROR: no TALYS channel ending in " << element
              << " reaches min_channel_fraction and min_channel_mb at any "
                 "selected strip."
              << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Particle *Simulator::NewProduct(const std::pair<Int_t, Int_t> &nucleus,
                                Float_t dEdxScale) {
  Particle *particle = nullptr;
  if (nucleus.second == 0) {
    particle = new Particle("gamma", 0.0, 0, /*SaveTrajectory=*/false);
  } else {
    const std::string name = NuF->GetName(nucleus.first, nucleus.second);
    particle = new Particle(
        name, NuF->GetMass(nucleus.first, nucleus.second, "MeV/c^2"),
        nucleus.first, /*SaveTrajectory=*/false);
    particle->SetMedium(&gas_, physics_, dEdxScale);
  }
  ownedParticles_.push_back(particle);
  return particle;
}

// Gives each class channel its own particles per step and, per strip, the
// lowest excitation from which each step's chain can finish; the step slots
// start on idle particles.
Bool_t Simulator::BuildChannelParticles() {
  Int_t longest = 0;
  for (size_t ch = 0; ch < channels_.size(); ++ch) {
    for (size_t k = 0; k < channels_[ch].res.size(); ++k)
      if (NuF->GetName(channels_[ch].res[k].first, channels_[ch].res[k].second)
              .empty()) {
        std::cerr << "musicsim ERROR: TALYS channel '" << channels_[ch].emitted
                  << "' leaves a nuclide (Z = " << channels_[ch].res[k].first
                  << ", A = " << channels_[ch].res[k].second
                  << ") that is not in the mass table." << std::endl;
        return kFALSE;
      }
    longest = std::max(longest, static_cast<Int_t>(channels_[ch].res.size()));
  }
  if (longest > maxEvaporations) {
    std::cerr << "musicsim ERROR: a TALYS channel has more than "
              << maxEvaporations << " steps." << std::endl;
    return kFALSE;
  }
  for (size_t ch = 0; ch < channels_.size(); ++ch) {
    ReactionChannel &channel = channels_[ch];
    const Int_t steps = static_cast<Int_t>(channel.res.size());
    for (Int_t k = 0; k < steps; ++k) {
      channel.evapParticle.push_back(
          NewProduct(channel.evap[k], ctf.classEvapDedxScale));
      channel.resParticle.push_back(
          NewProduct(channel.res[k], ctf.classResidueDedxScale));
    }
    if (channel.emitted.empty())
      continue;
    channel.stripFloor.assign(kReadoutStrips,
                              std::vector<Double_t>(steps, 0.0));
    channel.stripOpen.assign(kReadoutStrips, 0.0);
    for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
      if (!(channel.stripSigma[strip] > 0.0))
        continue;
      const std::vector<music::ExcitationTable> &tables =
          channel.stripSteps[strip];
      std::vector<Double_t> &floor = channel.stripFloor[strip];
      Bool_t reachable = kTRUE;
      for (Int_t k = steps - 2; reachable && k >= 0; --k) {
        // Step k + 1 needs its separation energy plus its lowest state.
        const Double_t separation = channel.evapParticle[k + 1]->Mass +
                                    channel.resParticle[k + 1]->Mass -
                                    channel.resParticle[k]->Mass;
        Double_t lowest = 0.0;
        reachable = tables[k + 1].LowestAbove(floor[k + 1], lowest);
        floor[k] = std::max(0.0, separation + lowest);
      }
      if (reachable)
        channel.stripOpen[strip] = tables[0].AllowedWeight(
            floor[0], std::numeric_limits<Double_t>::max());
      if (!(channel.stripOpen[strip] > 0.0)) {
        warnings_.push_back(StripWhere(talysRequests_[strip]) +
                            ": TALYS channel '" + channel.emitted +
                            "' has no state from which its chain can "
                            "finish, so it is left out there");
        channel.stripSigma[strip] = 0.0;
      }
    }
  }
  for (Int_t k = 0; k < longest; ++k) {
    EvaP[k] = new Particle("idle", 0.0, 0, /*SaveTrajectory=*/false);
    EvaR[k] = new Particle("idle", 0.0, 0, /*SaveTrajectory=*/false);
    ownedParticles_.push_back(EvaP[k]);
    ownedParticles_.push_back(EvaR[k]);
  }
  idleSlots_.clear();
  for (Int_t k = 0; k < longest; ++k)
    idleSlots_.push_back(std::make_pair(EvaP[k], EvaR[k]));
  numEvaporations = longest;
  channelOpen_.assign(channels_.size(), 0.0);
  return kTRUE;
}

// Draws by cross section x the open fraction of step 0's table at this CM
// energy, as SelectBranch does for weighted branches.
Int_t Simulator::SelectChannel(Double_t sqrtS) {
  Double_t total = 0.0;
  for (size_t ch = 0; ch < channels_.size(); ++ch) {
    const ReactionChannel &channel = channels_[ch];
    const Double_t sigma = channel.stripSigma[simStrip_];
    channelOpen_[ch] = 0.0;
    if (!(sigma > 0.0))
      continue;
    if (channel.emitted.empty()) {
      channelOpen_[ch] = sigma;
    } else {
      const Double_t eneAvail =
          sqrtS - channel.evapParticle[0]->Mass - channel.resParticle[0]->Mass;
      channelOpen_[ch] = sigma *
                         channel.stripSteps[simStrip_][0].AllowedWeight(
                             channel.stripFloor[simStrip_][0], eneAvail) /
                         channel.stripOpen[simStrip_];
    }
    total += channelOpen_[ch];
  }
  if (!(total > 0.0))
    return -1;
  Double_t pick = Rdm->Uniform(0.0, total);
  Int_t lastOpen = -1;
  for (size_t ch = 0; ch < channels_.size(); ++ch) {
    if (channelOpen_[ch] <= 0.0)
      continue;
    if (pick < channelOpen_[ch])
      return static_cast<Int_t>(ch);
    pick -= channelOpen_[ch];
    lastOpen = static_cast<Int_t>(ch);
  }
  // Rounding can leave pick at the very top of the cumulative weight.
  return lastOpen;
}

Bool_t Simulator::ChannelProducts(std::vector<std::string> &ions) {
  ions.clear();
  if (ctf.talysClass.IsNull())
    return kTRUE;
  if (!SelectChannels())
    return kFALSE;
  std::set<std::string> names;
  for (size_t ch = 0; ch < channels_.size(); ++ch) {
    for (size_t k = 0; k < channels_[ch].res.size(); ++k) {
      if (channels_[ch].evap[k].first > 0)
        names.insert(NuF->GetName(channels_[ch].evap[k].first,
                                  channels_[ch].evap[k].second));
      names.insert(NuF->GetName(channels_[ch].res[k].first,
                                channels_[ch].res[k].second));
    }
  }
  ions.assign(names.begin(), names.end());
  return kTRUE;
}

Bool_t Simulator::SrimTableSpec(std::string &dir, std::string &gas,
                                Double_t &pressure,
                                Double_t &temperature) const {
  dir = physics_.srimDir;
  gas = ctf.gas.Data();
  pressure = ctf.pressure;
  temperature = ctf.temperature;
  return physics_.stoppingModel != 0;
}
