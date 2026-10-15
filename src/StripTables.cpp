#include "Simulator.hpp"

#include <cmath>
#include <cstdio>

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

// Step-0 talys-cache tables: one TALYS input per strip at its centre energy;
// unless planning, reads each table and sets the strip's weights from it.
Bool_t Simulator::PrepareStripTables() {
  talysRequests_.clear();
  warnings_.clear();
  if (!ctf.computedWeights[0] || ctf.stripFirst < 0)
    return kTRUE;

  const Int_t beamZ = NuF->GetZ(ctf.beamName.Data());
  const Int_t beamA = NuF->GetA(ctf.beamName.Data());
  const Int_t targetZ = NuF->GetZ(ctf.target.Data());
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

  std::vector<ExcitationBranch> &branches = ctf.branches[0];
  for (size_t b = 0; b < branches.size(); ++b) {
    branches[b].stripTable.assign(kReadoutStrips, music::ExcitationTable());
    branches[b].stripWeight.assign(kReadoutStrips, 0.0);
    branches[b].stripElastic.assign(kReadoutStrips, music::ElasticTable());
  }
  for (Int_t strip = ctf.stripFirst; strip <= ctf.stripLast; ++strip) {
    const music::TalysRequest &request = talysRequests_[strip];
    char where[96];
    std::snprintf(where, sizeof(where), "strip %d (E_cm = %.4g MeV)", strip,
                  request.ecm);
    if (request.stopped) {
      // No TALYS table: any beam that straggles this far scatters elastically.
      for (size_t b = 0; b < branches.size(); ++b)
        branches[b].stripWeight[strip] = branches[b].perStrip ? 0.0 : 1.0;
      warnings_.push_back("strip " + std::to_string(strip) +
                          ": the mean beam stops before the strip centre, so "
                          "it has no TALYS table; any reaction there is "
                          "elastic, with pure Rutherford angles");
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
    for (size_t b = 0; b < branches.size(); ++b) {
      ExcitationBranch &branch = branches[b];
      std::string error;
      const Bool_t parsed =
          branch.perStrip
              ? branch.stripTable[strip].Parse(text, request.path, error, kTRUE)
              : branch.stripElastic[strip].Parse(text, request.path, error);
      if (!parsed) {
        std::cerr << "musicsim ERROR: " << error << std::endl;
        return kFALSE;
      }
      branch.stripWeight[strip] =
          branch.perStrip
              ? branch.stripTable[strip].sum
              : branch.stripElastic[strip].AboveMb(beamZ, targetZ, request.ecm,
                                                   branch.thetaCmMinDeg);
      if (branch.perStrip && branch.stripWeight[strip] == 0.0)
        warnings_.push_back(std::string(where) +
                            ": TALYS has no (a,a') channel open (below the "
                            "first excited level), so every reaction in "
                            "this strip is elastic");
      if (!branch.perStrip && branch.stripElastic[strip].theta.empty())
        warnings_.push_back(std::string(where) +
                            ": TALYS gave no elastic angular distribution, "
                            "so elastic is pure Rutherford there");
    }
    levelTableText_[request.path] = text;
  }
  for (size_t w = 0; workerId_ == 0 && w < warnings_.size(); ++w)
    std::cerr << "musicsim WARNING: " << warnings_[w] << std::endl;
  return kTRUE;
}
