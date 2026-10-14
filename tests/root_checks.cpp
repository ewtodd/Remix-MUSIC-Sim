#include "NuclideFinder.hpp"

#include <TFile.h>
#include <TFriendElement.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TList.h>
#include <TMath.h>
#include <TObjString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// Simulator::kAngleUnset: the angle sentinel ResetBranches writes.
constexpr Float_t kAngleUnset = -1000.0f;

struct RootFile {
  std::unique_ptr<TFile> file;
  TTree *events = nullptr;
  TTree *truth = nullptr;
};

bool Fail(const std::string &message) {
  std::cerr << "root_checks: " << message << '\n';
  return false;
}

RootFile Open(const char *path) {
  RootFile result;
  result.file.reset(TFile::Open(path, "READ"));
  if (!result.file || result.file->IsZombie()) {
    Fail(std::string("cannot open ") + path);
    return result;
  }
  result.events = dynamic_cast<TTree *>(result.file->Get("events_MeV"));
  result.truth = dynamic_cast<TTree *>(result.file->Get("MC"));
  if (!result.events || !result.truth)
    Fail(std::string("missing trees in ") + path);
  return result;
}

int KeyCycles(TFile *file, const char *name) {
  int count = 0;
  TIter next(file->GetListOfKeys());
  while (const auto *key = dynamic_cast<TKey *>(next()))
    if (std::string(key->GetName()) == name)
      ++count;
  return count;
}

bool Standard(const char *path, Long64_t expectedEntries) {
  RootFile root = Open(path);
  if (!root.events || !root.truth)
    return false;
  if (root.events->GetEntries() != expectedEntries ||
      root.truth->GetEntries() != expectedEntries)
    return Fail(std::string(path) + " has an unexpected entry count");
  if (KeyCycles(root.file.get(), "events_MeV") != 1 ||
      KeyCycles(root.file.get(), "MC") != 1)
    return Fail(std::string(path) + " has duplicate tree cycles");
  if (!root.events->GetListOfFriends() ||
      root.events->GetListOfFriends()->GetEntries() != 1)
    return Fail(std::string(path) + " is missing its MC friend");
  for (const char *branch :
       {"requested_strip", "reaction_strip", "event_index", "beam_termination",
        "evap_termination", "residue_termination"})
    if (!root.truth->GetBranch(branch))
      return Fail(std::string(path) + " is missing branch " + branch);
  for (const char *record :
       {"musicsim_version", "control_file_path", "control_file_toml",
        "resolved_configuration", "termination_codes", "straggling_sampler"})
    if (!root.file->Get((std::string("metadata/") + record).c_str()))
      return Fail(std::string(path) + " is missing metadata/" + record);
  const auto *resolved = dynamic_cast<TObjString *>(
      root.file->Get("metadata/resolved_configuration"));
  if (!resolved || !resolved->GetString().Contains("schema_version = 4\n"))
    return Fail(std::string(path) + " does not declare output schema 4");
  // Defaults change between builds, so a run is reproducible only if the
  // resolved values and the exact code revisions are recorded.
  const Int_t nResolvedKeys = 11;
  const char *resolvedKeys[nResolvedKeys] = {
      "git_revision = \"",       "catima_revision = \"",
      "catima_global = ",        "z_effective = \"",
      "low_energy = \"",         "straggling_z_effective = \"",
      "sim_step_cm = ",          "max_time_ns = ",
      "beam_energy_fwhm_mev = ", "total_entries = ",
      "musicsim_version = \""};
  for (Int_t i = 0; i < nResolvedKeys; ++i)
    if (!resolved->GetString().Contains(resolvedKeys[i]))
      return Fail(std::string(path) + " resolved configuration lacks " +
                  resolvedKeys[i]);
  if (resolved->GetString().Contains("git_revision = \"unknown\"") ||
      resolved->GetString().Contains("catima_revision = \"unknown\""))
    return Fail(std::string(path) + " does not record its code revisions");
  return true;
}

bool Conservation(const char *path) {
  RootFile root = Open(path);
  if (!root.events || !root.truth || root.events->GetEntries() != 1)
    return false;
  Float_t cathode = 0.0f;
  Float_t deadUpstream = 0.0f;
  Float_t deadDownstream = 0.0f;
  Float_t beamIn = 0.0f;
  Float_t beamOut = 0.0f;
  Int_t termination = 0;
  root.events->SetBranchAddress("Cathode", &cathode);
  root.truth->SetBranchAddress("DeadUS_dE", &deadUpstream);
  root.truth->SetBranchAddress("DeadDS_dE", &deadDownstream);
  root.truth->SetBranchAddress("beam_energy_gas", &beamIn);
  root.truth->SetBranchAddress("beam_energy_exit", &beamOut);
  root.truth->SetBranchAddress("beam_termination", &termination);
  root.events->GetEntry(0);
  root.truth->GetEntry(0);
  if (termination != 2 || beamOut != -1.0f)
    return Fail("stopping fixture did not stop in the gas");
  const double accounted = cathode + deadUpstream + deadDownstream;
  if (std::fabs(accounted - static_cast<double>(beamIn)) > 2e-5)
    return Fail(
        "stopping fixture lost kinetic energy: in=" + std::to_string(beamIn) +
        " deposited=" + std::to_string(accounted));
  return true;
}

// Every per-step product field and the residue state must still hold the
// ResetBranches sentinels. Reads entry 0 with its own branch addresses.
bool ProductTruthAtSentinels(TTree *truth, const std::string &label) {
  Int_t steps = 0;
  Float_t evapEnergy[10]{};
  Float_t residueEnergy[10]{};
  Float_t evapExitEnergy[10]{};
  Float_t residueExitEnergy[10]{};
  Float_t theta[10]{};
  Float_t phi[10]{};
  Float_t evapTheta[10]{};
  Float_t evapPhi[10]{};
  Float_t residueTheta[10]{};
  Float_t residuePhi[10]{};
  Int_t evapTermination[10]{};
  Int_t residueStep = 0;
  Int_t residueTermination = 0;
  truth->ResetBranchAddresses();
  truth->SetBranchAddress("n_steps", &steps);
  truth->SetBranchAddress("evap_energy", evapEnergy);
  truth->SetBranchAddress("residue_energy", residueEnergy);
  truth->SetBranchAddress("evap_energy_exit", evapExitEnergy);
  truth->SetBranchAddress("residue_energy_exit", residueExitEnergy);
  truth->SetBranchAddress("theta_cm", theta);
  truth->SetBranchAddress("phi_cm", phi);
  truth->SetBranchAddress("evap_theta", evapTheta);
  truth->SetBranchAddress("evap_phi", evapPhi);
  truth->SetBranchAddress("residue_theta", residueTheta);
  truth->SetBranchAddress("residue_phi", residuePhi);
  truth->SetBranchAddress("evap_termination", evapTermination);
  truth->SetBranchAddress("residue_step", &residueStep);
  truth->SetBranchAddress("residue_termination", &residueTermination);
  truth->GetEntry(0);
  // The locals above go out of scope; leave no branch pointing at them.
  truth->ResetBranchAddresses();
  for (Int_t step = 0; step < steps; ++step)
    if (evapEnergy[step] != -2.0f || residueEnergy[step] != -2.0f ||
        evapExitEnergy[step] != -2.0f || residueExitEnergy[step] != -2.0f ||
        theta[step] != kAngleUnset || phi[step] != kAngleUnset ||
        evapTheta[step] != kAngleUnset || evapPhi[step] != kAngleUnset ||
        residueTheta[step] != kAngleUnset || residuePhi[step] != kAngleUnset ||
        evapTermination[step] != 0)
      return Fail(label + " retained partial-chain truth");
  if (residueStep != -1 || residueTermination != 0)
    return Fail(label + " retained a residue transport state");
  return true;
}

bool Forbidden(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t requested = -2;
  Int_t reacted = -2;
  Int_t termination = 0;
  Float_t reactionEnergy = 0.0f;
  Float_t vertex = 0.0f;
  root.truth->SetBranchAddress("requested_strip", &requested);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("beam_termination", &termination);
  root.truth->SetBranchAddress("beam_energy_reaction", &reactionEnergy);
  root.truth->SetBranchAddress("vertex_z", &vertex);
  root.truth->GetEntry(0);
  if (requested != 0 || reacted != -1 || termination != 3 ||
      reactionEnergy != -2.0f || vertex != -1000.0f)
    return Fail("forbidden reaction left an inconsistent event history");
  return ProductTruthAtSentinels(root.truth, "forbidden reaction");
}

// Simulator::StripAtZ for the hardcoded anode in src/Geometry.cpp: a 3.590 cm
// upstream dead layer, then readout strips 0..17 of 1.578 cm each.
Int_t StripAtZ(Double_t z) {
  if (z < 0.0)
    return -1;
  Double_t rowEnd = 3.590;
  if (z < rowEnd)
    return -1;
  for (Int_t strip = 0; strip <= 17; ++strip) {
    rowEnd += 1.578;
    if (z < rowEnd)
      return strip;
  }
  return -1;
}

bool StoppedBeforeVertex(const char *path, Int_t requestedStrip) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t requested = -2;
  Int_t reacted = -2;
  Int_t termination = 0;
  Int_t stopStrip = -3;
  Float_t reactionEnergy = 0.0f;
  Float_t vertex = 0.0f;
  Float_t exitEnergy = 0.0f;
  Float_t stopZ = 0.0f;
  root.truth->SetBranchAddress("requested_strip", &requested);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("beam_termination", &termination);
  root.truth->SetBranchAddress("beam_stop_strip", &stopStrip);
  root.truth->SetBranchAddress("beam_energy_reaction", &reactionEnergy);
  root.truth->SetBranchAddress("vertex_z", &vertex);
  root.truth->SetBranchAddress("beam_energy_exit", &exitEnergy);
  root.truth->SetBranchAddress("beam_stop_z", &stopZ);
  root.truth->GetEntry(0);
  if (termination != 2)
    return Fail("stopped-before-vertex beam did not stop in the gas");
  if (requested != requestedStrip || reacted != -1 || reactionEnergy != -2.0f ||
      vertex != -1000.0f || exitEnergy != -1.0f)
    return Fail("stopped-before-vertex beam left reaction or exit truth");
  const Int_t expectedStrip = StripAtZ(static_cast<Double_t>(stopZ));
  if (stopStrip != expectedStrip || stopStrip < 0 ||
      stopStrip >= requestedStrip)
    return Fail("stopped-before-vertex beam_stop_strip " +
                std::to_string(stopStrip) + " does not match stop z " +
                std::to_string(stopZ) + " upstream of strip " +
                std::to_string(requestedStrip));
  return ProductTruthAtSentinels(root.truth, "stopped-before-vertex beam");
}

bool SideExit(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t reacted = -1;
  Int_t evapTermination[10]{};
  Int_t evapStopStrip[10]{};
  Float_t evapExitEnergy[10]{};
  Float_t evapStopX[10]{};
  Float_t evapStopY[10]{};
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("evap_termination", evapTermination);
  root.truth->SetBranchAddress("evap_stop_strip", evapStopStrip);
  root.truth->SetBranchAddress("evap_energy_exit", evapExitEnergy);
  root.truth->SetBranchAddress("evap_stop_x", evapStopX);
  root.truth->SetBranchAddress("evap_stop_y", evapStopY);
  root.truth->GetEntry(0);
  if (reacted < 0)
    return Fail("side-exit fixture did not react");
  if (evapTermination[0] != 5 || evapExitEnergy[0] != -1.0f ||
      evapStopStrip[0] != -1)
    return Fail("side-exit ejectile termination " +
                std::to_string(evapTermination[0]) +
                " did not record a side exit");
  // Half of AnodeLength (9 cm) and AnodeHeight (10 cm) from src/Geometry.cpp;
  // transport lands within 1e-10 cm of the face, far below this tolerance.
  const Double_t halfLength = 4.5;
  const Double_t halfHeight = 5.0;
  const Double_t tolerance = 1e-5;
  const Double_t x = std::fabs(static_cast<Double_t>(evapStopX[0]));
  const Double_t y = std::fabs(static_cast<Double_t>(evapStopY[0]));
  const Bool_t inside =
      x <= halfLength + tolerance && y <= halfHeight + tolerance;
  const Bool_t onSideFace = std::fabs(x - halfLength) <= tolerance ||
                            std::fabs(y - halfHeight) <= tolerance;
  if (!inside || !onSideFace)
    return Fail("side-exit ejectile stopped at x=" + std::to_string(x) +
                " y=" + std::to_string(y) + ", not on a side face");
  return true;
}

bool MultiStep(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t steps = 0;
  Int_t reacted = -1;
  Float_t theta[10]{};
  Float_t phi[10]{};
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("phi_cm", phi);
  root.truth->GetEntry(0);
  if (steps != 2 || reacted != 0)
    return Fail("multi-step fixture did not realize both reaction steps");
  for (int i = 0; i < steps; ++i)
    if (!std::isfinite(theta[i]) || theta[i] < 0.0f || theta[i] > 180.0f ||
        !std::isfinite(phi[i]) || phi[i] < -180.0f || phi[i] > 180.0f)
      return Fail("multi-step fixture stored an invalid CM angle");
  if (theta[0] == theta[1] && phi[0] == phi[1])
    return Fail("multi-step fixture repeated the final CM angle");
  return true;
}

struct FourMomentum {
  Double_t e = 0.0;
  Double_t px = 0.0;
  Double_t py = 0.0;
  Double_t pz = 0.0;
};

Double_t MomentumMagnitude(const FourMomentum &p) {
  return std::sqrt(p.px * p.px + p.py * p.py + p.pz * p.pz);
}

FourMomentum FromKineticEnergy(Double_t mass, Double_t kinetic,
                               Double_t thetaDeg, Double_t phiDeg) {
  const Double_t momentum = std::sqrt(kinetic * kinetic + 2.0 * mass * kinetic);
  const Double_t theta = thetaDeg * TMath::DegToRad();
  const Double_t phi = phiDeg * TMath::DegToRad();
  FourMomentum result;
  result.e = mass + kinetic;
  result.px = momentum * std::sin(theta) * std::cos(phi);
  result.py = momentum * std::sin(theta) * std::sin(phi);
  result.pz = momentum * std::cos(theta);
  return result;
}

// Same pure boost as FourVector::Boost: into the frame moving at beta.
FourMomentum ToFrameMovingAt(const FourMomentum &p, Double_t betaX,
                             Double_t betaY, Double_t betaZ) {
  const Double_t beta2 = betaX * betaX + betaY * betaY + betaZ * betaZ;
  if (beta2 <= 0.0)
    return p;
  const Double_t gamma = 1.0 / std::sqrt(1.0 - beta2);
  const Double_t betaDotP = betaX * p.px + betaY * p.py + betaZ * p.pz;
  const Double_t alongBeta = (gamma - 1.0) * betaDotP / beta2 - gamma * p.e;
  FourMomentum result;
  result.e = gamma * (p.e - betaDotP);
  result.px = p.px + alongBeta * betaX;
  result.py = p.py + alongBeta * betaY;
  result.pz = p.pz + alongBeta * betaZ;
  return result;
}

// Parent k = parent k-1 minus its ejectile, so no residue excitation is needed.
// In the parent frame the residue leaves opposite the ejectile, along theta_cm.
bool MultiStepReconstruction(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t steps = 0;
  Float_t reactionEnergy = 0.0f;
  Float_t evapEnergy[10]{};
  Float_t evapTheta[10]{};
  Float_t evapPhi[10]{};
  Float_t theta[10]{};
  Float_t phi[10]{};
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("beam_energy_reaction", &reactionEnergy);
  root.truth->SetBranchAddress("evap_energy", evapEnergy);
  root.truth->SetBranchAddress("evap_theta", evapTheta);
  root.truth->SetBranchAddress("evap_phi", evapPhi);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("phi_cm", phi);
  root.truth->GetEntry(0);
  if (steps != 2 || reactionEnergy <= 0.0f)
    return Fail("multi-step fixture has no realized two-step reaction");
  for (Int_t step = 0; step < steps; ++step)
    if (theta[step] == kAngleUnset || phi[step] == kAngleUnset ||
        evapTheta[step] == kAngleUnset || evapPhi[step] == kAngleUnset)
      return Fail("multi-step fixture left an angle unset");

  // Nuclides of tests/controls/multistep.toml, with the simulator's masses.
  NuclideFinder finder;
  const Double_t beamMass = finder.GetMass("98Zr", "MeV/c^2");
  const Double_t targetMass = finder.GetMass("4He", "MeV/c^2");
  const Double_t neutronMass = finder.GetMass("1n", "MeV/c^2");

  // Float_t truth is rounded to a relative epsilon. Each lab angle is at most
  // 2*pi in radians, so it carries an absolute error of epsilon * pi.
  const Double_t epsilon =
      static_cast<Double_t>(std::numeric_limits<Float_t>::epsilon());
  const Double_t beamKinetic = static_cast<Double_t>(reactionEnergy);
  const Double_t beamMomentum =
      std::sqrt(beamKinetic * beamKinetic + 2.0 * beamMass * beamKinetic);
  FourMomentum parent;
  parent.e = beamMass + beamKinetic + targetMass;
  parent.pz = beamMomentum;
  Double_t parentError =
      epsilon * beamKinetic * (beamMass + beamKinetic) / beamMomentum;

  for (Int_t step = 0; step < steps; ++step) {
    const std::string label = "multi-step step " + std::to_string(step);
    const Double_t betaX = parent.px / parent.e;
    const Double_t betaY = parent.py / parent.e;
    const Double_t betaZ = parent.pz / parent.e;
    const Double_t beta =
        std::sqrt(betaX * betaX + betaY * betaY + betaZ * betaZ);
    // Without a real boost the check could not tell lab from CM angles.
    if (beta < 0.05)
      return Fail(label + " parent frame is not appreciably boosted");
    const Double_t gamma = 1.0 / std::sqrt(1.0 - beta * beta);

    const Double_t evapKinetic = static_cast<Double_t>(evapEnergy[step]);
    const FourMomentum evapLab = FromKineticEnergy(
        neutronMass, evapKinetic, static_cast<Double_t>(evapTheta[step]),
        static_cast<Double_t>(evapPhi[step]));
    const FourMomentum evapRest = ToFrameMovingAt(evapLab, betaX, betaY, betaZ);
    const Double_t restMomentum = MomentumMagnitude(evapRest);
    if (restMomentum < 1e-3)
      return Fail(label + " has a degenerate CM momentum");

    const Double_t residueX = -evapRest.px / restMomentum;
    const Double_t residueY = -evapRest.py / restMomentum;
    const Double_t residueZ = -evapRest.pz / restMomentum;
    const Double_t storedTheta =
        static_cast<Double_t>(theta[step]) * TMath::DegToRad();
    const Double_t storedPhi =
        static_cast<Double_t>(phi[step]) * TMath::DegToRad();
    const Double_t storedX = std::sin(storedTheta) * std::cos(storedPhi);
    const Double_t storedY = std::sin(storedTheta) * std::sin(storedPhi);
    const Double_t storedZ = std::cos(storedTheta);
    // Chord between unit vectors: equals the angle when small, and phi drops
    // out on its own as sin(theta) goes to zero.
    const Double_t mismatch =
        std::sqrt((residueX - storedX) * (residueX - storedX) +
                  (residueY - storedY) * (residueY - storedY) +
                  (residueZ - storedZ) * (residueZ - storedZ));

    // Bound from rounded lab KE/angles and parent velocity carried through the
    // boost, plus the rounded CM angles; the tolerance is ten times that.
    const Double_t labMomentum = MomentumMagnitude(evapLab);
    const Double_t labError = epsilon * (evapKinetic * evapLab.e / labMomentum +
                                         2.0 * TMath::Pi() * labMomentum);
    const Double_t betaError = (1.0 + beta) * parentError / parent.e;
    const Double_t restError =
        gamma * (1.0 + beta) * (labError + evapLab.e * betaError);
    const Double_t tolerance =
        10.0 * (2.0 * TMath::Pi() * epsilon + restError / restMomentum);
    if (mismatch > tolerance)
      return Fail(label + " stored CM angle is not the one used: mismatch " +
                  std::to_string(mismatch) + " rad, tolerance " +
                  std::to_string(tolerance) + " rad");

    parent.e -= evapLab.e;
    parent.px -= evapLab.px;
    parent.py -= evapLab.py;
    parent.pz -= evapLab.pz;
    parentError += labError;
  }
  return true;
}

struct BoundaryRecord {
  Float_t exitEnergy = 0.0f;
  Float_t z = 0.0f;
  Int_t termination = 0;
};

std::optional<BoundaryRecord> Boundary(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return std::nullopt;
  BoundaryRecord value;
  root.truth->SetBranchAddress("beam_energy_exit", &value.exitEnergy);
  root.truth->SetBranchAddress("beam_stop_z", &value.z);
  root.truth->SetBranchAddress("beam_termination", &value.termination);
  root.truth->GetEntry(0);
  if (value.termination != 3 || std::fabs(value.z - 35.516f) > 2e-5f ||
      value.exitEnergy <= 0.0f) {
    Fail(std::string(path) +
         " did not terminate exactly on the downstream face");
    return std::nullopt;
  }
  return value;
}

using EventKey = std::pair<Int_t, ULong64_t>;
using EventRows = std::map<EventKey, std::vector<double>>;

std::optional<EventRows> Rows(const char *path) {
  RootFile root = Open(path);
  if (!root.events || !root.truth)
    return std::nullopt;
  std::vector<TLeaf *> leaves;
  for (TTree *tree : {root.events, root.truth}) {
    TIter next(tree->GetListOfLeaves());
    while (auto *leaf = dynamic_cast<TLeaf *>(next()))
      leaves.push_back(leaf);
  }
  EventRows rows;
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.events->GetEntry(entry);
    root.truth->GetEntry(entry);
    const auto requested = static_cast<Int_t>(
        root.truth->GetLeaf("requested_strip")->GetValueLong64());
    const auto eventIndex = static_cast<ULong64_t>(
        root.truth->GetLeaf("event_index")->GetValueLong64());
    std::vector<double> values;
    for (TLeaf *leaf : leaves)
      for (Int_t element = 0; element < leaf->GetLen(); ++element)
        values.push_back(leaf->GetValue(element));
    if (!rows.emplace(EventKey{requested, eventIndex}, std::move(values))
             .second) {
      Fail(std::string(path) + " has duplicate event identities");
      return std::nullopt;
    }
  }
  return rows;
}

// Closed-form lab kinetic energies of ejectile 3 at lab angle theta, for a beam
// of kinetic energy beamKinetic on a target at rest. From m4^2 = (E - E3)^2 -
// |P - p3|^2, which is a quadratic in p3; each root is -1 when not physical.
void TwoBodyEnergies(Double_t beamKinetic, Double_t theta, Double_t beamMass,
                     Double_t targetMass, Double_t ejectileMass,
                     Double_t residueMass, Double_t &lowRoot,
                     Double_t &highRoot) {
  lowRoot = -1.0;
  highRoot = -1.0;
  const Double_t total = beamKinetic + beamMass + targetMass;
  const Double_t momentum =
      std::sqrt(beamKinetic * beamKinetic + 2.0 * beamMass * beamKinetic);
  const Double_t cosTheta = std::cos(theta);
  const Double_t a = total * total - momentum * momentum +
                     ejectileMass * ejectileMass - residueMass * residueMass;
  const Double_t denominator =
      total * total - momentum * momentum * cosTheta * cosTheta;
  const Double_t discriminant =
      a * a - 4.0 * ejectileMass * ejectileMass * denominator;
  if (discriminant < 0.0 || denominator <= 0.0)
    return;
  const Double_t spread = total * std::sqrt(discriminant);
  const Double_t along = a * momentum * cosTheta;
  const Double_t roots[2] = {(along - spread) / (2.0 * denominator),
                             (along + spread) / (2.0 * denominator)};
  Double_t energies[2] = {-1.0, -1.0};
  for (Int_t i = 0; i < 2; ++i)
    // Squaring admits roots with E3 < 0 (E * E3 = a/2 + P p3 cos(theta)).
    if (roots[i] >= 0.0 && 0.5 * a + momentum * roots[i] * cosTheta > 0.0)
      energies[i] =
          std::sqrt(roots[i] * roots[i] + ejectileMass * ejectileMass) -
          ejectileMass;
  lowRoot = energies[0];
  highRoot = energies[1];
}

// Every event of tests/controls/two-body-ground.toml must put its proton at the
// energy two-body kinematics fixes for its lab angle and vertex beam energy.
// This shares no code with the simulator's CM sampling and boost.
bool TwoBodyGround(const char *path) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  Int_t steps = 0;
  Int_t reacted = -1;
  Float_t reactionEnergy = 0.0f;
  Float_t evapEnergy[10]{};
  Float_t evapTheta[10]{};
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("beam_energy_reaction", &reactionEnergy);
  root.truth->SetBranchAddress("evap_energy", evapEnergy);
  root.truth->SetBranchAddress("evap_theta", evapTheta);

  NuclideFinder finder;
  const Double_t beamMass = finder.GetMass("87Rb", "MeV/c^2");
  const Double_t targetMass = finder.GetMass("4He", "MeV/c^2");
  const Double_t protonMass = finder.GetMass("1H", "MeV/c^2");
  const Double_t residueMass = finder.GetMass("90Sr", "MeV/c^2");
  const Double_t epsilon =
      static_cast<Double_t>(std::numeric_limits<Float_t>::epsilon());

  Long64_t checked = 0;
  Double_t worstRatio = 0.0;
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    if (reacted < 0)
      continue;
    if (steps != 1)
      return Fail("two-body fixture has an unexpected step count");
    const Double_t beamKinetic = static_cast<Double_t>(reactionEnergy);
    const Double_t theta =
        static_cast<Double_t>(evapTheta[0]) * TMath::DegToRad();
    const Double_t stored = static_cast<Double_t>(evapEnergy[0]);
    Double_t lowRoot = -1.0;
    Double_t highRoot = -1.0;
    TwoBodyEnergies(beamKinetic, theta, beamMass, targetMass, protonMass,
                    residueMass, lowRoot, highRoot);
    // Rounding budget: the Float_t truth carries relative epsilon in energy
    // and absolute epsilon * pi in angle; propagate both through the formula.
    Double_t lowShifted = -1.0;
    Double_t highShifted = -1.0;
    TwoBodyEnergies(beamKinetic * (1.0 + epsilon),
                    theta + epsilon * TMath::Pi(), beamMass, targetMass,
                    protonMass, residueMass, lowShifted, highShifted);
    Double_t best = std::numeric_limits<Double_t>::infinity();
    Double_t bestTolerance = 0.0;
    const Double_t candidates[2] = {lowRoot, highRoot};
    const Double_t shifted[2] = {lowShifted, highShifted};
    for (Int_t i = 0; i < 2; ++i) {
      if (candidates[i] < 0.0)
        continue;
      const Double_t tolerance =
          10.0 * (std::fabs(shifted[i] - candidates[i]) + epsilon * stored);
      if (std::fabs(stored - candidates[i]) < best) {
        best = std::fabs(stored - candidates[i]);
        bestTolerance = tolerance;
      }
    }
    if (!std::isfinite(best) || best > bestTolerance)
      return Fail("two-body event " + std::to_string(entry) + " has proton " +
                  std::to_string(stored) + " MeV at " +
                  std::to_string(evapTheta[0]) + " deg; closed form gives " +
                  std::to_string(lowRoot) + " / " + std::to_string(highRoot));
    worstRatio = std::max(worstRatio, best / bestTolerance);
    ++checked;
  }
  if (checked < 150)
    return Fail("two-body fixture realized only " + std::to_string(checked) +
                " reactions");
  std::cout << "two-body: " << checked << " events, worst deviation "
            << worstRatio << " of tolerance" << std::endl;
  return true;
}

// Make a file look as if an older build with another default wrote it, so the
// next update run must refuse to append.
bool AlterResolved(const char *path) {
  std::unique_ptr<TFile> file(TFile::Open(path, "UPDATE"));
  if (!file || file->IsZombie())
    return Fail(std::string("cannot open ") + path);
  const TObjString *stored =
      dynamic_cast<TObjString *>(file->Get("metadata/resolved_configuration"));
  if (!stored)
    return Fail(std::string(path) + " has no resolved configuration");
  TString altered = stored->GetString();
  altered.ReplaceAll("z_effective = \"winger\"", "z_effective = \"atima14\"");
  if (altered == stored->GetString())
    return Fail(std::string(path) + " does not record z_effective = winger");
  file->cd("metadata");
  TObjString replacement(altered.Data());
  if (replacement.Write("resolved_configuration", TObject::kOverwrite) <= 0)
    return Fail(std::string("cannot rewrite metadata in ") + path);
  return true;
}

// The stored control text must be byte-identical to the control that ran.
bool ControlText(const char *rootPath, const char *controlPath) {
  std::unique_ptr<TFile> file(TFile::Open(rootPath, "READ"));
  if (!file || file->IsZombie())
    return Fail(std::string("cannot open ") + rootPath);
  const TObjString *stored =
      dynamic_cast<TObjString *>(file->Get("metadata/control_file_toml"));
  if (!stored)
    return Fail(std::string(rootPath) + " has no metadata/control_file_toml");
  std::ifstream input(controlPath, std::ios::binary);
  if (!input)
    return Fail(std::string("cannot read ") + controlPath);
  std::ostringstream expected;
  expected << input.rdbuf();
  if (expected.str().empty() ||
      std::string(stored->GetString().Data()) != expected.str())
    return Fail(std::string(rootPath) +
                " stores control text that differs "
                "from " +
                controlPath);
  return true;
}

bool Reproducible(const char *left, const char *right) {
  const auto a = Rows(left);
  const auto b = Rows(right);
  if (!a || !b || a->size() != b->size())
    return Fail(std::string(left) + " and " + right +
                " have different row counts");
  auto ai = a->begin();
  auto bi = b->begin();
  for (; ai != a->end(); ++ai, ++bi) {
    if (ai->first != bi->first || ai->second != bi->second)
      return Fail(std::string(left) + " and " + right +
                  " have different event data");
  }
  return true;
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    std::cerr << "usage: root_checks MODE FILE [ARGS...]\n";
    return 2;
  }
  const std::string mode = argv[1];
  bool ok = false;
  if (mode == "standard" && argc == 4)
    ok = Standard(argv[2], std::stoll(argv[3]));
  else if (mode == "conservation" && argc == 3)
    ok = Conservation(argv[2]);
  else if (mode == "forbidden" && argc == 3)
    ok = Forbidden(argv[2]);
  else if (mode == "multistep" && argc == 3)
    ok = MultiStep(argv[2]);
  else if (mode == "reconstruct" && argc == 3)
    ok = MultiStepReconstruction(argv[2]);
  else if (mode == "stopped" && argc == 4)
    ok = StoppedBeforeVertex(argv[2], std::stoi(argv[3]));
  else if (mode == "side-exit" && argc == 3)
    ok = SideExit(argv[2]);
  else if (mode == "boundary" && argc == 4) {
    const auto coarse = Boundary(argv[2]);
    const auto fine = Boundary(argv[3]);
    ok = coarse && fine &&
         std::fabs(coarse->exitEnergy - fine->exitEnergy) < 0.03f;
    if (!ok)
      Fail("coarse/fine transport does not converge at the boundary");
  } else if (mode == "compare" && argc == 4)
    ok = Reproducible(argv[2], argv[3]);
  else if (mode == "control-text" && argc == 4)
    ok = ControlText(argv[2], argv[3]);
  else if (mode == "alter-resolved" && argc == 3)
    ok = AlterResolved(argv[2]);
  else if (mode == "two-body" && argc == 3)
    ok = TwoBodyGround(argv[2]);
  else if (mode == "unique" && argc == 3) {
    const auto rows = Rows(argv[2]);
    ok = rows.has_value();
  } else {
    std::cerr << "root_checks: invalid mode or argument count\n";
    return 2;
  }
  return ok ? 0 : 1;
}
