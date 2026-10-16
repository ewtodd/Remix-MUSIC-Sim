#include "NuclideFinder.hpp"

#include <TFile.h>
#include <TFriendElement.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TList.h>
#include <TMD5.h>
#include <TMath.h>
#include <TObjString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
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

// True when `text` holds `line` as one whole line.
bool HasLine(const TString &text, const std::string &line) {
  std::istringstream lines(text.Data());
  std::string current;
  while (std::getline(lines, current))
    if (current == line)
      return true;
  return false;
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
  if (!resolved || !HasLine(resolved->GetString(), "schema_version = 6"))
    return Fail(std::string(path) + " does not declare output schema 6");
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
  Float_t excitation[10]{};
  Int_t branch[10]{};
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
  truth->SetBranchAddress("residue_excitation", excitation);
  truth->SetBranchAddress("reaction_branch", branch);
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
        evapTermination[step] != 0 || excitation[step] != -2.0f ||
        branch[step] != -1)
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

// branches*.toml: an (a,a') step with an elastic branch 0 (ground, Rutherford)
// and a Talys branch 1 (isotropic), drawn with `expectedFraction` for branch 1.
// The elastic branch is open at every vertex, so every event must react.
bool Branches(const char *path, double expectedFraction) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  Int_t reacted = -1;
  Float_t beamEnergy = 0.0f;
  Float_t evapEnergy[10]{};
  Float_t residueEnergy[10]{};
  Float_t excitation[10]{};
  Float_t theta[10]{};
  Int_t branch[10]{};
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("beam_energy_reaction", &beamEnergy);
  root.truth->SetBranchAddress("evap_energy", evapEnergy);
  root.truth->SetBranchAddress("residue_energy", residueEnergy);
  root.truth->SetBranchAddress("residue_excitation", excitation);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("reaction_branch", branch);

  Long64_t counts[2] = {0, 0};
  double thetaSum[2] = {0.0, 0.0};
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    if (reacted < 0)
      return Fail("an event was rejected though its elastic branch was open");
    if (branch[0] != 0 && branch[0] != 1)
      return Fail("a reacted event has no valid branch");
    if (branch[0] == 0 && excitation[0] != 0.0f)
      return Fail("the elastic branch left the residue excited");
    if (branch[0] == 1 && !(excitation[0] > 0.0f))
      return Fail("the Talys branch left the residue in its ground state");
    // Elastic or inelastic scattering has Q = 0: beam KE = KE out + E_ex.
    const double balance = static_cast<double>(beamEnergy) -
                           static_cast<double>(evapEnergy[0]) -
                           static_cast<double>(residueEnergy[0]) -
                           static_cast<double>(excitation[0]);
    if (std::fabs(balance) > 1e-3)
      return Fail("(a,a') step does not conserve energy: off by " +
                  std::to_string(balance) + " MeV");
    ++counts[branch[0]];
    thetaSum[branch[0]] += static_cast<double>(theta[0]);
  }

  const double total = static_cast<double>(counts[0] + counts[1]);
  if (counts[0] == 0 || total < 100.0)
    return Fail("too few events on the elastic branch");
  // At an expected fraction of 0, sigma is 0: a closed branch is never drawn.
  const double fraction = static_cast<double>(counts[1]) / total;
  const double sigma =
      std::sqrt(expectedFraction * (1.0 - expectedFraction) / total);
  if (std::fabs(fraction - expectedFraction) > 5.0 * sigma)
    return Fail("branch 1 drawn in " + std::to_string(fraction) +
                " of events, expected " + std::to_string(expectedFraction));
  // Rutherford from 1 degree is strongly forward; isotropic averages 90.
  const double elasticTheta = thetaSum[0] / static_cast<double>(counts[0]);
  const double talysTheta =
      counts[1] > 0 ? thetaSum[1] / static_cast<double>(counts[1]) : 90.0;
  if (elasticTheta > 45.0 || std::fabs(talysTheta - 90.0) > 10.0)
    return Fail("branch angular models were not applied per branch");
  return true;
}

// talys-step.toml: a step-level talys model on every reacted event, with
// E_ex inside the fixture's support (0, 6] MeV.
bool Talys(const char *path) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  Int_t reacted = -1;
  Float_t excitation[10]{};
  Int_t branch[10]{};
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("residue_excitation", excitation);
  root.truth->SetBranchAddress("reaction_branch", branch);
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    if (reacted < 0)
      return Fail("a step-level talys event did not react");
    if (branch[0] != 0 || !(excitation[0] > 0.0f) || excitation[0] > 6.0f)
      return Fail("step-level talys excitation is off the table support");
  }
  return root.truth->GetEntries() > 0;
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

// Prints each event's per-step CM angles, ordered by event identity. They are
// pure functions of the RNG draws, so they pin the reaction's draw order.
bool CmAngles(const char *path) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  Int_t steps = 0;
  Float_t theta[10]{};
  Float_t phi[10]{};
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("phi_cm", phi);
  std::map<EventKey, std::vector<std::string>> events;
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    const EventKey key{
        static_cast<Int_t>(
            root.truth->GetLeaf("requested_strip")->GetValueLong64()),
        static_cast<ULong64_t>(
            root.truth->GetLeaf("event_index")->GetValueLong64())};
    std::vector<std::string> &lines = events[key];
    for (Int_t step = 0; step < steps; ++step) {
      std::ostringstream line;
      line << std::setprecision(9) << key.first << ' ' << key.second << ' '
           << step << ' ' << static_cast<double>(theta[step]) << ' '
           << static_cast<double>(phi[step]);
      lines.push_back(line.str());
    }
  }
  for (std::map<EventKey, std::vector<std::string>>::const_iterator event =
           events.begin();
       event != events.end(); ++event)
    for (size_t i = 0; i < event->second.size(); ++i)
      std::cout << event->second[i] << std::endl;
  return !events.empty();
}

// Level table 0 in the metadata must hold the bytes of `tablePath` under the
// resolved `expectedPath`, and the resolved configuration must carry its MD5.
bool LevelTable(const char *rootPath, const char *expectedPath,
                const char *tablePath) {
  std::unique_ptr<TFile> file(TFile::Open(rootPath, "READ"));
  if (!file || file->IsZombie())
    return Fail(std::string("cannot open ") + rootPath);
  const TObjString *text =
      dynamic_cast<TObjString *>(file->Get("metadata/level_table_0"));
  const TObjString *path =
      dynamic_cast<TObjString *>(file->Get("metadata/level_table_0_path"));
  const TObjString *resolved =
      dynamic_cast<TObjString *>(file->Get("metadata/resolved_configuration"));
  if (!text || !path || !resolved || file->Get("metadata/level_table_1"))
    return Fail(std::string(rootPath) + " does not record exactly one table");
  std::ifstream input(tablePath, std::ios::binary);
  if (!input)
    return Fail(std::string("cannot read ") + tablePath);
  std::ostringstream expected;
  expected << input.rdbuf();
  if (std::string(text->GetString().Data()) != expected.str())
    return Fail(std::string(rootPath) +
                " stores table text that differs from " + tablePath);
  if (std::string(path->GetString().Data()) != expectedPath)
    return Fail(std::string(rootPath) + " records table path " +
                path->GetString().Data() + ", expected " + expectedPath);
  TMD5 md5;
  md5.Update(reinterpret_cast<const UChar_t *>(expected.str().data()),
             static_cast<UInt_t>(expected.str().size()));
  md5.Final();
  const std::string line =
      std::string(expectedPath) + "; md5 = " + md5.AsString();
  if (!resolved->GetString().Contains(line.c_str()))
    return Fail(std::string(rootPath) + " resolved configuration lacks '" +
                line + "'");
  return true;
}

// The number after `key` in `line`, e.g. "ecm_mev = " -> 8.11.
double NumberAfter(const std::string &line, const std::string &key) {
  const size_t at = line.find(key);
  return at == std::string::npos ? std::nan("")
                                 : std::atof(line.c_str() + at + key.size());
}

// Sum of the level and bin rows in step 0 of channel `name` of a talys-cache
// table, which must carry the channel's cross section; 0 if absent.
double ChannelStepSum(const std::string &text, const std::string &name) {
  std::istringstream lines(text);
  std::string line, channel;
  int step = -1;
  double sum = 0.0;
  while (std::getline(lines, line)) {
    std::istringstream fields(line);
    std::string tag;
    double a = 0.0, b = 0.0, c = 0.0;
    if (!(fields >> tag) || tag[0] == '#')
      continue;
    if (tag == "channel" && (fields >> channel))
      step = -1;
    else if (tag == "step")
      fields >> step;
    else if (channel != name || step != 0)
      continue;
    else if (tag == "level" && (fields >> a >> b))
      sum += b;
    else if (tag == "bin" && (fields >> a >> b >> c))
      sum += c;
  }
  return sum;
}

// A table's elastic rows integrated against Rutherford above thetaMinDeg, in
// units of (pi/4) k^2: cot^2(thetaMin/2) when the table has no elastic rows.
double ElasticIntegral(const std::string &text, double thetaMinDeg) {
  std::vector<double> theta, ratio;
  std::istringstream lines(text);
  std::string line;
  while (std::getline(lines, line)) {
    std::istringstream fields(line);
    std::string tag;
    double angle = 0.0, value = 0.0;
    if ((fields >> tag) && tag == "elastic" && (fields >> angle >> value)) {
      theta.push_back(angle);
      ratio.push_back(value);
    }
  }
  if (theta.empty()) {
    theta = {0.0, 180.0};
    ratio = {1.0, 1.0};
  }
  const double deg = TMath::Pi() / 180.0;
  double sum = 0.0;
  for (size_t i = 0; i + 1 < theta.size(); ++i) {
    double a = theta[i];
    double ra = ratio[i];
    if (theta[i + 1] <= thetaMinDeg)
      continue;
    if (a < thetaMinDeg) {
      ra += (ratio[i + 1] - ra) * (thetaMinDeg - a) / (theta[i + 1] - a);
      a = thetaMinDeg;
    }
    const double ua = std::pow(std::sin(0.5 * a * deg), 2);
    const double ub = std::pow(std::sin(0.5 * theta[i + 1] * deg), 2);
    sum += (1.0 / ua - 1.0 / ub) * 0.5 * (ra + ratio[i + 1]);
  }
  return sum;
}

// An aa class run, per strip: the elastic and (a,a') weights are the stored
// table's elastic integral and level sum, and channels, elastic angles past
// 150 deg, and vertex energies follow them.
bool ComputedWeights(const char *path, int chargeProduct, double thetaMinDeg) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  const TObjString *resolved = dynamic_cast<TObjString *>(
      root.file->Get("metadata/resolved_configuration"));
  if (!resolved)
    return Fail(std::string(path) + " has no resolved configuration");
  std::map<std::string, std::string> tables;
  for (int i = 0;; ++i) {
    const std::string name = "metadata/level_table_" + std::to_string(i);
    const TObjString *text =
        dynamic_cast<TObjString *>(root.file->Get(name.c_str()));
    const TObjString *tablePath =
        dynamic_cast<TObjString *>(root.file->Get((name + "_path").c_str()));
    if (!text || !tablePath)
      break;
    tables[tablePath->GetString().Data()] = text->GetString().Data();
  }

  const std::string prefix = "talys_strip_";
  std::map<int, double> inelasticShare;
  std::map<int, double> centreEnergy;
  std::map<int, double> elasticTail; // P(theta_cm > 150 deg | elastic)
  std::istringstream lines(resolved->GetString().Data());
  std::string line;
  while (std::getline(lines, line)) {
    if (line.rfind(prefix, 0) != 0)
      continue;
    const int strip = std::atoi(line.c_str() + prefix.size());
    const double ecm = NumberAfter(line, "ecm_mev = ");
    const size_t weightsAt = line.find("channels_mb = ");
    double elastic = std::nan("");
    double inelastic = std::nan("");
    if (weightsAt == std::string::npos ||
        std::sscanf(line.c_str() + weightsAt + 14, "%lf, %lf", &elastic,
                    &inelastic) != 2)
      return Fail("cannot read the weights in '" + line + "'");
    const size_t tableAt = line.find("table = ") + 8;
    const size_t md5At = line.find("; md5 = ");
    const std::string table = line.substr(tableAt, md5At - tableAt);
    const std::string md5Text = line.substr(md5At + 8, 32);
    const std::map<std::string, std::string>::const_iterator stored =
        tables.find(table);
    if (stored == tables.end())
      return Fail("strip " + std::to_string(strip) + " table " + table +
                  " is not stored in the metadata");
    const double k = chargeProduct * 1.439964548 / ecm;
    const double expected = 10.0 * 0.25 * TMath::Pi() * k * k *
                            ElasticIntegral(stored->second, thetaMinDeg);
    if (!(std::fabs(elastic - expected) <= 1e-9 * expected))
      return Fail("strip " + std::to_string(strip) + " elastic weight " +
                  std::to_string(elastic) +
                  " mb is not its table's elastic integral " +
                  std::to_string(expected));
    elasticTail[strip] = ElasticIntegral(stored->second, 150.0) /
                         ElasticIntegral(stored->second, thetaMinDeg);
    TMD5 md5;
    md5.Update(reinterpret_cast<const UChar_t *>(stored->second.data()),
               static_cast<UInt_t>(stored->second.size()));
    md5.Final();
    if (md5Text != md5.AsString())
      return Fail("strip " + std::to_string(strip) + " records the wrong MD5");
    const double sum = ChannelStepSum(stored->second, "a");
    if (!(std::fabs(inelastic - sum) <= 1e-7 * sum))
      return Fail("strip " + std::to_string(strip) + " inelastic weight " +
                  std::to_string(inelastic) +
                  " mb is not its (a,a') channel's sum " + std::to_string(sum));
    inelasticShare[strip] = inelastic / (elastic + inelastic);
    centreEnergy[strip] = NumberAfter(line, "beam_energy_mev = ");
  }
  if (inelasticShare.empty())
    return Fail(std::string(path) + " records no per-strip weights");

  Int_t requested = -1;
  Int_t reacted = -1;
  Int_t channel = -1;
  root.truth->SetBranchAddress("requested_strip", &requested);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("reaction_channel", &channel);
  Float_t beamEnergy = 0.0f;
  root.truth->SetBranchAddress("beam_energy_reaction", &beamEnergy);
  Float_t theta[10]{};
  root.truth->SetBranchAddress("theta_cm", theta);
  std::map<int, std::pair<long, long>> elasticCounts; // {events, past 150}
  std::map<int, double> energySum;
  std::map<int, std::pair<long, long>> counts; // strip -> {events, (a,a')}
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    if (reacted < 0)
      return Fail("an event was rejected though its elastic channel was open");
    ++counts[requested].first;
    energySum[requested] += static_cast<double>(beamEnergy);
    if (channel == 0) {
      ++elasticCounts[requested].first;
      if (theta[0] > 150.0f)
        ++elasticCounts[requested].second;
    }
    if (channel == 1)
      ++counts[requested].second;
  }
  for (std::map<int, double>::const_iterator share = inelasticShare.begin();
       share != inelasticShare.end(); ++share) {
    const std::pair<long, long> &count = counts[share->first];
    const double n = static_cast<double>(count.first);
    const double p = share->second;
    const double fraction = static_cast<double>(count.second) / n;
    if (count.first < 100 ||
        std::fabs(fraction - p) > 5.0 * std::sqrt(p * (1.0 - p) / n))
      return Fail("strip " + std::to_string(share->first) + " drew (a,a') in " +
                  std::to_string(fraction) + " of events, expected " +
                  std::to_string(p));
    // The vertex is uniform in the strip; one strip spans ~4 MeV of beam.
    const double meanEnergy = energySum[share->first] / n;
    if (std::fabs(meanEnergy - centreEnergy[share->first]) > 0.15)
      return Fail("strip " + std::to_string(share->first) +
                  " table is for beam energy " +
                  std::to_string(centreEnergy[share->first]) +
                  " MeV, but its vertices average " +
                  std::to_string(meanEnergy));
  }
  // Elastic angles follow each strip's table, pooled over the strips.
  double observed = 0.0, expected = 0.0, variance = 0.0;
  for (std::map<int, double>::const_iterator tail = elasticTail.begin();
       tail != elasticTail.end(); ++tail) {
    const std::pair<long, long> &count = elasticCounts[tail->first];
    observed += static_cast<double>(count.second);
    expected += static_cast<double>(count.first) * tail->second;
    variance +=
        static_cast<double>(count.first) * tail->second * (1.0 - tail->second);
  }
  if (std::fabs(observed - expected) > 5.0 * std::sqrt(variance))
    return Fail(std::to_string(observed) +
                " elastic events past 150 deg, expected " +
                std::to_string(expected));
  return true;
}

// A class channel's steps, read from its description, e.g.
// "37Cl(a,2n)39K: 1n + 40K, then 1n + 39K" -> {{"1n", "40K"}, {"1n", "39K"}}.
std::vector<std::pair<std::string, std::string>>
ChainSteps(const std::string &description) {
  std::vector<std::pair<std::string, std::string>> steps;
  std::string rest = description.substr(description.find(": ") + 2);
  while (!rest.empty()) {
    const size_t then = rest.find(", then ");
    const std::string step = rest.substr(0, then);
    const size_t plus = step.find(" + ");
    steps.push_back(
        std::make_pair(step.substr(0, plus), step.substr(plus + 3)));
    rest = then == std::string::npos ? "" : rest.substr(then + 7);
  }
  return steps;
}

double ProductMass(NuclideFinder &finder, const std::string &name) {
  return name == "gamma" ? 0.0 : finder.GetMass(name, "MeV/c^2");
}

// A talys_class run: per strip, channels are drawn in proportion to their
// recorded cross sections; each event runs its channel's whole chain, every
// step but the last leaving enough excitation for the next emission, and the
// slots past the chain hold the sentinels. With unreactedOk, a beam may stop
// short, but every strip must still have reactions.
bool ClassChannels(const char *path, bool unreactedOk) {
  RootFile root = Open(path);
  if (!root.truth)
    return false;
  const TObjString *resolved = dynamic_cast<TObjString *>(
      root.file->Get("metadata/resolved_configuration"));
  const TObjString *legend =
      dynamic_cast<TObjString *>(root.file->Get("metadata/channels"));
  if (!resolved || !legend)
    return Fail(std::string(path) + " lacks its channel metadata");
  NuclideFinder finder;
  std::vector<std::vector<std::pair<std::string, std::string>>> chains;
  std::vector<std::vector<double>> thresholds; // per step, but the last
  size_t longest = 0;
  std::map<int, std::vector<double>> sigma; // strip -> per channel [mb]
  std::istringstream lines(resolved->GetString().Data());
  std::string line;
  while (std::getline(lines, line)) {
    if (line.rfind("reaction_channel_", 0) == 0) {
      const size_t open = line.find('"');
      const std::string description =
          line.substr(open + 1, line.rfind('"') - open - 1);
      if (!legend->GetString().Contains(
              (std::to_string(chains.size()) + " = " + description).c_str()))
        return Fail("metadata/channels lacks " + description);
      chains.push_back(ChainSteps(description));
      const std::vector<std::pair<std::string, std::string>> &chain =
          chains.back();
      longest = std::max(longest, chain.size());
      thresholds.push_back(std::vector<double>());
      for (size_t k = 0; k + 1 < chain.size(); ++k)
        thresholds.back().push_back(ProductMass(finder, chain[k + 1].first) +
                                    ProductMass(finder, chain[k + 1].second) -
                                    ProductMass(finder, chain[k].second));
    }
    const size_t at = line.find("channels_mb = ");
    if (line.rfind("talys_strip_", 0) != 0 || at == std::string::npos)
      continue;
    const int strip = std::atoi(line.c_str() + 12);
    std::istringstream values(line.substr(at + 14));
    double value = 0.0;
    char comma = ',';
    while (comma == ',' && (values >> value)) {
      sigma[strip].push_back(value);
      values >> comma;
    }
  }
  if (chains.empty() || sigma.empty())
    return Fail(std::string(path) + " records no class channels");

  Int_t steps = 0, requested = -1, reacted = -1, channel = -1, residueStep = 0;
  Float_t evapEnergy[20]{}, residueEnergy[20]{}, excitation[20]{}, theta[20]{};
  Int_t branch[20]{}, evapTermination[20]{};
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("requested_strip", &requested);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("reaction_channel", &channel);
  root.truth->SetBranchAddress("residue_step", &residueStep);
  root.truth->SetBranchAddress("evap_energy", evapEnergy);
  root.truth->SetBranchAddress("residue_energy", residueEnergy);
  root.truth->SetBranchAddress("residue_excitation", excitation);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("reaction_branch", branch);
  root.truth->SetBranchAddress("evap_termination", evapTermination);
  std::map<int, std::map<int, long>> counts;
  std::map<int, long> events;
  for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
    root.truth->GetEntry(entry);
    const std::string where = "entry " + std::to_string(entry);
    if (static_cast<size_t>(steps) != longest)
      return Fail(where + " has n_steps " + std::to_string(steps) +
                  ", not the longest chain's " + std::to_string(longest));
    // A beam that stops short still puts its strip in `events`.
    if (unreactedOk && reacted < 0 && channel == -1) {
      events.insert(std::make_pair(static_cast<int>(requested), 0L));
      continue;
    }
    if (reacted < 0 || channel < 0 ||
        channel >= static_cast<Int_t>(chains.size()))
      return Fail(where + " did not react through a listed channel");
    ++events[requested];
    ++counts[requested][channel];
    const size_t length = chains[channel].size();
    if (residueStep != static_cast<Int_t>(length) - 1)
      return Fail(where + " transports the residue of step " +
                  std::to_string(residueStep));
    for (size_t k = 0; k < longest; ++k) {
      const bool active = k < length;
      if (active && (evapEnergy[k] < 0.0f || residueEnergy[k] < 0.0f ||
                     excitation[k] < 0.0f || branch[k] != 0 ||
                     theta[k] == kAngleUnset || evapTermination[k] == 0))
        return Fail(where + " step " + std::to_string(k) + " did not happen");
      if (!active && (evapEnergy[k] != -2.0f || residueEnergy[k] != -2.0f ||
                      excitation[k] != -2.0f || branch[k] != -1 ||
                      theta[k] != kAngleUnset || evapTermination[k] != 0))
        return Fail(where + " step " + std::to_string(k) +
                    " past its chain is not at the sentinels");
      if (active && k + 1 < length &&
          static_cast<double>(excitation[k]) < thresholds[channel][k] - 1e-5)
        return Fail(where + " step " + std::to_string(k) + " left " +
                    std::to_string(excitation[k]) +
                    " MeV, below the next emission's " +
                    std::to_string(thresholds[channel][k]));
    }
  }
  for (std::map<int, long>::const_iterator strip = events.begin();
       strip != events.end(); ++strip)
    if (strip->second == 0)
      return Fail("strip " + std::to_string(strip->first) +
                  " has no reactions");
  for (std::map<int, std::vector<double>>::const_iterator strip = sigma.begin();
       strip != sigma.end(); ++strip) {
    double total = 0.0;
    for (size_t c = 0; c < strip->second.size(); ++c)
      total += strip->second[c];
    const double n = static_cast<double>(events[strip->first]);
    for (size_t c = 0; c < strip->second.size(); ++c) {
      const double p = strip->second[c] / total;
      const double fraction =
          static_cast<double>(counts[strip->first][static_cast<int>(c)]) / n;
      if (std::fabs(fraction - p) > 5.0 * std::sqrt(p * (1.0 - p) / n) + 1e-12)
        return Fail("strip " + std::to_string(strip->first) + " drew channel " +
                    std::to_string(c) + " in " + std::to_string(fraction) +
                    " of events, expected " + std::to_string(p));
    }
  }
  return true;
}

// metadata/<name> must exist and contain `text`.
bool MetadataContains(const char *path, const char *name, const char *text) {
  std::unique_ptr<TFile> file(TFile::Open(path, "READ"));
  const TObjString *value =
      file && !file->IsZombie()
          ? dynamic_cast<TObjString *>(
                file->Get((std::string("metadata/") + name).c_str()))
          : nullptr;
  if (!value || !value->GetString().Contains(text))
    return Fail(std::string(path) + " metadata/" + name + " lacks '" + text +
                "'");
  return true;
}

// Runs of one control with different seeds must not replay each other's
// events: at most 1% of `right`'s may share a strip, vertex, and CM angle.
bool Independent(const char *left, const char *right) {
  std::set<std::tuple<int, float, float>> seen;
  int shared = 0, total = 0;
  for (int side = 0; side < 2; ++side) {
    RootFile root = Open(side == 0 ? left : right);
    if (!root.truth)
      return false;
    Int_t strip = -1;
    Float_t vertex = 0.0f, theta[20]{};
    root.truth->SetBranchAddress("requested_strip", &strip);
    root.truth->SetBranchAddress("vertex_z", &vertex);
    root.truth->SetBranchAddress("theta_cm", theta);
    for (Long64_t entry = 0; entry < root.truth->GetEntries(); ++entry) {
      root.truth->GetEntry(entry);
      const std::tuple<int, float, float> event(strip, vertex, theta[0]);
      if (side == 0) {
        seen.insert(event);
      } else {
        ++total;
        shared += seen.count(event) ? 1 : 0;
      }
    }
  }
  if (total == 0 || shared > total / 100)
    return Fail(std::to_string(shared) + " of " + std::to_string(total) +
                " events of " + right + " replay events of " + left);
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
  else if (mode == "branches" && argc == 4)
    ok = Branches(argv[2], std::stod(argv[3]));
  else if (mode == "talys" && argc == 3)
    ok = Talys(argv[2]);
  else if (mode == "boundary" && argc == 4) {
    const auto coarse = Boundary(argv[2]);
    const auto fine = Boundary(argv[3]);
    ok = coarse && fine &&
         std::fabs(coarse->exitEnergy - fine->exitEnergy) < 0.03f;
    if (!ok)
      Fail("coarse/fine transport does not converge at the boundary");
  } else if (mode == "compare" && argc == 4)
    ok = Reproducible(argv[2], argv[3]);
  else if (mode == "independent" && argc == 4)
    ok = Independent(argv[2], argv[3]);
  else if (mode == "control-text" && argc == 4)
    ok = ControlText(argv[2], argv[3]);
  else if (mode == "cm-angles" && argc == 3)
    ok = CmAngles(argv[2]);
  else if (mode == "level-table" && argc == 5)
    ok = LevelTable(argv[2], argv[3], argv[4]);
  else if (mode == "computed-weights" && argc == 5)
    ok = ComputedWeights(argv[2], std::stoi(argv[3]), std::stod(argv[4]));
  else if (mode == "class-channels" && argc == 3)
    ok = ClassChannels(argv[2], false);
  else if (mode == "class-channels" && argc == 4 &&
           std::string(argv[3]) == "unreacted-ok")
    ok = ClassChannels(argv[2], true);
  else if (mode == "metadata-contains" && argc == 5)
    ok = MetadataContains(argv[2], argv[3], argv[4]);
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
