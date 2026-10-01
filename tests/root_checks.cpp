#include <TFile.h>
#include <TFriendElement.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TList.h>
#include <TObjString.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace {

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
  if (!resolved || !resolved->GetString().Contains("schema_version = 3\n"))
    return Fail(std::string(path) + " does not declare output schema 3");
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

bool Forbidden(const char *path) {
  RootFile root = Open(path);
  if (!root.truth || root.truth->GetEntries() != 1)
    return false;
  Int_t steps = 0;
  Int_t requested = -2;
  Int_t reacted = -2;
  Int_t termination = 0;
  Float_t reactionEnergy = 0.0f;
  Float_t vertex = 0.0f;
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
  root.truth->SetBranchAddress("n_steps", &steps);
  root.truth->SetBranchAddress("requested_strip", &requested);
  root.truth->SetBranchAddress("reaction_strip", &reacted);
  root.truth->SetBranchAddress("beam_termination", &termination);
  root.truth->SetBranchAddress("beam_energy_reaction", &reactionEnergy);
  root.truth->SetBranchAddress("vertex_z", &vertex);
  root.truth->SetBranchAddress("evap_energy", evapEnergy);
  root.truth->SetBranchAddress("residue_energy", residueEnergy);
  root.truth->SetBranchAddress("evap_energy_exit", evapExitEnergy);
  root.truth->SetBranchAddress("residue_energy_exit", residueExitEnergy);
  root.truth->SetBranchAddress("theta_cm", theta);
  root.truth->SetBranchAddress("phi_cm", phi);
  root.truth->SetBranchAddress("evap_theta", evapTheta);
  root.truth->SetBranchAddress("evap_phi", evapPhi);
  root.truth->SetBranchAddress("residue_theta", residueTheta);
  root.truth->SetBranchAddress("residue_phi", residuePhi);
  root.truth->SetBranchAddress("evap_termination", evapTermination);
  root.truth->SetBranchAddress("residue_step", &residueStep);
  root.truth->SetBranchAddress("residue_termination", &residueTermination);
  root.truth->GetEntry(0);
  if (requested != 0 || reacted != -1 || termination != 3 ||
      reactionEnergy != -2.0f || vertex != -1000.0f)
    return Fail("forbidden reaction left an inconsistent event history");
  for (Int_t step = 0; step < steps; ++step)
    if (evapEnergy[step] != -2.0f || residueEnergy[step] != -2.0f ||
        evapExitEnergy[step] != -2.0f || residueExitEnergy[step] != -2.0f ||
        theta[step] != -1.0f || phi[step] != -1.0f ||
        evapTheta[step] != -1.0f || evapPhi[step] != -1.0f ||
        residueTheta[step] != -1.0f || residuePhi[step] != -1.0f ||
        evapTermination[step] != 0)
      return Fail("forbidden reaction retained partial-chain truth");
  if (residueStep != -1 || residueTermination != 0)
    return Fail("forbidden reaction retained a residue transport state");
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

bool Reproducible(const char *left, const char *right) {
  const auto a = Rows(left);
  const auto b = Rows(right);
  if (!a || !b || a->size() != b->size())
    return Fail("reproducibility files have different row counts");
  auto ai = a->begin();
  auto bi = b->begin();
  for (; ai != a->end(); ++ai, ++bi) {
    if (ai->first != bi->first || ai->second != bi->second)
      return Fail("one-thread and multi-thread event data differ");
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
  else if (mode == "boundary" && argc == 4) {
    const auto coarse = Boundary(argv[2]);
    const auto fine = Boundary(argv[3]);
    ok = coarse && fine &&
         std::fabs(coarse->exitEnergy - fine->exitEnergy) < 0.03f;
    if (!ok)
      Fail("coarse/fine transport does not converge at the boundary");
  } else if (mode == "compare" && argc == 4)
    ok = Reproducible(argv[2], argv[3]);
  else if (mode == "unique" && argc == 3) {
    const auto rows = Rows(argv[2]);
    ok = rows.has_value();
  } else {
    std::cerr << "root_checks: invalid mode or argument count\n";
    return 2;
  }
  return ok ? 0 : 1;
}
