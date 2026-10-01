#include "Simulator.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TObjString.h>

#ifndef MUSICSIM_VERSION
#define MUSICSIM_VERSION "dev"
#endif

// Initialize the "events_MeV" tree (detector-level output; branch layout
// mirrors the experimental "events" tree from the MUSIC EventBuilder) and the
// friended "MC" tree (truth-only).
// Energies are Float_t in MeV (the experimental data file uses ADC counts);
// analysis is expected to apply per-channel calibration to compare data
// against sim in MeV.
TTree *Simulator::InitTree(TFile *ROOTfile, TString FileOpt) {
  if (!ROOTfile || ROOTfile->IsZombie() || !ROOTfile->IsOpen()) {
    std::cerr << "musicsim ERROR: ROOT output file is not usable." << std::endl;
    return nullptr;
  }

  TTree *tree = nullptr;
  const Bool_t update = (FileOpt == "update" || FileOpt == "UPDATE");
  if (update) {
    tree = dynamic_cast<TTree *>(ROOTfile->Get("events_MeV"));
    MCTree = dynamic_cast<TTree *>(ROOTfile->Get("MC"));
    if (!tree || !MCTree) {
      std::cerr << "musicsim ERROR: update output must contain both "
                   "'events_MeV' and 'MC' trees."
                << std::endl;
      return nullptr;
    }
    if (tree->GetEntries() != MCTree->GetEntries()) {
      std::cerr << "musicsim ERROR: update trees have different entry counts."
                << std::endl;
      return nullptr;
    }
    if (MCTree->GetMinimum("n_steps") < 0.0 ||
        MCTree->GetMaximum("n_steps") > maxEvaporations) {
      std::cerr << "musicsim ERROR: update output has an unsafe n_steps value."
                << std::endl;
      return nullptr;
    }

    using Field = std::pair<const char *, const char *>;
    const std::vector<Field> eventFields = {{"LeftdE", "Float_t"},
                                            {"RightdE", "Float_t"},
                                            {"Strip0dE", "Float_t"},
                                            {"Strip17dE", "Float_t"},
                                            {"Cathode", "Float_t"}};
    const std::vector<Field> mcFields = {{"n_steps", "Int_t"},
                                         {"requested_strip", "Int_t"},
                                         {"reaction_strip", "Int_t"},
                                         {"event_index", "ULong64_t"},
                                         {"beam_energy_accel", "Float_t"},
                                         {"beam_energy_gas", "Float_t"},
                                         {"beam_energy_reaction", "Float_t"},
                                         {"beam_energy_exit", "Float_t"},
                                         {"beam_stop_x", "Float_t"},
                                         {"beam_stop_y", "Float_t"},
                                         {"beam_stop_z", "Float_t"},
                                         {"beam_stop_strip", "Int_t"},
                                         {"beam_termination", "Int_t"},
                                         {"vertex_x", "Float_t"},
                                         {"vertex_y", "Float_t"},
                                         {"vertex_z", "Float_t"},
                                         {"DeadUS_dE", "Float_t"},
                                         {"DeadDS_dE", "Float_t"},
                                         {"evap_energy", "Float_t"},
                                         {"residue_energy", "Float_t"},
                                         {"evap_energy_exit", "Float_t"},
                                         {"residue_energy_exit", "Float_t"},
                                         {"theta_cm", "Float_t"},
                                         {"phi_cm", "Float_t"},
                                         {"evap_theta", "Float_t"},
                                         {"evap_phi", "Float_t"},
                                         {"residue_theta", "Float_t"},
                                         {"residue_phi", "Float_t"},
                                         {"evap_stop_x", "Float_t"},
                                         {"evap_stop_y", "Float_t"},
                                         {"evap_stop_z", "Float_t"},
                                         {"evap_stop_strip", "Int_t"},
                                         {"evap_termination", "Int_t"},
                                         {"residue_stop_x", "Float_t"},
                                         {"residue_stop_y", "Float_t"},
                                         {"residue_stop_z", "Float_t"},
                                         {"residue_stop_strip", "Int_t"},
                                         {"residue_step", "Int_t"},
                                         {"residue_termination", "Int_t"}};
    auto validateSchema = [](TTree *candidate,
                             const std::vector<Field> &expected) {
      if (!candidate || candidate->GetListOfBranches()->GetEntries() !=
                            static_cast<Int_t>(expected.size()))
        return kFALSE;
      for (const auto &[name, type] : expected) {
        TBranch *branch = candidate->GetBranch(name);
        TLeaf *leaf = branch ? branch->GetLeaf(name) : nullptr;
        if (!leaf || std::string(leaf->GetTypeName()) != type)
          return kFALSE;
      }
      return kTRUE;
    };
    if (!validateSchema(tree, eventFields) ||
        !validateSchema(MCTree, mcFields)) {
      std::cerr << "musicsim ERROR: update output has an incompatible tree "
                   "schema."
                << std::endl;
      return nullptr;
    }
    if (tree->GetLeaf("LeftdE")->GetLenStatic() != N_SEG_STRIPS ||
        tree->GetLeaf("RightdE")->GetLenStatic() != N_SEG_STRIPS) {
      std::cerr << "musicsim ERROR: update output has incompatible detector "
                   "array lengths."
                << std::endl;
      return nullptr;
    }
    for (const char *name :
         {"evap_energy", "residue_energy", "evap_energy_exit",
          "residue_energy_exit", "theta_cm", "phi_cm", "evap_theta", "evap_phi",
          "residue_theta", "residue_phi", "evap_stop_x", "evap_stop_y",
          "evap_stop_z", "evap_stop_strip", "evap_termination"}) {
      TLeaf *leaf = MCTree->GetLeaf(name);
      if (!leaf->GetLeafCount() ||
          std::string(leaf->GetLeafCount()->GetName()) != "n_steps") {
        std::cerr << "musicsim ERROR: update branch '" << name
                  << "' has an incompatible array length." << std::endl;
        return nullptr;
      }
    }

    auto bind = [](TTree *candidate, const char *name, auto *address) {
      if (candidate->SetBranchAddress(name, address) < 0) {
        std::cerr << "musicsim ERROR: could not bind update branch '" << name
                  << "'." << std::endl;
        return kFALSE;
      }
      return kTRUE;
    };
    Bool_t bound = kTRUE;
    bound &= bind(tree, "LeftdE", LeftdE);
    bound &= bind(tree, "RightdE", RightdE);
    bound &= bind(tree, "Strip0dE", &Strip0dE);
    bound &= bind(tree, "Strip17dE", &Strip17dE);
    bound &= bind(tree, "Cathode", &Cathode);
    bound &= bind(MCTree, "n_steps", &n_steps);
    bound &= bind(MCTree, "requested_strip", &requested_strip);
    bound &= bind(MCTree, "reaction_strip", &reaction_strip);
    bound &= bind(MCTree, "event_index", &event_index);
    bound &= bind(MCTree, "beam_energy_accel", &beam_energy_accel);
    bound &= bind(MCTree, "beam_energy_gas", &beam_energy_gas);
    bound &= bind(MCTree, "beam_energy_reaction", &beam_energy_reaction);
    bound &= bind(MCTree, "beam_energy_exit", &beam_energy_exit);
    bound &= bind(MCTree, "beam_stop_x", &beam_stop_x);
    bound &= bind(MCTree, "beam_stop_y", &beam_stop_y);
    bound &= bind(MCTree, "beam_stop_z", &beam_stop_z);
    bound &= bind(MCTree, "beam_stop_strip", &beam_stop_strip);
    bound &= bind(MCTree, "beam_termination", &beam_termination);
    bound &= bind(MCTree, "vertex_x", &vertex_x);
    bound &= bind(MCTree, "vertex_y", &vertex_y);
    bound &= bind(MCTree, "vertex_z", &vertex_z);
    bound &= bind(MCTree, "DeadUS_dE", &DeadUS_dE);
    bound &= bind(MCTree, "DeadDS_dE", &DeadDS_dE);
    bound &= bind(MCTree, "evap_energy", evap_energy);
    bound &= bind(MCTree, "residue_energy", residue_energy);
    bound &= bind(MCTree, "evap_energy_exit", evap_energy_exit);
    bound &= bind(MCTree, "residue_energy_exit", residue_energy_exit);
    bound &= bind(MCTree, "theta_cm", theta_cm);
    bound &= bind(MCTree, "phi_cm", phi_cm);
    bound &= bind(MCTree, "evap_theta", evap_theta);
    bound &= bind(MCTree, "evap_phi", evap_phi);
    bound &= bind(MCTree, "residue_theta", residue_theta);
    bound &= bind(MCTree, "residue_phi", residue_phi);
    bound &= bind(MCTree, "evap_stop_x", evap_stop_x);
    bound &= bind(MCTree, "evap_stop_y", evap_stop_y);
    bound &= bind(MCTree, "evap_stop_z", evap_stop_z);
    bound &= bind(MCTree, "evap_stop_strip", evap_stop_strip);
    bound &= bind(MCTree, "evap_termination", evap_termination);
    bound &= bind(MCTree, "residue_stop_x", &residue_stop_x);
    bound &= bind(MCTree, "residue_stop_y", &residue_stop_y);
    bound &= bind(MCTree, "residue_stop_z", &residue_stop_z);
    bound &= bind(MCTree, "residue_stop_strip", &residue_stop_strip);
    bound &= bind(MCTree, "residue_step", &residue_step);
    bound &= bind(MCTree, "residue_termination", &residue_termination);
    if (!bound)
      return nullptr;
  } else {
    tree = new TTree("events_MeV", "Simulated MUSIC events (energies in MeV)");
    tree->Branch("LeftdE", LeftdE, Form("LeftdE[%d]/F", N_SEG_STRIPS));
    tree->Branch("RightdE", RightdE, Form("RightdE[%d]/F", N_SEG_STRIPS));
    tree->Branch("Strip0dE", &Strip0dE, "Strip0dE/F");
    tree->Branch("Strip17dE", &Strip17dE, "Strip17dE/F");
    tree->Branch("Cathode", &Cathode, "Cathode/F");

    MCTree = new TTree("MC", "Truth-level MUSIC simulation");
    // n_steps is the on-disk length of the per-step arrays, so its branch
    // must be defined before any of them.
    MCTree->Branch("n_steps", &n_steps, "n_steps/I");
    MCTree->Branch("requested_strip", &requested_strip, "requested_strip/I");
    MCTree->Branch("reaction_strip", &reaction_strip, "reaction_strip/I");
    MCTree->Branch("event_index", &event_index, "event_index/l");
    MCTree->Branch("beam_energy_accel", &beam_energy_accel,
                   "beam_energy_accel/F");
    MCTree->Branch("beam_energy_gas", &beam_energy_gas, "beam_energy_gas/F");
    MCTree->Branch("beam_energy_reaction", &beam_energy_reaction,
                   "beam_energy_reaction/F");
    MCTree->Branch("beam_energy_exit", &beam_energy_exit, "beam_energy_exit/F");
    MCTree->Branch("beam_stop_x", &beam_stop_x, "beam_stop_x/F");
    MCTree->Branch("beam_stop_y", &beam_stop_y, "beam_stop_y/F");
    MCTree->Branch("beam_stop_z", &beam_stop_z, "beam_stop_z/F");
    MCTree->Branch("beam_stop_strip", &beam_stop_strip, "beam_stop_strip/I");
    MCTree->Branch("beam_termination", &beam_termination, "beam_termination/I");
    MCTree->Branch("vertex_x", &vertex_x, "vertex_x/F");
    MCTree->Branch("vertex_y", &vertex_y, "vertex_y/F");
    MCTree->Branch("vertex_z", &vertex_z, "vertex_z/F");
    MCTree->Branch("DeadUS_dE", &DeadUS_dE, "DeadUS_dE/F");
    MCTree->Branch("DeadDS_dE", &DeadDS_dE, "DeadDS_dE/F");
    MCTree->Branch("evap_energy", evap_energy, "evap_energy[n_steps]/F");
    MCTree->Branch("residue_energy", residue_energy,
                   "residue_energy[n_steps]/F");
    MCTree->Branch("evap_energy_exit", evap_energy_exit,
                   "evap_energy_exit[n_steps]/F");
    MCTree->Branch("residue_energy_exit", residue_energy_exit,
                   "residue_energy_exit[n_steps]/F");
    MCTree->Branch("theta_cm", theta_cm, "theta_cm[n_steps]/F");
    MCTree->Branch("phi_cm", phi_cm, "phi_cm[n_steps]/F");
    MCTree->Branch("evap_theta", evap_theta, "evap_theta[n_steps]/F");
    MCTree->Branch("evap_phi", evap_phi, "evap_phi[n_steps]/F");
    MCTree->Branch("residue_theta", residue_theta, "residue_theta[n_steps]/F");
    MCTree->Branch("residue_phi", residue_phi, "residue_phi[n_steps]/F");
    MCTree->Branch("evap_stop_x", evap_stop_x, "evap_stop_x[n_steps]/F");
    MCTree->Branch("evap_stop_y", evap_stop_y, "evap_stop_y[n_steps]/F");
    MCTree->Branch("evap_stop_z", evap_stop_z, "evap_stop_z[n_steps]/F");
    MCTree->Branch("evap_stop_strip", evap_stop_strip,
                   "evap_stop_strip[n_steps]/I");
    MCTree->Branch("evap_termination", evap_termination,
                   "evap_termination[n_steps]/I");
    MCTree->Branch("residue_stop_x", &residue_stop_x, "residue_stop_x/F");
    MCTree->Branch("residue_stop_y", &residue_stop_y, "residue_stop_y/F");
    MCTree->Branch("residue_stop_z", &residue_stop_z, "residue_stop_z/F");
    MCTree->Branch("residue_stop_strip", &residue_stop_strip,
                   "residue_stop_strip/I");
    MCTree->Branch("residue_step", &residue_step, "residue_step/I");
    MCTree->Branch("residue_termination", &residue_termination,
                   "residue_termination/I");
    // Friended so users can `events_MeV->Draw("beam_energy_reaction:Cathode")`
    // without manually loading MCTree.
    tree->AddFriend(MCTree);
  }
  ResetBranches();
  if (PrintLevel > 0) {
    tree->Print();
    if (MCTree)
      MCTree->Print();
  }
  return tree;
}

void Simulator::ResetBranches() {
  for (Int_t k = 0; k < N_SEG_STRIPS; ++k) {
    LeftdE[k] = RightdE[k] = 0;
  }
  Strip0dE = Strip17dE = 0;
  Cathode = 0;
  n_steps = numEvaporations;
  requested_strip = -1;
  reaction_strip = -1;
  event_index = 0;
  beam_energy_gas = beam_energy_reaction = 0;
  beam_energy_exit = -2.0f; // N/A unless overwritten on unreacted-beam events
  beam_stop_x = beam_stop_y = 0;
  beam_stop_z = -1000;
  beam_stop_strip = -2;
  beam_termination = static_cast<Int_t>(TerminationReason::NotPropagated);
  vertex_x = vertex_y = 0;
  vertex_z = -1000;
  DeadUS_dE = DeadDS_dE = 0.0f;
  for (Int_t er = 0; er < maxEvaporations; er++) {
    // Same -2 = N/A sentinel on creation and exit energies, so unused slots
    // and disallowed-step slots aren't confused with "particle has KE=0".
    evap_energy[er] = residue_energy[er] = -2.0f;
    evap_energy_exit[er] = residue_energy_exit[er] = -2.0f;
    theta_cm[er] = phi_cm[er] = -1;
    evap_theta[er] = evap_phi[er] = -1;
    residue_theta[er] = residue_phi[er] = -1;
    evap_stop_x[er] = evap_stop_y[er] = 0;
    evap_stop_z[er] = -1000;
    evap_stop_strip[er] = -2;
    evap_termination[er] = static_cast<Int_t>(TerminationReason::NotPropagated);
  }
  residue_stop_x = residue_stop_y = 0;
  residue_stop_z = -1000;
  residue_stop_strip = -2;
  residue_step = -1;
  residue_termination = static_cast<Int_t>(TerminationReason::NotPropagated);
}

Bool_t Simulator::WriteRunMetadata(TFile *ROOTfile) {
  if (!ROOTfile || ROOTfile->IsZombie() || !ROOTfile->IsWritable())
    return kFALSE;
  TDirectory *metadata = ROOTfile->GetDirectory("metadata");
  if (!metadata)
    metadata = ROOTfile->mkdir("metadata");
  if (!metadata)
    return kFALSE;
  metadata->cd();

  std::ifstream input(ctrlFilePath_.Data());
  std::ostringstream source;
  if (input)
    source << input.rdbuf();

  auto excitationName = [](Int_t model) {
    return model == 0 ? "forced" : model == 1 ? "ground" : "uniform";
  };
  auto angularName = [](Int_t model) {
    return model == 1 ? "rutherford" : "isotropic";
  };
  std::ostringstream resolved;
  TTree *events = dynamic_cast<TTree *>(ROOTfile->Get("events_MeV"));
  resolved.precision(17);
  resolved << "schema_version = 3\n"
           << "master_seed = " << ctf.Seed << "\n"
           << "events_per_strip = " << ctf.NEvents << "\n"
           << "strip_first = " << ctf.stripFirst << "\n"
           << "strip_last = " << ctf.stripLast << "\n"
           << "total_entries = " << (events ? events->GetEntries() : 0) << "\n"
           << "beam_species = \"" << ctf.beamName << "\"\n"
           << "beam_energy_mev = " << ctf.BeamEnergy << "\n"
           << "gas_species = \"" << ctf.gas << "\"\n"
           << "gas_pressure_torr = " << ctf.pressure << "\n"
           << "gas_temperature_k = " << ctf.temperature << "\n"
           << "straggling = " << (physics_.stragglingEnabled ? "true" : "false")
           << "\n"
           << "stopping_model = \""
           << (physics_.stoppingModel == 0   ? "catima"
               : physics_.stoppingModel == 1 ? "srim"
                                             : "mean")
           << "\"\n";
  for (Int_t step = 0; step < ctf.NumEvapPart; ++step)
    resolved << "reaction_step_" << step << " = \"" << ctf.evap[step] << " + "
             << ctf.res[step] << "; " << excitationName(ctf.residueExc[step])
             << "; " << angularName(ctf.angularDist[step]) << "\"\n";

  const std::string terminationCodes =
      "0=not_propagated, 1=reached_vertex, 2=stopped, "
      "3=downstream_exit, 4=upstream_exit, 5=side_exit, 6=timeout, "
      "7=invalid_state";
  const std::vector<std::pair<const char *, std::string>> records = {
      {"musicsim_version", MUSICSIM_VERSION},
      {"control_file_path", ctrlFilePath_.Data()},
      {"control_file_toml", source.str()},
      {"resolved_configuration", resolved.str()},
      {"termination_codes", terminationCodes},
      {"straggling_sampler",
       "ROOT::Math::VavilovAccurate quantile table; Gaussian for kappa>10"}};
  for (const auto &[name, value] : records) {
    TObjString object(value.c_str());
    if (object.Write(name, TObject::kOverwrite) <= 0) {
      ROOTfile->cd();
      return kFALSE;
    }
  }
  ROOTfile->cd();
  return kTRUE;
}

void Simulator::CreateTracesAndTrajectories() {
  // Pick column colors from the first strip that has a complete set defined;
  // otherwise fall back to a generic palette.
  std::vector<Short_t> Chroma(static_cast<size_t>(AnodeCols));
  for (Int_t col = 0; col < AnodeCols; col++)
    Chroma[col] = static_cast<Short_t>(7 - col);
  for (Int_t stp = 0; stp < AnodeRows; stp++) {
    Int_t NotWhiteColumns = 0;
    for (Int_t col = 0; col < AnodeCols; col++)
      if (AnodeColor[stp][col] != kWhite)
        NotWhiteColumns++;
    if (NotWhiteColumns == AnodeCols) {
      for (Int_t col = 0; col < AnodeCols; col++)
        Chroma[col] = AnodeColor[stp][col];
      break;
    }
  }

  // Detector traces (one per column + a combined trace at index AnodeCols).
  Trace = new TGraph *[AnodeCols + 1];
  for (Int_t col = 0; col < AnodeCols + 1; col++) {
    Trace[col] = new TGraph();
    if (col == AnodeCols) {
      Trace[col]->SetName("full trace");
      Trace[col]->SetLineColor(kBlack);
      Trace[col]->SetLineWidth(2);
    } else {
      Trace[col]->SetName(Form("trace col %d", col));
      Trace[col]->SetLineColor(Chroma[col]);
      Trace[col]->SetLineStyle(2);
      Trace[col]->SetLineWidth(2);
    }
  }

  // Unreacted-beam traces.
  TraceUB = new TGraph *[AnodeCols + 1];
  for (Int_t col = 0; col < AnodeCols + 1; col++) {
    TraceUB[col] = new TGraph();
    if (col == AnodeCols) {
      TraceUB[col]->SetName("Beam trace");
      TraceUB[col]->SetLineColor(kGray);
      TraceUB[col]->SetLineWidth(3);
    } else {
      TraceUB[col]->SetName(Form("Beam trace col %d", col));
      TraceUB[col]->SetLineColor(kGray);
      TraceUB[col]->SetLineStyle(2);
      TraceUB[col]->SetLineWidth(2);
    }
  }

  // Beam traces (for the reacted-beam path).
  TraceB = new TGraph *[AnodeCols + 1];
  for (Int_t col = 0; col < AnodeCols + 1; col++) {
    TraceB[col] = new TGraph();
    if (col == AnodeCols) {
      TraceB[col]->SetName("Beam trace");
      TraceB[col]->SetLineColor(kGray + 2);
      TraceB[col]->SetLineWidth(3);
    } else {
      TraceB[col]->SetName(Form("Beam trace col %d", col));
      TraceB[col]->SetLineColor(kGray);
      TraceB[col]->SetLineStyle(2);
      TraceB[col]->SetLineWidth(2);
    }
  }

  // Evaporation-residue traces.
  TraceER = new TGraph **[numEvaporations];
  for (Int_t er = 0; er < numEvaporations; er++) {
    TraceER[er] = new TGraph *[AnodeCols + 1];
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      TraceER[er][col] = new TGraph();
      if (col == AnodeCols) {
        TraceER[er][col]->SetName(Form("%s trace", EvaR[er]->Name.Data()));
        TraceER[er][col]->SetLineColor(EvaR[er]->GetColor());
        TraceER[er][col]->SetLineWidth(3);
      } else {
        TraceER[er][col]->SetName(
            Form("%s trace col %d", EvaR[er]->Name.Data(), col));
        TraceER[er][col]->SetLineColor(
            static_cast<Color_t>(EvaR[er]->GetColor() - er - 1));
        TraceER[er][col]->SetLineStyle(2);
        TraceER[er][col]->SetLineWidth(2);
      }
    }
  }

  // Evaporated-particle traces.
  TraceEP = new TGraph **[numEvaporations];
  for (Int_t er = 0; er < numEvaporations; er++) {
    TraceEP[er] = new TGraph *[AnodeCols + 1];
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      TraceEP[er][col] = new TGraph();
      if (col == AnodeCols) {
        TraceEP[er][col]->SetName(Form("%s trace", EvaP[er]->Name.Data()));
        TraceEP[er][col]->SetLineColor(EvaP[er]->GetColor());
        TraceEP[er][col]->SetLineWidth(3);
      } else {
        TraceEP[er][col]->SetName(
            Form("%s trace col %d", EvaP[er]->Name.Data(), col));
        TraceEP[er][col]->SetLineColor(
            static_cast<Color_t>(EvaP[er]->GetColor() - er - 1));
        TraceEP[er][col]->SetLineStyle(2);
        TraceEP[er][col]->SetLineWidth(2);
      }
    }
  }
  tracesCreated = true;
}

void Simulator::UpdateVisuals(Int_t evt, Double_t Kbr, Double_t zr,
                              Double_t TOF, Int_t Wait) {
  if (PrintLevel > 0)
    Log << "Update visuals: evt=" << evt << " Kbr=" << Kbr << " MeV zr=" << zr
        << " cm TOF=" << TOF << " ns Wait=" << Wait << "\n3D stuff ..."
        << std::endl;

  if (Wait) {
    Double_t tracklength = TrackBeam->GetVector().Mag();
    TrackBeam->SetTubeR(static_cast<Float_t>(0.1 / tracklength));
    TrackBeam->ElementChanged();

    Short_t C, S, W;
    for (Int_t er = 0; er < numEvaporations; er++) {
      if (EvaP[er] && !EvaP[er]->DoNotPropagate) {
        EvaP[er]->GetTrajectoryAtt(C, S, W);
        tracklength = TrackEvaP[er]->GetVector().Mag();
        if (tracklength > 0) {
          TrackEvaP[er]->SetTubeR(static_cast<Float_t>(0.1 / tracklength));
          TrackEvaP[er]->ElementChanged();
        }
      }
      if (EvaR[er] && !EvaR[er]->DoNotPropagate) {
        tracklength = TrackEvaR[er]->GetVector().Mag();
        if (tracklength > 0) {
          EvaR[er]->GetTrajectoryAtt(C, S, W);
          TrackEvaR[er]->SetTubeR(static_cast<Float_t>(0.1 / tracklength));
          TrackEvaR[er]->ElementChanged();
        }
      }
    }
    Eve->Redraw3D();
  }

  if (PrintLevel > 0)
    Log << "2D stuff..." << std::endl;

  if (tracesCreated) {
    TraceCan->cd(1);
    LabelKine->Draw();
    TraceUB[AnodeCols]->Draw("l same");
    for (Int_t col = 0; col < AnodeCols; col++)
      Trace[col]->Draw("l same");
    Trace[AnodeCols]->Draw("*l same");
    if (LegCol->GetNRows() == 0) {
      LegCol->AddEntry(Trace[AnodeCols], "All columns", "l");
      for (Int_t col = 0; col < AnodeCols; col++)
        LegCol->AddEntry(Trace[col], Form("Column %d", col), "l");
      LegCol->Draw();
    }

    TraceCan->cd(2);
    TraceUB[AnodeCols]->Draw("l same");
    TraceB[AnodeCols]->Draw("l same");
    for (Int_t er = 0; er < numEvaporations; er++) {
      if (TraceER[er][AnodeCols]->GetN() > 0)
        TraceER[er][AnodeCols]->Draw("l same");
      if (TraceEP[er][AnodeCols]->GetN() > 0)
        TraceEP[er][AnodeCols]->Draw("l same");
    }
    Trace[AnodeCols]->Draw("*l same");
    if (LegPart->GetNRows() == 0) {
      LegPart->AddEntry(Trace[AnodeCols], "All particles", "l");
      LegPart->AddEntry(TraceB[AnodeCols], "beam", "l");
      for (Int_t er = 0; er < numEvaporations; er++) {
        if (TraceEP[er][AnodeCols]->GetN() > 0)
          LegPart->AddEntry(TraceEP[er][AnodeCols], EvaP[er]->Name.Data(), "l");
        if (TraceER[er][AnodeCols]->GetN() > 0)
          LegPart->AddEntry(TraceER[er][AnodeCols], EvaR[er]->Name.Data(), "l");
      }
    }
    LegPart->Draw();

    TraceCan->Update();
    if (Wait == 1)
      TraceCan->WaitPrimitive();
  }
}
