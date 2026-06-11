#include "Simulator.hpp"

// Initialize the "events_MeV" tree (detector-level output; branch layout
// mirrors the experimental "events" tree from the MUSIC EventBuilder) and the
// friended "MC" tree (truth-only).
// Energies are Float_t in MeV (the experimental data file uses ADC counts);
// analysis is expected to apply per-channel calibration to compare data
// against sim in MeV.
TTree *Simulator::InitTree(TFile *ROOTfile, TString FileOpt) {
  TTree *tree;
  const Bool_t update = (FileOpt == "update" || FileOpt == "UPDATE");
  if (ROOTfile && update) {
    tree = (TTree *)ROOTfile->Get("events_MeV");
    tree->SetBranchAddress("Left_0_17_dE", Left_0_17_dE);
    tree->SetBranchAddress("RightdE", RightdE);
    tree->SetBranchAddress("Cathode", &Cathode);
    MCTree = (TTree *)ROOTfile->Get("MC");
    MCTree->SetBranchAddress("n_steps", &n_steps);
    MCTree->SetBranchAddress("reaction_strip", &reaction_strip);
    MCTree->SetBranchAddress("beam_energy_accel", &beam_energy_accel);
    MCTree->SetBranchAddress("beam_energy_gas", &beam_energy_gas);
    MCTree->SetBranchAddress("beam_energy_reaction", &beam_energy_reaction);
    MCTree->SetBranchAddress("beam_energy_exit", &beam_energy_exit);
    MCTree->SetBranchAddress("beam_stop_x", &beam_stop_x);
    MCTree->SetBranchAddress("beam_stop_y", &beam_stop_y);
    MCTree->SetBranchAddress("beam_stop_z", &beam_stop_z);
    MCTree->SetBranchAddress("beam_stop_strip", &beam_stop_strip);
    MCTree->SetBranchAddress("vertex_x", &vertex_x);
    MCTree->SetBranchAddress("vertex_y", &vertex_y);
    MCTree->SetBranchAddress("vertex_z", &vertex_z);
    MCTree->SetBranchAddress("DeadUS_dE", &DeadUS_dE);
    MCTree->SetBranchAddress("DeadDS_dE", &DeadDS_dE);
    MCTree->SetBranchAddress("evap_energy", evap_energy);
    MCTree->SetBranchAddress("residue_energy", residue_energy);
    MCTree->SetBranchAddress("evap_energy_exit", evap_energy_exit);
    MCTree->SetBranchAddress("residue_energy_exit", residue_energy_exit);
    MCTree->SetBranchAddress("theta_cm", theta_cm);
    MCTree->SetBranchAddress("phi_cm", phi_cm);
    MCTree->SetBranchAddress("evap_theta", evap_theta);
    MCTree->SetBranchAddress("evap_phi", evap_phi);
    MCTree->SetBranchAddress("residue_theta", residue_theta);
    MCTree->SetBranchAddress("residue_phi", residue_phi);
    MCTree->SetBranchAddress("evap_stop_x", evap_stop_x);
    MCTree->SetBranchAddress("evap_stop_y", evap_stop_y);
    MCTree->SetBranchAddress("evap_stop_z", evap_stop_z);
    MCTree->SetBranchAddress("evap_stop_strip", evap_stop_strip);
    MCTree->SetBranchAddress("residue_stop_x", &residue_stop_x);
    MCTree->SetBranchAddress("residue_stop_y", &residue_stop_y);
    MCTree->SetBranchAddress("residue_stop_z", &residue_stop_z);
    MCTree->SetBranchAddress("residue_stop_strip", &residue_stop_strip);
    MCTree->SetBranchAddress("residue_step", &residue_step);
  } else {
    tree = new TTree("events_MeV", "Simulated MUSIC events (energies in MeV)");
    tree->Branch("Left_0_17_dE", Left_0_17_dE,
                 Form("Left_0_17_dE[%d]/F", N_STRIPS));
    tree->Branch("RightdE", RightdE, Form("RightdE[%d]/F", N_STRIPS));
    tree->Branch("Cathode", &Cathode, "Cathode/F");

    MCTree = new TTree("MC", "Truth-level MUSIC simulation");
    // n_steps is the on-disk length of the per-step arrays, so its branch
    // must be defined before any of them.
    MCTree->Branch("n_steps", &n_steps, "n_steps/I");
    MCTree->Branch("reaction_strip", &reaction_strip, "reaction_strip/I");
    MCTree->Branch("beam_energy_accel", &beam_energy_accel,
                   "beam_energy_accel/F");
    MCTree->Branch("beam_energy_gas", &beam_energy_gas, "beam_energy_gas/F");
    MCTree->Branch("beam_energy_reaction", &beam_energy_reaction,
                   "beam_energy_reaction/F");
    MCTree->Branch("beam_energy_exit", &beam_energy_exit,
                   "beam_energy_exit/F");
    MCTree->Branch("beam_stop_x", &beam_stop_x, "beam_stop_x/F");
    MCTree->Branch("beam_stop_y", &beam_stop_y, "beam_stop_y/F");
    MCTree->Branch("beam_stop_z", &beam_stop_z, "beam_stop_z/F");
    MCTree->Branch("beam_stop_strip", &beam_stop_strip, "beam_stop_strip/I");
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
    MCTree->Branch("residue_stop_x", &residue_stop_x, "residue_stop_x/F");
    MCTree->Branch("residue_stop_y", &residue_stop_y, "residue_stop_y/F");
    MCTree->Branch("residue_stop_z", &residue_stop_z, "residue_stop_z/F");
    MCTree->Branch("residue_stop_strip", &residue_stop_strip,
                   "residue_stop_strip/I");
    MCTree->Branch("residue_step", &residue_step, "residue_step/I");
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
  for (Int_t s = 0; s < N_STRIPS; ++s) {
    Left_0_17_dE[s] = RightdE[s] = 0;
  }
  Cathode = 0;
  n_steps = numEvaporations;
  reaction_strip = -1;
  beam_energy_gas = beam_energy_reaction = 0;
  beam_energy_exit = -2.0f; // N/A unless overwritten on unreacted-beam events
  beam_stop_x = beam_stop_y = 0;
  beam_stop_z = -1000;
  beam_stop_strip = -2;
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
  }
  residue_stop_x = residue_stop_y = 0;
  residue_stop_z = -1000;
  residue_stop_strip = -2;
  residue_step = -1;
}

void Simulator::CreateTracesAndTrajectories() {
  // Pick column colors from the first strip that has a complete set defined;
  // otherwise fall back to a generic palette.
  Short_t *Chroma = new Short_t[AnodeCols];
  for (Int_t col = 0; col < AnodeCols; col++)
    Chroma[col] = 7 - col;
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
        TraceER[er][col]->SetLineColor(EvaR[er]->GetColor() - er - 1);
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
        TraceEP[er][col]->SetLineColor(EvaP[er]->GetColor() - er - 1);
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
    TrackBeam->SetTubeR(0.1 / tracklength);
    TrackBeam->ElementChanged();

    Short_t C, S, W;
    for (Int_t er = 0; er < numEvaporations; er++) {
      if (EvaP[er] && !EvaP[er]->DoNotPropagate) {
        EvaP[er]->GetTrajectoryAtt(C, S, W);
        tracklength = TrackEvaP[er]->GetVector().Mag();
        if (tracklength > 0) {
          TrackEvaP[er]->SetTubeR(0.1 / tracklength);
          TrackEvaP[er]->ElementChanged();
        }
      }
      if (EvaR[er] && !EvaR[er]->DoNotPropagate) {
        tracklength = TrackEvaR[er]->GetVector().Mag();
        if (tracklength > 0) {
          EvaR[er]->GetTrajectoryAtt(C, S, W);
          TrackEvaR[er]->SetTubeR(0.1 / tracklength);
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

void Simulator::WriteTraces(char * /*FileName*/) {
  // Vestigial — traces are written per-event by ComputeDetectorResponse. Kept
  // as a no-op so the public API stays source-compatible with downstream code.
}
