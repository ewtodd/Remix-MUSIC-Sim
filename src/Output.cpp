#include "Simulator.hpp"

// Initialize the "events_MeV" tree (detector-level output, layout matches the
// upstream EventBuilderNearestGrid) and the friended "MC" tree (truth-only).
// Energies are Float_t in MeV (the upstream data file uses Int_t ADC counts);
// analysis is expected to apply per-channel calibration to compare data
// against sim in MeV.
TTree *Simulator::InitTree(TFile *ROOTfile, TString FileOpt) {
  TTree *tree;
  const Bool_t update = (FileOpt == "update" || FileOpt == "UPDATE");
  if (ROOTfile && update) {
    tree = (TTree *)ROOTfile->Get("events_MeV");
    tree->SetBranchAddress("LeftdE", LeftdE);
    tree->SetBranchAddress("RightdE", RightdE);
    tree->SetBranchAddress("TotaldE", TotaldE);
    tree->SetBranchAddress("Cathode", &Cathode);
    MCTree = (TTree *)ROOTfile->Get("MC");
    MCTree->SetBranchAddress("reacStp", &reacStp);
    MCTree->SetBranchAddress("BeamEnergyAccel", &BeamEnergyAccel);
    MCTree->SetBranchAddress("Kbi", &Kbi);
    MCTree->SetBranchAddress("Kbr", &Kbr);
    MCTree->SetBranchAddress("Kbeam_exit", &Kbeam_exit);
    MCTree->SetBranchAddress("DeadUS_dE", &DeadUS_dE);
    MCTree->SetBranchAddress("DeadDS_dE", &DeadDS_dE);
    MCTree->SetBranchAddress("Kl", Kl);
    MCTree->SetBranchAddress("Kh", Kh);
    MCTree->SetBranchAddress("Kl_exit", Kl_exit);
    MCTree->SetBranchAddress("Kh_exit", Kh_exit);
    MCTree->SetBranchAddress("theta_CM", theta_CM);
    MCTree->SetBranchAddress("theta_l", theta_l);
    MCTree->SetBranchAddress("theta_h", theta_h);
    MCTree->SetBranchAddress("phi_l", phi_l);
    MCTree->SetBranchAddress("phi_h", phi_h);
    MCTree->SetBranchAddress("xfl", xfl);
    MCTree->SetBranchAddress("yfl", yfl);
    MCTree->SetBranchAddress("zfl", zfl);
    MCTree->SetBranchAddress("xr", &xr);
    MCTree->SetBranchAddress("yr", &yr);
    MCTree->SetBranchAddress("zr", &zr);
    MCTree->SetBranchAddress("xfe", &xfe);
    MCTree->SetBranchAddress("yfe", &yfe);
    MCTree->SetBranchAddress("zfe", &zfe);
    MCTree->SetBranchAddress("resID", &resID);
  } else {
    tree = new TTree("events_MeV", "Simulated MUSIC events (energies in MeV)");
    tree->Branch("LeftdE", LeftdE, Form("LeftdE[%d]/F", N_STRIPS));
    tree->Branch("RightdE", RightdE, Form("RightdE[%d]/F", N_STRIPS));
    tree->Branch("TotaldE", TotaldE, Form("TotaldE[%d]/F", N_STRIPS));
    tree->Branch("Cathode", &Cathode, "Cathode/F");

    MCTree = new TTree("MC", "Truth-level MUSIC simulation");
    MCTree->Branch("reacStp", &reacStp, "reacStp/I");
    MCTree->Branch("BeamEnergyAccel", &BeamEnergyAccel, "BeamEnergyAccel/F");
    MCTree->Branch("Kbi", &Kbi, "Kbi/F");
    MCTree->Branch("Kbr", &Kbr, "Kbr/F");
    MCTree->Branch("Kbeam_exit", &Kbeam_exit, "Kbeam_exit/F");
    MCTree->Branch("DeadUS_dE", &DeadUS_dE, "DeadUS_dE/F");
    MCTree->Branch("DeadDS_dE", &DeadDS_dE, "DeadDS_dE/F");
    MCTree->Branch("Kl", Kl, Form("Kl[%d]/F", maxEvaporations));
    MCTree->Branch("Kh", Kh, Form("Kh[%d]/F", maxEvaporations));
    MCTree->Branch("Kl_exit", Kl_exit, Form("Kl_exit[%d]/F", maxEvaporations));
    MCTree->Branch("Kh_exit", Kh_exit, Form("Kh_exit[%d]/F", maxEvaporations));
    MCTree->Branch("theta_CM", theta_CM,
                   Form("theta_CM[%d]/F", maxEvaporations));
    MCTree->Branch("theta_l", theta_l, Form("theta_l[%d]/F", maxEvaporations));
    MCTree->Branch("theta_h", theta_h, Form("theta_h[%d]/F", maxEvaporations));
    MCTree->Branch("phi_l", phi_l, Form("phi_l[%d]/F", maxEvaporations));
    MCTree->Branch("phi_h", phi_h, Form("phi_h[%d]/F", maxEvaporations));
    MCTree->Branch("xfl", xfl, Form("xfl[%d]/F", maxEvaporations));
    MCTree->Branch("yfl", yfl, Form("yfl[%d]/F", maxEvaporations));
    MCTree->Branch("zfl", zfl, Form("zfl[%d]/F", maxEvaporations));
    MCTree->Branch("xr", &xr, "xr/F");
    MCTree->Branch("yr", &yr, "yr/F");
    MCTree->Branch("zr", &zr, "zr/F");
    MCTree->Branch("xfe", &xfe, "xfe/F");
    MCTree->Branch("yfe", &yfe, "yfe/F");
    MCTree->Branch("zfe", &zfe, "zfe/F");
    MCTree->Branch("resID", &resID, "resID/I");
    // Friended so users can `events_MeV->Draw("Kbr:Cathode")` without manually
    // loading MCTree.
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
    LeftdE[s] = RightdE[s] = TotaldE[s] = 0;
  }
  Cathode = 0;
  reacStp = -1;
  Kbi = Kbr = 0;
  Kbeam_exit = -2.0f; // N/A unless overwritten on unreacted-beam events
  DeadUS_dE = DeadDS_dE = 0.0f;
  for (Int_t er = 0; er < maxEvaporations; er++) {
    // Same -2 = N/A sentinel as Kl_exit/Kh_exit, so unused slots and
    // disallowed-step slots aren't confused with "outgoing particle has KE=0".
    Kl[er] = Kh[er] = -2.0f;
    Kl_exit[er] = Kh_exit[er] = -2.0f;
    phi_CM[er] = theta_CM[er] = -1;
    phi_l[er] = theta_l[er] = -1;
    phi_h[er] = theta_h[er] = -1;
    xfl[er] = yfl[er] = 0;
    zfl[er] = -1000;
  }
  xr = yr = 0;
  zr = -1000;
  xfe = yfe = 0;
  zfe = -1000;
  resID = -1;
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
