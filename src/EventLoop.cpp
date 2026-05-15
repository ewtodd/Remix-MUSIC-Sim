#include "Simulator.hpp"

// Main per-event simulation. For each event:
//   1. sample beam initial energy (accelerator FWHM → degrader → entrance
//   window)
//   2. pick a reaction vertex uniformly inside StpID (StpID == -1: unreacted
//   beam)
//   3. forward-propagate the beam to the vertex (per-step Vavilov straggling)
//   4. set reaction kinematics at the vertex
//   5. propagate outgoing particles
//   6. if reaction not allowed, continue propagating the beam to the exit
//   7. compute detector response and fill the tree
//   8. update visuals (interactive mode only)
void Simulator::Simulate(Int_t StpID, Int_t NEvents, Double_t MaxTime,
                         Double_t UserStep, Int_t UpdateVis, Int_t Wait,
                         TFile * /*ROOTfile*/) {
  if (VolAnode == 0) {
    std::cout << "Anode geometry not specified. Use SetAnode method."
              << std::endl;
    return;
  }

  this->NEvents = NEvents;

  TStopwatch StpWatch;
  LongDouble_t Frac[11] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.5,
                           0.6,  0.7,  0.8, 0.9, 1.0};
  Int_t FIndex = 0;

  if (verbose_)
    std::cout << "Simulating " << NEvents << " MUSIC traces for strip " << StpID
              << " ... " << std::endl;

  // Idealized beam (no straggling) used for the kinematic-range estimate and
  // the unreacted-beam reference trace.
  SetInitialKinematics(Kb_at_gas);

  Particle *BeamInit = new Particle("beam init");
  BeamInit->Copy(Beam);
  Particle *BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel > 0)
    BeamCopy->Print();
  PropagateParticle(BeamCopy, 0, MaxTime, UserStep, DeltaEB_ave);
  if (tracesCreated)
    for (Int_t stp = 0; stp < AnodeRows; stp++)
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);

  // Beam-energy bounds in the chosen strip (beam parallel to z).
  Double_t Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (Int_t stp = 1; stp < AnodeRows; stp++) {
    MinZ += AnodeDZ[stp - 1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0] == StpID)
      break;
  }
  Kb_max = Beam->GetFinalEnergy(0, Kb_at_gas, MinZ);
  MinT = Beam->GetTimeOfFlight(0);
  Kb_min = Beam->GetFinalEnergy(0, Kb_at_gas, MaxZ);
  MaxT = Beam->GetTimeOfFlight(0);
  if (PrintLevel > 0) {
    Log << "|---- Kinematic constraints for strip ";
    Log.width(3);
    Log << StpID;
    Log << " ---------\n"
        << "|     |   In    |   Out   |  Units  |\n"
        << "| zr  |";
    Log.width(9);
    Log << MinZ;
    Log << "|";
    Log.width(9);
    Log << MaxZ;
    Log << "|";
    Log.width(9);
    Log << "cm";
    Log << "|\n";
    Log << "| tof |";
    Log.width(9);
    Log << MinT;
    Log << "|";
    Log.width(9);
    Log << MaxT;
    Log << "|";
    Log.width(9);
    Log << "ns";
    Log << "|\n";
    Log << "| Kb  |";
    Log.width(9);
    Log << Kb_max;
    Log << "|";
    Log.width(9);
    Log << Kb_min;
    Log << "|";
    Log.width(9);
    Log << "MeV";
    Log << "|\n";
    Log << "|--------------------------------------------------" << std::endl;
  }

  CalculateCMEnergyRange();

  Log << "Initiating event for-loop" << std::endl;
  if (verbose_)
    std::cout << "\nInitiating event for-loop" << std::endl;

  Double_t ti, xi, yi, zi, tf, xf, yf, zf;
  for (Int_t evt = 0; evt < NEvents; evt++) {
    if (evt % 1000 == 0 && CheckMemoryUsage() == 0) {
      Log << "Exiting musicsim (memory limit exceeded)." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (PrintLevel > 0) {
      Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
          << "!!!       EVENT " << evt
          << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
          << std::endl;
    }
    ResetBranches();

    if (UpdateVis) {
      TraceCan->cd(1);
      HCT->Draw();
      TraceCan->cd(2);
      HPT->Draw();
    }

    // Reset per-event detector response.
    for (Int_t stp = 0; stp < AnodeRows; stp++)
      for (Int_t col = 0; col < AnodeCols + 1; col++) {
        DeltaEB[stp][col] = 0;
        DeltaEL[stp][col] = 0;
        DeltaEH[stp][col] = 0;
        for (Int_t er = 0; er < numEvaporations; er++) {
          DeltaE_EvaP[er][stp][col] = 0;
          DeltaE_EvaR[er][stp][col] = 0;
        }
      }

    // Per-event beam chain (physically correct order):
    //   accelerator KE [+ KbFWHM]  ->  degrader  ->  entrance window  ->  Kbi.
    Double_t Ebeam = ctf.BeamEnergy;
    if (ctf.KbFWHM > 0.0)
      Ebeam += Rdm->Gaus(0.0, ctf.KbFWHM / 2.355);
    {
      const Double_t amu_MeV = 931.49410242;
      Int_t A_beam =
          (Beam->Mass > 0) ? Int_t(std::round(Beam->Mass / amu_MeV)) : 0;
      if (hasDegrader_)
        Ebeam = EnergyThroughWithStraggling(A_beam, Beam->Z, Ebeam, degrader_);
      if (entranceWindowEnabled_)
        Ebeam = EnergyThroughWithStraggling(A_beam, Beam->Z, Ebeam,
                                            entranceWindow_);
    }
    this->Kbi = Ebeam;
    SetInitialKinematics(Ebeam);

    Int_t ReacAllowed = 0;
    this->Kbr = 0;
    Double_t TOF = 0;
    if (StpID > -1) {
      // Pick the vertex uniformly inside the chosen strip; forward-propagate
      // the beam to z=zr with per-step straggling so Kbr / TOF inherit the
      // sampled fluctuations (the reaction kinematics see the actually-
      // sampled energy at the vertex, not the deterministic mean).
      this->zr = Rdm->Uniform(MinZ, MaxZ);
      Beam->GetX(ti, xi, yi, zi); // entrance origin (kept for visualization)
      PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB,
                        /*endZ=*/this->zr);
      Double_t t_at_vertex, x_at_vertex, y_at_vertex, z_at_vertex;
      Beam->GetX(t_at_vertex, x_at_vertex, y_at_vertex, z_at_vertex);
      this->Kbr = (Float_t)Beam->GetKE();
      TOF = t_at_vertex;
      if (tracesCreated)
        for (Int_t stp = 0; stp < AnodeRows; stp++)
          for (Int_t col = 0; col < AnodeCols + 1; col++)
            TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      if (UpdateVis) {
        TrackBeam->SetOrigin(xi, yi, zi);
        TrackBeam->SetVector(x_at_vertex - xi, y_at_vertex - yi,
                             z_at_vertex - zi);
      }
      if (PrintLevel > 0)
        Log << "Kbr = " << this->Kbr << "  zr = " << this->zr
            << "  tof = " << TOF << std::endl;

      // SetReactionKinematics overwrites Beam->X / Beam->P to put the beam
      // at (TOF, 0, 0, zr) with KE=Kbr; the propagated state we already
      // captured into DeltaEB / TrackBeam is preserved.
      ReacAllowed = SetReactionKinematics(this->Kbr, this->zr, TOF);
      if (PrintLevel > 0) {
        Log << "Conservation of 4-momentum at reaction point (zr)" << std::endl;
        FourVector Pi("initial 4-momentum (lab)", 0, 0, 0, 0);
        Pi += Beam->GetP() + Target->GetP();
        FourVector Pf("final 4-momentum (lab)", 0, 0, 0, 0);
        for (Int_t er = 0; er < numEvaporations; er++) {
          if (!EvaR[er]->DoNotPropagate)
            Pf += EvaR[er]->GetP();
          if (!EvaP[er]->DoNotPropagate)
            Pf += EvaP[er]->GetP();
        }
        Pi.Print(Log);
        Pf.Print(Log);
      }
    } else {
      if (UpdateVis) {
        LabelKine->Clear();
        LabelKine->AddText("Kinematics");
        LabelKine->AddText(Form("beam: K=%.2f MeV", Ebeam));
      }
    }

    if (ReacAllowed) {
      // Beam was propagated entrance→vertex above. Propagate outgoing
      // particles next.
      for (Int_t er = 0; er < numEvaporations; er++) {
        EvaP[er]->GetX(ti, xi, yi, zi);
        PropagateParticle(EvaP[er], evt, MaxTime, UserStep, DeltaE_EvaP[er]);
        if (tracesCreated)
          for (Int_t stp = 0; stp < AnodeRows; stp++)
            for (Int_t col = 0; col < AnodeCols + 1; col++)
              TraceEP[er][col]->SetPoint(stp, stp, DeltaE_EvaP[er][stp][col]);
        EvaP[er]->GetX(tf, xf, yf, zf);
        if (UpdateVis) {
          TrackEvaP[er]->SetOrigin(xi, yi, zi);
          TrackEvaP[er]->SetVector(xf - xi, yf - yi, zf - zi);
        }
        xfl[er] = xf;
        yfl[er] = yf;
        zfl[er] = zf;

        EvaR[er]->GetX(ti, xi, yi, zi);
        PropagateParticle(EvaR[er], evt, MaxTime, UserStep, DeltaE_EvaR[er]);
        if (tracesCreated)
          for (Int_t stp = 0; stp < AnodeRows; stp++)
            for (Int_t col = 0; col < AnodeCols + 1; col++)
              TraceER[er][col]->SetPoint(stp, stp, DeltaE_EvaR[er][stp][col]);
        EvaR[er]->GetX(tf, xf, yf, zf);
        if (!EvaR[er]->DoNotPropagate) {
          xfe = xf;
          yfe = yf;
          zfe = zf;
          resID = er;
        }
        if (UpdateVis) {
          TrackEvaR[er]->SetOrigin(xi, yi, zi);
          TrackEvaR[er]->SetVector(xf - xi, yf - yi, zf - zi);
        }
      }
    } else {
      // Reaction forbidden, or StpID == -1 (unreacted beam). Sweep the beam
      // from its current position to the exit:
      //   StpID == -1: Beam is at the entrance (no vertex pick); reset_DE=true.
      //   StpID >  -1, forbidden: Beam was swept entrance→zr above (DeltaEB
      //     partially filled), then SetReactionKinematics overwrote its
      //     position to (TOF, 0, 0, zr). Continue forward with reset_DE=false
      //     so the entrance→zr deposits are kept.
      const Bool_t resume = (StpID > -1);
      Beam->GetX(ti, xi, yi, zi);
      PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB,
                        /*endZ=*/-1.0, /*reset_DE=*/!resume);
      Beam->GetX(tf, xf, yf, zf);
      if (UpdateVis) {
        TrackBeam->SetOrigin(xi, yi, zi);
        TrackBeam->SetVector(xf - xi, yf - yi, zf - zi);
      }
      if (tracesCreated)
        for (Int_t stp = 0; stp < AnodeRows; stp++)
          for (Int_t col = 0; col < AnodeCols + 1; col++)
            TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      if (StpID > -1) {
        std::cout << "Warning: reaction energetically not allowed for event "
                  << evt << " (Kbr= " << this->Kbr << " MeV)." << std::endl;
        Log << "Warning: reaction energetically not allowed for event " << evt
            << " (Kbr= " << this->Kbr << " MeV)." << std::endl;
      }
    }

    ComputeDetectorResponse(evt, StpID, UpdateVis);
    Log << "\tDone computing detector response." << std::endl;
    if (SimTree != 0) {
      FinalizeEvent(evt);
      SimTree->Fill();
      if (MCTree)
        MCTree->Fill();
    }
    if (UpdateVis)
      UpdateVisuals(evt, this->Kbr, this->zr, TOF, Wait);

    if (NEvents > 99 && (LongDouble_t)(evt) >= Frac[FIndex] * NEvents) {
      if (verbose_)
        std::cout << "\t" << Frac[FIndex] * 100 << "% processed ("
                  << StpWatch.RealTime() << " s)" << std::endl;
      StpWatch.Start(kFALSE);
      FIndex++;
    }
    NTraces++;
  }
  StpWatch.Stop();
  if (verbose_)
    StpWatch.Print();
  if (verbose_)
    std::cout << "Event for-loop concluded." << std::endl;
  CheckMemoryUsage(1);
}

// Angular-scan variant: instead of sampling random angles per event, walks a
// θ_CM × φ_CM grid for each anode strip.
void Simulator::Simulate(Int_t StpID, Double_t ThCMMin, Double_t ThCMMax,
                         Int_t ThSteps, Double_t PhiCMMin, Double_t PhiCMMax,
                         Int_t PhiSteps, Double_t MaxTime, Double_t UserStep,
                         Int_t Wait) {
  if (VolAnode == 0) {
    std::cout << "Anode geometry not specified. Use SetAnode method."
              << std::endl;
    return;
  }

  TFile *ROOTfile = new TFile("sim.root", "recreate");
  SimTree = InitTree(ROOTfile, "recreate");

  NEvents = PhiSteps * ThSteps * AnodeRows;

  CreateTracesAndTrajectories();
  if (verbose_)
    std::cout << "Simulating MUSIC traces ... " << std::endl;

  SetInitialKinematics(Kb_at_gas);

  Particle *BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  if (PrintLevel > 0)
    BeamCopy->Print();
  PropagateParticle(BeamCopy, Kb_at_gas, MaxTime, UserStep, DeltaEB_ave);
  for (Int_t stp = 0; stp < AnodeRows; stp++)
    for (Int_t col = 0; col < AnodeCols + 1; col++)
      TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);

  Double_t Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  Int_t evt = 0;
  Double_t theta_min = ThCMMin * pi / 180;
  Double_t theta_max = ThCMMax * pi / 180;
  Double_t phi_min = PhiCMMin * pi / 180;
  Double_t phi_max = PhiCMMax * pi / 180;

  MinZ = 0;
  MaxZ = AnodeDZ[0][0];
  for (Int_t stp = 1; stp < AnodeRows; stp++) {
    MinZ += AnodeDZ[stp - 1][0];
    MaxZ += AnodeDZ[stp][0];
    if (AnodeStpID[stp][0] == StpID)
      break;
  }
  Kb_max = Beam->GetFinalEnergy(0, Kb_at_gas, MinZ);
  MinT = Beam->GetTimeOfFlight(0);
  Kb_min = Beam->GetFinalEnergy(0, Kb_at_gas, MaxZ);
  MaxT = Beam->GetTimeOfFlight(0);
  if (PrintLevel > 0) {
    Log << "|---- Kinematic constraints for strip ";
    Log.width(3);
    Log << StpID;
    Log << " ---------\n"
        << "|     |   In    |   Out   |  Units  |\n"
        << "| zr  |";
    Log.width(9);
    Log << MinZ;
    Log << "|";
    Log.width(9);
    Log << MaxZ;
    Log << "|";
    Log.width(9);
    Log << "cm";
    Log << "|\n";
    Log << "| tof |";
    Log.width(9);
    Log << MinT;
    Log << "|";
    Log.width(9);
    Log << MaxT;
    Log << "|";
    Log.width(9);
    Log << "ns";
    Log << "|\n";
    Log << "| Kb  |";
    Log.width(9);
    Log << Kb_max;
    Log << "|";
    Log.width(9);
    Log << Kb_min;
    Log << "|";
    Log.width(9);
    Log << "MeV";
    Log << "|\n";
    Log << "|--------------------------------------------------" << std::endl;
  }

  Double_t theta, phi;
  for (Int_t ths = 0; ths < ThSteps; ths++) {
    theta = ths * (theta_max - theta_min) / (ThSteps - 1) + theta_min;
    for (Int_t phs = 0; phs < PhiSteps; phs++) {
      phi = phs * (phi_max - phi_min) / (PhiSteps - 1) + phi_min;

      if (PrintLevel > 0)
        Log << "\n***************** Event " << evt << "\n" << std::endl;

      TraceCan->cd(1);
      HCT->Draw();

      for (Int_t stp = 0; stp < AnodeRows; stp++)
        for (Int_t col = 0; col < AnodeCols + 1; col++) {
          DeltaEB[stp][col] = DeltaEL[stp][col] = DeltaEH[stp][col] = 0;
        }

      SetInitialKinematics(Kb_at_gas);

      this->zr = Rdm->Uniform(MinZ, MaxZ);
      this->Kbr = Beam->GetFinalEnergy(0, Kb_at_gas, this->zr);
      Double_t TOF = Beam->GetTimeOfFlight(0);
      if (PrintLevel > 0)
        Log << "Kbr = " << this->Kbr << "  zr = " << this->zr
            << "  tof = " << TOF << std::endl;

      Int_t ReacAllowed =
          SetReactionKinematics(this->Kbr, this->zr, TOF, theta, phi);
      if (ReacAllowed == 0) {
        std::cout << "Warning: reaction energetically not allowed for event "
                  << evt << " (Kbr= " << this->Kbr << " MeV)." << std::endl;
        continue;
      }

      // Beam back-propagated to entrance.
      PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
      PropagateParticle(Heavy, evt, MaxTime, UserStep, DeltaEH);
      PropagateParticle(Light, evt, MaxTime, UserStep, DeltaEL);

      ComputeDetectorResponse(evt, StpID, 1);
      FinalizeEvent(evt);
      SimTree->Fill();
      if (MCTree)
        MCTree->Fill();

      UpdateVisuals(evt, this->Kbr, this->zr, TOF, Wait);

      NTraces++;
      evt++;
    }
  }

  if (ROOTfile) {
    ROOTfile->cd();
    SimTree->Write("", TObject::kSingleKey);
    if (MCTree)
      MCTree->Write("", TObject::kSingleKey);
    ROOTfile->Close();
  }
}

// Trace-database generator: like the second Simulate() but writes its own file
// and iterates over every readout strip.
void Simulator::GenerateTraceDatabase(TString FileName, Double_t ThCMMin,
                                      Double_t ThCMMax, Int_t ThSteps,
                                      Double_t PhiCMMin, Double_t PhiCMMax,
                                      Int_t PhiSteps, Double_t MaxTime,
                                      Double_t UserStep, Int_t UpdateVis,
                                      Int_t Wait) {
  Double_t ti, xi, yi, zi, tf, xf, yf, zf;
  if (VolAnode == 0) {
    std::cout << "Anode geometry not specified. Use SetAnode method."
              << std::endl;
    return;
  }

  TStopwatch StpWatch;
  Long_t EvtsProcessed = 0;
  LongDouble_t Frac[6] = {0.01, 0.25, 0.5, 0.75, 0.9, 1.0};
  Int_t FIndex = 0;

  TFile *TDB = new TFile(FileName.Data(), "recreate");
  SimTree = InitTree(TDB, "recreate");

  NEvents = PhiSteps * ThSteps * AnodeRows;
  NTraces = 0;

  CreateTracesAndTrajectories();
  if (verbose_)
    std::cout << "Generating " << NEvents << " MUSIC traces ..." << std::endl;

  SetInitialKinematics(Kb_at_gas);

  Particle *BeamInit = new Particle("beam init");
  BeamInit->Copy(Beam);
  Particle *BeamCopy = new Particle("beam copy");
  BeamCopy->Copy(Beam);
  PropagateParticle(BeamCopy, 0, MaxTime, UserStep, DeltaEB_ave);
  if (tracesCreated)
    for (Int_t stp = 0; stp < AnodeRows; stp++)
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        TraceUB[col]->SetPoint(stp, stp, DeltaEB_ave[stp][col]);

  Double_t Kb_min, Kb_max, MinZ, MaxZ, MinT, MaxT;
  Int_t evt = 0;
  Double_t theta_min = ThCMMin * pi / 180;
  Double_t theta_max = ThCMMax * pi / 180;
  Double_t phi_min = PhiCMMin * pi / 180;
  Double_t phi_max = PhiCMMax * pi / 180;

  for (Int_t stp_base = 0; stp_base < AnodeRows; stp_base++) {
    if (AnodeStpID[stp_base][0] < 0)
      continue;

    MinZ = 0;
    MaxZ = AnodeDZ[0][0];
    for (Int_t stp = 1; stp < AnodeRows; stp++) {
      MinZ += AnodeDZ[stp - 1][0];
      MaxZ += AnodeDZ[stp][0];
      if (AnodeStpID[stp][0] == AnodeStpID[stp_base][0])
        break;
    }
    Kb_max = Beam->GetFinalEnergy(0, Kb_at_gas, MinZ);
    MinT = Beam->GetTimeOfFlight(0);
    Kb_min = Beam->GetFinalEnergy(0, Kb_at_gas, MaxZ);
    MaxT = Beam->GetTimeOfFlight(0);
    if (PrintLevel > 0) {
      Log << "|---- Kinematic constraints for strip ";
      Log.width(3);
      Log << AnodeStpID[stp_base][0];
      Log << " ---------\n"
          << "|     |   In    |   Out   |  Units  |\n"
          << "| zr  |";
      Log.width(9);
      Log << MinZ;
      Log << "|";
      Log.width(9);
      Log << MaxZ;
      Log << "|";
      Log.width(9);
      Log << "cm";
      Log << "|\n";
      Log << "| tof |";
      Log.width(9);
      Log << MinT;
      Log << "|";
      Log.width(9);
      Log << MaxT;
      Log << "|";
      Log.width(9);
      Log << "ns";
      Log << "|\n";
      Log << "| Kb  |";
      Log.width(9);
      Log << Kb_max;
      Log << "|";
      Log.width(9);
      Log << Kb_min;
      Log << "|";
      Log.width(9);
      Log << "MeV";
      Log << "|\n";
      Log << "|--------------------------------------------------" << std::endl;
    }
    Double_t theta, phi;
    for (Int_t ths = 0; ths < ThSteps; ths++) {
      theta = ths * (theta_max - theta_min) / (ThSteps - 1) + theta_min;
      for (Int_t phs = 0; phs < PhiSteps; phs++) {
        phi = phs * (phi_max - phi_min) / (PhiSteps - 1) + phi_min;

        if (PrintLevel > 0) {
          Log << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
              << "!!!       EVENT " << evt
              << "\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"
              << std::endl;
        }

        TraceCan->cd(1);
        HCT->Draw();

        for (Int_t stp = 0; stp < AnodeRows; stp++)
          for (Int_t col = 0; col < AnodeCols + 1; col++) {
            DeltaEB[stp][col] = DeltaEL[stp][col] = DeltaEH[stp][col] = 0;
            for (Int_t er = 0; er < numEvaporations; er++) {
              DeltaE_EvaP[er][stp][col] = 0;
              DeltaE_EvaR[er][stp][col] = 0;
            }
          }

        SetInitialKinematics(Kb_at_gas);

        this->zr = Rdm->Uniform(MinZ, MaxZ);
        this->Kbr = Beam->GetFinalEnergy(0, Kb_at_gas, this->zr);
        Double_t TOF = Beam->GetTimeOfFlight(0);
        if (PrintLevel > 0)
          Log << "Kbr = " << this->Kbr << "  zr = " << this->zr
              << "  tof = " << TOF << std::endl;

        Int_t ReacAllowed =
            SetReactionKinematics(this->Kbr, this->zr, TOF, theta, phi);
        if (PrintLevel > 0) {
          Log << "Conservation of 4-momentum at reaction point (zr)"
              << std::endl;
          FourVector Pi("initial 4-momentum (lab)", 0, 0, 0, 0);
          Pi += Beam->GetP() + Target->GetP();
          FourVector Pf("final 4-momentum (lab)", 0, 0, 0, 0);
          for (Int_t er = 0; er < numEvaporations; er++) {
            if (!EvaR[er]->DoNotPropagate)
              Pf += EvaR[er]->GetP();
            if (!EvaP[er]->DoNotPropagate)
              Pf += EvaP[er]->GetP();
          }
          Pi.Print(Log);
          Pf.Print(Log);
        }

        if (ReacAllowed) {
          Beam->GetX(tf, xf, yf, zf);
          PropagateParticle(Beam, evt, MaxTime, -UserStep, DeltaEB);
          Beam->GetX(ti, xi, yi,
                     zi); // post-propagate position is the new origin
          TrackBeam->SetOrigin(xi, yi, zi);
          TrackBeam->SetVector(xf - xi, yf - yi, zf - zi);

          for (Int_t er = 0; er < numEvaporations; er++) {
            EvaP[er]->GetX(ti, xi, yi, zi);
            PropagateParticle(EvaP[er], evt, MaxTime, UserStep,
                              DeltaE_EvaP[er]);
            EvaP[er]->GetX(tf, xf, yf, zf);
            TrackEvaP[er]->SetOrigin(xi, yi, zi);
            TrackEvaP[er]->SetVector(xf - xi, yf - yi, zf - zi);

            EvaR[er]->GetX(ti, xi, yi, zi);
            PropagateParticle(EvaR[er], evt, MaxTime, UserStep,
                              DeltaE_EvaR[er]);
            EvaR[er]->GetX(tf, xf, yf, zf);
            if (!EvaR[er]->DoNotPropagate) {
              xfe = xf;
              yfe = yf;
              zfe = zf;
              resID = er;
            }
            TrackEvaR[er]->SetOrigin(xi, yi, zi);
            TrackEvaR[er]->SetVector(xf - xi, yf - yi, zf - zi);
          }
        } else {
          Beam->Copy(BeamInit);
          PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB);
          std::cout << "Warning: reaction energetically not allowed for event "
                    << evt << " (Kbr= " << this->Kbr << " MeV)." << std::endl;
        }

        ComputeDetectorResponse(evt, stp_base, UpdateVis);
        if (SimTree != 0) {
          FinalizeEvent(evt);
          SimTree->Fill();
          if (MCTree)
            MCTree->Fill();
        }
        if (UpdateVis)
          UpdateVisuals(evt, this->Kbr, this->zr, TOF, Wait);

        EvtsProcessed++;
        if (NEvents > 99 &&
            (LongDouble_t)(EvtsProcessed) >= Frac[FIndex] * NEvents) {
          if (verbose_)
            std::cout << "\t" << Frac[FIndex] * 100 << "% processed ("
                      << StpWatch.RealTime() << " s)" << std::endl;
          StpWatch.Start(kFALSE);
          FIndex++;
        }
        NTraces++;
        evt++;
      }
    }
  }
  if (verbose_)
    std::cout << "Saving traces ..." << std::endl;
  if (verbose_)
    StpWatch.Print();

  TDB->cd();
  SimTree->Write("", TObject::kSingleKey);
  if (MCTree)
    MCTree->Write("", TObject::kSingleKey);
  TDB->Close();
  StpWatch.Stop();
  if (verbose_)
    StpWatch.Print();
}
