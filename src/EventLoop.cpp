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
void Simulator::Simulate(Int_t StpID, Int_t eventCount, Double_t MaxTime,
                         Double_t UserStep, Int_t UpdateVis, Int_t Wait,
                         TFile * /*ROOTfile*/) {
  if (VolAnode == 0) {
    std::cout << "Anode geometry not specified. Use SetAnode method."
              << std::endl;
    return;
  }

  this->NEvents = eventCount;

  TStopwatch StpWatch;
  LongDouble_t Frac[11] = {0.01, 0.05, 0.1, 0.2, 0.3, 0.5,
                           0.6,  0.7,  0.8, 0.9, 1.0};
  Int_t FIndex = 0;

  if (verbose_)
    std::cout << "Simulating " << eventCount << " MUSIC traces for strip "
              << StpID << " ... " << std::endl;

  // Idealized beam (no straggling) used for the kinematic-range estimate and
  // the unreacted-beam reference trace.
  SetInitialKinematics(Kb_at_gas);

  Particle beamCopy("beam copy");
  beamCopy.Copy(Beam);
  if (PrintLevel > 0)
    beamCopy.Print();
  PropagateParticle(&beamCopy, 0, MaxTime, UserStep, DeltaEB_ave);
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
  for (Int_t evt = 0; evt < eventCount; evt++) {
    SeedRandom(EventSeed(StpID, eventOffset_ + static_cast<ULong64_t>(evt)));
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
    eventReacted_ = kFALSE;
    beamTermination_ = TerminationReason::NotPropagated;
    evapTermination_.fill(TerminationReason::NotPropagated);
    residueTermination_.fill(TerminationReason::NotPropagated);

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
        for (Int_t er = 0; er < numEvaporations; er++) {
          DeltaE_EvaP[er][stp][col] = 0;
          DeltaE_EvaR[er][stp][col] = 0;
        }
      }

    // Per-event beam chain (physically correct order):
    //   accelerator KE [+ KbFWHM]  ->  degrader  ->  entrance window  ->  Kbi.
    Double_t Ebeam = ctf.BeamEnergy;
    if (ctf.KbFWHM > 0.0) {
      do {
        Ebeam = ctf.BeamEnergy + Rdm->Gaus(0.0, ctf.KbFWHM / 2.355);
      } while (Ebeam <= 0.0);
    }
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
    beam_energy_gas = static_cast<Float_t>(Ebeam);
    SetInitialKinematics(Ebeam);

    Int_t ReacAllowed = 0;
    beam_energy_reaction = -2.0f;
    Double_t attemptedReactionEnergy = -2.0;
    Double_t TOF = 0;
    if (StpID > -1) {
      // Pick the vertex uniformly inside the chosen strip; forward-propagate
      // the beam to the vertex with per-step straggling so the reaction
      // energy / TOF inherit the sampled fluctuations (the reaction
      // kinematics see the actually-sampled energy at the vertex, not the
      // deterministic mean).
      const Double_t requestedVertex = Rdm->Uniform(MinZ, MaxZ);
      Beam->GetX(ti, xi, yi, zi); // entrance origin (kept for visualization)
      const PropagationResult toVertex =
          PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB,
                            /*endZ=*/requestedVertex);
      beamTermination_ = toVertex.reason;
      Double_t t_at_vertex, x_at_vertex, y_at_vertex, z_at_vertex;
      Beam->GetX(t_at_vertex, x_at_vertex, y_at_vertex, z_at_vertex);
      TOF = t_at_vertex;
      if (tracesCreated)
        for (Int_t stp = 0; stp < AnodeRows; stp++)
          for (Int_t col = 0; col < AnodeCols + 1; col++)
            TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      if (UpdateVis) {
        TrackBeam->SetOrigin(static_cast<Float_t>(xi), static_cast<Float_t>(yi),
                             static_cast<Float_t>(zi));
        TrackBeam->SetVector(static_cast<Float_t>(x_at_vertex - xi),
                             static_cast<Float_t>(y_at_vertex - yi),
                             static_cast<Float_t>(z_at_vertex - zi));
      }
      if (toVertex.reachedVertex()) {
        attemptedReactionEnergy = Beam->GetKE();
        ReacAllowed = SetReactionKinematics(attemptedReactionEnergy,
                                            requestedVertex, TOF);
        if (ReacAllowed) {
          eventReacted_ = kTRUE;
          beam_energy_reaction = static_cast<Float_t>(attemptedReactionEnergy);
          vertex_x = static_cast<Float_t>(x_at_vertex);
          vertex_y = static_cast<Float_t>(y_at_vertex);
          vertex_z = static_cast<Float_t>(requestedVertex);
        }
      } else {
        if (PrintLevel > 0)
          Log << "Beam did not reach requested reaction vertex "
              << requestedVertex << " cm; termination reason "
              << static_cast<Int_t>(toVertex.reason) << std::endl;
      }
      if (PrintLevel > 0 && toVertex.reachedVertex())
        Log << "Kbr = " << attemptedReactionEnergy
            << "  zr = " << requestedVertex << "  tof = " << TOF << std::endl;

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
        evapTermination_[er] =
            PropagateParticle(EvaP[er], evt, MaxTime, UserStep, DeltaE_EvaP[er])
                .reason;
        if (tracesCreated)
          for (Int_t stp = 0; stp < AnodeRows; stp++)
            for (Int_t col = 0; col < AnodeCols + 1; col++)
              TraceEP[er][col]->SetPoint(stp, stp, DeltaE_EvaP[er][stp][col]);
        EvaP[er]->GetX(tf, xf, yf, zf);
        if (UpdateVis) {
          TrackEvaP[er]->SetOrigin(static_cast<Float_t>(xi),
                                   static_cast<Float_t>(yi),
                                   static_cast<Float_t>(zi));
          TrackEvaP[er]->SetVector(static_cast<Float_t>(xf - xi),
                                   static_cast<Float_t>(yf - yi),
                                   static_cast<Float_t>(zf - zi));
        }

        EvaR[er]->GetX(ti, xi, yi, zi);
        residueTermination_[er] =
            PropagateParticle(EvaR[er], evt, MaxTime, UserStep, DeltaE_EvaR[er])
                .reason;
        if (tracesCreated)
          for (Int_t stp = 0; stp < AnodeRows; stp++)
            for (Int_t col = 0; col < AnodeCols + 1; col++)
              TraceER[er][col]->SetPoint(stp, stp, DeltaE_EvaR[er][stp][col]);
        EvaR[er]->GetX(tf, xf, yf, zf);
        if (UpdateVis) {
          TrackEvaR[er]->SetOrigin(static_cast<Float_t>(xi),
                                   static_cast<Float_t>(yi),
                                   static_cast<Float_t>(zi));
          TrackEvaR[er]->SetVector(static_cast<Float_t>(xf - xi),
                                   static_cast<Float_t>(yf - yi),
                                   static_cast<Float_t>(zf - zi));
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
      const Bool_t resume =
          (StpID > -1 && beamTermination_ == TerminationReason::ReachedVertex);
      Beam->GetX(ti, xi, yi, zi);
      if (StpID == -1 || resume)
        beamTermination_ =
            PropagateParticle(Beam, evt, MaxTime, UserStep, DeltaEB,
                              /*endZ=*/-1.0, /*reset_DE=*/!resume)
                .reason;
      Beam->GetX(tf, xf, yf, zf);
      if (UpdateVis) {
        TrackBeam->SetOrigin(static_cast<Float_t>(xi), static_cast<Float_t>(yi),
                             static_cast<Float_t>(zi));
        TrackBeam->SetVector(static_cast<Float_t>(xf - xi),
                             static_cast<Float_t>(yf - yi),
                             static_cast<Float_t>(zf - zi));
      }
      if (tracesCreated)
        for (Int_t stp = 0; stp < AnodeRows; stp++)
          for (Int_t col = 0; col < AnodeCols + 1; col++)
            TraceB[col]->SetPoint(stp, stp, DeltaEB[stp][col]);
      if (StpID > -1 && resume) {
        std::cout << "Warning: reaction energetically not allowed for event "
                  << evt << " (Kbr= " << attemptedReactionEnergy << " MeV)."
                  << std::endl;
        Log << "Warning: reaction energetically not allowed for event " << evt
            << " (Kbr= " << attemptedReactionEnergy << " MeV)." << std::endl;
      }
    }

    ComputeDetectorResponse(evt, StpID, UpdateVis);
    Log << "\tDone computing detector response." << std::endl;
    if (SimTree != 0) {
      FinalizeEvent(evt);
      if (SimTree->Fill() < 0)
        ioFailed_ = kTRUE;
      if (MCTree && MCTree->Fill() < 0)
        ioFailed_ = kTRUE;
    }
    if (UpdateVis)
      UpdateVisuals(evt, beam_energy_reaction, vertex_z, TOF, Wait);

    if (eventCount > 99 &&
        static_cast<LongDouble_t>(evt) >= Frac[FIndex] * eventCount) {
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
