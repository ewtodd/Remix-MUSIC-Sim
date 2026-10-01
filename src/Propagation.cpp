#include "Simulator.hpp"

// FWHM = kFwhmPerSigma * sigma for a Gaussian; converts the % FWHM noise
// figures in detector.resolution_fwhm_percent to the Gaussian sigma.
static const Double_t kFwhmPerSigma = 2.0 * std::sqrt(2.0 * std::log(2.0));

// Transport a particle through the gas while ending every step at the first
// cell, detector, time, or requested-vertex boundary it intersects. Every
// traversed fraction is charged to the cell it actually crosses.
Simulator::PropagationResult
Simulator::PropagateParticle(Particle *PO, Int_t Event, Double_t MaxTime,
                             Double_t UserStep, Double_t **DE, Double_t endZ,
                             Bool_t reset_DE) {
  PropagationResult result;
  if (reset_DE) {
    for (Int_t stp = 0; stp < AnodeRows; stp++)
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        DE[stp][col] = 0;
  }
  auto rebuildTotals = [&]() {
    for (Int_t row = 0; row < AnodeRows; ++row) {
      Double_t total = 0.0;
      for (Int_t col = 0; col < AnodeCols; ++col)
        total += DE[row][col];
      DE[row][AnodeCols] = total;
    }
  };
  if (!PO || PO->DoNotPropagate) {
    result.reason = TerminationReason::NotPropagated;
    rebuildTotals();
    return result;
  }
  if (!std::isfinite(UserStep) || UserStep <= 0.0 || !std::isfinite(MaxTime) ||
      MaxTime <= 0.0) {
    result.reason = TerminationReason::InvalidState;
    rebuildTotals();
    return result;
  }

  Double_t time, x, y, z;
  PO->GetX(time, x, y, z);
  Double_t energy, px, py, pz;
  PO->GetP(energy, px, py, pz);
  Double_t kinetic = PO->GetKE();
  const Double_t restEnergy = std::max(0.0, energy - kinetic);
  const Double_t momentum0 = std::sqrt(px * px + py * py + pz * pz);
  if (!std::isfinite(momentum0) || momentum0 <= 0.0 || kinetic <= 0.0) {
    PO->SetP(restEnergy, 0.0, 0.0, 0.0);
    result.reason = TerminationReason::Stopped;
    rebuildTotals();
    return result;
  }
  const Double_t ux = px / momentum0;
  const Double_t uy = py / momentum0;
  const Double_t uz = pz / momentum0;

  PO->ResetTrace();
  if (PO->SaveTrajectory && PO->Trajectory) {
    PO->Trajectory->RemoveElements();
    PO->Trajectory->SetName(Form("%s evt %d", PO->Name.Data(), Event));
    PO->SetTracePoint(static_cast<Float_t>(time), static_cast<Float_t>(x),
                      static_cast<Float_t>(y), static_cast<Float_t>(z),
                      static_cast<Float_t>(kinetic));
  }

  constexpr Double_t eps = 1e-10;
  auto classifyOutside = [&]() {
    if (z >= AnodeDepth - eps)
      return TerminationReason::DownstreamExit;
    if (z <= eps)
      return TerminationReason::UpstreamExit;
    return TerminationReason::SideExit;
  };
  auto setFinalState = [&]() {
    const Double_t momentum =
        kinetic > 0.0
            ? std::sqrt(std::max(0.0, kinetic * (kinetic + 2.0 * restEnergy)))
            : 0.0;
    PO->SetX(time, x, y, z);
    PO->SetP(restEnergy + kinetic, momentum * ux, momentum * uy, momentum * uz);
    rebuildTotals();
  };

  if (endZ >= 0.0 && std::fabs(z - endZ) <= eps) {
    result.reason = TerminationReason::ReachedVertex;
    setFinalState();
    return result;
  }
  if (z < -eps || z > AnodeDepth + eps || x < -AnodeLength / 2.0 - eps ||
      x > AnodeLength / 2.0 + eps || y < -AnodeHeight / 2.0 - eps ||
      y > AnodeHeight / 2.0 + eps) {
    result.reason = classifyOutside();
    setFinalState();
    return result;
  }

  Int_t step = 0;
  Double_t distanceSinceTrace = 0.0;
  Double_t traceX = x, traceY = y, traceZ = z;
  while (step++ < 10000000) {
    if (time >= MaxTime - eps) {
      result.reason = TerminationReason::Timeout;
      break;
    }

    const Double_t momentum =
        std::sqrt(std::max(0.0, kinetic * (kinetic + 2.0 * restEnergy)));
    const Double_t totalEnergy = restEnergy + kinetic;
    const Double_t velocity =
        totalEnergy > 0.0 ? c * momentum / totalEnergy : 0.0;
    if (!std::isfinite(velocity) || velocity <= 0.0) {
      result.reason = TerminationReason::Stopped;
      kinetic = 0.0;
      break;
    }

    // Probe just ahead so a point exactly on an internal boundary is assigned
    // to the cell it is entering, not the one it just left.
    const Double_t probeX = x + ux * eps;
    const Double_t probeZ = z + uz * eps;
    Int_t row = -1;
    Int_t col = -1;
    Double_t rowLow = 0.0;
    Double_t rowHigh = 0.0;
    for (Int_t candidate = 0; candidate < AnodeRows; ++candidate) {
      rowHigh = rowLow + AnodeDZ[candidate][0];
      if (probeZ >= rowLow && probeZ < rowHigh) {
        row = candidate;
        break;
      }
      rowLow = rowHigh;
    }
    Double_t colLow = -AnodeLength / 2.0;
    Double_t colHigh = colLow;
    if (row >= 0) {
      for (Int_t candidate = 0; candidate < AnodeCols; ++candidate) {
        colHigh = colLow + AnodeDX[row][candidate];
        if (AnodeDX[row][candidate] > 0.0 && probeX >= colLow &&
            probeX < colHigh) {
          col = candidate;
          break;
        }
        colLow = colHigh;
      }
    }
    if (row < 0 || col < 0) {
      result.reason = classifyOutside();
      break;
    }

    Double_t distance = UserStep;
    TerminationReason boundaryReason = TerminationReason::InvalidState;
    Bool_t endsAtBoundary = kFALSE;
    auto consider = [&](Double_t candidate, TerminationReason reason,
                        Bool_t terminal) {
      if (candidate >= -eps && candidate < distance - eps) {
        distance = std::max(0.0, candidate);
        boundaryReason = reason;
        endsAtBoundary = terminal;
      }
    };

    consider(velocity * (MaxTime - time), TerminationReason::Timeout, kTRUE);
    if (endZ >= 0.0 && std::fabs(uz) > eps && (endZ - z) * uz >= -eps)
      consider((endZ - z) / uz, TerminationReason::ReachedVertex, kTRUE);
    if (uz > eps) {
      consider((AnodeDepth - z) / uz, TerminationReason::DownstreamExit, kTRUE);
      consider((rowHigh - z) / uz, TerminationReason::InvalidState, kFALSE);
    } else if (uz < -eps) {
      consider((0.0 - z) / uz, TerminationReason::UpstreamExit, kTRUE);
      consider((rowLow - z) / uz, TerminationReason::InvalidState, kFALSE);
    }
    if (ux > eps) {
      consider((AnodeLength / 2.0 - x) / ux, TerminationReason::SideExit,
               kTRUE);
      consider((colHigh - x) / ux, TerminationReason::InvalidState, kFALSE);
    } else if (ux < -eps) {
      consider((-AnodeLength / 2.0 - x) / ux, TerminationReason::SideExit,
               kTRUE);
      consider((colLow - x) / ux, TerminationReason::InvalidState, kFALSE);
    }
    if (uy > eps)
      consider((AnodeHeight / 2.0 - y) / uy, TerminationReason::SideExit,
               kTRUE);
    else if (uy < -eps)
      consider((-AnodeHeight / 2.0 - y) / uy, TerminationReason::SideExit,
               kTRUE);

    if (!std::isfinite(distance) || distance <= eps) {
      if (endsAtBoundary) {
        result.reason = boundaryReason;
        break;
      }
      result.reason = TerminationReason::InvalidState;
      break;
    }

    music::EnergyStepResult energyStep;
    if (gasEnabled_)
      energyStep = PO->TransportStep(0, kinetic, distance, Rdm);
    else {
      energyStep.finalEnergy = kinetic;
      energyStep.distance = distance;
    }
    const Double_t travelled = std::clamp(energyStep.distance, 0.0, distance);
    const Double_t dt = travelled / velocity;
    time += dt;
    x += ux * travelled;
    y += uy * travelled;
    z += uz * travelled;
    kinetic = std::max(0.0, energyStep.finalEnergy);
    DE[row][col] += energyStep.depositedEnergy;
    result.distance += travelled;
    result.depositedEnergy += energyStep.depositedEnergy;

    distanceSinceTrace += travelled;
    if (PO->SaveTrajectory && PO->Trajectory &&
        distanceSinceTrace >= 10.0 * UserStep) {
      PO->Trajectory->AddLine(
          static_cast<Float_t>(traceX), static_cast<Float_t>(traceY),
          static_cast<Float_t>(traceZ), static_cast<Float_t>(x),
          static_cast<Float_t>(y), static_cast<Float_t>(z));
      PO->SetTracePoint(static_cast<Float_t>(time), static_cast<Float_t>(x),
                        static_cast<Float_t>(y), static_cast<Float_t>(z),
                        static_cast<Float_t>(kinetic));
      traceX = x;
      traceY = y;
      traceZ = z;
      distanceSinceTrace = 0.0;
    }

    if (energyStep.stopped || kinetic <= 0.0) {
      result.reason = TerminationReason::Stopped;
      kinetic = 0.0;
      break;
    }
    if (endsAtBoundary && travelled >= distance - eps) {
      result.reason = boundaryReason;
      break;
    }
  }
  if (step >= 10000000)
    result.reason = TerminationReason::InvalidState;
  setFinalState();
  if (PrintLevel > 0)
    Log << "musicsim::PropagateParticle END (reason "
        << static_cast<Int_t>(result.reason) << ")" << std::endl;
  return result;
}

void Simulator::ComputeDetectorResponse(Int_t evt, Int_t reacStp,
                                        Int_t /*UpdateVis*/) {
  if (PrintLevel > 0)
    Log << "Compute detector response evt " << evt << std::endl;

  if (SimTree != 0) {
    requested_strip = reacStp;
    reaction_strip = eventReacted_ ? reacStp : -1;
  }

  auto addDeposit = [](Float_t &destination, Double_t deposit) {
    destination =
        static_cast<Float_t>(static_cast<Double_t>(destination) + deposit);
  };

  for (Int_t row = 0; row < AnodeRows; row++) {
    // Strip ID for this row. AnodeStpID[row][1] is unset (-1) for the
    // single-column rows (S0, S17, dead layers); col 0 carries the real id.
    Int_t rowStpid = AnodeStpID[row][0];

    // First pass: assemble noiseless dE per column from all particles, and
    // draw per-electrode noise for real anode electrodes (col 0..AnodeCols-1
    // on readout-strip rows). Dead layers get no noise — we have no
    // experimental info to anchor a sigma. The col == AnodeCols sum column
    // is a derived view used only for trace visualization, not a real
    // electrode, so it carries no independent noise either.
    std::vector<Double_t> baseDE(static_cast<size_t>(AnodeCols + 1), 0.0);
    std::vector<Double_t> noisedDE(static_cast<size_t>(AnodeCols + 1), 0.0);
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      Double_t DeltaE = DeltaEB[row][col];
      for (Int_t er = 0; er < numEvaporations; er++) {
        DeltaE += DeltaE_EvaP[er][row][col];
        DeltaE += DeltaE_EvaR[er][row][col];
      }
      baseDE[col] = DeltaE;
      noisedDE[col] = DeltaE;
      if (col < AnodeCols && rowStpid >= 0 && rowStpid <= 17) {
        // Relative resolution is % FWHM of the deposit, so its Gaussian sigma
        // scales with this event's energy on the electrode. Absolute noise is
        // already configured as sigma in MeV.
        const Double_t configured = ctf.Noise[ElectrodeIndex(rowStpid, col)];
        if (configured > 0.0 && DeltaE > 0.0) {
          const Double_t sigma =
              ctf.noiseMode == 1 ? configured
                                 : configured / 100.0 * DeltaE / kFwhmPerSigma;
          noisedDE[col] += Rdm->Gaus(0.0, sigma);
        }
      }
    }

    // Output-tree accumulation. Dead-layer rows go to their own scalars
    // (noiseless). Readout-strip rows split into per-electrode L/R, plus a
    // physically-summed Cathode contribution (one big plate).
    if (SimTree != 0) {
      if (rowStpid == -1) {
        addDeposit(DeadUS_dE, baseDE[0]);
      } else if (rowStpid == -2) {
        addDeposit(DeadDS_dE, baseDE[0]);
      } else if (rowStpid >= 0 && rowStpid <= 17) {
        for (Int_t col = 0; col < AnodeCols; col++) {
          // Skip the unused col=1 slot on full-width rows (S0, S17).
          if (AnodeStpID[row][col] != rowStpid)
            continue;
          // Cathode sees ions regardless of which anode finger sits above;
          // it gets the noiseless dE summed across electrodes, and a single
          // independent Gaussian is added after the row loop.
          addDeposit(Cathode, baseDE[col]);
          if (rowStpid == 0) {
            // Unsegmented strips 0 and 17: one scalar each, as in the
            // experimental events tree.
            addDeposit(Strip0dE, noisedDE[col]);
          } else if (rowStpid == 17) {
            addDeposit(Strip17dE, noisedDE[col]);
          } else {
            // Segmented strips 1..16 live in arrays of 16 indexed by strip - 1.
            if (col == 0)
              addDeposit(RightdE[rowStpid - 1], noisedDE[col]);
            else
              addDeposit(LeftdE[rowStpid - 1], noisedDE[col]);
          }
        }
      }
    }

    // Trace TGraphs (interactive visualization only): one point per readout
    // strip, indexed by stpid on the x-axis. Dead-layer rows are skipped.
    if (tracesCreated && rowStpid >= 0 && rowStpid <= 17) {
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        Trace[col]->SetPoint(rowStpid, rowStpid, noisedDE[col]);
    }
  }

  // Single independent Gaussian for the cathode readout channel (one plate,
  // one electronics chain). Applied after the row loop so the per-electrode
  // anode noise does not leak in. The cathode uses the same explicitly
  // selected absolute-sigma or relative-FWHM convention as the anodes.
  if (SimTree != 0 && ctf.NoiseCathode > 0.0 &&
      static_cast<Double_t>(Cathode) > 0.0) {
    const Double_t sigma =
        ctf.noiseMode == 1 ? ctf.NoiseCathode
                           : ctf.NoiseCathode / 100.0 *
                                 static_cast<Double_t>(Cathode) / kFwhmPerSigma;
    addDeposit(Cathode, Rdm->Gaus(0.0, sigma));
  }

  if (tracesCreated) {
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      if (col == AnodeCols)
        Trace[col]->Write(Form("Trace_s%d_e%d", reacStp, evt),
                          TObject::kOverwrite);
      else
        Trace[col]->Write(Form("Trace_s%d_c%d_e%d", reacStp, col, evt),
                          TObject::kOverwrite);
    }
  }
}

// Fill the exit-energy and stop-location branches from each particle's final
// state. Energy sentinels (set in ResetBranches): -1.0 = transported but did
// not reach the downstream exit, -2.0 = N/A. Stop-strip sentinels: -1 = did
// not stop in a readout strip (exited, or dead layer), -2 = N/A.
void Simulator::ComputeExitEnergies() {
  const Double_t amu_MeV = 931.49410242;
  auto Aof = [&](Particle *P) -> Int_t {
    return (P && P->Mass > 0) ? Int_t(std::round(P->Mass / amu_MeV)) : 0;
  };
  auto throughExitWindow = [&](Particle *particle, Double_t Kf) -> Float_t {
    if (!exitWindowEnabled_)
      return static_cast<Float_t>(Kf);
    const Double_t cosine = std::fabs(std::cos(particle->GetTheta()));
    if (cosine <= 1e-9)
      return 0.0f;
    catima::Material obliqueWindow = exitWindow_;
    obliqueWindow.thickness(exitWindow_.thickness() / cosine);
    return static_cast<Float_t>(
        EnergyOutOfMaterial(Aof(particle), particle->Z, Kf, obliqueWindow));
  };
  // Record the particle's final position and stop strip; return the exit
  // energy (through the exit window) or the non-downstream sentinel. Only call
  // for particles that were actually transported.
  auto recordExit = [&](Particle *P, TerminationReason reason, Float_t &sx,
                        Float_t &sy, Float_t &sz, Int_t &sstrip) -> Float_t {
    Double_t t, x, y, z;
    P->GetX(t, x, y, z);
    sx = static_cast<Float_t>(x);
    sy = static_cast<Float_t>(y);
    sz = static_cast<Float_t>(z);
    if (reason == TerminationReason::DownstreamExit) {
      sstrip = -1; // left out the back; sz holds the crossing point
      return throughExitWindow(P, P->GetKE());
    }
    sstrip = reason == TerminationReason::Stopped ? StripAtZ(z) : -1;
    return -1.0f;
  };
  auto transported = [](Particle *particle,
                        TerminationReason reason) -> Bool_t {
    return particle && !particle->DoNotPropagate &&
           reason != TerminationReason::NotPropagated;
  };

  // The chain's surviving residue: each allowed step stops the propagation of
  // the previous step's residue, so exactly the last one is active after a
  // successful reaction. Select it independently of the transport outcome:
  // an InvalidState termination must not make a real reaction look unreacted.
  if (eventReacted_)
    for (Int_t er = 0; er < numEvaporations; ++er)
      if (EvaR[er] && !EvaR[er]->DoNotPropagate)
        residue_step = er;

  // The beam survives only on unreacted events (no vertex picked, or the
  // sampled reaction was energetically disallowed and the beam swept on to
  // the exit). On reacted events it was consumed at the vertex: N/A.
  if (Beam && !eventReacted_)
    beam_energy_exit = recordExit(Beam, beamTermination_, beam_stop_x,
                                  beam_stop_y, beam_stop_z, beam_stop_strip);

  for (Int_t er = 0; er < numEvaporations; ++er)
    if (transported(EvaP[er], evapTermination_[er]))
      evap_energy_exit[er] =
          recordExit(EvaP[er], evapTermination_[er], evap_stop_x[er],
                     evap_stop_y[er], evap_stop_z[er], evap_stop_strip[er]);
  // Superseded residues decayed in place — their exit slots stay -2 (N/A).
  if (residue_step >= 0 &&
      transported(EvaR[residue_step], residueTermination_[residue_step]))
    residue_energy_exit[residue_step] = recordExit(
        EvaR[residue_step], residueTermination_[residue_step], residue_stop_x,
        residue_stop_y, residue_stop_z, residue_stop_strip);
}

void Simulator::FinalizeEvent(Int_t eventIndexInSlice) {
  event_index = eventOffset_ + static_cast<ULong64_t>(eventIndexInSlice);
  beam_termination = static_cast<Int_t>(beamTermination_);
  for (Int_t step = 0; step < numEvaporations; ++step)
    evap_termination[step] = static_cast<Int_t>(evapTermination_[step]);
  residue_termination =
      residue_step >= 0 ? static_cast<Int_t>(residueTermination_[residue_step])
                        : static_cast<Int_t>(TerminationReason::NotPropagated);
  ComputeExitEnergies();
  // ComputeExitEnergies resolves residue_step, so update the scalar after it.
  residue_termination =
      residue_step >= 0 ? static_cast<Int_t>(residueTermination_[residue_step])
                        : static_cast<Int_t>(TerminationReason::NotPropagated);
}
