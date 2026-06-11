#include "Simulator.hpp"

// FWHM = kFwhmPerSigma * sigma for a Gaussian; converts the % FWHM noise
// figures in [detector.eres] to the sigma fed to the random draw.
static const Double_t kFwhmPerSigma = 2.0 * std::sqrt(2.0 * std::log(2.0));

// Transport a particle step-by-step through the gas, accumulating per-strip
// energy deposits in DE. Stops when MaxTime is hit, the particle stops, the
// particle leaves the active volume, or it crosses endZ (when > 0). The last
// column of DE (index AnodeCols) holds the per-strip total.
Int_t Simulator::PropagateParticle(Particle *PO, Int_t Event, Double_t MaxTime,
                                   Double_t UserStep, Double_t **DE,
                                   Double_t endZ, Bool_t reset_DE) {
  if (reset_DE) {
    for (Int_t stp = 0; stp < AnodeRows; stp++)
      for (Int_t col = 0; col < AnodeCols + 1; col++)
        DE[stp][col] = 0;
  }
  const Double_t active_max_z =
      (endZ > 0.0) ? std::min(endZ, AnodeDepth) : AnodeDepth;

  if (PrintLevel > 0) {
    Log << "\nmusicsim::PropagateParticle START *********************************\n"
        << PO->Name << std::endl;
  }
  if (PO->DoNotPropagate) {
    if (PrintLevel > 0)
      Log << "Not propagating!" << std::endl;
    return 0;
  }

  Int_t step = 0;
  Double_t m = PO->Mass;
  Double_t tf = 0, xf = 0, yf = 0, zf = 0;
  Double_t ti, xi, yi, zi;
  PO->GetX(ti, xi, yi, zi);
  Double_t Ene, px, py, pz;
  PO->GetP(Ene, px, py, pz);
  Double_t p_mag = std::sqrt(px * px + py * py + pz * pz);
  Double_t Ki = PO->GetKE();

  PO->ResetTrace();
  Double_t tt = ti, xt = xi, yt = yi, zt = zi, Kt = Ki;
  PO->SetTracePoint((Float_t)tt, (Float_t)xt, (Float_t)yt, (Float_t)zt,
                    (Float_t)Kt);
  // Trace points are emitted every kTraceStepFactor simulation steps so we
  // stay well under Particle::MaxPoints even for full traversals.
  const Int_t kTraceStepFactor = 10;
  Double_t dist_since_trace = 0;
  PO->Trajectory->RemoveElements();
  PO->Trajectory->SetName(Form("%s evt %d", PO->Name.Data(), Event));

  if (PrintLevel > 0)
    Log << "MaxTime=" << MaxTime << " ns\nInitial time ti=" << ti << " ns"
        << std::endl;

  Double_t phi = PO->GetPhi();
  Double_t theta = PO->GetTheta();

  while (ti < MaxTime) {
    // Velocity in cm/ns. v = c·p/E (with c=30 cm/ns in our units).
    px = p_mag * std::cos(phi) * std::sin(theta);
    py = p_mag * std::sin(phi) * std::sin(theta);
    pz = p_mag * std::cos(theta);
    Double_t vx = c * px / Ene;
    Double_t vy = c * py / Ene;
    Double_t vz = c * pz / Ene;
    Double_t vel = std::sqrt(vx * vx + vy * vy + vz * vz);

    Double_t Dt = 0;
    if (vel > 0)
      Dt = UserStep / vel;

    tf = ti + Dt;
    xf = xi + vx * Dt;
    yf = yi + vy * Dt;
    zf = zi + vz * Dt;

    if (PrintLevel > 0 && (PrintLevel > 1 || step == 0)) {
      Log << "step " << step << " (";
      Log.precision(7);
      Log << tf << "ns, " << xf << "cm, " << yf << "cm, " << zf
          << "cm)  Dt=" << Dt << "ns" << std::endl;
    }

    // Exit if the particle leaves the active volume (or crosses endZ).
    if (zf > active_max_z || zf < 0 || xf > AnodeLength / 2 ||
        xf < -AnodeLength / 2 || yf > AnodeHeight / 2 ||
        yf < -AnodeHeight / 2) {
      if (PrintLevel > 0)
        Log << "Particle reached end of active volume." << std::endl;
      break;
    }

    // Find the anode cell containing (xf, yf, zf) by walking the AnodeDX /
    // AnodeDZ tables directly. The original code called Geo->FindNode but
    // that requires a TGeoManager (not available to MT workers) and the
    // geometry is a regular rectangular grid anyway.
    Int_t stp = -1;
    Int_t col = -1;
    {
      Double_t zacc = 0;
      for (Int_t r = 0; r < AnodeRows; ++r) {
        Double_t dz = AnodeDZ[r][0];
        if (zf >= zacc && zf < zacc + dz) {
          stp = r;
          break;
        }
        zacc += dz;
      }
      if (stp >= 0) {
        Double_t xacc = -AnodeLength / 2.0;
        for (Int_t c = 0; c < AnodeCols; ++c) {
          Double_t dx = AnodeDX[stp][c];
          if (dx <= 0)
            continue;
          if (xf >= xacc && xf < xacc + dx) {
            col = c;
            break;
          }
          xacc += dx;
        }
      }
    }
    if (stp == -1 || col == -1) {
      if (PrintLevel > 0)
        Log << "WARNING: anode segment not found." << std::endl;
      break;
    }

    Double_t dist = std::sqrt(std::pow(xf - xi, 2) + std::pow(yf - yi, 2) +
                              std::pow(zf - zi, 2));
    Double_t Kf = 0;
    if (!gasEnabled_) {
      // Pressure <= 0: particle propagates with no energy loss.
      Kf = Ki;
    } else if (UserStep > 0) {
      // Forward physics step: sample per-step Vavilov straggling.
      Kf = PO->GetFinalEnergyStraggled(0, Ki, dist, Rdm);
    } else {
      // Backward kinematic propagation (beam reconstructed back to entrance).
      // Inverse straggling isn't well-defined; use the mean.
      Kf = PO->GetInitialEnergy(0, Ki, dist);
    }

    if (Kf < 0.001) {
      if (PrintLevel > 0) {
        Log << "Particle stops inside active volume (Kf<1 keV)." << std::endl;
        if (Kf < 0)
          Log << "WARNING: Less than ZERO K.E.! " << Kf << std::endl;
      }
      break;
    }

    DE[stp][col] += std::fabs(Ki - Kf);

    if (PrintLevel > 0 && (PrintLevel > 1 || step == 0))
      Log << "d=" << dist << " cm \tKf=" << Kf << " \tDE=" << DE[stp][col]
          << std::endl;

    p_mag = std::sqrt(2 * m * Kf * (1 + Kf / 2 / m));
    // Signed kinetic-energy update handles both directions:
    //   forward (UserStep > 0): Kf ≤ Ki (clamped), so Ene decreases;
    //   backward (UserStep ≤ 0): Kf > Ki (deterministic), so Ene increases.
    // The sign of (Ki - Kf) carries the direction so we don't need fabs() here.
    Ene -= (Ki - Kf);

    ti = tf;
    xi = xf;
    yi = yf;
    zi = zf;
    Ki = Kf;
    step++;

    // Original code compared (tf - tt) > UserStep, which was a units mismatch
    // (ns vs cm) and effectively never fired. Drop a trace point every
    // kTraceStepFactor steps instead.
    dist_since_trace += dist;
    if (dist_since_trace >= kTraceStepFactor * std::fabs(UserStep)) {
      if (PO->SaveTrajectory)
        PO->Trajectory->AddLine(xt, yt, zt, xf, yf, zf);
      tt = ti;
      xt = xi;
      yt = yi;
      zt = zi;
      Kt = Ki;
      PO->SetTracePoint((Float_t)tt, (Float_t)xt, (Float_t)yt, (Float_t)zt,
                        (Float_t)Kt);
      dist_since_trace = 0;
    }
  }

  for (Int_t stp = 0; stp < AnodeRows; stp++)
    for (Int_t col = 0; col < AnodeCols; col++)
      if (DE[stp][col] > 0)
        DE[stp][AnodeCols] += DE[stp][col];

  PO->SetX(tf, xf, yf, zf);
  PO->SetP(Ene, p_mag * std::cos(phi) * std::sin(theta),
           p_mag * std::sin(phi) * std::sin(theta), p_mag * std::cos(theta));
  Log << "musicsim::PropagateParticle END ***********************************"
      << std::endl;
  return 1;
}

void Simulator::ComputeDetectorResponse(Int_t evt, Int_t reacStp,
                                        Int_t /*UpdateVis*/) {
  if (PrintLevel > 0)
    Log << "Compute detector response evt " << evt << std::endl;

  if (SimTree != 0)
    reaction_strip = reacStp;

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
    Double_t baseDE[3] = {0, 0, 0};
    Double_t noisedDE[3] = {0, 0, 0};
    for (Int_t col = 0; col < AnodeCols + 1; col++) {
      Double_t DeltaE = DeltaEB[row][col];
      for (Int_t er = 0; er < numEvaporations; er++) {
        DeltaE += DeltaE_EvaP[er][row][col];
        DeltaE += DeltaE_EvaR[er][row][col];
      }
      baseDE[col] = DeltaE;
      noisedDE[col] = DeltaE;
      if (col < AnodeCols && rowStpid >= 0 && rowStpid <= 17) {
        // Eres is a relative resolution in % FWHM of the deposit, so the
        // Gaussian sigma scales with this event's energy on the electrode.
        Double_t fwhmPct = ctf.Eres[ElectrodeIndex(rowStpid, col)];
        if (fwhmPct > 0.0 && DeltaE > 0.0) {
          Double_t sigma = fwhmPct / 100.0 * DeltaE / kFwhmPerSigma;
          noisedDE[col] += Rdm->Gaus(0.0, sigma);
        }
      }
    }

    // Output-tree accumulation. Dead-layer rows go to their own scalars
    // (noiseless). Readout-strip rows split into per-electrode L/R, plus a
    // physically-summed Cathode contribution (one big plate).
    if (SimTree != 0) {
      if (rowStpid == -1) {
        DeadUS_dE += baseDE[0];
      } else if (rowStpid == -2) {
        DeadDS_dE += baseDE[0];
      } else if (rowStpid >= 0 && rowStpid <= 17) {
        for (Int_t col = 0; col < AnodeCols; col++) {
          // Skip the unused col=1 slot on full-width rows (S0, S17).
          if (AnodeStpID[row][col] != rowStpid)
            continue;
          // Cathode sees ions regardless of which anode finger sits above;
          // it gets the noiseless dE summed across electrodes, and a single
          // independent Gaussian is added after the row loop.
          Cathode += baseDE[col];
          if (rowStpid == 0 || rowStpid == 17) {
            // Single-ended guard strips: full energy in the left slot, RightdE
            // stays 0 (matches the experimental Left_0_17_dE convention).
            Left_0_17_dE[rowStpid] += noisedDE[col];
          } else {
            if (col == 0)
              RightdE[rowStpid] += noisedDE[col];
            else
              Left_0_17_dE[rowStpid] += noisedDE[col];
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
  // anode noise does not leak in. EresCathode is % FWHM of the summed
  // cathode energy, same convention as the anode Eres entries.
  if (SimTree != 0 && ctf.EresCathode > 0.0 && Cathode > 0.0)
    Cathode += Rdm->Gaus(0.0, ctf.EresCathode / 100.0 * Cathode / kFwhmPerSigma);

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
// state. Energy sentinels (set in ResetBranches): -1.0 = stopped in the gas,
// -2.0 = N/A. Stop-strip sentinels: -1 = did not stop in a readout strip
// (exited, or dead layer), -2 = N/A.
void Simulator::ComputeExitEnergies() {
  const Double_t amu_MeV = 931.49410242;
  auto Aof = [&](Particle *P) -> Int_t {
    return (P && P->Mass > 0) ? Int_t(std::round(P->Mass / amu_MeV)) : 0;
  };
  auto throughExitWindow = [&](Int_t A, Int_t Z, Double_t Kf) -> Float_t {
    return exitWindowEnabled_
               ? (Float_t)EnergyOutOfMaterial(A, Z, Kf, exitWindow_)
               : (Float_t)Kf;
  };
  // Record the particle's final position and stop strip; return the exit
  // energy (through the exit window) or the stopped-in-gas sentinel. Only
  // call for particles that were actually transported.
  auto recordExit = [&](Particle *P, Float_t &sx, Float_t &sy, Float_t &sz,
                        Int_t &sstrip) -> Float_t {
    Double_t t, x, y, z;
    P->GetX(t, x, y, z);
    sx = (Float_t)x;
    sy = (Float_t)y;
    sz = (Float_t)z;
    if (z >= AnodeDepth) {
      sstrip = -1; // left out the back; sz holds the crossing point
      return throughExitWindow(Aof(P), P->Z, P->GetKE());
    }
    sstrip = StripAtZ(z);
    return -1.0f; // stopped inside the gas
  };
  // A particle was transported iff it isn't flagged DoNotPropagate and is not
  // at rest. Disallowed-step products are parked at rest (KE = 0); a particle
  // that ran out of energy mid-gas keeps the KE of its last completed step
  // (> 1 keV), so this never misclassifies a stopped particle.
  auto transported = [](Particle *P) -> Bool_t {
    return P && !P->DoNotPropagate && P->GetKE() > 0.001;
  };

  // The chain's surviving residue: each allowed step stops the propagation of
  // the previous step's residue, so at most the last one is transported.
  for (Int_t er = 0; er < numEvaporations; ++er)
    if (transported(EvaR[er]))
      residue_step = er;
  const Bool_t reacted = (residue_step >= 0);

  // The beam survives only on unreacted events (no vertex picked, or the
  // sampled reaction was energetically disallowed and the beam swept on to
  // the exit). On reacted events it was consumed at the vertex: N/A.
  if (Beam && !reacted)
    beam_energy_exit = recordExit(Beam, beam_stop_x, beam_stop_y, beam_stop_z,
                                  beam_stop_strip);

  for (Int_t er = 0; er < numEvaporations; ++er)
    if (transported(EvaP[er]))
      evap_energy_exit[er] =
          recordExit(EvaP[er], evap_stop_x[er], evap_stop_y[er],
                     evap_stop_z[er], evap_stop_strip[er]);
  // Superseded residues decayed in place — their exit slots stay -2 (N/A).
  if (reacted)
    residue_energy_exit[residue_step] =
        recordExit(EvaR[residue_step], residue_stop_x, residue_stop_y,
                   residue_stop_z, residue_stop_strip);
}

void Simulator::FinalizeEvent(Int_t eventIndex) { ComputeExitEnergies(); }
