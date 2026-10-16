#include "Simulator.hpp"

namespace {
constexpr Double_t thresholdTolerance = 1e-9;

// Share of a branch's excitation weight at or below EneAvail. Branches cannot
// be forced, so only a talys table can be partly or wholly closed.
Double_t OpenFraction(const Simulator::ExcitationBranch &branch,
                      Double_t EneAvail) {
  return branch.residueExc == 3 ? branch.table.AllowedWeight(EneAvail) : 1.0;
}
} // namespace

void Simulator::SetInitialKinematics(Double_t Kbi) {
  Double_t mb = Beam->Mass;
  Double_t pb = std::sqrt(2 * mb * Kbi * (1 + Kbi / (2 * mb)));
  Double_t theta_b = 0;
  Double_t phi_b = 0;
  Double_t Eb = std::sqrt(mb * mb + pb * pb);

  Beam->SetP(Eb, pb * std::sin(theta_b) * std::cos(phi_b),
             pb * std::sin(theta_b) * std::sin(phi_b), pb * std::cos(theta_b));
  Beam->SetX(0, 0, 0, 0);
  if (PrintLevel > 0) {
    Log << "musicsim::SetInitialKinematics ************************************"
        << std::endl;
    Beam->Print(Log);
  }
}

// Establish the kinematics of the particles at the reaction point. If theta_CM
// and phi_CM are both -1, they're assigned randomly.
//
// VALIDATION: 2023-12-20 MLA and DSG compared musicsim.log output with the
// LISE++ kinematics calculator; results were consistent.
Int_t Simulator::SetReactionKinematics(Double_t Kbr, Double_t zr, Double_t tof,
                                       Double_t theta_CM, Double_t phi_CM) {
  if (numEvaporations <= 0)
    return 0;
  Bool_t reactionAllowed = kTRUE;
  if (PrintLevel > 0) {
    Log << "musicsim::SetReactionKinematics ***********************************"
        << std::endl;
  }
  Double_t mb = Beam->Mass;
  Double_t mt = Target->Mass;
  Double_t mc = Compound->Mass;

  Double_t pb = std::sqrt(2 * mb * Kbr * (1 + Kbr / (2 * mb)));
  Double_t theta_b = 0;
  Double_t phi_b = 0;
  Double_t Eb = std::sqrt(mb * mb + pb * pb);

  Double_t BetaX = pb * std::sin(theta_b) * std::cos(phi_b) / (Eb + mt);
  Double_t BetaY = pb * std::sin(theta_b) * std::sin(phi_b) / (Eb + mt);
  Double_t BetaZ = pb * std::cos(theta_b) / (Eb + mt);
  if (PrintLevel > 0) {
    Log << "Center-of-mass velocity (v/c):" << std::endl;
    Log << "BetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ
        << std::endl;
    Log << "--- Beam and target particles --------------------------------------"
        << std::endl;
  }

  Beam->SetP(Eb, pb * std::sin(theta_b) * std::cos(phi_b),
             pb * std::sin(theta_b) * std::sin(phi_b), pb * std::cos(theta_b));
  Beam->SetX(tof, 0, 0, zr);
  if (PrintLevel > 0)
    Beam->Print(Log);

  Target->SetP(mt, 0, 0, 0);
  Target->SetX(tof, 0, 0, zr);
  if (PrintLevel > 0)
    Target->Print(Log);

  FourVector Ptot = Beam->GetP() + Target->GetP();

  Compound->SetP(Ptot);
  Compound->SetX(tof, 0, 0, zr);
  Compound->SetExcEnergy(std::sqrt(Ptot * Ptot) - mc);
  if (PrintLevel > 0) {
    Log << "--- Compound particle ----------------------------------------------"
        << std::endl;
    Compound->Print(Log);
  }

  // A class event draws its channel first, then runs that channel's chain
  // in the leading slots; the rest hold idle particles.
  Int_t steps = numEvaporations;
  Int_t chosen = -1;
  const ReactionChannel *channel = nullptr;
  if (!channels_.empty()) {
    chosen = SelectChannel(std::sqrt(Ptot * Ptot));
    if (chosen < 0) {
      if (PrintLevel > 0)
        Log << "No channel of the class is open." << std::endl;
      reactionAllowed = kFALSE;
      steps = 0;
    } else {
      channel = &channels_[chosen];
      steps = static_cast<Int_t>(channel->res.size());
      for (Int_t k = 0; k < numEvaporations; ++k) {
        EvaP[k] = k < steps ? channel->evapParticle[k] : idleSlots_[k].first;
        EvaR[k] = k < steps ? channel->resParticle[k] : idleSlots_[k].second;
      }
    }
  }

  // Build the chain transactionally. Products stay non-propagating until every
  // configured step succeeds, so a late failure cannot leave a half-reaction
  // alongside a resumed incident beam.
  for (Int_t er = 0; er < numEvaporations; er++) {
    EvaP[er]->DoNotPropagate = true;
    EvaP[er]->ResetKinematics();
    EvaR[er]->DoNotPropagate = true;
    EvaR[er]->ResetKinematics();
  }

  TString reacstr = Beam->Name + "(" + Target->Name + "," + EvaP[0]->Name +
                    ")" + EvaR[0]->Name;
  Log << "numEvaporations = " << numEvaporations << std::endl;
  for (Int_t er = 0; er < steps; er++) {
    Double_t ml = EvaP[er]->Mass;
    Double_t mh = EvaR[er]->Mass;
    Double_t Q0 =
        (er == 0) ? (ml + mh - mb - mt) : (ml + mh - EvaR[er - 1]->Mass);

    Double_t EneAvail = std::sqrt(Ptot * Ptot) - ml - mh;

    if (PrintLevel > 0) {
      if (er > 0) {
        Log << "\n--- Secondary Reaction (" << er
            << ") -------------------------------------------" << std::endl;
        reacstr =
            EvaR[er - 1]->Name + "->" + EvaP[er]->Name + "+" + EvaR[er]->Name;
      } else {
        Log << "\n--- Primary Reaction ------------------------------------------------"
            << std::endl;
      }
      Log << reacstr << std::endl;
      Log << "Q0=" << Q0 << " MeV" << std::endl;
      Log << "Max energy avail.=" << EneAvail << " MeV" << std::endl;
    }

    if (EneAvail < -thresholdTolerance) {
      if (PrintLevel > 0)
        Log << "Reaction step " << er << " is below threshold." << std::endl;
      reactionAllowed = kFALSE;
      break;
    }
    EneAvail = std::max(0.0, EneAvail);

    // The branch is drawn before any excitation or angle draw. A class step
    // samples its own table above the floor that lets the chain finish.
    Int_t branchIndex = 0;
    Int_t angular = 0; // 0 isotropic, 1 Rutherford, 2 TALYS elastic
    Double_t thetaMinDeg = 1.0;
    Double_t Ex = 0;
    if (channel && channel->emitted.empty()) {
      angular = 2;
      thetaMinDeg = ctf.classThetaCmMinDeg;
    } else if (channel) {
      Bool_t sampled = kTRUE;
      Ex = channel->stripSteps[simStrip_][er].Sample(
          channel->stripFloor[simStrip_][er], EneAvail, Rdm, sampled);
      if (!sampled) {
        if (PrintLevel > 0)
          Log << "Reaction step " << er << " has no open state." << std::endl;
        reactionAllowed = kFALSE;
        break;
      }
    } else {
      branchIndex = SelectBranch(er, EneAvail);
      if (branchIndex < 0) {
        if (PrintLevel > 0)
          Log << "Reaction step " << er << " has no open branch." << std::endl;
        reactionAllowed = kFALSE;
        break;
      }
      const ExcitationBranch &branch = ctf.branches[er][branchIndex];
      if (!SampleExcitation(branch, er, EneAvail, Ex)) {
        reactionAllowed = kFALSE;
        break;
      }
      angular = branch.angularDist;
      thetaMinDeg = branch.thetaCmMinDeg;
    }

    EvaR[er]->SetExcEnergy(Ex);
    reaction_branch[er] = branchIndex;
    residue_excitation[er] = static_cast<Float_t>(Ex);

    // For the first step, honour user-specified angles; for later steps,
    // randomise. Uniform on the unit sphere → cos(θ) uniform on [-1, 1].
    Double_t stepTheta = theta_CM;
    Double_t stepPhi = phi_CM;
    if ((stepTheta == -1 && stepPhi == -1) || er > 0) {
      if (angular == 2) {
        // This strip's TALYS elastic distribution, also with one draw.
        stepTheta = stripElastic_[simStrip_].Sample(thetaMinDeg, Rdm);
      } else if (angular == 1) {
        // Rutherford. With u = sin^2(theta/2), dsigma ~ du/u^2, so sampling u
        // on [u_min, 1] with that weight inverts in closed form:
        //     1/u = 1/u_min - r (1/u_min - 1)
        Double_t hmin = std::sin(0.5 * thetaMinDeg * pi / 180.0);
        Double_t umin = hmin * hmin;
        Double_t r = Rdm->Uniform(0.0, 1.0);
        Double_t inv = 1.0 / umin - r * (1.0 / umin - 1.0);
        Double_t u = 1.0 / inv;
        if (u > 1.0)
          u = 1.0;
        stepTheta = 2.0 * std::asin(std::sqrt(u));
      } else {
        stepTheta = std::acos(Rdm->Uniform(-1.0, 1.0));
      }
      stepPhi = Rdm->Uniform(-pi, pi);
    }
    theta_cm[er] = static_cast<Float_t>(stepTheta * 180.0 / pi);
    phi_cm[er] = static_cast<Float_t>(stepPhi * 180.0 / pi);

    if (PrintLevel > 0) {
      Log << "branch=" << branchIndex << std::endl;
      Log << "Ex(" << EvaR[er]->Name << ")=" << Ex
          << " MeV\ntheta_cm=" << stepTheta * 180 / pi
          << "\nphi_cm=" << stepPhi * 180 / pi << std::endl;
      Log << "--- Outgoing particles (evap res = " << er
          << ") -------------------------------" << std::endl;
    }

    if (reactionAllowed) {
      // pf_CM is the outgoing momentum magnitude in the CM. Ptot² is
      // Lorentz-invariant.
      const Double_t invariantMass2 = Ptot * Ptot;
      const Double_t radicand = (invariantMass2 - std::pow(ml + mh + Ex, 2)) *
                                (invariantMass2 - std::pow(ml - mh - Ex, 2)) /
                                (4.0 * invariantMass2);
      if (radicand < -thresholdTolerance || !std::isfinite(radicand)) {
        reactionAllowed = kFALSE;
        break;
      }
      Double_t pf_CM = std::sqrt(std::max(0.0, radicand));

      Double_t plxCM = -pf_CM * std::sin(stepTheta) * std::cos(stepPhi);
      Double_t plyCM = -pf_CM * std::sin(stepTheta) * std::sin(stepPhi);
      Double_t plzCM = -pf_CM * std::cos(stepTheta);
      Double_t ElCM = std::sqrt(ml * ml + pf_CM * pf_CM);
      EvaP[er]->SetP(ElCM, plxCM, plyCM, plzCM);

      Double_t phxCM = pf_CM * std::sin(stepTheta) * std::cos(stepPhi);
      Double_t phyCM = pf_CM * std::sin(stepTheta) * std::sin(stepPhi);
      Double_t phzCM = pf_CM * std::cos(stepTheta);
      Double_t EhCM = std::sqrt((mh + Ex) * (mh + Ex) + pf_CM * pf_CM);
      EvaR[er]->SetP(EhCM, phxCM, phyCM, phzCM);

      if (PrintLevel > 0) {
        Log << "(((((((((( Before lorentz boost ))))))))))" << std::endl;
        EvaR[er]->Print(Log);
        EvaP[er]->Print(Log);
      }

      // Lorentz boost into the lab frame. Sign of -Beta is correct.
      EvaP[er]->Boost(-BetaX, -BetaY, -BetaZ);
      EvaR[er]->Boost(-BetaX, -BetaY, -BetaZ);

      // The residue's 4-momentum is the new Ptot for the next step (so a
      // secondary decay's "reaction" is just a decay of this residue).
      Ptot = EvaR[er]->GetP();

      // Light particle starts at the vertex; it will be propagated.
      EvaP[er]->SetX(tof, 0, 0, zr);
      // Residue starts at the vertex; we don't yet know if it'll be propagated
      // (decided in the next loop iteration).
      EvaR[er]->SetX(tof, 0, 0, zr);

      if (PrintLevel > 0) {
        Log << ")))))))))) After lorentz boost ((((((((((" << std::endl;
        EvaR[er]->Print(Log);
        EvaP[er]->Print(Log);
      }

      EvaR[er]->GetBeta(BetaX, BetaY, BetaZ);
      if (PrintLevel > 0) {
        Log << "Evap residue beta (v/c):" << std::endl;
        Log << "\tBetaX=" << BetaX << "  BetaY=" << BetaY << "  BetaZ=" << BetaZ
            << std::endl;
      }
    }
  }

  if (!reactionAllowed) {
    for (Int_t er = 0; er < numEvaporations; ++er) {
      EvaP[er]->ResetKinematics();
      EvaP[er]->SetX(tof, 0.0, 0.0, zr);
      EvaP[er]->DoNotPropagate = true;
      EvaR[er]->ResetKinematics();
      EvaR[er]->SetX(tof, 0.0, 0.0, zr);
      EvaR[er]->DoNotPropagate = true;
      // A rejected chain is wholly absent from truth output. In particular,
      // erase angles already sampled by an earlier successful step before a
      // later step failed its threshold check.
      theta_cm[er] = phi_cm[er] = kAngleUnset;
      residue_excitation[er] = -2.0f;
      reaction_branch[er] = -1;
      evap_energy[er] = residue_energy[er] = -2.0f;
      evap_theta[er] = evap_phi[er] = kAngleUnset;
      residue_theta[er] = residue_phi[er] = kAngleUnset;
    }
    return 0;
  }

  // Slots past the chain stay idle; only the last residue is transported.
  for (Int_t er = 0; er < numEvaporations; ++er) {
    EvaP[er]->DoNotPropagate = er >= steps;
    EvaR[er]->DoNotPropagate = (er != steps - 1);
  }
  reaction_channel = channel ? chosen : 0;

  // Fill the reaction-kinematics branches.
  Kbr = Beam->GetKE();
  for (Int_t er = 0; er < steps; er++) {
    residue_energy[er] = static_cast<Float_t>(EvaR[er]->GetKE());
    evap_energy[er] = static_cast<Float_t>(EvaP[er]->GetKE());
    evap_theta[er] = static_cast<Float_t>(EvaP[er]->GetTheta() * 180.0 / pi);
    evap_phi[er] = static_cast<Float_t>(EvaP[er]->GetPhi() * 180.0 / pi);
    residue_theta[er] = static_cast<Float_t>(EvaR[er]->GetTheta() * 180.0 / pi);
    residue_phi[er] = static_cast<Float_t>(EvaR[er]->GetPhi() * 180.0 / pi);
  }

  if (PrintLevel > 0) {
    Log << Form("beam: K=%.2f MeV  z_{r}=%.2f cm  tof=%.1f ns", Kbr, zr, tof)
        << std::endl;
    for (Int_t er = 0; er < numEvaporations; er++) {
      if (EvaP[er] && !EvaP[er]->DoNotPropagate)
        Log << Form(
                   "%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
                   EvaP[er]->Name.Data(), EvaP[er]->GetKE(),
                   EvaP[er]->GetTheta() * 180 / pi,
                   EvaP[er]->GetPhi() * 180 / pi)
            << std::endl;
      if (EvaR[er] && !EvaR[er]->DoNotPropagate)
        Log << Form(
                   "%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
                   EvaR[er]->Name.Data(), EvaR[er]->GetKE(),
                   EvaR[er]->GetTheta() * 180 / pi,
                   EvaR[er]->GetPhi() * 180 / pi)
            << std::endl;
    }
  }
  if (LabelKine) {
    LabelKine->Clear();
    LabelKine->AddText("Kinematics");
    Log << "*** Kinematics ***" << std::endl;
    LabelKine->AddText(
        Form("beam: K=%.2f MeV  z_{r}=%.2f cm  tof=%.1f ns", Kbr, zr, tof));
    for (Int_t er = 0; er < numEvaporations; er++) {
      if (EvaP[er] && !EvaP[er]->DoNotPropagate)
        LabelKine->AddText(Form(
            "%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
            EvaP[er]->Name.Data(), EvaP[er]->GetKE(),
            EvaP[er]->GetTheta() * 180 / pi, EvaP[er]->GetPhi() * 180 / pi));
      if (EvaR[er] && !EvaR[er]->DoNotPropagate)
        LabelKine->AddText(Form(
            "%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
            EvaR[er]->Name.Data(), EvaR[er]->GetKE(),
            EvaR[er]->GetTheta() * 180 / pi, EvaR[er]->GetPhi() * 180 / pi));
    }
    LabelKine->AddText(
        Form("#theta_{c.m.}=%.1f deg", static_cast<Double_t>(theta_cm[0])));
  }

  return 1;
}

// Draws the branch for step `er` by weight x open fraction at this EneAvail;
// -1 when every branch is closed. A single-branch step consumes no RNG draw,
// so branch-free controls keep their pre-branch event streams.
Int_t Simulator::SelectBranch(Int_t er, Double_t EneAvail) {
  const std::vector<ExcitationBranch> &branches = ctf.branches[er];
  const Int_t count = static_cast<Int_t>(branches.size());
  if (count <= 1)
    return 0;
  Double_t open[kMaxBranchesPerStep];
  Double_t totalWeight = 0.0;
  for (Int_t i = 0; i < count; ++i) {
    open[i] = branches[i].weight * OpenFraction(branches[i], EneAvail);
    totalWeight += open[i];
  }
  if (!(totalWeight > 0.0))
    return -1;
  Double_t pick = Rdm->Uniform(0.0, totalWeight);
  Int_t lastOpen = 0;
  for (Int_t i = 0; i < count; ++i) {
    if (open[i] <= 0.0)
      continue;
    if (pick < open[i])
      return i;
    pick -= open[i];
    lastOpen = i;
  }
  // Rounding can leave pick at the very top of the cumulative weight.
  return lastOpen;
}

// Samples the step-`er` residue excitation for the selected branch. Returns
// kFALSE when the branch admits no excitation at this EneAvail.
Bool_t Simulator::SampleExcitation(const ExcitationBranch &branch, Int_t er,
                                   Double_t EneAvail, Double_t &Ex) {
  Bool_t ok = kTRUE;
  switch (branch.residueExc) {
  case 1:
    Ex = 0.0;
    break;
  case 2:
    Ex = Rdm->Uniform(0.0, EneAvail);
    break;
  case 0:
    if (er + 1 >= numEvaporations ||
        EneAvail + thresholdTolerance < minEx[er + 1])
      return kFALSE;
    // A forced intermediate residue carries at least the invariant-mass
    // excess required by the rest of the configured decay chain.
    Ex = minEx[er + 1] >= EneAvail ? EneAvail
                                   : Rdm->Uniform(minEx[er + 1], EneAvail);
    break;
  case 3:
    Ex = branch.table.Sample(EneAvail, Rdm, ok);
    break;
  default:
    return kFALSE;
  }
  return ok;
}

// Non-relativistic CM-energy ranges, per strip. Used to estimate the kinematic
// reach of a configuration before the event loop runs. When EnergeticsLog is
// open and energeticsWritten_ is false, writes the full table to
// energetics.log; otherwise only computes CMEMax/CMEMin silently.
void Simulator::CalculateCMEnergyRange() {
  Double_t mb = Beam->Mass;
  Double_t mt = Target->Mass;
  Double_t Kb = Kb_at_gas;
  Double_t TotalLength = 0;
  for (Int_t i = 0; i < AnodeRows; i++)
    TotalLength += AnodeDZ[i][0];

  Double_t CME_beg = Kb * mt / (mt + mb);
  CMEMax = CME_beg;

  Double_t Kb_min = MeanBeamEnergyAfter(Kb, TotalLength);
  Double_t CME_end = Kb_min * mt / (mt + mb);
  CMEMin = CME_end;

  // Write energetics table to file on first call only.
  if (EnergeticsLog.is_open() && !energeticsWritten_) {
    EnergeticsLog << "Center-of-mass energy range covered in " << TotalLength
                  << "cm (MUSIC length):" << std::endl;
    EnergeticsLog << "  Ecom(initial) = " << CME_beg << " MeV" << std::endl;
    EnergeticsLog << "   Ecom(final) = " << CME_end << " MeV" << std::endl;
    EnergeticsLog << "Energetics for each segment:" << std::endl;
    EnergeticsLog.fill(' ');
    EnergeticsLog.width(2);
    EnergeticsLog << "i";
    EnergeticsLog.width(5);
    EnergeticsLog << "stp";
    EnergeticsLog.width(6);
    EnergeticsLog << "L[cm]";
    EnergeticsLog.width(10);
    EnergeticsLog << "Ecm_in";
    EnergeticsLog.width(10);
    EnergeticsLog << "DeltaEcm";
    EnergeticsLog.width(10);
    EnergeticsLog << "Kb_in";
    EnergeticsLog.width(10);
    EnergeticsLog << "DeltaKb" << std::endl;

    Kb = Kb_at_gas;
    CME_beg = CMEMax;
    for (Int_t i = 0; i < AnodeRows; i++) {
      Double_t Kb_in = Kb;
      Kb = MeanBeamEnergyAfter(Kb, AnodeDZ[i][0]);
      Double_t Kb_out = Kb;
      CME_end = Kb * mt / (mt + mb);
      EnergeticsLog.fill(' ');
      EnergeticsLog.width(2);
      EnergeticsLog << i;
      EnergeticsLog.width(5);
      EnergeticsLog << AnodeStpID[i][0];
      EnergeticsLog.width(6);
      EnergeticsLog << AnodeDZ[i][0];
      EnergeticsLog.precision(5);
      EnergeticsLog.width(10);
      EnergeticsLog << CME_beg;
      EnergeticsLog.width(10);
      EnergeticsLog << CME_beg - CME_end;
      EnergeticsLog.width(10);
      EnergeticsLog << Kb_in;
      EnergeticsLog.width(10);
      EnergeticsLog << Kb_in - Kb_out << std::endl;
      if (i + 1 < AnodeRows)
        CME_beg = CME_end;
    }
    energeticsWritten_ = true;
  }
}
