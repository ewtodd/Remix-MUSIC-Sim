#include "Simulator.hpp"

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
  Int_t ReactionAllowed = 1;
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

  // Assume all the particles will be propagated; reset their 4-vectors and Eex.
  for (Int_t er = 0; er < numEvaporations; er++) {
    EvaP[er]->DoNotPropagate = false;
    EvaP[er]->ResetKinematics();
    EvaR[er]->DoNotPropagate = false;
    EvaR[er]->ResetKinematics();
  }

  TString reacstr = Beam->Name + "(" + Target->Name + "," + EvaP[0]->Name +
                    ")" + EvaR[0]->Name;
  Log << "numEvaporations = " << numEvaporations << std::endl;
  for (Int_t er = 0; er < numEvaporations; er++) {
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

    if (EneAvail < 0) {
      // Step is not energetically allowed; mark this and all later steps as
      // non-propagating and stop.
      if (PrintLevel > 0)
        Log << "Negative EneAvail!\nThe following particles will NOT be propagated:"
            << std::endl;
      for (Int_t i = er; i < numEvaporations; i++) {
        EvaP[i]->DoNotPropagate = true;
        EvaR[i]->DoNotPropagate = true;
        if (PrintLevel > 0)
          Log << i << " " << EvaP[i]->Name << ", " << EvaR[i]->Name
              << std::endl;
      }
      break;
    } else if (er > 0) {
      // Reaction at this step is allowed; stop propagating the previous
      // residue.
      EvaR[er - 1]->DoNotPropagate = true;
    }

    Double_t Ex = 0;
    if (EneAvail > minEx[er])
      // Force highly excited states (favors energetically-allowed chains).
      Ex = Rdm->Uniform(2 * EneAvail / 3, EneAvail);
    else
      ReactionAllowed = 0;

    EvaR[er]->SetExcEnergy(Ex);

    // For the first step, honour user-specified angles; for later steps,
    // randomise. Uniform on the unit sphere → cos(θ) uniform on [-1, 1].
    if ((theta_CM == -1 && phi_CM == -1) || er > 0) {
      theta_CM = std::acos(Rdm->Uniform(-1.0, 1.0));
      phi_CM = Rdm->Uniform(-pi, pi);
    }

    if (PrintLevel > 0) {
      Log << "Ex(" << EvaR[er]->Name << ")=" << Ex
          << " MeV\ntheta_cm=" << theta_CM * 180 / pi
          << "\nphi_cm=" << phi_CM * 180 / pi << std::endl;
      Log << "--- Outgoing particles (evap res = " << er
          << ") -------------------------------" << std::endl;
    }

    if (ReactionAllowed) {
      // pf_CM is the outgoing momentum magnitude in the CM. Ptot² is
      // Lorentz-invariant.
      Double_t pf_CM = std::sqrt((Ptot * Ptot - std::pow(ml + mh + Ex, 2)) *
                                 (Ptot * Ptot - std::pow(ml - mh - Ex, 2)) /
                                 (4 * (Ptot * Ptot)));

      Double_t plxCM = -pf_CM * std::sin(theta_CM) * std::cos(phi_CM);
      Double_t plyCM = -pf_CM * std::sin(theta_CM) * std::sin(phi_CM);
      Double_t plzCM = -pf_CM * std::cos(theta_CM);
      Double_t ElCM = std::sqrt(ml * ml + pf_CM * pf_CM);
      EvaP[er]->SetP(ElCM, plxCM, plyCM, plzCM);

      Double_t phxCM = pf_CM * std::sin(theta_CM) * std::cos(phi_CM);
      Double_t phyCM = pf_CM * std::sin(theta_CM) * std::sin(phi_CM);
      Double_t phzCM = pf_CM * std::cos(theta_CM);
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

      theta_l[er] = (EvaP[er]->GetTheta()) * 180 / pi;
      phi_l[er] = (EvaP[er]->GetPhi()) * 180 / pi;

      // Light particle starts at the vertex; it will be propagated.
      EvaP[er]->SetX(tof, 0, 0, zr);
      // Residue starts at the vertex; we don't yet know if it'll be propagated
      // (decided in the next loop iteration).
      EvaR[er]->SetX(tof, 0, 0, zr);

      theta_h[er] = (EvaR[er]->GetTheta()) * 180 / pi;
      phi_h[er] = (EvaR[er]->GetPhi()) * 180 / pi;

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
    } else {
      // Not allowed: park the products at rest at the vertex and mark them
      // as non-propagating.
      EvaP[er]->SetP(ml, 0, 0, 0);
      EvaP[er]->SetX(tof, 0, 0, zr);
      EvaP[er]->DoNotPropagate = true;
      EvaR[er]->SetP(mh + Ex, 0, 0, 0);
      EvaR[er]->SetX(tof, 0, 0, zr);
      EvaR[er]->DoNotPropagate = true;
    }
  }

  // Fill the reaction-kinematics branches.
  Kbr = Beam->GetKE();
  for (Int_t er = 0; er < numEvaporations; er++) {
    this->theta_CM[er] = theta_CM * 180 / pi;
    this->phi_CM[er] = phi_CM * 180 / pi;
    // -2 sentinel = "step not physically realised" (energetically disallowed,
    // DoNotPropagate set). Same convention as Kl_exit/Kh_exit so analysis can
    // mask both arrays the same way.
    Kh[er] = EvaR[er]->DoNotPropagate ? -2.0f : (Float_t)EvaR[er]->GetKE();
    Kl[er] = EvaP[er]->DoNotPropagate ? -2.0f : (Float_t)EvaP[er]->GetKE();
    theta_l[er] = (EvaP[er]->GetTheta()) * 180 / pi;
    phi_l[er] = (EvaP[er]->GetPhi()) * 180 / pi;
    theta_h[er] = (EvaR[er]->GetTheta()) * 180 / pi;
    phi_h[er] = (EvaR[er]->GetPhi()) * 180 / pi;
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
    if (DeDau1 && DeDau2) {
      LabelKine->AddText(
          Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
               DeDau1->Name.Data(), DeDau1->GetKE(),
               DeDau1->GetTheta() * 180 / pi, DeDau1->GetPhi() * 180 / pi));
      LabelKine->AddText(
          Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
               DeDau2->Name.Data(), DeDau2->GetKE(),
               DeDau2->GetTheta() * 180 / pi, DeDau2->GetPhi() * 180 / pi));
    } else if (Heavy) {
      LabelKine->AddText(
          Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
               Heavy->Name.Data(), Heavy->GetKE(), Heavy->GetTheta() * 180 / pi,
               Heavy->GetPhi() * 180 / pi));
    }
    if (Light) {
      LabelKine->AddText(
          Form("%s: K=%.2f MeV  #theta_{lab}=%.1f deg  #phi_{lab}=%.1f deg",
               Light->Name.Data(), Light->GetKE(), Light->GetTheta() * 180 / pi,
               Light->GetPhi() * 180 / pi));
    }
    LabelKine->AddText(Form("#theta_{c.m.}=%.1f deg", this->theta_CM[0]));
  }

  return ReactionAllowed;
}

// Non-relativistic CM-energy ranges, per strip. Used to estimate the kinematic
// reach of a configuration before the event loop runs.
void Simulator::CalculateCMEnergyRange() {
  Double_t mb = Beam->Mass;
  Double_t mt = Target->Mass;
  Double_t Kb = Kb_at_gas;
  Float_t TotalLength = 0;
  for (Int_t i = 0; i < AnodeRows; i++)
    TotalLength += AnodeDZ[i][0];

  Double_t CME_beg = Kb * mt / (mt + mb);
  CMEMax = CME_beg;
  if (verbose_) {
    std::cout << "Center-of-mass energy range covered in " << TotalLength
              << "cm (MUSIC length):" << std::endl;
    std::cout << "  Ecom(initial) = " << CME_beg << " MeV" << std::endl;
  }
  Double_t Kb_min = Beam->GetFinalEnergy(0, Kb, TotalLength);
  Double_t CME_end = Kb_min * mt / (mt + mb);
  CMEMin = CME_end;
  if (verbose_) {
    std::cout << "   Ecom(final) = " << CME_end << " MeV" << std::endl;
    std::cout << "Energetics for each segment:" << std::endl;
    std::cout.fill(' ');
    std::cout.width(2);
    std::cout << "i";
    std::cout.width(5);
    std::cout << "stp";
    std::cout.width(6);
    std::cout << "L[cm]";
    std::cout.width(10);
    std::cout << "Ecm_in";
    std::cout.width(10);
    std::cout << "DeltaEcm";
    std::cout.width(10);
    std::cout << "Kb_in";
    std::cout.width(10);
    std::cout << "DeltaKb" << std::endl;
  }

  for (Int_t i = 0; i < AnodeRows; i++) {
    Double_t Kb_in = Kb;
    Kb = Beam->GetFinalEnergy(0, Kb, AnodeDZ[i][0]);
    Double_t Kb_out = Kb;
    CME_end = Kb * mt / (mt + mb);
    if (verbose_) {
      std::cout.fill(' ');
      std::cout.width(2);
      std::cout << i;
      std::cout.width(5);
      std::cout << AnodeStpID[i][0];
      std::cout.width(6);
      std::cout << AnodeDZ[i][0];
      std::cout.precision(5);
      std::cout.width(10);
      std::cout << CME_beg;
      std::cout.width(10);
      std::cout << CME_beg - CME_end;
      std::cout.width(10);
      std::cout << Kb_in;
      std::cout.width(10);
      std::cout << Kb_in - Kb_out << std::endl;
    }
    if (i + 1 < AnodeRows)
      CME_beg = CME_end;
  }
}

void Simulator::CalculateExcEnergyRange() {
  Double_t mb = Beam->Mass;
  Double_t mt = Target->Mass;
  Double_t mf = Compound->Mass;
  FourVector Pb("Pb");
  FourVector Pt("Pt", mt, 0, 0, 0);
  FourVector Ptot("Total four-mom. in the lab");
  Double_t Kb = Kb_at_gas;
  Float_t TotalLength = 0;
  for (Int_t i = 0; i < AnodeRows; i++)
    TotalLength += AnodeDZ[i][0];

  Double_t pb = std::sqrt(2 * mb * Kb * (1 + Kb / (2 * mb)));
  Double_t Eb = std::sqrt(mb * mb + pb * pb);
  Pb.SetCoords(Eb, 0, 0, pb);
  Ptot = Pb + Pt;
  EexcMax = std::sqrt(Ptot * Ptot) - mf;
  Double_t Eexc_beg = EexcMax;
  std::cout << "Full excitation energy range covered in " << TotalLength
            << "cm:\n Eexc(beg) = " << Eexc_beg << " MeV";

  Double_t Kb_min = Beam->GetFinalEnergy(0, Kb, TotalLength);
  pb = std::sqrt(2 * mb * Kb_min * (1 + Kb_min / (2 * mb)));
  Eb = std::sqrt(mb * mb + pb * pb);
  Pb.SetCoords(Eb, 0, 0, pb);
  Ptot = Pb + Pt;
  EexcMin = std::sqrt(Ptot * Ptot) - mf;
  Double_t Eexc_end = EexcMin;
  std::cout << "   Eexc(end) = " << Eexc_end << " MeV" << std::endl;
  std::cout << "Exc. energy range in each segment:" << std::endl;

  for (Int_t i = 0; i < AnodeRows; i++) {
    Kb = Beam->GetFinalEnergy(0, Kb, AnodeDZ[i][0]);
    pb = std::sqrt(2 * mb * Kb * (1 + Kb / (2 * mb)));
    Eb = std::sqrt(mb * mb + pb * pb);
    Pb.SetCoords(Eb, 0, 0, pb);
    Ptot = Pb + Pt;
    Eexc_end = std::sqrt(Ptot * Ptot) - mf;
    SegEexcRange[i] = Eexc_beg - Eexc_end;
    std::cout << i << "\t" << AnodeDZ[i][0] << " cm \t" << SegEexcRange[i]
              << " MeV" << std::endl;
    if (i + 1 < AnodeRows)
      Eexc_beg = Eexc_end;
  }
}

void Simulator::PrintEnergetics(Double_t Kb, Double_t **DeltaEB) {
  Double_t mb = Beam->Mass;
  Double_t mt = Target->Mass;
  Double_t mc = Compound->Mass;

  std::cout.width(5);
  std::cout << "Stp";
  std::cout.width(15);
  std::cout << "ExMax[MeV]";
  std::cout.width(15);
  std::cout << "ExMin[MeV]";
  std::cout.width(15);
  std::cout << "Exmax[MeV]";
  std::cout.width(15);
  std::cout << "Exmin[MeV]";
  std::cout.width(15);
  std::cout << "ECMmax[MeV]";
  std::cout.width(15);
  std::cout << "ECMmin[MeV]" << std::endl;

  for (Int_t stp = 0; stp < AnodeRows + 1; stp++) {
    std::cout.width(5);
    std::cout << stp;
    Double_t pb = std::sqrt(2 * mb * Kb * (1 + Kb / 2 / mb));
    Double_t Eb = std::sqrt(mb * mb + pb * pb);
    FourVector Ptot;
    Ptot.SetCoords(Eb + mt, 0, 0, pb);
    std::cout.width(15);
    std::cout.precision(4);
    std::cout << std::sqrt(Ptot * Ptot) - mc;
    Double_t pCM_max =
        std::sqrt((Ptot * Ptot - std::pow(mt + mb, 2)) *
                  (Ptot * Ptot - std::pow(mb - mt, 2)) / (4 * (Ptot * Ptot)));

    Kb -= DeltaEB[stp][AnodeCols];
    if (Kb > 0) {
      pb = std::sqrt(2 * mb * Kb * (1 + Kb / 2 / mb));
      Eb = std::sqrt(mb * mb + pb * pb);
      Ptot.SetCoords(Eb + mt, 0, 0, pb);
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << std::sqrt(Ptot * Ptot) - mc;

      Double_t Ecm = mt * Kb / (mb + mt);
      Double_t Kt = Ecm * (mt + mb) / mb;
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << Kt;
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << std::sqrt(Ptot * Ptot) - mc;

      Double_t pCM_min =
          std::sqrt((Ptot * Ptot - std::pow(mt + mb, 2)) *
                    (Ptot * Ptot - std::pow(mb - mt, 2)) / (4 * (Ptot * Ptot)));
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << std::sqrt(mt * mt + pCM_max * pCM_max) - mt;
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << std::sqrt(mt * mt + pCM_min * pCM_min) - mt;
      std::cout.width(15);
      std::cout.precision(4);
      std::cout << Kb * mt / (mb + mt) << std::endl;
    } else {
      std::cout.width(15);
      std::cout << "0";
      std::cout.width(15);
      std::cout << "0";
      std::cout.width(15);
      std::cout << "0" << std::endl;
    }
  }
}
