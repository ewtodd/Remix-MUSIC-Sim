#include "Simulator.hpp"

void Simulator::SetBeamParticle(TString Name, Int_t Color, Float_t dEdxScale) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  Beam = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  Beam->SetTrajectoryAtt((Short_t)Color);
  Beam->SetMedium(&gas_, dEdxScale);
  if (ctf.Update) {
    TrackBeam->SetName(Name.Data());
    TrackBeam->SetMainColor(Color);
    TrackBeam->SetPickable(kTRUE);
    Eve->AddElement(TrackBeam);
  }
}

void Simulator::SetTargetParticle(TString Name) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  Target = new Particle(Name, m, Z);
}

void Simulator::SetCompoundParticle(TString Name) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  Compound = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  // Once beam and target are configured, log the maximum compound excitation.
  if (Beam && Target) {
    Double_t mb = Beam->Mass;
    Double_t Kb = Kb_at_gas;
    Double_t pb = std::sqrt(2 * mb * Kb * (1 + Kb / (2 * mb)));
    Double_t Eb = std::sqrt(mb * mb + pb * pb);
    FourVector Pb("Pb", Eb, 0, 0, pb);
    FourVector Pt("Pt", Target->Mass, 0, 0, 0);
    FourVector Ptot("Total four-mom. in the lab");
    Ptot = Pb + Pt;
    Double_t ExMax = std::sqrt(Ptot * Ptot) - Compound->Mass;
    if (verbose_)
      std::cout << "Maximum excitation energy of " << Name
                << " (compound) = " << ExMax << " MeV" << std::endl;
  }
}

void Simulator::SetDecayDaughter1(TString Name, Int_t Color) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  DeDau1 = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  DeDau1->SetTrajectoryAtt((Short_t)Color);
  DeDau1->SetMedium(&gas_);
}

void Simulator::SetDecayDaughter2(TString Name, Int_t Color) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  DeDau2 = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  DeDau2->SetTrajectoryAtt((Short_t)Color);
  DeDau2->SetMedium(&gas_);
}

void Simulator::SetEvapResAndPart(TString ResName, Int_t ResColor,
                                  TString ParName, Int_t ParColor,
                                  Float_t dEdxScaleRes, Float_t dEdxScalePar) {
  if (numEvaporations >= maxEvaporations) {
    std::cout << "Warning: No more than " << maxEvaporations
              << " evaporation particles allowed." << std::endl;
    return;
  }

  // Evaporated particle (p, n, α, …).
  Double_t mp = NuF->GetMass(ParName.Data(), "MeV/c^2");
  Int_t Zp = NuF->GetZ(ParName.Data());
  ParName += std::to_string(numEvaporations);
  EvaP[numEvaporations] =
      new Particle(ParName, mp, Zp, /*SaveTrajectory=*/false);
  EvaP[numEvaporations]->SetTrajectoryAtt((Short_t)ParColor);
  EvaP[numEvaporations]->SetMedium(&gas_, dEdxScalePar);
  if (verbose_)
    EvaP[numEvaporations]->Print();
  if (ctf.Update) {
    TrackEvaP[numEvaporations]->SetName(ParName.Data());
    TrackEvaP[numEvaporations]->SetMainColor(ParColor);
    TrackEvaP[numEvaporations]->SetPickable(kTRUE);
    Eve->AddElement(TrackEvaP[numEvaporations]);
  }
  // Evaporation residue (heavy product).
  Double_t mr = NuF->GetMass(ResName.Data(), "MeV/c^2");
  Int_t Zr = NuF->GetZ(ResName.Data());
  EvaR[numEvaporations] =
      new Particle(ResName, mr, Zr, /*SaveTrajectory=*/false);
  EvaR[numEvaporations]->SetTrajectoryAtt((Short_t)ResColor);
  EvaR[numEvaporations]->SetMedium(&gas_, dEdxScaleRes);
  if (verbose_)
    EvaR[numEvaporations]->Print();
  if (ctf.Update) {
    TrackEvaR[numEvaporations]->SetName(ResName.Data());
    TrackEvaR[numEvaporations]->SetMainColor(ResColor);
    TrackEvaR[numEvaporations]->SetPickable(kTRUE);
    Eve->AddElement(TrackEvaR[numEvaporations]);
  }
  numEvaporations++;
}

void Simulator::SetHeavyParticle(TString Name, Int_t Color, Int_t NEexc,
                                 Double_t *Eexc) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  Heavy = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  Heavy->SetTrajectoryAtt((Short_t)Color);
  Heavy->SetMedium(&gas_);
  Heavy->SetExcEnergies(NEexc, Eexc);
}

void Simulator::SetLightParticle(TString Name, Int_t Color) {
  Double_t m = NuF->GetMass(Name.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(Name.Data());
  Light = new Particle(Name, m, Z, /*SaveTrajectory=*/false);
  Light->SetTrajectoryAtt((Short_t)Color);
  Light->SetMedium(&gas_);
}
