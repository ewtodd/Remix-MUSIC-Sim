#include "Simulator.hpp"

void Simulator::SetBeamParticle(TString particleName, Int_t color,
                                Float_t dEdxScale) {
  Double_t m = NuF->GetMass(particleName.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(particleName.Data());
  Beam = new Particle(particleName, m, Z, /*SaveTrajectory=*/false);
  Beam->SetTrajectoryAtt(static_cast<Short_t>(color));
  Beam->SetMedium(&gas_, physics_, dEdxScale);
  if (ctf.Update) {
    TrackBeam->SetName(particleName.Data());
    TrackBeam->SetMainColor(static_cast<Color_t>(color));
    TrackBeam->SetPickable(kTRUE);
    Eve->AddElement(TrackBeam);
  }
}

void Simulator::SetTargetParticle(TString particleName) {
  Double_t m = NuF->GetMass(particleName.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(particleName.Data());
  Target = new Particle(particleName, m, Z);
}

void Simulator::SetCompoundParticle(TString particleName) {
  Double_t m = NuF->GetMass(particleName.Data(), "MeV/c^2");
  Int_t Z = NuF->GetZ(particleName.Data());
  Compound = new Particle(particleName, m, Z, /*SaveTrajectory=*/false);
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
      std::cout << "Maximum excitation energy of " << particleName
                << " (compound) = " << ExMax << " MeV" << std::endl;
  }
}

void Simulator::SetEvapResAndPart(TString residueName, Int_t residueColor,
                                  TString particleName, Int_t particleColor,
                                  Float_t dEdxScaleRes, Float_t dEdxScalePar) {
  if (numEvaporations >= maxEvaporations) {
    std::cout << "Warning: No more than " << maxEvaporations
              << " evaporation particles allowed." << std::endl;
    return;
  }

  // Evaporated particle (p, n, α, …).
  Double_t mp = NuF->GetMass(particleName.Data(), "MeV/c^2");
  Int_t Zp = NuF->GetZ(particleName.Data());
  particleName += std::to_string(numEvaporations);
  EvaP[numEvaporations] =
      new Particle(particleName, mp, Zp, /*SaveTrajectory=*/false);
  EvaP[numEvaporations]->SetTrajectoryAtt(static_cast<Short_t>(particleColor));
  EvaP[numEvaporations]->SetMedium(&gas_, physics_, dEdxScalePar);
  if (verbose_)
    EvaP[numEvaporations]->Print();
  if (ctf.Update) {
    TrackEvaP[numEvaporations]->SetName(particleName.Data());
    TrackEvaP[numEvaporations]->SetMainColor(
        static_cast<Color_t>(particleColor));
    TrackEvaP[numEvaporations]->SetPickable(kTRUE);
    Eve->AddElement(TrackEvaP[numEvaporations]);
  }
  // Evaporation residue (heavy product).
  Double_t mr = NuF->GetMass(residueName.Data(), "MeV/c^2");
  Int_t Zr = NuF->GetZ(residueName.Data());
  EvaR[numEvaporations] =
      new Particle(residueName, mr, Zr, /*SaveTrajectory=*/false);
  EvaR[numEvaporations]->SetTrajectoryAtt(static_cast<Short_t>(residueColor));
  EvaR[numEvaporations]->SetMedium(&gas_, physics_, dEdxScaleRes);
  if (verbose_)
    EvaR[numEvaporations]->Print();
  if (ctf.Update) {
    TrackEvaR[numEvaporations]->SetName(residueName.Data());
    TrackEvaR[numEvaporations]->SetMainColor(
        static_cast<Color_t>(residueColor));
    TrackEvaR[numEvaporations]->SetPickable(kTRUE);
    Eve->AddElement(TrackEvaR[numEvaporations]);
  }
  numEvaporations++;
}
