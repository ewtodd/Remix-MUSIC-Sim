#include "Particle.hpp"

Particle::Particle(TString particleName, Double_t mass, Int_t charge,
                   Bool_t saveTrajectory) {
  Name = particleName;
  Q = charge;
  Mass = mass;
  A = 0;
  Z = charge;
  NEexc = 1;
  SaveTrajectory = saveTrajectory;
  DoNotPropagate = false;

  AttColor = 1;
  AttStyle = 1;
  AttWidth = 2;
  CurrentExcState = 0;
  NumMedia = 0;
  P.SetName("P_" + particleName);
  P.SetCoords(mass, 0, 0, 0);
  RI = 0;
  TrPts = 0;
  X.SetName("X_" + particleName);
  X.SetCoords(0, 0, 0, 0);

  Eexc.assign(1, 0.0);
  IonInMedium.resize(MaxMedia);
  gas_ = nullptr;
  dEdxScale_ = 1.0;
  if (saveTrajectory) {
    TrT.resize(MaxPoints);
    TrX.resize(MaxPoints);
    TrY.resize(MaxPoints);
    TrZ.resize(MaxPoints);
    TrK.resize(MaxPoints);
    Trajectory = new TEveStraightLineSet();
  } else {
    Trajectory = nullptr;
  }
}

Particle::~Particle() { delete Trajectory; }

void Particle::Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ) {
  X.Boost(BetaX, BetaY, BetaZ);
  P.Boost(BetaX, BetaY, BetaZ);
}

// Note: does NOT modify Name.
void Particle::Copy(Particle *rhs) {
  if (this != rhs) {
    A = rhs->A;
    NEexc = rhs->NEexc;
    Eexc = rhs->Eexc;
    ProbExc = rhs->ProbExc;
    CurrentExcState = rhs->CurrentExcState;
    Mass = rhs->Mass;
    Q = rhs->Q;
    Z = rhs->Z;
    X = rhs->X;
    P = rhs->P;
    if (rhs->gas_)
      SetMedium(rhs->gas_, rhs->physics_, rhs->dEdxScale_);
  }
}

void Particle::CopyTrace(Int_t &NumPts, Float_t *t, Float_t *x, Float_t *y,
                         Float_t *z, Float_t *K) {
  if (!SaveTrajectory) {
    NumPts = 0;
    return;
  }
  Int_t TotPoints = TrPts;
  if (TrPts > MaxPoints)
    TotPoints = MaxPoints;
  for (Int_t p = 0; p < TotPoints; p++) {
    t[p] = TrT[p];
    x[p] = TrX[p];
    y[p] = TrY[p];
    z[p] = TrZ[p];
    K[p] = TrK[p];
  }
  NumPts = TotPoints;
}

void Particle::GetBeta(Double_t &BetaX, Double_t &BetaY, Double_t &BetaZ) {
  Double_t E = P.GetX0();
  BetaX = P.GetX1() / E;
  BetaY = P.GetX2() / E;
  BetaZ = P.GetX3() / E;
}

Int_t Particle::GetCurrentExcState() { return CurrentExcState; }

Double_t Particle::GetEexc() const {
  if (CurrentExcState >= 0 && CurrentExcState < NEexc)
    return Eexc[CurrentExcState];
  return 0.0;
}

Double_t Particle::GetEexc(Int_t ExcState) const {
  if (ExcState >= 0 && ExcState < NEexc)
    return Eexc[ExcState];
  return 0;
}

Double_t Particle::GetFinalEnergy(Int_t MediumID, Double_t InitE,
                                  Double_t PathLength) {
  if (MediumID >= 0 && MediumID < NumMedia) {
    Double_t FinalE = IonInMedium[MediumID]->GetFinalEnergy(InitE, PathLength);
    return (FinalE < 0) ? 0 : FinalE;
  }
  return InitE;
}

Double_t Particle::GetFinalEnergyStraggled(Int_t MediumID, Double_t InitE,
                                           Double_t PathLength, TRandom *rng) {
  if (MediumID >= 0 && MediumID < NumMedia) {
    Double_t FinalE =
        IonInMedium[MediumID]->GetFinalEnergyStraggled(InitE, PathLength, rng);
    return (FinalE < 0) ? 0 : FinalE;
  }
  return InitE;
}

music::EnergyStepResult Particle::TransportStep(Int_t MediumID, Double_t InitE,
                                                Double_t PathLength,
                                                TRandom *rng) {
  if (MediumID >= 0 && MediumID < NumMedia)
    return IonInMedium[MediumID]->TransportStep(InitE, PathLength, rng);
  music::EnergyStepResult result;
  result.finalEnergy = InitE;
  result.distance = std::max(0.0, PathLength);
  return result;
}

Double_t Particle::GetInitialEnergy(Int_t MediumID, Double_t FinalE,
                                    Double_t PathLength) {
  if (MediumID >= 0 && MediumID < NumMedia) {
    Double_t InitE =
        IonInMedium[MediumID]->GetInitialEnergy(FinalE, PathLength);
    return (InitE < 0) ? 0 : InitE;
  }
  return FinalE;
}

Double_t Particle::GetKE() const {
  Double_t KE = P.GetX0() - Mass - GetEexc();
  if (KE < 0.0) {
    KE = 0.0;
    std::cout << "Warning(" << Name << "): unrealistic negative KE -> set to 0."
              << std::endl;
  }
  return KE;
}

Double_t Particle::GetOptimumStepSize(Int_t MediumID, Double_t Energy) {
  if (MediumID >= 0 && MediumID < NumMedia)
    return 0.01 * IonInMedium[MediumID]->GetOptimumStepSize(Energy);
  return 0.1;
}

FourVector Particle::GetP() const { return P; }

void Particle::GetP(Double_t &P0, Double_t &P1, Double_t &P2,
                    Double_t &P3) const {
  P0 = P.GetX0();
  P1 = P.GetX1();
  P2 = P.GetX2();
  P3 = P.GetX3();
}

Double_t Particle::GetPathLength(Int_t MediumID, Double_t InitE,
                                 Double_t FinalE, Double_t DeltaT) {
  if (MediumID >= 0 && MediumID < NumMedia)
    return IonInMedium[MediumID]->GetPathLength(InitE, FinalE, DeltaT);
  return 0;
}

Double_t Particle::GetPhi() const {
  Double_t phi = std::atan2(P.GetX2(), P.GetX1());
  if (phi < 0)
    phi += 2 * M_PI;
  return phi;
}

Double_t Particle::GetPhiX() {
  Double_t phi = std::atan2(X.GetX2(), X.GetX1());
  if (phi < 0)
    phi += 2 * M_PI;
  return phi;
}

Double_t Particle::GetTheta() const {
  Double_t px = P.GetX1();
  Double_t py = P.GetX2();
  return std::atan2(std::sqrt(px * px + py * py), P.GetX3());
}

Double_t Particle::GetThetaX() {
  Double_t x = X.GetX1();
  Double_t y = X.GetX2();
  return std::atan2(std::sqrt(x * x + y * y), X.GetX3());
}

Double_t Particle::GetTimeOfFlight(Int_t MediumID) {
  if (MediumID >= 0 && MediumID < NumMedia)
    return IonInMedium[MediumID]->GetTimeOfFlight();
  return 0;
}

Double_t Particle::GetTimeOfFlight(Int_t MediumID, Float_t InitialEnergy,
                                   Float_t PathLength) {
  if (MediumID >= 0 && MediumID < NumMedia)
    return IonInMedium[MediumID]->GetTimeOfFlight(InitialEnergy, PathLength);
  return 0;
}

void Particle::GetTrajectoryAtt(Short_t &Color, Short_t &Style,
                                Short_t &Width) {
  Color = AttColor;
  Style = AttStyle;
  Width = AttWidth;
}

void Particle::GetX(Double_t &X0, Double_t &X1, Double_t &X2,
                    Double_t &X3) const {
  X0 = X.GetX0();
  X1 = X.GetX1();
  X2 = X.GetX2();
  X3 = X.GetX3();
}

void Particle::Print(std::ostream &log) {
  log << "|== Particle " << Name << " =================================|"
      << std::endl;
  log << "| mass = " << Mass << " MeV/c^2   Z = " << Q << " e\n"
      << "| Eexc = ";
  if (NEexc > 0) {
    for (Int_t n = 0; n < NEexc; n++) {
      log << Eexc[n];
      log << ((n < NEexc - 1) ? ", " : "\n");
    }
  } else {
    log << " 0 MeV" << std::endl;
  }
  log << "| KE = " << GetKE() << " MeV" << std::endl;
  log << "| ";
  X.Print(log);
  log << "| ";
  P.Print(log);
  if (SaveTrajectory) {
    log << "| Trajectory: " << Trajectory << " C=" << Trajectory->GetLineColor()
        << " W=" << Trajectory->GetLineWidth()
        << " S=" << Trajectory->GetLineStyle() << std::endl;
  } else {
    log << "| Trajectory object not saved." << std::endl;
  }
  if (NumMedia > 0 && gas_) {
    log << "| Stopping-power via "
        << (physics_.stoppingModel == 1   ? "SRIM tables"
            : physics_.stoppingModel == 2 ? "mean(catima, SRIM)"
                                          : "catima")
        << "; gas density = " << gas_->density()
        << " g/cm³, dE/dx scale = " << dEdxScale_ << std::endl;
  }
  log << "|==================================================|" << std::endl;
}

void Particle::ResetTrace() {
  if (!SaveTrajectory)
    return;
  Int_t TotPoints = TrPts;
  if (TrPts == 0)
    TotPoints = MaxPoints;
  for (Int_t p = 0; p < TotPoints; p++) {
    TrT[p] = -1000;
    TrX[p] = 0;
    TrY[p] = 0;
    TrZ[p] = -1000;
    TrK[p] = 0;
  }
  TrPts = 0;
}

void Particle::ResetKinematics() {
  P.SetCoords(Mass, 0, 0, 0);
  X.SetCoords(0, 0, 0, 0);
  Eexc.assign(1, 0.0);
  ProbExc.clear();
  NEexc = 1;
  CurrentExcState = 0;
}

void Particle::SetCurrentExcState(Int_t ExcState) {
  CurrentExcState = ExcState;
}

// If Prob is provided, it's the (un-normalised) selection weight per state;
// the cumulative probability ProbExc is built and normalised. If Prob=0, all
// excited states are equally likely.
void Particle::SetExcEnergies(Int_t count, Double_t *energies,
                              Double_t *probabilities) {
  if (count <= 0 || energies == nullptr)
    return;
  NEexc = count;
  Eexc.assign(energies, energies + count);
  ProbExc.assign(static_cast<size_t>(count + 1), 0.0);

  Double_t Norm = 0;
  for (Int_t n = 0; n < count; n++) {
    Norm += probabilities != nullptr ? probabilities[n] : 1.0 / count;
  }
  ProbExc[0] = 0.0;
  Double_t Cumulative = 0;
  for (Int_t n = 0; n < count; n++) {
    Cumulative += probabilities != nullptr ? probabilities[n] : 1.0 / count;
    ProbExc[n + 1] = Cumulative;
  }
  for (Int_t n = 0; n < count + 1; n++)
    ProbExc[n] /= Norm;

  CurrentExcState = 0;
}

void Particle::SampleExcitation(TRandom *rng) {
  if (!rng || NEexc <= 1 || ProbExc.size() != size_t(NEexc + 1)) {
    CurrentExcState = 0;
    return;
  }
  const Double_t draw = rng->Uniform();
  const auto upper = std::upper_bound(ProbExc.begin(), ProbExc.end(), draw);
  CurrentExcState =
      std::clamp(Int_t(upper - ProbExc.begin()) - 1, 0, NEexc - 1);
}

void Particle::SetExcEnergy(Double_t Ex) {
  Eexc.assign(1, Ex);
  ProbExc.clear();
  NEexc = 1;
  CurrentExcState = 0;
}

// Mass number A is derived from Mass (MeV/c²) / atomic mass unit.
void Particle::SetMedium(const catima::Material *gas,
                         const music::PhysicsConfig &physics,
                         Float_t dEdxScale) {
  if (gas == nullptr)
    return;
  if (Z <= 0) {
    // Neutral particles (e.g. neutrons) propagate via ExitWindow / kinematics
    // paths, not catima. Skip silently.
    return;
  }
  const Double_t amu_MeV = 931.49410242;
  Int_t A_derived = (Mass > 0.0) ? Int_t(std::round(Mass / amu_MeV)) : 0;
  if (A_derived <= 0) {
    std::cout << Name << ": cannot derive mass number from Mass=" << Mass
              << " MeV/c²" << std::endl;
    return;
  }
  if (A == 0)
    A = A_derived;
  NumMedia = 1;
  gas_ = gas;
  physics_ = physics;
  dEdxScale_ = dEdxScale;
  IonInMedium[0] =
      std::make_unique<EnergyLoss>(A_derived, Z, Mass, gas, physics, dEdxScale);
}

void Particle::SetP(FourVector V) {
  P.SetCoords(V.GetX0(), V.GetX1(), V.GetX2(), V.GetX3());
}

void Particle::SetP(Double_t P0, Double_t P1, Double_t P2, Double_t P3) {
  P.SetCoords(P0, P1, P2, P3);
}

void Particle::SetReactionIndex(Int_t reactionIndex) { RI = reactionIndex; }

void Particle::SetTracePoint(Float_t t, Float_t x, Float_t y, Float_t z,
                             Float_t K) {
  if (!SaveTrajectory)
    return;
  Int_t p = TrPts;
  if (p < MaxPoints) {
    TrT[p] = t;
    TrX[p] = x;
    TrY[p] = y;
    TrZ[p] = z;
    TrK[p] = K;
  } else {
    std::cout << "Warning: " << Name
              << " reached maximum number of trace points." << std::endl;
  }
  TrPts++;
}

// TEveStraightLineSet attributes can only be set after lines have been added
// (otherwise segfaults), so cache them here and apply later.
void Particle::SetTrajectoryAtt(Short_t Color, Short_t Style, Short_t Width) {
  AttColor = Color;
  AttStyle = Style;
  AttWidth = Width;
}

void Particle::SetX(Double_t X0, Double_t X1, Double_t X2, Double_t X3) {
  X.SetCoords(X0, X1, X2, X3);
}
