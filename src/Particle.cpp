#include "Particle.hpp"

Particle::Particle(TString Name, Double_t M, Int_t Q, Bool_t SaveTrajectory) {
  this->Name = Name;
  this->Q = Q;
  Mass = M;
  A = 0;
  Z = Q;
  NEexc = 1;
  this->SaveTrajectory = SaveTrajectory;
  DoNotPropagate = false;

  AttColor = 1;
  AttStyle = 1;
  AttWidth = 2;
  CurrentExcState = 0;
  NumMedia = 0;
  P.SetName("P_" + Name);
  P.SetCoords(M, 0, 0, 0);
  RI = 0;
  TrPts = 0;
  X.SetName("X_" + Name);
  X.SetCoords(0, 0, 0, 0);

  Eexc = new Double_t[NEexc];
  Eexc[0] = 0.0;
  PDF = 0;
  ProbExc = 0;
  IonInMedium = new EnergyLoss *[MaxMedia];
  for (Int_t m = 0; m < MaxMedia; m++)
    IonInMedium[m] = 0;
  gas_ = nullptr;
  dEdxScale_ = 1.0;
  TrT = new Float_t[MaxPoints];
  TrX = new Float_t[MaxPoints];
  TrY = new Float_t[MaxPoints];
  TrZ = new Float_t[MaxPoints];
  TrK = new Float_t[MaxPoints];
  Trajectory = new TEveStraightLineSet();
}

void Particle::Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ) {
  X.Boost(BetaX, BetaY, BetaZ);
  P.Boost(BetaX, BetaY, BetaZ);
}

// Note: does NOT modify Name.
void Particle::Copy(Particle *rhs) {
  if (this != rhs) {
    A = rhs->A;
    SetExcEnergies(rhs->NEexc, rhs->Eexc);
    Mass = rhs->Mass;
    Q = rhs->Q;
    Z = rhs->Z;
    X = rhs->X;
    P = rhs->P;
    if (rhs->gas_)
      SetMedium(rhs->gas_, rhs->dEdxScale_);
  }
}

void Particle::CopyTrace(Int_t &NumPts, Float_t *t, Float_t *x, Float_t *y,
                         Float_t *z, Float_t *K) {
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

Double_t Particle::GetEexc() {
  Double_t Eexc = 0;
  if (NEexc == 1) {
    Eexc = this->Eexc[0];
  } else if (PDF != 0 && ProbExc != 0 && this->Eexc != 0) {
    Double_t Rdm = PDF->Uniform();
    for (Int_t ne = 0; ne < NEexc; ne++) {
      if (Rdm >= ProbExc[ne] && Rdm < ProbExc[ne + 1]) {
        CurrentExcState = ne;
        Eexc = this->Eexc[ne];
        break;
      }
    }
  }
  return Eexc;
}

Double_t Particle::GetEexc(Int_t ExcState) {
  if (this->Eexc != 0 && ExcState >= 0 && ExcState < NEexc)
    return this->Eexc[ExcState];
  return 0;
}

// OBSOLETE — body retained as a no-op for legacy callers.
Double_t Particle::GetEnergyLoss(Int_t /*MediumID*/, Double_t /*InitE*/,
                                 Double_t /*PathLength*/) {
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

Double_t Particle::GetInitialEnergy(Int_t MediumID, Double_t FinalE,
                                    Double_t PathLength) {
  if (MediumID >= 0 && MediumID < NumMedia) {
    Double_t InitE =
        IonInMedium[MediumID]->GetInitialEnergy(FinalE, PathLength);
    return (InitE < 0) ? 0 : InitE;
  }
  return FinalE;
}

Double_t Particle::GetKE() {
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

FourVector Particle::GetP() { return P; }

void Particle::GetP(Double_t &P0, Double_t &P1, Double_t &P2, Double_t &P3) {
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

Double_t Particle::GetPhi() {
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

Double_t Particle::GetTheta() {
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

void Particle::GetX(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3) {
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
    log << "| Stopping-power via catima; gas density = " << gas_->density()
        << " g/cm³, dE/dx scale = " << dEdxScale_ << std::endl;
  }
  log << "|==================================================|" << std::endl;
}

void Particle::ResetTrace() {
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
  Eexc[0] = 0.0;
}

void Particle::SetCurrentExcState(Int_t ExcState) {
  CurrentExcState = ExcState;
}

// If Prob is provided, it's the (un-normalised) selection weight per state;
// the cumulative probability ProbExc is built and normalised. If Prob=0, all
// excited states are equally likely.
void Particle::SetExcEnergies(Int_t N, Double_t *Eexc, Double_t *Prob) {
  if (N <= 0 || Eexc == 0)
    return;
  NEexc = N;
  delete[] this->Eexc;
  delete[] ProbExc;
  delete PDF;
  this->Eexc = new Double_t[N];
  ProbExc = new Double_t[N + 1];

  Double_t Norm = 0;
  for (Int_t n = 0; n < N; n++) {
    this->Eexc[n] = Eexc[n];
    Norm += (Prob != 0) ? Prob[n] : 1.0 / N;
  }
  ProbExc[0] = 0.0;
  Double_t Cumulative = 0;
  for (Int_t n = 0; n < N; n++) {
    Cumulative += (Prob != 0) ? Prob[n] : 1.0 / N;
    ProbExc[n + 1] = Cumulative;
  }
  for (Int_t n = 0; n < N + 1; n++)
    ProbExc[n] /= Norm;

  PDF = new TRandom3();
  PDF->SetSeed();
}

void Particle::SetExcEnergy(Double_t Ex) { Eexc[0] = Ex; }

// Mass number A is derived from Mass (MeV/c²) / atomic mass unit.
void Particle::SetMedium(const catima::Material *gas, Float_t dEdxScale) {
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
  dEdxScale_ = dEdxScale;
  if (IonInMedium[0])
    delete IonInMedium[0];
  IonInMedium[0] = new EnergyLoss(A_derived, Z, Mass, gas, dEdxScale);
}

void Particle::SetP(FourVector V) {
  P.SetCoords(V.GetX0(), V.GetX1(), V.GetX2(), V.GetX3());
}

void Particle::SetP(Double_t P0, Double_t P1, Double_t P2, Double_t P3) {
  P.SetCoords(P0, P1, P2, P3);
}

void Particle::SetReactionIndex(Int_t RI) { this->RI = RI; }

void Particle::SetTracePoint(Float_t t, Float_t x, Float_t y, Float_t z,
                             Float_t K) {
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
