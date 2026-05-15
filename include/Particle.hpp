#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <cmath>
#include <iostream>

#include <TEveStraightLineSet.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TString.h>

#include "EnergyLoss.hpp"
#include "FourVector.hpp"

class Particle {
public:
  Particle(TString Name, Double_t M = 0, Int_t Q = 0,
           Bool_t SaveTrajectory = false);

  void Boost(Double_t BetaX, Double_t BetaY, Double_t BetaZ);
  void Copy(Particle *Other);
  void CopyTrace(Int_t &NumPts, Float_t *t, Float_t *x, Float_t *y, Float_t *z,
                 Float_t *K);
  void GetBeta(Double_t &BetaX, Double_t &BetaY, Double_t &BetaZ);
  Short_t GetColor() { return AttColor; }
  Int_t GetCurrentExcState();
  Double_t GetEexc();
  Double_t GetEexc(Int_t ExcState);
  Double_t GetEnergyLoss(Int_t MediumID, Double_t InitE, Double_t PathLength);
  Double_t GetFinalEnergy(Int_t MediumID, Double_t InitE, Double_t PathLength);
  Double_t GetFinalEnergyStraggled(Int_t MediumID, Double_t InitE,
                                   Double_t PathLength, TRandom *rng);
  Double_t GetInitialEnergy(Int_t MediumID, Double_t FinalE,
                            Double_t PathLength);
  Double_t GetKE();
  Double_t GetOptimumStepSize(Int_t MediumID, Double_t Energy);
  FourVector GetP();
  void GetP(Double_t &P0, Double_t &P1, Double_t &P2, Double_t &P3);
  Double_t GetPathLength(Int_t MediumID, Double_t InitE, Double_t FinalE,
                         Double_t DeltaT);
  Double_t GetPhi();
  Double_t GetPhiX();
  Double_t GetTheta();
  Double_t GetThetaX();
  Double_t GetTimeOfFlight(Int_t MediumID);
  Double_t GetTimeOfFlight(Int_t MediumID, Float_t InitialEnergy,
                           Float_t PathLength);
  void GetTrajectoryAtt(Short_t &Color, Short_t &Style, Short_t &Width);
  void GetX(Double_t &X0, Double_t &X1, Double_t &X2, Double_t &X3);
  void Print(std::ostream &log = std::cout);
  void ResetTrace();
  void ResetKinematics();
  void SetCurrentExcState(Int_t ExcState);
  void SetExcEnergies(Int_t N, Double_t *Eexc, Double_t *Prob = 0);
  void SetExcEnergy(Double_t Ex);
  void SetMedium(const catima::Material *gas, Float_t dEdxScale = 1.0);
  void SetP(FourVector P);
  void SetP(Double_t P0, Double_t P1, Double_t P2, Double_t P3);
  void SetReactionIndex(Int_t RI);
  void SetTracePoint(Float_t t, Float_t x, Float_t y, Float_t z, Float_t K);
  void SetX(Double_t X0, Double_t X1, Double_t X2, Double_t X3);
  void SetTrajectoryAtt(Short_t Color, Short_t Style = 1, Short_t Width = 1);

  Int_t A;
  Double_t Mass;
  static const Int_t MaxPoints = 50000;
  TString Name;
  Int_t NEexc;
  Int_t Q;
  Int_t RI;
  Bool_t SaveTrajectory;
  Bool_t DoNotPropagate;
  TEveStraightLineSet *Trajectory;
  Int_t Z;

private:
  // Trace storage (time, position, kinetic energy at each sampled point).
  Float_t *TrT;
  Float_t *TrX;
  Float_t *TrY;
  Float_t *TrZ;
  Float_t *TrK;
  Int_t TrPts;

  FourVector P;
  FourVector X;

  Double_t *Eexc;
  Double_t *ProbExc;
  TRandom3 *PDF;

  const catima::Material *gas_;
  Float_t dEdxScale_;
  EnergyLoss **IonInMedium;
  Int_t NumMedia;
  static const Int_t MaxMedia = 40;

  Int_t CurrentExcState;

  Short_t AttColor;
  Short_t AttStyle;
  Short_t AttWidth;
};

#endif
