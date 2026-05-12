#ifndef EnergyLoss_hpp_INCLUDED
#define EnergyLoss_hpp_INCLUDED

#include <string>
#include "catima/catima.h"

class EnergyLoss {
public:
  EnergyLoss(int A, int Z, double IonMass_MeV_per_c2,
             const catima::Material* gas, float dEdxScale = 1.0);
  ~EnergyLoss() = default;

  double GetFinalEnergy(double InitialEnergy, double PathLength, double StepSize = 0.0);
  double GetInitialEnergy(double FinalEnergy, double PathLength, double StepSize = 0.0);
  double GetEnergyLoss(double InitialEnergy, double PathLength);

  double GetOptimumStepSize(double Energy);
  double GetPathLength(double InitialEnergy, double FinalEnergy, double DeltaT);
  double GetTimeOfFlight();
  double GetTimeOfFlight(double InitialEnergy, double PathLength, double StepSize);

  void SetIonMass(double IonMass) { IonMass_ = IonMass; }

  bool GoodELossFile;
  std::string FileName;

private:
  catima::Material LayerWithThickness(double PathLength_cm) const;

  catima::Projectile proj_;
  const catima::Material* gas_;
  int A_;
  int Z_;
  double IonMass_;
  double dEdxScale_;
  double TOF_;

  static constexpr double c_cm_ns = 29.9792458;
};

#endif
