#ifndef ENERGYLOSS_HPP
#define ENERGYLOSS_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>

#include "VavilovSampler.hpp"
#include "catima/catima.h"

namespace music {
// catima's atima14 z_effective model — used for the mean dE/dx because it
// matches LISE++ — returns sigma_E = 0: that code path never populates the
// energy-loss straggling variance. So the straggling magnitude (for the gas,
// windows, and degrader alike) is read from this separate Config, whose
// z_effective is a straggling-capable model (pierce_blann by default,
// overridable via the [physics] straggling_z_effective key). sigma_E is only
// weakly model-dependent, so pairing it with the atima14 mean is well-behaved.
// Set in Simulator::loadCtrlFile.
extern catima::Config gStragglingConfig;
// Master switch for energy-loss straggling (gas Vavilov sampling and the
// Gaussian window/degrader smearing alike). Set from the [physics] straggling
// key; default on. Off means every energy loss is the catima mean.
extern Bool_t gStragglingEnabled;
} // namespace music

class EnergyLoss {
public:
  EnergyLoss(Int_t A, Int_t Z, Double_t IonMass_MeV_per_c2,
             const catima::Material *gas, Float_t dEdxScale = 1.0);
  ~EnergyLoss() = default;

  Double_t GetFinalEnergy(Double_t InitialEnergy, Double_t PathLength);
  // Forward propagation with per-step Vavilov-sampled straggling. At our κ ≈
  // 0.1 the Bohr Gaussian catima exposes is the wrong shape (it can sample Kf >
  // Ki, unphysical). We sample from Vavilov via music::VavilovSampler (Yi & Han
  // convolution), falling back to Gaussian outside the band [1e-3, 10]. Final
  // energy is clamped to [0, Ki]. dEdxScale multiplies the mean only.
  // Pass rng=nullptr for the mean-only result.
  Double_t GetFinalEnergyStraggled(Double_t InitialEnergy, Double_t PathLength,
                                   TRandom *rng);
  Double_t GetInitialEnergy(Double_t FinalEnergy, Double_t PathLength);
  Double_t GetEnergyLoss(Double_t InitialEnergy, Double_t PathLength);

  Double_t GetOptimumStepSize(Double_t Energy);
  Double_t GetPathLength(Double_t InitialEnergy, Double_t FinalEnergy,
                         Double_t DeltaT);
  Double_t GetTimeOfFlight();
  Double_t GetTimeOfFlight(Double_t InitialEnergy, Double_t PathLength);

  void SetIonMass(Double_t IonMass) { IonMass_ = IonMass; }

  Bool_t GoodELossFile;
  TString FileName;

private:
  catima::Material LayerWithThickness(Double_t PathLength_cm) const;
  Double_t GasZoverA();
  void BuildTables();
  Double_t InterpAt(const std::vector<Double_t> &vals, Double_t Ki) const;

  catima::Projectile proj_;
  const catima::Material *gas_;
  Int_t A_;
  Int_t Z_;
  Double_t IonMass_;
  Double_t dEdxScale_;
  Double_t TOF_;

  Double_t gas_ZoverA_ = -1.0;

  // Pre-tabulated stopping power and straggling so the hot path skips catima.
  static constexpr Int_t kNTable = 256;
  std::vector<Double_t> log_ki_grid_;
  std::vector<Double_t> eloss_per_cm_;
  std::vector<Double_t> sigma2_per_cm_;

  static constexpr Double_t c_cm_ns = 29.9792458;
};

#endif
