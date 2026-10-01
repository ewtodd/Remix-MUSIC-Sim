#ifndef ENERGYLOSS_HPP
#define ENERGYLOSS_HPP

#include <algorithm>
#include <cmath>
#include <mutex>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>

#include "VavilovSampler.hpp"
#include "catima/catima.h"

namespace music {
// Immutable, instance-owned physics settings. Keeping these values out of
// process globals makes worker setup deterministic and removes configuration
// data races. catima itself still owns a process-wide cache, so all calls into
// it are serialized through CatimaMutex(). Energy-loss tables make that lock a
// setup-time cost rather than a hot-loop cost.
struct PhysicsConfig {
  catima::Config mean;
  catima::Config straggling;
  Bool_t stragglingEnabled = kTRUE;
  Int_t stoppingModel = 0; // 0 = catima, 1 = SRIM, 2 = arithmetic mean
  std::string srimDir = ".";
  std::string srimGasTag;
};

std::mutex &CatimaMutex();

struct EnergyStepResult {
  Double_t finalEnergy = 0.0;
  Double_t distance = 0.0;
  Double_t depositedEnergy = 0.0;
  Bool_t stopped = kFALSE;
};
} // namespace music

class EnergyLoss {
public:
  EnergyLoss(Int_t A, Int_t Z, Double_t IonMass_MeV_per_c2,
             const catima::Material *gas, const music::PhysicsConfig &physics,
             Float_t dEdxScale = 1.0);
  ~EnergyLoss() = default;

  Double_t GetFinalEnergy(Double_t InitialEnergy, Double_t PathLength);
  // Forward propagation with per-step Vavilov-sampled straggling. At our κ ≈
  // 0.1 the Bohr Gaussian catima exposes is the wrong shape (it can sample Kf >
  // Ki, unphysical). We sample from Vavilov via music::VavilovSampler (Yi & Han
  // table): below the ROOT-supported band the lowest-kappa, Landau-like shape
  // is retained; above it the Gaussian limit is used. Final energy is clamped
  // to [0, Ki]. dEdxScale multiplies the mean only.
  // Pass rng=nullptr for the mean-only result.
  Double_t GetFinalEnergyStraggled(Double_t InitialEnergy, Double_t PathLength,
                                   TRandom *rng);
  music::EnergyStepResult TransportStep(Double_t InitialEnergy,
                                        Double_t PathLength, TRandom *rng);
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
  music::PhysicsConfig physics_;

  Double_t gas_ZoverA_ = -1.0;

  // Pre-tabulated stopping power and straggling so the hot path skips catima.
  static constexpr Int_t kNTable = 256;
  std::vector<Double_t> log_ki_grid_;
  std::vector<Double_t> eloss_per_cm_;
  std::vector<Double_t> sigma2_per_cm_;

  static constexpr Double_t c_cm_ns = 29.9792458;
};

#endif
