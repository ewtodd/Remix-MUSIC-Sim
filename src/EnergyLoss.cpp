#include "EnergyLoss.hpp"
#include "VavilovSampler.hpp"

#include <TRandom.h>
#include <algorithm>
#include <cmath>
#include <iostream>

EnergyLoss::EnergyLoss(int A, int Z, double IonMass_MeV_per_c2,
                       const catima::Material* gas, float dEdxScale)
    : GoodELossFile(true),
      proj_(double(A), double(Z)),
      gas_(gas),
      A_(A),
      Z_(Z),
      IonMass_(IonMass_MeV_per_c2),
      dEdxScale_(dEdxScale),
      TOF_(0.0) {}

catima::Material EnergyLoss::LayerWithThickness(double PathLength_cm) const {
  catima::Material layer = *gas_;
  layer.thickness_cm(PathLength_cm);
  return layer;
}

// Kinetic energy (MeV total, not MeV/u) after traversing PathLength of gas.
// StepSize is ignored — catima integrates internally.
double EnergyLoss::GetFinalEnergy(double InitialEnergy, double PathLength, double /*StepSize*/) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return InitialEnergy;
  catima::Material layer = LayerWithThickness(PathLength);
  proj_.T = InitialEnergy / A_;
  double Eout_per_u = catima::energy_out(proj_, layer);
  double Eout = Eout_per_u * A_;
  double Eloss = (InitialEnergy - Eout) * dEdxScale_;
  double Efinal = InitialEnergy - Eloss;
  if (Efinal < 0.0)
    Efinal = 0.0;
  return Efinal;
}

// <Z/A> weighted by mass fraction of each constituent. Constant for a fixed
// gas, so cached on first use.
double EnergyLoss::GasZoverA() {
  if (gas_ZoverA_ > 0.0) return gas_ZoverA_;
  if (gas_ == nullptr) return 0.0;
  double zoa = 0.0;
  // Cast away constness — catima exposes ncomponents/get_element/weight_fraction
  // as non-const accessors despite reading only. Local pointer keeps the
  // member const-correct.
  catima::Material* mat = const_cast<catima::Material*>(gas_);
  for (int i = 0; i < mat->ncomponents(); ++i) {
    catima::Target t = mat->get_element(i);
    if (t.A > 0)
      zoa += mat->weight_fraction(i) * (t.Z / t.A);
  }
  gas_ZoverA_ = zoa;
  return gas_ZoverA_;
}

// Same as GetFinalEnergy but with per-step energy straggling. See header
// comment for the physics rationale (Bohr Gaussian wrong shape for our κ;
// switch to Vavilov via Yi & Han's convolution sampler).
double EnergyLoss::GetFinalEnergyStraggled(double InitialEnergy, double PathLength, TRandom* rng) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return InitialEnergy;
  catima::Material layer = LayerWithThickness(PathLength);
  proj_.T = InitialEnergy / A_;
  catima::Result r = catima::calculate(proj_, layer);
  const double Eout    = r.Eout    * A_;
  const double sigma_E = r.sigma_E * A_;
  double Eloss = (InitialEnergy - Eout) * dEdxScale_;

  if (sigma_E > 0.0 && rng) {
    // Compute κ and β² for Vavilov. Standard definitions (Yi & Han Eqs. 1–3):
    //   ξ        = (K/2) · z² · <Z/A> · ρ·t / β²     (Landau scale, MeV)
    //   ε_max    = 2 m_e c² β² γ² / (1 + 2γ m_e/M + (m_e/M)²)
    //   κ        = ξ / ε_max
    constexpr double K       = 0.307075;   // MeV cm²/g (Bethe-Bloch)
    constexpr double me_MeV  = 0.510998950;
    const double M     = (IonMass_ > 0.0) ? IonMass_ : (A_ * 931.49410242);
    const double gamma = 1.0 + InitialEnergy / M;
    const double beta2 = std::max(1e-9, 1.0 - 1.0 / (gamma * gamma));
    const double rho   = (gas_ != nullptr) ? gas_->density() : 0.0;  // g/cm³
    const double rho_t = rho * PathLength;                            // g/cm²
    const double zoa   = GasZoverA();
    double standardized;
    if (rho_t <= 0.0 || zoa <= 0.0) {
      standardized = rng->Gaus(0.0, 1.0);
    } else {
      const double meM   = me_MeV / M;
      const double xi    = 0.5 * K * Z_ * Z_ * zoa * rho_t / beta2;
      const double emax  = 2.0 * me_MeV * beta2 * gamma * gamma
                           / (1.0 + 2.0 * gamma * meM + meM * meM);
      const double kappa = (emax > 0.0) ? (xi / emax) : 0.0;
      standardized = music::VavilovSampler::Instance()
                       .SampleStandardized(kappa, beta2, rng);
    }
    Eloss += sigma_E * standardized;
  }

  // Clamp to physical loss range: 0 ≤ Eloss ≤ InitialEnergy. The lower clamp
  // rejects the residual unphysical "energy gain" tail; the upper clamp
  // handles cases where the sample would drive the particle through zero.
  if (Eloss < 0.0) Eloss = 0.0;
  if (Eloss > InitialEnergy) Eloss = InitialEnergy;
  return InitialEnergy - Eloss;
}

// Inverse of GetFinalEnergy via bisection (energy_out is monotonic in InitialEnergy).
double EnergyLoss::GetInitialEnergy(double FinalEnergy, double PathLength, double /*StepSize*/) {
  if (PathLength <= 0.0)
    return FinalEnergy;
  double lo = FinalEnergy;
  double hi = FinalEnergy + 1.0;
  while (GetFinalEnergy(hi, PathLength, 0.0) < FinalEnergy && hi < 1e6) {
    lo = hi;
    hi *= 2.0;
  }
  for (int i = 0; i < 60; ++i) {
    double mid = 0.5 * (lo + hi);
    double Ef = GetFinalEnergy(mid, PathLength, 0.0);
    if (Ef < FinalEnergy)
      lo = mid;
    else
      hi = mid;
    if (hi - lo < 1e-6)
      break;
  }
  return 0.5 * (lo + hi);
}

double EnergyLoss::GetEnergyLoss(double InitialEnergy, double PathLength) {
  return InitialEnergy - GetFinalEnergy(InitialEnergy, PathLength, 0.0);
}

double EnergyLoss::GetOptimumStepSize(double Energy) {
  if (Energy <= 0.0 || gas_ == nullptr)
    return 0.01;
  proj_.T = Energy / A_;
  double dEdx_per_gcm2 = catima::dedx(proj_, *gas_);
  double rho = gas_->density();
  if (dEdx_per_gcm2 <= 0.0 || rho <= 0.0)
    return 0.01;
  double dEdx_MeV_per_cm = dEdx_per_gcm2 * A_ * rho;
  return 0.01 * Energy / dEdx_MeV_per_cm;
}

double EnergyLoss::GetPathLength(double InitialEnergy, double FinalEnergy, double DeltaT) {
  if (FinalEnergy >= InitialEnergy)
    return 0.0;
  double step = GetOptimumStepSize(InitialEnergy);
  if (step <= 0.0)
    return 0.0;
  double E = InitialEnergy;
  double L = 0.0;
  double tof_ns = 0.0;
  while (E > FinalEnergy && L < 1e5) {
    double Enew = GetFinalEnergy(E, step, 0.0);
    if (Enew <= 0.0) {
      L += step;
      break;
    }
    double Emid = 0.5 * (E + Enew);
    double beta = std::sqrt(1.0 - std::pow(IonMass_ / (Emid + IonMass_), 2));
    double v = beta * c_cm_ns;
    tof_ns += step / v;
    L += step;
    E = Enew;
  }
  TOF_ = tof_ns;
  if (DeltaT > 0.0 && tof_ns > DeltaT) {
    L *= DeltaT / tof_ns;
    TOF_ = DeltaT;
  }
  return L;
}

double EnergyLoss::GetTimeOfFlight() { return TOF_; }

double EnergyLoss::GetTimeOfFlight(double InitialEnergy, double PathLength, double StepSize) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return 0.0;
  if (StepSize <= 0.0)
    StepSize = GetOptimumStepSize(InitialEnergy);
  double E = InitialEnergy;
  double remaining = PathLength;
  double tof_ns = 0.0;
  while (remaining > 0.0 && E > 0.0) {
    double dx = std::min(StepSize, remaining);
    double Enew = GetFinalEnergy(E, dx, 0.0);
    double Emid = 0.5 * (E + Enew);
    double beta = std::sqrt(1.0 - std::pow(IonMass_ / (Emid + IonMass_), 2));
    double v = beta * c_cm_ns;
    if (v <= 0.0) break;
    tof_ns += dx / v;
    E = Enew;
    remaining -= dx;
  }
  TOF_ = tof_ns;
  return tof_ns;
}
