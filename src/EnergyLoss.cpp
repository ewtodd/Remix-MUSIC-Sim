#include "EnergyLoss.hpp"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <mutex>
#include <set>
#include <sstream>
#include <vector>

namespace music {
// Defaults to catima's compiled default config; Simulator::loadCtrlFile
// overwrites it (inheriting low_energy from the main config, forcing a
// straggling-capable z_effective) once a control file is parsed.
catima::Config gStragglingConfig;
Bool_t gStragglingEnabled = kTRUE;
Int_t gStoppingModel = 0;
std::string gSrimDir;
std::string gSrimGasTag;
} // namespace music

namespace {
// Element symbols, index = Z. Only needed to name a SRIM table file.
const char *kElem[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U"};
const Int_t kNElem = Int_t(sizeof(kElem) / sizeof(kElem[0]));

// Parse a SRIM output table into (kinetic energy [MeV], dE/dx [MeV/cm]).
// The file lists electronic and nuclear stopping separately in the units named
// by its "Stopping Units" header line; both are summed, since the anode sees
// the total deposit.
Bool_t ReadSrimTable(const std::string &path, std::vector<Double_t> &E_MeV,
                     std::vector<Double_t> &dEdx_MeV_per_cm) {
  std::ifstream in(path);
  if (!in)
    return kFALSE;
  Double_t unit_to_MeV_per_cm = 0.0;
  Bool_t in_table = kFALSE;
  std::string line;
  while (std::getline(in, line)) {
    if (!in_table) {
      if (line.find("Stopping Units") != std::string::npos) {
        if (line.find("MeV/mm") != std::string::npos)
          unit_to_MeV_per_cm = 10.0;
        else if (line.find("MeV/cm") != std::string::npos)
          unit_to_MeV_per_cm = 1.0;
        else if (line.find("keV/um") != std::string::npos)
          unit_to_MeV_per_cm = 10.0; // keV/um == MeV/mm
        in_table = (unit_to_MeV_per_cm > 0.0);
      }
      continue;
    }
    std::istringstream ss(line);
    Double_t e = 0, el = 0, nu = 0;
    std::string unit;
    if (!(ss >> e >> unit >> el >> nu))
      continue;
    Double_t scale = 0.0;
    if (unit == "eV")
      scale = 1e-6;
    else if (unit == "keV")
      scale = 1e-3;
    else if (unit == "MeV")
      scale = 1.0;
    else if (unit == "GeV")
      scale = 1e3;
    else
      continue;
    E_MeV.push_back(e * scale);
    dEdx_MeV_per_cm.push_back((el + nu) * unit_to_MeV_per_cm);
  }
  return E_MeV.size() > 2;
}
} // namespace

EnergyLoss::EnergyLoss(Int_t A, Int_t Z, Double_t IonMass_MeV_per_c2,
                       const catima::Material *gas, Float_t dEdxScale)
    : GoodELossFile(true), proj_(Double_t(A), Double_t(Z)), gas_(gas), A_(A),
      Z_(Z), IonMass_(IonMass_MeV_per_c2), dEdxScale_(dEdxScale), TOF_(0.0) {
  BuildTables();
}

catima::Material EnergyLoss::LayerWithThickness(Double_t PathLength_cm) const {
  catima::Material layer = *gas_;
  layer.thickness_cm(PathLength_cm);
  return layer;
}

// Pre-tabulate (mean energy loss / cm) and (straggling variance / cm) as a
// function of total kinetic energy. Each grid point is one catima::calculate
// call on a thin reference layer; the per-cm rate is recovered by dividing
// out the layer thickness. Range covers every projectile/energy combination
// we plausibly simulate; outside that the hot-path interpolation clamps to
// the endpoints.
void EnergyLoss::BuildTables() {
  if (gas_ == nullptr)
    return;
  log_ki_grid_.resize(kNTable);
  eloss_per_cm_.resize(kNTable);
  sigma2_per_cm_.resize(kNTable);

  const Double_t Ki_min = 1e-4; // MeV total
  const Double_t Ki_max = 1e4;
  const Double_t log_min = std::log(Ki_min);
  const Double_t log_max = std::log(Ki_max);
  // Reference thickness small enough that the per-cm rate is well approximated
  // by (Ki − Eout)/dx, large enough that catima doesn't return numerical
  // garbage. 1e-5 cm of gas at ~mg/cm² density is well below any physical
  // step we take.
  const Double_t dx_ref = 1e-5;
  catima::Material layer = *gas_;
  layer.thickness_cm(dx_ref);

  for (Int_t i = 0; i < kNTable; ++i) {
    const Double_t logKi = log_min + (log_max - log_min) * i / (kNTable - 1);
    log_ki_grid_[i] = logKi;
    const Double_t Ki = std::exp(logKi);
    proj_.T = Ki / A_;
    // Mean energy loss from the default (atima14) config; straggling variance
    // from the straggling config, since atima14 returns sigma_E = 0.
    catima::Result rMean = catima::calculate(proj_, layer);
    catima::Result rStr =
        catima::calculate(proj_, layer, music::gStragglingConfig);
    const Double_t Eout = rMean.Eout * A_;
    const Double_t sigma_E = rStr.sigma_E * A_;
    eloss_per_cm_[i] = std::max(0.0, (Ki - Eout) / dx_ref);
    sigma2_per_cm_[i] = std::max(0.0, sigma_E * sigma_E / dx_ref);
  }

  // SRIM ("srim") replaces, and "mean" averages with, only the MEAN dE/dx; the
  // straggling variance above stays, because SRIM tables carry no variance. A
  // missing table is fatal rather than a silent fall back to catima: the
  // models differ by ~10% in helium, so quietly substituting one would
  // invalidate a dedx_scale calibrated on the other.
  if (music::gStoppingModel == 0)
    return;
  const std::string ion =
      (Z_ > 0 && Z_ < kNElem) ? (std::to_string(A_) + kElem[Z_]) : "";
  if (ion.empty()) {
    std::cerr << "musicsim ERROR: no element symbol for Z=" << Z_
              << "; cannot name a SRIM table" << std::endl;
    std::exit(1);
  }
  const std::string path =
      music::gSrimDir + "/" + ion + "_in_" + music::gSrimGasTag + ".srim";
  std::vector<Double_t> tE, tS;
  if (!ReadSrimTable(path, tE, tS)) {
    std::cerr << "musicsim ERROR: physics.stopping = srim but no usable table "
              << "at " << path << "\n  generate it with srim-cache, or set "
              << "physics.stopping = catima" << std::endl;
    std::exit(1);
  }
  for (Int_t i = 0; i < kNTable; ++i) {
    const Double_t Ki = std::exp(log_ki_grid_[i]);
    // Clamp outside the tabulated range; the grid deliberately spans far more
    // than SRIM tabulates.
    Double_t v;
    if (Ki <= tE.front())
      v = tS.front();
    else if (Ki >= tE.back())
      v = tS.back();
    else {
      size_t k =
          size_t(std::lower_bound(tE.begin(), tE.end(), Ki) - tE.begin());
      const Double_t f = (Ki - tE[k - 1]) / (tE[k] - tE[k - 1]);
      v = tS[k - 1] + f * (tS[k] - tS[k - 1]);
    }
    // "srim" takes the table value; "mean" averages it with the catima value
    // already sitting in the slot, per ApJ 983:142 sec 2.2.
    eloss_per_cm_[i] = std::max(
        0.0, music::gStoppingModel == 2 ? 0.5 * (eloss_per_cm_[i] + v) : v);
  }
  // One line per distinct table: EnergyLoss is constructed per particle per
  // worker, so an unguarded message would repeat hundreds of times.
  {
    static std::mutex m;
    static std::set<std::string> seen;
    std::lock_guard<std::mutex> lk(m);
    if (seen.insert(path).second)
      std::cout << "  [stopping] "
                << (music::gStoppingModel == 2 ? "mean(catima, SRIM)" : "SRIM")
                << " table for " << ion << ": " << tE.size() << " points from "
                << path << std::endl;
  }
}

Double_t EnergyLoss::InterpAt(const std::vector<Double_t> &vals,
                              Double_t Ki) const {
  if (log_ki_grid_.empty())
    return 0.0;
  const Double_t x = std::log(std::max(Ki, 1e-30));
  if (x <= log_ki_grid_.front())
    return vals.front();
  if (x >= log_ki_grid_.back())
    return vals.back();
  const Double_t t =
      (x - log_ki_grid_.front()) / (log_ki_grid_.back() - log_ki_grid_.front());
  const Double_t fpos = t * (kNTable - 1);
  const Int_t i = static_cast<Int_t>(fpos);
  const Double_t frac = fpos - i;
  return vals[i] * (1.0 - frac) + vals[i + 1] * frac;
}

Double_t EnergyLoss::GetFinalEnergy(Double_t InitialEnergy,
                                    Double_t PathLength) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return InitialEnergy;
  Double_t Eloss =
      InterpAt(eloss_per_cm_, InitialEnergy) * PathLength * dEdxScale_;
  if (Eloss < 0.0)
    Eloss = 0.0;
  if (Eloss > InitialEnergy)
    Eloss = InitialEnergy;
  return InitialEnergy - Eloss;
}

// <Z/A> weighted by mass fraction. Constant for a fixed gas; cached.
Double_t EnergyLoss::GasZoverA() {
  if (gas_ZoverA_ > 0.0)
    return gas_ZoverA_;
  if (gas_ == nullptr)
    return 0.0;
  Double_t zoa = 0.0;
  // catima exposes ncomponents/get_element/weight_fraction as non-const
  // accessors despite reading only — cast away constness locally so the
  // member stays const-correct.
  catima::Material *mat = const_cast<catima::Material *>(gas_);
  for (Int_t i = 0; i < mat->ncomponents(); ++i) {
    catima::Target t = mat->get_element(i);
    if (t.A > 0)
      zoa += mat->weight_fraction(i) * (t.Z / t.A);
  }
  gas_ZoverA_ = zoa;
  return gas_ZoverA_;
}

Double_t EnergyLoss::GetFinalEnergyStraggled(Double_t InitialEnergy,
                                             Double_t PathLength,
                                             TRandom *rng) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return InitialEnergy;
  // sigma_E for the step scales as √(dsigma²/dx · dx). dEdxScale only
  // multiplies the mean; straggling magnitude stays at catima's value.
  const Double_t dedx_per_cm = InterpAt(eloss_per_cm_, InitialEnergy);
  const Double_t sigma2_per_cm = InterpAt(sigma2_per_cm_, InitialEnergy);
  const Double_t sigma_E = std::sqrt(std::max(0.0, sigma2_per_cm * PathLength));
  Double_t Eloss = dedx_per_cm * PathLength * dEdxScale_;

  if (music::gStragglingEnabled && sigma_E > 0.0 && rng) {
    // κ and β² for Vavilov (Yi & Han Eqs. 1–3):
    //   ξ    = (K/2) · z² · <Z/A> · ρ·t / β²     (Landau scale, MeV)
    //   εmax = 2 mₑc² β² γ² / (1 + 2γ mₑ/M + (mₑ/M)²)
    //   κ    = ξ / εmax
    constexpr Double_t K = 0.307075; // MeV cm²/g (Bethe–Bloch)
    constexpr Double_t me_MeV = 0.510998950;
    const Double_t M = (IonMass_ > 0.0) ? IonMass_ : (A_ * 931.49410242);
    const Double_t gamma = 1.0 + InitialEnergy / M;
    const Double_t beta2 = std::max(1e-9, 1.0 - 1.0 / (gamma * gamma));
    const Double_t rho = (gas_ != nullptr) ? gas_->density() : 0.0;
    const Double_t rho_t = rho * PathLength;
    const Double_t zoa = GasZoverA();
    Double_t standardized;
    if (rho_t <= 0.0 || zoa <= 0.0) {
      standardized = rng->Gaus(0.0, 1.0);
    } else {
      const Double_t meM = me_MeV / M;
      const Double_t xi = 0.5 * K * Z_ * Z_ * zoa * rho_t / beta2;
      const Double_t emax = 2.0 * me_MeV * beta2 * gamma * gamma /
                            (1.0 + 2.0 * gamma * meM + meM * meM);
      const Double_t kappa = (emax > 0.0) ? (xi / emax) : 0.0;
      standardized = music::VavilovSampler::Instance().SampleStandardized(
          kappa, beta2, rng);
    }
    Eloss += sigma_E * standardized;
  }

  // Clamp to 0 ≤ Eloss ≤ InitialEnergy: reject the residual unphysical
  // "energy gain" tail and the through-zero tail.
  if (Eloss < 0.0)
    Eloss = 0.0;
  if (Eloss > InitialEnergy)
    Eloss = InitialEnergy;
  return InitialEnergy - Eloss;
}

// Inverse of GetFinalEnergy via bisection (energy_out is monotonic in Ki).
Double_t EnergyLoss::GetInitialEnergy(Double_t FinalEnergy,
                                      Double_t PathLength) {
  if (PathLength <= 0.0)
    return FinalEnergy;
  Double_t lo = FinalEnergy;
  Double_t hi = FinalEnergy + 1.0;
  while (GetFinalEnergy(hi, PathLength) < FinalEnergy && hi < 1e6) {
    lo = hi;
    hi *= 2.0;
  }
  for (Int_t i = 0; i < 60; ++i) {
    Double_t mid = 0.5 * (lo + hi);
    Double_t Ef = GetFinalEnergy(mid, PathLength);
    if (Ef < FinalEnergy)
      lo = mid;
    else
      hi = mid;
    if (hi - lo < 1e-6)
      break;
  }
  return 0.5 * (lo + hi);
}

Double_t EnergyLoss::GetEnergyLoss(Double_t InitialEnergy,
                                   Double_t PathLength) {
  return InitialEnergy - GetFinalEnergy(InitialEnergy, PathLength);
}

Double_t EnergyLoss::GetOptimumStepSize(Double_t Energy) {
  if (Energy <= 0.0 || gas_ == nullptr)
    return 0.01;
  proj_.T = Energy / A_;
  Double_t dEdx_per_gcm2 = catima::dedx(proj_, *gas_);
  Double_t rho = gas_->density();
  if (dEdx_per_gcm2 <= 0.0 || rho <= 0.0)
    return 0.01;
  // catima::dedx returns the TOTAL ion dE/dx in MeV/(g/cm^2), not a
  // per-nucleon value, so there is no A_ factor here.
  Double_t dEdx_MeV_per_cm = dEdx_per_gcm2 * rho;
  return 0.01 * Energy / dEdx_MeV_per_cm;
}

Double_t EnergyLoss::GetPathLength(Double_t InitialEnergy, Double_t FinalEnergy,
                                   Double_t DeltaT) {
  if (FinalEnergy >= InitialEnergy)
    return 0.0;
  Double_t step = GetOptimumStepSize(InitialEnergy);
  if (step <= 0.0)
    return 0.0;
  Double_t E = InitialEnergy;
  Double_t L = 0.0;
  Double_t tof_ns = 0.0;
  while (E > FinalEnergy && L < 1e5) {
    Double_t Enew = GetFinalEnergy(E, step);
    if (Enew <= 0.0) {
      L += step;
      break;
    }
    Double_t Emid = 0.5 * (E + Enew);
    Double_t beta = std::sqrt(1.0 - std::pow(IonMass_ / (Emid + IonMass_), 2));
    Double_t v = beta * c_cm_ns;
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

Double_t EnergyLoss::GetTimeOfFlight() { return TOF_; }

Double_t EnergyLoss::GetTimeOfFlight(Double_t InitialEnergy,
                                     Double_t PathLength) {
  if (PathLength <= 0.0 || InitialEnergy <= 0.0)
    return 0.0;
  // Pick an integration step locally — GetOptimumStepSize keeps it small
  // enough relative to dE/dx that per-slice dE stays in the linear regime of
  // the tabulated stopping power.
  const Double_t dxInt = GetOptimumStepSize(InitialEnergy);
  Double_t E = InitialEnergy;
  Double_t remaining = PathLength;
  Double_t tof_ns = 0.0;
  while (remaining > 0.0 && E > 0.0) {
    Double_t dx = std::min(dxInt, remaining);
    Double_t Enew = GetFinalEnergy(E, dx);
    Double_t Emid = 0.5 * (E + Enew);
    Double_t beta = std::sqrt(1.0 - std::pow(IonMass_ / (Emid + IonMass_), 2));
    Double_t v = beta * c_cm_ns;
    if (v <= 0.0)
      break;
    tof_ns += dx / v;
    E = Enew;
    remaining -= dx;
  }
  TOF_ = tof_ns;
  return tof_ns;
}
