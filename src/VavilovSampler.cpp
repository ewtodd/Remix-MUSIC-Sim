#include "VavilovSampler.hpp"

#include <Math/Vavilov.h>

namespace music {

namespace {

constexpr Int_t kNKappa = 64;
constexpr Int_t kNBeta2 = 17;
constexpr Int_t kNU = 513;

// Return the lower interpolation index and fractional offset. Quantile-axis
// extrapolation is deliberate: it avoids artificial point masses caused by
// clamping uniform draws into an interior probability grid.
void LerpIndex(const std::vector<Double_t> &grid, Double_t x, Int_t &i,
               Double_t &t, Bool_t extrapolate = kFALSE) {
  const Int_t n = static_cast<Int_t>(grid.size());
  if (x <= grid.front()) {
    i = 0;
    t = extrapolate ? (x - grid[0]) / (grid[1] - grid[0]) : 0.0;
    return;
  }
  if (x >= grid.back()) {
    i = n - 2;
    t = extrapolate ? (x - grid[n - 2]) / (grid[n - 1] - grid[n - 2]) : 1.0;
    return;
  }
  const auto upper = std::upper_bound(grid.begin(), grid.end(), x);
  i = static_cast<Int_t>(upper - grid.begin()) - 1;
  t = (x - grid[i]) / (grid[i + 1] - grid[i]);
}

} // namespace

const VavilovSampler &VavilovSampler::Instance() {
  static const VavilovSampler singleton;
  return singleton;
}

VavilovSampler::VavilovSampler() {
  log_kappa_grid_.resize(kNKappa);
  beta2_grid_.resize(kNBeta2);
  u_grid_.resize(kNU);

  const Double_t logKmin = std::log(kKappaMin);
  const Double_t logKmax = std::log(kKappaMax);
  for (Int_t i = 0; i < kNKappa; ++i)
    log_kappa_grid_[i] = logKmin + (logKmax - logKmin) * i / (kNKappa - 1);
  for (Int_t i = 0; i < kNBeta2; ++i)
    beta2_grid_[i] = Double_t(i) / (kNBeta2 - 1);
  for (Int_t i = 0; i < kNU; ++i)
    u_grid_[i] = (i + 0.5) / kNU;

  quantiles_.resize(size_t(kNKappa) * kNBeta2 * kNU);
  ROOT::Math::VavilovAccurate vavilov;
  for (Int_t ik = 0; ik < kNKappa; ++ik) {
    const Double_t kappa = std::exp(log_kappa_grid_[ik]);
    for (Int_t ib = 0; ib < kNBeta2; ++ib) {
      const Double_t beta2 = beta2_grid_[ib];
      vavilov.SetKappaBeta2(kappa, beta2);
      const Double_t mean = ROOT::Math::Vavilov::Mean(kappa, beta2);
      const Double_t sigma =
          std::sqrt(ROOT::Math::Vavilov::Variance(kappa, beta2));
      for (Int_t iu = 0; iu < kNU; ++iu) {
        const size_t index =
            (size_t(ik) * kNBeta2 + size_t(ib)) * kNU + size_t(iu);
        quantiles_[index] = (vavilov.Quantile(u_grid_[iu]) - mean) / sigma;
      }
    }
  }
}

Double_t VavilovSampler::StandardizedQuantile(Double_t kappa, Double_t beta2,
                                              Double_t u) const {
  if (!std::isfinite(kappa) || !std::isfinite(beta2) || !std::isfinite(u))
    return 0.0;
  kappa = std::clamp(kappa, kKappaMin, kKappaMax);
  beta2 = std::clamp(beta2, 0.0, 1.0);
  u = std::clamp(u, 0.0, 1.0);
  Int_t ik = 0;
  Int_t ib = 0;
  Int_t iu = 0;
  Double_t tk = 0.0;
  Double_t tb = 0.0;
  Double_t tu = 0.0;
  LerpIndex(log_kappa_grid_, std::log(kappa), ik, tk);
  LerpIndex(beta2_grid_, beta2, ib, tb);
  LerpIndex(u_grid_, u, iu, tu, kTRUE);

  auto at = [&](Int_t k, Int_t b, Int_t q) {
    const size_t index = (size_t(k) * kNBeta2 + size_t(b)) * kNU + size_t(q);
    return quantiles_[index];
  };
  auto alongU = [&](Int_t k, Int_t b) {
    return at(k, b, iu) * (1.0 - tu) + at(k, b, iu + 1) * tu;
  };
  const Double_t q00 = alongU(ik, ib);
  const Double_t q01 = alongU(ik, ib + 1);
  const Double_t q10 = alongU(ik + 1, ib);
  const Double_t q11 = alongU(ik + 1, ib + 1);
  const Double_t q0 = q00 * (1.0 - tb) + q01 * tb;
  const Double_t q1 = q10 * (1.0 - tb) + q11 * tb;
  return q0 * (1.0 - tk) + q1 * tk;
}

Double_t VavilovSampler::SampleStandardized(Double_t kappa, Double_t beta2,
                                            TRandom *rng) const {
  if (!rng || !std::isfinite(kappa) || !std::isfinite(beta2))
    return 0.0;
  if (kappa > kKappaMax)
    return rng->Gaus(0.0, 1.0);

  // ROOT's accurate implementation supports kappa >= 1e-3. For still thinner
  // steps retain its lowest-kappa (strongly Landau-like) shape instead of
  // substituting a symmetric Gaussian.
  return StandardizedQuantile(kappa, beta2, rng->Uniform());
}

} // namespace music
