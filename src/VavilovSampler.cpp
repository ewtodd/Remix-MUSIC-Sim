#include "VavilovSampler.hpp"

namespace music {

namespace {

// 64 log-spaced κ points covers [1e-3, 10] at ~0.07-decade resolution; 257
// CDF points gives ~0.4% quantile resolution. Storage ≈ 66 kB.
constexpr Int_t kNKappa = 64;
constexpr Int_t kNU = 257; // odd so we include 0.5 exactly; endpoints excluded

// Convert λ_L (Landau parameter used by VavilovAccurate) to λ_V:
//     λ_V = κ · (λ_L + ln κ)
// Yi & Han's convolution holds in λ_V space, so we store quantiles in λ_V.
inline Double_t LtoV(Double_t lambdaL, Double_t kappa) {
  return kappa * (lambdaL + std::log(kappa));
}

// Return the highest grid index i such that grid[i] <= x, plus the
// fractional offset t ∈ [0,1] placing x between grid[i] and grid[i+1].
void LerpIndex(const std::vector<Double_t> &grid, Double_t x, Int_t &i,
               Double_t &t) {
  const Int_t n = static_cast<Int_t>(grid.size());
  if (x <= grid.front()) {
    i = 0;
    t = 0.0;
    return;
  }
  if (x >= grid.back()) {
    i = n - 2;
    t = 1.0;
    return;
  }
  Int_t lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    Int_t mid = (lo + hi) / 2;
    if (grid[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  i = lo;
  t = (x - grid[lo]) / (grid[lo + 1] - grid[lo]);
}

} // namespace

const VavilovSampler &VavilovSampler::Instance() {
  static const VavilovSampler kSingleton;
  return kSingleton;
}

VavilovSampler::VavilovSampler() {
  log_kappa_grid_.resize(kNKappa);
  const Double_t logKmin = std::log(kKappaMin);
  const Double_t logKmax = std::log(kKappaMax);
  for (Int_t i = 0; i < kNKappa; ++i)
    log_kappa_grid_[i] = logKmin + (logKmax - logKmin) * i / (kNKappa - 1);

  // Interior CDF grid u_j = (j+1)/(kNU+1), strictly inside (0, 1) so
  // VavilovAccurate.Quantile is well-defined.
  u_grid_.resize(kNU);
  for (Int_t j = 0; j < kNU; ++j)
    u_grid_[j] = (j + 1.0) / (kNU + 1.0);

  // VavilovAccurate's valid β² range is (0, 1]; use a small ε for β² = 0.
  BuildTable(1e-6, q_b0_, mean_b0_, var_b0_);
  BuildTable(1.0, q_b1_, mean_b1_, var_b1_);
}

void VavilovSampler::BuildTable(Double_t beta2,
                                std::vector<std::vector<Double_t>> &q,
                                std::vector<Double_t> &means,
                                std::vector<Double_t> &variances) {
  q.assign(kNKappa, std::vector<Double_t>(kNU));
  means.assign(kNKappa, 0.0);
  variances.assign(kNKappa, 0.0);

  ROOT::Math::VavilovAccurate vav;
  for (Int_t i = 0; i < kNKappa; ++i) {
    const Double_t kappa = std::exp(log_kappa_grid_[i]);
    vav.SetKappaBeta2(kappa, beta2);
    Double_t m1 = 0.0;
    Double_t m2 = 0.0;
    for (Int_t j = 0; j < kNU; ++j) {
      const Double_t lambdaL = vav.Quantile(u_grid_[j]);
      const Double_t lambdaV = LtoV(lambdaL, kappa);
      q[i][j] = lambdaV;
      m1 += lambdaV;
      m2 += lambdaV * lambdaV;
    }
    const Double_t mean = m1 / kNU;
    means[i] = mean;
    variances[i] = std::max(1e-30, m2 / kNU - mean * mean);
  }
}

Double_t VavilovSampler::QuantileAt(const std::vector<std::vector<Double_t>> &q,
                                    Double_t kappa, Double_t u) const {
  Int_t iK;
  Double_t tK;
  Int_t iU;
  Double_t tU;
  LerpIndex(log_kappa_grid_, std::log(kappa), iK, tK);
  LerpIndex(u_grid_, u, iU, tU);
  const Double_t q00 = q[iK][iU];
  const Double_t q01 = q[iK][iU + 1];
  const Double_t q10 = q[iK + 1][iU];
  const Double_t q11 = q[iK + 1][iU + 1];
  const Double_t qK0 = q00 * (1 - tU) + q01 * tU;
  const Double_t qK1 = q10 * (1 - tU) + q11 * tU;
  return qK0 * (1 - tK) + qK1 * tK;
}

Double_t VavilovSampler::ScalarAt(const std::vector<Double_t> &s,
                                  Double_t kappa) const {
  Int_t iK;
  Double_t tK;
  LerpIndex(log_kappa_grid_, std::log(kappa), iK, tK);
  return s[iK] * (1 - tK) + s[iK + 1] * tK;
}

Double_t VavilovSampler::SampleStandardized(Double_t kappa, Double_t beta2,
                                            TRandom *rng) const {
  // Outside the Vavilov band: fall back to a standard normal. Caller scales
  // by σ_E and adds the catima mean, recovering Bohr (κ ≫ 1) or
  // Landau-truncated-by-clamp (κ ≪ 1) — the relevant physics at the limits.
  if (!rng)
    return 0.0;
  if (kappa <= kKappaMin || kappa >= kKappaMax)
    return rng->Gaus(0.0, 1.0);
  beta2 = std::clamp(beta2, 0.0, 1.0);

  // Yi & Han Eq. (16): Φ(κ, β²) = Φ((1−β²)κ, 0) ⋆ Φ(β²κ, 1). Independent
  // samples summed in λ_V space.
  const Double_t k1 = (1.0 - beta2) * kappa;
  const Double_t k2 = beta2 * kappa;

  Double_t sample_lV = 0.0;
  Double_t mean_lV = 0.0;
  Double_t var_lV = 0.0;

  if (k1 > kKappaMin) {
    const Double_t u = rng->Uniform();
    sample_lV += QuantileAt(q_b0_, k1, u);
    mean_lV += ScalarAt(mean_b0_, k1);
    var_lV += ScalarAt(var_b0_, k1);
  }
  if (k2 > kKappaMin) {
    const Double_t u = rng->Uniform();
    sample_lV += QuantileAt(q_b1_, k2, u);
    mean_lV += ScalarAt(mean_b1_, k2);
    var_lV += ScalarAt(var_b1_, k2);
  }

  if (var_lV <= 0.0)
    return 0.0;
  return (sample_lV - mean_lV) / std::sqrt(var_lV);
}

} // namespace music
