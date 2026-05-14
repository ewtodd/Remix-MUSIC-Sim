// Vavilov energy-loss straggling sampler.
//
// Background: per-step straggling in MUSIC sits at κ ≈ 0.08 (Vavilov regime),
// where the Bohr Gaussian that catima exposes is the wrong distribution. The
// correct physics is asymmetric — a long high-loss tail from δ-electron
// transfers and no negative-loss left tail.
//
// We use the convolution decomposition of Yi & Han (NIM B 149, 1999) to keep
// the hot path cheap: Φ(λ_V; κ, β²) factors into the convolution of two
// distributions, one at β²=0 and one at β²=1. So instead of a 2-D (κ, β²)
// sampler we precompute two 1-D κ-indexed quantile tables (filled from
// ROOT::Math::VavilovAccurate at startup) and sum two table lookups at
// sample time. No per-step VavilovAccurate reconfiguration; no per-(κ,β²)
// cache invalidation logic.

#include "VavilovSampler.hpp"

#include <algorithm>
#include <cmath>

#include <Math/VavilovAccurate.h>
#include <TRandom.h>

namespace music {

namespace {

// Grid sizes. 64 log-spaced κ points covers the [1e-3, 10] Vavilov band
// with ≈0.07-decade resolution; 257 CDF points gives ≈0.4% quantile
// resolution. Storage is ~66 kB total — trivial.
constexpr int kNKappa = 64;
constexpr int kNU = 257;  // odd so we include 0.5 exactly; endpoints excluded

// Convert λ_L (the Landau parameter VavilovAccurate works in) to λ_V.
//     λ_V = κ · (λ_L + ln κ)
// Yi & Han's convolution holds in λ_V space, so we store quantiles in λ_V.
inline double LtoV(double lambdaL, double kappa) {
  return kappa * (lambdaL + std::log(kappa));
}

// Binary search for the highest grid index i such that grid[i] <= x, with
// linear extrapolation handled by the caller via the returned fractional
// position. Returns i ∈ [0, grid.size()-2]; the fractional offset t in [0,1]
// places x between grid[i] and grid[i+1].
void LerpIndex(const std::vector<double>& grid, double x, int& i, double& t) {
  const int n = static_cast<int>(grid.size());
  if (x <= grid.front()) { i = 0;        t = 0.0; return; }
  if (x >= grid.back())  { i = n - 2;    t = 1.0; return; }
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) {
    int mid = (lo + hi) / 2;
    if (grid[mid] <= x) lo = mid; else hi = mid;
  }
  i = lo;
  t = (x - grid[lo]) / (grid[lo + 1] - grid[lo]);
}

}  // namespace

const VavilovSampler& VavilovSampler::Instance() {
  static const VavilovSampler kSingleton;
  return kSingleton;
}

VavilovSampler::VavilovSampler() {
  // log-spaced κ grid in [kKappaMin, kKappaMax]
  log_kappa_grid_.resize(kNKappa);
  const double logKmin = std::log(kKappaMin);
  const double logKmax = std::log(kKappaMax);
  for (int i = 0; i < kNKappa; ++i) {
    log_kappa_grid_[i] = logKmin + (logKmax - logKmin) * i / (kNKappa - 1);
  }

  // Interior CDF grid: u_j = (j+1)/(kNU+1), j = 0 .. kNU-1.
  // Strictly inside (0,1) so VavilovAccurate.Quantile is well-defined.
  u_grid_.resize(kNU);
  for (int j = 0; j < kNU; ++j) {
    u_grid_[j] = (j + 1.0) / (kNU + 1.0);
  }

  // β² = 0 corresponds to the non-relativistic limit. VavilovAccurate's
  // valid β² range is (0,1]; we use a small ε to stay inside.
  BuildTable(1e-6, q_b0_, mean_b0_, var_b0_);
  BuildTable(1.0,  q_b1_, mean_b1_, var_b1_);
}

void VavilovSampler::BuildTable(double beta2,
                                std::vector<std::vector<double>>& q,
                                std::vector<double>& means,
                                std::vector<double>& variances) {
  q.assign(kNKappa, std::vector<double>(kNU));
  means.assign(kNKappa, 0.0);
  variances.assign(kNKappa, 0.0);

  ROOT::Math::VavilovAccurate vav;
  for (int i = 0; i < kNKappa; ++i) {
    const double kappa = std::exp(log_kappa_grid_[i]);
    vav.SetKappaBeta2(kappa, beta2);
    double m1 = 0.0;
    double m2 = 0.0;
    for (int j = 0; j < kNU; ++j) {
      // VavilovAccurate works in λ_L; convert to λ_V for the convolution.
      const double lambdaL = vav.Quantile(u_grid_[j]);
      const double lambdaV = LtoV(lambdaL, kappa);
      q[i][j] = lambdaV;
      m1 += lambdaV;
      m2 += lambdaV * lambdaV;
    }
    const double mean = m1 / kNU;
    means[i] = mean;
    variances[i] = std::max(1e-30, m2 / kNU - mean * mean);
  }
}

double VavilovSampler::QuantileAt(const std::vector<std::vector<double>>& q,
                                  double kappa, double u) const {
  int iK; double tK;
  int iU; double tU;
  LerpIndex(log_kappa_grid_, std::log(kappa), iK, tK);
  LerpIndex(u_grid_, u, iU, tU);
  const double q00 = q[iK    ][iU    ];
  const double q01 = q[iK    ][iU + 1];
  const double q10 = q[iK + 1][iU    ];
  const double q11 = q[iK + 1][iU + 1];
  const double qK0 = q00 * (1 - tU) + q01 * tU;
  const double qK1 = q10 * (1 - tU) + q11 * tU;
  return qK0 * (1 - tK) + qK1 * tK;
}

double VavilovSampler::ScalarAt(const std::vector<double>& s, double kappa) const {
  int iK; double tK;
  LerpIndex(log_kappa_grid_, std::log(kappa), iK, tK);
  return s[iK] * (1 - tK) + s[iK + 1] * tK;
}

double VavilovSampler::SampleStandardized(double kappa, double beta2,
                                          TRandom* rng) const {
  // Outside the Vavilov band: fall back to a standard normal. Caller will
  // scale by σ_E and add the catima mean, recovering the Bohr Gaussian
  // (κ ≫ 1) or the Landau-truncated-by-clamp (κ ≪ 1) behaviour the limits
  // naturally produce.
  if (!rng) return 0.0;
  if (kappa <= kKappaMin || kappa >= kKappaMax) {
    return rng->Gaus(0.0, 1.0);
  }
  beta2 = std::clamp(beta2, 0.0, 1.0);

  // Yi & Han Eq. (16): Φ(κ, β²) is the convolution of Φ((1−β²)κ, 0) and
  // Φ(β²κ, 1). Independent samples; sum them in λ_V space.
  const double k1 = (1.0 - beta2) * kappa;
  const double k2 = beta2          * kappa;

  double sample_lV = 0.0;
  double mean_lV   = 0.0;
  double var_lV    = 0.0;

  if (k1 > kKappaMin) {
    const double u = rng->Uniform();
    sample_lV += QuantileAt(q_b0_, k1, u);
    mean_lV   += ScalarAt(mean_b0_, k1);
    var_lV    += ScalarAt(var_b0_,  k1);
  }
  if (k2 > kKappaMin) {
    const double u = rng->Uniform();
    sample_lV += QuantileAt(q_b1_, k2, u);
    mean_lV   += ScalarAt(mean_b1_, k2);
    var_lV    += ScalarAt(var_b1_,  k2);
  }

  if (var_lV <= 0.0) return 0.0;
  return (sample_lV - mean_lV) / std::sqrt(var_lV);
}

}  // namespace music
