#ifndef VavilovSampler_hpp_INCLUDED
#define VavilovSampler_hpp_INCLUDED

// Vavilov energy-loss straggling sampler using the convolution decomposition
// of Yi & Han (NIM B 149 (1999) 263–271): for any β² ∈ [0,1] and κ,
//     Φ(λ_V; κ, β²) = Φ(λ_V; (1−β²)κ, β²=0) ⋆ Φ(λ_V; β²κ, β²=1).
// A sample from Φ(κ, β²) is then the sum of two independent samples from
// the 1-D κ-indexed tables at β²=0 and β²=1. This reduces a 2-D (κ, β²)
// sampler to two 1-D lookups, and we never pay the per-step cost of
// VavilovAccurate construction or SetKappaBeta2.
//
// The 1-D tables themselves are pre-filled at first use from
// ROOT::Math::VavilovAccurate (≈3 ms total of startup work), so we're
// effectively riding ROOT's accurate Vavilov implementation while taking
// the Yi & Han convolution trick to avoid hot-path overhead.

#include <vector>

class TRandom;

namespace music {

class VavilovSampler {
public:
  // Singleton — the tables are heavy enough that we want exactly one copy.
  static const VavilovSampler& Instance();

  // Returns a standardised Vavilov deviate (zero mean, unit variance under
  // the precomputed table moments) for the given κ and β². Caller multiplies
  // by σ_E and adds the catima mean to get an energy-loss sample.
  // For κ outside the Vavilov-applicable band [κ_min, κ_max], falls back to
  // a Gaussian draw — there the Bohr limit (κ ≫ 1) or the Landau limit
  // (κ ≪ 1) is the relevant physics and Vavilov degenerates anyway.
  double SampleStandardized(double kappa, double beta2, TRandom* rng) const;

  static constexpr double kKappaMin = 1e-3;
  static constexpr double kKappaMax = 10.0;

private:
  VavilovSampler();

  // Fill the per-β² 1-D tables from VavilovAccurate. Internal use.
  void BuildTable(double beta2,
                  std::vector<std::vector<double>>& quantiles,
                  std::vector<double>& means,
                  std::vector<double>& variances);

  // Bilinear lookup: λ_V quantile at (κ, u).
  double QuantileAt(const std::vector<std::vector<double>>& q,
                    double kappa, double u) const;
  // Linear lookup of a scalar (mean / variance) at κ.
  double ScalarAt(const std::vector<double>& s, double kappa) const;

  // log-spaced κ grid in [kKappaMin, kKappaMax]
  std::vector<double> log_kappa_grid_;
  // Linearly spaced CDF grid in (0, 1), interior only (1/(N+1) … N/(N+1))
  std::vector<double> u_grid_;

  // Tables indexed as [κ_idx][u_idx]: λ_V quantile for that (κ, u). Two
  // tables: one at β² = 0, one at β² = 1. Yi & Han's decomposition only
  // ever needs these two endpoints.
  std::vector<std::vector<double>> q_b0_;
  std::vector<std::vector<double>> q_b1_;
  std::vector<double> mean_b0_, mean_b1_;
  std::vector<double> var_b0_,  var_b1_;
};

}  // namespace music

#endif
