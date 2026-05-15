#ifndef VAVILOVSAMPLER_HPP
#define VAVILOVSAMPLER_HPP

// Vavilov energy-loss straggling sampler using Yi & Han's convolution
// decomposition (NIM B 149 (1999) 263–271):
//     Φ(λ_V; κ, β²) = Φ(λ_V; (1−β²)κ, 0) ⋆ Φ(λ_V; β²κ, 1)
// A 2-D (κ, β²) sampler becomes two 1-D κ-indexed table lookups. The 1-D
// tables are pre-filled at first use from ROOT::Math::VavilovAccurate.

#include <algorithm>
#include <cmath>
#include <vector>

#include <Math/VavilovAccurate.h>
#include <TRandom.h>

namespace music {

class VavilovSampler {
public:
  static const VavilovSampler &Instance();

  // Returns a standardised Vavilov deviate (zero mean, unit variance under
  // the precomputed table moments). Caller multiplies by σ_E and adds the
  // catima mean to get an energy-loss sample. Outside [kKappaMin, kKappaMax]
  // falls back to a Gaussian draw — the Bohr (κ ≫ 1) or Landau (κ ≪ 1) limit
  // is the relevant physics there and Vavilov degenerates anyway.
  Double_t SampleStandardized(Double_t kappa, Double_t beta2,
                              TRandom *rng) const;

  static constexpr Double_t kKappaMin = 1e-3;
  static constexpr Double_t kKappaMax = 10.0;

private:
  VavilovSampler();

  void BuildTable(Double_t beta2, std::vector<std::vector<Double_t>> &quantiles,
                  std::vector<Double_t> &means,
                  std::vector<Double_t> &variances);

  Double_t QuantileAt(const std::vector<std::vector<Double_t>> &q,
                      Double_t kappa, Double_t u) const;
  Double_t ScalarAt(const std::vector<Double_t> &s, Double_t kappa) const;

  std::vector<Double_t> log_kappa_grid_;
  std::vector<Double_t> u_grid_;

  // [κ_idx][u_idx] → λ_V quantile. Two tables: β² = 0 and β² = 1 (the only
  // endpoints Yi & Han's decomposition ever needs).
  std::vector<std::vector<Double_t>> q_b0_;
  std::vector<std::vector<Double_t>> q_b1_;
  std::vector<Double_t> mean_b0_, mean_b1_;
  std::vector<Double_t> var_b0_, var_b1_;
};

} // namespace music

#endif
