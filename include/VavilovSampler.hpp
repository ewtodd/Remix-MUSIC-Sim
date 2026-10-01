#ifndef VAVILOVSAMPLER_HPP
#define VAVILOVSAMPLER_HPP

// Direct Vavilov energy-loss sampler. A standardized quantile table is built
// from ROOT::Math::VavilovAccurate over (κ, β², u), then interpolated in
// the hot path. Direct tabulation avoids the small-component discontinuity of
// the former convolution approximation and uses ROOT's analytic moments.

#include <algorithm>
#include <cmath>
#include <vector>

#include <Math/VavilovAccurate.h>
#include <TRandom.h>

namespace music {

class VavilovSampler {
public:
  static const VavilovSampler &Instance();

  // Returns a standardized Vavilov deviate. Below the accurate table range,
  // the thinnest supported Vavilov shape is used (Landau-like); above it, the
  // Gaussian thick-absorber limit is used.
  Double_t SampleStandardized(Double_t kappa, Double_t beta2,
                              TRandom *rng) const;

  // Deterministic quantile counterpart used by validation code. Kappa below
  // the ROOT-supported range uses the lowest supported Landau-like shape.
  Double_t StandardizedQuantile(Double_t kappa, Double_t beta2,
                                Double_t probability) const;

  static constexpr Double_t kKappaMin = 1e-3;
  static constexpr Double_t kKappaMax = 10.0;

private:
  VavilovSampler();

  std::vector<Double_t> log_kappa_grid_;
  std::vector<Double_t> beta2_grid_;
  std::vector<Double_t> u_grid_;
  // Flat [kappa][beta2][u] standardized quantile table.
  std::vector<Double_t> quantiles_;
};

} // namespace music

#endif
