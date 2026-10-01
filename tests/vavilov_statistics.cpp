#include "VavilovSampler.hpp"

#include <Math/Vavilov.h>
#include <Math/VavilovAccurate.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <iostream>

int main() {
  const auto &sampler = music::VavilovSampler::Instance();
  constexpr std::array<double, 6> kappas = {1e-5, 1e-3, 0.01, 0.1, 1.0, 10.0};
  constexpr std::array<double, 3> beta2s = {0.0, 0.5, 1.0};
  constexpr std::array<double, 5> probabilities = {0.01, 0.1, 0.5, 0.9, 0.99};

  ROOT::Math::VavilovAccurate reference;
  for (double requestedKappa : kappas) {
    const double kappa = std::max(requestedKappa, 1e-3);
    for (double beta2 : beta2s) {
      reference.SetKappaBeta2(kappa, beta2);
      const double mean = ROOT::Math::Vavilov::Mean(kappa, beta2);
      const double sigma =
          std::sqrt(ROOT::Math::Vavilov::Variance(kappa, beta2));
      for (double probability : probabilities) {
        const double expected =
            (reference.Quantile(probability) - mean) / sigma;
        const double actual =
            sampler.StandardizedQuantile(requestedKappa, beta2, probability);
        const double tolerance =
            probability == 0.01 || probability == 0.99 ? 0.35 : 0.12;
        if (!std::isfinite(actual) ||
            std::fabs(actual - expected) > tolerance) {
          std::cerr << "Vavilov quantile mismatch: kappa=" << requestedKappa
                    << " beta2=" << beta2 << " p=" << probability
                    << " expected=" << expected << " actual=" << actual << '\n';
          return 1;
        }
      }
    }
  }

  // The high-kappa branch must be a standardized Gaussian.
  TRandom3 random(18731);
  constexpr int samples = 200000;
  double sum = 0.0;
  double sumSquares = 0.0;
  for (int i = 0; i < samples; ++i) {
    const double value = sampler.SampleStandardized(11.0, 0.5, &random);
    sum += value;
    sumSquares += value * value;
  }
  const double mean = sum / samples;
  const double variance = sumSquares / samples - mean * mean;
  if (std::fabs(mean) > 0.015 || std::fabs(variance - 1.0) > 0.025) {
    std::cerr << "Gaussian-limit moments failed: mean=" << mean
              << " variance=" << variance << '\n';
    return 1;
  }

  // Below ROOT's supported range, preserve the lowest-kappa Vavilov shape.
  TRandom3 lowA(991);
  TRandom3 lowB(991);
  for (int i = 0; i < 1000; ++i)
    if (sampler.SampleStandardized(1e-8, 0.5, &lowA) !=
        sampler.SampleStandardized(1e-3, 0.5, &lowB)) {
      std::cerr << "low-kappa clamp is not deterministic\n";
      return 1;
    }
  return 0;
}
