#include "VavilovSampler.hpp"

#include <Math/Vavilov.h>
#include <Math/VavilovAccurate.h>
#include <TRandom3.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

int main() {
  const music::VavilovSampler &sampler = music::VavilovSampler::Instance();
  // 0.0137, 0.6, 2.2 and beta2 0.01, 0.27 fall between table nodes, so every
  // interpolation axis carries a fractional weight. MUSIC ions sit at
  // beta2 ~ 0.004-0.03.
  const Int_t nKappa = 9;
  const Double_t kappas[nKappa] = {1e-5, 1e-3, 0.01, 0.0137, 0.1,
                                   0.6,  1.0,  2.2,  10.0};
  const Int_t nBeta2 = 5;
  const Double_t beta2s[nBeta2] = {0.0, 0.01, 0.27, 0.5, 1.0};
  const Int_t nProbability = 5;
  const Double_t probabilities[nProbability] = {0.01, 0.1, 0.5, 0.9, 0.99};
  // Measured worst off-node errors are 0.0007 sigma (core) and 0.024 sigma
  // (1%/99% tails).
  const Double_t coreTolerance = 0.005;
  const Double_t tailTolerance = 0.06;

  ROOT::Math::VavilovAccurate reference;
  for (Int_t k = 0; k < nKappa; ++k) {
    const Double_t requestedKappa = kappas[k];
    const Double_t kappa = std::max(requestedKappa, 1e-3);
    for (Int_t b = 0; b < nBeta2; ++b) {
      const Double_t beta2 = beta2s[b];
      reference.SetKappaBeta2(kappa, beta2);
      const Double_t mean = ROOT::Math::Vavilov::Mean(kappa, beta2);
      const Double_t sigma =
          std::sqrt(ROOT::Math::Vavilov::Variance(kappa, beta2));
      for (Int_t p = 0; p < nProbability; ++p) {
        const Double_t probability = probabilities[p];
        const Double_t expected =
            (reference.Quantile(probability) - mean) / sigma;
        const Double_t actual =
            sampler.StandardizedQuantile(requestedKappa, beta2, probability);
        const Bool_t tail = probability < 0.05 || probability > 0.95;
        const Double_t tolerance = tail ? tailTolerance : coreTolerance;
        if (!std::isfinite(actual) ||
            std::fabs(actual - expected) > tolerance) {
          std::cerr << "Vavilov quantile mismatch: kappa=" << requestedKappa
                    << " beta2=" << beta2 << " p=" << probability
                    << " expected=" << expected << " actual=" << actual
                    << std::endl;
          return 1;
        }
      }
    }
  }

  // Sampled draws must stay standardized across the heavy-ion step regime.
  // Over 20 seeds at 400k draws the worst deviations were 0.0043 (mean) and
  // 0.025 (variance, kappa 0.0137).
  const Int_t nSampledKappa = 3;
  const Double_t sampledKappas[nSampledKappa] = {0.0137, 0.137, 2.2};
  const Int_t momentSamples = 400000;
  TRandom3 momentRandom(5150);
  for (Int_t k = 0; k < nSampledKappa; ++k) {
    Double_t sum = 0.0;
    Double_t sumSquares = 0.0;
    for (Int_t i = 0; i < momentSamples; ++i) {
      const Double_t value =
          sampler.SampleStandardized(sampledKappas[k], 0.01, &momentRandom);
      sum += value;
      sumSquares += value * value;
    }
    const Double_t mean = sum / momentSamples;
    const Double_t variance = sumSquares / momentSamples - mean * mean;
    if (std::fabs(mean) > 0.01 || std::fabs(variance - 1.0) > 0.05) {
      std::cerr << "sampled Vavilov moments failed: kappa=" << sampledKappas[k]
                << " mean=" << mean << " variance=" << variance << std::endl;
      return 1;
    }
  }

  // Whole sampled distributions against ROOT's CDF (Kolmogorov-Smirnov) at the
  // production regime: 10 um steps in 400-555 Torr He, light ions at the
  // kappa floor, 37Cl/87Rb at 0.03-1. The bound is the 0.1% critical value.
  const Int_t nShape = 4;
  const Double_t shapeKappas[nShape] = {1e-3, 0.03, 0.3, 1.0};
  const Double_t shapeBeta2s[nShape] = {0.006, 0.006, 0.01, 0.006};
  const Int_t shapeSamples = 100000;
  const Double_t ksBound =
      1.95 / std::sqrt(static_cast<Double_t>(shapeSamples));
  TRandom3 shapeRandom(777);
  std::vector<Double_t> cdf(shapeSamples);
  for (Int_t k = 0; k < nShape; ++k) {
    reference.SetKappaBeta2(shapeKappas[k], shapeBeta2s[k]);
    const Double_t mean =
        ROOT::Math::Vavilov::Mean(shapeKappas[k], shapeBeta2s[k]);
    const Double_t sigma = std::sqrt(
        ROOT::Math::Vavilov::Variance(shapeKappas[k], shapeBeta2s[k]));
    for (Int_t i = 0; i < shapeSamples; ++i)
      cdf[i] = reference.Cdf(
          mean + sigma * sampler.SampleStandardized(
                             shapeKappas[k], shapeBeta2s[k], &shapeRandom));
    std::sort(cdf.begin(), cdf.end());
    Double_t distance = 0.0;
    for (Int_t i = 0; i < shapeSamples; ++i) {
      const Double_t above = (i + 1.0) / shapeSamples - cdf[i];
      const Double_t below = cdf[i] - static_cast<Double_t>(i) / shapeSamples;
      distance = std::max(distance, std::max(above, below));
    }
    if (distance > ksBound) {
      std::cerr << "sampled Vavilov shape differs from ROOT: kappa="
                << shapeKappas[k] << " beta2=" << shapeBeta2s[k]
                << " KS distance=" << distance << " bound=" << ksBound
                << std::endl;
      return 1;
    }
  }

  // The high-kappa branch must be a standardized Gaussian.
  TRandom3 random(18731);
  const Int_t samples = 200000;
  Double_t sum = 0.0;
  Double_t sumSquares = 0.0;
  for (Int_t i = 0; i < samples; ++i) {
    const Double_t value = sampler.SampleStandardized(11.0, 0.5, &random);
    sum += value;
    sumSquares += value * value;
  }
  const Double_t mean = sum / samples;
  const Double_t variance = sumSquares / samples - mean * mean;
  if (std::fabs(mean) > 0.015 || std::fabs(variance - 1.0) > 0.025) {
    std::cerr << "Gaussian-limit moments failed: mean=" << mean
              << " variance=" << variance << std::endl;
    return 1;
  }

  // Below ROOT's supported range, preserve the lowest-kappa Vavilov shape.
  TRandom3 lowA(991);
  TRandom3 lowB(991);
  for (Int_t i = 0; i < 1000; ++i)
    if (sampler.SampleStandardized(1e-8, 0.5, &lowA) !=
        sampler.SampleStandardized(1e-3, 0.5, &lowB)) {
      std::cerr << "low-kappa clamp is not deterministic" << std::endl;
      return 1;
    }
  return 0;
}
