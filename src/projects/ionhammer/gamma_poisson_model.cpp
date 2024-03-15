//
// Created by Vasiliy Ershov on 08/11/2016.
//

#include "gamma_poisson_model.hpp"

using namespace n_gamma_poisson_model;

std::array<double, 100000> PoissonGammaDistribution::log_gamma_integer_cache_ =
    []() -> std::array<double, 100000> {
  std::array<double, 100000> cache;
  for (size_t i = 0; i < cache.size(); ++i) {
    cache[i] = boost::math::lgamma(i + 1);
  }
  return cache;
}();
