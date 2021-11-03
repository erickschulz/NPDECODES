/**
 * @file parametricfiniteelements_main.cc
 * @brief NPDE homework ParametricFiniteElements code
 * @author Amélie Loher
 * @date 04.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include "parametricfiniteelements.h"

using namespace ParametricFiniteElements;

int main() {
  unsigned int n = 3;
  Eigen::VectorXd mu((n + 1) * (n + 1));

  // Psi \in C^1([0,1]), Psi > 0
  auto Psi = [](double x) -> double { return x * x + 1.0; };

  // 1 <= alpha(x) <= 2
  auto alpha = [](Eigen::Vector2d x) -> double { return 3.0 / 2.0; };

  mu = geoThermSolve(n, alpha, Psi);

  // Surface Integral over Gamma_S of expansion coefficient vector mu
  double val = geoThermSurfInt(n, Psi, mu);

  std::cout
      << "Basis Expansion Coefficient vector of solution Variational Problem: "
      << mu << std::endl;
  std::cout << "Value of expansion coefficient vector mu integrated over "
               "Surface Gamma_S : "
            << val << std::endl;

  return 0;
}
