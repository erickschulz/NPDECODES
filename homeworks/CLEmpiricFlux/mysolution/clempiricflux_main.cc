/**
 * @file clempiricflux_main.cc
 * @brief NPDE exam problem summer 2019 "CLEmpiricFlux" code
 * @author Oliver Rietmann
 * @date 18.07.2019
 * @copyright Developed at ETH Zurich
 */

#include <iostream>

#include "clempiricflux.h"
#include "solvecauchyproblem.h"
#include "uniformcubicspline.h"

using namespace CLEmpiricFlux;

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision,
                                       Eigen::DontAlignCols, ", ", "\n");

int main() {
  // Burgers equation
  double a = -1.1;
  double b = 1.1;
  auto f_lambda = [](double u) { return 0.5 * u * u; };
  auto M_lambda = [](double u) { return 1.0; };
  int n = 10;
  Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(n, a, b);
  UniformCubicSpline f(a, b, u.unaryExpr(f_lambda), u.unaryExpr(M_lambda));

  // initial data and its support
  auto u0 = [](double x) { return x < 0.0 ? 1.0 : 0.0; };

  double h = 0.05;  // spacial meshwidth
  double T = 1.0;   // final time

  // compute solution
  Eigen::VectorXd mu0 = computeInitVec(f, u0, h, T);
  Eigen::VectorXd muT = solveCauchyProblem(f, mu0, h, T);

  // print solution
  std::cout << "Your solution at time T = " << T << ":" << std::endl;
  std::cout << muT.transpose().format(CSVFormat) << std::endl;

  return 0;
}
