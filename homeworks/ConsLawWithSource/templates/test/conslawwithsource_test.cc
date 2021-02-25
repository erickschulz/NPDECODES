/**
 * @file conslawwithsource_test.cc
 * @brief NPDE exam TEST FILE
 * @author Oliver Rietmann
 * @date 20.07.2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../conslawwithsource.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <cmath>

namespace ConsLawWithSource::test {

TEST(ConsLawWithSource, godnfn) {
  int N = 10;

  // My own Godunov flux
  Eigen::VectorXd v(N);
  v << -3.0, -2.2, -5.0, 1.0, 0.0, 3.4, 3.4, 0.1, 3.6, -0.5;
  Eigen::VectorXd w(N);
  w << 5.0, -2.4, -2.0, 0.0, -0.5, -3.4, 0.5, 0.1, 4.5, -0.5;
  Eigen::VectorXd Fvw = v.binaryExpr(w, &godnfn);

  // Reference Godunov flux
  Eigen::VectorXd Fvw_ref(N);
  Fvw_ref << 1.0, 2.49071795328941, 2.13533528323661, 1.71828182845905,
      1.10653065971263, 26.564100047397, 26.564100047397, 1.00517091807565,
      32.998234443678, 1.10653065971263;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (Fvw_ref - Fvw).lpNorm<Eigen::Infinity>(), tol);
}

TEST(ConsLawWithSource, fluxdiffsource) {
  int N = 10;

  // Compute vector mu of cell averages
  const double PI = 3.14159265358979323846;
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(N, -PI, 2.0 * PI);
  Eigen::VectorXd mu = x.unaryExpr([](double x) { return std::sin(x); });

  // Use Burgers equation with a central flux function
  auto f = [](double u) { return 0.5 * u * u; };
  auto F = [f](double v, double w) { return 0.5 * (f(v) + f(w)); };
  auto s = [](double u) { return -u; };

  // Cell size
  double h = (x(N - 1) - x(0)) / N;

  // My own fluxdiffsource
  Eigen::VectorXd fluxdiffsource_values = fluxdiffsource(mu, F, s, h);

  // Reference fluxdiffsource
  Eigen::VectorXd fluxdiffsource_values_ref(N);
  fluxdiffsource_values_ref << 0.1875, -0.628709713905398, -1.0037097139054,
      0.0, 1.0037097139054, 0.628709713905398, 0.0, -0.628709713905398,
      -1.0037097139054, -0.1875;

  double tol = 1.0e-8;
  Eigen::VectorXd difference =
      fluxdiffsource_values_ref - fluxdiffsource_values;
  ASSERT_NEAR(0.0, difference.lpNorm<Eigen::Infinity>(), tol);
}

TEST(ConsLawWithSource, traceMass) {
  // Compute my total mass at different times
  auto u0 = [](double x) { return (0.0 <= x && x < 1.0) ? 1.0 : 0.0; };
  unsigned int N = 1500;
  Eigen::VectorXd m = ConsLawWithSource::traceMass(u0, N);

  // Pick three time steps from my total mass vector
  int M = m.size() - 1;
  Eigen::Vector3d m_values = {m(0), m(M / 2), m(M)};

  // The corresponding time steps from the analytical solution
  Eigen::Vector3d m_values_ref = {1.0, 0.22182064707654, 0.0492043994694548};

  double tol = 1.0e-2;
  ASSERT_NEAR(0.0, (m_values_ref - m_values).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace ConsLawWithSource::test
