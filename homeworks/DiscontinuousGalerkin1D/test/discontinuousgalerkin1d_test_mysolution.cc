/**
 * @file discontinuousgalerkin1d_test_mysolution.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include "../mysolution/discontinuousgalerkin1d.h"

#include <Eigen/Core>

namespace DiscontinuousGalerkin1D::test {

TEST(DiscontinuousGalerkin1D, compBmat) {
  int Ml = 1;
  int Mr = 2;
  double h = 0.55;

  // to test
  Eigen::MatrixXd B(compBmat(Ml, Mr, h));

  // reference
  int N = 2 * (Ml + Mr + 1);
  Eigen::VectorXd d(N);
  for (int i = 0; i < N; i += 2) {
    d[i] = 0.55;
    d[i + 1] = 0.013864583333333;
  }
  Eigen::MatrixXd B_ref = d.asDiagonal();

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (B - B_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, G) {
  Eigen::VectorXd mu(6);
  mu << 0.1, 0.2, 0.3, 0.4, 0.3, 0.2;
  auto f = [](double x) { return 0.5 * x * x; };
  int Ml = 1;
  int Mr = 1;
  double h = 0.55;
  double tol = 1.0e-8;

  Eigen::VectorXd Gvec;
  Eigen::VectorXd Gvec_ref(6);

  auto Fv = [](double v, double w) { return v; };
  Gvec = G(mu, f, Fv, Ml, Mr, h);
  Gvec_ref << 0.155, 0.0395977083333333, 0.255, 0.129515833333333, -0.055,
      0.185347708333333;
  ASSERT_NEAR(0.0, (Gvec - Gvec_ref).lpNorm<Eigen::Infinity>(), tol);

  auto Fw = [](double v, double w) { return w; };
  Gvec = G(mu, f, Fw, Ml, Mr, h);
  Gvec_ref << 0.145, 0.0615977083333333, 0.055, 0.0937658333333333, -0.245,
      0.0423477083333333;
  ASSERT_NEAR(0.0, (Gvec - Gvec_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, dgcl) {
  Eigen::VectorXd mu0(6);
  mu0 << 0.1, 0.2, 0.3, 0.4, 0.3, 0.2;
  auto f = [](double x) { return 0.5 * x * x; };
  auto F = [](double v, double w) { return v; };
  int Ml = 1;
  int Mr = 1;
  double h = 0.55;
  double T = 1.0;
  unsigned int m = 2;

  Eigen::VectorXd mu =
      DiscontinuousGalerkin1D::dgcl(mu0, f, F, T, Ml, Mr, h, m);

  Eigen::VectorXd mu_ref(6);
  mu_ref << 0.524398178802972, 5.22582383256893, 0.788213539122826,
      17.3887446322133, -0.540429592450715, 17.1803412500935;

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (mu - mu_ref).lpNorm<Eigen::Infinity>(), tol);
}

TEST(DiscontinuousGalerkin1D, Feo) {
  Eigen::Vector4d v = {1.0, -1.0, 1.0, 0.3};
  Eigen::Vector4d w = {2.0, 2.0, 1.0, -2.0};

  Eigen::Vector4d F = v.binaryExpr(w, &Feo);
  Eigen::Vector4d F_ref = {-2.0, -4.25, 0.0, 0.21};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (F - F_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace DiscontinuousGalerkin1D::test
