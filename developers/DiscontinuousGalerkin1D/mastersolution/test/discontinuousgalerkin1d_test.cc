/**
 * @file discontinuousgalerkin1d_test_mastersolution.cc
 * @brief NPDE homework "DiscontinuousGalerkin1D" code
 * @author Oliver Rietmann
 * @date 22.05.2019
 * @copyright Developed at ETH Zurich
 */

#include "../discontinuousgalerkin1d.h"

#include <gtest/gtest.h>

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

TEST(DiscontinuousGalerkin1D, solveTrafficFlow) {
  int Ml = 40;
  int Mr = 40;
  int N_half = Mr + Ml + 1;

  Eigen::VectorXd x_ref = Eigen::VectorXd::LinSpaced(N_half, -2.0, 2.0);

  Eigen::VectorXd u_ref(N_half);
  u_ref << -4.11013912654857e-18, -4.72660905295286e-18, 8.00581225662214e-18,
      4.11294391613077e-18, -2.2298012370052e-18, -1.05808440058776e-17,
      -5.3814646573403e-18, -5.00378603987006e-18, 2.6736532253624e-18,
      2.84094167941902e-18, -1.2523525927762e-18, 2.49080935940472e-18,
      -4.20177737702227e-18, 5.70615548094113e-18, 5.0685770325926e-18,
      -9.15307859643567e-19, -1.22119029795752e-18, 7.46818780493834e-19,
      -5.40276403793583e-18, 3.07749725287703e-18, 4.41889301544291e-18,
      -2.21988392021262e-18, 1.57077336758211e-18, -4.44237387106423e-18,
      8.62069924053091e-19, 9.25970198093325e-18, -1.03135121722185e-17,
      1.655713505636e-18, 4.06697517341438e-18, -8.24864055819056e-18,
      4.52246236255391e-18, 5.02435700190905e-18, -4.80750538913644e-18,
      9.15939210382569e-18, -1.1633673052399e-17, 4.26142989729195e-18,
      -1.71746476628415e-18, 1.34457582404414e-18, -6.79037234015735e-19,
      0.0219503971136849, 0.995324251295099, 0.972296207390145,
      0.950832864802787, 0.928466219406586, 0.905577311693144,
      0.882347925306832, 0.858875149626509, 0.835215768767552,
      0.811405161855874, 0.787466190828235, 0.763413733301244,
      0.739257109012437, 0.715001364837435, 0.690647841631483,
      0.666194170542805, 0.641633639164249, 0.616953578035746,
      0.592131736882971, 0.567127387298308, 0.541853584294137,
      0.516028406912738, 0.483971593087262, 0.458146415705863,
      0.432872612701692, 0.407868263117029, 0.383046421964254,
      0.358366360835751, 0.333805829457195, 0.309352158368517,
      0.284998635162565, 0.260742890987563, 0.236586266698756,
      0.212533809171765, 0.188594838144126, 0.164784231232448,
      0.141124850373491, 0.117652074693168, 0.0944226883068558,
      0.0715337805934138, 0.0491671351972135, 0.0277037926098554;

  DiscontinuousGalerkin1D::Solution solution = solveTrafficFlow();

  // test for correct length
  ASSERT_EQ(solution.x_.size(), solution.u_.size());
  ASSERT_EQ(solution.x_.size(), N_half);

  // test for correct values
  if (solution.x_.size() == N_half && solution.u_.size() == N_half) {
    double tol = 1.0e-5;
    ASSERT_NEAR(0.0, (solution.x_ - x_ref).lpNorm<Eigen::Infinity>(), tol);
    ASSERT_NEAR(0.0, (solution.u_ - u_ref).lpNorm<Eigen::Infinity>(), tol);
  }
}

}  // namespace DiscontinuousGalerkin1D::test
