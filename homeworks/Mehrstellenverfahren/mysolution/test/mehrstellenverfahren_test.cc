/**
 * @file
 * @brief NPDE homework TEMPLATE MAIN FILE
 * @author Tobas Rohner
 * @date 25-03-2022
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../mehrstellenverfahren.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <iostream>

namespace mehrstellenverfahren::test {

TEST(mehrstellenerfahren, compMehrstellenA) {
  const unsigned int M = 4;
  Eigen::MatrixXd A = compMehrstellenA(M);
  Eigen::MatrixXd A_ref(M * M, M * M);
  A_ref << 20, -4, 0, 0, -4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 20, -4, 0,
      -1, -4, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 20, -4, 0, -1, -4, -1, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, -4, 20, 0, 0, -1, -4, 0, 0, 0, 0, 0, 0, 0, 0, -4,
      -1, 0, 0, 20, -4, 0, 0, -4, -1, 0, 0, 0, 0, 0, 0, -1, -4, -1, 0, -4, 20,
      -4, 0, -1, -4, -1, 0, 0, 0, 0, 0, 0, -1, -4, -1, 0, -4, 20, -4, 0, -1, -4,
      -1, 0, 0, 0, 0, 0, 0, -1, -4, 0, 0, -4, 20, 0, 0, -1, -4, 0, 0, 0, 0, 0,
      0, 0, 0, -4, -1, 0, 0, 20, -4, 0, 0, -4, -1, 0, 0, 0, 0, 0, 0, -1, -4, -1,
      0, -4, 20, -4, 0, -1, -4, -1, 0, 0, 0, 0, 0, 0, -1, -4, -1, 0, -4, 20, -4,
      0, -1, -4, -1, 0, 0, 0, 0, 0, 0, -1, -4, 0, 0, -4, 20, 0, 0, -1, -4, 0, 0,
      0, 0, 0, 0, 0, 0, -4, -1, 0, 0, 20, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1,
      -4, -1, 0, -4, 20, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -4, -1, 0, -4,
      20, -4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -4, 0, 0, -4, 20;
  A_ref /= 6;
  ASSERT_TRUE(A.isApprox(A_ref));
}

TEST(mehrstellenverfahren, compMehrstellenf) {
  const unsigned int M = 4;
  const auto f = [](double x, double y) -> double {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  };
  Eigen::VectorXd phi = compMehrstellenf(f, M);
  Eigen::VectorXd phi_ref(M * M);
  phi_ref << 0.012940, 0.020937, 0.020937, 0.012940, 0.020937, 0.033877,
      0.033877, 0.020937, 0.020937, 0.033877, 0.033877, 0.020937, 0.012940,
      0.020937, 0.020937, 0.012940;
  const bool is_close = (phi - phi_ref).lpNorm<Eigen::Infinity>() < 1e-6;
  ASSERT_TRUE(is_close);
}

TEST(mehrstellenverfahren, solveMehrstellen) {
  const unsigned int M = 32;
  const auto f = [](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  };
  const auto u_exact = [](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) / (2 * M_PI * M_PI);
  };
  const Eigen::VectorXd mu = solveMehrstellen(f, M);
  Eigen::VectorXd mu_ref(M * M);
  const double h = 1. / (M + 1);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      const double x = (i + 1) * h;
      const double y = (j + 1) * h;
      const int k = i * M + j;
      mu_ref[k] = u_exact(x, y);
    }
  }
  const double diff_infty_norm = (mu - mu_ref).lpNorm<Eigen::Infinity>();
  ASSERT_LT(diff_infty_norm, 1e-6);
}

TEST(mehrstellenverfahren, compgriderr) {
  Eigen::VectorXi Ms(4);
  Ms << 4, 8, 16, 32;
  Eigen::VectorXd errs_ref(4);
  errs_ref << 0.000019333, 0.0000020112, 0.00000016239, 0.000000011526;
  for (int i = 0; i < Ms.size(); ++i) {
    const unsigned int M = Ms[i];
    const double err = compgriderr(M);
    const double err_ref = errs_ref[i];
    const double rel = err / err_ref;
    ASSERT_NEAR(rel, 1, 1e-4);
  }
}

}  // namespace mehrstellenverfahren::test
