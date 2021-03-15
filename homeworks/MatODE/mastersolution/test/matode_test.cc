/**
 * @file matode_test.cc
 * @brief NPDE homework MatODE code
 * @author Oliver Rietmann
 * @date 14.03.2021
 * @copyright Developed at ETH Zurich
 */

#include "../matode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <vector>

namespace MatODE::test {


TEST(MatODE, eeulstep) {
  // Tools
  Eigen::VectorXd h(5);
  h << 0.05, 0.1, 0.3, 0.5, 1.0;
  Eigen::Matrix3d A;
  A << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0;
  const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();

  std::vector<Eigen::Matrix3d> result(5);
  for (int i = 0; i < 5; ++i){
  result[i] = MatODE::eeulstep(A, Id, h[i]);
}

  std::vector<Eigen::Matrix3d> reference(5);
  reference[0] << 1.0, 0.05, 0.0, 0.05, 1.0, 0.05, 0.05, 0.05, 1.0;
  reference[1] << 1.0, 0.1, 0.0, 0.1, 1.0, 0.1, 0.1, 0.1, 1.0;
  reference[2] << 1.0, 0.3, 0.0, 0.3, 1.0, 0.3, 0.3, 0.3, 1.0;
  reference[3] << 1.0, 0.5, 0.0, 0.5, 1.0, 0.5, 0.5, 0.5, 1.0;
  reference[4] << 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;

  Eigen::VectorXd error(5);
  for (int i = 0; i < 5; ++i) {
    error[i] = (reference[i] - result[i]).lpNorm<Eigen::Infinity>();
  }

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

TEST(MatODE, ieulstep) {
  // Tools
  Eigen::VectorXd h(5);
  h << 0.05, 0.1, 0.3, 0.5, 1.0;
  Eigen::Matrix3d A;
  A << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0;
  const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();

  std::vector<Eigen::Matrix3d> result(5);
  for (int i = 0; i < 5; ++i){
  result[i] = MatODE::ieulstep(A, Id, h[i]);
}

  std::vector<Eigen::Matrix3d> reference(5);
  reference[0] << 1.00263852242744, 0.0502575700464882, 0.00251287850232441,
      0.0527704485488127, 1.00515140092976, 0.0502575700464882,
      0.0527704485488127, 0.0527704485488127, 1.00263852242744;
  reference[1] << 1.01123595505618, 0.102145045965271, 0.0102145045965271,
      0.112359550561798, 1.02145045965271, 0.102145045965271, 0.112359550561798,
      0.112359550561798, 1.01123595505618;
  reference[2] << 1.14754098360656, 0.378310214375788, 0.113493064312736,
      0.491803278688525, 1.26103404791929, 0.378310214375788, 0.491803278688525,
      0.491803278688525, 1.14754098360656;
  reference[3] << 2.0, 1.33333333333333, 0.666666666666667, 2.0,
      2.66666666666667, 1.33333333333333, 2.0, 2.0, 2.0;
  reference[4] << 0.0, -0.5, -0.5, -1.0, -0.5, -0.5, -1.0, -1.0, -0.0;

   Eigen::VectorXd error(5);
  for (int i = 0; i < 5; ++i) {
    error[i] = (reference[i] - result[i]).lpNorm<Eigen::Infinity>();
  }

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

TEST(MatODE, impstep) {
  // Tools
  Eigen::VectorXd h(5);
  h << 0.05, 0.1, 0.3, 0.5, 1.0;
  Eigen::Matrix3d A;
  A << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0;
  const Eigen::Matrix3d Id = Eigen::Matrix3d::Identity();

  std::vector<Eigen::Matrix3d> result(5);
  for (int i = 0; i < 5; ++i){
  result[i] = MatODE::impstep(A, Id, h[i]);
}

  std::vector<Eigen::Matrix3d> reference(5);
  reference[0] << 1.00128287363695, 0.0500633614418248, 0.00125158403604562,
      0.0513149454778704, 1.00253445767299, 0.0500633614418248,
      0.0513149454778704, 0.0513149454778704, 1.00128287363695;
  reference[1] << 1.00527704485488, 0.100515140092977, 0.00502575700464883,
      0.105540897097625, 1.01030280185953, 0.100515140092976, 0.105540897097625,
      0.105540897097625, 1.00527704485488;
  reference[2] << 1.05438066465257, 0.315250229869959, 0.0472875344804939,
      0.362537764350453, 1.10166819913306, 0.315250229869959, 0.362537764350453,
      0.362537764350453, 1.05438066465257;
  reference[3] << 1.18181818181818, 0.581818181818182, 0.145454545454545,
      0.727272727272727, 1.32727272727273, 0.581818181818182, 0.727272727272727,
      0.727272727272727, 1.18181818181818;
  reference[4] << 3.0, 2.66666666666667, 1.33333333333333, 4.0,
      4.33333333333333, 2.66666666666667, 4.0, 4.0, 3.0;

   Eigen::VectorXd error(5);
  for (int i = 0; i < 5; ++i) {
    error[i] = (reference[i] - result[i]).lpNorm<Eigen::Infinity>();
  }

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, error.lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace MatODE::test
