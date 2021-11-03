/**
 * @file ipdgfem_test_mastersolution.cc
 * @brief NPDE homework IPDGFEM code
 * @author Philippe Peter
 * @date 22.11.2019
 * @copyright Developed at ETH Zurich
 */

// HACK:
#undef SOLUTION
#define SOLUTION 1

#include "../ipdgfem.h"

#include <gtest/gtest.h>

#include <Eigen/Core>

namespace IPDGFEM::test {

TEST(NewProblem, dummyFunction) {
  double x = 0.0;
  int n = 0;

  Eigen::Vector2d v = IPDGFEM::dummyFunction(x, n);

  Eigen::Vector2d v_ref = {1.0, 1.0};

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, (v - v_ref).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace IPDGFEM::test
