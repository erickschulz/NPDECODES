/**
 * @file newproblem_test_mastersolution.cc
 * @brief NPDE homework NewProblem code
 * @author Oliver Rietmann, Erick Schulz
 * @date 01.01.2020
 * @copyright Developed at ETH Zurich
 */

// HACK:
#undef SOLUTION
#define SOLUTION 1

#include <gtest/gtest.h>

#if SOLUTION
#include "../mastersolution/electrostaticforce.h"
#else
#include "../mysolution/electrostaticforce.h"
#endif

#include <Eigen/Core>

namespace ElectrostaticForce::test {

TEST(ElectrostaticForce, test1) {

  double tol = 1.0e-8;
  ASSERT_NEAR(0.0, 0.0, tol);
}

}  // namespace ElectrostaticForce::test
