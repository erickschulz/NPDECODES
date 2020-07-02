/**
 * @file radauthreetimesteppingode_test.cc
 * @brief NPDE homework "RadauThreeTimestepping" code
 * @author Tobias Rohner
 * @date 16.03.2020
 * @copyright Developed at ETH Zurich
 */

#include <vector>

#include <gtest/gtest.h>

#include "../radauthreetimesteppingode.h"

namespace RadauThreeTimestepping::test {

TEST(RadauThreeTimestepping, twoStageRadauTimeSteppingLinScalODE) {
  // Solve the ODE using a single step and check whether the factor is correct
  const std::vector<double> sol =
      RadauThreeTimestepping::twoStageRadauTimesteppingLinScalODE(1);
  ASSERT_TRUE(sol.size() == 2)
      << "Have you forgotten to add the initial state to the vector?";
  ASSERT_NEAR(sol[0], 1.0, 1e-10) << "Initial value should be 1.0";
  ASSERT_NEAR(sol[1], -4. / 51, 1e-10);
}

}  // end namespace RadauThreeTimestepping::test
