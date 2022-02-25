/**
 * @file
 * @brief NPDE homework ConvBLFMatrixProvider
 * @author Ralf Hiptmair
 * @date May 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../convblfmatrixprovider.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <iostream>

namespace cblfdemo::test {

TEST(ConvBLFMatrixProvider, Test) {
  /* Macros available in the Google test framework:
     EXPECT_EQ(x,y), EXPECT_NE(x,y), EXPECT_LT(x,y), EXPECT_LE(x,y),
     EXPECT_GT(x,y), EXPECT_GE(x,y) EXPECT_STREQ(x,y), EXPECT_STRNE(x,y) -> for
     C-strings only ! EXPECT_NEAR(x,y,abs_tol) All testing macros can output a
     message by a trailing << ....
   */
  // Obtain a purely triangular mesh from the collection of LehrFEM++'s
  // built-in meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p{
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3)};
  // vectors for testing
  const Eigen::Vector2d a(1.0, 2.0);
  const Eigen::Vector2d b(3.0, 2.0);
  // Run test computation
  const double itg = cblfdemo::testCDBLF(mesh_p, a, b);
  EXPECT_NEAR(itg, 63, 1E-6);
}
}  // namespace cblfdemo::test
