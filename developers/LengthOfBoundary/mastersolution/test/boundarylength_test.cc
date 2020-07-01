/**
 * @ file master_tests.cc
 * @ brief NPDE homework LengthOfBoundary code
 * @ author Christian Mitsch
 * @ date 03.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "../boundarylength.h"

#include <gtest/gtest.h>
#include <lf/mesh/test_utils/test_meshes.h>

namespace LengthOfBoundary::test {

TEST(BoundaryLength, area_test) {
  // define the test mesh
  // "auto" = std::shared_ptr<lf::mesh::Mesh>
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_NEAR(LengthOfBoundary::volumeOfDomain(mesh_p), 9.0, 1e-12);
}

TEST(BoundaryLength, length_test) {
  // define the test mesh
  // "auto" = std::shared_ptr<lf::mesh::Mesh>
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh();

  EXPECT_NEAR(LengthOfBoundary::lengthOfBoundary(mesh_p), 12.0, 1e-12);
}

}  // namespace LengthOfBoundary::test
