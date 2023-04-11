/**
 * @file MinimalGraphSurface_test.cc
 * @brief NPDE homework MinimalGraphSurface code
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../minimalgraphsurface.h"

#include <gtest/gtest.h>
#include <lf/fe/fe_tools.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/mesh_function_global.h>

#include <Eigen/Core>

namespace MinimalGraphSurface::test {

TEST(MinimalGraphSurface, computeGraphArea) {
  // A triangular mesh of the unit square
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // A globally linear FE function
  lf::mesh::utils::MeshFunctionGlobal mf_lin(
      [](Eigen::Vector2d x) -> double { return x[0] + 2 * x[1]; });
  Eigen::VectorXd lin_vec = lf::fe::NodalProjection(*fe_space_p, mf_lin);
  // Compute area of graph
  double area = computeGraphArea(fe_space_p, lin_vec);
  EXPECT_NEAR(area, std::sqrt(6), 1.0E-10);
}

}  // namespace MinimalGraphSurface::test

TEST(MinimalGraphSurface, coeffTensorA) {}
TEST(MinimalGraphSurface, coeffScalarc) {}
TEST(MinimalGraphSurface, computeNewtonCorrection) {}
