/**
 * @file
 * @brief NPDE homework ResidualErrorEstimator Tests
 * @author Ralf Hiptmair
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../residualerrorestimator.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <iostream>

namespace REE::test {

TEST(REE, TwoTriangleMesh) {
  /* Macros available in the Google test framework:
     EXPECT_EQ(x,y), EXPECT_NE(x,y), EXPECT_LT(x,y), EXPECT_LE(x,y),
     EXPECT_GT(x,y), EXPECT_GE(x,y) EXPECT_STREQ(x,y), EXPECT_STRNE(x,y) -> for
     C-strings only ! EXPECT_NEAR(x,y,abs_tol) All testing macros can output a
     message by a trailing << ....
   */
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;

  // Build a mesh with two triangles
  // Create helper object: mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Generate nodes of the mesh
  // clang-format off
  std::array<std::array<double, 2>, 4> node_coord{
  std::array<double, 2>({0 , 0 }),
  std::array<double, 2>({0.5 , 0 }),
  std::array<double, 2>({0.5 , 0.5 }),
  std::array<double, 2>({0 , 0.5 })};
  // clang-format on
  // Add nodes to the mesh via the MeshFactory object
  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(coord_t({node[0], node[1]}));
  }
  // Add plain triangles to the mesh, defined by their vertex nodes.
  // Since no particular geometry is specified, the triangles are assumed to
  // have straght edges.
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({0, 1, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({0, 2, 3}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Get a pointer to the coarsest mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();

  // Print information about the coarsest mesh
  std::cout << "\t Coarsest mesh for demonstration run\n";
  lf::mesh::utils::PrintInfo(std::cout, *mesh_p, 100);

  std::function<double(Eigen::Vector2d)> alpha =
      [](Eigen::Vector2d x) -> double {
    return ((x[0] - x[1]) > 0) ? 1.0 / 3.0 : 1.0 / 7.0;
  };

  auto f = [](Eigen::Vector2d x) -> double { return (1.0); };

  // Defining the discretized boundary value problem including the
  // finite-element space
  const dataDiscreteBVP disc_bvp(mesh_p, alpha, f);

  // Fix a finite element function for testing
  Eigen::VectorXd mu(4);
  mu(0) = 1;
  mu(1) = 0;
  mu(2) = 1;
  mu(3) = 0;

  // Compute edge contributions
  lf::mesh::utils::CodimMeshDataSet<double> ed_res{edgeResiduals(disc_bvp, mu)};

  // Values of the jump residual for the five edges of the mesh
  std::array<double, 5> vals({0.0, 2.721088435374150, 0.0, 0.0, 0.0});
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    EXPECT_NEAR(ed_res(*edge), vals[mesh_p->Index(*edge)], 1.0E-6);
  }
}

TEST(REE, ZeroEdgeResioduals) {
  // Obtain test mesh
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Defining the discretized boundary value problem including the
  // finite-element space
  std::function<double(Eigen::Vector2d)> alpha =
      [](Eigen::Vector2d x) -> double { return 1.0; };
  auto f = [](Eigen::Vector2d) -> double { return 1.0; };
  const dataDiscreteBVP disc_bvp(mesh_p, alpha, f);
  // A simple linear function
  auto lin_fun = [](Eigen::Vector2d x) -> double {
    return 2.0 * x[0] + x[1] + 1.0;
  };
  // Sample nodal values
  lf::mesh::utils::MeshFunctionGlobal mf_lin(lin_fun);
  Eigen::VectorXd mu =
      lf::fe::NodalProjection(*disc_bvp.pwlinfespace_p_, mf_lin);
  // Compute edge contributions
  lf::mesh::utils::CodimMeshDataSet<double> ed_res{edgeResiduals(disc_bvp, mu)};
  // All edge residual should vanish
  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    EXPECT_NEAR(ed_res(*edge), 0.0, 1.0E-6);
  }
}

}  // namespace REE::test
