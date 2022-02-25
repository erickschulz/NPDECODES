/**
 * @file
 * @brief NPDE homework TEMPLATE MAIN FILE
 * @author
 * @date
 * @copyright Developed at SAM, ETH Zurich
 */

#include "../nitschemethod.h"

#include <gtest/gtest.h>

#include <cmath>

namespace NitscheMethod::test {

TEST(NitscheMethod, EvalTest) {
  /* Macros available in the Google test framework:
     EXPECT_EQ(x,y), EXPECT_NE(x,y), EXPECT_LT(x,y), EXPECT_LE(x,y),
     EXPECT_GT(x,y), EXPECT_GE(x,y) EXPECT_STREQ(x,y), EXPECT_STRNE(x,y) -> for
     C-strings only ! EXPECT_NEAR(x,y,abs_tol) All testing macros can output a
     message by a trailing << ....
   */
  using size_type = lf::base::size_type;
  using glb_idx_t = lf::assemble::glb_idx_t;
  using coord_t = Eigen::Vector2d;
  // Build a "mesh" containing the reference triangle only
  // Create helper object: mesh factory
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
  // Set coordinates of the nodes of the mesh
  // clang-format off
  std::array<std::array<double, 2>, 3> node_coord{
  std::array<double, 2>({0 , 0 }),
  std::array<double, 2>({1 , 0 }),
  std::array<double, 2>({0 , 1})};
  // clang-format on
  // Add nodes to the mesh via the MeshFactory object
  for (const auto& node : node_coord) {
    mesh_factory_ptr->AddPoint(coord_t({node[0], node[1]}));
  }
  // Define the single triangle of the mesh
  mesh_factory_ptr->AddEntity(lf::base::RefEl::kTria(),
                              std::vector<size_type>({0, 1, 2}),
                              std::unique_ptr<lf::geometry::Geometry>(nullptr));
  // Generate mesh and obtain a pointer to it
  std::shared_ptr<lf::mesh::Mesh> mesh_p = mesh_factory_ptr->Build();
  // Print information about the mesh
  std::cout << "\t Coarsest mesh for demonstration run\n";
  lf::mesh::utils::PrintInfo(std::cout, *mesh_p, 100);
  /*
    Node 0: [ 0 0 ],  Node 1: [ 1 0 ],  Node 2: [ 0 1 ]
    Edge 0: 0 -> 1, Edge 1: 2 -> 0, Edge 2: 1 -> 2
    TRIA: Nodes = 0 1 2, Edges = 0 2 1
   */
  // Set edge #2 as a boundary edge
  const int single_bd_ed_idx = 2;
  lf::mesh::utils::CodimMeshDataSet<bool> edge_marker(mesh_p, 1, false);
  edge_marker(*mesh_p->EntityByIndex(1, single_bd_ed_idx)) = true;

  // Compute element matrix for the single triangle of the mesh
  const double c = std::sqrt(2) / 2.0;
  const lf::mesh::Entity& K = *mesh_p->EntityByIndex(0, 0);
  NitscheMethod::LinearFENitscheElementMatrix elmat_builder(edge_marker, c);
  const Eigen::Matrix3d elmat = elmat_builder.Eval(K);
  std::cout << "Element matrix = " << std::endl << elmat << std::endl;

  // Reference result
  const double sq2 = std::sqrt(2);
  Eigen::Matrix3d A =
      (Eigen::Matrix3d() << 1, 0.5, 0.5, 0.5, -0.5 + c / 3.0 * sq2,
       -1.0 + c / 6.0 * sq2, 0.5, -1.0 + c / 6.0 * sq2, -0.5 + c / 3.0 * sq2)
          .finished();
  EXPECT_NEAR((elmat - A).norm(), 0.0, 1E-10);
}
}  // namespace NitscheMethod::test
