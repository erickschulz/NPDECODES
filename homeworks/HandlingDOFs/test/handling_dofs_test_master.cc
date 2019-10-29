#include "../mastersolution/handling_dofs.h"
#include <gtest/gtest.h>
#include <string>

/** @brief Generate a mesh simply consisting of a single triangle
 *  (0,1)
 *    |   \
 *    |    \
 *  (0,0)--(1,0)
 * @return A shared ptr to const lf::mesh::Mesh with above triangle
 */
std::shared_ptr<const lf::mesh::Mesh> singleTriagMesh() {
  lf::mesh::hybrid2d::MeshFactory builder(2);

  // Add points
  builder.AddPoint(Eigen::Vector2d{0, 0});  // (0)
  builder.AddPoint(Eigen::Vector2d{1, 0});  // (1)
  builder.AddPoint(Eigen::Vector2d{0, 1});  // (2)

  // Add the triangle
  // First set the coordinates of its nodes:
  Eigen::MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 0, 1, 0, 0, 0, 1;
  builder.AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      {0, 1, 2},                 // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  std::shared_ptr<const lf::mesh::Mesh> mesh = builder.Build();

  return mesh;
}

TEST(Homework_2_9, SimpleNumEntityDofs) {
  // Create a triangular mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh = singleTriagMesh();

  // Create a dofhandler for the linear lagrangian FE space
  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  std::array<std::size_t, 3> entityDofs =
      HandlingDOFs::countEntityDofs(lin_dofh);
  EXPECT_EQ(entityDofs[0], 0);  // 0 dofs in cells
  EXPECT_EQ(entityDofs[1], 0);  // 0 dofs on edges
  EXPECT_EQ(entityDofs[2], 3);  // 3 dofs on nodes

  // Create a dofhandler for the quadratic lagrangian FE space
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 1},  // 1 dof on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  entityDofs = HandlingDOFs::countEntityDofs(quad_dofh);
  EXPECT_EQ(entityDofs[0], 0);  // 0 dofs in cells
  EXPECT_EQ(entityDofs[1], 3);  // 3 dofs on edges
  EXPECT_EQ(entityDofs[2], 3);  // 3 dofs on nodes
}

TEST(Homework_2_9, SimpleNumBoundaryDofs) {
  // Create a triangular mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh = singleTriagMesh();

  // Create a dofhandler for the linear lagrangian FE space
  // For the S_1^0 space we have 1 dof per node, and 0 for all other entities:
  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof per node
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // 3 nodes on boundary
  EXPECT_EQ(HandlingDOFs::countBoundaryDofs(lin_dofh), 3);

  // Create a dofhandler for the quadratic lagrangian FE space
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof per node
             {lf::base::RefEl::kSegment(), 1},  // 1 dof per edge
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // 3 nodes and 3 edges on boundary
  EXPECT_EQ(HandlingDOFs::countBoundaryDofs(quad_dofh), 6);
}

TEST(Homework_2_9, NumEntityDofs) {
  // Create a triangular mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Create a dofhandler for the linear lagrangian FE space
  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  std::array<std::size_t, 3> entityDofs =
      HandlingDOFs::countEntityDofs(lin_dofh);
  EXPECT_EQ(entityDofs[0], 0);   // 0 dofs in cells
  EXPECT_EQ(entityDofs[1], 0);   // 0 dofs on edges
  EXPECT_EQ(entityDofs[2], 13);  // 13 dofs on nodes
}

TEST(Homework_2_9, NumBoundaryDofs) {
  // Create a triangular mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Create a dofhandler for the linear lagrangian FE space
  // For the S_1^0 space we have 1 dof per node, and 0 for all other entities:
  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // 9 nodes on boundary
  EXPECT_EQ(HandlingDOFs::countBoundaryDofs(lin_dofh), 9); 
}

TEST(Homework_2_9, SimpleIntegration) {
  // Generate a mesh simply consisting of a single triangle
  //  (0,1)
  //    |   \
  //    |    \
  //  (0,0)--(1,0)
  // And integrate f(x) = x_1*x_2 over it. For linear integration,
  // where only the values on the nodes are taken into account we expect
  // to get 0, since f(nodes) = 0.
  // For quadratic integration only the edge between (1,0) and (0,1) has
  // a non-zero contribution, at the midpoint of that edge the value of
  // f is f([0.5 0.5]) = 0.25, hence I = 0.25*|K|/3 = 0.125/3
  // Create a dofhandler for the linear lagrangian FE space
  // For the S_1^0 space we have 1 dof per node, and 0 for all other entities:
  std::shared_ptr<const lf::mesh::Mesh> mesh = singleTriagMesh();

  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 1},  // 1 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // integrating a linear function
  auto f = [](const Eigen::Vector2d& x) { return x[0] * x[1]; };
  Eigen::VectorXd mu = HandlingDOFs::buildCoefVector(f, lin_dofh);
  Eigen::VectorXd zeta = HandlingDOFs::buildCoefVector(f, quad_dofh);

  EXPECT_EQ(HandlingDOFs::integrateLinearFEFunction(lin_dofh, mu), 0);
  const double precision = 1e-15;
  EXPECT_NEAR(HandlingDOFs::integrateQuadraticFEFunction(quad_dofh, zeta),
              0.125 / 3, precision);
}

TEST(Homework_2_9, Integration) {
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 1},  // 1 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // integrating a linear function
  auto f = [](const Eigen::Vector2d& x) { return x[0] * x[1]; };
  Eigen::VectorXd mu = HandlingDOFs::buildCoefVector(f, lin_dofh);
  Eigen::VectorXd zeta = HandlingDOFs::buildCoefVector(f, quad_dofh);

  const double precision = 1e-15;
  EXPECT_NEAR(HandlingDOFs::integrateLinearFEFunction(lin_dofh, mu), 
              487.0/24.0, precision);
  EXPECT_NEAR(HandlingDOFs::integrateQuadraticFEFunction(quad_dofh, zeta),
              20.25, precision);
}

TEST(Homework_2_9, ConstructingZetaFromMu) {
  std::shared_ptr<const lf::mesh::Mesh> mesh =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);

  lf::assemble::UniformFEDofHandler lin_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 0},  // 0 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals
  lf::assemble::UniformFEDofHandler quad_dofh(
      mesh, {{lf::base::RefEl::kPoint(), 1},    // 1 dof on nodes
             {lf::base::RefEl::kSegment(), 1},  // 1 dofs on edges
             {lf::base::RefEl::kTria(), 0},     // 0 dofs in triangles
             {lf::base::RefEl::kQuad(), 0}});   // 0 dofs in quadrilaterals

  // integrating a linear function
  auto f = [](const Eigen::Vector2d& x) { return x[0] * x[1]; };
  Eigen::VectorXd mu = HandlingDOFs::buildCoefVector(f, lin_dofh);
  Eigen::VectorXd zeta = HandlingDOFs::buildCoefVector(f, quad_dofh);
  Eigen::VectorXd zeta_constructed = HandlingDOFs::convertDOFsLinearQuadratic(lin_dofh, quad_dofh, mu);

  double rel_difference = (zeta - zeta_constructed).norm()/zeta.norm();
  double tol = 0.05;

  EXPECT_NEAR(rel_difference, 0, tol);
}
