#include <gtest/gtest.h>
#include "../SimpleLinearFEM2D.h"
#include "../local_assembler.h"

TEST(SimpleLinearFEM2DTest, TestMatrixAssembler) {
  TriaMesh2D mesh("test/meshes/test_msh.txt");
  MatrixAssembler mat_assembler(ElementMatrix_Lapl_LFE);
  Eigen::MatrixXd A = mat_assembler.Assemble(mesh);
  EXPECT_NEAR(A.sum(), 0., 0.0001);
  EXPECT_NEAR(A.maxCoeff(), 1., 0.0001);
}

TEST(SimpleLinearFEM2DTest, TestSlowMatrixAssembler) {
  TriaMesh2D mesh("test/meshes/test_msh.txt");
  SlowMatrixAssembler mat_assembler(ElementMatrix_Lapl_LFE);
  auto A = mat_assembler.Assemble(mesh);
  EXPECT_NEAR(A.sum(), 0., 0.0001);
}

TEST(SimpleLinearFEM2DTest, TestVectorAssembler) {
  TriaMesh2D mesh("test/meshes/test_msh.txt");
  auto f = [](const Eigen::Vector2d &x) { return 1.; };
  VectorAssembler vec_assembler(localLoadLFE, f);
  auto phi = vec_assembler.Assemble(mesh);
  EXPECT_NEAR(phi.sum(), 0.5, 0.0001);
}

TEST(SimpleLinearFEM2DTest, TestSolver) {
  TriaMesh2D mesh("test/meshes/test_msh.txt");
  auto f = [](const Eigen::Vector2d &x) { return 1.; };
  FESolver solver(f);
  auto mu = solver.Solve(mesh);
  EXPECT_NEAR(mu.sum(), 3, 0.0001);
}
