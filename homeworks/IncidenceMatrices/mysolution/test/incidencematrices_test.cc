#include "../incidencematrices.h"

#include <gtest/gtest.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <memory>

namespace IncidenceMatrices::test {

// Create demo mesh from exercise sheet and test if edge vertex incidence
// matrix is the same as the one calculated by hand
TEST(Homework_2_6, EdgeVertexIncidenceMatrix) {
  std::shared_ptr<lf::mesh::Mesh> demoMesh =
      IncidenceMatrices::createDemoMesh();
  Eigen::SparseMatrix<int> G_sp =
      IncidenceMatrices::computeEdgeVertexIncidenceMatrix(*demoMesh);
  Eigen::MatrixXi G(G_sp);

  Eigen::MatrixXi G_expected(6, 5);
  // clang-format off
  G_expected << 1, -1,  0,  0,  0, 
               -1,  0,  0,  1,  0, 
                0,  1, -1,  0,  0, 
                0, -1,  0,  0,  1,
                0,  0,  1,  0, -1, 
                0,  0,  0, -1,  1;
  // clang-format on
  EXPECT_EQ(G.rows(), G_expected.rows());
  EXPECT_EQ(G.cols(), G_expected.cols());
  if (G.rows() == G_expected.rows() && G.cols() == G_expected.cols())
    EXPECT_EQ(G, G_expected);
}

// Create demo mesh from exercise sheet and test if cell edge incidence
// matrix is the same as the one calculated by hand
TEST(Homework_2_6, CellEdgeIncidenceMatrix) {
  std::shared_ptr<lf::mesh::Mesh> demoMesh =
      IncidenceMatrices::createDemoMesh();
  Eigen::SparseMatrix<int> D_sp =
      IncidenceMatrices::computeCellEdgeIncidenceMatrix(*demoMesh);
  Eigen::MatrixXi D(D_sp);

  Eigen::MatrixXi D_expected(2, 6);
  // clang-format off
  D_expected << 0, 0, 1,  1, 1, 0,
                1, 1, 0, -1, 0, 1;
  // clang-format on
  EXPECT_EQ(D.rows(), D_expected.rows());
  EXPECT_EQ(D.cols(), D_expected.cols());
  if (D.rows() == D_expected.rows() && D.cols() == D_expected.cols())
    EXPECT_EQ(D, D_expected);
}

// Test co-chain complex property (D*G = 0) at the example of the mesh
// in the exercise sheet
TEST(Homework_2_6, CoChainComplexProperty) {
  std::shared_ptr<lf::mesh::Mesh> demoMesh =
      IncidenceMatrices::createDemoMesh();
  EXPECT_TRUE(IncidenceMatrices::testZeroIncidenceMatrixProduct(*demoMesh));
}

}  // namespace IncidenceMatrices::test
