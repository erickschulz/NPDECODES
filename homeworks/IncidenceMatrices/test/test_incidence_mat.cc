#include <gtest/gtest.h>
#include "../mysolution/incidence_mat.h"  

// Create demo mesh from exercise sheet and test if edge vertex incidence
// matrix is the same as the one calculated by hand
TEST(Homework_2_6, EdgeVertexIncidenceMatrix) {
  std::shared_ptr<Mesh> demoMesh = IncidenceMatrices::createDemoMesh();
  SparseMatrix<int> G_sp =
      IncidenceMatrices::computeEdgeVertexIncidenceMatrix(*demoMesh);
  MatrixXi G(G_sp);

  MatrixXi G_expected(6, 5);
  // clang-format off
  G_expected << 1, -1,  0,  0,  0, 
               -1,  0,  0,  1,  0, 
                0,  1, -1,  0,  0, 
                0, -1,  0,  0,  1,
                0,  0,  1,  0, -1, 
                0,  0,  0, -1,  1;
  // clang-format on

  EXPECT_EQ(G, G_expected);
}

// Create demo mesh from exercise sheet and test if cell edge incidence
// matrix is the same as the one calculated by hand
TEST(Homework_2_6, CellEdgeIncidenceMatrix) {
  std::shared_ptr<Mesh> demoMesh = IncidenceMatrices::createDemoMesh();
  SparseMatrix<int> D_sp =
      IncidenceMatrices::computeCellEdgeIncidenceMatrix(*demoMesh);
  MatrixXi D(D_sp);

  MatrixXi D_expected(2, 6);
  // clang-format off
  D_expected << 0, 0, 1,  1, 1, 0, 
                1, 1, 0, -1, 0, 1;
  // clang-format on

  EXPECT_EQ(D, D_expected);
}

// Test co-chain complex property (D*G = 0) at the example of the mesh
// in the exercise sheet
TEST(Homework_2_6, CoChainComplexProperty) {
  std::shared_ptr<Mesh> demoMesh = IncidenceMatrices::createDemoMesh();
  EXPECT_TRUE(IncidenceMatrices::testZeroIncidenceMatrixProduct(*demoMesh));
}
