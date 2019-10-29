#include "incidence_mat.h"

using namespace Eigen;
using lf::mesh::Mesh;

namespace IncidenceMatrices {

// @brief Create the mesh consisting of a triangle and quadrilateral
//        from the exercise sheet.
// @return Shared pointer to the hybrid2d mesh.
std::shared_ptr<Mesh> createDemoMesh() {
  // data type for a hybrid mesh in a world of dimension 2
  lf::mesh::hybrid2d::MeshFactory builder(2);

  // Add points
  builder.AddPoint(Vector2d{0, 0});    // (0)
  builder.AddPoint(Vector2d{1, 0});    // (1)
  builder.AddPoint(Vector2d{1, 1});    // (2)
  builder.AddPoint(Vector2d{0, 1});    // (3)
  builder.AddPoint(Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  builder.AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      {1, 2, 4},                 // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  builder.AddEntity(lf::base::RefEl::kQuad(), {0, 1, 4, 3},
                    std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<Mesh> demoMesh = builder.Build();

  return demoMesh;
}

// @brief Compute the edge-vertex incidence matrix G for a given mesh
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
SparseMatrix<int> computeEdgeVertexIncidenceMatrix(const Mesh& mesh) {
  // Store edge-vertex incidence matrix here
  SparseMatrix<int, RowMajor> G;

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here */
  /* END_SOLUTION */

  return G;
}

// @brief Compute the cell-edge incidence matrix D for a given mesh
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
SparseMatrix<int> computeCellEdgeIncidenceMatrix(const Mesh& mesh) {
  // Store cell-edge incidence matrix here
  SparseMatrix<int, RowMajor> D;

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here */
  /* END_SOLUTION */

  return D;
}

// @brief For a given mesh test if the product of cell-edge and edge-vertex
//        incidence matrix is zero: D*G == 0?
// @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
//             such as lf::mesh::hybrid2d::Mesh)
// @return true, if the product is zero and false otherwise
bool testZeroIncidenceMatrixProduct(const Mesh& mesh) {
  bool isZero = false;

  /* BEGIN_SOLUTION */
  /* TODO Your implementation goes here */
  /* END_SOLUTION */

  return isZero;
}

}  // namespace IncidenceMatrices
