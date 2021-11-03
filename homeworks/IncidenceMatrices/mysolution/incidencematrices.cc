#include "incidencematrices.h"

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <array>
#include <memory>

namespace IncidenceMatrices {

/** @brief Create the mesh consisting of a triangle and quadrilateral
 *         from the exercise sheet.
 * @return Shared pointer to the hybrid2d mesh.
 */
std::shared_ptr<lf::mesh::Mesh> createDemoMesh() {
  // builder for a hybrid mesh in a world of dimension 2
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
      std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);

  // Add points
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 0});    // (0)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 0});    // (1)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{1, 1});    // (2)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0, 1});    // (3)
  mesh_factory_ptr->AddPoint(Eigen::Vector2d{0.5, 1});  // (4)

  // Add the triangle
  // First set the coordinates of its nodes:
  Eigen::MatrixXd nodesOfTria(2, 3);
  nodesOfTria << 1, 1, 0.5, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kTria(),  // we want a triangle
      std::array<lf::mesh::Mesh::size_type, 3>{
          {1, 2, 4}},  // indices of the nodes
      std::make_unique<lf::geometry::TriaO1>(nodesOfTria));  // node coords

  // Add the quadrilateral
  Eigen::MatrixXd nodesOfQuad(2, 4);
  nodesOfQuad << 0, 1, 0.5, 0, 0, 0, 1, 1;
  mesh_factory_ptr->AddEntity(
      lf::base::RefEl::kQuad(),
      std::array<lf::mesh::Mesh::size_type, 4>{{0, 1, 4, 3}},
      std::make_unique<lf::geometry::QuadO1>(nodesOfQuad));

  std::shared_ptr<lf::mesh::Mesh> demoMesh_p = mesh_factory_ptr->Build();

  return demoMesh_p;
}

/** @brief Compute the edge-vertex incidence matrix G for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The edge-vertex incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<int> computeEdgeVertexIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store edge-vertex incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> G;

  //====================
  // Your code goes here
  //====================

  return G;
}
/* SAM_LISTING_END_1 */

/** @brief Compute the cell-edge incidence matrix D for a given mesh
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *        such as lf::mesh::hybrid2d::Mesh)
 * @return The cell-edge incidence matrix as Eigen::SparseMatrix<int>
 */
/* SAM_LISTING_BEGIN_2 */
Eigen::SparseMatrix<int> computeCellEdgeIncidenceMatrix(
    const lf::mesh::Mesh &mesh) {
  // Store cell-edge incidence matrix here
  Eigen::SparseMatrix<int, Eigen::RowMajor> D;

  //====================
  // Your code goes here
  //====================

  return D;
}
/* SAM_LISTING_END_2 */

/** @brief For a given mesh test if the product of cell-edge and edge-vertex
 *        incidence matrix is zero: D*G == 0?
 * @param mesh The input mesh of type lf::mesh::Mesh (or of derived type,
 *             such as lf::mesh::hybrid2d::Mesh)
 * @return true, if the product is zero and false otherwise
 */
/* SAM_LISTING_BEGIN_3 */
bool testZeroIncidenceMatrixProduct(const lf::mesh::Mesh &mesh) {
  bool isZero = false;

  //====================
  // Your code goes here
  //====================
  return isZero;
}
/* SAM_LISTING_END_3 */

}  // namespace IncidenceMatrices
