/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "MyLinearLoadVector.h"

namespace ElementMatrixComputation {

MyLinearLoadVector::ElemVec MyLinearLoadVector::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};
  const lf::base::size_type num_nodes{ref_el.NumNodes()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  return computeLoadVector(vertices, f_);
}

MyLinearLoadVector::ElemVec computeLoadVector(
    Eigen::MatrixXd vertices, const MyLinearLoadVector::FHandle_t f) {
  // Number of nodes of the element: triangles = 3, rectangles = 4
  const int num_nodes = vertices.cols();

  // Vector for returning element vector
  MyLinearLoadVector::elem_vec_t elem_vec =
      MyLinearLoadVector::elem_vec_t::Zero();

  /* BEGIN_SOLUTION */
  // Your implementation goes here!
  /* END_SOLUTION */

  return elem_vec;
}

}  // namespace ElementMatrixComputation
