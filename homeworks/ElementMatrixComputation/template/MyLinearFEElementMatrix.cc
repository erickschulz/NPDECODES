/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 03.03.2019
 * @copyright Developed at ETH Zurich
 */

#include "MyLinearFEElementMatrix.h"

namespace ElementMatrixComputation {

MyLinearFEElementMatrix::ElemMat MyLinearFEElementMatrix::Eval(
    const lf::mesh::Entity &cell) {
  // Topological type of the cell
  const lf::base::RefEl ref_el{cell.RefEl()};

  // Obtain the vertex coordinates of the cell, which completely
  // describe its shape.
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();

  // Matrix storing corner coordinates in its columns
  auto vertices = geo_ptr->Global(ref_el.NodeCoords());

  // Matrix for returning element matrix
  MyLinearFEElementMatrix::elem_mat_t elem_mat;

  /* BEGIN_SOLUTION */
  // Your implementation goes here!
  /* END_SOLUTION */

  return elem_mat;
}
}  // namespace ElementMatrixComputation