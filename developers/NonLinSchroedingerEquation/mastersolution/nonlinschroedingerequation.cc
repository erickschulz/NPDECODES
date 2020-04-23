/**
 * @file nonlinschroedingerequation.cc
 * @brief NPDE homework NonLinSchroedingerEquation code
 * @author Oliver Rietmann
 * @date 22.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "nonlinschroedingerequation.h"

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>

namespace NonLinSchroedingerEquation {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix3d MassElementMatrixProvider::Eval(const lf::mesh::Entity &cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(), "Unsupported cell type " << cell.RefEl());
  Eigen::Matrix3d element_matrix;
#if SOLUTION
  double area = lf::geometry::Volume(*(cell.Geometry()));
  element_matrix = area / 3.0 * Eigen::Matrix3d::Identity();
#else
  //====================
  // Your code goes here
  //====================
#endif
  return element_matrix;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d StiffnessElementMatrixProvider::Eval(const lf::mesh::Entity &cell) {
  LF_VERIFY_MSG(cell.RefEl() == lf::base::RefEl::kTria(), "Unsupported cell type " << cell.RefEl());
  Eigen::Matrix3d element_matrix;
#if SOLUTION
  const lf::geometry::Geometry *geo_ptr = cell.Geometry();
  Eigen::Matrix<double, 2, 3> a = lf::geometry::Corners(*geo_ptr);
  double area = lf::geometry::Volume(*geo_ptr);
  Eigen::Matrix<double, 2, 3> g;
  g << a(1, 1) - a(1, 2), a(1, 2) - a(1, 0), a(1, 0) - a(1, 1),
       a(0, 2) - a(0, 1), a(0, 0) - a(0, 2), a(0, 1) - a(0, 0);
  g *= 1.0 / (2.0 * area);
  element_matrix = area * g.transpose() * g;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return element_matrix;
}
/* SAM_LISTING_END_2 */

/*
  Eigen::Matrix3d S33;
  S33 << 0.0, 1.0, -1.0, -1.0, 0.0, 1.0, 1.0, -1.0, 0.0;
  Eigen::Matrix2d S22;
  S22 << 0.0, -1.0, 1.0, 0.0;
  Eigen::Matrix<double, 2, 3> g = 1.0 / (2.0 * area) * S22 * a * S33;
*/

}  // namespace NonLinSchroedingerEquation
