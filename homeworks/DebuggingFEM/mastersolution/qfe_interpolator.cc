/**
 * @file qfe_interpolator.cc
 * @brief NPDE homework DebuggingFEM code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "qfe_interpolator.h"

namespace DebuggingFEM {

Eigen::Vector2d globalCoordinate(int idx, const lf::mesh::Entity &cell) {
  auto geom = cell.Geometry();
  Eigen::Vector2d result;
  Eigen::MatrixXd corners(2, 3);
  corners << 0., 1., 0., 0., 0., 1.;
  /* SOLUTION_BEGIN */
  switch (idx) {
    case (0):
      result = geom->Global(corners.col(0));
      break;
    case (1):
      result = geom->Global(corners.col(1));
      break;
    case (2):
      result = geom->Global(corners.col(2));
      break;
    case (3):
      result = geom->Global((corners.col(0) + corners.col(1)) / 2.);
      break;
    case (4):
      result = geom->Global((corners.col(1) + corners.col(2)) / 2.);
      break;
    case (5):
      result = geom->Global((corners.col(2) + corners.col(0)) / 2.);
      break;
    default:
      throw std::invalid_argument("idx needs to be in range [0,5]");
  }
  /* SOLUTION_END */
  return result;
}

}  // namespace DebuggingFEM
