/**
 * @file qfeinterpolator.cc
 * @brief NPDE homework DebuggingFEM code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "qfeinterpolator.h"

#include <Eigen/Core>

#include <lf/mesh/mesh.h>

namespace DebuggingFEM {

Eigen::Vector2d globalCoordinate(int idx, const lf::mesh::Entity &cell) {
  auto geom = cell.Geometry();
  Eigen::Vector2d result;
  Eigen::MatrixXd corners(2, 3);
  corners << 0., 1., 0., 0., 0., 1.;
  //====================
  // Your code goes here
  //====================
  return result;
}

}  // namespace DebuggingFEM
