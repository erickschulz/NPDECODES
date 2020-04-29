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

/* SAM_LISTING_BEGIN_1 */
Eigen::Vector2d globalCoordinate(int idx, const lf::mesh::Entity &cell) {
  // Consistency check for arguments
  LF_ASSERT_MSG(cell.RefEl() == lf::base::RefEl::kTria(),
                "Implemented for triangles only");
  // Fetch pointer to asscoiated geometry object
  lf::geometry::Geometry *geom = cell.Geometry();
  // For returning the global coordinates of the interpolation node
  Eigen::Vector2d result;
  // Reference coordinates of the vertices of the triangle
  Eigen::Matrix<double, 2, 3> corners(2, 3);
  corners << 0., 1., 0., 0., 0., 1.;
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

} // namespace DebuggingFEM
