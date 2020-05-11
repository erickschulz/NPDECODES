/**
 * @file qfeinterpolator.cc
 * @brief NPDE homework DebuggingFEM code
 * @author Simon Meierhans
 * @date 27/03/2019
 * @copyright Developed at ETH Zurich
 */

#include "qfeinterpolator.h"

#include <lf/mesh/mesh.h>

#include <Eigen/Core>

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
#if SOLUTION
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
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

}  // namespace DebuggingFEM
