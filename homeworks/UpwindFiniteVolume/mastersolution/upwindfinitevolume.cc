/**
 * @file upwindfinitevolume.cc
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include "upwindfinitevolume.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>

namespace UpwindFiniteVolume {

/* SAM_LISTING_BEGIN_1 */
Eigen::Matrix<double, 2, 3> gradbarycoordinates(
    const Eigen::Matrix<double, 2, 3> &triangle) {
  Eigen::Matrix3d X;
  // Solve for the coefficients of the barycentric coordinate functions
  X.block<3, 1>(0, 0) = Eigen::Vector3d::Ones();
  X.block<3, 2>(0, 1) = triangle.transpose();
  return X.inverse().block<2, 3>(1, 0);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Vector2d computeCircumcenters(const Eigen::Vector2d &a1,
                                     const Eigen::Vector2d &a2,
                                     const Eigen::Vector2d &a3) {

  Eigen::Vector2d mp1 = 0.5 * (a1 + a2);
  Eigen::Vector2d mp2 = 0.5 * (a2 + a3);

  Eigen::Vector2d dir1 = a2 - a1;
  Eigen::Vector2d dir2 = a3 - a2;

  Eigen::Matrix2d dir;
  dir << dir1(1), dir2(1), -dir1(0), -dir2(0);

  Eigen::Vector2d p = dir.colPivHouseholderQr().solve(mp2 - mp1);

  Eigen::Vector2d center = mp1 + p(0) * dir.col(0);

  // Barycentric coordinates from
  // Christer Ericson’s ’Real-Time Collision Detection’

  Eigen::Vector2d v0 = a2 - a1, v1 = a3 - a1, v2 = center - a1;
  double d00 = v0.dot(v0);
  double d01 = v0.dot(v1);
  double d11 = v1.dot(v1);
  double d20 = v2.dot(v0);
  double d21 = v2.dot(v1);

  double denom = d00 * d11 - d01 * d01;
  double v = (d11 * d20 - d01 * d21) / denom;
  double w = (d00 * d21 - d01 * d20) / denom;
  double u = 1.0 - v - w;

  if (v <= 0 || w <= 0 || u <= 0) {
    throw std::runtime_error("Obtused triangle!");
    return Eigen::Vector2d::Zero();
  } else {
    return center;
  }

  // Alternative:
  /*
  // Calculate the midpoint of the edges
  const Eigen::Vector2d midpoint1 = (a1 + a2) * 0.5;
  const Eigen::Vector2d midpoint2 = (a2 + a3) * 0.5;
  const Eigen::Vector2d midpoint3 = (a1 + a3) * 0.5;

  // Calculate the slope of the line perpendicular to the edges
  Eigen::Matrix<double, 2, 3> corners;
  corners.col(0) = a1;
  corners.col(1) = a2;
  corners.col(2) = a3;

  Eigen::MatrixXd barycenters = gradbarycoordinates(corners);

  // Reorder the vectors corresponding to the edges
  Eigen::Vector2d n1 = barycenters.col(2).normalized();
  Eigen::Vector2d n2 = barycenters.col(0).normalized();
  Eigen::Vector2d n3 = barycenters.col(1).normalized();

  // Check obtuse triangle
  double alpha1 = acos(n1.dot(-n3));
  double alpha2 = acos(n1.dot(-n2));
  double alpha3 = acos(n2.dot(-n3));

  double rad_90deg = M_PI / 2.0;
  if (alpha1 >= rad_90deg || alpha2 >=  rad_90deg || alpha3 >=  rad_90deg) {
    std::cout << "Obtused triangle!" << std::endl;
    std::cout << "alpha1 " << alpha1 * 180.0 / M_PI << "degree" << std::endl;
    std::cout << "alpha2 " << alpha2 * 180.0 / M_PI << "degree" << std::endl;
    std::cout << "alpha3 " << alpha3 * 180.0 / M_PI << "degree" << std::endl;
    throw std::runtime_error("Obtused triangle!");
  }

  // Compute intersection the two vectors
  // midpoint1 + n1 * t1 = midpoint2 + n2 * t2
  // (midpoint1x - midpoint2x) = t2 * n2x - t1 * n1x
  // (midpoint1y - midpoint2y) = t2 * n2y - t1 * n1y
  // solving for t1:
  double t1 = ((midpoint1[0] - midpoint2[0]) * n2[1] -
               (midpoint1[1] - midpoint2[1]) * n2[0]) /
              (n1[1] * n2[0] - n1[0] * n2[1]);

  Eigen::Vector2d intersection = midpoint1 + n1 * t1;

  return intersection;
  */
}
/* SAM_LISTING_END_2 */

}  // namespace UpwindFiniteVolume
