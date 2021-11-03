/**
 * @file upwindquadrature.cc
 * @brief NPDE homework template
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include "upwindquadrature.h"

namespace UpwindQuadrature {

/**
 * @brief Computes for all corners a^j of a triangle the direction of the vector
 * field -v(a^j)
 * @param geo Geometry object describing the triangle
 * @param velocities values of the vector field v, evaluated at the corners of
 * the triangle described by geo, stored in a 2x3 matrix.
 * @return A vector, containing for all 3 corners the corresponding direction of
 * -v(a^j)
 */
std::vector<Direction> opposite_velocity_directions(
    const lf::geometry::Geometry &geo, const Eigen::MatrixXd &velocities) {
  std::vector<Direction> res(3);
  const Eigen::MatrixXd corners = lf::geometry::Corners(geo);

#if SOLUTION
  // Observe that -v(a^j) points into the triangle K
  // iff y = a^j - v(a^j) fullfills the two constraints of K
  // which are tight at a^j.
  // In order to simplify the computations we can transform y to the unit
  // triangle and verify that yhat fullfills the two constraints tight at the
  // corresponding corner of the unit triangle.

  // transform a^j - v(a^j) back to the unit triangle
  Eigen::MatrixXd InvJacobian = geo.Jacobian(corners.col(0)).inverse();
  Eigen::MatrixXd y = corners - velocities;
  Eigen::MatrixXd yhat(2, 3);
  for (int i = 0; i < 3; ++i) {
    yhat.col(i) = InvJacobian * (y.col(i) - corners.col(0));
  }

  // verify the constraints at the first corner of the unit triangle.
  Eigen::Vector2d yhat0 = yhat.col(0);
  if (yhat0(0) >= 0 && yhat0(1) >= 0) {
    if (yhat0(0) == 0 || yhat0(1) == 0) {
      res[0] = Direction::ALONG_EDGE;
    } else {
      res[0] = Direction::INWARDS;
    }
  } else {
    res[0] = Direction::OUTWARDS;
  }

  // verify the constraints at the second corner of the unit triangle.
  Eigen::Vector2d yhat1 = yhat.col(1);
  if (yhat1(1) >= 0 && yhat1(1) <= 1 - yhat1(0)) {
    if (yhat1(1) == 0.0 || yhat1(1) == 1.0 - yhat1(0)) {
      res[1] = Direction::ALONG_EDGE;
    } else {
      res[1] = Direction::INWARDS;
    }
  } else {
    res[1] = Direction::OUTWARDS;
  }

  // verify the constraints at the third corner of the unit triangle.
  Eigen::Vector2d yhat2 = yhat.col(2);
  if (yhat2(0) >= 0 && yhat2(1) <= 1 - yhat2(0)) {
    if (yhat2(0) == 0 || yhat2(1) == 1 - yhat2(0)) {
      res[2] = Direction::ALONG_EDGE;
    } else {
      res[2] = Direction::INWARDS;
    }
  } else {
    res[2] = Direction::OUTWARDS;
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return res;
}

/**
 * @brief Computes the masses m(p) of all vertices of the mesh
 * @param mesh_p pointer to the mesh.
 * @return Datastructure containing the masses m(p) for all vertices p of the
 * mesh represented by mesh_p.
 */
lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p) {
  lf::mesh::utils::CodimMeshDataSet<double> masses(mesh_p, 2, 0.0);

#if SOLUTION
  // compute masses using a cell-based approach.
  for (const lf::mesh::Entity *entity : mesh_p->Entities(0)) {
    const lf::geometry::Geometry *geo_ptr = entity->Geometry();
    double area = lf::geometry::Volume(*geo_ptr);
    for (const lf::mesh::Entity *corner : entity->SubEntities(2)) {
      masses(*corner) += area / 3.0;
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return masses;
}

}  // namespace UpwindQuadrature
