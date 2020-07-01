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

  //====================
  // Your code goes here
  //====================
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

  //====================
  // Your code goes here
  //====================
  return masses;
}

}  // namespace UpwindQuadrature
