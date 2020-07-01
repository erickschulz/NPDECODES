/**
 * @file upwindquadrature.h
 * @brief NPDE homework template
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <memory>
#include <vector>

#include <Eigen/Core>

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

namespace UpwindQuadrature {

/**
 * @headerfile upwindquadrature.h
 * @brief Computes the local matrices for the convection term based on linear
 * finite elements and the trapezoidal rule (standard Galerkin approach).
 *
 * @tparam FUNCTOR function that defines the vector valued velocity
 * coefficient v.
 */
template <typename FUNCTOR>
class ConvectionElementMatrixProvider {
 public:
  /**
   * @brief
   * @param v functor for the velocity field
   */
  explicit ConvectionElementMatrixProvider(FUNCTOR v) : v_(v) {}

  /**
   * @brief main routine for the computation of element matrices
   * @param enttity reference to the TRIANGULAR cell for which the element
   * matrix should be computed.
   * @return a 3x3 matrix containing the element matrix.
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);

  /** @brief Default implementation: all cells are active */
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }

 private:
  FUNCTOR v_;  // functor for the velocity field.
};

template <typename FUNCTOR>
Eigen::Matrix3d ConvectionElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);

  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  // evaluate velocity field at the vertices:
  Eigen::MatrixXd velocities(2, 3);
  velocities << v_(corners.col(0)), v_(corners.col(1)), v_(corners.col(2));

  // compute local matrix using local trapezoidal rule:
  loc_mat = lf::geometry::Volume(*geo_ptr) / 3.0 * velocities.transpose() *
            grad_basis;

  return loc_mat;
}

// At any vertex a^j of a triangle, the vector field -v(a^j) may either point
// into the triangle, along an edge of the triangle or outside the triangle.
enum class Direction { INWARDS, OUTWARDS, ALONG_EDGE };

/**
 * @brief Computes for all corners a^j of a triangle the direction of the vector
 * field -v(a^j)
 * @param geo Geometry object describing the triangle
 * @param velocities values of the vector field v, evaluated at the corners of
 * the triangle described by geo, stored in a 2x3 matrix.
 * @return A vector containing for all 3 corners the corresponding direction of
 * -v(a^j)
 */
std::vector<Direction> opposite_velocity_directions(
    const lf::geometry::Geometry &geo, const Eigen::MatrixXd &velocities);

/**
 * @headerfile upwindquadrature.h
 * @brief Computes the local matrices for the convection term based on the
 * upwind quadrature scheme introced in the exercise.
 * @tparam FUNCTOR function that defines the vector valued velocity
 * coefficient v.
 */
template <typename FUNCTOR>
class UpwindConvectionElementMatrixProvider {
 public:
  /**
   * @brief
   * @param v functor for the velocity field
   * @param masses data structure storing the masses m(a^j) for all vertices of
   * the mesh.
   */
  explicit UpwindConvectionElementMatrixProvider(
      FUNCTOR v, lf::mesh::utils::CodimMeshDataSet<double> masses)
      : v_(v), masses_(masses) {}

  /**
   * @brief main routine for the computation of element matrices.
   * @param entity reference to the TRIANGULAR cell for which the element
   * matrix should be computed.
   * @return a 3x3 matrix containing the element matrix.
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);

  /** @brief Default implementation: all cells are active */
  bool isActive(const lf::mesh::Entity & /*entity*/) const { return true; }

 private:
  FUNCTOR v_;  // velocity field
  lf::mesh::utils::CodimMeshDataSet<double>
      masses_;  // masses of all vertices of the mesh.
};

template <typename FUNCTOR>
Eigen::Matrix3d UpwindConvectionElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);
  Eigen::Matrix3d loc_mat;

  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  // compute velocities at the vertices:
  Eigen::MatrixXd velocities(2, 3);
  velocities << v_(corners.col(0)), v_(corners.col(1)), v_(corners.col(2));

  // determine the direction of -v at the corners of the triangle
  std::vector<Direction> directions =
      opposite_velocity_directions(*geo_ptr, velocities);

  // extract  masses of the corners.
  std::vector<double> local_masses;
  for (const lf::mesh::Entity *sub_ent : entity.SubEntities(2)) {
    local_masses.push_back(masses_(*sub_ent));
  }

  // compute rows of the local matrix according to the upwind quadrature scheme.
  for (int i = 0; i < 3; ++i) {
    switch (directions[i]) {
      case Direction::INWARDS:
        loc_mat.row(i) =
            local_masses[i] * velocities.transpose().row(i) * grad_basis;
        break;
      case Direction::OUTWARDS:
        loc_mat.row(i) = Eigen::Vector3d::Zero();
        break;
      case Direction::ALONG_EDGE:
        loc_mat.row(i) =
            0.5 * local_masses[i] * velocities.transpose().row(i) * grad_basis;
        break;
    }
  }

  return loc_mat;
}

/**
 * @brief Computes the masses m(p) of all vertices of the mesh
 * @param mesh_p pointer to the mesh.
 * @return data structure containing the masses m(p) for all vertices p of the
 * mesh represented by mesh_p.
 */
lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

}  // namespace UpwindQuadrature
