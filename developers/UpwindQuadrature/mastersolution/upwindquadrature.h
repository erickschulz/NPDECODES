#ifndef UPWIND_QUADRATURE_H
#define UPWIND_QUADRATURE_H

/**
 * @file upwindquadrature.h
 * @brief NPDE homework template
 * @author Philippe Peter
 * @date June 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <memory>
#include <vector>

namespace UpwindQuadrature {

/**
 * @brief Computes the masses m(p) of all vertices of the mesh
 * @param mesh_p pointer to the mesh.
 * @return data structure containing the masses m(p) for all vertices p of the
 * mesh represented by mesh_p.
 */
lf::mesh::utils::CodimMeshDataSet<double> initializeMasses(
    std::shared_ptr<const lf::mesh::Mesh> mesh_p);

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

/* SAM_LISTING_BEGIN_1 */
template <typename FUNCTOR>
Eigen::Matrix3d UpwindConvectionElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);
  Eigen::Matrix3d loc_mat;

#if SOLUTION
  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  // extract  masses of the corners.
  std::vector<double> local_masses;
  for (const lf::mesh::Entity *sub_ent : entity.SubEntities(2)) {
    local_masses.push_back(masses_(*sub_ent));
  }

  // compute velocities at the vertices:
  Eigen::MatrixXd velocities(2, 3);
  velocities << v_(corners.col(0)), v_(corners.col(1)), v_(corners.col(2));

  // compute rows of the local matrix according to the upwind quadrature scheme.
  // -v(a^j) points into the triangle K (K is the upwind triangle)
  // iff y = a^j - v(a^j) fullfills the two constraints of K
  // that are tight at a^j.
  // In order to simplify the computations we can transform y to the unit
  // triangle and verify that yhat fullfills the two constraints tight at the
  // corresponding corner of the unit triangle.

  // transform a^j - v(a^j) back to the unit triangle
  Eigen::MatrixXd InvJacobian = geo_ptr->Jacobian(corners.col(0)).inverse();
  Eigen::MatrixXd y = corners - velocities;
  Eigen::MatrixXd yhat(2, 3);
  for (int i = 0; i < 3; ++i) {
    yhat.col(i) = InvJacobian * (y.col(i) - corners.col(0));
  }

  // verify the constraints at the first corner of the unit triangle.
  Eigen::Vector2d yhat0 = yhat.col(0);
  if (yhat0(0) >= 0 && yhat0(1) >= 0) {
    if (yhat0(0) == 0 || yhat0(1) == 0) {
      loc_mat.row(0) =
          0.5 * local_masses[0] * velocities.transpose().row(0) * grad_basis;
    } else {
      loc_mat.row(0) =
          local_masses[0] * velocities.transpose().row(0) * grad_basis;
    }
  } else {
    loc_mat.row(0) = Eigen::Vector3d::Zero();
  }

  // verify the constraints at the second corner of the unit triangle.
  Eigen::Vector2d yhat1 = yhat.col(1);
  if (yhat1(1) >= 0 && yhat1(1) <= 1 - yhat1(0)) {
    if (yhat1(1) == 0.0 || yhat1(1) == 1.0 - yhat1(0)) {
      loc_mat.row(1) =
          0.5 * local_masses[1] * velocities.transpose().row(1) * grad_basis;
    } else {
      loc_mat.row(1) =
          local_masses[1] * velocities.transpose().row(1) * grad_basis;
    }
  } else {
    loc_mat.row(1) = Eigen::Vector3d::Zero();
  }

  // verify the constraints at the third corner of the unit triangle.
  Eigen::Vector2d yhat2 = yhat.col(2);
  if (yhat2(0) >= 0 && yhat2(1) <= 1 - yhat2(0)) {
    if (yhat2(0) == 0 || yhat2(1) == 1 - yhat2(0)) {
      loc_mat.row(2) =
          0.5 * local_masses[2] * velocities.transpose().row(2) * grad_basis;
    } else {
      loc_mat.row(2) =
          local_masses[2] * velocities.transpose().row(2) * grad_basis;
    }
  } else {
    loc_mat.row(2) = Eigen::Vector3d::Zero();
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return loc_mat;
}
/* SAM_LISTING_END_1 */

}  // namespace UpwindQuadrature

#endif  // UPWIND_QUADRATURE_H
