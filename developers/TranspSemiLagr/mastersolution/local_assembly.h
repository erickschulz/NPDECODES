/**
 * @file loacl_assembly.h
 * @brief NPDE homework TranspSemiLagr local assembly code
 * @author Philippe Peter
 * @date November 2020
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>

#include <Eigen/Core>
#include <memory>

#ifndef LOCAL_ASSEMBLY_H
#define LOCAL_ASSEMBLY_H

namespace TranspSemiLagr {

/**
 * @headerfile local_assembly.h
 * @brief computes the local vectors that arrise in the fully discreteized first
 * order Semi-Lagrangian scheme
 *
 * @tparam FUNCTOR function that defines the vector vallued velocity field v
 * @tparam MESH_FUNCTION meshfunction that defines the solution u0 at the
 * previous time step
 *
 * The local vector is given by
 *
 * v^K_i = 1/3*|K|*u0(p_i - \tau*v(p_i)) , i=1,...,3
 *
 * where p_i is the i-th node of the triangle K.
 */
template <typename FUNCTOR, typename MESH_FUNCTION>
class UpwindLagrangianElementVectorProvider {
  static_assert(lf::mesh::utils::isMeshFunction<MESH_FUNCTION>);

 public:
  /**
   * @brief
   * @param v velocity field
   * @param tau  step size
   * @param mesh_p pointer to uthe underlying mesh
   * @param U0 solution at previous time step
   */
  UpwindLagrangianElementVectorProvider(
      FUNCTOR v, double tau, std::shared_ptr<const lf::mesh::Mesh> mesh_p,
      MESH_FUNCTION U0)
      : v_(v), tau_(tau), mesh_p_(std::move(mesh_p)), U0_(U0) {}

  /**
   * @brief actual computation of the element vector
   * @param entity reference to the triangle for which the matrix is needed
   * @return element vector
   */
  Eigen::Vector3d Eval(const lf::mesh::Entity& entity);

  /**
   * @brief If true, then an element is takein into account during assembly
   */
  bool isActive(const lf::mesh::Entity& /* entity */) const { return true; }

 private:
  /**
   * @brief computes the bullback of the identity function on K,
   * i.e. the local coordinates of a point in K
   * @param x point at which the pullback is evaluated (not necessarily in K)
   * @return (\Phi^* Id)(x)
   */
  Eigen::Vector2d pull_back_(Eigen::Vector2d x,
                             const lf::geometry::Geometry& geo);

  double tau_;        // time step size
  FUNCTOR v_;         // velocity field
  MESH_FUNCTION U0_;  // solution at previous time step
  std::shared_ptr<const lf::mesh::Mesh> mesh_p_;  // pointer to underlying mesh
};

template <typename FUNCTOR, typename MESH_FUNCTION>
Eigen::Vector3d
UpwindLagrangianElementVectorProvider<FUNCTOR, MESH_FUNCTION>::Eval(
    const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  // initialize element vector
  Eigen::Vector3d result;
  result.setZero();

  // retrieve geometric information
  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);

  // evaluate the three components of the element vector
  for (unsigned int i = 0; i < 3; ++i) {
    // compute evaluation point for i-th coordinate
    Eigen::Vector2d y = corners.col(i) - tau_ * v_(corners.col(i));

    // Since U0 is given as a mesh function, we need to find the triangle
    // in which the evaluation point y lies.
    for (const lf::mesh::Entity* cell : mesh_p_->Entities(0)) {
      // compute the local coordinates of y wrt. cell
      Eigen::Vector2d yhat = pull_back_(y, *(cell->Geometry()));
      // check if the local coordinates  yhat of y wrt. cell lie in the
      // reference triangle. This is the case if and only if y lies in cell
      if (yhat(0) >= 0 && yhat(1) >= 0 && yhat.sum() <= 1) {
        // evaluate mesh function at the correct local coordinates.
        result(i) = U0_(*cell, yhat)[0] * area / 3.0;
      }
    }
  }
  return result;
}

template <typename FUNCTOR, typename MESH_FUNCTION>
Eigen::Vector2d
UpwindLagrangianElementVectorProvider<FUNCTOR, MESH_FUNCTION>::pull_back_(
    Eigen::Vector2d x, const lf::geometry::Geometry& geo) {
  const Eigen::MatrixXd corners = lf::geometry::Corners(geo);
  const Eigen::MatrixXd InvJacobian = geo.Jacobian(corners.col(0)).inverse();
  return InvJacobian * (x - corners.col(0));
}

/**
 * @brief local_assembly.h
 * @brief computes the local mass matrix on a triangle  on the 2d trapezoidal
 * rule on the triangle (mass-lumping).
 *
 * @tparam FUNCTOR function that defines the scalar valued coefficient c
 *
 *
 * This helper class corresponds to the local matrix for the local biliniear
 * form \int_K c(x)*b^i(x)*b^j(x)dx
 *
 * where the integral is evalauted using the 2d trapezoidal rule on the
 * triangle., i.e.
 *
 * \int_K f(x) dx \approx |K|/3*(f(p_1) + f(p_2) + f(p_3))
 *
 */
template <typename FUNCTOR>
class LumpedMassElementMatrixProvider {
 public:
  /**
   * @brief
   * @param c scalar valued coefficient function
   */
  LumpedMassElementMatrixProvider(FUNCTOR c) : c_(c) {}

  /**
   * @brief actual computation of the element matrix
   * @param entity reference to the triangle for which the matrix is needed.
   * @return a 3x3 dense matrix, containing the element matrix.
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);

  /**
   * @brief If true, an element is taken into account during assembly.
   */
  bool isActive(const lf::mesh::Entity& /* entity */) const { return true; }

 private:
  FUNCTOR c_;  // coefficient function
};

template <typename FUNCTOR>
Eigen::Matrix3d LumpedMassElementMatrixProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  Eigen::Matrix3d result = Eigen::Matrix3d::Zero();

  // evaluation points of the quadrature rule
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const double area = lf::geometry::Volume(*geo_ptr);

  // evalauate element matrix based on the 2d trapezoidal rule
  for (int i = 0; i < 3; ++i) {
    result(i, i) = area / 3.0 * c_(corners.col(i));
  }

  return result;
}

}  // namespace TranspSemiLagr

#endif