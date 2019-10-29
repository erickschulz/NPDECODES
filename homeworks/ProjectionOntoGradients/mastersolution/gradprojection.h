#ifndef __GRADPROJECTION_H
#define __GRADPROJECTION_H

#include <iostream>
#include <memory>

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include "lf/mesh/test_utils/test_meshes.h"

namespace ProjectionOntoGradients {

// Start of sub-problem e)
/* SAM_LISTING_BEGIN_1 */
class ElementMatrixProvider {
 private:
  using coord_t = Eigen::Vector2d;

 public:
  Eigen::Matrix3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
};
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
Eigen::Matrix3d ElementMatrixProvider::Eval(const lf::mesh::Entity &entity) {
  const lf::geometry::Geometry *geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;

  /* BEGIN_SOLUTION */
  // get area of the entity
  const double area = lf::geometry::Volume(*geo_ptr);

  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.block(0, 0, 3, 1) = Eigen::Vector3d::Ones();
  grad_helper.block(0, 1, 3, 2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().block(1, 0, 2, 3);

  loc_mat = area * (grad_basis.transpose() * grad_basis);
  /* END_SOLUTION */
  return loc_mat;
}
/* SAM_LISTING_END_2 */
// End of sub-problem e)

// Start of sub-problem g)
/* SAM_LISTING_BEGIN_3 */
template <typename FUNCTOR>
class GradProjRhsProvider {
 private:
  FUNCTOR f_;
  using coord_t = Eigen::Vector2d;

 public:
  explicit GradProjRhsProvider(FUNCTOR f) : f_(f) {}
  Eigen::Vector3d Eval(const lf::mesh::Entity &entity);
  bool isActive(const lf::mesh::Entity &entity) const { return true; }
};

/* SAM_LISTING_END_3 */

/* SAM_LISTING_BEGIN_4 */
template <typename FUNCTOR>
Eigen::Vector3d GradProjRhsProvider<FUNCTOR>::Eval(
    const lf::mesh::Entity &entity) {
  Eigen::Vector3d loc_vec;

  const lf::geometry::Geometry *geo_ptr = entity.Geometry();

  /* BEGIN_SOLUTION */
  // get area of the entity
  const double area = lf::geometry::Volume(*geo_ptr);

  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");

  // calculate center of mass
  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);
  const coord_t c = corners.rowwise().sum() / 3.;

  // value of the functor f at center of mass of the entity
  const Eigen::Vector2d func_value = f_(c);

  // calculate the gradients of the basis functions
  Eigen::Matrix3d grad_helper;
  grad_helper.block(0, 0, 3, 1) = Eigen::Vector3d::Ones();
  grad_helper.block(0, 1, 3, 2) = corners.transpose();
  // vector of the gradients of the basis function evaluated at c
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().block(1, 0, 2, 3);

  loc_vec = area * (grad_basis.transpose() * func_value);
  /* END_SOLUTION */

  return loc_vec;
}
/* SAM_LISTING_END_4 */
// End of sub-problem e)

// Start of sub-problem h)
template <typename FUNCTOR>
Eigen::VectorXd projectOntoGradients(const lf::assemble::DofHandler &dofh,
                                     FUNCTOR f) {
  const std::size_t N_dofs = dofh.NoDofs();
  Eigen::VectorXd sol_vec;

  /* SAM_LISTING_BEGIN_5 */
  // ASSEMBLE GLOBAL MATRIX
  /* BEGIN_SOLUTION */
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  ProjectionOntoGradients::ElementMatrixProvider my_mat_provider;
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, my_mat_provider, A);
  /* END_SOLUTION */
  // ASSEMBLE GLOBAL RHS VECOR
  /* BEGIN_SOLUTION */
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  ProjectionOntoGradients::GradProjRhsProvider my_vec_provider(f);
  // co-dimension 0 because we locally assemble on cells
  lf::assemble::AssembleVectorLocally(0, dofh, my_vec_provider, phi);
  /* END_SOLUTION */
  // ENFORCE HOMOGENEOUS DIRICHLET BOUNDARY CONDITIONS
  /* BEGIN_SOLUTION */
  // We do this by selecting the DOFs on the boundary and setting the
  // values to zero. Note, that we could also use this to set the boundary
  // to any other value (if bc were not homogeneous)

  // To select the right DOFs we need a selector:
  // * for all DOFs on the boundary return true and the value on the boundary
  // * for all other DOFs return false (and whatever as value)
  const double boundary_val = 0;
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(dofh.Mesh(), 2)};
  auto my_selector = [&dofh, &bd_flags, &boundary_val](unsigned int dof_idx) {
    if (bd_flags(dofh.Entity(dof_idx))) {
      return (std::pair<bool, double>(true, boundary_val));
    } else {
      // interior node: the value we return here does not matter
      return (std::pair<bool, double>(false, 42.0));
    }
  };
  // Since we know the values on the boundary we know the solution on these
  // DOFs and we can write the Galerkin LSE in block format and solve only
  // for the unknown coefficients. This modification is done by the following
  // function fix\_flagged\_solution\_components(). We use the selector we have
  // defined above.
  // See \lref{sec:essbdc} for explanations.
  lf::assemble::fix_flagged_solution_components<double>(my_selector, A, phi);
  /* END_SOLUTION */
  /* SAM_LISTING_END_5 */
  // CALCULATE THE SOLUTION
  /* BEGIN_SOLUTION */
  // Convert from triplet to CRS format
  const Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  sol_vec = solver.solve(phi);
  /* END_SOLUTION */
  return sol_vec;
}
// End of sub-problem h)

}  // namespace ProjectionOntoGradients

#endif  // define __GRADPROJECTION_H
