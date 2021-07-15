#ifndef SUPG_H
#define SUPG_H

/**
 * @file supg.h
 * @brief Solves the CD BVP based on the SUPG method
 * @author Philippe Peter
 * @date July 2021
 * @copyright Developed at SAM, ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <memory>

#include "cd_tools.h"

namespace ConvectionDiffusion {

/**
 * @brief Computes the local control parameter delta for the SUPG method
 */
double Delta(const lf::mesh::Entity& entity, double eps, Eigen::Vector2d v) {
  double h = Diameter(entity);
  double v_norm = v.lpNorm<Eigen::Infinity>();

  if (v_norm * h / (2 * eps) <= 1.0) {
    return h * h / eps;
  } else {
    return h;
  }
}

/**
 * @brief Computes the local matrices for the stabilization term in the SUPG
 * method based on linear finite elements and the trapezoidal rule .
 */
template <typename FUNCTOR_EPS, typename FUNCTOR_V, typename FUNCTOR_F>
class SupgStabilizationEMP {
 public:
  /**
   * @brief
   * @param eps  diffusion coefficient
   * @param v  velocity field
   * @param f source function
   */
  explicit SupgStabilizationEMP(FUNCTOR_EPS eps, FUNCTOR_V v, FUNCTOR_F f)
      : eps_(eps), v_(v), f_(f) {}

  /**
   * @brief main routine, computes the element matrix on a TRIANGULAR entity
   */
  Eigen::Matrix3d Eval(const lf::mesh::Entity& entity);

  /** @brief Default implementation: all cells are active */
  bool isActive(const lf::mesh::Entity& cell /*entity*/) const { return true; }

 private:
  FUNCTOR_EPS eps_;  // diffusion coefficient
  FUNCTOR_V v_;      // velocity field
  FUNCTOR_F f_;      // source function
};

template <typename FUNCTOR_EPS, typename FUNCTOR_V, typename FUNCTOR_F>
Eigen::Matrix3d SupgStabilizationEMP<FUNCTOR_EPS, FUNCTOR_V, FUNCTOR_F>::Eval(
    const lf::mesh::Entity& entity) {
  LF_ASSERT_MSG(lf::base::RefEl::kTria() == entity.RefEl(),
                "Function only defined for triangular cells");
  const lf::geometry::Geometry* geo_ptr = entity.Geometry();
  Eigen::Matrix3d loc_mat;
  loc_mat.setZero();

  const Eigen::MatrixXd corners = lf::geometry::Corners(*geo_ptr);

  // calculate the gradients of the basis functions.
  // See \lref{cpp:gradbarycoords}, \lref{mc:ElementMatrixLaplLFE} for details.
  Eigen::Matrix3d grad_helper;
  grad_helper.col(0) = Eigen::Vector3d::Ones();
  grad_helper.rightCols(2) = corners.transpose();
  // Matrix with gradients of the local shape functions in its columns
  const Eigen::MatrixXd grad_basis = grad_helper.inverse().bottomRows(2);

  // Local control parameter
  const double delta = Delta(entity, eps_(corners.col(0)), v_(corners.col(0)));

  // weight for the local trapezoidal rule:
  const double weight = delta * lf::geometry::Volume(*geo_ptr) / 3.0;

  // compute local matrix based on local trapezoidal rule:
  for (int i = 0; i < 3; ++i) {
    Eigen::Vector2d corner = corners.col(i);
    Eigen::MatrixXd trial_eval =
        v_(corner).transpose() * grad_basis -
        f_(corner) * Eigen::Vector3d::Ones().transpose();
    Eigen::MatrixXd test_eval = v_(corner).transpose() * grad_basis;

    loc_mat += weight * test_eval.transpose() * trial_eval;
  }

  return loc_mat;
}

/**
 * @brief Solves the Convection-Diffusion BVP with nonhomogeneous Dirichlet
 * boundary conditions using the SUPG method
 */
template <typename DIFFUSION_COEFF, typename CONVECTION_COEFF,
          typename FUNCTOR_F, typename FUNCTOR_G>
Eigen::VectorXd SolveCDBVPSupg(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>>& fe_space,
    DIFFUSION_COEFF eps, CONVECTION_COEFF v, FUNCTOR_F f, FUNCTOR_G g) {
  // Wrap functions into mesh functions
  lf::mesh::utils::MeshFunctionGlobal mf_g{g};
  lf::mesh::utils::MeshFunctionGlobal mf_eps{eps};
  lf::mesh::utils::MeshFunctionGlobal mf_f{f};

  // mesh and dofhanlder
  auto mesh_p = fe_space->Mesh();
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(dofh.NumDofs(), dofh.NumDofs());

  // ASSEMBLE GALERKIN MATRIX
  // First the part corresponding to the laplacian
  lf::uscalfe::ReactionDiffusionElementMatrixProvider laplacian_provider(
      fe_space, mf_eps, lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, laplacian_provider, A);

  // Next part corresponding to the convection term:
  ConvectionElementMatrixProvider convection_provider(v);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, convection_provider, A);

  // Stabilization term:
  SupgStabilizationEMP stabilization_provider(eps, v, f);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stabilization_provider, A);

  // RIGHT-HAND SIDE VECTOR
  Eigen::VectorXd phi(dofh.NumDofs());
  phi.setZero();
  lf::uscalfe::ScalarLoadElementVectorProvider elvec_provider(fe_space, mf_f);
  lf::assemble::AssembleVectorLocally(0, dofh, elvec_provider, phi);
  // IMPOSE DIRICHLET BC
  // Obtain specification for shape functions on edges
  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // Create a dataset of boolean flags indicating edges on the boundary of the
  // mesh
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Fetch flags and values for degrees of freedom located on Dirichlet
  // boundary.
  auto ess_bdc_flags_values{
      lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_g)};

  // Eliminate Dirichlet dofs from linear system
  lf::assemble::FixFlaggedSolutionComponents<double>(
      [&ess_bdc_flags_values](lf::uscalfe::glb_idx_t gdof_idx) {
        return ess_bdc_flags_values[gdof_idx];
      },
      A, phi);

  // SOLVE LINEAR SYSTEM
  Eigen::SparseMatrix A_crs = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_crs);
  Eigen::VectorXd sol_vec = solver.solve(phi);

  return sol_vec;
}

}  // namespace ConvectionDiffusion

#endif  // SUPG_H