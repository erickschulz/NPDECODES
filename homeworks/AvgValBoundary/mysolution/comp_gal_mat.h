/**
 * @ file comp_gal_mat.h
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

#include <boost/filesystem.hpp>

namespace AvgValBoundary {

/**
 * @brief Assembly of general Galerkin matrix for 2nd-order elliptic
 *        boundary value problem
 *
 * @tparam FUNC_ALPHA functor type for diffusion coefficient
 * @tparam FUNC_GAMMA functor type for reaction coefficient
 * @tparam FUNC_BETA functor type for impedance coefficient
 *
 * @param lfe_dofh DofHandler object providing mesh a local-to-global
 *        index mapping for global shape functions.
 *
 * This function computes the finite element Galerkin matrix for the
 * lowest-order Lagrangian FEM and the bilinear form on \f$H^1(\Omega)\f$
 * \f[
     (u,v) \mapsto \int\limits_{\Omega} \alpha(\mathbf{x})
        \mathbf{grad}\,u\cdot\mathfb{grad}\,v +
        \gamma(\mathbf{x})\,u\,v\,\mathrm{d}\mathbf{x} +
        \int\limits_{\Omega}\beta(\mathbf{x})\,u\,v\,\mathrm{d}S(\mathbf{x})
    \f]
 */
template <typename FUNC_ALPHA, typename FUNC_GAMMA, typename FUNC_BETA>
Eigen::SparseMatrix<double> compGalerkinMatrix(
    const lf::assemble::DofHandler& lfe_dofh, FUNC_ALPHA alpha,
    FUNC_GAMMA gamma, FUNC_BETA beta) {
  // obtain mesh and set up fe_space (p.w. linear Lagrangian FEM)
  auto mesh = lfe_dofh.Mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);

  // get the number of degrees of freedom = dimension of FE space
  const lf::base::size_type N_dofs(lfe_dofh.NoDofs());
  // Set up an empty sparse matrix to hold the Galerkin matrix
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Initialize ELEMENT_MATRIX_PROVIDER object
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<double, decltype(alpha),
                                                      decltype(gamma)>
      elmat_builder(fe_space, alpha, gamma);
  // Cell-oriented assembly over the whole computational domain
  lf::assemble::AssembleMatrixLocally(0, lfe_dofh, lfe_dofh, elmat_builder, A);

  // Add contributions of boundary term in the bilinear form using
  // a LehrFEM++ built-in high-level ENTITY\_MATRIX\_PROVIDER class
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh, 1)};
  lf::uscalfe::MassEdgeMatrixProvider<double, decltype(beta),
                                      decltype(bd_flags)>
      edgemat_builder(fe_space, beta, bd_flags);
  lf::assemble::AssembleMatrixLocally(1, lfe_dofh, lfe_dofh, edgemat_builder,
                                      A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();
  return A_crs;
}

double compH1seminorm(const lf::assemble::DofHandler& dofh,
                      const Eigen::VectorXd& u);

/**
 * @brief computes boundary functional as in exercise c)
 * @param dofh DofHandler of FEspace.
 *        u coefficient vector
 */
template <typename FUNCTION>
double compBoundaryFunctional(const lf::assemble::DofHandler& dofh_lfe,
                              const Eigen::VectorXd& u, const FUNCTION& w) {
  double result = 0.;
  /* SOLUTION_BEGIN */
  /* TODO Your implementation goes here! */
  /* SOLUTION_END */
  return result;
}

Eigen::VectorXd solveTestProblem(const lf::assemble::DofHandler& dofh);

std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L);

}  // namespace AvgValBoundary
