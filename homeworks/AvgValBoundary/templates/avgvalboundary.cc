/**
 * @ file avgvalboundary.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "avgvalboundary.h"

#include <memory>
#include <utility>
#include <vector>

#include <Eigen/SparseCore>

#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace AvgValBoundary {

/**
 * @brief computes H1 seminorm over the computational domain
 * @param dofh DofHandler of FEspace.
 *        u coefficient vector
 */
/* SAM_LISTING_BEGIN_1 */
double compH1seminorm(const lf::assemble::DofHandler &dofh,
                      const Eigen::VectorXd &u) {
  double result = 0.;
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_1 */

/**
 * @brief solves pde for some simple test problem with
 *        alpha = beta = gamma := 1. and load f := 1.
 * @param dofh DofHandler of FEspace.
 */
Eigen::VectorXd solveTestProblem(const lf::assemble::DofHandler &dofh) {
  // constant identity mesh function
  lf::mesh::utils::MeshFunctionConstant mf_identity{1.};

  // obtain Galerkin matrix for alpha = beta = gamma := 1.
  auto A = AvgValBoundary::compGalerkinMatrix(dofh, mf_identity, mf_identity,
                                              mf_identity);

  // Set up load vector
  auto mesh = dofh.Mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::base::size_type N_dofs(dofh.NumDofs());
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_identity)>
      elvec_builder(fe_space, mf_identity);
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // solve system of linear equations
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  Eigen::VectorXd mu = solver.solve(phi);

  return mu;
}

/** @brief generate sequence of nested triangular meshes with L+1 levels */
std::shared_ptr<lf::refinement::MeshHierarchy> generateTestMeshSequence(
    unsigned int L) {
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3.);
  std::shared_ptr<lf::refinement::MeshHierarchy> meshes =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, L);
  return meshes;
}

/**
 * @brief Compute the boundary functional values for a sequence of L meshes
 * @param L The number of meshes in the sequence
 * @returns A vector of pairs containing the number of DOFs and value of the
 *	    boundary functional for each level
 */
/* SAM_LISTING_BEGIN_5 */
std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L) {
  std::vector<std::pair<unsigned int, double>> result{};
  //====================
  // Your code goes here
  //====================
  return result;
}
/* SAM_LISTING_END_5 */

}  // namespace AvgValBoundary
