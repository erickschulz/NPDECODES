/**
 * @ file avg_val_boundary.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "avgvalboundary.h"
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/refinement/refinement.h>

namespace AvgValBoundary
{

/**
 * @brief computes H1 seminorm over the computational domain
 * @param dofh DofHandler of FEspace.
 *        u coefficient vector
 */
/* SAM_LISTING_BEGIN_1 */
double compH1seminorm(const lf::assemble::DofHandler &dofh,
                      const Eigen::VectorXd &u)
{
  double result = 0.;
#if SOLUTION
  // constant identity mesh function
  lf::mesh::utils::MeshFunctionConstant mf_identity{1.};
  // constant zero mesh function
  lf::mesh::utils::MeshFunctionConstant mf_zero{0.};

  // compute Stiffness matrix
  // alpha := 1, beta, gamma := 0
  auto A = compGalerkinMatrix(dofh, mf_identity, mf_zero, mf_zero);

  // compute seminorm using the Stiffness matrix
  result = u.transpose() * A * u;
  result = std::sqrt(result);
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_1 */

/**
 * @brief solves pde for some simple test problem with
 *        alpha = beta = gamma := 1. and load f := 1.
 * @param dofh DofHandler of FEspace.
 */
Eigen::VectorXd solveTestProblem(const lf::assemble::DofHandler &dofh)
{
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
    unsigned int L)
{
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3.);
  std::shared_ptr<lf::refinement::MeshHierarchy> meshes =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, L);
  return meshes;
}

/**
 * @brief Compute the boundary functional values for a sequence of `L+1` meshes
 * @param L The number of refinements to perform to get the sequence of meshes
 * @returns A vector of pairs containing the number of DOFs and value of the
 *	    boundary functional for each level
 */
/* SAM_LISTING_BEGIN_5 */
std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L)
{
  std::vector<std::pair<unsigned int, double>> result{};
#if SOLUTION
  auto meshes = generateTestMeshSequence(L - 1);
  int num_meshes = meshes->NumLevels();
  for (int level = 0; level < num_meshes; ++level)
  {
    auto mesh_p = meshes->getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs(dofh.NumDofs());
    lf::mesh::utils::MeshFunctionConstant mf_identity{1.};
    // compute galerkin matrix
    auto A = AvgValBoundary::compGalerkinMatrix(dofh, mf_identity, mf_identity,
                                                mf_identity);
    // compute load vector for f(x) = x.norm()
    auto f = [](Eigen::Vector2d x) -> double { return x.norm(); };
    lf::mesh::utils::MeshFunctionGlobal mf_f{f};
    Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
    phi.setZero();
    lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_identity)>
        elvec_builder(fe_space, mf_identity);
    AssembleVectorLocally(0, dofh, elvec_builder, phi);

    // solve system of equations
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    Eigen::VectorXd u = solver.solve(phi);

    // set up weight function
    auto w = [](Eigen::Vector2d x) -> double { return 1.; };
    lf::mesh::utils::MeshFunctionGlobal mf_w{w};
    double functional_value = compBoundaryFunctional(dofh, u, mf_w);

    result.push_back({N_dofs, functional_value});
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return result;
}
/* SAM_LISTING_END_5 */

} // namespace AvgValBoundary
