/**
 * @ file avg_val_boundary.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "comp_gal_mat.h"

namespace AvgValBoundary {

/**
 * @brief computes H1 seminorm over the computational domain
 * @param dofh DofHandler of FEspace.
 *        u coefficient vector
 */
/* SAM_LISTING_BEGIN_1 */
double compH1seminorm(const lf::assemble::DofHandler& dofh,
                      const Eigen::VectorXd& u) {
  double result = 0.;
  /* SOLUTION_BEGIN */
  // constant identity mesh function
  lf::uscalfe::MeshFunctionConstant mf_identity{1.};
  // constant zero mesh function
  lf::uscalfe::MeshFunctionConstant mf_zero{0.};

  // compute Stiffness matrix
  // alpha := 1, beta, gamma := 0
  auto A = compGalerkinMatrix(dofh, mf_identity, mf_zero, mf_zero);

  // compute seminorm using the Stiffness matrix
  result = u.transpose() * A * u;
  result = std::sqrt(result);
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_1 */

/**
 * @brief solves pde for some simple test problem with
 *        alpha = beta = gamma := 1. and load f := 1.
 * @param dofh DofHandler of FEspace.
 */
Eigen::VectorXd solveTestProblem(const lf::assemble::DofHandler& dofh) {
  // constant identity mesh function
  lf::uscalfe::MeshFunctionConstant mf_identity{1.};

  // obtain Galerkin matrix for alpha = beta = gamma := 1.
  auto A = AvgValBoundary::compGalerkinMatrix(dofh, mf_identity, mf_identity,
                                              mf_identity);

  // Set up load vector
  auto mesh = dofh.Mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::base::size_type N_dofs(dofh.NoDofs());
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

/** @brief generate sequence of nested triangular meshes with L levels */
std::shared_ptr<lf::refinement::MeshHierarchy> generateTestMeshSequence(
    unsigned int L) {
  auto mesh = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3.);
  std::shared_ptr<lf::refinement::MeshHierarchy> meshes =
      lf::refinement::GenerateMeshHierarchyByUniformRefinemnt(mesh, L);
  return meshes;
}

/* SAM_LISTING_BEGIN_5 */
std::vector<std::pair<unsigned int, double>> approxBoundaryFunctionalValues(
    unsigned int L) {
  std::vector<std::pair<unsigned int, double>> result{};
  /* SOLUTION_BEGIN */
  auto meshes = generateTestMeshSequence(L - 1);
  int num_meshes = meshes->NumLevels();
  for (int level = 0; level < num_meshes; ++level) {
    auto mesh_p = meshes->getMesh(level);
    // Set up global FE space; lowest order Lagrangian finite elements
    auto fe_space =
        std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
    // Obtain local->global index mapping for current finite element space
    const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
    const lf::base::size_type N_dofs(dofh.NoDofs());
    lf::uscalfe::MeshFunctionConstant mf_identity{1.};
    // compute galerkin matrix
    auto A = AvgValBoundary::compGalerkinMatrix(dofh, mf_identity, mf_identity,
                                                mf_identity);
    // compute load vector for f(x) = x.norm()
    auto f = [](Eigen::Vector2d x) -> double { return x.norm(); };
    lf::uscalfe::MeshFunctionGlobal mf_f{f};
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
    lf::uscalfe::MeshFunctionGlobal mf_w{w};
    double functional_value = compBoundaryFunctional(dofh, u, mf_w);

    result.push_back({N_dofs, functional_value});
  }
  /* SOLUTION_END */
  return result;
}
/* SAM_LISTING_END_5 */

}  // namespace AvgValBoundary
