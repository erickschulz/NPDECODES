/**
 * @ file main.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */

#include "comp_gal_mat.h"

/* SAM_LISTING_BEGIN_1 */
int main() {
  // read in mesh and set up finite element space
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();
  // obtain dofh for lagrangian finite element space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Solve test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute H1 seminorm of the solution
  double h1s_norm = AvgValBoundary::compH1seminorm(dofh, mu);
  // constant identity mesh function
  lf::uscalfe::MeshFunctionConstant mf_identity{1.};
  // compute boundary functional
  double boundary_functional =
      AvgValBoundary::compBoundaryFunctional(dofh, mu, mf_identity);

  std::cout << "H1s-norm: " << h1s_norm << "\n";
  std::cout << "F: " << boundary_functional << "\n";

  /* SOLUTION_BEGIN */
  auto results = AvgValBoundary::approxBoundaryFunctionalValues(7);
  double ground_truth = results[6].second;
  std::cout << "N_Dofs    Error\n";
  for (int i = 0; i < 6; i++) {
    std::cout << results[i].first << "  "
              << std::abs(results[i].second - ground_truth) << "\n";
  }
  /* SOLUTION_END */
  return 0;
}
/* SAM_LISTING_END_1 */
