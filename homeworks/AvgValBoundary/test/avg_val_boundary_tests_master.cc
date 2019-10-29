#include <gtest/gtest.h>
#include <boost/filesystem.hpp>
#include "../mastersolution/comp_gal_mat.h"

namespace AvgValBoundary::test {
TEST(AvgValBoundary, TestH1SemiNorm) {
  // read in mesh and set up finite element space
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";

  // constant identity mesh function
  lf::uscalfe::MeshFunctionConstant mf_identity{1.};

  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute H1 seminorm
  double h1s_norm = AvgValBoundary::compH1seminorm(dofh, mu);

  ASSERT_NEAR(h1s_norm, 0.151178, 0.00001);
}

TEST(AvgValBoundary, TestBoundaryFunctional) {
  // read in mesh and set up finite element space
  boost::filesystem::path here = __FILE__;
  auto square_path = here.parent_path().parent_path() / "meshes/square.msh";

  // constant identity mesh function
  lf::uscalfe::MeshFunctionConstant mf_identity{1.};

  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), square_path.string());
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute boundary functional
  double boundary_functional =
      AvgValBoundary::compBoundaryFunctional(dofh, mu, mf_identity);
  std::cout << boundary_functional;
  ASSERT_NEAR(boundary_functional, 0.880602, 0.00001);
}
}  // namespace AvgValBoundary::test
