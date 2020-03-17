/**
 * @ file avgvalboundary_test.cc
 * @ brief NPDE homework AvgValBoundary code
 * @ author Simon Meierhans, edited by Oliver Rietmann
 * @ date 11.03.2019
 * @ copyright Developed at ETH Zurich
 */


#include <gtest/gtest.h>

#include <memory>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/io/io.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../avgvalboundary.h"

namespace AvgValBoundary::test {

constexpr char mesh_file[] = CURRENT_SOURCE_DIR "/../../meshes/square.msh";

TEST(AvgValBoundary, TestH1SemiNorm) {
  // constant identity mesh function
  lf::mesh::utils::MeshFunctionConstant mf_identity{1.0};

  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute H1 seminorm
  double h1s_norm = AvgValBoundary::compH1seminorm(dofh, mu);

  ASSERT_NEAR(h1s_norm, 0.151178, 0.00001);
}

TEST(AvgValBoundary, TestBoundaryFunctional) {
  // constant identity mesh function
  lf::mesh::utils::MeshFunctionConstant mf_identity{1.0};

  // obtain dofh for lagrangian finite element space
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh = reader.mesh();
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh);
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};

  // get solution of test problem
  Eigen::VectorXd mu = AvgValBoundary::solveTestProblem(dofh);
  // compute boundary functional
  double boundary_functional =
      AvgValBoundary::compBoundaryFunctional(dofh, mu, mf_identity);
  ASSERT_NEAR(boundary_functional, 0.880602, 0.00001);
}
} // namespace AvgValBoundary::test
