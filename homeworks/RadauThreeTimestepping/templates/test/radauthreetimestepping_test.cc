/**
 * @file radauthreetimestepping_test.cc
 * @brief NPDE homework "RadauThreeTimestepping" code
 * @author Tobias Rohner
 * @date 16.03.2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../radauthreetimestepping.h"

namespace RadauThreeTimestepping::test {

TEST(RadauThreeTimestepping, TrapRuleLinFEElemVecProvider) {
  // Get some triangular test mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3);
  // Define some easy functions to test the provider with
  auto f1 = [](const Eigen::Vector2d &x) -> double { return 0; };
  auto f2 = [](const Eigen::Vector2d &x) -> double { return 1; };
  auto f3 = [](const Eigen::Vector2d &x) -> double { return x[0]; };
  // Check the element vector for each triangle
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f1p(f1);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f2p(f2);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider f3p(f3);
  for (const auto tria : mesh_p->Entities(0)) {
    const auto geom = tria->Geometry();
    auto ev1 = f1p.Eval(*tria);
    auto ev2 = f2p.Eval(*tria);
    auto ev3 = f3p.Eval(*tria);
    ASSERT_TRUE(ev1.isApprox(Eigen::Vector3d::Zero()));
    ASSERT_TRUE(ev2.isApprox(
        Eigen::Vector3d::Constant(lf::geometry::Volume(*geom) / 3)));
    Eigen::Vector3d e3_correct =
        lf::geometry::Volume(*geom) / 3 * lf::geometry::Corners(*geom).row(0);
    ASSERT_TRUE(ev3.isApprox(e3_correct));
  }
}

TEST(RadauThreeTimestepping, rhsVectorheatSource) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto &dofh = fespace.LocGlobMap();
  // Assemble the vectors for t=0 and t=0.5
  const Eigen::VectorXd rhs0 =
      RadauThreeTimestepping::rhsVectorheatSource(dofh, 0);
  const Eigen::VectorXd rhs1 =
      RadauThreeTimestepping::rhsVectorheatSource(dofh, 0.5);
  // Get the DOFs on the boundary
  const auto boundary = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);

  // Create a functional for time t=0 and t=0.5
  auto f0 = [](const Eigen::Vector2d &x) -> double {
    if ((x[0] - 0.5) * (x[0] - 0.5) + x[1] * x[1] < 0.25) {
      return 1;
    } else {
      return 0;
    }
  };
  auto f1 = [](const Eigen::Vector2d &x) -> double {
    if (x[0] * x[0] + (x[1] - 0.5) * (x[0] * 0.5) < 0.25) {
      return 1;
    } else {
      return 0;
    }
  };

  // Assume TrapRuleLinFEElemVecProvider works correctly,
  // as it is tested above
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider provider0(f0);
  RadauThreeTimestepping::TrapRuleLinFEElemVecProvider provider1(f1);
  Eigen::VectorXd rhs0_test(dofh.NumDofs());
  Eigen::VectorXd rhs1_test(dofh.NumDofs());
  lf::assemble::AssembleVectorLocally(0, dofh, provider0, rhs0_test);
  lf::assemble::AssembleVectorLocally(0, dofh, provider1, rhs1_test);

  // Make sure the TrapRuleLinFEElemVecProvider is indeed implemented already
  ASSERT_TRUE(rhs0_test.norm() > 1e-10);
  ASSERT_TRUE(rhs1_test.norm() > 1e-10);

  for (int i = 0; i < dofh.NumDofs(); ++i) {
    if (boundary(dofh.Entity(i))) {
      // Check whether boundary DOFs are set to zero
      ASSERT_DOUBLE_EQ(rhs0[i], 0)
          << "Have you forgotten to set the boundary values to zero?";
      ASSERT_DOUBLE_EQ(rhs1[i], 0)
          << "Have you forgotten to set the boundary values to zero?";
    } else {
      // Check whether the value coincides with the one
      // computed previously
      ASSERT_NEAR(rhs0[i], rhs0_test[i], 1e-10);
      ASSERT_NEAR(rhs1[i], rhs1_test[i], 1e-10);
    }
  }
}

TEST(RadauThreeTimestepping, LinFEMassMatrixProvider) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto &dofh = fespace.LocGlobMap();

  // Compare the element matrices for every cell
  RadauThreeTimestepping::LinFEMassMatrixProvider provider;
  for (const auto cell : mesh_p->Entities(0)) {
    const auto geom = cell->Geometry();
    // Compute the correct element matrix
    Eigen::Matrix3d em_correct;
    em_correct.setConstant(lf::geometry::Volume(*geom) / 12);
    em_correct.diagonal() *= 2;
    // Get the element matrix from the element matrix provider
    Eigen::Matrix3d em_prov = provider.Eval(*cell);
    // Compare the matrices
    ASSERT_TRUE(em_correct.isApprox(em_prov));
  }
}

TEST(RadauThreeTimestepping, discreteEvolutionOperator) {
  // Generate a triangular test mesh on [0,1]^2
  const auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1. / 3);
  // Create a DOF handler
  const lf::uscalfe::FeSpaceLagrangeO1<double> fespace(mesh_p);
  const auto &dofh = fespace.LocGlobMap();

  // Compute A and M
  lf::uscalfe::LinearFELaplaceElementMatrix A_provider;
  RadauThreeTimestepping::LinFEMassMatrixProvider M_provider;
  lf::assemble::COOMatrix<double> A_COO(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::COOMatrix<double> M_COO(dofh.NumDofs(), dofh.NumDofs());
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A_provider, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_provider, M_COO);
  const auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);
  const auto selector = [&](unsigned idx) { return bd_flags(dofh.Entity(idx)); };
  RadauThreeTimestepping::dropMatrixRowsColumns(selector, A_COO);
  RadauThreeTimestepping::dropMatrixRowsColumns(selector, M_COO);
  Eigen::SparseMatrix<double> A = A_COO.makeSparse();
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();

  const double dt = 0.0001;
  const Eigen::VectorXd phi = RadauThreeTimestepping::rhsVectorheatSource(dofh, dt);
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver(M + dt*A);
  RadauThreeTimestepping::Radau3MOLTimestepper timestepper(dofh);
  // Test for different mu
  for (unsigned i = 0 ; i < dofh.NumDofs() ; ++i) {
      Eigen::VectorXd mu0(dofh.NumDofs());
      mu0.setZero();
      mu0[i] = 1;
      // Compute mu at the next timestep using discreteEvolutionOperator
      const Eigen::VectorXd mu_dEO = timestepper.discreteEvolutionOperator(0, dt, mu0);
      // Compute mu at the next timestep using implicit euler
      const Eigen::VectorXd mu_iE = solver.solve(M*mu0 + dt*phi);
      // Compare the two mu
      ASSERT_TRUE((mu_dEO-mu_iE).array().abs().maxCoeff() < 1e-4);
  }
}

}  // end namespace RadauThreeTimestepping::test
