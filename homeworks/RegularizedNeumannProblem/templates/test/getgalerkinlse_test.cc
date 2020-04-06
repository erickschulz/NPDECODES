/**
 * @file getgalerkinlse_test.cc
 * @brief NPDE homework "RegularizedNeumannProblem" code
 * @author Philipp Lindenberger
 * @date 03.04.2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../getgalerkinlse.h"

namespace RegularizedNeumann::test {

// Test for getGalerkinLSE method (not used in problem, but example code in PDF)
TEST(RegularizedNeumannProblem, getGalerkinLSE) {
  // We test on a triangular mesh
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4);
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto& dofh = fe_space_p->LocGlobMap();

  // Initialize constant mesh functions
  auto mf_f0 = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto mf_f2 = lf::mesh::utils::MeshFunctionConstant(2.0);
  auto mf_h0 = lf::mesh::utils::MeshFunctionConstant(0.0);
  auto mf_h3 = lf::mesh::utils::MeshFunctionConstant(3.0);

  // Get boundary edges
  auto bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1);

  // Compute exact solution with f = 2.0, h=0.0
  Eigen::VectorXd f2_h0_solution(dofh.NumDofs());
  f2_h0_solution.setZero();

  for (const lf::mesh::Entity* cell : mesh_p->Entities(0)) {
    auto geo_ptr = cell->Geometry();
    auto area = lf::geometry::Volume(*geo_ptr);
    auto glob_dof_indices = dofh.GlobalDofIndices(*cell);
    int local_dofs = dofh.NumLocalDofs(*cell);
    for (int i = 0; i < local_dofs; i++) {
      f2_h0_solution(glob_dof_indices[i]) += area / 3.0;
    }
  }

  f2_h0_solution *= 2.0;

  // Compute exact solution with f = 0.0, h=3.0
  Eigen::VectorXd f0_h3_solution(dofh.NumDofs());
  f0_h3_solution.setZero();

  for (const lf::mesh::Entity* edge : mesh_p->Entities(1)) {
    if (bd_flags(*edge)) {
      auto geo_ptr = edge->Geometry();
      auto area = lf::geometry::Volume(*geo_ptr);
      auto glob_dof_indices = dofh.GlobalDofIndices(*edge);
      int local_dofs = dofh.NumLocalDofs(*edge);
      for (int i = 0; i < local_dofs; i++) {
        f0_h3_solution(glob_dof_indices[i]) += area / 2.0;
      }
    }
  }

  f0_h3_solution *= 3.0;

  // Get FEM Solution from getGalerkinLSE
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f2, mf_h0);
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair_2 =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f0, mf_h3);
  std::pair<Eigen::SparseMatrix<double>, Eigen::VectorXd> sol_pair_3 =
      RegularizedNeumannProblem::getGalerkinLSE(fe_space_p, mf_f2, mf_h3);

  // We only test rhs since A is simple laplacian
  auto approx_rhs_f2_h0 = sol_pair.second;
  auto approx_rhs_f0_h3 = sol_pair_2.second;
  auto approx_rhs_f2_h3 = sol_pair_3.second;

  EXPECT_NEAR((f2_h0_solution - approx_rhs_f2_h0).norm(), 0.0, 1.0e-8);
  EXPECT_NEAR((f0_h3_solution - approx_rhs_f0_h3).norm(), 0.0, 1.0e-8);

  // Test solution of f=2.0, h=3.0 -> Should be equal to adding up our two exact
  // solutions
  EXPECT_NEAR((f0_h3_solution + f2_h0_solution - approx_rhs_f2_h3).norm(), 0.0,
              1.0e-8);
}

}  // namespace RegularizedNeumann::test