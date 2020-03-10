/**
 * @file regularizedneumannproblem_test.cc
 * @brief NPDE homework RegularizedNeumannProblem code
 * @author Christian Mitsch, Philippe Peter
 * @date March 2020
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <memory>
#include <utility>

#include <Eigen/Core>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../regularizedneumannproblem.h"

namespace RegularizedNeumann::test {

// Test for sub-exercise c with functions h and f being constant
TEST(RegularizedNeumannProblem, solution_test_dropDof_const) {

  // read mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c = RegularizedNeumannProblem::getGalerkinLSE_dropDof(fe_space, f, h);

  // Now we compare the returned values to hard coded results
  const double eps = 1e-10;

  // Check Matrix
  Eigen::MatrixXd solution_mat_c(5, 5);
  solution_mat_c << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 1,
      -1, 0, -1, -1, -1, 4;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_c - result_c.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_c(5);
  solution_vec_c << 0, 1.1666666667, 1.1666666667, 1.1666666667, 0.33333333333;
  // Compare with expected results
  EXPECT_NEAR((solution_vec_c - result_c.second).norm(), 0.0, eps);
}

// Test for sub-exercise c with functions h and f not being constant
TEST(RegularizedNeumannProblem, solution_test_dropDof_gen) {

  // read  mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_c = RegularizedNeumannProblem::getGalerkinLSE_dropDof(fe_space, f, h);

  // Now we compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_c(5, 5);
  solution_mat_c << 1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 1,
      -1, 0, -1, -1, -1, 4;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_c - result_c.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_c(5);
  solution_vec_c << 0, 1.1666666667, 1.91667, 1.1666666667, 0.33333333333;
  // Compare with expected results
  EXPECT_NEAR((solution_vec_c - result_c.second).norm(), 0.0, eps);
}

// Test for sub-exercise f with functions h and f being constant
TEST(RegularizedNeumannProblem, solution_test_augment_const) {

  // read mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_f = RegularizedNeumannProblem::getGalerkinLSE_augment(fe_space, f, h);

  // Compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_f(6, 6);
  solution_mat_f << 1, 0, 0, 0, -1, 0.1666667, 0, 1, 0, 0, -1, 0.1666667, 0, 0,
      1, 0, -1, 0.1666667, 0, 0, 0, 1, -1, 0.1666667, -1, -1, -1, -1, 4,
      0.3333333, 0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_f - result_f.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_f(6);
  solution_vec_f << 1.1666666667, 1.1666666667, 1.1666666667, 1.1666666667,
      0.33333333333, 0;
  EXPECT_NEAR((solution_vec_f - result_f.second).norm(), 0.0, eps);
}

// Test for sub-exercise f with functions h and f not being constant
TEST(RegularizedNeumannProblem, solution_test_augment_gen) {

  // read  mesh file
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/meshes/test.msh");
  auto mesh_p = reader.mesh();

  // source and boundary functions for testing
  const auto f = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });
  const auto h = lf::mesh::utils::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return x(0) + x(1); });

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // Compute solution
  auto result_f = RegularizedNeumannProblem::getGalerkinLSE_augment(fe_space, f, h);

  // Compare the returned values to hard coded results
  const double eps = 1e-5;

  // Check Matrix
  Eigen::MatrixXd solution_mat_f(6, 6);
  solution_mat_f << 1, 0, 0, 0, -1, 0.1666667, 0, 1, 0, 0, -1, 0.1666667, 0, 0,
      1, 0, -1, 0.1666667, 0, 0, 0, 1, -1, 0.1666667, -1, -1, -1, -1, 4,
      0.3333333, 0.166667, 0.166667, 0.166667, 0.166667, 0.3333333, 0;
  // Compare with expected results
  EXPECT_NEAR((solution_mat_f - result_f.first).norm(), 0.0, eps);

  // Check rhs vector
  Eigen::VectorXd solution_vec_f(6);
  solution_vec_f << 0.416667, 1.1666666667, 1.91667, 1.1666666667,
      0.33333333333, 0;
  EXPECT_NEAR((solution_vec_f - result_f.second).norm(), 0.0, eps);
}

} // namespace RegularizedNeumannProblem::test
