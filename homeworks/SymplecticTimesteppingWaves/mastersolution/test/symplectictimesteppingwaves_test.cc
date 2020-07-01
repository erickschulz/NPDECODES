/**
 * @file symplectictimesteppingwaves_test.cc
 * @brief NPDE homework "SymplecticTimesteppingWaves" code
 * @author Amélie Loher
 * @date 09.04.2020
 * @copyright Developed at ETH Zurich
 */

#include "../symplectictimesteppingwaves.h"
#include "../symplectictimesteppingwaves_assemble.h"
#include "../symplectictimesteppingwaves_ode.h"

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <Eigen/SparseLU>

#include <utility>

#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/uscalfe/uscalfe.h>

namespace SymplecticTimesteppingWaves::test {

TEST(SymplecticTimesteppingWaves, sympTimesteppingHarmonicOscillatorODE) {
  unsigned int m = 10;

  Eigen::Vector2d sol = sympTimesteppingHarmonicOscillatorODE(m);

  Eigen::Vector2d ref_sol;
  ref_sol << 0.000675997, 1;

  double tol = 1.0e-4;
  ASSERT_NEAR((sol - ref_sol).lpNorm<Eigen::Infinity>(), 0.0, tol);
}

TEST(SymplecticTimesteppingWaves, assembleGalerkinMatrix) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(
      std::move(mesh_factory), CURRENT_SOURCE_DIR "/../../meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto one = [](Eigen::Vector2d x) -> double { return 1.0; };

  Eigen::SparseMatrix<double> galMat =
      assembleGalerkinMatrix(fe_space, one, one, one);

  Eigen::MatrixXd sol(galMat);

  Eigen::MatrixXd ref_sol(5, 5);
  ref_sol << 1.75, 0.1875, 0, 0.1875, -0.958333, 0.1875, 1.75, 0.1875, 0,
      -0.958333, 0, 0.1875, 1.75, 0.1875, -0.958333, 0.1875, 0, 0.1875, 1.75,
      -0.958333, -0.958333, -0.958333, -0.958333, -0.958333, 4.16667;

  double tol = 1.0e-4;
  ASSERT_NEAR((sol - ref_sol).lpNorm<Eigen::Infinity>(), 0.0, tol);
}

TEST(SymplecticTimesteppingWaves, computeEnergies) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(
      std::move(mesh_factory), CURRENT_SOURCE_DIR "/../../meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto c = [](Eigen::Vector2d x) -> double { return 1.0 + x.dot(x); };

  SympTimestepWaveEq stepper(fe_space, c);

  Eigen::VectorXd p = Eigen::VectorXd::Ones(5);
  Eigen::VectorXd q = 2 * Eigen::VectorXd::Ones(5);

  double sol = stepper.computeEnergies(p, q);

  double ref_sol = 3.83333;

  double tol = 1.0e-4;
  ASSERT_NEAR((sol - ref_sol), 0.0, tol);
}

TEST(SymplecticTimesteppingWaves, solvewave) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(
      std::move(mesh_factory), CURRENT_SOURCE_DIR "/../../meshes/simple.msh");
  auto mesh_p = reader.mesh();

  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto c = [](Eigen::Vector2d x) -> double { return 1.0 + x.dot(x); };

  Eigen::VectorXd u0 = Eigen::VectorXd::Ones(5);
  Eigen::VectorXd v0 = 2 * Eigen::VectorXd::Ones(5);

  double T = 1.0;
  unsigned int m = 10;

  std::pair<Eigen::VectorXd, Eigen::VectorXd> sol =
      solvewave(fe_space, c, u0, v0, T, m);

  Eigen::VectorXd ref_sol_first(5);
  ref_sol_first << 2.01346, 1.76133, 1.53227, 1.76133, 1.77866;

  Eigen::VectorXd ref_sol_second(11);
  ref_sol_second << 2.83333, 2.83333, 2.83334, 2.83336, 2.83338, 2.8334,
      2.83342, 2.83346, 2.8335, 2.83355, 2.83359;

  double tol = 1.0e-4;
  ASSERT_NEAR((sol.first - ref_sol_first).lpNorm<Eigen::Infinity>(), 0.0, tol);
  ASSERT_NEAR((sol.second - ref_sol_second).lpNorm<Eigen::Infinity>(), 0.0,
              tol);
}

} /* namespace SymplecticTimesteppingWaves::test */
