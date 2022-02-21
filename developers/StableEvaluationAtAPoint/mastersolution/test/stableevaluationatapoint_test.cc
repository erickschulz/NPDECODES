/**
 * @file stableevaluationatapoint_test.cc
 * @brief NPDE homework StableEvaluationAtAPoint
 * @author Amélie Loher
 * @date 29/04/2020
 * @copyright Developed at ETH Zurich
 */

#include "../stableevaluationatapoint.h"

#include <gtest/gtest.h>
#include <lf/fe/fe.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <cmath>
#include <memory>
#include <utility>

TEST(StableEvaluationAtAPoint, PSL) {
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR "/../../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  const auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };

  const Eigen::Vector2d x(0.3, 0.4);

  const double val = StableEvaluationAtAPoint::PSL(mesh_p, u, x);

  const double ref_val = 0.15525;

  double tol = 1.e-4;

  ASSERT_NEAR(std::abs(ref_val - val), 0.0, tol);
}

TEST(StableEvaluationAtAPoint, PDL) {
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR "/../../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  const auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };

  const Eigen::Vector2d x(0.3, 0.4);

  const double val = StableEvaluationAtAPoint::PDL(mesh_p, u, x);

  const double ref_val = -0.484226;

  double tol = 1.e-4;

  ASSERT_NEAR(std::abs(ref_val - val), 0.0, tol);
}

TEST(StableEvaluationAtAPoint, PointEval) {
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR "/../../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  double error = StableEvaluationAtAPoint::PointEval(mesh_p);

  double ref_error = 0.0784387;

  double tol = 1.e-4;

  ASSERT_NEAR(std::abs(ref_error - error), 0.0, tol);
}

TEST(StableEvaluationAtAPoint, Jstar) {
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR
                                 "/../../meshes/square7.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };

  lf::mesh::utils::MeshFunctionGlobal mf_u{u};
  Eigen::VectorXd uFE = lf::fe::NodalProjection(*fe_space, mf_u);

  const Eigen::Vector2d x(0.3, 0.4);

  double val = StableEvaluationAtAPoint::Jstar(fe_space, uFE, x);

  double ref_val = u(x);

  double tol = 1.e-2;

  ASSERT_NEAR(val, ref_val, tol);
}

/*
TEST(StableEvaluationAtAPoint, stab_pointEval) {
  auto mesh_factory_init = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader_init(std::move(mesh_factory_init),
                                 CURRENT_SOURCE_DIR "/../../meshes/square.msh");
  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader_init.mesh();

  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  const auto u = [](Eigen::Vector2d x) -> double {
    Eigen::Vector2d one(1.0, 0.0);
    return std::log((x + one).norm());
  };

  const Eigen::Vector2d x(0.3, 0.4);

  double val = StableEvaluationAtAPoint::StablePointEvaluation(fe_space, u, x);

  double ref_val = 0.0;

  double tol = 1.e-4;

  ASSERT_NEAR(std::abs(ref_val - val), 0.0, tol);
}
*/
