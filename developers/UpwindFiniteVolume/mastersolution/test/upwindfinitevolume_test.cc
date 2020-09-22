/**
 * @file upwindfinitevolume_test.cc
 * @brief NPDE homework UpwindFiniteVolume code
 * @author Philipp Egg
 * @date 08.09.2020
 * @copyright Developed at ETH Zurich
 */

#include "../upwindfinitevolume.h"

#include <gtest/gtest.h>
#include <lf/assemble/assemble.h>
#include <lf/io/gmsh_reader.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/mesh.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <memory>
#include <vector>

namespace UpwindFiniteVolume::test {

TEST(UpwindFiniteVolume, computeCircumcenters) {
  Eigen::Vector2d a1 = Eigen::Vector2d(0, 0);
  Eigen::Vector2d a2 = Eigen::Vector2d(1, 0);
  Eigen::Vector2d a3 = Eigen::Vector2d(0.3, 1);

  Eigen::Vector2d result01 =
      UpwindFiniteVolume::computeCircumcenters(a1, a2, a3);
  Eigen::Vector2d result02 =
      UpwindFiniteVolume::computeCircumcenters(a3, a1, a2);
  Eigen::Vector2d result03 =
      UpwindFiniteVolume::computeCircumcenters(a2, a3, a1);

  int throwed_error = 0;
  try {
    Eigen::Vector2d a4 = Eigen::Vector2d(1, 1);
    Eigen::Vector2d temp = UpwindFiniteVolume::computeCircumcenters(a1, a2, a4);
  } catch (const std::exception&) {
    throwed_error = 1;
  }

  double tol = 1.0e-6;
  Eigen::Vector2d ref_res = Eigen::Vector2d(0.5, 0.395);
  ASSERT_NEAR(0.0, (ref_res - result01).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (ref_res - result02).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_NEAR(0.0, (ref_res - result03).lpNorm<Eigen::Infinity>(), tol);
  ASSERT_EQ(1, throwed_error);
}

TEST(ElementMatrixProvider, Eval) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/mesh.msh");

  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  double eps = 1e-6;
  auto v = [](Eigen::Vector2d x) -> Eigen::Vector2d {
    return Eigen::Vector2d(1, 0.5 - x[1]);
  };

  const lf::assemble::UniformFEDofHandler cur_dofh(
      mesh_p, {{lf::base::RefEl::kPoint(), 1},
               {lf::base::RefEl::kSegment(), 0},
               {lf::base::RefEl::kTria(), 0},
               {lf::base::RefEl::kQuad(), 0}});
  int N_dofs = cur_dofh.NumDofs();

  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  UpwindFiniteVolume::ElementMatrixProvider my_mat_provider(v, eps);
  lf::assemble::AssembleMatrixLocally(0, cur_dofh, cur_dofh, my_mat_provider,
                                      A);

  const Eigen::MatrixXd A_crs = A.makeSparse();

  Eigen::MatrixXd A_ref(N_dofs, N_dofs);
  A_ref << 0, 0, 0, 0, 0.125, 0, 0, 0, 0, 0, 0, 0.0416667, 0, 0, 0, 0,
      0.0564236, 0, 0, 0, 0, 0, -0.15191, 0, 0, 0, 0, 0.0416667, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.15191, 0, 0, 0, 0, 0.0416667, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.125, 0.0416667,
      0, 0, 0, 0.0564236, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.125, 0.00584975, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0386423, 0.20752, 0, 0, 0, 0, 0.125, 0, 0, 0,
      -0.0888792, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0908909, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, -0.367792, 1.75493e-08, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 1.75493e-08, -0.367792, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.125, 0, 0, 0, 0, 0, -0.0888792, 0, 0, 0, 0, 0.0908909, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00584975, -0.125, 0, 0, 0, 0,
      0.0386423, 0, 0, 0, 0, 0, 0.20752, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      -0.0416667, 4.97828e-08, 0, 0, 0.179547, 0, 0, 0, 0.160075, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 4.97828e-08, -0.0416667, 0, 0, 0, 0, 0.179547, 0,
      0.160075, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.380219, 0.106647,
      0, 0.106647, 0, 0, 0, 0.0892857, 0, 0, 0, 0.0269097, 0, 0, 0, 0, 0.187769,
      0, 0, 0, 0, 0, -0.411068, 0, 0, 0, 0, 0, 0.12511, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0.125118, 0, -0.353218, 0, 0, 0, 0, 0, 0.156841, 0,
      0.0269097, 0, 0, 0, 0, 0.187769, 0, 0, 0, 0, 0, 0, 0, 0, -0.411068, 0, 0,
      0, 0.12511, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.125118, 0, 0, 0,
      -0.353218, 0.156841, 0, 0, 0, 0, 0, 0, 0, 0, 0.0830295, 0, 0, 0, 0, 0, 0,
      0.0146152, 0, 0, 0.21353, 0, -0.36436, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.100752, 0, 0.078606, 0, 0.078606, 0, -0.32015, 0, 0, 0, 0, 0,
      0, 0, 0, 0.138357, 0.138357, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.339506,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0830295, 0, 0, 0, 0.0146152, 0.21353, 0, 0,
      0, 0, 0, 0, -0.36436;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (A_ref - A_crs).lpNorm<Eigen::Infinity>(), tol);
}

TEST(ElementVectorProvider, Eval) {
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  lf::io::GmshReader reader(std::move(mesh_factory),
                            CURRENT_SOURCE_DIR "/mesh.msh");

  std::shared_ptr<lf::mesh::Mesh> mesh_p = reader.mesh();

  double eps = 1e-6;
  auto u = [](Eigen::Vector2d x) -> double {
    return std::sin(4 * M_PI * x[1]) * std::exp(-x[0]);
  };

  auto f = [&u, eps](Eigen::Vector2d x) -> double {
    return eps * (1.0 - 16.0 * M_PI * M_PI) * u(x) - 2 * u(x) +
           (0.5 - x[1]) * 4 * M_PI * std::cos(4.0 * M_PI * x[1]) *
               std::exp(-x[0]);
  };

  const lf::assemble::UniformFEDofHandler cur_dofh(
      mesh_p, {{lf::base::RefEl::kPoint(), 1},
               {lf::base::RefEl::kSegment(), 0},
               {lf::base::RefEl::kTria(), 0},
               {lf::base::RefEl::kQuad(), 0}});
  int N_dofs = cur_dofh.NumDofs();

  Eigen::VectorXd phi(N_dofs);
  phi.setZero();
  UpwindFiniteVolume::ElementVectorProvider my_vec_provider(f);
  lf::assemble::AssembleVectorLocally(0, cur_dofh, my_vec_provider, phi);

  Eigen::VectorXd phi_ref(N_dofs);
  phi_ref << 0.169079, 0.0622006, -0.0622006, -0.169079, 0.166834, 0.119542,
      0.00933809, -0.00933809, -0.119542, -0.166834, -0.0260237, 0.0260237,
      2.96479e-13, 0.0444807, 0.0668182, -0.0444807, -0.0668182, -0.139637,
      -3.23969e-13, 2.74695e-13, 0.139637;

  double tol = 1.0e-6;
  ASSERT_NEAR(0.0, (phi_ref - phi).lpNorm<Eigen::Infinity>(), tol);
}

}  // namespace UpwindFiniteVolume::test
