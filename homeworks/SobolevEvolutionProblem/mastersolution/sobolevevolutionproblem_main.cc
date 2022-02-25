/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author
 * @ date
 * @ copyright Developed at SAM, ETH Zurich
 */

#include <lf/mesh/test_utils/test_meshes.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <utility>

#include "sobolevevolutionproblem.h"

int main(int /*argc*/, char** /*argv*/) {
  // Obtain a triangular mesh of the unit square from the collection of
  // test meshes
  std::shared_ptr<lf::mesh::Mesh> mesh_p =
      lf::mesh::test_utils::GenerateHybrid2DTestMesh(3, 1.0 / 3.0);
  // Simplest Lagrangian finite-element space
  auto fe_space_p =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  const int N_dofs = fe_space_p->LocGlobMap().NumDofs();
  // Coefficient function beta
  auto beta = [](const Eigen::Vector2d x) -> double {
    return (1.0 + x.squaredNorm());
  };
  // Coefficient function alpha
  auto alpha = [](const Eigen::Vector2d x) -> double {
    return std::exp(x.norm());
  };
  // Wrap them into MeshFunctions
  lf::mesh::utils::MeshFunctionGlobal mf_beta{beta};
  lf::mesh::utils::MeshFunctionGlobal mf_alpha{alpha};

  // Butcher tableau for classical RK-SSM of order 4
  Eigen::MatrixXd RK_mat(4, 4);
  RK_mat.setZero();
  RK_mat(1, 0) = 0.5;
  RK_mat(2, 1) = 0.5;
  RK_mat(3, 2) = 1.0;
  Eigen::VectorXd b = ((Eigen::VectorXd(4) << 1, 2, 2, 1).finished()) / 6.0;
  // Initial coefficient vector
  Eigen::VectorXd mu0 = Eigen::VectorXd::Constant(N_dofs, 1.0);

  // Reference solution
  const Eigen::VectorXd mu_final = SobolevEVP::solveRKSobEvl(
      fe_space_p, mf_beta, mf_alpha, mu0, 1.0, RK_mat, b, 1000);

  // Repeat numerical integration with different numbers of timesteps
  int M = 10;
  for (int l = 0; l < 4; l++, M *= 2) {
    const Eigen::VectorXd mu = SobolevEVP::solveRKSobEvl(
        fe_space_p, mf_beta, mf_alpha, mu0, 1.0, RK_mat, b, M);
    std::cout << "M = " << M
              << ": timestepping error = " << (mu - mu_final).norm()
              << std::endl;
  }
  return 0;
}
