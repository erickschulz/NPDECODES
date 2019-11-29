/**
 * @file electrostaticforce.cc
 * @brief ElectrostaticForce
 * @author Erick Schulz
 * @date 27.11.2019
 * @copyright Developed at ETH Zurich
 */

// HACK:
#undef SOLUTION
#define SOLUTION 1

#include <gtest/gtest.h>

#include "../electrostaticforce.h"

#include <Eigen/Core>

namespace ElectrostaticForce::test {

TEST(ElectrostaticForce, computeExactForce) {
  Eigen::Vector2d exactForce = ElectrostaticForce::computeExactForce();

  double tol = 1.0e-4;
  std::cout << exactForce << std::endl;
  ASSERT_NEAR(13.0776, exactForce(0), tol);
  ASSERT_NEAR(0.0, exactForce(1), tol);
}

TEST(ElectrostaticForce, computeForceDomainFunctional) {
  // READ MESH INTO LEHRFEMPP
  // Load mesh into a Lehrfem++ object
  std::string mesh_file =
      CURRENT_SOURCE_DIR "/../../meshes/emforce" + std::to_string(4) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>
  // Finite element space
  auto fe_space_p = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // SOLVE POISSON DIRICHLET BVP
  Eigen::VectorXd approx_sol = ElectrostaticForce::solvePoissonBVP(fe_space_p);
  Eigen::Vector2d approx_force_domain_functional =
      ElectrostaticForce::computeForceDomainFunctional(fe_space_p, approx_sol);

  double tol = 1.0e-4;
  ASSERT_NEAR(13.0776, approx_force_domain_functional(0), tol);
  ASSERT_NEAR(0.0, approx_force_domain_functional(1), tol);
}

TEST(ElectrostaticForce, computeForceBoundaryFunctional) {
  // READ MESH INTO LEHRFEMPP
  // Load mesh into a Lehrfem++ object
  std::string mesh_file =
      CURRENT_SOURCE_DIR "/../../meshes/emforce" + std::to_string(4) + ".msh";
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory), mesh_file);
  auto mesh_p = reader.mesh();  // type shared_ptr< const lf::mesh::Mesh>
  // Finite element space
  auto fe_space_p = std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // SOLVE POISSON DIRICHLET BVP
  Eigen::VectorXd approx_sol = ElectrostaticForce::solvePoissonBVP(fe_space_p);
  Eigen::Vector2d approx_force_boundary_functional =
      ElectrostaticForce::computeForceBoundaryFunctional(fe_space_p, approx_sol);

  double tol = 1.0e-4;
  ASSERT_NEAR(13.0776, approx_force_boundary_functional(0), tol);
  ASSERT_NEAR(0.0, approx_force_boundary_functional(1), tol);
}

}  // namespace ElectrostaticForce::test
