/**
 * @file
 * @brief UNITESTS for NPDE homework ErrorEstimatesForTraces
 * @author Erick Schulz
 * @date 28/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <gtest/gtest.h>

#include <string>

#include "../teelaplrobinassembly.h"

namespace ErrorEstimatesForTraces::test {

TEST(ErrorEstimatesForTraces, TestResult) {
  // Load mesh into a Lehrfem++ object
  auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
  const lf::io::GmshReader reader(std::move(mesh_factory),
                                  CURRENT_SOURCE_DIR "/../../meshes/hex4.msh");
  auto mesh_p = reader.mesh();  // mesh_p has type shared_ptr<lf::mesh::Mesh>

  // Finite element space
  auto fe_space = std::make_shared<linear_lagrange>(mesh_p);

  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  // Solve the boundary value problem with Robin boundary conditions
  Eigen::VectorXd sol_vec = solveBVP(fe_space);

  // Integrate the solution sol_vec over the flagged edges
  double bd_functional_val = bdFunctionalEval(fe_space, sol_vec);

  ASSERT_NEAR(bd_functional_val, 2.081541059732923, 1e-08);
}

}  // namespace ErrorEstimatesForTraces::test
