/**
 * @file transpsemilagr_test.cc
 * @brief NPDE homework TranspSemiLagr test file
 * @author Philippe Peter
 * @date November 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <cmath>
#include <memory>

#include <gtest/gtest.h>

#include <Eigen/Core>

#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include "../transpsemilagr.h"

namespace TranspSemiLagr::test {

// verifies that the vector u of nodal values satisfies zero boundary conditions
void verify_zero_bc(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const Eigen::VectorXd& u) {
  auto mesh_p = fe_space->Mesh();
  const auto& dof_h = fe_space->LocGlobMap();
  auto boundary_nodes = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2);

  for (int i = 0; i < u.size(); ++i) {
    if (boundary_nodes(dof_h.Entity(i))) {
      EXPECT_NEAR(u[i], 0.0, 1.0E-6);
    }
  }
}

// The first set of test cases verifies, that all methods produce solutions that
// satisfy the zero dirichlet boundary conditions.

TEST(ReactionStep, boundary_conditions) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u_function));
  auto c_function = [](Eigen::VectorXd x) { return x(0) + x(1); };

  verify_zero_bc(fe_space, u0_vector);
  Eigen::VectorXd u_new = reaction_step(fe_space, u0_vector, c_function, 1.0);
  verify_zero_bc(fe_space, u_new);
}

TEST(SemiLagrStep, boundary_conditions) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u_function));

  auto v = [](Eigen::Vector2d x) {
    return (Eigen::Vector2d() << -x(1) + 3.0 * x(0) * x(0), x(0)).finished();
  };

  verify_zero_bc(fe_space, u0_vector);
  Eigen::VectorXd u_new = semiLagr_step(fe_space, u0_vector, v, 1.0);
  verify_zero_bc(fe_space, u_new);
}

TEST(solverot, boundary_conditions) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto u_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u_function));

  Eigen::VectorXd u_new = solverot(fe_space, u0_vector, 10, 1.0);
  verify_zero_bc(fe_space, u_new);
}

TEST(solvetrp, boundary_conditions) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  auto u_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };

  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u_function));
  Eigen::VectorXd u_new = solvetrp(fe_space, u0_vector, 10, 1.0);
  verify_zero_bc(fe_space, u_new);
}

// The following two test cases test the reaction_step function
// based on an exact solution u(x,t). Even though this exact solution
// is a polynomial of degree four in x, it is still the case that
// the approximate solution to the evolution problem is a good approximation
// of the nodal projection of the exact solution at time T. Since
// the reaction_step relies on a lumped mass matrix, the system of
// equations decouble and for each component ODE the RK-scheme is applied.

// u(x,t) = x0*(1-x0)*x1*(1-x1)*exp(2t)
TEST(ReactionStep, exact_solution_constant_c) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // solution function at initial time t=0
  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };

  // solution function at final time t=1
  auto u1_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0)) * std::exp(2 * 1.0);
  };

  // c function in the ODE:
  auto c_function = [](Eigen::VectorXd /*x*/) { return 2.0; };

  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));
  Eigen::VectorXd u1_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u1_function));

  // apply 100 reaction steps with time step 0.01
  for (int i = 0; i < 100; ++i) {
    u0_vector = reaction_step(fe_space, u0_vector, c_function, 0.01);
  }
  EXPECT_NEAR((u0_vector - u1_vector).norm(), 0.0, 1.0E-3);
}

// based on exact solution u(x,t) = x0*(1-x0)*x1*(1-x1)*exp(2t*x0 + t*x1))
TEST(ReactionStep, exact_solution_linear_c) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  // solution function at initial time t=0
  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };

  // solution function at final time t=1
  auto u1_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0)) * std::exp(2 * x(0) + x(1));
  };

  // c function in the ODE:
  auto c_function = [](Eigen::VectorXd x) { return 2.0 * x(0) + x(1); };

  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));
  Eigen::VectorXd u1_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u1_function));

  // apply 100 reaction steps with time step 0.01
  for (int i = 0; i < 100; ++i) {
    u0_vector = reaction_step(fe_space, u0_vector, c_function, 0.01);
  }
  EXPECT_NEAR((u0_vector - u1_vector).norm(), 0.0, 1.0E-3);
}

// The following few testcases verify that solverot (and solvetrp) satisfy the
// following: If solverot simulates N timesteps of the PDE on [0,T] with initial
// condition u_0 the result should be the same as simulating N/2 timesteps on
// [0,T/2] and then another N/2 timestpes on [0,T/2] with the updated initial
// conditin, since all coefficients are indepentent of time.
TEST(solverot, consistency_1) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solverot(fe_space, u0_vector, 2, 1.0);
  Eigen::VectorXd sol_2 =
      solverot(fe_space, solverot(fe_space, u0_vector, 1, 0.5), 1, 0.5);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

TEST(solverot, consistency_2) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solverot(fe_space, u0_vector, 10, 1.0);
  Eigen::VectorXd sol_2 =
      solverot(fe_space, solverot(fe_space, u0_vector, 5, 0.5), 5, 0.5);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

TEST(solverot, consistency_3) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solverot(fe_space, u0_vector, 10, 1.0);
  Eigen::VectorXd sol_2 =
      solverot(fe_space, solverot(fe_space, u0_vector, 1, 0.1), 9, 0.9);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

TEST(solvetrp, consistency_1) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solvetrp(fe_space, u0_vector, 2, 1.0);
  Eigen::VectorXd sol_2 =
      solvetrp(fe_space, solvetrp(fe_space, u0_vector, 1, 0.5), 1, 0.5);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

TEST(solvetrp, consistency_2) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solvetrp(fe_space, u0_vector, 10, 1.0);
  Eigen::VectorXd sol_2 =
      solvetrp(fe_space, solvetrp(fe_space, u0_vector, 5, 0.5), 5, 0.5);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

TEST(solvetrp, consistency_3) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd sol_1 = solvetrp(fe_space, u0_vector, 10, 1.0);
  Eigen::VectorXd sol_2 =
      solvetrp(fe_space, solvetrp(fe_space, u0_vector, 1, 0.1), 9, 0.9);

  EXPECT_NEAR((sol_1 - sol_2).norm(), 0.0, 1.0E-4);
}

// The final set of test cases checks the output of solverot and solvetrp
// against
// a numerical refernce solution.

TEST(solverot, reference_1) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd u_ref(u0_vector.size());
  u_ref << 0, 0, 0, 0, 0, 0.0159701, 0.0150475, 0, 0, 0.0171946, 0.0170436, 0,
      0, 0, 0, 0;
  Eigen::VectorXd u = solverot(fe_space, u0_vector, 1, 0.1);

  EXPECT_NEAR((u - u_ref).norm(), 0.0, 1.0E-5);
}

TEST(solverot, reference_2) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd u_ref(u0_vector.size());
  u_ref << 0, 0, 0, 0, 0, 0.00860437, 0.00827551, 0, 0, 0.00910899, 0.00901993,
      0, 0, 0, 0, 0;

  Eigen::VectorXd u = solverot(fe_space, u0_vector, 10, 0.1);

  EXPECT_NEAR((u - u_ref).norm(), 0.0, 1.0E-5);
}

TEST(solvetrp, reference_1) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd u_ref(u0_vector.size());
  u_ref << 0, 0, 0, 0, 0, 0.00993255, 0.00841276, 0, 0, 0.0105583, 0.00937419,
      0, 0, 0, 0, 0;
  Eigen::VectorXd u = solvetrp(fe_space, u0_vector, 1, 0.1);

  EXPECT_NEAR((u - u_ref).norm(), 0.0, 1.0E-5);
}

TEST(solvetrp, reference_2) {
  auto mesh_p = lf::mesh::test_utils::GenerateHybrid2DTestMesh(4, 1.0);
  auto fe_space =
      std::make_shared<const lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);

  auto u0_function = [](Eigen::VectorXd x) {
    return x(1) * x(0) * (1 - x(1)) * (1 - x(0));
  };
  Eigen::VectorXd u0_vector = lf::uscalfe::NodalProjection(
      *fe_space, lf::mesh::utils::MeshFunctionGlobal(u0_function));

  Eigen::VectorXd u_ref(u0_vector.size());
  u_ref << 0, 0, 0, 0, 0, 0.00670092, 0.00583407, 0, 0, 0.00707033, 0.00636342,
      0, 0, 0, 0, 0;
  Eigen::VectorXd u = solvetrp(fe_space, u0_vector, 10, 0.1);
  EXPECT_NEAR((u - u_ref).norm(), 0.0, 1.0E-5);
}

}  // namespace TranspSemiLagr::test