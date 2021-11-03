/**
 * @file transpsemilagr.h
 * @brief NPDE homework TranspSemiLagr code
 * @author Philippe Peter
 * @date November 2020
 * @copyright Developed at SAM, ETH Zurich
 */
#include <lf/assemble/assemble.h>
#include <lf/base/base.h>
#include <lf/fe/fe.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <memory>

#include "local_assembly.h"

namespace TranspSemiLagr {

/**
 * @brief encorce zero dirichlet boundary conditions
 * @param fe_space shared point to the fe spce on which the problem is defined.
 * @param A reference to the square coefficient matrix in COO format
 * @param b reference to the right-hand-side vector
 */
void enforce_zero_boundary_conditions(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    lf::assemble::COOMatrix<double>& A, Eigen::VectorXd& b) {
  lf::mesh::utils::MeshFunctionGlobal mf_zero{
      [](const Eigen::Vector2d& /*x*/) { return 0.0; }};
  const lf::fe::ScalarReferenceFiniteElement<double>* rsf_edge_p =
      fe_space->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(fe_space->Mesh(), 1)};
  auto flag_values{
      lf::fe::InitEssentialConditionFromFunction(*fe_space, bd_flags, mf_zero)};

  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&flag_values](lf::assemble::glb_idx_t dof_idx) {
        return flag_values[dof_idx];
      },
      A, b);
}

/**
 * @brief performs a semi lagrangian step according to the update
 * formula 7.3.4.13
 * @param fe_space finite element space on which the problem is solved
 * @param u0_vector nodal values of the solution at the previous time step
 * @param v velocity field (time independent)
 * @param tau time step size
 * @return nodal values of the approximated solution at current time step
 */
template <typename FUNCTOR>
Eigen::VectorXd semiLagr_step(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const Eigen::VectorXd& u0_vector, FUNCTOR v, double tau) {
  // Assemble left hand side A = A_lm + tau*A_s
  lf::assemble::COOMatrix<double> A(fe_space->LocGlobMap().NumDofs(),
                                    fe_space->LocGlobMap().NumDofs());

  // stiffness matrix tau*A_s
  lf::uscalfe::ReactionDiffusionElementMatrixProvider
      stiffness_element_matrix_provider(
          fe_space, lf::mesh::utils::MeshFunctionConstant(tau),
          lf::mesh::utils::MeshFunctionConstant(0.0));
  lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(),
                                      fe_space->LocGlobMap(),
                                      stiffness_element_matrix_provider, A);

  // lumped mass matrix A_lm
  LumpedMassElementMatrixProvider lumped_mass_element_matrix_provider(
      [](Eigen::Vector2d /*x*/) { return 1.0; });
  lf::assemble::AssembleMatrixLocally(0, fe_space->LocGlobMap(),
                                      fe_space->LocGlobMap(),
                                      lumped_mass_element_matrix_provider, A);

  // warp u0 into a mesh function (required by the Vector provider) & assemble
  // rhs.
  auto u0_mf = lf::fe::MeshFunctionFE(fe_space, u0_vector);
  UpwindLagrangianElementVectorProvider vector_provider(
      v, tau, fe_space->Mesh(), u0_mf);
  Eigen::VectorXd b(fe_space->LocGlobMap().NumDofs());
  b.setZero();
  lf::assemble::AssembleVectorLocally(0, fe_space->LocGlobMap(),
                                      vector_provider, b);

  enforce_zero_boundary_conditions(fe_space, A, b);

  // solve LSE
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  return solver.solve(b);
}

/**
 * @brief approximates the solution to the first model problem specified in the
 * exercise sheet based on N uniform timesteps of the Semi Lagrangian method
 * @param fe_space (linear) finite element space on which the solution is
 * approximated
 * @param u0_vector vector of nodal values of the initial condition
 * @param N number of time steps
 * @param T final time
 */
Eigen::VectorXd solverot(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::VectorXd u0_vector, int N, double T) {
  // time step size
  double tau = T / N;

  // velocity field of the model problem
  auto v = [](Eigen::Vector2d x) {
    return (Eigen::Vector2d() << -x(1) + 3.0 * x(0) * x(0), x(0)).finished();
  };

  // approximate solution based on N uniform time steps.
  for (int i = 0; i < N; ++i) {
    u0_vector = semiLagr_step(fe_space, u0_vector, v, tau);
  }

  return u0_vector;
}

/**
 * @param solves the variational evolution problem specified in the exercise
 * sheet over on time step.
 * @param fe_space finite element space on which the variational evolution
 * problem is solved
 * @param u0_vector nodal values of the solution at the previous time step
 * @param c coefficient function in the variational evolution problem
 * @param tau time step size
 */
template <typename FUNCTOR>
Eigen::VectorXd reaction_step(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    const Eigen::VectorXd& u0_vector, FUNCTOR c, double tau) {
  // The explicit midpoint rule computes the solution to y' = f(y) as
  // y* = y_n + tau/2*f(y_n)
  // y_{n+1} =   y_n + tau*f(y*)

  // In the evolution problem f(y) solves the discretized variational problem
  // A_m*f(y) = A_{m,c}*y
  // where A_m is the lumped mass matrix
  // and A_{m,c} is the lumped mass matrix with coefficient function c.

  // Assemble matrix A_{m,c}
  LumpedMassElementMatrixProvider mass_element_matrix_provider_c(c);
  lf::assemble::COOMatrix<double> mass_matrix_c(
      fe_space->LocGlobMap().NumDofs(), fe_space->LocGlobMap().NumDofs());
  lf::assemble::AssembleMatrixLocally(
      0, fe_space->LocGlobMap(), fe_space->LocGlobMap(),
      mass_element_matrix_provider_c, mass_matrix_c);
  Eigen::SparseMatrix<double> mass_matrix_c_sparse = mass_matrix_c.makeSparse();

  // Assemble matrix A_m
  LumpedMassElementMatrixProvider mass_element_matrix_provider_1(
      [](Eigen::Vector2d /*x*/) { return 1.0; });
  lf::assemble::COOMatrix<double> mass_matrix_1(
      fe_space->LocGlobMap().NumDofs(), fe_space->LocGlobMap().NumDofs());
  lf::assemble::AssembleMatrixLocally(
      0, fe_space->LocGlobMap(), fe_space->LocGlobMap(),
      mass_element_matrix_provider_1, mass_matrix_1);
  Eigen::SparseMatrix<double> mass_matrix_1_sparse = mass_matrix_1.makeSparse();

  // initialize sparse solver for mass matrix A_m
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(mass_matrix_1_sparse);

  // evaluate f(y_n)
  Eigen::VectorXd k1 = solver.solve(mass_matrix_c_sparse * u0_vector);

  // evaluate f(y*)
  Eigen::VectorXd k2 =
      solver.solve(mass_matrix_c_sparse * (u0_vector + 0.5 * tau * k1));

  // update to y_{n+1}
  return u0_vector + tau * k2;
}

/**
 * @brief approximates the solution to the second model problem specified in the
 * exercise sheet based on N uniform time steps of the Strang-splitting
 * split-step method
 * @param fe_space (linear) finite element space on which the solution is
 * approximated
 * @param u0_vector vector of nodal values of the initial condition
 * @param N number of time steps
 * @param T final time
 */
Eigen::VectorXd solvetrp(
    std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space,
    Eigen::VectorXd u0_vector, int N, double T) {
  // time step size
  double tau = T / N;

  // coefficient functiosn for the semiLagr and reaction step
  auto v = [](Eigen::Vector2d x) {
    return (Eigen::Vector2d() << -x(1) + 3.0 * x(0) * x(0), x(0)).finished();
  };
  auto c = [](Eigen::Vector2d x) { return -6.0 * x(0); };

  // Strang splitting scheme
  //-----------------------
  // first SemiLagr half step:
  u0_vector = semiLagr_step(fe_space, u0_vector, v, 0.5 * tau);

  // intermediate time steps: Combine two semiLagr half steps to one step
  for (int i = 0; i < N - 1; ++i) {
    u0_vector = reaction_step(fe_space, u0_vector, c, tau);
    u0_vector = semiLagr_step(fe_space, u0_vector, v, tau);
  }

  // final reaction step and semiLagr half step
  u0_vector = reaction_step(fe_space, u0_vector, c, tau);
  u0_vector = semiLagr_step(fe_space, u0_vector, v, tau * 0.5);

  return u0_vector;

}

}  // namespace TranspSemiLagr
