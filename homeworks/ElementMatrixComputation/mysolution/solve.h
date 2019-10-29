/**
 * @file
 * @brief NPDE homework ElementMatrixComputation code
 * @author Janik Schüttler
 * @date 06.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <lf/assemble/assemble.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/uscalfe/uscalfe.h>

#ifndef _SOLVEH_
#define _SOLVEH_

#include "../meshes/mesh.h"

#include "MyLinearFEElementMatrix.h"
#include "MyLinearLoadVector.h"

#endif

namespace ElementMatrixComputation {

/**
 * @brief      Given element builders, it creates FE space, assembles Galerkin
 * matrix and vector, and solves the LSE.
 *
 * @param[in]  mesh_p         A mesh pointer
 * @param[in]  elmat_builder  An element matrix builder object
 * @param[in]  elvec_builder  An element vector builder object
 *
 * @return     The solution vector
 */
template <class ELMAT_BUILDER, class ELVEC_BUILDER>
Eigen::VectorXd solve(ELMAT_BUILDER& elmat_builder,
                      ELVEC_BUILDER& elvec_builder) {
  // Use one of LehrFEM++'s default meshes. Try different meshes by changing the
  // function index parameter
  // std::shared_ptr<lf::mesh::Mesh> mesh_p =
  //     lf::mesh::test_utils::GenerateHybrid2DTestMesh(0, 1.0 / 3.0);
  auto mesh_p = Generate2DTestMesh();
  // We use a linear Lagrangian FE space
  auto fe_space =
      std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p);
  // Reference to current mesh, obtained from the FE space
  const lf::mesh::Mesh& mesh{*(fe_space->Mesh())};
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler& dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space`
  const lf::base::size_type N_dofs(dofh.NoDofs());

  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Invoke assembly on cells (co-dimension = 0). The element matrix builder is
  // passed as an argument
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
  Eigen::SparseMatrix<double> A_crs = A.makeSparse();

  // Right-hand side vector; has to be set to zero initially
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();
  // Invoke assembly on cells (codim == 0). The element vector builder is passed
  // as an argument
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // Define solution vector
  Eigen::VectorXd sol_vec = Eigen::VectorXd::Zero(N_dofs);

  /* BEGIN_SOLUTION */
  // Your implementation goes here!
  /* END_SOLUTION */

  double solver_error = (A_crs * sol_vec - phi).norm();
  double solver_relative_error = solver_error / phi.norm();

  // Postprocessing: Compute and output norms of the finite element solution
  // Helper class for L2 error computation
  lf::uscalfe::MeshFunctionL2NormDifference lc_L2(
      fe_space, lf::uscalfe::MeshFunctionConstant<double>(0.0), 2);
  // Helper class for H1 semi norm
  lf::uscalfe::MeshFunctionL2GradientDifference lc_H1(
      fe_space,
      lf::uscalfe::MeshFunctionConstant<Eigen::Vector2d>(
          Eigen::Vector2d(0.0, 0.0)),
      2);

  double L2norm = lf::uscalfe::NormOfDifference(dofh, lc_L2, sol_vec);
  double H1snorm = lf::uscalfe::NormOfDifference(dofh, lc_H1, sol_vec);
  std::cout << "Error of solver: Absolute error = " << solver_error
            << ", relative error = " << solver_relative_error << std::endl;
  std::cout << "Norms of FE solution: L2-norm = " << L2norm
            << ", H1-seminorm = " << H1snorm << std::endl;
  return sol_vec;
}

/**
 * @brief      Right hand side functional `f(x)`
 *
 * @param[in]  x     A 2D vector that provides the evaluation of the source
 *
 * @return     The value of `f` evaluated at point `x`
 */
double f(Eigen::Vector2d x) { return 1 + x(0) * x(0) + x(1) * x(1); };

/**
 * @brief Solve Poisson's equation -△u = f using LehrFEM++'s built in element
 * matrix and vector builders
 *
 * @return     The solution vector
 */
Eigen::VectorXd solvePoissonBVP() {
  // Convert tPoissonda function f to a LehrFEM++ mesh function object
  lf::uscalfe::MeshFunctionGlobal mf_f{f};

  // Define the solution vector
  Eigen::VectorXd solution;

  /* SAM_LISTING_BEGIN_2 */
  /* BEGIN_SOLUTION */
  // Your implementation goes here!
  /* END_SOLUTION */
  /* SAM_LISTING_END_2 */

  return solution;
}

/**
 * @brief Solve Neumann equation -△u + u = f where ɑ is a constant diffusion
 * coefficient using a custom implementation of element matrix and element
 * vector builders
 *
 * @return     The solution vector
 */
Eigen::VectorXd solveNeumannEq() {
  // Define the solution vector
  Eigen::VectorXd solution;

  /* BEGIN_SOLUTION */
  // Your implementation goes here!
  /* END_SOLUTION */

  return solution;
}
}  // namespace ElementMatrixComputation
