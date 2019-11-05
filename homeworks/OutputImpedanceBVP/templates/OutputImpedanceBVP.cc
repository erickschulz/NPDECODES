/** @file
 * @brief NPDE OutputImpedanceBVP
 * @author Erick Schulz
 * @date 12/07/2019
 * @copyright Developed at ETH Zurich
 */

#include "OutputImpedanceBVP.h"

#include <cassert>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>

namespace OutputImpedanceBVP
{

Eigen::VectorXd solveImpedanceBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    Eigen::Vector2d g)
{
  // Related implementations:
  // Homework problem ErrorEstimatesForTraces:
  // https://gitlab.math.ethz.ch/ralfh/npdecodes/tree/master/homeworks/ErrorEstimatesForTraces
  Eigen::VectorXd discrete_solution;
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());
  // Obtain specification for shape functions on edges
  std::shared_ptr<const lf::uscalfe::ScalarReferenceFiniteElement<double>>
      rsf_edge_p = fe_space_p->ShapeFunctionLayout(lf::base::RefEl::kSegment());

  // I : ASSEMBLY
  // Matrix in triplet format holding Galerkin matrix, zero initially.
  lf::assemble::COOMatrix<double> A(N_dofs, N_dofs);
  // Right hand side vector, must be initialized with 0!
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();

  // I.i : Computing volume matrix for negative Laplace operator
  /* SOLUTION_BEGIN */

  // WRITE YOUR CODE HERE ...

  /* SOLUTION_END */

  // I.ii : Computing mass edge matrix resulting from Robin B.C.
  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  /* SOLUTION_BEGIN */

  // WRITE YOUR CODE HERE ...

  /* SOLUTION_END */

  // I.iii : Computing right-hand side vector
  // Right-hand side source function f
  auto mf_f =
      lf::uscalfe::MeshFunctionGlobal([](coord_t x) -> double { return 0.0; });
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  // I.iv : Imposing essential boundary conditions
  // Dirichlet data
  auto mf_g = lf::uscalfe::MeshFunctionGlobal(
      [&g](coord_t x) -> double { return g.dot(x); });
  /* SOLUTION_BEGIN */

  // WRITE YOUR CODE HERE ...

  /* SOLUTION_END */

  // Assembly completed! Convert COO matrix A into CRS format using Eigen's
  // internal conversion routines.
  Eigen::SparseMatrix<double> A_sparse = A.makeSparse();

  // II : SOLVING  THE LINEAR SYSTEM
  // II.i : Setting up Eigen's sparse direct elimination
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  solver.compute(A_sparse);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "LU decomposition failed");
  // II.ii : Solving
  discrete_solution = solver.solve(phi);
  LF_VERIFY_MSG(solver.info() == Eigen::Success, "Solving LSE failed");

  return discrete_solution;
};

double computeBoundaryOutputFunctional(
    const Eigen::VectorXd eta,
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p,
    Eigen::Vector2d d)
{
  double func_val = 0.0;
  // Pointer to current mesh
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space_p->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space_p->LocGlobMap()};

  // Obtain an array of boolean flags for the edges of the mesh, 'true'
  // indicates that the edge lies on the boundary
  auto bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};

  /* SOLUTION_BEGIN */

  // WRITE YOUR CODE HERE ...

  /* SOLUTION_END */

  // Computing value of the functional
  for (const lf::mesh::Entity &edge : mesh_p->Entities(1))
  {
    /* SOLUTION_BEGIN */

    // WRITE YOUR CODE HERE ...

    /* SOLUTION_END */
  }
  return func_val;
};

} // namespace OutputImpedanceBVP
