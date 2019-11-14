/**
 * @file
 * @brief NPDE homework CoupledSecondOrderBVP
 * @author Erick Schulz
 * @date 13/11/2019
 * @copyright Developed at ETH Zurich
 */

#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/mesh/utils/utils.h>
#include <lf/refinement/refinement.h>
#include <lf/uscalfe/uscalfe.h>

namespace CoupledSecondOrderBVP {

// Function solving the coupled BVP
template <typename FUNCTOR>
Eigen::VectorXd solveCoupledBVP(
    std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space,
    double gamma, FUNCTOR &&f) {
  Eigen::VectorXd sol_vec;  // solution vector
  // Get pointer to current f32divf32x
  std::shared_ptr<const lf::mesh::Mesh> mesh_p = fe_space->Mesh();
  // Obtain local->global index mapping for current finite element space
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()};
  // Dimension of finite element space
  const lf::uscalfe::size_type N_dofs(dofh.NumDofs());

  // Data related to Dirichlet B.C.
  // Obtain an array of boolean flags for the nodes of the mesh, 'true'
  // indicates that the node lies on the boundary
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 0)};
  // Inspired by the example in the documentation of
  // InitEssentialConditionFromFunction()
  // https://craffael.github.io/lehrfempp/namespacelf_1_1uscalfe.html#a5afbd94919f0382cf3fb200c452797ac
  // Creating a predicate that will guarantee that the computations are carried
  // only on the exterior boundary nodes of the mesh using the boundary flags
  auto nodes_predicate_Dirichlet =
      [&nodes_bd_flags](const lf::mesh::Entity &node) -> bool {
    return nodes_bd_flags(node);
  };

  /* I : Creating coefficients as Lehrfem++ mesh functions */
  // Coefficients used in the class template
  // ReactionDiffusionElementMatrixProvider<SCALAR,DIFF_COEFF,REACTION_COEFF>
  auto const_one = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
  auto const_zero = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });

  /* II: Instantiating finite element matrices */
  // Matrices in triplet format holding  the Galerkin matrices
  // (set to zero initially)
  lf::assemble::COOMatrix<double> A0(N_dofs, N_dofs);  // upper left block
  lf::assemble::COOMatrix<double> A1(N_dofs, N_dofs);  // lower right block
  lf::assemble::COOMatrix<double> M(N_dofs, N_dofs);   // off-diag blocks
  lf::assemble::COOMatrix<double> L(2 * N_dofs, 2 * N_dofs);  // full matrix
  // Right hand side vector
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi(N_dofs);
  phi.setZero();  // has to be zero initially

  /* III : Computing the Galerkin matrices */
  // III.i Computing A0 : standard negative Laplace Galerkin Matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_zero)>
      A0_builder(fe_space, const_one, const_zero);
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A0_builder, A0);
  // III.ii Computing A01 : Laplacian with reaction term
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_one), decltype(const_one)>
      A1_builder(fe_space, const_one, const_one);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, A1_builder, A1);
  // III.iii Computing M : standard mass matrix
  lf::uscalfe::ReactionDiffusionElementMatrixProvider<
      double, decltype(const_zero), decltype(const_one)>
      M_builder(fe_space, const_zero, const_one);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, M_builder, M);
  // III.iv : Computing right-hand side vector
  // Wrap the right-hand side lambda source function f in a Lehrfem++ MeshFunction
  auto mf_f =
      lf::uscalfe::MeshFunctionGlobal(f);
  lf::uscalfe::ScalarLoadElementVectorProvider<double, decltype(mf_f)>
      elvec_builder(fe_space_p, mf_f);
  // Invoke assembly on cells (codim == 0)
  AssembleVectorLocally(0, dofh, elvec_builder, phi);

  /* IV : Assemble the full linear system matrix and right hand side vector */
  //                     _        _
  //    L (u  p)^T  :=  |  A0    M | (u  p)^T  = (f 0)^T
  //                    |_ M^T  A1_|
  //
  // IV.i Inserting A0 in L
  const std::vector<Eigen::Triplet<double>> A0_triplets_vec = A0.triplets();
  for (auto &triplet : A0_triplets_vec) {
    L.AddToEntry(triplet.row(), triplet.col(), triplet.value());
  }
  // IV.ii Inserting A1 in L
  const std::vector<Eigen::Triplet<double>> A1_triplets_vec = A1.triplets();
  for (auto &triplet : A1_triplets_vec) {
    L.AddToEntry(triplet.row() + N_dofs, triplet.col() + N_dofs,
                 triplet.value());
  }
  // IV.iii Inserting M in L
  // Notice the symmetry between the entries. The two indices of the entry
  // are flipped so that it is the tranpose of M that is added to the left lower
  // diagonal block of L.
  const std::vector<Eigen::Triplet<double>> M_triplets_vec = M.triplets();
  for (auto &triplet : M_triplets_vec) {
    L.AddToEntry(triplet.row(), triplet.col() + N_dofs,
                 triplet.value());  // for M in upper right block
    L.AddToEntry(triplet.col() + N_dofs, triplet.row(),
                 triplet.value());  // for M^T in lower left block
  }

  return sol_vec;
}

}  // namespace CoupledSecondOrderBVP
