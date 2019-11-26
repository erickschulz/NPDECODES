/** @file electrostaticforce.cc
 * @brief NPDE homework ElectrostaticForce code
 * @author Erick Schulz
 * @date 25/11/2019
 * @copyright Developed at ETH Zurich */

#include "electrostaticforce.h"

#include <math.h>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

namespace ElectrostaticForce {

Eigen::Vector2d computeExactForce() {
  Eigen::Vector2d force;  // return vector

  unsigned int N = 1e6;  // nb. of quadrature points

  // Data
  double r = 4.0 / 15.0;  // radius of the circle
  Eigen::Vector2d a(-16.0 / 15.0, 0.0);
  Eigen::Vector2d b(-1.0 / 15.0, 0.0);

  // Tools
  Eigen::Vector2d x;     // euclidean coordinates
  Eigen::Vector2d n;     // normal vector at the coordinates
  Eigen::Vector2d grad;  // gradient of the exact solution
  double angle;          // interior angle

  /* Compute the force using the trapezoidal rule */
  force << 0.0, 0.0;  // initially zero
  for (unsigned int i = 0; i < N; i++) {
    // I. Find the euclidean coordinates of integration node
    angle = i * (2 * M_PI) / N;
    x(0) = r * std::cos(angle);
    x(1) = r * std::sin(angle);
    // II. Evaluate the gradient at the coordinates
    grad = ((x - a) / (x - a).squaredNorm() - (x - b) / (x - b).squaredNorm()) /
           std::log(2.0);
    // III. Evaluate the normal vector at the coordinates
    n = -x / x.norm();
    // IV. Evaluate full integrand
    force += grad.dot(n) * grad;
  }
  // V. Scale the summation
  force *= r * M_PI / N;

  return force;
}

Eigen::VectorXd solvePoissonBVP(
    const std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> &fe_space_p) {
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
  // Computing volume matrix for negative Laplace operator
#if SOLUTION
  // Initialize object taking care of local mass (volume) computations.
  lf::uscalfe::LinearFELaplaceElementMatrix elmat_builder{};
  // Invoke assembly on cells (co-dimension = 0 as first argument)
  // Information about the mesh and the local-to-global map is passed through
  // a Dofhandler object, argument 'dofh'. This function call adds triplets to
  // the internal COO-format representation of the sparse matrix A.
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, elmat_builder, A);
#else
  //====================
  // Your code goes here
  //====================
#endif

  // II. IMPOSING ESSENTIAL BOUNDARY CONDITIONS
  // Creating Dirichlet data as mesh functions
  auto mf_zero_f = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 0.0; });
  auto mf_one_f = lf::uscalfe::MeshFunctionGlobal(
      [](Eigen::Vector2d x) -> double { return 1.0; });
#if SOLUTION
  // Obtain arrays of boolean flags for the edges and nodes of the mesh, 'true'
  // indicates that the edge or node lies on the boundary
  auto edges_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 1)};
  auto nodes_bd_flags{lf::mesh::utils::flagEntitiesOnBoundary(mesh_p, 2)};
  // Inspired by the example in the documentation of
  // InitEssentialConditionFromFunction()
  // https://craffael.github.io/lehrfempp/namespacelf_1_1uscalfe.html#a5afbd94919f0382cf3fb200c452797ac
  // Creating a predicate that will guarantee that the computations are carried
  // only on the proper edges of the mesh using the boundary flags
  auto edges_predicate_interior_bd =
      [&edges_bd_flags](const lf::mesh::Entity &edge) -> bool {
    if (edges_bd_flags(edge)) {
      auto endpoints = lf::geometry::Corners(*(edge.Geometry()));
      if (endpoints.col(0).norm() > 2.7) {
        return false;
      }
    }
    return true;
  };
  // II.i Imposing interior Dirichlet boundary condition
  // Determine the fixed dofs on the interior boundary and their values
  // Alternative: See lecturedemoDirichlet() in
  // https://github.com/craffael/lehrfempp/blob/master/examples/lecturedemos/lecturedemoassemble.cc
  auto edges_flag_interior_values{
      lf::uscalfe::InitEssentialConditionFromFunction(
          dofh, *rsf_edge_p, edges_predicate_interior_bd, mf_one_f)};
  // Eliminate interior Dirichlet dofs from the linear system
  lf::assemble::FixFlaggedSolutionCompAlt<double>(
      [&edges_flag_interior_values](lf::assemble::glb_idx_t dof_idx) {
        return edges_flag_interior_values[dof_idx];
      },
      A, phi);
  // II.ii Imposing exterior Dirichlet boundary condition
  // Index predicate for the selectvals FUNCTOR of dropMatrixRowsAndColumns
  auto nodes_exterior_bd_selector = [&nodes_bd_flags, &dofh](unsigned int idx) -> bool {
    if (nodes_bd_flags(dofh.Entity(idx))) {
      auto coordinates = lf::geometry::Corners(*(dofh.Entity(idx).Geometry()));
      if (coordinates.norm() > 0.27) {
        return true;
      }
    }
    return false;
  };
  dropMatrixRowsAndColumns(nodes_exterior_bd_selector, A);
  // Assigning zero to the exterior boundary values of phi
  for (unsigned int dof_idx = 0; dof_idx < N_dofs; ++dof_idx){
  if (nodes_exterior_bd_selector(dof_idx)){
      phi(dof_idx) = 0.0;
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif

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

}  // namespace ElectrostaticForce
