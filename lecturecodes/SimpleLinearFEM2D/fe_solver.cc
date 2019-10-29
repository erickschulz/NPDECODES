/* **********************************************************************
   Demo code for course "Numerical Methods for PDEs
   Section "Case Study: Triangular Linear FEM in Two Dimensions
   ********************************************************************** */

#include "SimpleLinearFEM2D.h"
#include "local_assembler.h"

Eigen::VectorXd FESolver::Solve(TriaMesh2D const& mesh) {
  Eigen::SparseMatrix<double> A;
  // Initialize Galerkin matrix assembler
  // ElementMatrix_LaplMass_LFE as defined in local_assembler.cc is used to
  // locally compute elment matrices for the bilinear form \int_{K}
  // grad(b_i)grad(b_j) + b_i b_j dx
  if (FESolver::inefficient) {
    SlowMatrixAssembler mass_stiffness_matrix_assembler(
        &ElementMatrix_LaplMass_LFE);
    A = mass_stiffness_matrix_assembler.Assemble(mesh);
  } else {
    MatrixAssembler mass_stiffness_matrix_assembler(
        &ElementMatrix_LaplMass_LFE);
    A = mass_stiffness_matrix_assembler.Assemble(mesh);
  }

  // Initialize load vector assembler
  VectorAssembler load_vector_assembler(&localLoadLFE,
                                        FESolver::sourceFunction);
  Eigen::VectorXd phi = load_vector_assembler.Assemble(mesh);

  // Solve A \mu = phi for \mu to obtain the coefficients of the basis functions
  Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >
      solver;
  solver.analyzePattern(A);
  solver.factorize(A);
  Eigen::VectorXd mu = solver.solve(phi);

  return mu;
}
