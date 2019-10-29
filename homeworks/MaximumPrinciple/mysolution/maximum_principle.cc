/**
 * @file  maximum_principle_main.cc
 * @brief NPDE homework "MaximumPrinciple" code
 * @author Oliver Rietmann
 * @date 25.03.2019
 * @copyright Developed at ETH Zurich
 */

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace MaximumPrinciple {

/**
 * @brief Assembly
 *
 * Compute the global Galerkin matrix from the local
 * element matrix.
 *
 * @param M Number of interior vertices in x and y direction.
 * @param B_K Local element matrix.
 * @return Global Galerkin matrix of size M^2 times M^2.
 */
/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> assemble(int M, const Eigen::Matrix3d &B_K) {
  int M2 = M * M;
  Eigen::SparseMatrix<double> A(M2, M2);
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return A;
}

Eigen::SparseMatrix<double> computeGalerkinMatrix(int M, double c) {
  Eigen::Matrix3d B_K;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return assemble(M, B_K);
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_4 */
Eigen::SparseMatrix<double> computeGalerkinMatrixTR(int M, double c) {
  Eigen::Matrix3d B_K;
  /* SOLUTION_BEGIN */
  /* SOLUTION_END */
  return assemble(M, B_K);
}
/* SAM_LISTING_END_4 */

}  // namespace MaximumPrinciple
