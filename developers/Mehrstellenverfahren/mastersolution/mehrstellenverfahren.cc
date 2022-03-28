/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Tobias Rohner
 * @ date 25-03-2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "mehrstellenverfahren.h"

#include <iomanip>
#include <iostream>
#include <vector>

namespace mehrstellenverfahren {

/* SAM_LISTING_BEGIN_1 */
Eigen::SparseMatrix<double> compMehrstellenA(unsigned int M) {
  // For the sake of efficiency the use of Eigen's sparse matrix data type is
  // essential. The matrix is stored in CCS format.
  Eigen::SparseMatrix<double> A(M * M, M * M);
#if SOLUTION
  // We already know that the matrix has at most $9M$ non-zero entries per row
  // and column. This information is passed to Eigen via the reserve() member
  // funtion.
  A.reserve(Eigen::VectorXi::Constant(M * M, 9));
  // Iterate over all interior nodes of the mesh and apply the stencil and
  // initialize the matrix in column-wise order, from top to bottom in every
  // column, which is most efficient for the CCS storage format.
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      // Index of the current node
      const int k = i * M + j;
      // Interaction term with the node below to the left
      if (i > 0 && j > 0) {
        A.insert(k - M - 1, k) = -1;
      }
      // Interaction term with the node below
      if (i > 0) {
        A.insert(k - M, k) = -4;
      }
      // Interaction term with the node below to the right
      if (i > 0 && j < M - 1) {
        A.insert(k - M + 1, k) = -1;
      }
      // Interaction term with the node to the left
      if (j > 0) {
        A.insert(k - 1, k) = -4;
      }
      // Interaction term with itself
      A.insert(k, k) = 20;
      // Interaction term with the node to the right
      if (j < M - 1) {
        A.insert(k + 1, k) = -4;
      }
      // Interaction term with the node above to the left
      if (i < M - 1 && j > 0) {
        A.insert(k + M - 1, k) = -1;
      }
      // Interaction term with the node above
      if (i < M - 1) {
        A.insert(k + M, k) = -4;
      }
      // Interaction term with the node above to th right
      if (i < M - 1 && j < M - 1) {
        A.insert(k + M + 1, k) = -1;
      }
    }
  }
  A.makeCompressed();
  A /= 6;
#else
  //====================
  // Your code goes here
  //====================
#endif
  return A;
}
/* SAM_LISTING_END_1 */

/* SAM_LISTING_BEGIN_2 */
double compgriderr(unsigned int M) {
  double err = 0;
#if SOLUTION
  // Obtain the Mehrstellen solution to the PDE
  const auto f = [](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y);
  };
  const Eigen::VectorXd mu = solveMehrstellen(f, M);
  // The exact solution u
  const auto u = [](double x, double y) {
    return std::sin(M_PI * x) * std::sin(M_PI * y) / (2 * M_PI * M_PI);
  };
  // Find the maximum difference between the discretized and exact solution at
  // the mesh nodes
  const double h = 1. / (M + 1);
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      const double x = (i + 1) * h;
      const double y = (j + 1) * h;
      // The index of the node at (i,j)
      const int k = i * M + j;
      // Update the error norm
      const double diff = std::fabs(u(x, y) - mu[k]);
      if (diff > err) {
        err = diff;
      }
    }
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  return err;
}
/* SAM_LISTING_END_2 */

/* SAM_LISTING_BEGIN_3 */
void tabulateMehrstellenError() {
  std::vector<double> errs;
  std::vector<unsigned int> Ms = {5, 10, 20, 40, 80, 160};
#if SOLUTION
  // Fill the errs vector with the errors at the different mesh resolutions
  for (unsigned int M : Ms) {
    errs.push_back(compgriderr(M));
  }
#else
  //====================
  // Your code goes here
  //====================
#endif
  // Print the collected data out in a table
  std::cout << "| M   | h       | error        |\n";
  std::cout << "-------------------------------|\n";
  for (int i = 0; i < Ms.size(); ++i) {
    const unsigned int M = Ms[i];
    const double h = 1. / (M + 1);
    const double err = errs[i];
    std::cout << "| " << std::setw(3) << M << " | " << std::fixed
              << std::setprecision(5) << h << " | " << std::scientific << err
              << " |\n";
  }
}
/* SAM_LISTING_END_3 */

}  // namespace mehrstellenverfahren
