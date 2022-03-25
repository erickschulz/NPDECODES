/**
 * @ file
 * @ brief NPDE homework TEMPLATE MAIN FILE
 * @ author Tobias Rohner
 * @ date 25-03-2022
 * @ copyright Developed at SAM, ETH Zurich
 */

#include "mehrstellenverfahren.h"

#include <cstdio>
#include <vector>

namespace mehrstellenverfahren {

Eigen::SparseMatrix<double> compMehrstellenA(unsigned int M) {
  Eigen::SparseMatrix<double> A(M * M, M * M);
#if SOLUTION
  // We have at most 9nnz per row of A and A has M^2 rows
  const int nnz = 9 * M * M;
  A.reserve(nnz);
  // Iterate over all interior nodes of the mesh and apply the stencil
  for (int i = 0; i < M; ++i) {
    for (int j = 0; j < M; ++j) {
      // Index of the current node
      const int k = i * M + j;
      // Interaction term with the node below to the left
      if (i > 0 && j > 0) {
        A.insert(k, k - M - 1) = -1;
      }
      // Interaction term with the node below
      if (i > 0) {
        A.insert(k, k - M) = -4;
      }
      // Interaction term with the node below to the right
      if (i > 0 && j < M - 1) {
        A.insert(k, k - M + 1) = -1;
      }
      // Interaction term with the node to the left
      if (j > 0) {
        A.insert(k, k - 1) = -4;
      }
      // Interaction term with itself
      A.insert(k, k) = 20;
      // Interaction term with the node to the right
      if (j < M - 1) {
        A.insert(k, k + 1) = -4;
      }
      // Interaction term with the node above to the left
      if (i < M - 1 && j > 0) {
        A.insert(k, k + M - 1) = -1;
      }
      // Interaction term with the node above
      if (i < M - 1) {
        A.insert(k, k + M) = -4;
      }
      // Interaction term with the node above to th right
      if (i < M - 1 && j < M - 1) {
        A.insert(k, k + M + 1) = -1;
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
  printf("| M   | h       | error        |\n");
  printf("-------------------------------|\n");
  for (int i = 0; i < Ms.size(); ++i) {
    const unsigned int M = Ms[i];
    const double h = 1. / (M + 1);
    const double err = errs[i];
    printf("| %3d | %.05f | %e |\n", M, h, err);
  }
}

}  // namespace mehrstellenverfahren
