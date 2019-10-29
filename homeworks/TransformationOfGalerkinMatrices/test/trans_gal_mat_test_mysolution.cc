#include <gtest/gtest.h>

#include "../mysolution/trans_gal_mat.h"

namespace TransformationOfGalerkinMatrices::test {

TEST(TransformationOfGalerkinMatrices, TestTransformation) {
  std::cout << "NPDE homework TransformationOfGalerkinMatrices: unit test"
            << std::endl;

  // Create change of basis matrix S
  typedef Eigen::SparseMatrix<double> SpMat;
  std::vector<triplet_t> S_triplets;

  int N = 4;
  SpMat S(2 * N, 2 * N);
  for (int i = 0; i < N; i++) {
    S_triplets.push_back(triplet_t(i, i * 2, 1));
    S_triplets.push_back(triplet_t(i, i * 2 + 1, 1));
    S_triplets.push_back(triplet_t(N + i, i * 2, 1));
    S_triplets.push_back(triplet_t(N + i, i * 2 + 1, -1));
  }
  S.setFromTriplets(S_triplets.begin(), S_triplets.end());

  SpMat sol_mat(2 * N, 2 * N);
  SpMat A_mat(2 * N, 2 * N);
  SpMat A_tilde_mat(2 * N, 2 * N);
  std::vector<triplet_t> A, A_tilde;

  // Case 1
  A.push_back(triplet_t(2 * N - 1, 2 * N - 1, 1));  // i, j even
  A.push_back(triplet_t(3 - 1, 3 - 1, 5));          // i, j odd
  A_mat.setFromTriplets(A.begin(), A.end());
  A_tilde = TransformationOfGalerkinMatrices::transformCOOmatrix(A);
  A_tilde_mat.setFromTriplets(A_tilde.begin(), A_tilde.end());
  sol_mat = S * A_mat * S.transpose();

  ASSERT_TRUE((sol_mat - A_tilde_mat).norm() == 0);

  // Case 2
  A.clear();
  A_tilde.clear();
  A.push_back(triplet_t(2 * N - 1, 2 * N - 1, 1));
  A.push_back(triplet_t(4 - 1, 3 - 1, 5));  // i even, j odd
  A_mat.setFromTriplets(A.begin(), A.end());
  A_tilde = TransformationOfGalerkinMatrices::transformCOOmatrix(A);
  A_tilde_mat.setFromTriplets(A_tilde.begin(), A_tilde.end());
  sol_mat = S * A_mat * S.transpose();

  ASSERT_TRUE((sol_mat - A_tilde_mat).norm() == 0);

  // Case 3
  A.clear();
  A_tilde.clear();
  A.push_back(triplet_t(2 * N - 1, 2 * N - 1, 1));
  A.push_back(triplet_t(3 - 1, 4 - 1, 5));  // i odd, j even
  A_mat.setFromTriplets(A.begin(), A.end());
  A_tilde = TransformationOfGalerkinMatrices::transformCOOmatrix(A);
  A_tilde_mat.setFromTriplets(A_tilde.begin(), A_tilde.end());
  sol_mat = S * A_mat * S.transpose();

  ASSERT_TRUE((sol_mat - A_tilde_mat).norm() == 0);
}

}  // namespace TransformationOfGalerkinMatrices::test
