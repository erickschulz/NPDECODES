/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices code
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <boost/assert.hpp>

#include "trans_gal_mat.h"

namespace TransformationOfGalerkinMatrices {

std::vector<triplet_t> transformCOOmatrix(const std::vector<triplet_t> &A) {
  std::vector<triplet_t> A_t{};  // return value

  // First step: find the size of the matrix by searching the maximal
  // indices. Depends on the assumption that no zero rows/columns occur.
  int rows_max_idx = 0, cols_max_idx = 0;
  for (const triplet_t &triplet : A) {
    rows_max_idx =
        (triplet.row() > rows_max_idx) ? triplet.row() : rows_max_idx;
    cols_max_idx =
        (triplet.col() > cols_max_idx) ? triplet.col() : cols_max_idx;
  }
  int n_rows = rows_max_idx + 1;
  int n_cols = cols_max_idx + 1;

  BOOST_ASSERT_MSG(n_rows == n_cols, "Matrix must be square");
  BOOST_ASSERT_MSG(n_cols % 2 == 0, "Matrix dimension must be even");

  int N = n_cols;      // Size of (square) matrix
  int M = n_cols / 2;  // Half the size

  /* SOLUTION_BEGIN */
  // clang-format off
  /* SAM_LISTING_BEGIN_1 */
  // Distribute entries of "old" matrix to new matrix
  for (const triplet_t &it : A) {
    // row and column indices of current A triplet
    const int I = it.row() + 1; // $\cob{k}$ in \prbeqref{D1}--\prbeqref{D4}
    const int J = it.col() + 1; // $\cob{\ell}$ in \prbeqref{D1}--\prbeqref{D4}

    // Distinguish different parities of indices
    if (I % 2 == 0 & J % 2 == 0) {
      // even, even, \prbeqref{D1}
      A_t.emplace_back(I / 2 - 1, J / 2 - 1, it.value());
      A_t.emplace_back(I / 2 + M - 1, J / 2 + M - 1, it.value());
      A_t.emplace_back(I / 2 + M - 1, J / 2 - 1, -it.value());
      A_t.emplace_back(I / 2 - 1, J / 2 + M - 1, -it.value());
    } else if (I % 2 != 0 & J % 2 != 0) {
      // odd, odd, see \prbeqref{D2}
      A_t.emplace_back((I + 1) / 2 - 1, (J + 1) / 2 - 1, it.value());
      A_t.emplace_back((I + 1) / 2 + M - 1, (J + 1) / 2 + M - 1, it.value());
      A_t.emplace_back((I + 1) / 2 - 1, (J + 1) / 2 + M - 1, it.value());
      A_t.emplace_back((I + 1) / 2 + M - 1, (J + 1) / 2 - 1, it.value());
    } else if (I % 2 == 0 & J % 2 != 0) {
      // even, odd, see \prbeqref{D3}
      A_t.emplace_back(I / 2 - 1, (J + 1) / 2 - 1, it.value());
      A_t.emplace_back(I / 2 - 1, (J + 1) / 2 + M - 1, it.value());
      A_t.emplace_back(I / 2 + M - 1, (J + 1) / 2 + M - 1, -it.value());
      A_t.emplace_back(I / 2 + M - 1, (J + 1) / 2 - 1, -it.value());
    } else if (I % 2 != 0 & J % 2 == 0) {
      // odd, even, see \prbeqref{D4}
      A_t.emplace_back((I + 1) / 2 - 1, J / 2 - 1, it.value());
      A_t.emplace_back((I + 1) / 2 + M - 1, J / 2 + M - 1, -it.value());
      A_t.emplace_back((I + 1) / 2 + M - 1, J / 2 - 1, it.value());
      A_t.emplace_back((I + 1) / 2 - 1, J / 2 + M - 1, -it.value());
    } else {
      BOOST_ASSERT_MSG(false, "Should never get there!");
    }
  }
  /* SAM_LISTING_END_1 */
  // clang-format on
  /* SOLUTION_END */
  return A_t;
}

}  // namespace TransformationOfGalerkinMatrices
