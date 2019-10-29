/**
 * @file
 * @brief NPDE homework TransformationOfGalerkinMatrices template
 * @author Erick Schulz
 * @date 01/03/2019
 * @copyright Developed at ETH Zurich
 */

#include <boost/assert.hpp>

#include "trans_gal_mat.h"

namespace TransformationOfGalerkinMatrices {

std::vector<triplet_t> transformCOOmatrix(const std::vector<triplet_t> &A) {
  std::vector<triplet_t> A_tilde{};  // return value

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
  /* SAM_LISTING_BEGIN_1 */

  /* SAM_LISTING_END_1 */
  /* SOLUTION_END */
  return A_tilde;
}

}  // namespace TransformationOfGalerkinMatrices
